from __future__ import division
import warnings
import sys
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen_repo') # (If one does not want to change $PYTHONPATH)
import unittest
import pymatgen
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Vasprun
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import *
from pymatgen.core.structure_modifier import StructureEditor
import numpy as np
import os
import fit_elas

__author__="Maarten de Jong"
__copyright__ = "Copyright 2011, The Materials Project"
__credits__ = "Mark Asta"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ ="Jan 24, 2012"

np.set_printoptions(precision=3)			# pretty printing of numpy matrices
np.set_printoptions(suppress=True)			# pretty printing of numpy matrices

class DeformGeometry(object):
	"""
	A class which takes a CIF-file, generates a pymatgen core structure 
	and applies a range of lattice deformations. The deformed lattices 
	are to be used as inputs to an ab initio code such as VASP to compute 
	the stress tensor and from this, the full elastic tensor of rank 4
	"""

	def __init__(self, path_to_cif):
		"""
		Constructor to initialize path to CIF-file
		
		Args: 
			path_to_cif - string with full path to CIF-file, e.g.
			'/home/MDEJONG1/pythonplayground/pymatgen/classes/7048.cif'
		"""
		self._path = path_to_cif
	
	def CIF2struct(self):
		"""
		Parse CIF to pymatgen core structure for easy manipulation
		
		Args: 
			none
		"""
		myCIF = CifParser(self._path).get_structures()[0]
		self.base_struct = myCIF

		return self.base_struct

	def get_residual_stress(self, residual_stress=np.zeros((3,3)), res_tresh_hold=2):
		"""
	    Parse the residual stress tensor 

		Args:
			residual_stress - 3x3 numpy matrix containing residual stresses, 
			defaults to zero stress tensor
			res_tresh_hold - what residual stresses (in kB) are acceptable
		"""		
		self._residual_stress = residual_stress
		
		if (np.abs(residual_stress)>res_tresh_hold).any() == True:	
			warnings.warn('***Residual stresses are quite large, structural relaxation might not have converged...***')
	
		return self._residual_stress

	def deform(self, nd=0.02, ns=0.02, m=6, n=6): # this class still requires some (optional) user input
		"""
		Take original geometry (from CIF-file) and apply a range of deformations. 
		Default values generally work well for metals (in my experience). However, one might need to
		make changes for material such as oxides.

		Note: deformed structures are stored in a dict, named defstructures. Also stored are
		deformation gradient tensor, strain tensor, and strain tensor-indices.
		defstructures[0] contains first deformed structure and its properties, defstructures[1]
		contains the second etc.

		Args: 
			nd - maximum amount of normal strain  
			ns - maximum amount of shear strain
			m - number of normal deformations used for structural deformations, even integer required
			n - number of shear deformations used for structural deformations, even integer required
		"""	
		self.msteps = np.int(m)
		self.nsteps = np.int(n)

		if m%2!=0:
			raise ValueError("m has to be even.")

		if n%2!=0:
			raise ValueError("n has to be even.")

		mystrains = np.zeros((3, 3, np.int(m*3) + np.int(n*3)))

		defs = np.linspace(-nd, nd, num=m+1)
		defs = np.delete(defs, np.int(m/2), 0)
		sheardef = np.linspace(-ns, ns, num=n+1)
		sheardef = np.delete(sheardef, np.int(n/2), 0)
	
		defstructures = dict()

		counter = 0
		# First apply non-shear deformations
		for i1 in range(0, 3):

			for i2 in range(0, len(defs)):

				s = StructureEditor(self.base_struct)
				F = np.identity(3)
				F[i1, i1] = F[i1, i1] + defs[i2]
				E = 0.5*(np.transpose(F)*F-np.eye(3))      # Green-Lagrange strain tensor
				s.apply_strain_transformation(F)           # let deformation gradient tensor act on undistorted lattice
				defstructures[counter] = [s.modified_structure, F, E, (i1, i1)]
				counter += 1

		# Now apply shear deformations #		
		F_index = [[0, 1], [0, 2], [1, 2]]
		for j1 in range(0, 3):
		
			for j2 in range(0, len(sheardef)):

				s = StructureEditor(self.base_struct)
				F = np.identity(3)
				F[F_index[j1][0], F_index[j1][1]] = F[F_index[j1][0], F_index[j1][1]] + sheardef[j2]
				F = np.matrix(F)						   # Not sure why, but without this statement, calculation of strain tensor gets messed up...
				E = 0.5*(np.transpose(F)*F-np.eye(3))      # Green-Lagrange strain tensor
				s.apply_strain_transformation(F)           # let deformation gradient tensor act on undistorted lattice
				defstructures[counter] = [s.modified_structure, F, E, (F_index[j1][0], F_index[j1][1])]
				counter += 1
		
		self.defstructures = defstructures

		return self.defstructures

	def append_stress_tensors(self, stress_tensor_dict):
		"""
		After the ab initio engine is run, call this method to append the stress tensors to the corresponding structures,
		stored in defstructures. Residual stresses are subtracted out for increased accuracy.

		Args:
			stress_tensor_dict - a dict with  3x3 numpy matrices, containing the stresses. Key should be structure number,
			value should be the computed stress tensor, corresponding to that specific structure.  
			
		"""
		if np.int(3*(self.msteps + self.nsteps)) != len(stress_tensor_dict):
			raise ValueError("Number of stress tensors should match number of strain tensors.")
		
		for i in range(0, len(stress_tensor_dict)):
			self.defstructures[i].append(stress_tensor_dict[i])

		return self.defstructures

	def fit_cij(self, tol1=1.0, sym=False, origin=False):
		"""
		 Use the strain tensors and computed stress tensors to fit the elastic constants.
		 By default, the elastic tensor is not symmetrized. For every elastic constant, we
		 choose the largest set of strains for fitting the stress-strain relationship that
		 obeys the "linearity criterion".

		Args:
			tol1 - tolerance used for comparing linear parts of linear and quadratic fits of elastic constants
			these may not differ by more than [tol1] GPa
			sym - specifies whether or not the elastic tensor is symmetrized
			origin - specifies whether or not the linear least squares-fit should 
			be forced to pass through the origin (zero stress-zero strain point)
		"""
		self.tol1 = tol1

		Cij = np.zeros((6,6))	# elastic tensor in Voigt notation
		count = 0

		for n1 in range(0,3):	# loop over normal modes
		
			eps_tot = []
			sig_tot = []

			for n2 in range(0, np.int(self.msteps/2)):	# loop over def. magnitudes

				start_index = np.int(self.msteps/2) - n2 - 1
				end_index = np.int(self.msteps/2) + n2 + 1

				eps_array = []
				sig_array = []
				
				for n3 in range(start_index, end_index):

					n3 = np.int(n3 + n1*self.msteps)
					eps_array.append(Q.defstructures[n3][2][Q.defstructures[n3][3]])                           
					
					sig_array.append([Q.defstructures[n3][4][0,0],Q.defstructures[n3][4][1,1],Q.defstructures[n3][4][2,2],  
					Q.defstructures[n3][4][0,1], Q.defstructures[n3][4][0,2], Q.defstructures[n3][4][1,2]]) 
				
				eps_tot.append(eps_array)
				sig_tot.append(sig_array)
			
			Cij_col = fit_elas.fit_elas(eps_tot, sig_tot, tol1, origin)
			
			for Cijn in range(0,6):

				Cij[Cijn, n1] = Cij_col[Cijn, 0]					# fill up elastic tensor

		for n1 in range(3,6):   # loop over shear modes

			eps_tot = []
			sig_tot = []

			for n2 in range(0, np.int(self.nsteps/2)):  # loop over def. magnitudes

				start_index = np.int(self.nsteps/2) - n2 - 1
				end_index = np.int(self.nsteps/2) + n2 + 1

				eps_array = []
				sig_array = []

				for n3 in range(start_index, end_index):

					n3 = np.int(n3 + n1*self.nsteps)
					eps_array.append(Q.defstructures[n3][2][Q.defstructures[n3][3]])

					sig_array.append([Q.defstructures[n3][4][0,0],Q.defstructures[n3][4][1,1],Q.defstructures[n3][4][2,2],
					Q.defstructures[n3][4][0,1], Q.defstructures[n3][4][0,2], Q.defstructures[n3][4][1,2]])

				eps_tot.append(eps_array)
				sig_tot.append(sig_array)

				Cij_col = fit_elas.fit_elas(eps_tot, sig_tot, tol1, origin)

				for Cijn in range(0,6):

					Cij[Cijn, n1] = 0.5*Cij_col[Cijn, 0]			# factor 0.5 is required for shear modes only

		if sym == True:
			Cij = 0.50*(Cij + np.transpose(Cij))

		self.Cij = Cij

	def check_tensor_symmetry(self):
		"""
		This method returns the theoretical symmetry of the elastic tensor, based on the space group. 
		For this, a symmetry analyzer needs to be called, such as phonopy or cctbx.
		
		Args:
			None
		"""

#		To be implemented		


	def born_huang(self):
		"""
		This method checks the Born-Huang criterions for the mechanical stability of the crystal.
		
		Args:
			None
		"""

#		To be implemented



##### Below shows example how to use this code #####

#Q = DeformGeometry('/home/MDEJONG1/pythonplayground/pymatgen/classes/7048.cif')
Q = DeformGeometry('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen_repo/pymatgen/phonons/aluminum.cif')
Q.CIF2struct()
Q.deform(0.015, 0.015, 6, 6)
Q.get_residual_stress(np.zeros((3,3)), 1)
stress_dict = dict()

## run VASP-calculation here ##

## continue once calculations have finished ##

for i in range(0, 36):

	A = Vasprun('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen_repo/pymatgen/phonons/test4/F'+str(i)+'/vasprun.xml')
	stress_dict[i] = A.ionic_steps[-1]['stress']

Q.append_stress_tensors(stress_dict)

Q.fit_cij()
print Q.Cij

#for c in range(0, len(Q.defstructures)):

#   w = Poscar(Q.defstructures[c][0])
#   w.write_file('POSCAR_' + str(c))

#print Q.__dict__.keys()

