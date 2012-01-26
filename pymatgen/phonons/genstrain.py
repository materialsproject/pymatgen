from __future__ import division
import unittest
import pymatgen
import sys
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

__author__="Maarten de Jong"
__copyright__ = "Copyright 2011, The Materials Project"
__credits__ = "Mark Asta"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ ="Jan 24, 2012"

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
			/home/MDEJONG1/pythonplayground/pymatgen/classes/7048.cif
		"""
		
		self._path = path_to_cif
	
	def CIF2struct(self):
		"""
		Parse CIF to pymatgen core structure for easy manipulation
		
		Args: 
			none
		"""

		coords = list()
		species = []
		myCIF = CifParser(self._path).get_structures()[0]
#		print myCIF
		s = StructureEditor(myCIF)
		


#		S2 = StructureEditor(myCIF).replace_species({'Al', 'O'})
#		print S2	


#		print type(myCIF)
#		print myCIF.__dict__.keys()
#		print myCIF._lattice.__dict__.keys()
#		print myCIF._lattice._mc2d
#		print myCIF._lattice._md2c
#		print myCIF._lattice._matrix

		

#		check this part, cifio-code does this already!

#		print myCIF


		for k in range(len(myCIF._sites)):
			coords.append(Lattice(np.linalg.inv(myCIF._lattice._mc2d)).get_fractional_coords(myCIF._sites[k]._coords))
			species.append(myCIF._sites[k]._species.keys()[0])

		lattice = Lattice(myCIF._lattice._matrix)
		struct = Structure(lattice,species,coords)

		

#		print struct

		self.base_struct = struct
#		self._coords = coords
#		self._lattice = myCIF._lattice._matrix
#		self._species = species

		return self.base_struct

	def get_residual_stress(self):
		"""
		Given a crystal structure, this method starts up a static calculation of the stress tensor. These
		are residual stresses which need to be subtracted later on. 
		
		Args:
			none
		"""
		
		self._residual_stress = np.zeros((3,3))


		return self._residual_stress

	def deform(self, nd=0.02, ns=0.02, m=6, n=6): # this class still requires some user input
		"""
		Take original geometry (from CIF-file) and apply a range of deformations. 
		Default values generally work well for metals. However, one might need to
		make changes for material such as oxides.

		Args: 
			nd - maximum amount of normal strain  
			ns - maximum amount of shear strain
			m - number of normal deformations used for fitting elastic constants per deformation mode, even integer required
			n - number of shear deformations used for fitting elastic constants per deformation mode, even integer required
		"""
		
		self.msteps = m
		self.nsteps = n

		if m%2!=0:
			raise ValueError("m has to be even.")

		if n%2!=0:
			raise ValueError("n has to be even.")

		mystrains = np.zeros((3, 3, np.int(m*3) + np.int(n*3)))

		defs = np.linspace(-nd, nd, num=m+1)
		defs = np.delete(defs, np.int(m/2), 0)
		defstructures = []

		counter = 0
		# First apply non-shear deformations
		# Depending on crystal symmetry, these may yield shear stresses, however
		for i1 in range(0, 3):

			for i2 in range(0, len(defs)):
				
				F = np.eye(3)
				F[i1, i1] = F[i1, i1] + defs[i2]
				C = np.transpose(F)*F
				E = 0.5*(C-np.eye(3)) # Green-Lagrange strain tensor
				Lp = self._lattice*F  # Deformation mapping of lattice vectors
				mystrains[:, :, counter] = E

				lat = Lattice(Lp)
				struc = Structure(lat,self._species,self._coords)
				defstructures.append(struc)
				counter += 1

		# Now apply shear deformations #		
		sheardef = np.linspace(-ns, ns, num=n+1)
		sheardef = np.delete(sheardef, np.int(n/2), 0)

		F_index = [[0, 1], [0, 2], [1, 2]]
		for j1 in range(0, 3):
		
			for j2 in range(0, len(sheardef)):

				F = np.zeros((3,3))
				F = np.eye(3)
				F[F_index[j1][0], F_index[j1][1]] = F[F_index[j1][0], F_index[j1][1]] + sheardef[j2]
				F = np.matrix(F)
		
				C = np.transpose(F)*F
				E = 0.5*(C-np.eye(3)) # Green-Lagrange strain tensor
				Lp = self._lattice*F  # Deformation mapping of lattice vectors
				mystrains[:, :, counter] = E
				
				lat = Lattice(Lp)
				struc = Structure(lat,self._species,self._coords)
				defstructures.append(struc)
				counter += 1
		
		self.defstructures = defstructures
		self.mystrains = mystrains

		return self.defstructures, self.mystrains


	def get_stress_tensors(self, path_list):
		"""
		This method takes a list of paths, which contain output files (currently only vasprun.xml's are supported)
		containing the stress tensor, and stores these in a 3D-matrix.
		
		Args:
			path_list - list of paths that contain vasprun.xml's which contain calculated stress tensor for every
			deformed geometry
		"""

		self.path_list = path_list

		if np.int((self.msteps + self.nsteps)*3) != len(path_list):
			raise ValueError("Number of stress tensors does not match number of strain tensors.")

		stress_tensors = np.zeros((3, 3, np.int(self.msteps*3) + np.int(self.nsteps*3)))
		
		counter = 0
		for stress_path in path_list:
			
			A = Vasprun(stress_path)
			stressmatrix = A.ionic_steps[-1]['stress']
			stress_tensors[:, :, counter] = stressmatrix - self._residual_stress
			counter += 1


		self.stress_tensors = stress_tensors


	def fit_cij(self, tol1=0.02, sym=False):
		"""
		Take stress and strain tensors and fit elastic constants. By default, the elastic tensor is not symmetrized.

		Args:
			tol1 - tolerance used for comparing linear and quadratic fits of elastic constants
			sym - specifies whether or not the elastic tensor is symmetrized
		"""



Q = DeformGeometry('/home/MDEJONG1/pythonplayground/pymatgen/classes/7048.cif')


## Calling sequence ##
Q.CIF2struct()



#Q.deform(0.02, 0.02, 4, 4)
#Q.get_residual_stress()

# Perform VASP-calculations in specified paths #
#Q.get_stress_tensors(['/home/MDEJONG1/pythonplayground/pymatgen/classes/vasprun.xml']*24)



#print Q.__dict__.keys()
#print Q.defstructures
#print Q.mystrains

#A = Vasprun('vasprun.xml')
#print A.ionic_steps[-1]['stress']




