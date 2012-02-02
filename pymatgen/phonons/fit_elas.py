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

__author__="Maarten de Jong"
__copyright__ = "Copyright 2011, The Materials Project"
__credits__ = "Mark Asta"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ ="Jan 29, 2012"

def fit_elas(strain, stress, tol, origin):
	"""
	Function for fitting stress-strain. Works correctly, but some additional safety checks should be added.
	"""

	Cij = np.zeros((6,1))					# represents one column in elastc tensor

	for k in range(0, 6):					# loop over Cij in every column in the elastic tensor

		C = 0
		
		for mag in range(0, len(strain)):	# loop over strain magnitudes, used for calculating each Cij

			linear_flag = False

			if len(stress[mag]) > 2:		# do not attempt to fit parabola for only 2 strain - stress points

				sig_fit = []
			
				for n in range(0, len(stress[mag])):

					sig_fit.append(stress[mag][n][k])					# create array containing stresses

				sig_fit.reverse()
				
				if origin == False:
					
					lincoeff = np.polyfit(strain[mag], sig_fit, 1)*0.1		# linear fit, convert to GPa
				
				else:														# force line through origin
					
					lincoeff = np.sum(np.array(strain[mag])*np.array(sig_fit))/np.sum(np.array(strain[mag])*np.array(strain[mag]))*0.1
				
				quadcoeff = np.polyfit(strain[mag], sig_fit, 2)*0.1		# quadratic fit, convert to GPa

				if mag == 1 and np.abs(lincoeff[0] - quadcoeff[1]) > tol:

					raise ValueError("Linearity criterion failed on smallest amount of strain (4 points), please use smaller strains.")

				if np.abs(lincoeff[0] - quadcoeff[1]) < tol:			# check linearity condition
		
					linear_flag = True
					C = lincoeff[0]

		Cij[k, 0] = C	
	
	return Cij

