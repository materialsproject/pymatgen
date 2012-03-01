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

class matrix_ops(object):
	
	def __init__(self, matrix=np.zeros((3,3))):
		self._matrix = matrix
#		self.trace = 3

	def apply_homo_def(self, scaling):
		self._matrix = scaling*self._matrix
		return matrix_ops(self._matrix)
	
#	def trace(self):
#		trace = self._matrix[0,0] + self._matrix[1,1] + self._matrix[2,2]
#		self.trace = trace
#		return self.trace
#		self.trace = trace
#		return(self.trace)
#		return trace

	@property
	def plus_mat(self):
		return self._matrix + 2



s = matrix_ops(np.identity(3))
print s.__dict__.keys()
#s.plus_mat



#print s._matrix
#print s.trace()
#print s.trace()

#s.apply_homo_def(2)
#print s._matrix
#print s.trace()



#print s.__dict__.keys()
#print s.matrix
#print s.trace

#s.apply_homo_def(2)
#print s.matrix
#print s.trace
