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
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ ="Jan 24, 2012"

class Strain(object):
    """
    Class that performs various operations on deformation and strain matrices
    """

    def __init__(self, deformation_matrix):
        self.dfm = deformation_matrix
        self.strain = np.zeros((3,3))

    # construct def. matrix from indices and amount
    @staticmethod
    def from_ind_amt_dfm(matrixpos, amt):
        F = np.identity(3)
        F[matrixpos] = F[matrixpos] + amt
        return Strain(F)

    # change def. matrix by scaling it
    def change_dfm(self, ScaleFactor):
        self.dfm = self.dfm*ScaleFactor
        self.get_strain_matrix
        return Strain(self.dfm*ScaleFactor)

    # return Green-Lagrange strain matrix
    @ property
    def get_strain_matrix(self):
        self.strain = 0.5*(np.matrix(self.dfm)*np.transpose(np.matrix(self.dfm)) - np.eye(3))
        return self.strain		

mat = np.eye(3)
mat[0,0] = 1.01
ind = (0,0)

obj = Strain(mat)
print obj.__dict__.keys()
print obj.dfm
print obj.strain

#obj.change_dfm(2)
#print obj.dfm
#print obj.strain



