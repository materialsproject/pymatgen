import warnings, sys, os
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/') # (If one does not want to change $PYTHONPATH)
import unittest
import pymatgen
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Vasprun
from pymatgen.io.cifio import CifWriter
from pymatgen.io.cifio import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import *
from pymatgen.core.structure_modifier import StructureEditor
import numpy as np

__author__="Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ ="March 22, 2012"


class Stress(object):
 
    def __init__(self, stress_matrix):
        self._sigma = stress_matrix

    # return a scaled version of this matrix
    def get_scaled(self, scale_factor):
        stress_matrix = self._sigma * scale_factor
        return Stress(stress_matrix)

    @property
    def issymmetric(self, tol=0.001):
        s= self._sigma
        st = np.transpose(self._sigma)

        if len (np.nonzero(np.abs(s-st)>tol)[0]) == 0:
            return True
        else:
            return False

    @property
    def stress_matrix(self):
        return self._sigma

    @property
    def value(self, i, j):         # put value in matrix method
        return self._sigma[i, j]

    
if __name__ == "__main__":

    mat = np.eye(3)
    mat[0,2] = 0.1
    mat[2,0] = 0.1
    s = Stress(mat)
    print s.issymmetric
    


