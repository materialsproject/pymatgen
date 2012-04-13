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

class Cij(object):
    """
    Class for doing useful general operation on *square* matrices, without 
    restrictions on what type of matrix (stress, elastic, strain etc.).
    An error is thrown when the class is initialized with non-square matrix.
    """

    def __init__(self, matrix):
        self._matrix = matrix
        if self._matrix.shape[0]!=self._matrix.shape[1]:
            raise ValueError("This class operates only on square matrices")
        
    def is_symmetric(self, tol=0.001):
        if len(np.nonzero(np.abs(self._matrix-np.transpose(self._matrix))>tol)[0]) == 0:
            return True
        else:
            return False

    def is_rotation_matrix(self, tol=0.001):
        
        arr1 =  np.nonzero(np.abs(np.linalg.inv(self._matrix)-np.transpose(self._matrix))>tol)[0]
        arr2 = np.linalg.det(self._matrix)
        arr2 = np.float(arr2 - 1)
        print arr1 ==0 and arr2 < tol


    def symmetrize(self):
        return 0.5*(self._matrix + np.transpose(self._matrix))

    def value(self):
        return self._matrix

class StressOps(Cij):
    """
    This class contains methods that are to operate on stress tensors specifically, hence
    matrices passed into this class need to be 3x3, otherise an error is thrown. Moreover,
    the stress tensor is required to be symmetric up to a tolerance tol
    """
    def __init__(self, stress_matrix):

        super(StressOps, self).__init__(stress_matrix)
        if super(StressOps, self).is_symmetric() == False:
            raise ValueError("The stress tensor needs to be symmetric up to a tolerance tol")

    def transform(self, rotation):
        
        super(StressOps, self).__init__(rotation)
        super(StressOps, self).is_rotation_matrix()

#        if super(StressOps, self).is_rotation_matrix() == False:
#            raise ValueError("Not a valid rotation matrix")
#        else:
#            print 'passed test'


class ElasticOps(Cij):
    """
    This class contains methods that are to operate on elastic tensors specifically, hence
    matrices passed into this class need to be 6x6, otherise error is thrown. 
    """
    
        

if __name__ == "__main__":

    sigma = Cij(np.random.randn(3,3))
    sigma = sigma.symmetrize()
    rot = np.identity(3)

    B = StressOps(sigma)
    B.transform(rot)

