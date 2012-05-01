import warnings, sys, os, math
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

class SQTensor(object):
    """
    Class for doing useful general operation on *square* matrices, without 
    restrictions on what type of matrix (stress, elastic, strain etc.).
    An error is thrown when the class is initialized with non-square matrix.
    TODO: change class names
    TODO: AJ asks if there is an existing matrix implementation to subclass instead of object, that for instance already implements 'rotate'
    """
    def __init__(self, matrix):
        self._matrix = matrix
        if self._matrix.shape[0]!=self._matrix.shape[1]:
            raise ValueError("SQTensor operates only on square matrices")
        
    def is_symmetric(self, tol=0.001):
        if len(np.nonzero(np.abs(self._matrix-np.transpose(self._matrix))>tol)[0]) == 0:
            return True
        else:
            return False

    def is_rotation_matrix(self, tol=0.001):
        arr1 = np.nonzero(np.abs(np.linalg.inv(self._matrix)-np.transpose(self._matrix))>tol)[0]
        arr2 = np.linalg.det(self._matrix)
        arr2 = np.float(arr2 - 1)

        if len(arr1) ==0 and arr2 < tol:
            return True
        else:
            return False

    def get_scaled(self, scale_factor):
        return self._matrix*scale_factor
        
    def symmetrize(self):
        return 0.5*(self._matrix + np.transpose(self._matrix))


    def rotate(self, rotation):
        if not self.is_rotation_matrix():
            raise ValueError("Specified rotation matrix is invalid")
        raise NotImplementedError("matrix rotations are not yet supported")

    def returntensor(self):
        return self._matrix

    def value(self, i, j):         
        return self._matrix[i, j]

    @property
    def PrincipalInvariants(self):
        if np.shape(self._matrix)[0] != 3:
            raise NotImplementedError("Principal Invariants are currently only supported for 3x3 matrices.")
        I1 = self._matrix[0,0] + self._matrix[1,1] + self._matrix[2,2]
        I2 = self._matrix[0,0]*self._matrix[1,1]+ \
             self._matrix[1,1]*self._matrix[2,2]+ \
             self._matrix[0,0]*self._matrix[2,2]- \
             self._matrix[0,1]**2 - self._matrix[1,2]**2 - self._matrix[2,0]**2
        I3 = np.linalg.det(self._matrix)
        I = [I1, I2, I3]
        return I

#    def rotate(self, rotation):
#        super(StressOps, self).__init__(rotation)
#        super(StressOps, self).is_rotation_matrix()

#        if super(StressOps, self).is_rotation_matrix() == False:
#            raise ValueError("Not a valid rotation matrix")

#        else:
#            return rotation*self._StressMatrix*np.transpose(rotation)

if __name__ == "__main__":

    eye = np.identity(3)
    sigma = SQTensor(np.random.randn(3,3))
    sigma.PrincipalInvariants

    

#    print sigma.is_rotation_matrix(np.matrix(eye))
#    print np.linalg.inv(eye)
        
#    sigma = SQTensor(np.identity(3))    
#    print sigma.is_symmetric()
#    eye = np.identity(3)   

#    print sigma.symmetrize()        
#    sigma.rotate(np.identity(3))

#    print sigma.rotate(eye)

#    sigma = sigma.symmetrize()
#    rot = np.identity(3)
#    B = StressOps(sigma)
#    print B.MeanStress
#    print B.DeviatorStress
#    print B.DeviatorPrincipalInvariants    
#    print B.VonMises


