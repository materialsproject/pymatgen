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

    @property
    def symmetrize(self):
        return 0.5*(self._matrix + np.transpose(self._matrix))


    def rotate(self, rotation):
        if not self.is_rotation_matrix():
            raise ValueError("Specified rotation matrix is invalid")
        raise NotImplementedError("matrix rotations are not yet supported")

    @property
    def tensor(self):
        return self._matrix

    @property
    def principal_invariants(self):
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

    @property
    def KG_average(self):
        if np.shape(self._matrix)[0] != 6 or np.shape(self._matrix)[1] != 6:
            raise ValueError("This method only makes sense for 6x6 elastic tensors")
        if self.is_symmetric()==False:
            raise ValueError("This method takes only symmetric tensors")

        K_Voigt = (self._matrix[0,0]+self._matrix[1,1]+self._matrix[2,2])*1.0/9 + (self._matrix[0,1]+self._matrix[0,2]+self._matrix[1,2])*2.0/9
        G_Voigt = (self._matrix[0,0]+self._matrix[1,1]+self._matrix[2,2]-self._matrix[0,1]-self._matrix[0,2]-self._matrix[1,2])*1.0/15 + (self._matrix[3,3]+self._matrix[4,4] + self._matrix[5,5])*1.0/5

        S = np.linalg.inv(self._matrix)

        K_Reuss = 1.0/(S[0,0]+S[1,1]+S[2,2]+2*S[0,1]+2*S[0,2]+2*S[1,2])
        G_Reuss = 15.0/(4*(S[0,0]+S[1,1]+S[2,2]-S[0,1]-S[0,2]-S[1,2])+3*(S[3,3]+S[4,4]+S[5,5]))

        K_Voigt_Reuss_Hill = 0.5*(K_Voigt+K_Reuss)
        G_Voigt_Reuss_Hill = 0.5*(G_Voigt+G_Reuss)

        average_Cij = [K_Voigt, G_Voigt, K_Reuss, G_Reuss, K_Voigt_Reuss_Hill, G_Voigt_Reuss_Hill]

        return average_Cij

    @property
    def universal_anisotropy(self):
        if self.is_symmetric()==False:
            raise ValueError("This method takes only symmetric tensors")
        average_Cij = self.KG_average
        ua = 5*average_Cij[1]/average_Cij[3]+average_Cij[0]/average_Cij[2]-6
        return ua

    @property
    def homogeneous_poisson(self):
        if self.is_symmetric()==False:
            raise ValueError("This method takes only symmetric tensors")
        average_Cij = self.KG_average
        nu = (1-2.0/3*average_Cij[5]/average_Cij[4])/(2+2.0/3*average_Cij[5]/average_Cij[4])






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
#    print sigma.PrincipalInvariants

    mat1 = np.random.randn(6,6)
    mat1 = mat1 + np.transpose(mat1)

    sigma2 = SQTensor(mat1)

    sigma2.universal_anisotropy


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


