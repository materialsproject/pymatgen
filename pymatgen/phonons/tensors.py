import warnings, sys, os, math
import unittest
import pymatgen
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp import Vasprun
from pymatgen.io.cif import CifWriter
from pymatgen.io.cif import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import *
import numpy as np

__author__="Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ ="March 22, 2012"

class SQTensor(np.matrix):
    """
    Class for doing useful general operations on *square* matrices, without 
    restrictions on what type of matrix (stress, elastic, strain etc.).
    An error is thrown when the class is initialized with non-square matrix.

    .. attribute::

        matrix corresponding to the tensor values
    """

    def __new__(cls, matrix):
        super(SQTensor,cls).__new__(matrix)

    def __init__(self, matrix):
        """
        Creates a SQTensor 

        Args:
            matrix (NxN array-like): input array-like with square shape.  If this
            argument is not square, throws an ValueError.
        """
        super(SQTensor,self).__init__(matrix)
        if self.shape[0] != self.shape[1]:
            raise ValueError("SQTensor operates only on square matrices")

    def __array_finalize__
    def __repr__(self):
        return "SQTensor({})".format(self.__str__())
    '''
    def __str__(self):
        return self
    '''
    @property
    def det(self):
        '''
        returns the determinant of the SQTensor
        '''
        return np.linalg.det(self)

    @property
    def is_symmetric(self, tol=0.001):
        """
        Test to see if tensor is symmetric to a user-defined tolerance.
        This is determined by subtracting the transpose; if any of the 
        resultant elements are above the specified tolerance, returns 
        False.  Otherwise returns true.

        Args:
            tol (float): tolerance to test whether the matrix is symmetric
        """
        # TODO: write a unit test
        # TODO: JM asks whether tolerance test is necessay
        return (np.abs(self - self.T) < tol).all()

        """
        return self == self
        if len(np.nonzero(np.abs(self-np.transpose(self))>tol)[0]) == 0:
            return True
        else:
            return False
        """

    @property
    def is_rotation(self, tol=0.001):
        """
        Test to see if tensor is a valid rotation matrix, performs a 
        test to check whether the inverse is equal to the transpose
        and if the determinant is equal to one within the specified
        tolerance

        Args:
            tol (float): tolerance to both tests of whether the 
                the determinant is zero and the inverse is equal
                to the transpose
        """

        return (np.abs(self.I - self.T) < tol).all() \
                and (np.linalg.det(self) - 1. < tol)
        
    @property
    def symmetrized(self):
        """
        Returns a symmetrized matrix from the input matrix, 
        calculated by taking the sum of the matrix and its
        transpose
        """
        return 0.5*(self + self.T)

    def rotate(self, rotation):
        if not self.is_rotation():
            raise ValueError("Specified rotation matrix is invalid")
        raise NotImplementedError("matrix rotations are not yet supported")
    
    def get_scaled(self, scale_factor):
        '''
        Scales the tensor by a certain multiplicative scale factor
        '''
        return SQTensor(self * scale_factor)

    @property
    def principal_invariants(self):
        '''
        returns a list of principal invariants for the tensor,
        which are the values of the coefficients of the characteristic 
        polynomial for the matrix
        '''
        # TODO: JM asks whether this fulfills the necessary sign conventions
        return np.poly(self)[1:]
 
    @property
    def KG_average(self):
        '''
        calculates the Voigt-Reuss-Hill average
        '''
        if self.shape[0] != 6:
            raise ValueError("This method only makes sense for 6x6 elastic tensors")
        if self.is_symmetric == False:
            raise ValueError("This method takes only symmetric tensors")

        #k_voigt = self[:3,:3].trace()/9. + self
        K_Voigt = (self[0,0]+self[1,1]+self[2,2])*1./9. \
                + (self[0,1]+self[0,2]+self[1,2])*2./9.
        G_Voigt = (self[0,0]+self[1,1]+self[2,2] - \
                   self[0,1]-self[0,2]-self[1,2])*1./15.\
                + (self[3,3]+self[4,4] + self[5,5])*1./5.

        S = self.I

        K_Reuss = 1.0/(S[0,0]+S[1,1]+S[2,2]+2*S[0,1]+2*S[0,2]+2*S[1,2])
        G_Reuss = 15.0/(4*(S[0,0]+S[1,1]+S[2,2]-S[0,1]-S[0,2]-S[1,2])+3*(S[3,3]+S[4,4]+S[5,5]))

        K_Voigt_Reuss_Hill = 0.5*(K_Voigt+K_Reuss)
        G_Voigt_Reuss_Hill = 0.5*(G_Voigt+G_Reuss)

        average_Cij = [K_Voigt, G_Voigt, K_Reuss, G_Reuss, 
                       K_Voigt_Reuss_Hill, G_Voigt_Reuss_Hill]

        return average_Cij

    @property
    def universal_anisotropy(self):
        '''
        calculates value for universal anisotropy, only valid for
        symmetric tensors, throws an error if not symmetric
        '''
        if self.is_symmetric == False:
            raise ValueError("This method takes only symmetric tensors")
        average_Cij = self.KG_average
        ua = 5*average_Cij[1]/average_Cij[3] + average_Cij[0]/average_Cij[2] - 6
        return ua

    @property
    def homogeneous_poisson(self):
        '''
        calculates homogeneous poisson 
        '''
        # TODO: JM asks should this be implemented?
        if self.is_symmetric == False:
            raise ValueError("This method takes only symmetric tensors")
        average_Cij = self.KG_average
        nu = (1 - 2./3. * average_Cij[5]/average_Cij[4]) / \
                (2 + 2./3. * average_Cij[5]/average_Cij[4])

        # TODO: JM asks should this be implemented?
#    def rotate(self, rotation):
#        super(StressOps, self).__init__(rotation)
#        super(StressOps, self).is_rotatio()

#        if super(StressOps, self).is_rotatio() == False:
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
    print sigma2

#    print sigma.is_rotatio(np.matrix(eye))
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


