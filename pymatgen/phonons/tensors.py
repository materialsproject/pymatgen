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
from scipy.linalg import polar
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
    Class for doing useful general operations on *square* rank 2 tensors, without 
    restrictions on what type (stress, elastic, strain etc.).
    Error is thrown when the class is initialized with non-square matrix.
    """

    def __new__(cls, input_matrix):
        obj = np.asmatrix(input_matrix).view(cls)
        if obj.shape[0] != obj.shape[1]:
            raise ValueError("SQTensor only takes square arrays as input")
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    def __repr__(self):
        return "SQTensor({})".format(self.__str__())
    
    @property
    def T(self):
        """
        shorthand for transpose on SQTensor
        """
        return SQTensor(np.transpose(self))

    @property
    def I(self):
        """
        shorthand for matrix inverse on SQTensor
        """
        return SQTensor(np.linalg.inv(self))

    @property
    def det(self):
        """
        shorthand for the determinant of the SQTensor
        """
        return np.linalg.det(self)

    # TODO: Note that the boolean properties have keywords,
    #           but those can't be specified, since they're 
    #           decorated as properties.  We could either add
    #           a tolerance property to use in the boolean properties
    #           or perhaps redecorate them as methods -JHM

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
        return (np.abs(self - self.T) < tol).all()

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

    # TODO: JHM asks should this be implemented?
    """def rotate(self, rotation):
        super(StressOps, self).__init__(rotation)
        super(StressOps, self).is_rotatio()

        if super(StressOps, self).is_rotatio() == False:
            raise ValueError("Not a valid rotation matrix")

        else:
            return rotation*self._StressMatrix*np.transpose(rotation)"""

    def rotate(self, rotation):
        if not self.is_rotation():
            raise ValueError("Specified rotation matrix is invalid")
        raise NotImplementedError("matrix rotations are not yet supported")
    
    def get_scaled(self, scale_factor):
        """
        Scales the tensor by a certain multiplicative scale factor
        """
        return SQTensor(self * scale_factor)

    @property
    def principal_invariants(self):
        """
        returns a list of principal invariants for the tensor,
        which are the values of the coefficients of the characteristic 
        polynomial for the matrix
        """
        # TODO: JM asks whether this fulfills the necessary sign conventions
        return np.poly(self)[1:] 


voigt_map = [[0,0],[1,1],[2,2],[1,2],[2,0],[0,1]]

class ElasticTensor(SQTensor):
    """
    This class extends SQTensor to describe the 6x6
    elastic tensor, C_{ij}.
    """

    def __new__(cls, input_matrix):
        obj = SQTensor(input_matrix).view(cls)
        if obj.shape[0]!=6:
            raise ValueError("Elastic tensor input must be\
                             a 6x6 array")
        if not obj.is_symmetric:
            raise ValueError("Elastic tensor input must be\
                             symmetric")
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    def __repr__(self):
        return "ElasticTensor({})".format(self.__str__())

    # TODO: JHM suggests might want to split this up
    @property
    def KG_average(self):
        """
        calculates the Voigt-Reuss-Hill average
        """
        K_voigt = self[:3,:3].mean() 
        G_voigt = (2.*self[:3,:3].trace() - np.triu(self[:3,:3]).sum() \
                    + 3*self[3:,3:].trace()) / 15.
        K_reuss = 1./self.I[:3,:3].sum()
        G_reuss = 15./(8.*self.I[:3,:3].trace() \
                    - 4.*np.triu(self.I[:3,:3]).sum() \
                    + 3.*self.I[3:,3:].trace())
        K_vrh = 0.5*(K_voigt + K_reuss)
        G_vrh = 0.5*(G_voigt + G_reuss)
        
        return [K_voigt, G_voigt, K_reuss, G_reuss, K_vrh, G_vrh]

    @property
    def universal_anisotropy(self):
        """
        calculates value for universal anisotropy
        """
        average_Cij = self.KG_average
        return 5.*average_Cij[1]/average_Cij[3] \
                + average_Cij[0]/average_Cij[2] - 6.
    
    @property
    def polar_decomposition(self,side='right'):
        """
        calculates matrices for polar decomposition
        """
        return polar(self,side=side)

    @property
    def homogeneous_poisson(self):
        """
        calculates homogeneous poisson
        """
        average_Cij = self.KG_average
        return (1. - 2./3. * average_Cij[5]/average_Cij[4]) / \
                (2. + 2./3. * average_Cij[5]/average_Cij[4])

    
    @property
    def full_tensor(self):
        """
        Returns the tensor in standard notation (i. e. 
        a 4th order 3-dimensional tensor, C_{ijkl}), which
        is represented in a np.array with shape (3,3,3,3)
        """
        c = np.zeros((3,3,3,3))
        for p in range(6):
            for q in range(6):
                i,j = voigt_map[p]
                k,l = voigt_map[q]
                c[i,j,k,l] = c[j,i,k,l] = c[i,j,l,k] = \
                        c[j,i,l,k] = c[k,l,i,j] = self[p,q]
        return c

    @classmethod
    def from_full_tensor(cls,c_ijkl):
        """
        Factory method to construct elastic tensor from fourth order
        tensor C_ijkl.  First tests for appropriate symmetries and then 
        constructs the 6x6 voigt notation tensor.

        Args:
            c_ijkl (3x3x3x3 array-like): input fourth order tensor corresponding
                to the full elastic tensor
        """
        # Test symmetry of elastic tensor
        c_ijkl = np.array(c_ijkl)
        tol = 1e-10
        if not ((c_ijkl - np.transpose(c_ijkl,(1,0,2,3)) < tol).all() and
                (c_ijkl - np.transpose(c_ijkl,(0,1,3,2)) < tol).all() and
                (c_ijkl - np.transpose(c_ijkl,(1,0,3,2)) < tol).all() and
                (c_ijkl - np.transpose(c_ijkl,(3,2,0,1)) < tol).all()):
            raise ValueError("Input elasticity tensor does \
                             not satisfy necessary symmetries")
        # Construct elastic tensor
        c_pq = np.zeros((6,6))
        for p in range(6):
            for q in range(6):
                i,j = voigt_map[p]
                k,l = voigt_map[q]
                c_pq[p,q] = c_ijkl[i,j,k,l]

        return cls(c_pq)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

    eye = np.identity(3)
    sigma = SQTensor(np.random.randn(3,3))

    mat1 = np.random.randn(6,6)
    mat1 = mat1 + np.transpose(mat1)

    sigma2 = ElasticTensor(mat1)
    sigma2.KG_average
    sigma2.universal_anisotropy
    sigma2.polar_decomposition
    this_tensor = sigma2.full_tensor
    #print sigma2
    sigma3 = ElasticTensor.from_full_tensor(this_tensor)
    this_tensor[0,1,0,1] += 0.1
    sigma4 = ElasticTensor.from_full_tensor(this_tensor)

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


