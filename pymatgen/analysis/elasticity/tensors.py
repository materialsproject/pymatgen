# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals
from __future__ import absolute_import

"""
This module provides a base class, for tensor-like objects, methods for
basic tensor manipulation.  It also provides a class, SquareTensor,
that provides basic methods for creating and manipulating rank 2 tensors
"""


__author__ = "Maarten de Jong, Joseph Montoya"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Wei Chen, Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ = "March 22, 2012"


from scipy.linalg import polar
import numpy as np
import itertools

class TensorBase(np.ndarray):
    """
    Base class for doing useful general operations on Nth order tensors,
    without restrictions on the type (stress, elastic, strain, piezo, etc.)
    """

    def __new__(cls, input_array):
        """
        Create a TensorBase object.  Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.

        Args:
            input_array: (3xN array-like): the 3xN array-like representing
                a tensor quantity
        """
        obj = np.asarray(input_array).view(cls)
        obj.rank = len(obj.shape)
        if not all([i == 3 for i in obj.shape]):
            raise ValueError("Pymatgen only supports 3-dimensional tensors")
        
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.rank = getattr(obj, 'rank', None)

    def __array_wrap__(self, obj):
        """
        Overrides __array_wrap__ methods in ndarray superclass to avoid errors
        associated with functions that return scalar values
        """

        if len(obj.shape) == 0:
            return obj[()]
        else:
            return np.ndarray.__array_wrap__(self, obj)
        
    def __hash__(self):
        """
        define a hash function, since numpy arrays
        have their own __eq__ method
        """
        return hash(self.tostring())

    def __repr__(cls):
        return "{}({})".format(cls.__class__.__name__,
                               cls.__str__())

    def zeroed(self, tol = 1e-3):
        """
        returns the matrix with all entries below a certain threshold
        (i.e. tol) set to zero
        """
        new_tensor = self.copy()
        new_tensor[abs(new_tensor) < tol] = 0
        return new_tensor
    
    @property
    def symmetrized(self):
        """
        Returns a symmetrized matrix from the input matrix,
        calculated by taking the sum of the matrix and its
        transpose
        """
        perms = list(itertools.permutations(range(self.rank)))
        return sum([np.transpose(self, ind) for ind in perms]) / len(perms)

    @property
    def is_symmetric(self, tol=1e-5):
        """
        Tests whether a tensor is symmetric or not based on the residual
        with its symmetric part, from self.symmetrized
        """
        return (self - self.symmetrized < tol).all()


class SquareTensor(TensorBase):
    """
    Base class for doing useful general operations on second rank tensors
    without restrictions on what type (stress, elastic, strain etc.).
    """

    def __new__(cls, input_array):
        """
        Create a SquareTensor object.  Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.  Error is thrown when the class is
        initialized with non-square matrix.

        Args:
            stress_matrix (3x3 array-like): the 3x3 array-like
                representing the Green-Lagrange strain
        """

        obj = TensorBase(input_array).view(cls)
        if not (len(obj.shape) == 2):
            raise ValueError("SquareTensor only takes 2-D "
                             "tensors as input")
        return obj
        
    @property
    def trans(self):
        """
        shorthand for transpose on SquareTensor
        """
        return SquareTensor(np.transpose(self))

    @property
    def inv(self):
        """
        shorthand for matrix inverse on SquareTensor
        """
        if self.det == 0:
            raise ValueError("SquareTensor is non-invertible")
        return SquareTensor(np.linalg.inv(self))

    @property
    def det(self):
        """
        shorthand for the determinant of the SquareTensor
        """
        return np.linalg.det(self)

    def is_symmetric(self, tol=1e-5):
        """
        Test to see if tensor is symmetric to a user-defined tolerance.
        This is determined by subtracting the transpose; if any of the
        resultant elements are above the specified tolerance, returns
        False.  Otherwise returns true.

        Args:
            tol (float): tolerance to symmetry test
        """
        return (np.abs(self - self.trans) < tol).all()

    def is_rotation(self, tol=1e-5):
        """
        Test to see if tensor is a valid rotation matrix, performs a
        test to check whether the inverse is equal to the transpose
        and if the determinant is equal to one within the specified
        tolerance

        Args:
            tol (float): tolerance to both tests of whether the
                the determinant is one and the inverse is equal
                to the transpose
        """

        return (np.abs(self.inv - self.trans) < tol).all() \
            and (np.linalg.det(self) - 1. < tol)
 
    def rotate(self, rotation):
        """
        Returns a rotated tensor based on input of a another
        rotation tensor.

        Args:
            rotation (3x3 array-like): rotation tensor, is tested
                for rotation properties and then operates on self
        """
        if self.shape != (3, 3):
            raise NotImplementedError("Rotations are only implemented for "
                                      "3x3 tensors.")
        rotation = SquareTensor(rotation)
        if not rotation.is_rotation():
            raise ValueError("Specified rotation matrix is invalid")
        return np.dot(rotation, np.dot(self, rotation.trans))

    def get_scaled(self, scale_factor):
        """
        Scales the tensor by a certain multiplicative scale factor

        Args:
            scale_factor (float): scalar multiplier to be applied to the
                SquareTensor object
        """
        return SquareTensor(self * scale_factor)

    @property
    def principal_invariants(self):
        """
        Returns a list of principal invariants for the tensor,
        which are the values of the coefficients of the characteristic
        polynomial for the matrix
        """
        if self.shape == (3, 3):
            return np.poly(self)[1:]*np.array([-1, 1, -1])
        else:
            raise ValueError("Principal invariants is only intended for use "
                             "with 3x3 SquareTensors")

    def polar_decomposition(self, side='right'):
        """
        calculates matrices for polar decomposition
        """
        return polar(self, side=side)

    
