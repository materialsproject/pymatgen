# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals
from __future__ import absolute_import

"""
This module provides classes for the Piezoelectric tensor
"""
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
import warnings
from six.moves import range

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Feb, 2016"

voigt_map = [(0, 0), (1, 1), (2, 2), (1, 2), (2, 0), (0, 1)]


class PiezoTensor(np.ndarray):
    """
    This class describes the 3 x 6 piezo tensor in Voigt-notation
    """

    def __new__(cls, input_array):
        """
        Create an PiezoTensor object.  The constructor throws an error if
        the shape of the input_matrix argument is not 3x6, i. e. in Voigt-
        notation. Note that the constructor uses __new__ rather than
        __init__ according to the standard method of subclassing numpy
        ndarrays.

        Args:
            input_matrix (3x6 array-like): the Voigt-notation 3x6 array-like
                representing the piezo tensor
        """
        obj = np.asarray(input_array).view(cls)
        if obj.shape != (3, 6):
            raise ValueError("Default piezo tensor constructor requires "
                             "input argument to be the Voigt-notation 3x6 "
                             "array.  To construct from a 3x3x3 array, use "
                             "PiezoTensor.from_full_tensor")
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    @property
    def full_tensor(self):
        """
        Returns the tensor in standard notation (i. e.
        a 3th order 3-dimensional tensor, C_{ijk}), which
        is represented in a np.array with shape (3,3,3)
        """
        c = np.zeros((3, 3, 3))
        for i in range(3):
            for p in range(6):
                j, k = voigt_map[p]
                c[i, j, k] = c[i, k, j] = self[i, p]
        return c

    @classmethod
    def from_full_tensor(cls, c_ijk, tol=1e-5):
        """
        Factory method to construct piezo tensor from third order
        tensor C_ijk.  First tests for appropriate symmetries and then
        constructs the 3x6 voigt notation tensor.

        Args:
            c_ijkl (3x3x3 array-like): third order tensor corresponding
                to the piezo tensor
            tol (float): tolerance for the symmetry test of the tensor
        """
        # Test symmetry of elastic tensor
        c_ijk = np.array(c_ijk)
        if not (c_ijk - np.transpose(c_ijk, (0, 2, 1)) < tol).all():
            raise ValueError("Input piezo tensor does "
                             "not satisfy necessary symmetries")
        # Construct elastic tensor
        c_ip = np.zeros((3, 6))
        for i in range(3):
            for p in range(6):
                j, k = voigt_map[p]
                c_ip[i, p] = c_ijk[i, j, k]

        return cls(c_ip)

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

    def transform(self, sym):
        """
        Returns a transformed tensor based on input of symmetry operation

        Args:
            symp (SymmOp): symmetry operation
        """
        return PiezoTensor.from_full_tensor(sym.transform_tensor(self.full_tensor))

    def is_valid(self, structure, symprec=0.1, tol=1e-3):
        """
        Checks the piezo tensor against the symmetries of the structure and
        determine if the piezoelectric tensor is valid for that structure

        Args:
            structure (Structure): structure to check against
            symprec (float): Tolerance for symmetry finding, default good for
                             as-made structures. Set to 0.1 for MP structures
            tol (float): tolerance for equality
        """

        sg = SpacegroupAnalyzer(structure, symprec)

        # Check every symmetry operation in the space group
        for symm in sg.get_symmetry_operations(cartesian=True):
            # does it produce the same tensor?
            diff = self.transform(symm) - self
            if not (np.abs(diff) < tol).all():
                print(diff)
                return False

        return True

    def symmeterize(self, structure, symprec=0.1, tol=1e-3):
        """
        Averages the piezo tensor based on the symmetries of the structure

        Args:
            structure (Structure): structure to check against
            symprec (float): Tolerance for symmetry finding, default good for
                             as-made structures. Set to 0.1 for MP structures
            tol (float): tolerance for zero'ing out indicies. The average procedure produces very small numbers rather then 0. This tolerance is used to zero out those values to make the tensor less messy.
        """

        sg = SpacegroupAnalyzer(structure, symprec)
        new_pt = PiezoTensor(self)

        for symm in sg.get_symmetry_operations(cartesian=True):
            new_pt = (new_pt + new_pt.transform(symm)) / 2

        low_values_indices = np.abs(new_pt) < tol
        new_pt [low_values_indices] = 0

        return new_pt
