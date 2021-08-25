# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module provides the Stress class used to create, manipulate, and
calculate relevant properties of the stress tensor.
"""

import math

import numpy as np

from pymatgen.core.tensors import SquareTensor

__author__ = "Joseph Montoya"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Maarten de Jong, Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Production"
__date__ = "July 24, 2018"


class Stress(SquareTensor):
    """
    This class extends SquareTensor as a representation of the
    stress
    """

    symbol = "s"

    def __new__(cls, stress_matrix):
        """
        Create a Stress object.  Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.

        Args:
            stress_matrix (3x3 array-like): the 3x3 array-like
                representing the stress
        """
        obj = super().__new__(cls, stress_matrix)
        return obj.view(cls)

    @property
    def dev_principal_invariants(self):
        """
        returns the principal invariants of the deviatoric stress tensor,
        which is calculated by finding the coefficients of the characteristic
        polynomial of the stress tensor minus the identity times the mean
        stress
        """
        return self.deviator_stress.principal_invariants * np.array([1, -1, 1])

    @property
    def von_mises(self):
        """
        returns the von mises stress
        """
        if not self.is_symmetric():
            raise ValueError(
                "The stress tensor is not symmetric, Von Mises " "stress is based on a symmetric stress tensor."
            )
        return math.sqrt(3 * self.dev_principal_invariants[1])

    @property
    def mean_stress(self):
        """
        returns the mean stress
        """
        return 1.0 / 3.0 * self.trace()

    @property
    def deviator_stress(self):
        """
        returns the deviatoric component of the stress
        """
        if not self.is_symmetric:
            raise ValueError("The stress tensor is not symmetric, so deviator stress will not be either")
        return self - self.mean_stress * np.eye(3)

    def piola_kirchoff_1(self, def_grad):
        """
        calculates the first Piola-Kirchoff stress

        Args:
            def_grad (3x3 array-like): deformation gradient tensor
        """
        if not self.is_symmetric:
            raise ValueError(
                "The stress tensor is not symmetric, \
                             PK stress is based on a symmetric stress tensor."
            )
        def_grad = SquareTensor(def_grad)
        return def_grad.det * np.dot(self, def_grad.inv.trans)

    def piola_kirchoff_2(self, def_grad):
        """
        calculates the second Piola-Kirchoff stress

        Args:
            def_grad (3x3 array-like): rate of deformation tensor
        """

        def_grad = SquareTensor(def_grad)
        if not self.is_symmetric:
            raise ValueError(
                "The stress tensor is not symmetric, \
                             PK stress is based on a symmetric stress tensor."
            )
        return def_grad.det * np.dot(np.dot(def_grad.inv, self), def_grad.inv.trans)
