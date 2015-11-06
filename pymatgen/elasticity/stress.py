from __future__ import absolute_import
from pymatgen.elasticity import voigt_map
from pymatgen.elasticity.tensors import SQTensor
import math
import numpy as np
import warnings

__author__ = "Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ = "March 22, 2012"


class Stress(SQTensor):
    """
    This class extends SQTensor as a representation of the 
    stress
    """
    def __new__(cls, stress_matrix):
        """
        Create a Stress object.  Note that the constructor uses __new__ 
        rather than __init__ according to the standard method of 
        subclassing numpy ndarrays.  

        Args:
            stress_matrix (3x3 array-like): the 3x3 array-like
                representing the stress
        """
        obj = SQTensor(stress_matrix).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    def __repr__(self):
        return "Stress({})".format(self.__str__())

    @property
    def dev_principal_invariants(self):
        # TODO: check sign convention
        """
        returns the principal invariants of the deviatoric stress tensor, 
        which is calculated by finding the coefficients of the characteristic
        polynomial of the stress tensor minus the identity times the mean
        stress
        """
        return self.deviator_stress.principal_invariants

    # TODO: fix this method, is there a physical meaning to 
    #           negative J1, and how should it be handled?
    @property
    def von_mises(self):
        """
        returns the von mises stress
        """
        if not self.is_symmetric():
            raise ValueError("The stress tensor is not symmetric, Von Mises "
                             "stress is based on a symmetric stress tensor.")
        return math.sqrt(3*self.dev_principal_invariants[1])

    @property
    def mean_stress(self):
        """
        returns the mean stress
        """
        return 1./3.*self.trace()

    @property
    def deviator_stress(self):
        """
        returns the deviatoric component of the stress
        """
        if not self.is_symmetric:
            raise warnings.warn("The stress tensor is not symmetric, "
                                "so deviator stress will not be either")
        return self - self.mean_stress*np.eye(3)

    def piola_kirchoff_1(self, f):
        """
        calculates the first Piola-Kirchoff stress

        Args:
            f (3x3 array-like): rate of deformation tensor
        """
        if not self.is_symmetric:
            raise ValueError("The stress tensor is not symmetric, \
                             PK stress is based on a symmetric stress tensor.")
        f = SQTensor(f)
        return f.det*np.dot(self, f.I.T)

    def piola_kirchoff_2(self, f):
        """
        calculates the second Piola-Kirchoff stress

        Args:
            f (3x3 array-like): rate of deformation tensor
        """

        f = SQTensor(f)
        if not self.is_symmetric:
            raise ValueError("The stress tensor is not symmetric, \
                             PK stress is based on a symmetric stress tensor.")
        return f.det*np.dot(np.dot(f.I, self),f.I.T)

    @property
    def voigt(self):
        """
        returns the vector representing to the stress tensor in voigt notation
        """
        if not self.is_symmetric(1e-2):
            raise ValueError("Conversion to voigt notation requires a "
                             "symmetric stress.")
        return [self[ind] for ind in voigt_map]
