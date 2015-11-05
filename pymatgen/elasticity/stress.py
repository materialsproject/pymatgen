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
    def __new__(cls, input_matrix):
        """
        Constructs the Stress object similarly to other constructors
        involving SQTensor subclasses
        """
        obj = SQTensor(input_matrix).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    def __repr__(self):
        return "Stress({})".format(self.__str__())

    @property
    def deviator_principal_invariants(self):
        # TODO: check sign convention
        """
        returns the principal invariants of the deviatoric stress tensor, 
        which is calculated by finding the coefficients of the characteristic
        polynomial of the stress tensor minus the identity times the mean
        stress
        """
        return self.deviator_stress.principal_invariants
        
    @property
    def von_mises(self):
        """
        calculates the von mises stress
        """
        if not self.is_symmetric():
            raise ValueError("The stress tensor is not symmetric, Von Mises "
                             "stress is based on a symmetric stress tensor.")
        return math.sqrt(3*self.deviator_principal_invariants[1])

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
        return f.det*self*f.I.T

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
        return f.det*f.I*self*f.I.T

    @property
    def voigt(self):
        """
        returns the vector representing to the stress tensor in voigt notation
        """
        return [self[ind] for ind in voigt_map]

if __name__ == "__main__":

    mat = np.eye(3)
    mat = np.random.randn(3, 3)
#    mat[0,2] = 0.1
#    mat[2,0] = 0.1
    s = Stress(mat)
    print s.principal_invariants
    
#    for property, value in vars(s).iteritems():
#            print property, ": ", value


#    print s.get_scaled(2.0)
#    s.get_scaled(1.2)
#    print s.is_symmetric()
#    print s.value(1,2)
#    print s.returntensor()
#    print s.stress_matrix
#    print s.PiolaKirchoff1(mat)
#    print s.PiolaKirchoff2(mat)
