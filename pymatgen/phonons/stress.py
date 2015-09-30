import warnings, sys, os
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
from pymatgen.phonons.tensors import SQTensor
import numpy as np

__author__="Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ ="March 22, 2012"

class Stress(SQTensor):
    """
    This class extends SQTensor
    """
    def __init__(self, stress_matrix):
        super(Stress, self).__init__(stress_matrix)

    @property
    def deviator_principal_invariants(self):
        # TODO: check sign convention
        # TODO: Also JM says this may be too abstract and might want
        #   to use the original convention here, old code left
        """
        returns the principal invariants of the deviatoric stress tensor, 
        which is calculated by finding the coefficients of the characteristic
        polynomial of the stress tensor minus the identity times the mean
        stress
        """
        return self.deviator_stress.principal_invariants
        """
        I = self.principal_invariants
        J1 = 0
        J2 = 1./3. * I[0]**2 - I[1]
        J3 = 2./27. * I[0]**3 - 1./3.*I[0]*I[1] + I[2]
        J = [J1, J2, J3]
        return J
        """

    @property
    def stress_matrix(self):
        return self._matrix

    @property
    def von_mises(self):
        """
        calculates the von mises stress
        """
        if self.is_symmetric() == False:
            raise ValueError("The stress tensor is not symmetric, \
                             VM stress is based on a symmetric stress tensor.")
        return math.sqrt(3*self.deviator_principal_invariants[1])

    @property
    def mean_stress(self):
        """
        returns the mean stress
        """
        return 1./3.*self._matrix.trace()

    @property
    def deviator_stress(self):
        """
        returns the deviatoric component of the stress
        """
        if self.is_symmetric() == False:
            raise ValueError("The stress tensor is not symmetric, \
                             so deviator stress will not be either")
        return SQTensor(self._matrix - self.mean_stress*np.eye(3))

    # TODO: JM asks is there a more descriptive way to distinguish these?
    # TODO: JM asks is the F argument here necessary or should it operate
    #   on the matrix attribute?
    def piola_kirchoff_1(self, F):
        if self.is_symmetric() == False:
            raise ValueError("The stress tensor is not symmetric, \
                             PK stress is based on a symmetric stress tensor.")
        return np.linalg.det(F)*self._matrix*np.transpose(np.linalg.inv(F))

    def piola_kirchoff_2(self, F):
        if self.is_symmetric() == False:
            raise ValueError("The stress tensor is not symmetric, \
                             PK stress is based on a symmetric stress tensor.")
        return np.linalg.det(F)*np.linalg.inv(F)*self._sigma*np.transpose(np.linalg.inv(F))


if __name__ == "__main__":

    mat = np.eye(3)
    mat = np.random.randn(3,3)
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


