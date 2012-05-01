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
    #TODO: AJ says method names are not PEP8. e.g., VonMises should be von_mises.
 
    def __init__(self, stress_matrix):
        self._sigma = stress_matrix
        super(Stress, self).__init__(self._sigma)
        #TODO: AJ says defining a sigma parameter is bad, just use the _matrix parameter from the superclass

    def get_scaled(self, scale_factor):
        return super(Stress, self).get_scaled(scale_factor)

    @property
    def DeviatorPrincipalInvariants(self):
        I = self.PrincipalInvariants
        J1 = 0
        J2 = 1.0/3*I[0]**2 - I[1]
        J3 = 2.0/27*I[0]**3 - 1.0/3*I[0]*I[1] + I[2]
        J = [J1, J2, J3]
        return J

    @property
    def VonMises(self):
        if self.is_symmetric() == False:
            raise ValueError("The stress tensor is not symmetric, VM stress is based on a symmetric stress tensor.")
        J = self.DeviatorPrincipalInvariants        
        sigma_mises = math.sqrt(3*J[1])
        return sigma_mises

    @property
    def MeanStress(self):
        if self.is_symmetric() == False:
            raise ValueError("The stress tensor is not symmetric, mean stress is based on a symmetric stress tensor.")
        return 1.0/3*(self._sigma[0,0] + self._sigma[1,1] + self._sigma[2,2])

    @property
    def DeviatorStress(self):
        if self.is_symmetric() == False:
            raise ValueError("The stress tensor is not symmetric, so deviator stress will not be either")
        return self._sigma - self.MeanStress

    def PiolaKirchoff1(self, F):
        if self.is_symmetric() == False:
            raise ValueError("The stress tensor is not symmetric, PK stress is based on a symmetric stress tensor.")
        return np.linalg.det(F)*self._sigma*np.transpose(np.linalg.inv(F))

    def PiolaKirchoff2(self, F):
        if self.is_symmetric() == False:
            raise ValueError("The stress tensor is not symmetric, PK stress is based on a symmetric stress tensor.")
        return np.linalg.det(F)*np.linalg.inv(F)*self._sigma*np.transpose(np.linalg.inv(F))


if __name__ == "__main__":

    mat = np.eye(3)
    mat = np.random.randn(3,3)
#    mat[0,2] = 0.1
#    mat[2,0] = 0.1
    s = Stress(mat)
    print s.PrincipalInvariants
    
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


