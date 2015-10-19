import warnings
import sys
import unittest
import pymatgen
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp import Vasprun
from pymatgen.io.cif import CifWriter
from pymatgen.io.cif import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.phonons.tensors import SQTensor
from pymatgen.phonons.deformation import Deformation
from pymatgen.transformations.standard_transformations import *
import warnings
import numpy as np
import os

__author__ = "Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ = "March 13, 2012"


class Strain(SQTensor):
    """
    Subclass of SQTensor that describes the strain tensor
    """
        
    def __new__(cls, strain_matrix, dfm=None):
        obj = SQTensor(strain_matrix).view(cls)
        obj._dfm = dfm
        if dfm == None:
            warnings.warn("Constructing a strain object without a deformation "\
                          "matrix makes many methods unusable.  Use "\
                          "Strain.from_deformation to construct a Strain object"\
                          " from a deformation gradient.")
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._dfm = getattr(obj, "_dfm", None)

    def __repr__(self):
        return "Strain({})".format(self.__str__())

    @classmethod
    def from_deformation(cls, deformation):
        """
        constructor that returns a Strain object from a deformation
        gradient

        Args:
            deformation (3x3 array-like):
        """
        dfm = Deformation(deformation)
        return cls(0.5*(dfm.T*dfm - np.eye(3)), dfm)

    @property
    def deformation_matrix(self):
        """
        returns the deformation matrix
        """
        return self._dfm
    
    @property
    def independent_deformation(self):
        """
        determines whether the deformation matrix represents an
        independent deformation, raises a value error if not.
        Returns the index of the deformation gradient corresponding
        to the independent deformation

        Args: tol
        """
        if self._dfm == None:
            raise ValueError("No deformation matrix supplied "\
                             "for this strain tensor.") 
        return self._dfm.check_independent()

    @property
    def voigt(self):
        """
        translates a strain tensor into a voigt notation vector
        """
        return [self[0,0],self[1,1],self[2,2],
                2.*self[1,2],2.*self[0,2],2.*self[0,1]]

class IndependentStrain(Strain):
    """
    Class for independent strains intended for use with old Materials Project
    elasticity workflow.  Note that the default constructor constructs from 
    a deformation matrix, rather than an array representing the strain, to 
    emulate the old workflow behavior.
    """
    def __new__(cls, deformation_matrix):
        obj = Strain.from_deformation(deformation_matrix).view(cls)
        (obj._i,obj._j) = obj.independent_deformation
        return obj

    @property
    def i(self):
        return self._i

    @property
    def j(self):
        return self._j


if __name__ == "__main__":

    mat = np.eye(3)
    mat[0,1] = 0.001
#    print mat

    my_strain = IndependentStrain.from_deformation(mat)
    #my_strain.check_F()


#    print my_strain._strain
    
    
    
#    print type(mat)

#    print my_strain.deformation_matrix
#    print my_strain.strain

#    my_strain2 = IndependentStrain(mat)
#    print my_strain2.__dict__.keys()
#    print my_strain2.__hash__()

#    print my_strain2._j
#    print my_strain2.check_F()
#    my_strain2.checkF
#    print my_strain.__dict__.keys()
#    print my_strain.deformation_matrix
#    print my_strain.strain
#    my_strain.index
#    my_scaled_strain = my_strain.get_scaled(1.05)
#    print my_scaled_strain.deformation_matrix
#    print my_scaled_strain.strain
#    print my_strain == my_scaled_strain
#    mat2 = np.eye(3)
#    mat2[0,0] = 1.01
#    my_strain2 = Strain(mat)
#    print my_strain == my_strain2


