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
from pymatgen.transformations.standard_transformations import *
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


class Deformation(SQTensor):
    """
    Subclass of SQTensor that describes the deformation gradient tensor
    """

    def __new__(cls, stress_matrix, dfm=None):
        obj = SQTensor(deformation_gradient).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    def __repr__(self):
        return "Deformation({})".format(self.__str__())
 
    @property
    def check_independent(self, tol=0.00001):
        """
        a check to determine whether the deformation matrix represents an
        independent deformation, raises a ValueError if not.  If so, returns 
        the indices of the deformation gradient entry representing the 
        independent component

        Args: tol
        """
        indices = zip(*np.asarray(self - np.eye(3)).nonzero())
        if len(indices) != 1:
            raise ValueError("One and only one independent deformation\
                             must be applied.")

        return indices[0]

    @property
    def euler_lagrange_strain(self):
        """
        calculates the euler-lagrange strain from
        the deformation gradient
        """
        return 0.5*self*self.T - np.eye(3)

    @classmethod
    def from_ind_amt_dfm(matrixpos, amt):
        """
        Factory method for constructing a Deformation object
        from a matrix position and amount
        """
        F = np.identity(3)
        F[matrixpos] = F[matrixpos] + amt
        return cls(F)


class DeformedStructureSet(object):
    """
    class that generates a set of deformed structures that
    can be used to fit the elastic tensor of a material
    """
    def __init__(self,rlxd_str, nd=0.01, ns=0.08, m=4, n=4, 
                 delete_center=True, symmetry=False):
        """
        constructs the deformed geometries of a structure.  Generates
        m + n deformed structures according to the supplied parameters.

        Args:
            rlxd_str (structure): structure to undergo deformation
            nd (float): maximum perturbation applied to normal deformation
            ns (float): maximum perturbation applied to shear deformation
            m (int): number of deformation structures to generate for 
                normal deformation
            n (int): number of deformation structures to generate for 
                shear deformation
        """

        # TODO: JHM suggests that we might want to think about the
        #           conventions specified here
        # TODO: JHM needs to figure out how to do symmetry

        if m%2 != 0:
            raise ValueError("m has to be even.")
        if n%2 != 0:
            raise ValueError("n has to be even.")
        
        self.msteps = np.int(m)
        self.nsteps = np.int(n)
        
        normal_defs = np.linspace(-nd, nd, num=m+1)
        if delete_center:
            normal_defs = np.delete(normal_defs, np.int(m/2), 0)
        shear_defs = np.linspace(-ns, ns, num=n+1)
        if delete_center:
            sheardef = np.delete(sheardef, np.int(n/2), 0)

        self.deformed_structures = {}

        # TODO: Make this section more pythonic
        # TODO: JHM asks whether indexing the dictionary
        #           on the strain object is efficient?

        if symmetry:
            raise NotImplementedError("Symmetry reduction of structure "\
                                      "generation is not yet implemented")
        else:
            # Apply normal deformations
            self.symmetry = None
            for i1 in range(0, 3):
                for i2 in range(0, len(defs)):
                    s=rlxd_str.copy()
                    F = np.identity(3)
                    F[i1, i1] = F[i1, i1] + normal_defs[i2]
                    StrainObject = IndependentStrain.from_deformation(F)
                    s.apply_deformation_gradient(F)
                    defstructures[StrainObject] = s

            # Apply shear deformations 
            F_index = [[0, 1], [0, 2], [1, 2]]
            for j1 in range(0, 3):
                for j2 in range(0, len(sheardef)):
                    s=rlxd_str.copy()
                    F = np.identity(3)
                    F[F_index[j1][0], F_index[j1][1]] = F[F_index[j1][0], F_index[j1][1]] + sheardef[j2]
                    StrainObject = IndependentStrain.from_deformation(F)
                    s.apply_deformation_gradient(F)
                    defstructures[StrainObject] = s

        return defstructures



if __name__ == "__main__":

    mat = np.eye(3)
    mat[0,1] = 0.001
#    print mat

    my_strain = IndependentStrain.from_deformation(mat)
    my_strain.check_F()


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


