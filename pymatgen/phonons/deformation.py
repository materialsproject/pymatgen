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
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.phonons.strain import IndependentStrain,Strain
from pymatgen.phonons.tensors import voigt_map
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

    def __new__(cls, deformation_gradient, dfm=None):
        obj = SQTensor(deformation_gradient).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    def __repr__(self):
        return "Deformation({})".format(self.__str__())
 
    @classmethod
    def check_independent(self):
        """
        a check to determine whether the deformation matrix represents an
        independent deformation, raises a ValueError if not.  If so, returns 
        the indices of the deformation gradient entry representing the 
        independent component

        Args: tol
        """
        indices = zip(*np.asarray(self - np.eye(3)).nonzero())
        if len(indices) != 1:
            raise ValueError("One and only one independent deformation"\
                             "must be applied.")
        return indices[0]

    @property
    def euler_lagrange_strain(self):
        """
        calculates the euler-lagrange strain from
        the deformation gradient
        """
        return Strain.from_deformation(self)

    @classmethod
    def apply_to_structure(self,structure):
        """
        Apply the deformation gradient to a structure.
        
        Args:
            structure (Structure object): the structure object to
                be modified by the deformation
        """
        def_struct = structure.copy()
        def_struct.modify_lattice(Lattice(np.dot(self._lattice.matrix,
                                                 np.array(self))))
        return def_struct

    @classmethod
    def from_index_amount(cls,matrixpos, amt):
        """
        Factory method for constructing a Deformation object
        from a matrix position and amount

        Args:
            matrixpos (tuple): tuple corresponding the matrix position to
                have a perturbation added
            amt (float): amount to add to the identity matrix at position 
                matrixpos
        """
        F = np.identity(3)
        F[matrixpos] += amt
        return cls(F)

class DeformedStructureSet(object):
    """
    class that generates a set of deformed structures that
    can be used to fit the elastic tensor of a material
    """

    def __init__(self,rlxd_str, nd=0.01, ns=0.08, num_norm=4, num_shear=4, symmetry=False):
        """
        constructs the deformed geometries of a structure.  Generates
        m + n deformed structures according to the supplied parameters.

        Args:
            rlxd_str (structure): structure to undergo deformation, if 
                fitting elastic tensor is desired, should be a geometry 
                optimized structure
            nd (float): maximum perturbation applied to normal deformation
            ns (float): maximum perturbation applied to shear deformation
            m (int): number of deformation structures to generate for 
                normal deformation
            n (int): number of deformation structures to generate for 
                shear deformation
        """

        if num_norm%2 != 0:
            raise ValueError("Number of normal deformations (num_norm)"\
                             " must be even.")
        if num_shear%2 != 0:
            raise ValueError("Number of shear deformations (num_shear)"\
                             " must be even.")
        
        self.msteps = np.int(m)
        self.nsteps = np.int(n)
        
        norm_deformations = np.linspace(-nd, nd, num=m+1)
        norm_deformations = norm_deformations[norm_deformations.nonzero()]
        shear_deformations = np.linspace(-ns, ns, num=n+1)
        shear_deformations = shear_deformations[shear_deformations.nonzero()]

        self.equilibrium_structure = rlxd_str
        self.deformations = []
        self.def_structs = []
        if symmetry:
            raise NotImplementedError("Symmetry reduction of structure "\
                                      "generation is not yet implemented") 
        else:
            self.symmetry = None
            # Determine normal deformation gradients
            # Apply normal deformations
            for ind in [(1,1),(2,2),(3,3)]:
                for amount in normal_deformations:
                    defo_tensor = Deformation.from_index_amount(ind,amount)
                    self.deformations.append(defo_tensor)
                    self.def_structs.append(defo_tensor.apply_to_structure(s))

            # Apply shear deformations 
            for ind in [(0,1),(0,2),(1,2)]:
                for amount in shear_deformations:
                    defo_tensor = Deformation.from_index_amount(ind,amount)
                    self.def_structs.append(defo_tensor.apply_to_structure(s))

        def as_stress_strain_dict(self):

if __name__ == "__main__":
    from pymatgen.matproj.rest import MPRester
    mpr = MPRester()
    Cu_struct = mpr.get_structures('Cu')[0]
    dss = DeformedStructureSet(Cu_struct)
