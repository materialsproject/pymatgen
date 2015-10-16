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
from pymatgen.phonons.strain import IndependentStrain
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
    def from_index_amount(matrixpos, amt):
        """
        Factory method for constructing a Deformation object
        from a matrix position and amount

        Args:
            matrixpos (tuple): 
        """
        F = np.identity(3)
        F[matrixpos] += amt
        return cls(F)

# TODO: should this be PMGSONABLE?
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
        normal_defs = np.delete(normal_defs, np.int(m/2), 0)
        shear_defs = np.linspace(-ns, ns, num=n+1)
        shear_defs = np.delete(shear_defs, np.int(n/2), 0)

        self.equilibrium_structure = rlxd_str
        self.deformed_structures = {}

        # TODO: Make this section more pythonic
        # TODO: JHM asks whether indexing the dictionary
        #           on the strain object is efficient?
        # TODO: Integrate new deformatinos class
        if symmetry:
            raise NotImplementedError("Symmetry reduction of structure "\
                                      "generation is not yet implemented")
            '''
            recp_lattice = rlxd_str.lattice.reciprocal_lattice_crystallographic
            recp_lattice = recp_lattice.scale(1)
            recp = Structure(recp_lattice, ["H"], [[0,0,0]])
            analyzer = SpacegroupAnalyzer(recp, symprec=0.001)
            symm_ops = analyzer.get_symmetry_operations()
            self.symmetry = symm_ops
            # generate a list of unique deformations
            unique_defs = []
            for i1 in range(0, 3):
                for i2 in range(0, len(normal_defs)):
                    F = np.eye(3)
                    F[i1, i1] = F[i1, i1] + normal_defs[i2]
                    #import pdb; pdb.set_trace()
                    #print F
                    if not [f for f in [op.operate_multi(F) for op in symm_ops]\
                            if f in unique_defs]:
                        unique_defs += [F]
            for j1 in [(0,1),(0,2),(1,2)]:
                for j2 in range(0, len(shear_defs)):
                    F = np.eye(3)
                    F[j1] = F[j1] + shear_defs[j2]
                    #import pdb; pdb.set_trace()
                    if not [f for f in [op.operate_multi(F) for op in symm_ops]\
                            if f in unique_defs]:
                        unique_defs += [F]
            for unique_def in unique_defs:
                s = rlxd_str.copy()
                StrainObject = IndependentStrain.from_deformation(F)
                s.apply_deformation_gradient(unique_def)
                self.deformed_structures[StrainObject] = s
            '''
        else:
            # Apply normal deformations
            self.symmetry = None
            for i1 in range(0, 3):
                for i2 in range(0, len(normal_defs)):
                    s=rlxd_str.copy()
                    F = np.identity(3)
                    F[i1, i1] += normal_defs[i2]
                    StrainObject = IndependentStrain.from_deformation(F)
                    s.apply_deformation_gradient(F)
                    self.deformed_structures[StrainObject] = s

            # Apply shear deformations 
            F_index = [[0, 1], [0, 2], [1, 2]]
            for j1 in [(0,1),(0,2),(1,2)]:
                for j2 in range(0, len(shear_defs)):
                    s=rlxd_str.copy()
                    F = np.identity(3)
                    F[j1] += shear_defs[j2]
                    StrainObject = IndependentStrain.from_deformation(F)
                    s.apply_deformation_gradient(F)
                    self.deformed_structures[StrainObject] = s

def list_intersect(list_1,list_2):
    return None 
    
if __name__ == "__main__":
    from pymatgen.matproj.rest import MPRester
    mpr = MPRester()
    Cu_struct = mpr.get_structures('Cu')[0]
    dss = DeformedStructureSet(Cu_struct)
