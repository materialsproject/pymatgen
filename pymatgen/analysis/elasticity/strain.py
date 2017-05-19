# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

"""
This module provides classes and methods used to describe deformations and
strains, including applying those deformations to structure objects and
generating deformed structure sets for further calculations.
"""

from pymatgen.core.lattice import Lattice
from pymatgen.analysis.elasticity.tensors import SquareTensor, symmetry_reduce
import numpy as np
from six.moves import zip
import collections

__author__ = "Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Joseph Montoya, Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Development"
__date__ = "March 13, 2012"


class Deformation(SquareTensor):
    """
    Subclass of SquareTensor that describes the deformation gradient tensor
    """

    def __new__(cls, deformation_gradient):
        """
        Create a Deformation object.  Note that the constructor uses __new__
        rather than __init__ according to the standard method of subclassing
        numpy ndarrays.

        Args:
            deformation_gradient (3x3 array-like): the 3x3 array-like
                representing the deformation gradient
        """
        obj = super(Deformation, cls).__new__(cls, deformation_gradient)
        return obj.view(cls)

    def is_independent(self, tol=1e-8):
        """
        checks to determine whether the deformation is independent
        """
        return len(self.get_perturbed_indices(tol)) == 1
    
    def get_perturbed_indices(self, tol=1e-8):
        """
        Gets indices of perturbed elements of the deformation gradient,
        i. e. those that differ from the identity
        """
        indices = list(zip(*np.where(abs(self - np.eye(3)) > tol)))
        return indices

    @property
    def green_lagrange_strain(self):
        """
        calculates the euler-lagrange strain from
        the deformation gradient
        """
        return Strain.from_deformation(self)

    def apply_to_structure(self, structure):
        """
        Apply the deformation gradient to a structure.

        Args:
            structure (Structure object): the structure object to
                be modified by the deformation
        """
        def_struct = structure.copy()
        def_struct.modify_lattice(Lattice(np.dot(def_struct.lattice.matrix,
                                                 self)))
        return def_struct

    @classmethod
    def from_index_amount(cls, matrixpos, amt):
        """
        Factory method for constructing a Deformation object
        from a matrix position and amount

        Args:
            matrixpos (tuple): tuple corresponding the matrix position to
                have a perturbation added
            amt (float): amount to add to the identity matrix at position
                matrixpos
        """
        f = np.identity(3)
        f[matrixpos] += amt
        return cls(f)


class DeformedStructureSet(collections.Sequence):
    """
    class that generates a set of independently deformed structures that
    can be used to calculate linear stress-strain response
    """

    def __init__(self, rlxd_str, nd=0.01, ns=0.08,
                 num_norm=4, num_shear=4, symmetry=False):
        """
        constructs the deformed geometries of a structure.  Generates
        m + n deformed structures according to the supplied parameters.

        Args:
            rlxd_str (structure): structure to undergo deformation, if
                fitting elastic tensor is desired, should be a geometry
                optimized structure
            nd (float): maximum perturbation applied to normal deformation
            ns (float): maximum perturbation applied to shear deformation
            num_norm (int): number of deformation structures to generate for
                normal deformation, must be even
            num_shear (int): number of deformation structures to generate for
                shear deformation, must be even
        """

        if num_norm % 2 != 0:
            raise ValueError("Number of normal deformations (num_norm)"
                             " must be even.")
        if num_shear % 2 != 0:
            raise ValueError("Number of shear deformations (num_shear)"
                             " must be even.")

        norm_deformations = np.linspace(-nd, nd, num=num_norm + 1)
        norm_deformations = norm_deformations[norm_deformations.nonzero()]
        shear_deformations = np.linspace(-ns, ns, num=num_shear + 1)
        shear_deformations = shear_deformations[shear_deformations.nonzero()]

        self.undeformed_structure = rlxd_str
        self.deformations = []
        self.def_structs = []

        # Generate deformations
        for ind in [(0, 0), (1, 1), (2, 2)]:
            for amount in norm_deformations:
                defo = Deformation.from_index_amount(ind, amount)
                self.deformations.append(defo)

        for ind in [(0, 1), (0, 2), (1, 2)]:
            for amount in shear_deformations:
                defo = Deformation.from_index_amount(ind, amount)
                self.deformations.append(defo)

        # Perform symmetry reduction if specified
        if symmetry:
            self.sym_dict = symmetry_reduce(self.deformations, self.undeformed_structure)
            self.deformations = list(self.sym_dict.keys())
        self.def_structs = [defo.apply_to_structure(rlxd_str)
                            for defo in self.deformations]

    def __iter__(self):
        return iter(self.def_structs)
    
    def __len__(self):
        return len(self.def_structs)

    def __getitem__(self, ind):
        return self.def_structs[ind]

    def as_strain_dict(self):
        """
        Returns dictionary of deformed structures indexed by independent
        strain objects in accordance with legacy behavior of elasticity
        package
        """
        strains = [IndependentStrain(defo) for defo in self.deformations]
        return dict(zip(strains, self.def_structs))


class Strain(SquareTensor):
    """
    Subclass of SquareTensor that describes the Green-Lagrange strain tensor.
    """

    def __new__(cls, strain_matrix, dfm=None):
        """
        Create a Strain object.  Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.  Note also that the default constructor
        does not include the deformation gradient

        Args:
            strain_matrix (3x3 array-like): the 3x3 array-like
                representing the Green-Lagrange strain
        """
        vscale = np.ones((6,))
        vscale[3:] *= 2
        obj = super(Strain, cls).__new__(cls, strain_matrix, vscale=vscale)
        if dfm is None:
            obj._dfm = convert_strain_to_deformation(obj)
        else:
            dfm = Deformation(dfm)
            gls_test = 0.5 * (np.dot(dfm.trans, dfm) - np.eye(3))
            if (gls_test - obj > 1e-10).any():
                raise ValueError("Strain and deformation gradients "
                                 "do not match!")
            obj._dfm = Deformation(dfm)

        if not obj.is_symmetric():
            raise ValueError("Strain objects must be initialized "
                             "with a symmetric array or a voigt-notation "
                             "vector with six entries.")
        return obj.view(cls)

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.rank = getattr(obj, "rank", None)
        self._dfm = getattr(obj, "_dfm", None)
        self._vscale = getattr(obj, "_vscale", None)

    @classmethod
    def from_deformation(cls, deformation):
        """
        Factory method that returns a Strain object from a deformation
        gradient

        Args:
            deformation (3x3 array-like):
        """
        dfm = Deformation(deformation)
        return cls(0.5 * (np.dot(dfm.trans, dfm) - np.eye(3)), dfm)

    @property
    def deformation_matrix(self):
        """
        returns the deformation matrix
        """
        return self._dfm

    @property
    def von_mises_strain(self):
        """
         Equivalent strain to Von Mises Stress
        """
        eps = self - 1/3*np.trace(self)*np.identity(3)

        return np.sqrt(np.dot(eps.voigt, eps.voigt)*2/3)


class IndependentStrain(Strain):
    """
    Class for independent strains intended for use with old Materials Project
    elasticity workflow.  Note that the default constructor constructs from
    a deformation matrix, rather than an array representing the strain, to
    emulate the legacy behavior.
    """

    def __new__(cls, deformation_gradient, tol=1e-8):
        """
        Create an Independent Strain object.  Note that the constructor uses
        __new__ rather than __init__ according to the standard method of
        subclassing numpy ndarrays.  Note also that, unlike the Strain class,
        the default constructor of IndependentStrain takes the deformation
        gradient as input, rather than an array representing the Green-Lagrange
        strain.

        Args:
            deformation_gradient (3x3 array-like): the 3x3 array-like
                representing the deformation gradient
        """
        deformation_gradient = Deformation(deformation_gradient)
        if not deformation_gradient.is_independent(tol):
            raise ValueError("IndependentStrain must be constructed from "
                             "an independent deformation gradient")
        obj = Strain.from_deformation(deformation_gradient).view(cls)
        (obj._i, obj._j) = obj._dfm.get_perturbed_indices(tol=tol)[0]
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._dfm = getattr(obj, "_dfm", None)
        self._i = getattr(obj, "_i", None)
        self._j = getattr(obj, "_j", None)

    @property
    def ij(self):
        """
        Convenience method to return independent indices
        """
        return self._i, self._j


def convert_strain_to_deformation(strain):
    """
    This function converts a strain to a deformation gradient that will
    produce that strain
    
    Args:
        strain (3x3 array-like): strain matrix
    """
    strain = SquareTensor(strain)
    ftdotf = 2*strain + np.eye(3)
    eigs, eigvecs = np.linalg.eigh(ftdotf)
    rotated = ftdotf.rotate(np.transpose(eigvecs))
    rotated = rotated.round(10)
    defo = Deformation(np.sqrt(rotated))
    result = defo.rotate(eigvecs)
    return result
