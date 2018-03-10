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
import scipy
import itertools
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
        old_latt = def_struct.lattice.matrix
        new_latt = np.transpose(np.dot(self, np.transpose(old_latt)))
        def_struct.modify_lattice(Lattice(new_latt)) 
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

    def __init__(self, structure, norm_strains=None, shear_strains=None,
                 symmetry=False):
        """
        constructs the deformed geometries of a structure.  Generates
        m + n deformed structures according to the supplied parameters.

        Args:
            structure (Structure): structure to undergo deformation
            norm_strains (list of floats): strain values to apply
                to each normal mode.
            shear_strains (list of floats): strain values to apply
                to each shear mode.
            symmetry (bool): whether or not to use symmetry reduction.
        """
        norm_strains = norm_strains or [-0.01, -0.005, 0.005, 0.01]
        shear_strains = shear_strains or [-0.06, -0.03, 0.03, 0.06]

        self.undeformed_structure = structure
        self.deformations = []
        self.def_structs = []

        # Generate deformations
        for ind in [(0, 0), (1, 1), (2, 2)]:
            for amount in norm_strains:
                strain = Strain.from_index_amount(ind, amount)
                self.deformations.append(strain.deformation_matrix)

        for ind in [(0, 1), (0, 2), (1, 2)]:
            for amount in shear_strains:
                strain = Strain.from_index_amount(ind, amount)
                self.deformations.append(strain.deformation_matrix)

        # Perform symmetry reduction if specified
        if symmetry:
            self.sym_dict = symmetry_reduce(self.deformations, structure)
            self.deformations = list(self.sym_dict.keys())
        self.deformed_structures = [defo.apply_to_structure(structure)
                                    for defo in self.deformations]

    def __iter__(self):
        return iter(self.deformed_structures)
    
    def __len__(self):
        return len(self.deformed_structures)

    def __getitem__(self, ind):
        return self.deformed_structures[ind]


class Strain(SquareTensor):
    """
    Subclass of SquareTensor that describes the Green-Lagrange strain tensor.
    """

    def __new__(cls, strain_matrix, dfm=None, dfm_shape="upper"):
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
            obj._dfm = convert_strain_to_deformation(obj, dfm_shape)
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

    @classmethod
    def from_index_amount(cls, idx, amount):
        """
        Like Deformation.from_index_amount, except generates
        a strain from the zero 3x3 tensor or voigt vector with
        the amount specified in the index location.  Ensures
        symmetric strain.

        Args:
            idx (tuple or integer): index to be perturbed, can be voigt or 
                full-tensor notation
            amount (float): amount to perturb selected index
        """
        if np.array(idx).ndim == 0:
            v = np.zeros(6)
            v[idx] = amount
            return cls.from_voigt(v)
        elif np.array(idx).ndim == 1:
            v = np.zeros((3, 3))
            for i in itertools.permutations(idx):
                v[i] = amount
            return cls(v)
        else:
            raise ValueError("Index must either be 2-tuple or integer "
                             "corresponding to full-tensor or voigt index")

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

        return np.sqrt(np.sum(eps * eps) * 2/3)


def convert_strain_to_deformation(strain, shape="upper"):
    """
    This function converts a strain to a deformation gradient that will
    produce that strain.  Supports three methods:
    
    Args:
        strain (3x3 array-like): strain matrix
        shape: (string): method for determining deformation, supports
            "upper" produces an upper triangular defo
            "lower" produces a lower triangular defo
            "symmetric" produces a symmetric defo
    """
    strain = SquareTensor(strain)
    ftdotf = 2*strain + np.eye(3)
    if shape == "upper":
        result = scipy.linalg.cholesky(ftdotf)
    elif shape == "symmetric":
        result = scipy.linalg.sqrtm(ftdotf)
    else:
        raise ValueError("shape must be \"upper\" or \"symmetric\"")
    return Deformation(result)
