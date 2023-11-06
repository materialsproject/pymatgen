"""
This module provides classes and methods used to describe deformations and
strains, including applying those deformations to structure objects and
generating deformed structure sets for further calculations.
"""

from __future__ import annotations

import collections
import itertools
from typing import TYPE_CHECKING, Literal

import numpy as np
import scipy

from pymatgen.core.lattice import Lattice
from pymatgen.core.tensors import SquareTensor, symmetry_reduce

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from pymatgen.core.structure import Structure

__author__ = "Joseph Montoya"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Maarten de Jong, Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Production"
__date__ = "July 24, 2018"


class Deformation(SquareTensor):
    """Subclass of SquareTensor that describes the deformation gradient tensor."""

    symbol = "d"

    def __new__(cls, deformation_gradient):
        """
        Create a Deformation object. Note that the constructor uses __new__ rather than
        __init__ according to the standard method of subclassing numpy ndarrays.

        Args:
            deformation_gradient (3x3 array-like): the 3x3 array-like
                representing the deformation gradient
        """
        obj = super().__new__(cls, deformation_gradient)
        return obj.view(cls)

    def is_independent(self, tol: float = 1e-8):
        """Checks to determine whether the deformation is independent."""
        return len(self.get_perturbed_indices(tol)) == 1

    def get_perturbed_indices(self, tol: float = 1e-8):
        """
        Gets indices of perturbed elements of the deformation gradient,
        i. e. those that differ from the identity.
        """
        return list(zip(*np.where(abs(self - np.eye(3)) > tol)))

    @property
    def green_lagrange_strain(self):
        """Calculates the Euler-Lagrange strain from the deformation gradient."""
        return Strain.from_deformation(self)

    def apply_to_structure(self, structure: Structure):
        """
        Apply the deformation gradient to a structure.

        Args:
            structure (Structure object): the structure object to
                be modified by the deformation
        """
        def_struct = structure.copy()
        old_latt = def_struct.lattice.matrix
        new_latt = np.transpose(np.dot(self, np.transpose(old_latt)))
        def_struct.lattice = Lattice(new_latt)
        return def_struct

    @classmethod
    def from_index_amount(cls, matrixpos, amt):
        """
        Factory method for constructing a Deformation object
        from a matrix position and amount.

        Args:
            matrixpos (tuple): tuple corresponding the matrix position to
                have a perturbation added
            amt (float): amount to add to the identity matrix at position
                matrixpos
        """
        f = np.identity(3)
        f[matrixpos] += amt
        return cls(f)


class DeformedStructureSet(collections.abc.Sequence):
    """
    class that generates a set of independently deformed structures that
    can be used to calculate linear stress-strain response.
    """

    def __init__(self, structure: Structure, norm_strains=None, shear_strains=None, symmetry=False):
        """
        Construct the deformed geometries of a structure. Generates m + n deformed structures
        according to the supplied parameters.

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
        self.deformations: list[Deformation] = []
        self.def_structs: list[Structure] = []

        # Generate deformations
        for ind in [(0, 0), (1, 1), (2, 2)]:
            for amount in norm_strains:
                strain = Strain.from_index_amount(ind, amount)
                self.deformations.append(strain.get_deformation_matrix())

        for ind in [(0, 1), (0, 2), (1, 2)]:
            for amount in shear_strains:
                strain = Strain.from_index_amount(ind, amount)
                self.deformations.append(strain.get_deformation_matrix())

        # Perform symmetry reduction if specified
        if symmetry:
            self.sym_dict = symmetry_reduce(self.deformations, structure)
            self.deformations = list(self.sym_dict)
        self.deformed_structures = [defo.apply_to_structure(structure) for defo in self.deformations]

    def __iter__(self):
        return iter(self.deformed_structures)

    def __len__(self):
        return len(self.deformed_structures)

    def __getitem__(self, ind):
        return self.deformed_structures[ind]


class Strain(SquareTensor):
    """Subclass of SquareTensor that describes the Green-Lagrange strain tensor."""

    symbol = "e"

    def __new__(cls, strain_matrix):
        """
        Create a Strain object. Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays. Note also that the default constructor
        does not include the deformation gradient.

        Args:
            strain_matrix (ArrayLike): 3x3 matrix or length-6 Voigt notation vector
                representing the Green-Lagrange strain
        """
        vscale = np.ones((6,))
        vscale[3:] *= 2
        obj = super().__new__(cls, strain_matrix, vscale=vscale)
        if not obj.is_symmetric():
            raise ValueError(
                "Strain must be initialized with a symmetric array or a Voigt-notation vector with six entries."
            )
        return obj.view(cls)

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.rank = getattr(obj, "rank", None)
        self._vscale = getattr(obj, "_vscale", None)

    @classmethod
    def from_deformation(cls, deformation: ArrayLike) -> Strain:
        """
        Factory method that returns a Strain object from a deformation
        gradient.

        Args:
            deformation (ArrayLike): 3x3 array defining the deformation
        """
        dfm = Deformation(deformation)
        return cls(0.5 * (np.dot(dfm.trans, dfm) - np.eye(3)))

    @classmethod
    def from_index_amount(cls, idx, amount):
        """
        Like Deformation.from_index_amount, except generates
        a strain from the zero 3x3 tensor or Voigt vector with
        the amount specified in the index location. Ensures
        symmetric strain.

        Args:
            idx (tuple or integer): index to be perturbed, can be Voigt or full-tensor notation
            amount (float): amount to perturb selected index
        """
        if np.array(idx).ndim == 0:
            v = np.zeros(6)
            v[idx] = amount
            return cls.from_voigt(v)
        if np.array(idx).ndim == 1:
            v = np.zeros((3, 3))
            for i in itertools.permutations(idx):
                v[i] = amount
            return cls(v)
        raise ValueError("Index must either be 2-tuple or integer corresponding to full-tensor or Voigt index")

    def get_deformation_matrix(self, shape: Literal["upper", "lower", "symmetric"] = "upper"):
        """
        Returns the deformation matrix.

        Args:
            shape ('upper' | 'lower' | 'symmetric'): method for determining deformation
                'upper' produces an upper triangular defo
                'lower' produces a lower triangular defo
                'symmetric' produces a symmetric defo
        """
        return convert_strain_to_deformation(self, shape=shape)

    @property
    def von_mises_strain(self):
        """Equivalent strain to Von Mises Stress."""
        eps = self - 1 / 3 * np.trace(self) * np.identity(3)

        return np.sqrt(np.sum(eps * eps) * 2 / 3)


def convert_strain_to_deformation(strain, shape: Literal["upper", "lower", "symmetric"]):
    """
    This function converts a strain to a deformation gradient that will
    produce that strain. Supports three methods:

    Args:
        strain (3x3 array-like): strain matrix
        shape: ('upper' | 'lower' | 'symmetric'): method for determining deformation
            'upper' produces an upper triangular defo
            'lower' produces a lower triangular defo
            'symmetric' produces a symmetric defo
    """
    strain = SquareTensor(strain)
    ft_dot_f = 2 * strain + np.eye(3)
    if shape == "upper":
        result = scipy.linalg.cholesky(ft_dot_f)
    elif shape == "symmetric":
        result = scipy.linalg.sqrtm(ft_dot_f)
    else:
        raise ValueError('shape must be "upper" or "symmetric"')
    return Deformation(result)
