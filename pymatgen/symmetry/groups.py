# coding: utf-8

from __future__ import division, unicode_literals

"""
Defines SymmetryGroup parent class and PointGroup and SpaceGroup classes.
Shyue Ping Ong thanks Marc De Graef for his generous sharing of his
SpaceGroup data as published in his textbook "Structure of Materials".
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2013, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "4/4/14"

import os
from itertools import product
from fractions import Fraction
import numpy as np
import json

from pymatgen.core.operations import SymmOp


with open(os.path.join(os.path.dirname(__file__), "symm_data.json")) as f:
    SYMM_DATA = json.load(f)

GENERATOR_MATRICES = SYMM_DATA["generator_matrices"]
POINT_GROUP_ENC = SYMM_DATA["point_group_encoding"]
SPACE_GROUP_ENC = SYMM_DATA["space_group_encoding"]
TRANSLATIONS = {k: Fraction(v) for k, v in SYMM_DATA["translations"].items()}


class SymmetryGroup(object):
    """
    This is simply a parent class for a symmetry group to make it easier to
    share code that is common between any group.
    """

    def get_orbit(self, p, tol=1e-5):
        """
        Returns the orbit for a point.

        Args:
            p: Point as a 3x1 array.
            tol: Tolerance for determining if sites are the same. 1e-5 should
                be sufficient for most purposes. Set to 0 for exact matching
                (and also needed for symbolic orbits).

        Returns:
            ([array]) Orbit for point.
        """
        orbit = []
        for o in self.symmetry_ops:
            pp = o.operate(p)
            pp = np.mod(pp, 1)
            if not in_array_list(orbit, pp, tol=tol):
                orbit.append(pp)
        return orbit


class PointGroup(SymmetryGroup):
    """
    Class representing a Point Group, with generators and symmetry operations.

    .. attribute:: symbol

        Full International or Hermann-Mauguin Symbol.

    .. attribute:: generators

        List of generator matrices. Note that 3x3 matrices are used for Point
        Groups.

    .. attribute:: symmetry_ops

        Full set of symmetry operations as matrices.
    """

    def __init__(self, int_symbol):
        """
        Initializes a Point Group from its international symbol.

        Args:
            int_symbol (str): International or Hermann-Mauguin Symbol.
        """
        self.symbol = int_symbol
        self.generators = [GENERATOR_MATRICES[c]
                           for c in POINT_GROUP_ENC[int_symbol]]
        self.symmetry_ops = [SymmOp.from_rotation_and_translation(m)
                             for m in self._generate_full_symmetry_ops()]
        self.order = len(self.symmetry_ops)

    def _generate_full_symmetry_ops(self):
        symm_ops = list(self.generators)
        new_ops = self.generators
        while len(new_ops) > 0:
            gen_ops = []
            for g1, g2 in product(new_ops, symm_ops):
                op = np.dot(g1, g2)
                if not in_array_list(symm_ops, op):
                    gen_ops.append(op)
                    symm_ops.append(op)
            new_ops = gen_ops
        return symm_ops


class SpaceGroup(SymmetryGroup):
    """
    Class representing a SpaceGroup.

    .. attribute:: symbol

        Full International or Hermann-Mauguin Symbol.

    .. attribute:: int_number

        International number

    .. attribute:: generators

        List of generator matrices. Note that 4x4 matrices are used for Space
        Groups.

    .. attribute:: order

        Order of Space Group
    """

    #Contains the entire list of supported Space Group symbols.
    SG_SYMBOLS = tuple(SPACE_GROUP_ENC.keys())

    def __init__(self, int_symbol):
        """
        Initializes a Space Group from its *full* international symbol.

        Args:
            int_symbol (str): Full International or Hermann-Mauguin Symbol.
                The notation is a LaTeX-like string, with screw axes being
                represented by an underscore. For example, "P6_3/mmc". Note
                that for rhomohedral cells, the hexagonal setting can be
                accessed by adding a "H", e.g., "R-3mH".
        """
        self.symbol = int_symbol
        # TODO: Support different origin choices.
        enc = list(SPACE_GROUP_ENC[int_symbol]["enc"])
        inversion = int(enc.pop(0))
        ngen = int(enc.pop(0))
        symm_ops = [np.eye(4)]
        if inversion:
            symm_ops.append(np.array(
                [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0],
                 [0, 0, 0, 1]]))
        for i in range(ngen):
            m = np.eye(4)
            m[:3, :3] = GENERATOR_MATRICES[enc.pop(0)]
            m[0, 3] = TRANSLATIONS[enc.pop(0)]
            m[1, 3] = TRANSLATIONS[enc.pop(0)]
            m[2, 3] = TRANSLATIONS[enc.pop(0)]
            symm_ops.append(m)
        self.generators = symm_ops
        self.int_number = SPACE_GROUP_ENC[int_symbol]["int_number"]
        self.order = SPACE_GROUP_ENC[int_symbol]["order"]
        self._symmetry_ops = None

    def _generate_full_symmetry_ops(self):
        symm_ops = np.array(self.generators)
        for op in symm_ops:
            op[0:3, 3] = np.mod(op[0:3, 3], 1)
        new_ops = symm_ops
        while len(new_ops) > 0 and len(symm_ops) < self.order:
            gen_ops = []
            for g in new_ops:
                temp_ops = np.einsum('ijk,kl', symm_ops, g)
                for op in temp_ops:
                    op[0:3, 3] = np.mod(op[0:3, 3], 1)
                    ind = np.where(np.abs(1 - op[0:3, 3]) < 1e-5)
                    op[ind, 3] = 0
                    if not in_array_list(symm_ops, op):
                        gen_ops.append(op)
                        symm_ops = np.append(symm_ops, [op], axis=0)
            new_ops = gen_ops
        assert len(symm_ops) == self.order
        return symm_ops

    @property
    def symmetry_ops(self):
        """
        Full set of symmetry operations as matrices. Lazily initialized as
        generation sometimes takes a bit of time.
        """
        if self._symmetry_ops is None:
            self._symmetry_ops = [
                SymmOp(m) for m in self._generate_full_symmetry_ops()]
        return self._symmetry_ops

    @property
    def crystal_system(self):
        i = self.int_number
        if i <= 2:
            return "triclinic"
        elif i <= 15:
            return "monoclinic"
        elif i <= 74:
            return "orthorhombic"
        elif i <= 142:
            return "tetragonal"
        elif i <= 167:
            return "trigonal"
        elif i <= 194:
            return "hexagonal"
        else:
            return "cubic"

    @classmethod
    def from_int_number(cls, int_number, hexagonal=True):
        """
        Obtains a SpaceGroup from its international number.

        Args:
            int_number (int): International number.
            hexagonal (bool): For rhombohedral groups, whether to return the
                hexagonal setting (default) or rhombohedral setting.

        Returns:
            (SpaceGroup)
        """
        return SpaceGroup(sg_symbol_from_int_number(int_number,
                                                    hexagonal=hexagonal))

    def __str__(self):
        return "Spacegroup %s with international number %d and order %d" % (
            self.symbol, self.int_number, len(self.symmetry_ops))


def sg_symbol_from_int_number(int_number, hexagonal=True):
    """
    Obtains a SpaceGroup name from its international number.

    Args:
        int_number (int): International number.
        hexagonal (bool): For rhombohedral groups, whether to return the
            hexagonal setting (default) or rhombohedral setting.

    Returns:
        (str) Spacegroup symbol
    """
    syms = []
    for n, v in SPACE_GROUP_ENC.items():
        if v["int_number"] == int_number:
            syms.append(n)
    if len(syms) == 0:
        raise ValueError("Invalid international number!")
    if len(syms) == 2:
        if hexagonal:
            syms = list(filter(lambda s: s.endswith("H"), syms))
        else:
            syms = list(filter(lambda s: not s.endswith("H"), syms))
    return syms.pop()


def in_array_list(array_list, a, tol=1e-5):
    """
    Extremely efficient nd-array comparison using numpy's broadcasting. This
    function checks if a particular array a, is present in a list of arrays.
    It works for arrays of any size, e.g., even matrix searches.

    Args:
        array_list ([array]): A list of arrays to compare to.
        a (array): The test array for comparison.
        tol (float): The tolerance. Defaults to 1e-5. If 0, an exact match is
            done.

    Returns:
        (bool)
    """
    if len(array_list) == 0:
        return False
    axes = tuple(range(1, a.ndim + 1))
    if not tol:
        return np.any(np.all(np.equal(array_list, a[None, :]), axes))
    else:
        return np.any(np.sum(np.abs(array_list - a[None, :]), axes) < tol)
