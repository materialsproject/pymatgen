# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes for non-standard space-group settings
"""

from __future__ import annotations

import re
from fractions import Fraction

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import MagSymmOp, SymmOp
from pymatgen.util.string import transformation_to_string

__author__ = "Matthew Horton"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matthew Horton"
__email__ = "mkhorton@lbl.gov"
__status__ = "Development"
__date__ = "Apr 2017"


class JonesFaithfulTransformation:
    """
    Transformation for space-groups defined in a non-standard setting
    """

    def __init__(self, P, p):
        """
        Transform between settings using matrix P and origin shift vector p,
        using same notation as reference.

        Should initialize using `from_transformation_string` in Jones
        faithful notation, given by a string specifying both a
        transformation matrix and an origin shift, with parts delimited
        by a semi-colon. Best shown by example:

        * `a,b,c;0,0,0` is the identity (no change)
        * `-b+c,a+c,-a+b+c;0,0,0` is R3:r to R3:h (rhombohedral to
          hexagonal setting)
        * `a,b,c;-1/4,-1/4,-1/4` is Pnnn:1 to Pnnn:2 (change in origin
          choice)
        * `b,c,a;-1/2,-1/2,-1/2` is Bbab:1 to Ccca:2 (change settin
          and origin)

        Can transform points (coords), lattices and symmetry operations.

        Used for transforming magnetic space groups since these are
        commonly used in multiple settings, due to needing to transform
        between magnetic and non-magnetic settings.

        See: International Tables for Crystallography (2016). Vol. A,
        Chapter 1.5, pp. 75-106.
        """
        # using capital letters in violation of PEP8 to
        # be consistent with variables in supplied reference,
        # for easier debugging in future
        self._P, self._p = P, p

    @classmethod
    def from_transformation_string(cls, transformation_string="a,b,c;0,0,0"):
        """
        Construct SpaceGroupTransformation from its transformation string.
        :param P: matrix
        :param p: origin shift vector
        :return:
        """
        P, p = JonesFaithfulTransformation.parse_transformation_string(transformation_string)
        return cls(P, p)

    @classmethod
    def from_origin_shift(cls, origin_shift="0,0,0"):
        """
        Construct SpaceGroupTransformation from its origin shift string.
        :param p: origin shift vector
        :return:
        """
        P = np.identity(3)
        p = [float(Fraction(x)) for x in origin_shift.split(",")]
        return cls(P, p)

    @staticmethod
    def parse_transformation_string(
        transformation_string: str = "a,b,c;0,0,0",
    ) -> tuple[list[list[float]] | np.ndarray, list[float]]:
        """
        Args:
            transformation_string (str, optional): Defaults to "a,b,c;0,0,0".

        Raises:
            ValueError: When transformation string fails to parse.

        Returns:
            tuple[list[list[float]] | np.ndarray, list[float]]: transformation matrix & vector
        """
        try:
            a = np.array([1, 0, 0])
            b = np.array([0, 1, 0])
            c = np.array([0, 0, 1])
            b_change, o_shift = transformation_string.split(";")
            basis_change = b_change.split(",")
            origin_shift = o_shift.split(",")
            # add implicit multiplication symbols
            basis_change = [
                re.sub(
                    r"(?<=\w|\))(?=\() | (?<=\))(?=\w) | (?<=(\d|a|b|c))(?=([abc]))",
                    r"*",
                    x,
                    flags=re.X,
                )
                for x in basis_change
            ]
            # should be fine to use eval here but be mindful for security
            # reasons
            # see http://lybniz2.sourceforge.net/safeeval.html
            # could replace with regex? or sympy expression?
            P = np.array([eval(x, {"__builtins__": None}, {"a": a, "b": b, "c": c}) for x in basis_change])
            P = P.transpose()  # by convention
            p = [float(Fraction(x)) for x in origin_shift]
            return P, p
        except Exception:
            raise ValueError("Failed to parse transformation string.")

    @property
    def P(self) -> list[list[float]]:
        """
        :return: transformation matrix
        """
        return self._P

    @property
    def p(self) -> list[float]:
        """
        :return: translation vector
        """
        return self._p

    @property
    def inverse(self) -> JonesFaithfulTransformation:
        """
        :return: JonesFaithfulTransformation
        """
        Q = np.linalg.inv(self.P)
        return JonesFaithfulTransformation(Q, -np.matmul(Q, self.p))

    @property
    def transformation_string(self) -> str:
        """
        :return: transformation string
        """
        return self._get_transformation_string_from_Pp(self.P, self.p)

    @staticmethod
    def _get_transformation_string_from_Pp(P: list[list[float]] | np.ndarray, p: list[float]) -> str:
        P = np.array(P).transpose()
        P_string = transformation_to_string(P, components=("a", "b", "c"))
        p_string = transformation_to_string(np.zeros((3, 3)), p)
        return P_string + ";" + p_string

    def transform_symmop(self, symmop: SymmOp | MagSymmOp) -> SymmOp | MagSymmOp:
        """
        Takes a symmetry operation and transforms it.
        :param symmop: SymmOp or MagSymmOp
        :return:
        """
        W_rot = symmop.rotation_matrix
        w_translation = symmop.translation_vector
        Q = np.linalg.inv(self.P)
        W_ = np.matmul(np.matmul(Q, W_rot), self.P)
        w_ = np.matmul(Q, (w_translation + np.matmul(W_rot - np.identity(3), self.p)))
        w_ = np.mod(w_, 1.0)
        if isinstance(symmop, MagSymmOp):
            return MagSymmOp.from_rotation_and_translation_and_time_reversal(
                rotation_matrix=W_,
                translation_vec=w_,
                time_reversal=symmop.time_reversal,
                tol=symmop.tol,
            )
        if isinstance(symmop, SymmOp):
            return SymmOp.from_rotation_and_translation(rotation_matrix=W_, translation_vec=w_, tol=symmop.tol)
        raise RuntimeError

    def transform_coords(self, coords: list[list[float]] | np.ndarray) -> list[list[float]]:
        """
        Takes a list of coordinates and transforms them.
        :param coords: List of coords
        :return:
        """
        new_coords = []
        for x in coords:
            x = np.array(x)
            Q = np.linalg.inv(self.P)
            x_ = np.matmul(Q, (x - self.p))
            new_coords.append(x_.tolist())
        return new_coords

    def transform_lattice(self, lattice: Lattice) -> Lattice:
        """
        Takes a Lattice object and transforms it.
        :param lattice: Lattice
        :return:
        """
        return Lattice(np.matmul(lattice.matrix, self.P))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return np.allclose(self.P, other.P) and np.allclose(self.p, other.p)

    def __str__(self):
        return str(JonesFaithfulTransformation.transformation_string)

    def __repr__(self):
        return f"JonesFaithfulTransformation with P:\n{self.P}\nand p:\n{self.p}"
