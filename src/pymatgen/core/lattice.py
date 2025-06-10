"""This module defines the Lattice class, the fundamental class for representing
periodic crystals. It is essentially a matrix with some extra methods and attributes.
"""

from __future__ import annotations

import itertools
import math
import operator
import warnings
from collections import defaultdict
from fractions import Fraction
from functools import cached_property, reduce
from typing import TYPE_CHECKING

import numpy as np
from monty.dev import deprecated
from monty.json import MSONable

from pymatgen.util.coord import pbc_shortest_vectors
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import Literal

    from numpy.typing import ArrayLike, NDArray
    from typing_extensions import Self

    from pymatgen.core.operations import SymmOp

__author__ = "Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2011, The Materials Project"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"


class Lattice(MSONable):
    """Essentially a matrix with conversion matrices. In general,
    it is assumed that lengths are in Angstrom and angles are in
    degrees unless otherwise stated.

    Properties lazily generated for efficiency.
    """

    def __init__(
        self,
        matrix: ArrayLike,
        pbc: tuple[bool, bool, bool] = (True, True, True),
    ) -> None:
        """Create a lattice from any sequence of 9 numbers. Note that the sequence
        is assumed to be read one row at a time. Each row represents one
        lattice vector.

        Args:
            matrix: Sequence of numbers in any form. Examples of acceptable
                input.
                i) An actual numpy array.
                ii) [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                iii) [1, 0, 0 , 0, 1, 0, 0, 0, 1]
                iv) (1, 0, 0, 0, 1, 0, 0, 0, 1)
                Each row should correspond to a lattice vector.
                e.g. [[10, 0, 0], [20, 10, 0], [0, 0, 30]] specifies a lattice
                with lattice vectors [10, 0, 0], [20, 10, 0] and [0, 0, 30].
            pbc: a tuple defining the periodic boundary conditions along the three
                axis of the lattice.
        """
        mat = np.array(matrix, dtype=np.float64).reshape((3, 3))
        mat.setflags(write=False)
        self._matrix: NDArray[np.float64] = mat
        self._inv_matrix: NDArray[np.float64] | None = None
        self._diags = None
        self._lll_matrix_mappings: dict[float, tuple[NDArray[np.float64], NDArray[np.float64]]] = {}
        self._lll_inverse = None

        self.pbc = pbc

    def __repr__(self) -> str:
        return "\n".join(
            [
                "Lattice",
                f"    abc : {' '.join(map(repr, self.lengths))}",
                f" angles : {' '.join(map(repr, self.angles))}",
                f" volume : {self.volume!r}",
                f"      A : {' '.join(map(repr, self._matrix[0]))}",
                f"      B : {' '.join(map(repr, self._matrix[1]))}",
                f"      C : {' '.join(map(repr, self._matrix[2]))}",
                f"    pbc : {' '.join(map(repr, self.pbc))}",
            ]
        )

    def __eq__(self, other: object) -> bool:
        """A lattice is considered to be equal to another if the internal matrix
        representation satisfies np.allclose(matrix1, matrix2) and
        share the same periodicity.
        """
        if not hasattr(other, "matrix") or not hasattr(other, "pbc"):
            return NotImplemented

        # Shortcut the np.allclose if the memory addresses are the same
        # (very common in Structure.from_sites)
        if self is other:
            return True

        return np.allclose(self.matrix, other.matrix) and self.pbc == other.pbc

    def __hash__(self) -> int:
        return hash((self.lengths, self.angles, self.pbc))

    def __str__(self) -> str:
        return "\n".join(" ".join([f"{i:.6f}" for i in row]) for row in self._matrix)

    def __format__(self, fmt_spec: str = "") -> str:
        """Support format printing.

        Supported fmt_spec (str) are:
        1. "l" for a list format that can be easily copied and pasted, e.g.
           ".3fl" prints something like
           "[[10.000, 0.000, 0.000], [0.000, 10.000, 0.000], [0.000, 0.000, 10.000]]"
        2. "p" for lattice parameters ".1fp" prints something like
           "{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}"
        3. Default will simply print a 3x3 matrix form. E.g.
           10 0 0
           0 10 0
           0 0 10
        """
        matrix = self._matrix.tolist()
        if fmt_spec.endswith("l"):
            fmt = "[[{}, {}, {}], [{}, {}, {}], [{}, {}, {}]]"
            fmt_spec = fmt_spec[:-1]
        elif fmt_spec.endswith("p"):
            fmt = "{{{}, {}, {}, {}, {}, {}}}"
            fmt_spec = fmt_spec[:-1]
            matrix = (self.lengths, self.angles)
        else:
            fmt = "{} {} {}\n{} {} {}\n{} {} {}"

        return fmt.format(*(format(c, fmt_spec) for row in matrix for c in row))

    @property
    def pbc(self) -> tuple[bool, bool, bool]:
        """Tuple defining the periodicity of the Lattice."""
        return self._pbc  # type:ignore[return-value]

    @pbc.setter
    def pbc(self, pbc: tuple[bool, bool, bool]) -> None:
        if len(pbc) != 3 or any(item not in {True, False} for item in pbc):
            raise ValueError(f"pbc must be a tuple of three True/False values, got {pbc}")
        self._pbc = tuple(bool(item) for item in pbc)

    @property
    def matrix(self) -> NDArray[np.float64]:
        """Copy of matrix representing the Lattice."""
        return self._matrix

    @cached_property
    def lengths(self) -> tuple[float, float, float]:
        """Lattice lengths.

        Returns:
            The lengths (a, b, c) of the lattice.
        """
        return tuple(np.sqrt(np.sum(self._matrix**2, axis=1)).tolist())  # type: ignore[return-value]

    @cached_property
    def angles(self) -> tuple[float, float, float]:
        """Lattice angles.

        Returns:
            The angles (alpha, beta, gamma) of the lattice.
        """
        matrix, lengths = self._matrix, self.lengths
        angles = np.zeros(3)
        for dim in range(3):
            jj = (dim + 1) % 3
            kk = (dim + 2) % 3
            angles[dim] = np.clip(np.dot(matrix[jj], matrix[kk]) / (lengths[jj] * lengths[kk]), -1, 1)
        angles = np.arccos(angles) * 180.0 / np.pi  # type: ignore[assignment]
        return tuple(angles.tolist())  # type: ignore[return-value]

    @cached_property
    def volume(self) -> float:
        """Volume of the unit cell in Angstrom^3."""
        matrix = self._matrix
        return float(abs(np.dot(np.cross(matrix[0], matrix[1]), matrix[2])))

    @property
    def is_orthogonal(self) -> bool:
        """Whether all angles are 90 degrees."""
        return all(abs(a - 90) < 1e-5 for a in self.angles)

    @property
    def is_3d_periodic(self) -> bool:
        """True if the Lattice is periodic in all directions."""
        return all(self.pbc)

    @property
    def inv_matrix(self) -> NDArray[np.float64]:
        """Inverse of lattice matrix."""
        if self._inv_matrix is None:
            self._inv_matrix = np.linalg.inv(self._matrix)  # type: ignore[assignment]
            self._inv_matrix.setflags(write=False)
        return self._inv_matrix  # type: ignore[return-value]

    @property
    def metric_tensor(self) -> NDArray[np.float64]:
        """The metric tensor of the lattice."""
        return np.dot(self._matrix, self._matrix.T)

    def copy(self) -> Self:
        """Make a copy of this lattice."""
        return type(self)(self.matrix.copy(), pbc=self.pbc)

    def get_cartesian_coords(self, fractional_coords: ArrayLike) -> NDArray[np.float64]:
        """Get the Cartesian coordinates given fractional coordinates.

        Args:
            fractional_coords (3x1 array): Fractional coords.

        Returns:
            Cartesian coordinates
        """
        return np.dot(fractional_coords, self._matrix)

    def get_fractional_coords(self, cart_coords: ArrayLike) -> NDArray[np.float64]:
        """Get the fractional coordinates given Cartesian coordinates.

        Args:
            cart_coords (3x1 array): Cartesian coords.

        Returns:
            Fractional coordinates.
        """
        return np.dot(cart_coords, self.inv_matrix)

    def get_vector_along_lattice_directions(
        self,
        cart_coords: ArrayLike,
    ) -> NDArray[np.float64]:
        """Get the coordinates along lattice directions given Cartesian coordinates.

        Note, this is different than a projection of the Cartesian vector along the
        lattice parameters. It is simply the fractional coordinates multiplied by the
        lattice vector magnitudes.

        For example, this method is helpful when analyzing the dipole moment (in
        units of electron Angstroms) of a ferroelectric crystal. See the `Polarization`
        class in `pymatgen.analysis.ferroelectricity.polarization`.

        Args:
            cart_coords (3x1 array): Cartesian coords.

        Returns:
            Lattice coordinates.
        """
        return self.get_fractional_coords(cart_coords) * self.lengths

    def d_hkl(self, miller_index: tuple[int, ...]) -> float:
        """Get the distance between the hkl plane and the origin.

        Args:
            miller_index (tuple[int, ...]): Miller index of plane

        Returns:
            float: distance between hkl plane and origin
        """
        g_star = self.reciprocal_lattice_crystallographic.metric_tensor
        hkl = np.array(miller_index)
        return 1 / ((np.dot(np.dot(hkl, g_star), hkl.T)) ** (1 / 2))

    @classmethod
    def cubic(cls, a: float, pbc: tuple[bool, bool, bool] = (True, True, True)) -> Self:
        """Convenience constructor for a cubic lattice.

        Args:
            a (float): The *a* lattice parameter of the cubic cell.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Cubic lattice of dimensions (a x a x a).
        """
        return cls([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]], pbc)

    @classmethod
    def tetragonal(
        cls,
        a: float,
        c: float,
        pbc: tuple[bool, bool, bool] = (True, True, True),
    ) -> Self:
        """Convenience constructor for a tetragonal lattice.

        Args:
            a (float): *a* lattice parameter of the tetragonal cell.
            c (float): *c* lattice parameter of the tetragonal cell.
            pbc (tuple): The periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Tetragonal lattice of dimensions (a x a x c).
        """
        return cls.from_parameters(a, a, c, 90, 90, 90, pbc=pbc)

    @classmethod
    def orthorhombic(
        cls,
        a: float,
        b: float,
        c: float,
        pbc: tuple[bool, bool, bool] = (True, True, True),
    ) -> Self:
        """Convenience constructor for an orthorhombic lattice.

        Args:
            a (float): *a* lattice parameter of the orthorhombic cell.
            b (float): *b* lattice parameter of the orthorhombic cell.
            c (float): *c* lattice parameter of the orthorhombic cell.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Orthorhombic lattice of dimensions (a x b x c).
        """
        return cls.from_parameters(a, b, c, 90, 90, 90, pbc=pbc)

    @classmethod
    def monoclinic(
        cls,
        a: float,
        b: float,
        c: float,
        beta: float,
        pbc: tuple[bool, bool, bool] = (True, True, True),
    ) -> Self:
        """Convenience constructor for a monoclinic lattice.

        Args:
            a (float): *a* lattice parameter of the monoclinic cell.
            b (float): *b* lattice parameter of the monoclinic cell.
            c (float): *c* lattice parameter of the monoclinic cell.
            beta (float): *beta* angle between lattice vectors b and c in
                degrees.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Monoclinic lattice of dimensions (a x b x c) with non
                right-angle beta between lattice vectors a and c.
        """
        return cls.from_parameters(a, b, c, 90, beta, 90, pbc=pbc)

    @classmethod
    def hexagonal(cls, a: float, c: float, pbc: tuple[bool, bool, bool] = (True, True, True)) -> Self:
        """Convenience constructor for a hexagonal lattice.

        Args:
            a (float): *a* lattice parameter of the hexagonal cell.
            c (float): *c* lattice parameter of the hexagonal cell.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Hexagonal lattice of dimensions (a x a x c).
        """
        return cls.from_parameters(a, a, c, 90, 90, 120, pbc=pbc)

    @classmethod
    def rhombohedral(
        cls,
        a: float,
        alpha: float,
        pbc: tuple[bool, bool, bool] = (True, True, True),
    ) -> Self:
        """Convenience constructor for a rhombohedral lattice.

        Args:
            a (float): *a* lattice parameter of the rhombohedral cell.
            alpha (float): Angle for the rhombohedral lattice in degrees.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Rhombohedral lattice of dimensions (a x a x a).
        """
        return cls.from_parameters(a, a, a, alpha, alpha, alpha, pbc=pbc)

    @classmethod
    def from_parameters(
        cls,
        a: float,
        b: float,
        c: float,
        alpha: float,
        beta: float,
        gamma: float,
        *,
        vesta: bool = False,
        pbc: tuple[bool, bool, bool] = (True, True, True),
    ) -> Self:
        """Create a Lattice using unit cell lengths (in Angstrom) and angles (in degrees).

        Args:
            a (float): *a* lattice parameter.
            b (float): *b* lattice parameter.
            c (float): *c* lattice parameter.
            alpha (float): *alpha* angle in degrees.
            beta (float): *beta* angle in degrees.
            gamma (float): *gamma* angle in degrees.
            vesta (bool): True if you import Cartesian coordinates from VESTA.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Lattice with the specified lattice parameters.
        """
        angles_r = np.radians([alpha, beta, gamma])
        cos_alpha, cos_beta, cos_gamma = np.cos(angles_r)
        sin_alpha, sin_beta, sin_gamma = np.sin(angles_r)

        if vesta:
            c1 = c * cos_beta
            c2 = (c * (cos_alpha - (cos_beta * cos_gamma))) / sin_gamma

            vector_a = [float(a), 0.0, 0.0]
            vector_b = [b * cos_gamma, b * sin_gamma, 0]
            vector_c = [c1, c2, math.sqrt(c**2 - c1**2 - c2**2)]

        else:
            val = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta)
            val = np.clip(val, -1, 1)  # rounding errors may cause values slightly > 1
            gamma_star = np.arccos(val)

            vector_a = [a * sin_beta, 0.0, a * cos_beta]
            vector_b = [
                -b * sin_alpha * np.cos(gamma_star),
                b * sin_alpha * np.sin(gamma_star),
                b * cos_alpha,
            ]
            vector_c = [0.0, 0.0, float(c)]

        return cls([vector_a, vector_b, vector_c], pbc)

    @classmethod
    def from_dict(
        cls,
        dct: dict,
        fmt: str | None = None,
        **kwargs,
    ) -> Self:
        """Create a Lattice from a dictionary.

        If fmt is None, the dict should contain the a, b, c,
        alpha, beta, and gamma parameters.

        If fmt == "abivars", the function build a `Lattice` object from a
        dictionary with the Abinit variables `acell` and `rprim` in Bohr.
        If acell is not given, the Abinit default is used i.e. [1,1,1] Bohr

        Example:
            Lattice.from_dict(fmt="abivars", acell=3*[10], rprim=np.eye(3))
        """
        if fmt == "abivars":
            from pymatgen.io.abinit.abiobjects import lattice_from_abivars

            kwargs |= dct
            return lattice_from_abivars(cls=cls, **kwargs)

        pbc = dct.get("pbc", (True, True, True))
        if "matrix" in dct:
            return cls(dct["matrix"], pbc=pbc)
        return cls.from_parameters(
            dct["a"],
            dct["b"],
            dct["c"],
            dct["alpha"],
            dct["beta"],
            dct["gamma"],
            pbc=pbc,
        )

    @property
    def a(self) -> float:
        """*a* lattice parameter."""
        return self.lengths[0]

    @property
    def b(self) -> float:
        """*b* lattice parameter."""
        return self.lengths[1]

    @property
    def c(self) -> float:
        """*c* lattice parameter."""
        return self.lengths[2]

    @property
    def abc(self) -> tuple[float, float, float]:
        """Lengths of the lattice vectors, i.e. (a, b, c)."""
        return self.lengths

    @property
    def alpha(self) -> float:
        """Angle alpha of lattice in degrees."""
        return self.angles[0]

    @property
    def beta(self) -> float:
        """Angle beta of lattice in degrees."""
        return self.angles[1]

    @property
    def gamma(self) -> float:
        """Angle gamma of lattice in degrees."""
        return self.angles[2]

    @property
    def parameters(self) -> tuple[float, float, float, float, float, float]:
        """6-tuple of floats (a, b, c, alpha, beta, gamma)."""
        return (*self.lengths, *self.angles)

    @property
    def params_dict(self) -> dict[str, float]:
        """Dictionary of lattice parameters."""
        return dict(zip(("a", "b", "c", "alpha", "beta", "gamma"), self.parameters, strict=True))

    @property
    def reciprocal_lattice(self) -> Self:
        """The reciprocal lattice. Note that this is the standard
        reciprocal lattice used for solid state physics with a factor of 2 *
        pi. If you are looking for the crystallographic reciprocal lattice,
        use the reciprocal_lattice_crystallographic property.
        The property is lazily generated for efficiency.
        """
        inv_mat = np.linalg.inv(self._matrix).T
        return type(self)(inv_mat * 2 * np.pi)

    @property
    def reciprocal_lattice_crystallographic(self) -> Self:
        """The *crystallographic* reciprocal lattice, i.e. no factor of 2 * pi."""
        return type(self)(self.reciprocal_lattice.matrix / (2 * np.pi))

    @property
    def lll_matrix(self) -> NDArray[np.float64]:
        """The matrix for LLL reduction."""
        if 0.75 not in self._lll_matrix_mappings:
            self._lll_matrix_mappings[0.75] = self._calculate_lll()
        return self._lll_matrix_mappings[0.75][0]

    @property
    def lll_mapping(self) -> NDArray[np.float64]:
        """The mapping between the LLL reduced lattice and the original lattice."""
        if 0.75 not in self._lll_matrix_mappings:
            self._lll_matrix_mappings[0.75] = self._calculate_lll()
        return self._lll_matrix_mappings[0.75][1]

    @property
    def lll_inverse(self) -> NDArray[np.float64]:
        """Inverse of self.lll_mapping."""
        return np.linalg.inv(self.lll_mapping)  # type: ignore[return-value]

    @property
    def selling_vector(self) -> NDArray[np.float64]:
        """The (1,6) array of Selling Scalars."""
        a, b, c = self.matrix
        d = -(a + b + c)
        tol = 1e-10

        selling_vector = np.array(
            [
                np.dot(b, c),
                np.dot(a, c),
                np.dot(a, b),
                np.dot(a, d),
                np.dot(b, d),
                np.dot(c, d),
            ]
        )
        selling_vector = np.array([s if abs(s) > tol else 0 for s in selling_vector])

        reduction_matrices = [
            [
                [-1, 0, 0, 0, 0, 0],
                [1, 1, 0, 0, 0, 0],
                [1, 0, 0, 0, 1, 0],
                [-1, 0, 0, 1, 0, 0],
                [1, 0, 1, 0, 0, 0],
                [1, 0, 0, 0, 0, 1],
            ],
            [
                [1, 1, 0, 0, 0, 0],
                [0, -1, 0, 0, 0, 0],
                [0, 1, 0, 1, 0, 0],
                [0, 1, 1, 0, 0, 0],
                [0, -1, 0, 0, 1, 0],
                [0, 1, 0, 0, 0, 1],
            ],
            [
                [1, 0, 1, 0, 0, 0],
                [0, 0, 1, 1, 0, 0],
                [0, 0, -1, 0, 0, 0],
                [0, 1, 1, 0, 0, 0],
                [0, 0, 1, 0, 1, 0],
                [0, 0, -1, 0, 0, 1],
            ],
            [
                [1, 0, 0, -1, 0, 0],
                [0, 0, 1, 1, 0, 0],
                [0, 1, 0, 1, 0, 0],
                [0, 0, 0, -1, 0, 0],
                [0, 0, 0, 1, 1, 0],
                [0, 0, 0, 1, 0, 1],
            ],
            [
                [0, 0, 1, 0, 1, 0],
                [0, 1, 0, 0, -1, 0],
                [1, 0, 0, 0, 1, 0],
                [0, 0, 0, 1, 1, 0],
                [0, 0, 0, 0, -1, 0],
                [0, 0, 0, 0, 1, 1],
            ],
            [
                [0, 1, 0, 0, 0, 1],
                [1, 0, 0, 0, 0, 1],
                [0, 0, 1, 0, 0, -1],
                [0, 0, 0, 1, 0, 1],
                [0, 0, 0, 0, 1, 1],
                [0, 0, 0, 0, 0, -1],
            ],
        ]

        while np.greater(np.max(selling_vector), 0):
            max_index = selling_vector.argmax()
            selling_vector = np.dot(reduction_matrices[max_index], selling_vector)

        return selling_vector

    def selling_dist(self, other: Self) -> float:
        """Get the minimum Selling distance between two lattices."""
        vcp_matrices = [
            [
                [-1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            [
                [1, 0, 0, 0, 0, 0],
                [0, -1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            [
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, -1, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            [
                [1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, -1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            [
                [0, 0, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, -1, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            [
                [0, 1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, -1],
            ],
        ]

        reflection_matrices = [
            [
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            [
                [1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 1, 0],
            ],
            [
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 1, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
            ],
            [
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
            ],
            [
                [0, 1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            [
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 1, 0, 0],
            ],
            [
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 1, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
            ],
            [
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 1, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
            ],
            [
                [0, 0, 1, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
            ],
            [
                [0, 0, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 1, 0, 0],
            ],
            [
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
            ],
            [
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
            ],
            [
                [0, 0, 0, 1, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 1, 0, 0, 0],
            ],
            [
                [0, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 1, 0, 0, 0, 0],
            ],
            [
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 1, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            [
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
            ],
            [
                [0, 0, 0, 0, 1, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 0],
            ],
            [
                [0, 0, 0, 0, 1, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
                [1, 0, 0, 0, 0, 0],
            ],
            [
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            [
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
            ],
            [
                [0, 0, 0, 0, 0, 1],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 1, 0, 0, 0, 0],
            ],
            [
                [0, 0, 0, 0, 0, 1],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [1, 0, 0, 0, 0, 0],
            ],
            [
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 1, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0],
            ],
            [
                [0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 1, 0],
                [1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
            ],
        ]

        selling1 = self.selling_vector
        selling2 = other.selling_vector

        vcps = np.dot(selling1, vcp_matrices)[0]

        all_reflections = []
        for vcp in vcps:
            for reflection_matrix in reflection_matrices:
                all_reflections.append(np.dot(vcp, reflection_matrix))

        for reflection_matrix in reflection_matrices:
            all_reflections.append(np.dot(selling1, reflection_matrix))

        return min(np.linalg.norm(reflection - selling2) for reflection in all_reflections)  # type: ignore[return-value, type-var]

    def as_dict(self, verbosity: Literal[0, 1] = 0) -> dict:
        """MSONable dict representation of the Lattice.

        Args:
            verbosity (0 | 1): Default of 0 only includes the matrix representation.
                Set to 1 to include the lattice parameters.
        """
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "matrix": self._matrix.tolist(),
            "pbc": self.pbc,
        }

        if verbosity not in {0, 1}:
            warnings.warn(
                f"`verbosity={verbosity}` is deprecated and will be disallowed in a future version. "
                "Please use 0 (silent) or 1 (verbose) explicitly.",
                DeprecationWarning,
                stacklevel=2,
            )
        if verbosity > 0:  # TODO: explicit check `verbosity == 1`
            dct |= self.params_dict
            dct["volume"] = self.volume

        return dct

    def find_all_mappings(
        self,
        other_lattice: Self,
        ltol: float = 1e-5,
        atol: float = 1,
        skip_rotation_matrix: bool = False,
    ) -> Iterator[tuple[Lattice, NDArray[np.float64] | None, NDArray[np.float64]]]:
        """Find all mappings between current lattice and another lattice.

        Args:
            other_lattice (Lattice): Another lattice that is equivalent to this one.
            ltol (float): Tolerance for matching lengths. Defaults to 1e-5.
            atol (float): Tolerance for matching angles. Defaults to 1.
            skip_rotation_matrix (bool): Whether to skip calculation of the
                rotation matrix

        Yields:
            (aligned_lattice, rotation_matrix, scale_matrix) if a mapping is
            found. aligned_lattice is a rotated version of other_lattice that
            has the same lattice parameters, but which is aligned in the
            coordinate system of this lattice so that translational points
            match up in 3D. rotation_matrix is the rotation that has to be
            applied to other_lattice to obtain aligned_lattice, i.e.,
            aligned_matrix = np.inner(other_lattice, rotation_matrix) and
            op = SymmOp.from_rotation_and_translation(rotation_matrix)
            aligned_matrix = op.operate_multi(latt.matrix)
            Finally, scale_matrix is the integer matrix that expresses
            aligned_matrix as a linear combination of this
            lattice, i.e., aligned_matrix = np.dot(scale_matrix, self.matrix)

            None is returned if no matches are found.
        """

        def get_angles(v1, v2, l1, l2):
            x = np.inner(v1, v2) / l1[:, None] / l2
            x[x > 1] = 1
            x[x < -1] = -1
            return np.arccos(x) * 180.0 / np.pi

        lengths = other_lattice.lengths
        alpha, beta, gamma = other_lattice.angles

        frac, dist, _, _ = self.get_points_in_sphere(  # type: ignore[misc]
            [[0, 0, 0]], [0, 0, 0], max(lengths) * (1 + ltol), zip_results=False
        )
        cart = self.get_cartesian_coords(frac)  # type: ignore[arg-type]
        # This can't be broadcast because they're different lengths
        inds = [np.logical_and(dist / ln < 1 + ltol, dist / ln > 1 / (1 + ltol)) for ln in lengths]  # type: ignore[operator]
        c_a, c_b, c_c = (cart[i] for i in inds)
        f_a, f_b, f_c = (frac[i] for i in inds)
        l_a, l_b, l_c = (np.sum(c**2, axis=-1) ** 0.5 for c in (c_a, c_b, c_c))

        alpha_b = np.isclose(get_angles(c_b, c_c, l_b, l_c), alpha, atol=atol, rtol=0)
        beta_b = np.isclose(get_angles(c_a, c_c, l_a, l_c), beta, atol=atol, rtol=0)
        gamma_b = np.isclose(get_angles(c_a, c_b, l_a, l_b), gamma, atol=atol, rtol=0)

        for idx, all_j in enumerate(gamma_b):
            inds = np.logical_and(all_j[:, None], np.logical_and(alpha_b, beta_b[idx][None, :]))
            for j, k in np.argwhere(inds):
                scale_m = np.array((f_a[idx], f_b[j], f_c[k]), dtype=np.int64)  # type: ignore[index]
                if abs(np.linalg.det(scale_m)) < 1e-8:
                    continue

                aligned_m = np.array((c_a[idx], c_b[j], c_c[k]))

                rotation_m = None if skip_rotation_matrix else np.linalg.solve(aligned_m, other_lattice.matrix)

                yield type(self)(aligned_m), rotation_m, scale_m

    def find_mapping(
        self,
        other_lattice: Self,
        ltol: float = 1e-5,
        atol: float = 1,
        skip_rotation_matrix: bool = False,
    ) -> tuple[Lattice, NDArray[np.float64] | None, NDArray[np.float64]] | None:
        """Find a mapping between current lattice and another lattice. There
        are an infinite number of choices of basis vectors for two entirely
        equivalent lattices. This method returns a mapping that maps
        other_lattice to this lattice.

        Args:
            other_lattice (Lattice): Another lattice that is equivalent to
                this one.
            ltol (float): Tolerance for matching lengths. Defaults to 1e-5.
            atol (float): Tolerance for matching angles. Defaults to 1.
            skip_rotation_matrix (bool): Whether to skip calculation of the rotation matrix.
                Defaults to False.

        Returns:
            tuple[Lattice, NDArray[np.float_], NDArray[np.float_]]: (aligned_lattice, rotation_matrix, scale_matrix)
            if a mapping is found. aligned_lattice is a rotated version of other_lattice that
            has the same lattice parameters, but which is aligned in the
            coordinate system of this lattice so that translational points
            match up in 3D. rotation_matrix is the rotation that has to be
            applied to other_lattice to obtain aligned_lattice, i.e.,
            aligned_matrix = np.inner(other_lattice, rotation_matrix) and
            op = SymmOp.from_rotation_and_translation(rotation_matrix)
            aligned_matrix = op.operate_multi(latt.matrix)
            Finally, scale_matrix is the integer matrix that expresses
            aligned_matrix as a linear combination of this
            lattice, i.e., aligned_matrix = np.dot(scale_matrix, self.matrix)

            None is returned if no matches are found.
        """
        return next(
            self.find_all_mappings(other_lattice, ltol, atol, skip_rotation_matrix),
            None,
        )

    def get_lll_reduced_lattice(self, delta: float = 0.75) -> Self:
        """Lenstra-Lenstra-Lovasz lattice basis reduction.

        Args:
            delta: Delta parameter.

        Returns:
            Lattice: LLL reduced
        """
        if delta not in self._lll_matrix_mappings:
            self._lll_matrix_mappings[delta] = self._calculate_lll()
        return type(self)(self._lll_matrix_mappings[delta][0])

    def _calculate_lll(self, delta: float = 0.75) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Perform a Lenstra-Lenstra-Lovasz lattice basis reduction to obtain a
        c-reduced basis. This method returns a basis which is as "good" as
        possible, with "good" defined by orthogonality of the lattice vectors.

        This basis is used for all the periodic boundary condition calculations.

        Args:
            delta (float): Reduction parameter. Default of 0.75 is usually fine.

        Returns:
            Reduced lattice matrix, mapping to get to that lattice.
        """
        # Transpose the lattice matrix first so that basis vectors are columns.
        # Makes life easier.

        a = self._matrix.copy().T

        b = np.zeros((3, 3))  # Vectors after the Gram-Schmidt process
        u = np.zeros((3, 3))  # Gram-Schmidt coefficients
        m = np.zeros(3)  # These are the norm squared of each vec

        b[:, 0] = a[:, 0]
        m[0] = np.dot(b[:, 0], b[:, 0])
        for i in range(1, 3):
            u[i, :i] = np.dot(a[:, i].T, b[:, :i]) / m[:i]
            b[:, i] = a[:, i] - np.dot(b[:, :i], u[i, :i].T)
            m[i] = np.dot(b[:, i], b[:, i])

        k = 2

        mapping = np.identity(3, dtype=np.double)
        while k <= 3:
            # Size reduction
            for i in range(k - 1, 0, -1):
                q = round(u[k - 1, i - 1])
                if q != 0:
                    # Reduce the k-th basis vector
                    a[:, k - 1] -= q * a[:, i - 1]
                    mapping[:, k - 1] -= q * mapping[:, i - 1]
                    uu = list(u[i - 1, 0 : (i - 1)])
                    uu.append(1)  # type:ignore[arg-type]
                    # Update the GS coefficients
                    u[k - 1, 0:i] -= q * np.array(uu)

            # Check the Lovasz condition
            if np.dot(b[:, k - 1], b[:, k - 1]) >= (delta - abs(u[k - 1, k - 2]) ** 2) * np.dot(
                b[:, (k - 2)], b[:, (k - 2)]
            ):
                # Increment k if the Lovasz condition holds
                k += 1
            else:
                # If the Lovasz condition fails, swap the k-th and (k-1)-th basis vector
                v = a[:, k - 1].copy()
                a[:, k - 1] = a[:, k - 2].copy()
                a[:, k - 2] = v

                v_m = mapping[:, k - 1].copy()
                mapping[:, k - 1] = mapping[:, k - 2].copy()
                mapping[:, k - 2] = v_m

                # Update the Gram-Schmidt coefficients
                for s in range(k - 1, k + 1):
                    u[s - 1, : (s - 1)] = np.dot(a[:, s - 1].T, b[:, : (s - 1)]) / m[: (s - 1)]
                    b[:, s - 1] = a[:, s - 1] - np.dot(b[:, : (s - 1)], u[s - 1, : (s - 1)].T)
                    m[s - 1] = np.dot(b[:, s - 1], b[:, s - 1])

                if k > 2:
                    k -= 1
                else:
                    # We have to do p/q, so do lstsq(q.T, p.T).T instead
                    p = np.dot(a[:, k:3].T, b[:, (k - 2) : k])
                    q = np.diag(m[(k - 2) : k])

                    result = np.linalg.lstsq(q.T, p.T, rcond=None)[0].T
                    u[k:3, (k - 2) : k] = result

        return a.T, mapping.T  # type: ignore[return-value]

    def get_lll_frac_coords(self, frac_coords: ArrayLike) -> NDArray[np.float64]:
        """Given fractional coordinates in the lattice basis, returns corresponding
        fractional coordinates in the lll basis.
        """
        return np.dot(frac_coords, self.lll_inverse)

    def get_frac_coords_from_lll(self, lll_frac_coords: ArrayLike) -> NDArray[np.float64]:
        """Given fractional coordinates in the lll basis, returns corresponding
        fractional coordinates in the lattice basis.
        """
        return np.dot(lll_frac_coords, self.lll_mapping)

    @due.dcite(
        Doi("10.1107/S010876730302186X"),
        description="Numerically stable algorithms for the computation of reduced unit cells",
    )
    def get_niggli_reduced_lattice(self, tol: float = 1e-5) -> Self:
        """Get the Niggli reduced lattice using the numerically stable algo
        proposed by R. W. Grosse-Kunstleve, N. K. Sauter, & P. D. Adams,
        Acta Crystallographica Section A Foundations of Crystallography, 2003,
        60(1), 1-6. doi:10.1107/S010876730302186X.

        Args:
            tol (float): The numerical tolerance. The default of 1e-5 should
                result in stable behavior for most cases.

        Returns:
            Lattice: Niggli-reduced lattice.
        """
        # lll reduction is more stable for skewed cells
        matrix = self.lll_matrix
        e = tol * self.volume ** (1 / 3)

        # Define metric tensor
        G = np.dot(matrix, matrix.T)

        # This sets an upper limit on the number of iterations.
        for _ in range(100):
            # The steps are labelled as Ax as per the labelling scheme in the
            # paper.
            A, B, C, E, N, Y = (
                G[0, 0],
                G[1, 1],
                G[2, 2],
                2 * G[1, 2],
                2 * G[0, 2],
                2 * G[0, 1],
            )

            if B + e < A or (abs(A - B) < e and abs(E) > abs(N) + e):
                # A1
                M = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, -1]])
                G = np.dot(np.transpose(M), np.dot(G, M))
                # update lattice parameters based on new G (gh-3657)
                A, B, C, E, N, Y = (
                    G[0, 0],
                    G[1, 1],
                    G[2, 2],
                    2 * G[1, 2],
                    2 * G[0, 2],
                    2 * G[0, 1],
                )

            if (C + e < B) or (abs(B - C) < e and abs(N) > abs(Y) + e):
                # A2
                M = np.array([[-1, 0, 0], [0, 0, -1], [0, -1, 0]])
                G = np.dot(np.transpose(M), np.dot(G, M))
                continue

            ll = 0 if abs(E) < e else E / abs(E)
            m = 0 if abs(N) < e else N / abs(N)
            n = 0 if abs(Y) < e else Y / abs(Y)
            if ll * m * n == 1:
                # A3
                i = -1 if ll == -1 else 1
                j = -1 if m == -1 else 1
                k = -1 if n == -1 else 1
                M = np.diag((i, j, k))
                G = np.dot(np.transpose(M), np.dot(G, M))
            elif ll * m * n in (0, -1):
                # A4
                i = -1 if ll == 1 else 1
                j = -1 if m == 1 else 1
                k = -1 if n == 1 else 1

                if i * j * k == -1:
                    if n == 0:
                        k = -1
                    elif m == 0:
                        j = -1
                    elif ll == 0:
                        i = -1
                M = np.diag((i, j, k))
                G = np.dot(np.transpose(M), np.dot(G, M))

            A, B, C, E, N, Y = (
                G[0, 0],
                G[1, 1],
                G[2, 2],
                2 * G[1, 2],
                2 * G[0, 2],
                2 * G[0, 1],
            )

            # A5
            if abs(E) > B + e or (abs(E - B) < e and Y - e > 2 * N) or (abs(E + B) < e and -e > Y):
                M = np.array([[1, 0, 0], [0, 1, -E / abs(E)], [0, 0, 1]])
                G = np.dot(np.transpose(M), np.dot(G, M))
                continue

            # A6
            if abs(N) > A + e or (abs(A - N) < e and Y - e > 2 * E) or (abs(A + N) < e and -e > Y):
                M = np.array([[1, 0, -N / abs(N)], [0, 1, 0], [0, 0, 1]])
                G = np.dot(np.transpose(M), np.dot(G, M))
                continue

            # A7
            if abs(Y) > A + e or (abs(A - Y) < e and N - e > 2 * E) or (abs(A + Y) < e and -e > N):
                M = np.array([[1, -Y / abs(Y), 0], [0, 1, 0], [0, 0, 1]])
                G = np.dot(np.transpose(M), np.dot(G, M))
                continue

            # A8
            if -e > E + N + Y + A + B or (abs(E + N + Y + A + B) < e < Y + (A + N) * 2):
                M = np.array([[1, 0, 1], [0, 1, 1], [0, 0, 1]])
                G = np.dot(np.transpose(M), np.dot(G, M))
                continue

            break

        A = G[0, 0]
        B = G[1, 1]
        C = G[2, 2]
        E = 2 * G[1, 2]
        N = 2 * G[0, 2]
        Y = 2 * G[0, 1]
        a = math.sqrt(A)
        b = math.sqrt(B)
        c = math.sqrt(C)
        alpha = math.acos(E / 2 / b / c) / math.pi * 180
        beta = math.acos(N / 2 / a / c) / math.pi * 180
        gamma = math.acos(Y / 2 / a / b) / math.pi * 180
        lattice = type(self).from_parameters(a, b, c, alpha, beta, gamma)

        mapped = self.find_mapping(lattice, e, skip_rotation_matrix=True)
        if mapped is not None:
            if np.linalg.det(mapped[0].matrix) > 0:
                return mapped[0]  # type:ignore[return-value]
            return type(self)(-mapped[0].matrix)

        raise ValueError("can't find niggli")

    def scale(self, new_volume: float) -> Self:
        """Return a new Lattice with volume new_volume by performing a
        scaling of the lattice vectors so that length proportions and angles
        are preserved.

        Args:
            new_volume:
                New volume to scale to.

        Returns:
            New lattice with desired volume.
        """
        versors = self.matrix / self.lengths

        geo_factor = abs(np.dot(np.cross(versors[0], versors[1]), versors[2]))

        ratios = np.array(self.lengths) / self.c

        new_c = (new_volume / (geo_factor * np.prod(ratios))) ** (1 / 3.0)

        return type(self)(versors * (new_c * ratios), pbc=self.pbc)

    def get_wigner_seitz_cell(self) -> list[list[NDArray[np.float64]]]:
        """Get the Wigner-Seitz cell for the given lattice.

        Returns:
            A list of list of coordinates.
            Each element in the list is a "facet" of the boundary of the
            Wigner Seitz cell. For instance, a list of four coordinates will
            represent a square facet.
        """
        vec1, vec2, vec3 = self._matrix

        list_k_points = []
        for ii, jj, kk in itertools.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
            list_k_points.append(ii * vec1 + jj * vec2 + kk * vec3)
        from scipy.spatial import Voronoi

        tess = Voronoi(list_k_points)
        out = []
        for r in tess.ridge_dict:
            if r[0] == 13 or r[1] == 13:
                out.append([tess.vertices[i] for i in tess.ridge_dict[r]])

        return out

    def get_brillouin_zone(self) -> list[list[NDArray[np.float64]]]:
        """Get the Wigner-Seitz cell for the reciprocal lattice, aka the
        Brillouin Zone.

        Returns:
            A list of list of coordinates.
            Each element in the list is a "facet" of the boundary of the
            Brillouin Zone. For instance, a list of four coordinates will
            represent a square facet.
        """
        return self.reciprocal_lattice.get_wigner_seitz_cell()

    def dot(
        self,
        coords_a: ArrayLike,
        coords_b: ArrayLike,
        frac_coords: bool = False,
    ) -> NDArray[np.float64]:
        """Compute the scalar product of vector(s).

        Args:
            coords_a: Array-like coordinates.
            coords_b: Array-like coordinates.
            frac_coords (bool): True if the vectors are fractional (as opposed to Cartesian) coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        coords_a, coords_b = np.reshape(coords_a, (-1, 3)), np.reshape(coords_b, (-1, 3))

        if len(coords_a) != len(coords_b):
            raise ValueError("Coordinates must have same length!")

        for coord in (coords_a, coords_b):
            if np.iscomplexobj(coord):
                raise TypeError(f"Complex array are not supported, got {coord=}")

        if not frac_coords:
            cart_a, cart_b = coords_a, coords_b
        else:
            cart_a = np.reshape([self.get_cartesian_coords(vec) for vec in coords_a], (-1, 3))
            cart_b = np.reshape([self.get_cartesian_coords(vec) for vec in coords_b], (-1, 3))

        return np.array(list(itertools.starmap(np.dot, zip(cart_a, cart_b, strict=True))))

    def norm(self, coords: ArrayLike, frac_coords: bool = True) -> NDArray[np.float64]:
        """Compute the norm of vector(s).

        Args:
            coords:
                Array-like object with the coordinates.
            frac_coords:
                Boolean stating whether the vector corresponds to fractional or
                Cartesian coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        return np.sqrt(self.dot(coords, coords, frac_coords=frac_coords))

    def get_points_in_sphere(
        self,
        frac_points: ArrayLike,
        center: ArrayLike,
        r: float,
        zip_results: bool = True,
    ) -> list[tuple[NDArray, float, int, NDArray]] | tuple[NDArray] | NDArray:
        """Find all points within a sphere from the point taking into account
        periodic boundary conditions. This includes sites in other periodic images.

        Algorithm:

        1. place sphere of radius r in crystal and determine minimum supercell
           (parallelepiped) which would contain a sphere of radius r. for this
           we need the projection of a_1 on a unit vector perpendicular
           to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
           determine how many a_1's it will take to contain the sphere.

           Nxmax = r * length_of_b_1 / (2 Pi)

        2. keep points falling within r.

        Args:
            frac_points: All points in the lattice in fractional coordinates.
            center: Cartesian coordinates of center of sphere.
            r: radius of sphere.
            zip_results (bool): Whether to zip the results together to group by
                point, or return the raw frac_coord, dist, index arrays

        Returns:
            if zip_results:
                [(frac_coord, dist, index, supercell_image) ...] since most of the time, subsequent
                processing requires the distance, index number of the atom, or index of the image
            else:
                frac_coords, dists, inds, image
        """
        try:
            from pymatgen.optimization.neighbors import find_points_in_spheres
        except ImportError:
            return self.get_points_in_sphere_py(frac_points=frac_points, center=center, r=r, zip_results=zip_results)  # type: ignore[return-value]

        else:
            frac_points = np.ascontiguousarray(frac_points, dtype=float)
            cart_coords = np.ascontiguousarray(self.get_cartesian_coords(frac_points), dtype=float)
            center_coords = np.ascontiguousarray([center], dtype=float)
            pbc = np.ascontiguousarray(self.pbc, dtype=np.int64)
            latt_matrix = np.ascontiguousarray(self.matrix, dtype=float)

            _, indices, images, distances = find_points_in_spheres(
                all_coords=cart_coords,
                center_coords=center_coords,
                r=float(r),
                pbc=pbc,
                lattice=latt_matrix,
                tol=1e-8,
            )

            if len(indices) < 1:
                # Return empty np.array (not list or tuple) to ensure consistent return type
                # whether sphere contains points or not
                return np.array([]) if zip_results else tuple(np.array([]) for _ in range(4))  # type: ignore[return-value]
            frac_coords = frac_points[indices] + images
            if zip_results:
                return tuple(zip(frac_coords, distances, indices, images, strict=True))  # type: ignore[return-value]
            return frac_coords, distances, indices, images  # type: ignore[return-value]

    def get_points_in_sphere_py(
        self,
        frac_points: ArrayLike,
        center: ArrayLike,
        r: float,
        zip_results: bool = True,
    ) -> list[tuple[NDArray[np.float64], float, int, NDArray[np.float64]]] | list[NDArray[np.float64]]:
        """Find all points within a sphere from the point taking into account
        periodic boundary conditions. This includes sites in other periodic
        images.

        Algorithm:

        1. place sphere of radius r in crystal and determine minimum supercell
           (parallelepiped) which would contain a sphere of radius r. for this
           we need the projection of a_1 on a unit vector perpendicular
           to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
           determine how many a_1's it will take to contain the sphere.

           Nxmax = r * length_of_b_1 / (2 Pi)

        2. keep points falling within r.

        Args:
            frac_points: All points in the lattice in fractional coordinates.
            center: Cartesian coordinates of center of sphere.
            r: radius of sphere.
            zip_results (bool): Whether to zip the results together to group by
                point, or return the raw frac_coord, dist, index arrays

        Returns:
            if zip_results:
                [(frac_coord, dist, index, supercell_image) ...] since most of the time, subsequent
                processing requires the distance, index number of the atom, or index of the image
            else:
                frac_coords, dists, inds, image
        """
        cart_coords = self.get_cartesian_coords(frac_points)
        neighbors = get_points_in_spheres(
            all_coords=cart_coords,
            center_coords=np.array([center]),
            r=r,
            pbc=self.pbc,
            numerical_tol=1e-8,
            lattice=self,
            return_fcoords=True,
        )[0]
        if len(neighbors) < 1:
            return [] if zip_results else [()] * 4  # type: ignore[return-value]
        if zip_results:
            return neighbors
        return list(map(np.array, zip(*neighbors, strict=True)))

    @deprecated(get_points_in_sphere, "This is retained purely for checking purposes.")
    def get_points_in_sphere_old(
        self,
        frac_points: ArrayLike,
        center: ArrayLike,
        r: float,
        zip_results=True,
    ) -> (
        list[tuple[NDArray[np.float64], float, int, NDArray[np.float64]]]
        | tuple[list[NDArray[np.float64]], list[float], list[int], list[NDArray[np.float64]]]
    ):
        """Find all points within a sphere from the point taking into account
        periodic boundary conditions. This includes sites in other periodic
        images. Does not support partial periodic boundary conditions.

        Algorithm:

        1. place sphere of radius r in crystal and determine minimum supercell
           (parallelepiped) which would contain a sphere of radius r. for this
           we need the projection of a_1 on a unit vector perpendicular
           to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
           determine how many a_1's it will take to contain the sphere.

           Nxmax = r * length_of_b_1 / (2 Pi)

        2. keep points falling within r.

        Args:
            frac_points: All points in the lattice in fractional coordinates.
            center: Cartesian coordinates of center of sphere.
            r: radius of sphere.
            zip_results (bool): Whether to zip the results together to group by
                point, or return the raw frac_coord, dist, index arrays

        Returns:
            if zip_results:
                [(frac_coord, dist, index, supercell_image) ...] since most of the time, subsequent
                processing requires the distance, index number of the atom, or index of the image
            else:
                frac_coords, dists, inds, image
        """
        if self.pbc != (True, True, True):
            raise RuntimeError("get_points_in_sphere_old does not support partial periodic boundary conditions")
        # TODO: refactor to use lll matrix (nmax will be smaller)
        # Determine the maximum number of supercells in each direction
        # required to contain a sphere of radius n
        recp_len = np.array(self.reciprocal_lattice.lengths) / (2 * np.pi)
        nmax = float(r) * recp_len + 0.01

        # Get the fractional coordinates of the center of the sphere
        pcoords = self.get_fractional_coords(center)
        center = np.array(center)

        # Prepare the list of output atoms
        n = len(frac_points)  # type: ignore[arg-type]
        frac_coords = np.array(frac_points) % 1
        indices = np.arange(n)

        # Generate all possible images that could be within `r` of `center`
        mins = np.floor(pcoords - nmax)
        maxes = np.ceil(pcoords + nmax)
        arange = np.arange(start=mins[0], stop=maxes[0], dtype=np.int64)
        brange = np.arange(start=mins[1], stop=maxes[1], dtype=np.int64)
        crange = np.arange(start=mins[2], stop=maxes[2], dtype=np.int64)
        arange = arange[:, None] * np.array([1, 0, 0], dtype=np.int64)[None, :]  # type:ignore[assignment]
        brange = brange[:, None] * np.array([0, 1, 0], dtype=np.int64)[None, :]  # type:ignore[assignment]
        crange = crange[:, None] * np.array([0, 0, 1], dtype=np.int64)[None, :]  # type:ignore[assignment]
        images = arange[:, None, None] + brange[None, :, None] + crange[None, None, :]

        # Generate the coordinates of all atoms within these images
        shifted_coords = frac_coords[:, None, None, None, :] + images[None, :, :, :, :]

        # Determine distance from `center`
        cart_coords = self.get_cartesian_coords(frac_coords)
        cart_images = self.get_cartesian_coords(images)
        coords = cart_coords[:, None, None, None, :] + cart_images[None, :, :, :, :]
        coords -= center[None, None, None, None, :]
        coords **= 2
        d_2 = np.sum(coords, axis=4)

        # Determine which points are within `r` of `center`
        within_r = np.where(d_2 <= r**2)
        #  `within_r` now contains the coordinates of each image that is
        #    inside of the cutoff distance. It has 4 coordinates:
        #   0 - index of the image within `frac_points`
        #   1,2,3 - index of the supercell which holds the images in the x, y, z directions

        if zip_results:
            return list(
                zip(
                    shifted_coords[within_r],
                    np.sqrt(d_2[within_r]),
                    indices[within_r[0]],
                    images[within_r[1:]],
                    strict=True,
                )
            )
        return (  # type: ignore[return-value]
            shifted_coords[within_r],
            np.sqrt(d_2[within_r]),
            indices[within_r[0]],
            images[within_r[1:]],
        )

    def get_all_distances(
        self,
        frac_coords1: ArrayLike,
        frac_coords2: ArrayLike,
    ) -> NDArray[np.float64]:
        """Get the distances between two lists of coordinates taking into
        account periodic boundary conditions and the lattice. Note that this
        computes an MxN array of distances (i.e. the distance between each
        point in frac_coords1 and every coordinate in frac_coords2). This is
        different functionality from pbc_diff.

        Args:
            frac_coords1: First set of fractional coordinates. e.g. [0.5, 0.6,
                0.7] or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
                coord or any array of coords.
            frac_coords2: Second set of fractional coordinates.

        Returns:
            2d array of Cartesian distances. E.g the distance between
            frac_coords1[i] and frac_coords2[j] is distances[i,j]
        """
        _v, d2 = pbc_shortest_vectors(self, frac_coords1, frac_coords2, return_d2=True)
        return np.sqrt(d2)

    def is_hexagonal(
        self,
        hex_angle_tol: float = 5,
        hex_length_tol: float = 0.01,
    ) -> bool:
        """
        Args:
            hex_angle_tol: Angle tolerance
            hex_length_tol: Length tolerance.

        Returns:
            Whether lattice corresponds to hexagonal lattice.
        """
        lengths = self.lengths
        angles = self.angles
        right_angles = [i for i in range(3) if math.isclose(angles[i], 90, abs_tol=hex_angle_tol, rel_tol=0)]
        hex_angles = [
            idx
            for idx in range(3)
            if math.isclose(angles[idx], 60, abs_tol=hex_angle_tol, rel_tol=0)
            or math.isclose(angles[idx], 120, abs_tol=hex_angle_tol, rel_tol=0)
        ]

        return (
            len(right_angles) == 2
            and len(hex_angles) == 1
            and math.isclose(lengths[right_angles[0]], lengths[right_angles[1]], abs_tol=hex_length_tol, rel_tol=0)
        )

    def get_distance_and_image(
        self,
        frac_coords1: ArrayLike,
        frac_coords2: ArrayLike,
        jimage: ArrayLike | None = None,
    ) -> tuple[float, NDArray[np.int_]]:
        """Get distance between two frac_coords assuming periodic boundary
        conditions. If the index jimage is not specified it selects the j
        image nearest to the i atom and returns the distance and jimage
        indices in terms of lattice vector translations. If the index jimage
        is specified it returns the distance between the frac_coords1 and
        the specified jimage of frac_coords2, and the given jimage is also
        returned.

        Args:
            frac_coords1 (3x1 array): Reference frac_coords to get distance from.
            frac_coords2 (3x1 array): frac_coords to get distance from.
            jimage (3x1 array): Specific periodic image in terms of
                lattice translations, e.g. [1,0,0] implies to take periodic
                image that is one a-lattice vector away. If jimage is None,
                the image that is nearest to the site is found.

        Returns:
            tuple[float, NDArray[np.int_]]: distance and periodic lattice translations (jimage)
                of the other site for which the distance applies. This means that
                the distance between frac_coords1 and (jimage + frac_coords2) is
                equal to distance.
        """
        if jimage is None:
            v, d2 = pbc_shortest_vectors(self, frac_coords1, frac_coords2, return_d2=True)
            fc = self.get_fractional_coords(v[0][0]) + frac_coords1 - frac_coords2  # type: ignore[operator]
            fc = np.array(np.round(fc), dtype=np.int64)
            return np.sqrt(d2[0, 0]), fc

        jimage = np.array(jimage)
        mapped_vec = self.get_cartesian_coords(jimage + frac_coords2 - frac_coords1)  # type: ignore[operator]
        return np.linalg.norm(mapped_vec), jimage  # type: ignore[return-value]

    def get_miller_index_from_coords(
        self,
        coords: ArrayLike,
        coords_are_cartesian: bool = True,
        round_dp: int = 4,
        verbose: bool = True,
    ) -> tuple[int, int, int]:
        """Get the Miller index of a plane from a list of site coordinates.

        A minimum of 3 sets of coordinates are required. If more than 3 sets of
        coordinates are given, the best plane that minimizes the distance to all
        Points will be calculated.

        Args:
            coords (iterable): A list or numpy array of coordinates. Can be
                Cartesian or fractional coordinates. If more than three sets of
                coordinates are provided, the best plane that minimises the
                distance to all sites will be calculated.
            coords_are_cartesian (bool, optional): Whether the coordinates are
                in Cartesian space. If using fractional coordinates set to
                False.
            round_dp (int, optional): The number of decimal places to round the
                miller index to.
            verbose (bool, optional): Whether to print warnings.

        Returns:
            tuple: The Miller index.
        """
        if coords_are_cartesian:
            coords = [self.get_fractional_coords(c) for c in coords]  # type: ignore[union-attr]

        coords = np.asarray(coords)
        g = coords.sum(axis=0) / coords.shape[0]

        # Run singular value decomposition
        _, _, vh = np.linalg.svd(coords - g)

        # Get unitary normal vector
        u_norm = vh[2, :]
        return get_integer_index(u_norm, round_dp=round_dp, verbose=verbose)  # type: ignore[return-value]

    def get_recp_symmetry_operation(self, symprec: float = 0.01) -> list[SymmOp]:
        """Find the symmetric operations of the reciprocal lattice,
        to be used for hkl transformations.

        Args:
            symprec: default is 0.001.
        """
        recp_lattice = self.reciprocal_lattice_crystallographic
        # Get symmetry operations from input conventional unit cell
        # Need to make sure recp lattice is big enough, otherwise symmetry
        # determination will fail. We set the overall volume to 1.
        recp_lattice = recp_lattice.scale(1)
        # Need a localized import of structure to build a
        # pseudo empty lattice for SpacegroupAnalyzer

        from pymatgen.core.structure import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        recp = Structure(recp_lattice, ["H"], [[0, 0, 0]])  # type:ignore[list-item]
        # Create a function that uses the symmetry operations in the
        # structure to find Miller indices that might give repetitive slabs
        analyzer = SpacegroupAnalyzer(recp, symprec=symprec)
        return analyzer.get_symmetry_operations()


def get_integer_index(
    miller_index: tuple[int, ...],
    round_dp: int = 4,
    verbose: bool = True,
) -> tuple[int, ...]:
    """Attempt to convert a vector of floats to whole numbers.

    Args:
        miller_index (tuple[int, ...]): Miller index.
        round_dp (int, optional): The number of decimal places to round the
            miller index to.
        verbose (bool, optional): Whether to print warnings.

    Returns:
        tuple[int, ...]: The Miller index.
    """
    mi = np.asarray(miller_index)
    # Deal with the case we have small irregular floats
    # that are all equal or factors of each other
    mi /= min(m for m in mi if m != 0)
    mi /= np.max(np.abs(mi))

    # Deal with the case where we have nice fractions
    md = [Fraction(n).limit_denominator(12).denominator for n in mi]
    mi *= reduce(operator.mul, md)
    int_miller_index = np.round(mi, 1).astype(int)
    mi /= np.abs(reduce(math.gcd, int_miller_index))

    # Round to a reasonable precision
    mi = np.array([round(h, round_dp) for h in mi])

    # Need to recalculate this after rounding as values may have changed
    int_miller_index = np.round(mi, 1).astype(int)
    if np.any(np.abs(mi - int_miller_index) > 1e-6) and verbose:
        warnings.warn("Non-integer encountered in Miller index", stacklevel=2)
    else:
        mi = int_miller_index

    # Minimise the number of negative indexes
    mi += 0  # converts -0 to 0

    def n_minus(index):
        return sum(h < 0 for h in index)

    if n_minus(mi) > n_minus(mi * -1):
        mi *= -1

    # If only one index is negative, make sure it is the smallest
    # e.g. (-2 1 0) -> (2 -1 0)
    if sum(mi != 0) == 2 and n_minus(mi) == 1 and abs(min(mi)) > max(mi):
        mi *= -1

    return tuple(mi)


def get_points_in_spheres(
    all_coords: NDArray[np.float64],
    center_coords: NDArray[np.float64],
    r: float,
    pbc: bool | list[bool] | tuple[bool, bool, bool] = True,
    numerical_tol: float = 1e-8,
    lattice: Lattice | None = None,
    return_fcoords: bool = False,
) -> list[list[tuple[NDArray[np.float64], float, int, NDArray[np.float64]]]]:
    """For each point in `center_coords`, get all the neighboring points
    in `all_coords` that are within the cutoff radius `r`.

    Args:
        all_coords: (list of Cartesian coordinates) all available points
        center_coords: (list of Cartesian coordinates) all centering points
        r: (float) cutoff radius
        pbc: (bool or a list of bool) whether to set periodic boundaries
        numerical_tol: (float) numerical tolerance
        lattice: (Lattice) lattice to consider when PBC is enabled
        return_fcoords: (bool) whether to return fractional coords when pbc is set.

    Returns:
        List[List[Tuple[coords, distance, index, image]]]
    """
    if isinstance(pbc, bool):
        pbc = [pbc] * 3
    _pbc = np.array(pbc, dtype=bool)
    if return_fcoords and lattice is None:
        raise ValueError("Lattice needs to be supplied to compute fractional coordinates")
    center_coords_min = np.min(center_coords, axis=0)
    center_coords_max = np.max(center_coords, axis=0)

    # The lower bound of all considered atom coords
    global_min = center_coords_min - r - numerical_tol
    global_max = center_coords_max + r + numerical_tol
    if np.any(_pbc):
        if lattice is None:
            raise ValueError("Lattice needs to be supplied when considering periodic boundary")
        recp_len = np.array(lattice.reciprocal_lattice.lengths)
        maxr = np.ceil((r + 0.15) * recp_len / (2 * math.pi))
        frac_coords = lattice.get_fractional_coords(center_coords)
        nmin_temp = np.floor(np.min(frac_coords, axis=0)) - maxr
        nmax_temp = np.ceil(np.max(frac_coords, axis=0)) + maxr
        nmin = np.zeros_like(nmin_temp)
        nmin[_pbc] = nmin_temp[_pbc]
        nmax = np.ones_like(nmax_temp)
        nmax[_pbc] = nmax_temp[_pbc]
        all_ranges = [np.arange(x, y, dtype="int64") for x, y in zip(nmin, nmax, strict=True)]
        matrix = lattice.matrix

        # Temporarily hold the fractional coordinates
        image_offsets = lattice.get_fractional_coords(all_coords)
        all_frac_coords = []

        # Only wrap periodic boundary
        for kk in range(3):
            if _pbc[kk]:
                all_frac_coords.append(np.mod(image_offsets[:, kk : kk + 1], 1))
            else:
                all_frac_coords.append(image_offsets[:, kk : kk + 1])
        all_frac_coords = np.concatenate(all_frac_coords, axis=1)
        image_offsets -= all_frac_coords
        coords_in_cell = np.dot(all_frac_coords, matrix)

        # Filter out those beyond max range
        valid_coords = []
        valid_images = []
        valid_indices = []
        for image in itertools.product(*all_ranges):
            coords = np.dot(image, matrix) + coords_in_cell
            valid_index_bool = np.all(
                np.bitwise_and(coords > global_min[None, :], coords < global_max[None, :]),
                axis=1,
            )
            ind = np.arange(len(all_coords))
            if np.any(valid_index_bool):
                valid_coords.append(coords[valid_index_bool])
                valid_images.append(np.tile(image, [np.sum(valid_index_bool), 1]) - image_offsets[valid_index_bool])
                valid_indices.extend([k for k in ind if valid_index_bool[k]])
        if not valid_coords:
            return [[]] * len(center_coords)
        valid_coords = np.concatenate(valid_coords, axis=0)
        valid_images = np.concatenate(valid_images, axis=0)  # type:ignore[assignment]

    else:
        valid_coords = all_coords  # type: ignore[assignment]
        valid_images = [[0, 0, 0]] * len(valid_coords)  # type: ignore[list-item]
        valid_indices = np.arange(len(valid_coords))  # type: ignore[assignment]

    # Divide the valid 3D space into cubes and compute the cube ids
    all_cube_index = _compute_cube_index(valid_coords, global_min, r)  # type: ignore[arg-type]
    nx, ny, nz = _compute_cube_index(global_max, global_min, r) + 1
    all_cube_index = _three_to_one(all_cube_index, ny, nz)
    site_cube_index = _three_to_one(_compute_cube_index(center_coords, global_min, r), ny, nz)

    # Create cube index to coordinates, images, and indices map
    cube_to_coords: dict[int, list] = defaultdict(list)
    cube_to_images: dict[int, list] = defaultdict(list)
    cube_to_indices: dict[int, list] = defaultdict(list)
    for ii, jj, kk, ll in zip(all_cube_index.ravel(), valid_coords, valid_images, valid_indices, strict=True):  # type: ignore[assignment]
        cube_to_coords[ii].append(jj)  # type: ignore[index]
        cube_to_images[ii].append(kk)  # type: ignore[index]
        cube_to_indices[ii].append(ll)  # type: ignore[index]

    # Find all neighboring cubes for each atom in the lattice cell
    site_neighbors = find_neighbors(site_cube_index, nx, ny, nz)
    neighbors: list[list[tuple[NDArray[np.float64], float, int, NDArray[np.float64]]]] = []

    for cc, jj in zip(center_coords, site_neighbors, strict=True):
        l1 = np.array(_three_to_one(jj, ny, nz), dtype=np.int64).ravel()
        # Use the cube index map to find the all the neighboring
        # coords, images, and indices
        ks = [k for k in l1 if k in cube_to_coords]
        if not ks:
            neighbors.append([])
            continue
        nn_coords = np.concatenate([cube_to_coords[k] for k in ks], axis=0)  # type:ignore[index]
        nn_images = itertools.chain(*(cube_to_images[k] for k in ks))  # type:ignore[index]
        nn_indices = itertools.chain(*(cube_to_indices[k] for k in ks))  # type:ignore[index]
        distances = np.linalg.norm(nn_coords - cc[None, :], axis=1)  # type:ignore[index]
        nns: list[tuple[NDArray[np.float64], float, int, NDArray[np.float64]]] = []
        for coord, index, image, dist in zip(nn_coords, nn_indices, nn_images, distances, strict=True):
            # Filtering out all sites that are beyond the cutoff
            # Here there is no filtering of overlapping sites
            if dist < r + numerical_tol:
                if return_fcoords and (lattice is not None):
                    coord = np.round(lattice.get_fractional_coords(coord), 10)
                nn = (coord, float(dist), int(index), image)
                nns.append(nn)  # type:ignore[arg-type]
        neighbors.append(nns)
    return neighbors


# The following internal functions are used in the get_points_in_sphere method
def _compute_cube_index(
    coords: NDArray[np.float64],
    global_min: float,
    radius: float,
) -> NDArray[np.int_]:
    """Compute the cube index from coordinates
    Args:
        coords: (nx3 array) atom coordinates
        global_min: (float) lower boundary of coordinates
        radius: (float) cutoff radius.

    Returns:
        NDArray[np.float_]: nx3 array int indices
    """
    return np.array(np.floor((coords - global_min) / radius), dtype=np.int_)


def _one_to_three(label1d: NDArray[np.int_], ny: int, nz: int) -> NDArray[np.int_]:
    """Convert a 1D index array to 3D index array.

    Args:
        label1d: (array) 1D index array
        ny: (int) number of cells in y direction
        nz: (int) number of cells in z direction

    Returns:
        NDArray[np.float_]: nx3 array int indices
    """
    last = np.mod(label1d, nz)
    second = np.mod((label1d - last) / nz, ny)
    first = (label1d - last - second * nz) / (ny * nz)
    return np.concatenate([first, second, last], axis=1)


def _three_to_one(label3d: NDArray[np.int_], ny: int, nz: int) -> NDArray[np.int_]:
    """The reverse of _one_to_three."""
    return np.array(label3d[:, 0] * ny * nz + label3d[:, 1] * nz + label3d[:, 2]).reshape((-1, 1))


def find_neighbors(label: NDArray[np.int_], nx: int, ny: int, nz: int) -> list[NDArray[np.int_]]:
    """Given a cube index, find the neighbor cube indices.

    Args:
        label: (array) (n,) or (n x 3) indice array
        nx: (int) number of cells in y direction
        ny: (int) number of cells in y direction
        nz: (int) number of cells in z direction

    Returns:
        Neighbor cell indices.
    """
    array = [[-1, 0, 1]] * 3
    neighbor_vectors = np.array(list(itertools.product(*array)), dtype=np.int64)
    label3d = _one_to_three(label, ny, nz) if np.shape(label)[1] == 1 else label
    all_labels = label3d[:, None, :] - neighbor_vectors[None, :, :]
    filtered_labels = []
    # Filter out out-of-bound labels i.e., label < 0
    for labels in all_labels:
        ind = (labels[:, 0] < nx) * (labels[:, 1] < ny) * (labels[:, 2] < nz) * np.all(labels > -1e-5, axis=1)
        filtered_labels.append(labels[ind])
    return filtered_labels
