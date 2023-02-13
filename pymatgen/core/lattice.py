# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Defines the classes relating to 3D lattices.
"""

from __future__ import annotations

import collections
import itertools
import math
import warnings
from fractions import Fraction
from functools import reduce
from typing import Iterator, Sequence

import numpy as np
from monty.dev import deprecated
from monty.json import MSONable
from numpy import dot, pi, transpose
from numpy.linalg import inv

from pymatgen.util.coord import pbc_shortest_vectors
from pymatgen.util.num import abs_cap
from pymatgen.util.typing import ArrayLike

__author__ = "Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2011, The Materials Project"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"


class Lattice(MSONable):
    """
    A lattice object. Essentially a matrix with conversion matrices. In
    general, it is assumed that length units are in Angstroms and angles are in
    degrees unless otherwise stated.
    """

    # Properties lazily generated for efficiency.

    def __init__(self, matrix: ArrayLike, pbc: tuple[bool, bool, bool] = (True, True, True)):
        """
        Create a lattice from any sequence of 9 numbers. Note that the sequence
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
                E.g., [[10, 0, 0], [20, 10, 0], [0, 0, 30]] specifies a lattice
                with lattice vectors [10, 0, 0], [20, 10, 0] and [0, 0, 30].
            pbc: a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.
        """
        m = np.array(matrix, dtype=np.float64).reshape((3, 3))
        m.setflags(write=False)
        self._matrix: np.ndarray = m
        self._inv_matrix: np.ndarray | None = None
        self._diags = None
        self._lll_matrix_mappings: dict[float, tuple[np.ndarray, np.ndarray]] = {}
        self._lll_inverse = None
        self._pbc = tuple(pbc)

    @property
    def lengths(self) -> tuple[float, float, float]:
        """
        Lattice lengths.

        :return: The lengths (a, b, c) of the lattice.
        """
        return tuple(np.sqrt(np.sum(self._matrix**2, axis=1)).tolist())  # type: ignore

    @property
    def angles(self) -> tuple[float, float, float]:
        """
        Lattice angles.

        :return: The angles (alpha, beta, gamma) of the lattice.
        """
        m = self._matrix
        lengths = self.lengths
        angles = np.zeros(3)
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            angles[i] = abs_cap(dot(m[j], m[k]) / (lengths[j] * lengths[k]))
        angles = np.arccos(angles) * 180.0 / pi
        return tuple(angles.tolist())  # type: ignore

    @property
    def is_orthogonal(self) -> bool:
        """
        :return: Whether all angles are 90 degrees.
        """
        return all(abs(a - 90) < 1e-5 for a in self.angles)

    def __format__(self, fmt_spec=""):
        """
        Support format printing. Supported formats are:

        1. "l" for a list format that can be easily copied and pasted, e.g.,
           ".3fl" prints something like
           "[[10.000, 0.000, 0.000], [0.000, 10.000, 0.000], [0.000, 0.000, 10.000]]"
        2. "p" for lattice parameters ".1fp" prints something like
           "{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}"
        3. Default will simply print a 3x3 matrix form. E.g.,
           10.000 0.000 0.000
           0.000 10.000 0.000
           0.000 0.000 10.000
        """
        m = self._matrix.tolist()
        if fmt_spec.endswith("l"):
            fmt = "[[{}, {}, {}], [{}, {}, {}], [{}, {}, {}]]"
            fmt_spec = fmt_spec[:-1]
        elif fmt_spec.endswith("p"):
            fmt = "{{{}, {}, {}, {}, {}, {}}}"
            fmt_spec = fmt_spec[:-1]
            m = (self.lengths, self.angles)
        else:
            fmt = "{} {} {}\n{} {} {}\n{} {} {}"
        return fmt.format(*(format(c, fmt_spec) for row in m for c in row))

    def copy(self):
        """Deep copy of self."""
        return self.__class__(self.matrix.copy(), pbc=self.pbc)

    @property
    def matrix(self) -> np.ndarray:
        """Copy of matrix representing the Lattice"""
        return self._matrix

    @property
    def pbc(self) -> tuple[bool, bool, bool]:
        """Tuple defining the periodicity of the Lattice"""
        return self._pbc  # type: ignore

    @property
    def is_3d_periodic(self) -> bool:
        """True if the Lattice is periodic in all directions"""
        return all(self._pbc)

    @property
    def inv_matrix(self) -> np.ndarray:
        """
        Inverse of lattice matrix.
        """
        if self._inv_matrix is None:
            self._inv_matrix = inv(self._matrix)
            self._inv_matrix.setflags(write=False)
        return self._inv_matrix

    @property
    def metric_tensor(self) -> np.ndarray:
        """
        The metric tensor of the lattice.
        """
        return dot(self._matrix, self._matrix.T)

    def get_cartesian_coords(self, fractional_coords: ArrayLike) -> np.ndarray:
        """
        Returns the Cartesian coordinates given fractional coordinates.

        Args:
            fractional_coords (3x1 array): Fractional coords.

        Returns:
            Cartesian coordinates
        """
        return dot(fractional_coords, self._matrix)

    def get_fractional_coords(self, cart_coords: ArrayLike) -> np.ndarray:
        """
        Returns the fractional coordinates given Cartesian coordinates.

        Args:
            cart_coords (3x1 array): Cartesian coords.

        Returns:
            Fractional coordinates.
        """
        return dot(cart_coords, self.inv_matrix)

    def get_vector_along_lattice_directions(self, cart_coords: ArrayLike) -> np.ndarray:
        """
        Returns the coordinates along lattice directions given Cartesian coordinates.

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
        return self.lengths * self.get_fractional_coords(cart_coords)  # type: ignore

    def d_hkl(self, miller_index: ArrayLike) -> float:
        """
        Returns the distance between the hkl plane and the origin

        Args:
            miller_index ([h,k,l]): Miller index of plane

        Returns:
            d_hkl (float)
        """
        gstar = self.reciprocal_lattice_crystallographic.metric_tensor
        hkl = np.array(miller_index)
        return 1 / ((dot(dot(hkl, gstar), hkl.T)) ** (1 / 2))

    @staticmethod
    def cubic(a: float, pbc: tuple[bool, bool, bool] = (True, True, True)) -> Lattice:
        """
        Convenience constructor for a cubic lattice.

        Args:
            a (float): The *a* lattice parameter of the cubic cell.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Cubic lattice of dimensions a x a x a.
        """
        return Lattice([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]], pbc)

    @staticmethod
    def tetragonal(a: float, c: float, pbc: tuple[bool, bool, bool] = (True, True, True)) -> Lattice:
        """
        Convenience constructor for a tetragonal lattice.

        Args:
            a (float): *a* lattice parameter of the tetragonal cell.
            c (float): *c* lattice parameter of the tetragonal cell.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Tetragonal lattice of dimensions a x a x c.
        """
        return Lattice.from_parameters(a, a, c, 90, 90, 90, pbc=pbc)

    @staticmethod
    def orthorhombic(a: float, b: float, c: float, pbc: tuple[bool, bool, bool] = (True, True, True)) -> Lattice:
        """
        Convenience constructor for an orthorhombic lattice.

        Args:
            a (float): *a* lattice parameter of the orthorhombic cell.
            b (float): *b* lattice parameter of the orthorhombic cell.
            c (float): *c* lattice parameter of the orthorhombic cell.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Orthorhombic lattice of dimensions a x b x c.
        """
        return Lattice.from_parameters(a, b, c, 90, 90, 90, pbc=pbc)

    @staticmethod
    def monoclinic(
        a: float, b: float, c: float, beta: float, pbc: tuple[bool, bool, bool] = (True, True, True)
    ) -> Lattice:
        """
        Convenience constructor for a monoclinic lattice.

        Args:
            a (float): *a* lattice parameter of the monoclinc cell.
            b (float): *b* lattice parameter of the monoclinc cell.
            c (float): *c* lattice parameter of the monoclinc cell.
            beta (float): *beta* angle between lattice vectors b and c in
                degrees.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Monoclinic lattice of dimensions a x b x c with non right-angle
            beta between lattice vectors a and c.
        """
        return Lattice.from_parameters(a, b, c, 90, beta, 90, pbc=pbc)

    @staticmethod
    def hexagonal(a: float, c: float, pbc: tuple[bool, bool, bool] = (True, True, True)) -> Lattice:
        """
        Convenience constructor for a hexagonal lattice.

        Args:
            a (float): *a* lattice parameter of the hexagonal cell.
            c (float): *c* lattice parameter of the hexagonal cell.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Hexagonal lattice of dimensions a x a x c.
        """
        return Lattice.from_parameters(a, a, c, 90, 90, 120, pbc=pbc)

    @staticmethod
    def rhombohedral(a: float, alpha: float, pbc: tuple[bool, bool, bool] = (True, True, True)) -> Lattice:
        """
        Convenience constructor for a rhombohedral lattice.

        Args:
            a (float): *a* lattice parameter of the rhombohedral cell.
            alpha (float): Angle for the rhombohedral lattice in degrees.
            pbc (tuple): a tuple defining the periodic boundary conditions along the three
                axis of the lattice. If None periodic in all directions.

        Returns:
            Rhombohedral lattice of dimensions a x a x a.
        """
        return Lattice.from_parameters(a, a, a, alpha, alpha, alpha, pbc=pbc)

    @classmethod
    def from_parameters(
        cls,
        a: float,
        b: float,
        c: float,
        alpha: float,
        beta: float,
        gamma: float,
        vesta: bool = False,
        pbc: tuple[bool, bool, bool] = (True, True, True),
    ):
        """
        Create a Lattice using unit cell lengths (in Angstrom) and angles (in degrees).

        Args:
            a (float): *a* lattice parameter.
            b (float): *b* lattice parameter.
            c (float): *c* lattice parameter.
            alpha (float): *alpha* angle in degrees.
            beta (float): *beta* angle in degrees.
            gamma (float): *gamma* angle in degrees.
            vesta: True if you import Cartesian coordinates from VESTA.
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
            # Sometimes rounding errors result in values slightly > 1.
            val = abs_cap(val)
            gamma_star = np.arccos(val)

            vector_a = [a * sin_beta, 0.0, a * cos_beta]
            vector_b = [
                -b * sin_alpha * np.cos(gamma_star),
                b * sin_alpha * np.sin(gamma_star),
                b * cos_alpha,
            ]
            vector_c = [0.0, 0.0, float(c)]

        return Lattice([vector_a, vector_b, vector_c], pbc)

    @classmethod
    def from_dict(cls, d: dict, fmt: str | None = None, **kwargs):
        """
        Create a Lattice from a dictionary containing the a, b, c, alpha, beta,
        and gamma parameters if fmt is None.

        If fmt == "abivars", the function build a `Lattice` object from a
        dictionary with the Abinit variables `acell` and `rprim` in Bohr.
        If acell is not given, the Abinit default is used i.e. [1,1,1] Bohr

        Example:
            Lattice.from_dict(fmt="abivars", acell=3*[10], rprim=np.eye(3))
        """
        if fmt == "abivars":
            # pylint: disable=C0415
            from pymatgen.io.abinit.abiobjects import lattice_from_abivars

            kwargs.update(d)
            return lattice_from_abivars(cls=cls, **kwargs)

        pbc = d.get("pbc", (True, True, True))
        if "matrix" in d:
            return cls(d["matrix"], pbc=pbc)
        return cls.from_parameters(d["a"], d["b"], d["c"], d["alpha"], d["beta"], d["gamma"], pbc=pbc)

    @property
    def a(self) -> float:
        """
        *a* lattice parameter.
        """
        return self.lengths[0]

    @property
    def b(self) -> float:
        """
        *b* lattice parameter.
        """
        return self.lengths[1]

    @property
    def c(self) -> float:
        """
        *c* lattice parameter.
        """
        return self.lengths[2]

    @property
    def abc(self) -> tuple[float, float, float]:
        """
        Lengths of the lattice vectors, i.e. (a, b, c)
        """
        return self.lengths

    @property
    def alpha(self) -> float:
        """
        Angle alpha of lattice in degrees.
        """
        return self.angles[0]

    @property
    def beta(self) -> float:
        """
        Angle beta of lattice in degrees.
        """
        return self.angles[1]

    @property
    def gamma(self) -> float:
        """
        Angle gamma of lattice in degrees.
        """
        return self.angles[2]

    @property
    def volume(self) -> float:
        """
        Volume of the unit cell.
        """
        m = self._matrix
        return float(abs(dot(np.cross(m[0], m[1]), m[2])))

    @property
    def parameters(self) -> tuple[float, float, float, float, float, float]:
        """
        Returns: (a, b, c, alpha, beta, gamma).
        """
        return (*self.lengths, *self.angles)

    @property
    def reciprocal_lattice(self) -> Lattice:
        """
        Return the reciprocal lattice. Note that this is the standard
        reciprocal lattice used for solid state physics with a factor of 2 *
        pi. If you are looking for the crystallographic reciprocal lattice,
        use the reciprocal_lattice_crystallographic property.
        The property is lazily generated for efficiency.
        """
        v = np.linalg.inv(self._matrix).T
        return Lattice(v * 2 * np.pi)

    @property
    def reciprocal_lattice_crystallographic(self) -> Lattice:
        """
        Returns the *crystallographic* reciprocal lattice, i.e., no factor of
        2 * pi.
        """
        return Lattice(self.reciprocal_lattice.matrix / (2 * np.pi))

    @property
    def lll_matrix(self) -> np.ndarray:
        """
        :return: The matrix for LLL reduction
        """
        if 0.75 not in self._lll_matrix_mappings:
            self._lll_matrix_mappings[0.75] = self._calculate_lll()
        return self._lll_matrix_mappings[0.75][0]

    @property
    def lll_mapping(self) -> np.ndarray:
        """
        :return: The mapping between the LLL reduced lattice and the original
            lattice.
        """
        if 0.75 not in self._lll_matrix_mappings:
            self._lll_matrix_mappings[0.75] = self._calculate_lll()
        return self._lll_matrix_mappings[0.75][1]

    @property
    def lll_inverse(self) -> np.ndarray:
        """
        :return: Inverse of self.lll_mapping.
        """
        return np.linalg.inv(self.lll_mapping)

    @property
    def selling_vector(self) -> np.ndarray:
        """
        Returns the (1,6) array of Selling Scalars.
        """
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

        reduction_matrices = np.array(
            [
                np.array(
                    [
                        [-1, 0, 0, 0, 0, 0],
                        [1, 1, 0, 0, 0, 0],
                        [1, 0, 0, 0, 1, 0],
                        [-1, 0, 0, 1, 0, 0],
                        [1, 0, 1, 0, 0, 0],
                        [1, 0, 0, 0, 0, 1],
                    ]
                ),
                np.array(
                    [
                        [1, 1, 0, 0, 0, 0],
                        [0, -1, 0, 0, 0, 0],
                        [0, 1, 0, 1, 0, 0],
                        [0, 1, 1, 0, 0, 0],
                        [0, -1, 0, 0, 1, 0],
                        [0, 1, 0, 0, 0, 1],
                    ]
                ),
                np.array(
                    [
                        [1, 0, 1, 0, 0, 0],
                        [0, 0, 1, 1, 0, 0],
                        [0, 0, -1, 0, 0, 0],
                        [0, 1, 1, 0, 0, 0],
                        [0, 0, 1, 0, 1, 0],
                        [0, 0, -1, 0, 0, 1],
                    ]
                ),
                np.array(
                    [
                        [1, 0, 0, -1, 0, 0],
                        [0, 0, 1, 1, 0, 0],
                        [0, 1, 0, 1, 0, 0],
                        [0, 0, 0, -1, 0, 0],
                        [0, 0, 0, 1, 1, 0],
                        [0, 0, 0, 1, 0, 1],
                    ]
                ),
                np.array(
                    [
                        [0, 0, 1, 0, 1, 0],
                        [0, 1, 0, 0, -1, 0],
                        [1, 0, 0, 0, 1, 0],
                        [0, 0, 0, 1, 1, 0],
                        [0, 0, 0, 0, -1, 0],
                        [0, 0, 0, 0, 1, 1],
                    ]
                ),
                np.array(
                    [
                        [0, 1, 0, 0, 0, 1],
                        [1, 0, 0, 0, 0, 1],
                        [0, 0, 1, 0, 0, -1],
                        [0, 0, 0, 1, 0, 1],
                        [0, 0, 0, 0, 1, 1],
                        [0, 0, 0, 0, 0, -1],
                    ]
                ),
            ]
        )
        while np.greater(np.max(selling_vector), 0):
            max_index = selling_vector.argmax()
            selling_vector = np.dot(reduction_matrices[max_index], selling_vector)

        return selling_vector

    def selling_dist(self, other):
        """
        Returns the minimum Selling distance between two lattices.
        """
        vcp_matrices = [
            np.array(
                [
                    [-1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                ]
            ),
            np.array(
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, -1, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                ]
            ),
            np.array(
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, -1, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                ]
            ),
            np.array(
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, -1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                ]
            ),
            np.array(
                [
                    [0, 0, 1, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, -1, 0],
                    [0, 0, 0, 0, 0, 1],
                ]
            ),
            np.array(
                [
                    [0, 1, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, -1],
                ]
            ),
        ]

        reflection_matrices = [
            np.array(
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                ]
            ),
            np.array(
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 1, 0],
                ]
            ),
            np.array(
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 1, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 1, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                ]
            ),
            np.array(
                [
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 1, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 1, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 1, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 1, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 1, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 1, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 1, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 1, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 1, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 1, 0, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 1, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 1, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 0, 1, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                    [1, 0, 0, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 0, 0, 1],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 0, 0, 1],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [1, 0, 0, 0, 0, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 1, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                ]
            ),
            np.array(
                [
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 1, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                ]
            ),
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

        return min(np.linalg.norm(reflection - selling2) for reflection in all_reflections)

    def __repr__(self):
        outs = [
            "Lattice",
            "    abc : " + " ".join(map(repr, self.lengths)),
            " angles : " + " ".join(map(repr, self.angles)),
            " volume : " + repr(self.volume),
            "      A : " + " ".join(map(repr, self._matrix[0])),
            "      B : " + " ".join(map(repr, self._matrix[1])),
            "      C : " + " ".join(map(repr, self._matrix[2])),
            "    pbc : " + " ".join(map(repr, self._pbc)),
        ]
        return "\n".join(outs)

    def __eq__(self, other: object) -> bool:
        """
        A lattice is considered to be equal to another if the internal matrix
        representation satisfies np.allclose(matrix1, matrix2) to be True and
        share the same periodicity.
        """
        if not hasattr(other, "matrix") or not hasattr(other, "pbc"):
            return NotImplemented
        # shortcut the np.allclose if the memory addresses are the same
        # (very common in Structure.from_sites)
        if self is other:
            return True
        return np.allclose(self.matrix, other.matrix) and self.pbc == other.pbc  # type: ignore

    def __hash__(self):
        return 7

    def __str__(self):
        return "\n".join(" ".join([f"{i:.6f}" for i in row]) for row in self._matrix)

    def as_dict(self, verbosity: int = 0) -> dict:
        """
        Json-serialization dict representation of the Lattice.

        Args:
            verbosity (int): Verbosity level. Default of 0 only includes the
                matrix representation. Set to 1 for more details.
        """
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "matrix": self._matrix.tolist(),
            "pbc": self._pbc,
        }
        a, b, c, alpha, beta, gamma = self.parameters
        if verbosity > 0:
            d.update(
                {
                    "a": a,
                    "b": b,
                    "c": c,
                    "alpha": alpha,
                    "beta": beta,
                    "gamma": gamma,
                    "volume": self.volume,
                }
            )

        return d

    def find_all_mappings(
        self,
        other_lattice: Lattice,
        ltol: float = 1e-5,
        atol: float = 1,
        skip_rotation_matrix: bool = False,
    ) -> Iterator[tuple[Lattice, np.ndarray | None, np.ndarray]]:
        """
        Finds all mappings between current lattice and another lattice.

        Args:
            other_lattice (Lattice): Another lattice that is equivalent to
                this one.
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
        lengths = other_lattice.lengths
        (alpha, beta, gamma) = other_lattice.angles

        frac, dist, _, _ = self.get_points_in_sphere(
            [[0, 0, 0]], [0, 0, 0], max(lengths) * (1 + ltol), zip_results=False
        )
        cart = self.get_cartesian_coords(frac)  # type: ignore
        # this can't be broadcast because they're different lengths
        inds = [np.logical_and(dist / l < 1 + ltol, dist / l > 1 / (1 + ltol)) for l in lengths]  # type: ignore
        c_a, c_b, c_c = (cart[i] for i in inds)
        f_a, f_b, f_c = (frac[i] for i in inds)  # type: ignore
        l_a, l_b, l_c = (np.sum(c**2, axis=-1) ** 0.5 for c in (c_a, c_b, c_c))

        def get_angles(v1, v2, l1, l2):
            x = np.inner(v1, v2) / l1[:, None] / l2
            x[x > 1] = 1
            x[x < -1] = -1
            angles = np.arccos(x) * 180.0 / pi
            return angles

        alphab = np.abs(get_angles(c_b, c_c, l_b, l_c) - alpha) < atol
        betab = np.abs(get_angles(c_a, c_c, l_a, l_c) - beta) < atol
        gammab = np.abs(get_angles(c_a, c_b, l_a, l_b) - gamma) < atol

        for i, all_j in enumerate(gammab):
            inds = np.logical_and(all_j[:, None], np.logical_and(alphab, betab[i][None, :]))
            for j, k in np.argwhere(inds):
                scale_m = np.array((f_a[i], f_b[j], f_c[k]), dtype=int)  # type: ignore
                if abs(np.linalg.det(scale_m)) < 1e-8:  # type: ignore
                    continue

                aligned_m = np.array((c_a[i], c_b[j], c_c[k]))

                if skip_rotation_matrix:
                    rotation_m = None
                else:
                    rotation_m = np.linalg.solve(aligned_m, other_lattice.matrix)

                yield Lattice(aligned_m), rotation_m, scale_m

    def find_mapping(
        self,
        other_lattice: Lattice,
        ltol: float = 1e-5,
        atol: float = 1,
        skip_rotation_matrix: bool = False,
    ) -> tuple[Lattice, np.ndarray | None, np.ndarray] | None:
        """
        Finds a mapping between current lattice and another lattice. There
        are an infinite number of choices of basis vectors for two entirely
        equivalent lattices. This method returns a mapping that maps
        other_lattice to this lattice.

        Args:
            other_lattice (Lattice): Another lattice that is equivalent to
                this one.
            ltol (float): Tolerance for matching lengths. Defaults to 1e-5.
            atol (float): Tolerance for matching angles. Defaults to 1.

        Returns:
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
        for x in self.find_all_mappings(other_lattice, ltol, atol, skip_rotation_matrix=skip_rotation_matrix):
            return x
        return None

    def get_lll_reduced_lattice(self, delta: float = 0.75) -> Lattice:
        """
        :param delta: Delta parameter.
        :return: LLL reduced Lattice.
        """
        if delta not in self._lll_matrix_mappings:
            self._lll_matrix_mappings[delta] = self._calculate_lll()
        return Lattice(self._lll_matrix_mappings[delta][0])

    def _calculate_lll(self, delta: float = 0.75) -> tuple[np.ndarray, np.ndarray]:
        """
        Performs a Lenstra-Lenstra-Lovasz lattice basis reduction to obtain a
        c-reduced basis. This method returns a basis which is as "good" as
        possible, with "good" defined by orthongonality of the lattice vectors.

        This basis is used for all the periodic boundary condition calculations.

        Args:
            delta (float): Reduction parameter. Default of 0.75 is usually
                fine.

        Returns:
            Reduced lattice matrix, mapping to get to that lattice.
        """
        # Transpose the lattice matrix first so that basis vectors are columns.
        # Makes life easier.
        # pylint: disable=E1136,E1137,E1126
        a = self._matrix.copy().T

        b = np.zeros((3, 3))  # Vectors after the Gram-Schmidt process
        u = np.zeros((3, 3))  # Gram-Schmidt coefficients
        m = np.zeros(3)  # These are the norm squared of each vec.

        b[:, 0] = a[:, 0]
        m[0] = dot(b[:, 0], b[:, 0])
        for i in range(1, 3):
            u[i, 0:i] = dot(a[:, i].T, b[:, 0:i]) / m[0:i]
            b[:, i] = a[:, i] - dot(b[:, 0:i], u[i, 0:i].T)
            m[i] = dot(b[:, i], b[:, i])

        k = 2

        mapping = np.identity(3, dtype=np.double)
        while k <= 3:
            # Size reduction.
            for i in range(k - 1, 0, -1):
                q = round(u[k - 1, i - 1])
                if q != 0:
                    # Reduce the k-th basis vector.
                    a[:, k - 1] = a[:, k - 1] - q * a[:, i - 1]
                    mapping[:, k - 1] = mapping[:, k - 1] - q * mapping[:, i - 1]
                    uu = list(u[i - 1, 0 : (i - 1)])
                    uu.append(1)
                    # Update the GS coefficients.
                    u[k - 1, 0:i] = u[k - 1, 0:i] - q * np.array(uu)

            # Check the Lovasz condition.
            if dot(b[:, k - 1], b[:, k - 1]) >= (delta - abs(u[k - 1, k - 2]) ** 2) * dot(b[:, (k - 2)], b[:, (k - 2)]):
                # Increment k if the Lovasz condition holds.
                k += 1
            else:
                # If the Lovasz condition fails,
                # swap the k-th and (k-1)-th basis vector
                v = a[:, k - 1].copy()
                a[:, k - 1] = a[:, k - 2].copy()
                a[:, k - 2] = v

                v_m = mapping[:, k - 1].copy()
                mapping[:, k - 1] = mapping[:, k - 2].copy()
                mapping[:, k - 2] = v_m

                # Update the Gram-Schmidt coefficients
                for s in range(k - 1, k + 1):
                    u[s - 1, 0 : (s - 1)] = dot(a[:, s - 1].T, b[:, 0 : (s - 1)]) / m[0 : (s - 1)]
                    b[:, s - 1] = a[:, s - 1] - dot(b[:, 0 : (s - 1)], u[s - 1, 0 : (s - 1)].T)
                    m[s - 1] = dot(b[:, s - 1], b[:, s - 1])

                if k > 2:
                    k -= 1
                else:
                    # We have to do p/q, so do lstsq(q.T, p.T).T instead.
                    p = dot(a[:, k:3].T, b[:, (k - 2) : k])
                    q = np.diag(m[(k - 2) : k])  # type: ignore
                    # pylint: disable=E1101
                    result = np.linalg.lstsq(q.T, p.T, rcond=None)[0].T  # type: ignore
                    u[k:3, (k - 2) : k] = result

        return a.T, mapping.T

    def get_lll_frac_coords(self, frac_coords: ArrayLike) -> np.ndarray:
        """
        Given fractional coordinates in the lattice basis, returns corresponding
        fractional coordinates in the lll basis.
        """
        return dot(frac_coords, self.lll_inverse)

    def get_frac_coords_from_lll(self, lll_frac_coords: ArrayLike) -> np.ndarray:
        """
        Given fractional coordinates in the lll basis, returns corresponding
        fractional coordinates in the lattice basis.
        """
        return dot(lll_frac_coords, self.lll_mapping)

    def get_niggli_reduced_lattice(self, tol: float = 1e-5) -> Lattice:
        """
        Get the Niggli reduced lattice using the numerically stable algo
        proposed by R. W. Grosse-Kunstleve, N. K. Sauter, & P. D. Adams,
        Acta Crystallographica Section A Foundations of Crystallography, 2003,
        60(1), 1-6. doi:10.1107/S010876730302186X

        Args:
            tol (float): The numerical tolerance. The default of 1e-5 should
                result in stable behavior for most cases.

        Returns:
            Niggli-reduced lattice.
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
            (A, B, C, E, N, Y) = (
                G[0, 0],
                G[1, 1],
                G[2, 2],
                2 * G[1, 2],
                2 * G[0, 2],
                2 * G[0, 1],
            )

            if A > B + e or (abs(A - B) < e and abs(E) > abs(N) + e):
                # A1
                M = [[0, -1, 0], [-1, 0, 0], [0, 0, -1]]
                G = dot(transpose(M), dot(G, M))
            if (B > C + e) or (abs(B - C) < e and abs(N) > abs(Y) + e):
                # A2
                M = [[-1, 0, 0], [0, 0, -1], [0, -1, 0]]
                G = dot(transpose(M), dot(G, M))
                continue

            l = 0 if abs(E) < e else E / abs(E)
            m = 0 if abs(N) < e else N / abs(N)
            n = 0 if abs(Y) < e else Y / abs(Y)
            if l * m * n == 1:
                # A3
                i = -1 if l == -1 else 1
                j = -1 if m == -1 else 1
                k = -1 if n == -1 else 1
                M = [[i, 0, 0], [0, j, 0], [0, 0, k]]
                G = dot(transpose(M), dot(G, M))
            elif l * m * n == 0 or l * m * n == -1:
                # A4
                i = -1 if l == 1 else 1
                j = -1 if m == 1 else 1
                k = -1 if n == 1 else 1

                if i * j * k == -1:
                    if n == 0:
                        k = -1
                    elif m == 0:
                        j = -1
                    elif l == 0:
                        i = -1
                M = [[i, 0, 0], [0, j, 0], [0, 0, k]]
                G = dot(transpose(M), dot(G, M))

            (A, B, C, E, N, Y) = (
                G[0, 0],
                G[1, 1],
                G[2, 2],
                2 * G[1, 2],
                2 * G[0, 2],
                2 * G[0, 1],
            )

            # A5
            if abs(E) > B + e or (abs(E - B) < e and 2 * N < Y - e) or (abs(E + B) < e and Y < -e):
                M = [[1, 0, 0], [0, 1, -E / abs(E)], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
                continue

            # A6
            if abs(N) > A + e or (abs(A - N) < e and 2 * E < Y - e) or (abs(A + N) < e and Y < -e):
                M = [[1, 0, -N / abs(N)], [0, 1, 0], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
                continue

            # A7
            if abs(Y) > A + e or (abs(A - Y) < e and 2 * E < N - e) or (abs(A + Y) < e and N < -e):
                M = [[1, -Y / abs(Y), 0], [0, 1, 0], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
                continue

            # A8
            if E + N + Y + A + B < -e or (abs(E + N + Y + A + B) < e < Y + (A + N) * 2):
                M = [[1, 0, 1], [0, 1, 1], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
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

        latt = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        mapped = self.find_mapping(latt, e, skip_rotation_matrix=True)
        if mapped is not None:
            if np.linalg.det(mapped[0].matrix) > 0:
                return mapped[0]
            return Lattice(-mapped[0].matrix)

        raise ValueError("can't find niggli")

    def scale(self, new_volume: float) -> Lattice:
        """
        Return a new Lattice with volume new_volume by performing a
        scaling of the lattice vectors so that length proportions and angles
        are preserved.

        Args:
            new_volume:
                New volume to scale to.

        Returns:
            New lattice with desired volume.
        """
        versors = self.matrix / self.abc

        geo_factor = abs(dot(np.cross(versors[0], versors[1]), versors[2]))

        ratios = np.array(self.abc) / self.c

        new_c = (new_volume / (geo_factor * np.prod(ratios))) ** (1 / 3.0)

        return Lattice(versors * (new_c * ratios), pbc=self.pbc)

    def get_wigner_seitz_cell(self) -> list[list[np.ndarray]]:
        """
        Returns the Wigner-Seitz cell for the given lattice.

        Returns:
            A list of list of coordinates.
            Each element in the list is a "facet" of the boundary of the
            Wigner Seitz cell. For instance, a list of four coordinates will
            represent a square facet.
        """
        vec1 = self._matrix[0]
        vec2 = self._matrix[1]
        vec3 = self._matrix[2]

        list_k_points = []
        for i, j, k in itertools.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
            list_k_points.append(i * vec1 + j * vec2 + k * vec3)
        # pylint: disable=C0415
        from scipy.spatial import Voronoi

        tess = Voronoi(list_k_points)
        to_return = []
        for r in tess.ridge_dict:
            if r[0] == 13 or r[1] == 13:
                to_return.append([tess.vertices[i] for i in tess.ridge_dict[r]])

        return to_return

    def get_brillouin_zone(self) -> list[list[np.ndarray]]:
        """
        Returns the Wigner-Seitz cell for the reciprocal lattice, aka the
        Brillouin Zone.

        Returns:
            A list of list of coordinates.
            Each element in the list is a "facet" of the boundary of the
            Brillouin Zone. For instance, a list of four coordinates will
            represent a square facet.
        """
        return self.reciprocal_lattice.get_wigner_seitz_cell()

    def dot(self, coords_a: ArrayLike, coords_b: ArrayLike, frac_coords: bool = False) -> np.ndarray:
        """
        Compute the scalar product of vector(s).

        Args:
            coords_a, coords_b: Array-like objects with the coordinates.
            frac_coords (bool): Boolean stating whether the vector
                corresponds to fractional or Cartesian coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        coords_a, coords_b = (
            np.reshape(coords_a, (-1, 3)),
            np.reshape(coords_b, (-1, 3)),
        )

        if len(coords_a) != len(coords_b):
            raise ValueError("")

        if np.iscomplexobj(coords_a) or np.iscomplexobj(coords_b):
            raise TypeError("Complex array!")

        if not frac_coords:
            cart_a, cart_b = coords_a, coords_b
        else:
            cart_a = np.reshape([self.get_cartesian_coords(vec) for vec in coords_a], (-1, 3))
            cart_b = np.reshape([self.get_cartesian_coords(vec) for vec in coords_b], (-1, 3))

        return np.array([dot(a, b) for a, b in zip(cart_a, cart_b)])

    def norm(self, coords: ArrayLike, frac_coords: bool = True) -> np.ndarray:
        """
        Compute the norm of vector(s).

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
        zip_results=True,
    ) -> list[tuple[np.ndarray, float, int, np.ndarray]] | list[np.ndarray] | list:
        """
        Find all points within a sphere from the point taking into account
        periodic boundary conditions. This includes sites in other periodic
        images.

        Algorithm:

        1. place sphere of radius r in crystal and determine minimum supercell
           (parallelpiped) which would contain a sphere of radius r. for this
           we need the projection of a_1 on a unit vector perpendicular
           to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
           determine how many a_1"s it will take to contain the sphere.

           Nxmax = r * length_of_b_1 / (2 Pi)

        2. keep points falling within r.

        Args:
            frac_points: All points in the lattice in fractional coordinates.
            center: Cartesian coordinates of center of sphere.
            r: radius of sphere.
            zip_results (bool): Whether to zip the results together to group by
                 point, or return the raw fcoord, dist, index arrays

        Returns:
            if zip_results:
                [(fcoord, dist, index, supercell_image) ...] since most of the time, subsequent
                processing requires the distance, index number of the atom, or index of the image
            else:
                fcoords, dists, inds, image
        """
        try:
            # pylint: disable=C0415
            from pymatgen.optimization.neighbors import find_points_in_spheres
        except ImportError:
            return self.get_points_in_sphere_py(frac_points=frac_points, center=center, r=r, zip_results=zip_results)
        else:
            frac_points = np.ascontiguousarray(frac_points, dtype=np.float_)
            r = float(r)
            lattice_matrix = np.array(self.matrix)
            lattice_matrix = np.ascontiguousarray(lattice_matrix)
            cart_coords = self.get_cartesian_coords(frac_points)
            _, indices, images, distances = find_points_in_spheres(
                all_coords=cart_coords,
                center_coords=np.ascontiguousarray([center], dtype=float),
                r=r,
                pbc=np.array(self.pbc, dtype=int),
                lattice=lattice_matrix,
                tol=1e-8,
            )
            if len(indices) < 1:
                return [] if zip_results else [()] * 4
            fcoords = frac_points[indices] + images
            if zip_results:
                return list(
                    zip(
                        fcoords,
                        distances,
                        indices,
                        images,
                    )
                )
            return [
                fcoords,
                distances,
                indices,
                images,
            ]

    def get_points_in_sphere_py(
        self,
        frac_points: ArrayLike,
        center: ArrayLike,
        r: float,
        zip_results=True,
    ) -> list[tuple[np.ndarray, float, int, np.ndarray]] | list[np.ndarray]:
        """
        Find all points within a sphere from the point taking into account
        periodic boundary conditions. This includes sites in other periodic
        images.

        Algorithm:

        1. place sphere of radius r in crystal and determine minimum supercell
           (parallelpiped) which would contain a sphere of radius r. for this
           we need the projection of a_1 on a unit vector perpendicular
           to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
           determine how many a_1"s it will take to contain the sphere.

           Nxmax = r * length_of_b_1 / (2 Pi)

        2. keep points falling within r.

        Args:
            frac_points: All points in the lattice in fractional coordinates.
            center: Cartesian coordinates of center of sphere.
            r: radius of sphere.
            zip_results (bool): Whether to zip the results together to group by
                 point, or return the raw fcoord, dist, index arrays

        Returns:
            if zip_results:
                [(fcoord, dist, index, supercell_image) ...] since most of the time, subsequent
                processing requires the distance, index number of the atom, or index of the image
            else:
                fcoords, dists, inds, image
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
            return [] if zip_results else [()] * 4  # type: ignore
        if zip_results:
            return neighbors
        return [np.array(i) for i in list(zip(*neighbors))]

    @deprecated(get_points_in_sphere, "This is retained purely for checking purposes.")
    def get_points_in_sphere_old(
        self,
        frac_points: ArrayLike,
        center: ArrayLike,
        r: float,
        zip_results=True,
    ) -> (
        list[tuple[np.ndarray, float, int, np.ndarray]]
        | tuple[list[np.ndarray], list[float], list[int], list[np.ndarray]]
    ):
        """
        Find all points within a sphere from the point taking into account
        periodic boundary conditions. This includes sites in other periodic
        images. Does not support partial periodic boundary conditions.

        Algorithm:

        1. place sphere of radius r in crystal and determine minimum supercell
           (parallelpiped) which would contain a sphere of radius r. for this
           we need the projection of a_1 on a unit vector perpendicular
           to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
           determine how many a_1"s it will take to contain the sphere.

           Nxmax = r * length_of_b_1 / (2 Pi)

        2. keep points falling within r.

        Args:
            frac_points: All points in the lattice in fractional coordinates.
            center: Cartesian coordinates of center of sphere.
            r: radius of sphere.
            zip_results (bool): Whether to zip the results together to group by
                 point, or return the raw fcoord, dist, index arrays

        Returns:
            if zip_results:
                [(fcoord, dist, index, supercell_image) ...] since most of the time, subsequent
                processing requires the distance, index number of the atom, or index of the image
            else:
                fcoords, dists, inds, image
        """
        if self.pbc != (True, True, True):
            raise RuntimeError("get_points_in_sphere_old does not support partial periodic boundary conditions")
        # TODO: refactor to use lll matrix (nmax will be smaller)
        # Determine the maximum number of supercells in each direction
        #  required to contain a sphere of radius n
        recp_len = np.array(self.reciprocal_lattice.abc) / (2 * pi)
        nmax = float(r) * recp_len + 0.01

        # Get the fractional coordinates of the center of the sphere
        pcoords = self.get_fractional_coords(center)
        center = np.array(center)

        # Prepare the list of output atoms
        n = len(frac_points)  # type: ignore
        fcoords = np.array(frac_points) % 1
        indices = np.arange(n)

        # Generate all possible images that could be within `r` of `center`
        mins = np.floor(pcoords - nmax)
        maxes = np.ceil(pcoords + nmax)
        arange = np.arange(start=mins[0], stop=maxes[0], dtype=int)
        brange = np.arange(start=mins[1], stop=maxes[1], dtype=int)
        crange = np.arange(start=mins[2], stop=maxes[2], dtype=int)
        arange = arange[:, None] * np.array([1, 0, 0], dtype=int)[None, :]
        brange = brange[:, None] * np.array([0, 1, 0], dtype=int)[None, :]
        crange = crange[:, None] * np.array([0, 0, 1], dtype=int)[None, :]
        images = arange[:, None, None] + brange[None, :, None] + crange[None, None, :]

        # Generate the coordinates of all atoms within these images
        shifted_coords = fcoords[:, None, None, None, :] + images[None, :, :, :, :]

        # Determine distance from `center`
        cart_coords = self.get_cartesian_coords(fcoords)
        cart_images = self.get_cartesian_coords(images)
        coords = cart_coords[:, None, None, None, :] + cart_images[None, :, :, :, :]  # pylint: disable=E1126
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
                )
            )
        return (  # type: ignore
            shifted_coords[within_r],
            np.sqrt(d_2[within_r]),
            indices[within_r[0]],
            images[within_r[1:]],
        )

    def get_all_distances(
        self,
        fcoords1: ArrayLike,
        fcoords2: ArrayLike,
    ) -> np.ndarray:
        """
        Returns the distances between two lists of coordinates taking into
        account periodic boundary conditions and the lattice. Note that this
        computes an MxN array of distances (i.e. the distance between each
        point in fcoords1 and every coordinate in fcoords2). This is
        different functionality from pbc_diff.

        Args:
            fcoords1: First set of fractional coordinates. e.g., [0.5, 0.6,
                0.7] or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
                coord or any array of coords.
            fcoords2: Second set of fractional coordinates.

        Returns:
            2d array of Cartesian distances. E.g the distance between
            fcoords1[i] and fcoords2[j] is distances[i,j]
        """
        v, d2 = pbc_shortest_vectors(self, fcoords1, fcoords2, return_d2=True)
        return np.sqrt(d2)

    def is_hexagonal(self, hex_angle_tol: float = 5, hex_length_tol: float = 0.01) -> bool:
        """
        :param hex_angle_tol: Angle tolerance
        :param hex_length_tol: Length tolerance
        :return: Whether lattice corresponds to hexagonal lattice.
        """
        lengths = self.lengths
        angles = self.angles
        right_angles = [i for i in range(3) if abs(angles[i] - 90) < hex_angle_tol]
        hex_angles = [
            i for i in range(3) if abs(angles[i] - 60) < hex_angle_tol or abs(angles[i] - 120) < hex_angle_tol
        ]

        return (
            len(right_angles) == 2
            and len(hex_angles) == 1
            and abs(lengths[right_angles[0]] - lengths[right_angles[1]]) < hex_length_tol
        )

    def get_distance_and_image(
        self,
        frac_coords1: ArrayLike,
        frac_coords2: ArrayLike,
        jimage: ArrayLike | None = None,
    ) -> tuple[float, np.ndarray]:
        """
        Gets distance between two frac_coords assuming periodic boundary
        conditions. If the index jimage is not specified it selects the j
        image nearest to the i atom and returns the distance and jimage
        indices in terms of lattice vector translations. If the index jimage
        is specified it returns the distance between the frac_coords1 and
        the specified jimage of frac_coords2, and the given jimage is also
        returned.

        Args:
            frac_coords1 (3x1 array): Reference fcoords to get distance from.
            frac_coords2 (3x1 array): fcoords to get distance from.
            jimage (3x1 array): Specific periodic image in terms of
                lattice translations, e.g., [1,0,0] implies to take periodic
                image that is one a-lattice vector away. If jimage is None,
                the image that is nearest to the site is found.

        Returns:
            (distance, jimage): distance and periodic lattice translations
            of the other site for which the distance applies. This means that
            the distance between frac_coords1 and (jimage + frac_coords2) is
            equal to distance.
        """
        if jimage is None:
            v, d2 = pbc_shortest_vectors(self, frac_coords1, frac_coords2, return_d2=True)
            fc = self.get_fractional_coords(v[0][0]) + frac_coords1 - frac_coords2  # type: ignore
            fc = np.array(np.round(fc), dtype=int)
            return np.sqrt(d2[0, 0]), fc

        jimage = np.array(jimage)
        mapped_vec = self.get_cartesian_coords(jimage + frac_coords2 - frac_coords1)  # type: ignore
        return np.linalg.norm(mapped_vec), jimage  # type: ignore

    def get_miller_index_from_coords(
        self,
        coords: ArrayLike,
        coords_are_cartesian: bool = True,
        round_dp: int = 4,
        verbose: bool = True,
    ) -> tuple[int, int, int]:
        """
        Get the Miller index of a plane from a list of site coordinates.

        A minimum of 3 sets of coordinates are required. If more than 3 sets of
        coordinates are given, the best plane that minimises the distance to all
        points will be calculated.

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
            (tuple): The Miller index.
        """
        if coords_are_cartesian:
            coords = [self.get_fractional_coords(c) for c in coords]  # type: ignore

        coords = np.asarray(coords)
        g = coords.sum(axis=0) / coords.shape[0]

        # run singular value decomposition
        _, _, vh = np.linalg.svd(coords - g)

        # get unitary normal vector
        u_norm = vh[2, :]
        return get_integer_index(u_norm, round_dp=round_dp, verbose=verbose)

    def get_recp_symmetry_operation(self, symprec: float = 0.01) -> list:
        """
        Find the symmetric operations of the reciprocal lattice,
        to be used for hkl transformations
        Args:
            symprec: default is 0.001
        """
        recp_lattice = self.reciprocal_lattice_crystallographic
        # get symmetry operations from input conventional unit cell
        # Need to make sure recp lattice is big enough, otherwise symmetry
        # determination will fail. We set the overall volume to 1.
        recp_lattice = recp_lattice.scale(1)
        # need a localized import of structure to build a
        # pseudo empty lattice for SpacegroupAnalyzer
        # pylint: disable=C0415
        from pymatgen.core.structure import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        recp = Structure(recp_lattice, ["H"], [[0, 0, 0]])
        # Creates a function that uses the symmetry operations in the
        # structure to find Miller indices that might give repetitive slabs
        analyzer = SpacegroupAnalyzer(recp, symprec=symprec)
        recp_symmops = analyzer.get_symmetry_operations()

        return recp_symmops


def get_integer_index(miller_index: Sequence[float], round_dp: int = 4, verbose: bool = True) -> tuple[int, int, int]:
    """
    Attempt to convert a vector of floats to whole numbers.

    Args:
        miller_index (list of float): A list miller indexes.
        round_dp (int, optional): The number of decimal places to round the
            miller index to.
        verbose (bool, optional): Whether to print warnings.

    Returns:
        (tuple): The Miller index.
    """
    mi = np.asarray(miller_index)
    # deal with the case we have small irregular floats
    # that are all equal or factors of each other
    mi /= min(m for m in mi if m != 0)
    mi /= np.max(np.abs(mi))

    # deal with the case we have nice fractions
    md = [Fraction(n).limit_denominator(12).denominator for n in mi]
    mi *= reduce(lambda x, y: x * y, md)
    int_miller_index = np.round(mi, 1).astype(int)
    mi /= np.abs(reduce(math.gcd, int_miller_index))

    # round to a reasonable precision
    mi = np.array([round(h, round_dp) for h in mi])

    # need to recalculate this after rounding as values may have changed
    int_miller_index = np.round(mi, 1).astype(int)
    if np.any(np.abs(mi - int_miller_index) > 1e-6) and verbose:
        warnings.warn("Non-integer encountered in Miller index")
    else:
        mi = int_miller_index

    # minimise the number of negative indexes
    mi += 0  # converts -0 to 0

    def n_minus(index):
        return len([h for h in index if h < 0])

    if n_minus(mi) > n_minus(mi * -1):
        mi *= -1

    # if only one index is negative, make sure it is the smallest
    # e.g. (-2 1 0) -> (2 -1 0)
    if sum(mi != 0) == 2 and n_minus(mi) == 1 and abs(min(mi)) > max(mi):
        mi *= -1

    return tuple(mi)  # type: ignore


def get_points_in_spheres(
    all_coords: np.ndarray,
    center_coords: np.ndarray,
    r: float,
    pbc: bool | list[bool] | tuple[bool, bool, bool] = True,
    numerical_tol: float = 1e-8,
    lattice: Lattice | None = None,
    return_fcoords: bool = False,
) -> list[list[tuple[np.ndarray, float, int, np.ndarray]]]:
    """
    For each point in `center_coords`, get all the neighboring points in `all_coords` that are within the
    cutoff radius `r`.

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
    pbc = np.array(pbc, dtype=bool)  # type: ignore
    if return_fcoords and lattice is None:
        raise ValueError("Lattice needs to be supplied to compute fractional coordinates")
    center_coords_min = np.min(center_coords, axis=0)
    center_coords_max = np.max(center_coords, axis=0)
    # The lower bound of all considered atom coords
    global_min = center_coords_min - r - numerical_tol
    global_max = center_coords_max + r + numerical_tol
    if np.any(pbc):
        if lattice is None:
            raise ValueError("Lattice needs to be supplied when considering periodic boundary")
        recp_len = np.array(lattice.reciprocal_lattice.abc)
        maxr = np.ceil((r + 0.15) * recp_len / (2 * math.pi))
        frac_coords = lattice.get_fractional_coords(center_coords)
        nmin_temp = np.floor(np.min(frac_coords, axis=0)) - maxr
        nmax_temp = np.ceil(np.max(frac_coords, axis=0)) + maxr
        nmin = np.zeros_like(nmin_temp)
        nmin[pbc] = nmin_temp[pbc]
        nmax = np.ones_like(nmax_temp)
        nmax[pbc] = nmax_temp[pbc]
        all_ranges = [np.arange(x, y, dtype="int64") for x, y in zip(nmin, nmax)]
        matrix = lattice.matrix
        # temporarily hold the fractional coordinates
        image_offsets = lattice.get_fractional_coords(all_coords)
        all_fcoords = []
        # only wrap periodic boundary
        for k in range(3):
            if pbc[k]:  # type: ignore
                all_fcoords.append(np.mod(image_offsets[:, k : k + 1], 1))
            else:
                all_fcoords.append(image_offsets[:, k : k + 1])
        all_fcoords = np.concatenate(all_fcoords, axis=1)
        image_offsets = image_offsets - all_fcoords
        coords_in_cell = np.dot(all_fcoords, matrix)
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
        if len(valid_coords) < 1:
            return [[]] * len(center_coords)
        valid_coords = np.concatenate(valid_coords, axis=0)
        valid_images = np.concatenate(valid_images, axis=0)

    else:
        valid_coords = all_coords  # type: ignore
        valid_images = [[0, 0, 0]] * len(valid_coords)
        valid_indices = np.arange(len(valid_coords))  # type: ignore

    # Divide the valid 3D space into cubes and compute the cube ids
    all_cube_index = _compute_cube_index(valid_coords, global_min, r)  # type: ignore
    nx, ny, nz = _compute_cube_index(global_max, global_min, r) + 1
    all_cube_index = _three_to_one(all_cube_index, ny, nz)
    site_cube_index = _three_to_one(_compute_cube_index(center_coords, global_min, r), ny, nz)
    # create cube index to coordinates, images, and indices map
    cube_to_coords: dict[int, list] = collections.defaultdict(list)
    cube_to_images: dict[int, list] = collections.defaultdict(list)
    cube_to_indices: dict[int, list] = collections.defaultdict(list)
    for i, j, k, l in zip(all_cube_index.ravel(), valid_coords, valid_images, valid_indices):
        cube_to_coords[i].append(j)
        cube_to_images[i].append(k)
        cube_to_indices[i].append(l)

    # find all neighboring cubes for each atom in the lattice cell
    site_neighbors = find_neighbors(site_cube_index, nx, ny, nz)
    neighbors: list[list[tuple[np.ndarray, float, int, np.ndarray]]] = []

    for i, j in zip(center_coords, site_neighbors):
        l1 = np.array(_three_to_one(j, ny, nz), dtype=int).ravel()
        # use the cube index map to find the all the neighboring
        # coords, images, and indices
        ks = [k for k in l1 if k in cube_to_coords]
        if not ks:
            neighbors.append([])
            continue
        nn_coords = np.concatenate([cube_to_coords[k] for k in ks], axis=0)
        nn_images = itertools.chain(*(cube_to_images[k] for k in ks))
        nn_indices = itertools.chain(*(cube_to_indices[k] for k in ks))
        dist = np.linalg.norm(nn_coords - i[None, :], axis=1)
        nns: list[tuple[np.ndarray, float, int, np.ndarray]] = []
        for coord, index, image, d in zip(nn_coords, nn_indices, nn_images, dist):
            # filtering out all sites that are beyond the cutoff
            # Here there is no filtering of overlapping sites
            if d < r + numerical_tol:
                if return_fcoords and (lattice is not None):
                    coord = np.round(lattice.get_fractional_coords(coord), 10)
                nn = (coord, float(d), int(index), image)
                nns.append(nn)
        neighbors.append(nns)
    return neighbors


# The following internal methods are used in the get_points_in_sphere method.
def _compute_cube_index(coords: np.ndarray, global_min: float, radius: float) -> np.ndarray:
    """
    Compute the cube index from coordinates
    Args:
        coords: (nx3 array) atom coordinates
        global_min: (float) lower boundary of coordinates
        radius: (float) cutoff radius

    Returns: (nx3 array) int indices

    """
    return np.array(np.floor((coords - global_min) / radius), dtype=int)


def _one_to_three(label1d: np.ndarray, ny: int, nz: int) -> np.ndarray:
    """
    Convert a 1D index array to 3D index array

    Args:
        label1d: (array) 1D index array
        ny: (int) number of cells in y direction
        nz: (int) number of cells in z direction

    Returns: (nx3) int array of index

    """
    last = np.mod(label1d, nz)
    second = np.mod((label1d - last) / nz, ny)
    first = (label1d - last - second * nz) / (ny * nz)
    return np.concatenate([first, second, last], axis=1)


def _three_to_one(label3d: np.ndarray, ny: int, nz: int) -> np.ndarray:
    """
    The reverse of _one_to_three
    """
    return np.array(label3d[:, 0] * ny * nz + label3d[:, 1] * nz + label3d[:, 2]).reshape((-1, 1))


def find_neighbors(label: np.ndarray, nx: int, ny: int, nz: int) -> list[np.ndarray]:
    """
    Given a cube index, find the neighbor cube indices.

    Args:
        label: (array) (n,) or (n x 3) indice array
        nx: (int) number of cells in y direction
        ny: (int) number of cells in y direction
        nz: (int) number of cells in z direction

    Returns:
        Neighbor cell indices.
    """
    array = [[-1, 0, 1]] * 3
    neighbor_vectors = np.array(list(itertools.product(*array)), dtype=int)
    if np.shape(label)[1] == 1:
        label3d = _one_to_three(label, ny, nz)
    else:
        label3d = label
    all_labels = label3d[:, None, :] - neighbor_vectors[None, :, :]
    filtered_labels = []
    # filter out out-of-bound labels i.e., label < 0
    for labels in all_labels:
        ind = (labels[:, 0] < nx) * (labels[:, 1] < ny) * (labels[:, 2] < nz) * np.all(labels > -1e-5, axis=1)
        filtered_labels.append(labels[ind])
    return filtered_labels
