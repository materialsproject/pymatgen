# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import math
import itertools
import warnings

from six.moves import map, zip

import numpy as np
from numpy.linalg import inv
from numpy import pi, dot, transpose, radians

from monty.json import MSONable
from monty.dev import deprecated
from pymatgen.util.coord_utils import pbc_shortest_vectors
from pymatgen.util.num_utils import abs_cap

"""
This module defines the classes relating to 3D lattices.
"""


__author__ = "Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"



class Lattice(MSONable):
    """
    A lattice object.  Essentially a matrix with conversion matrices. In
    general, it is assumed that length units are in Angstroms and angles are in
    degrees unless otherwise stated.
    """

    # Properties lazily generated for efficiency.


    def __init__(self, matrix):
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
        """
        m = np.array(matrix, dtype=np.float64).reshape((3, 3))
        lengths = np.sqrt(np.sum(m ** 2, axis=1))
        angles = np.zeros(3)
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            angles[i] = abs_cap(dot(m[j], m[k]) / (lengths[j] * lengths[k]))

        self._angles = np.arccos(angles) * 180. / pi
        self._lengths = lengths
        self._matrix = m
        self._inv_matrix = None
        self._metric_tensor = None
        self._diags = None
        self._lll_matrix_mappings = {}
        self._lll_inverse = None
        self.is_orthogonal = all([abs(a - 90) < 1e-5 for a in self._angles])

    def __format__(self, fmt_spec=''):
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
        m = self.matrix.tolist()
        if fmt_spec.endswith("l"):
            fmt = "[[{}, {}, {}], [{}, {}, {}], [{}, {}, {}]]"
            fmt_spec = fmt_spec[:-1]
        elif fmt_spec.endswith("p"):
            fmt = "{{{}, {}, {}, {}, {}, {}}}"
            fmt_spec = fmt_spec[:-1]
            m = self.lengths_and_angles
        else:
            fmt = "{} {} {}\n{} {} {}\n{} {} {}"
        return fmt.format(*[format(c, fmt_spec) for row in m
                            for c in row])

    def copy(self):
        """Deep copy of self."""
        return self.__class__(self.matrix.copy())

    @property
    def matrix(self):
        """Copy of matrix representing the Lattice"""
        return np.copy(self._matrix)

    @property
    def inv_matrix(self):
        """
        Inverse of lattice matrix.
        """
        if self._inv_matrix is None:
            self._inv_matrix = inv(self._matrix)
        return self._inv_matrix

    @property
    def metric_tensor(self):
        """
        The metric tensor of the lattice.
        """
        if self._metric_tensor is None:
            self._metric_tensor = np.dot(self._matrix, self._matrix.T)
        return self._metric_tensor

    def get_cartesian_coords(self, fractional_coords):
        """
        Returns the cartesian coordinates given fractional coordinates.

        Args:
            fractional_coords (3x1 array): Fractional coords.

        Returns:
            Cartesian coordinates
        """
        return dot(fractional_coords, self._matrix)

    def get_fractional_coords(self, cart_coords):
        """
        Returns the fractional coordinates given cartesian coordinates.

        Args:
            cart_coords (3x1 array): Cartesian coords.

        Returns:
            Fractional coordinates.
        """
        return dot(cart_coords, self.inv_matrix)

    @staticmethod
    def cubic(a):
        """
        Convenience constructor for a cubic lattice.

        Args:
            a (float): The *a* lattice parameter of the cubic cell.

        Returns:
            Cubic lattice of dimensions a x a x a.
        """
        return Lattice([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]])

    @staticmethod
    def tetragonal(a, c):
        """
        Convenience constructor for a tetragonal lattice.

        Args:
            a (float): *a* lattice parameter of the tetragonal cell.
            c (float): *c* lattice parameter of the tetragonal cell.

        Returns:
            Tetragonal lattice of dimensions a x a x c.
        """
        return Lattice.from_parameters(a, a, c, 90, 90, 90)

    @staticmethod
    def orthorhombic(a, b, c):
        """
        Convenience constructor for an orthorhombic lattice.

        Args:
            a (float): *a* lattice parameter of the orthorhombic cell.
            b (float): *b* lattice parameter of the orthorhombic cell.
            c (float): *c* lattice parameter of the orthorhombic cell.

        Returns:
            Orthorhombic lattice of dimensions a x b x c.
        """
        return Lattice.from_parameters(a, b, c, 90, 90, 90)

    @staticmethod
    def monoclinic(a, b, c, beta):
        """
        Convenience constructor for a monoclinic lattice.

        Args:
            a (float): *a* lattice parameter of the monoclinc cell.
            b (float): *b* lattice parameter of the monoclinc cell.
            c (float): *c* lattice parameter of the monoclinc cell.
            beta (float): *beta* angle between lattice vectors b and c in
                degrees.

        Returns:
            Monoclinic lattice of dimensions a x b x c with non right-angle
            beta between lattice vectors a and c.
        """
        return Lattice.from_parameters(a, b, c, 90, beta, 90)

    @staticmethod
    def hexagonal(a, c):
        """
        Convenience constructor for a hexagonal lattice.

        Args:
            a (float): *a* lattice parameter of the hexagonal cell.
            c (float): *c* lattice parameter of the hexagonal cell.

        Returns:
            Hexagonal lattice of dimensions a x a x c.
        """
        return Lattice.from_parameters(a, a, c, 90, 90, 120)

    @staticmethod
    def rhombohedral(a, alpha):
        """
        Convenience constructor for a rhombohedral lattice.

        Args:
            a (float): *a* lattice parameter of the rhombohedral cell.
            alpha (float): Angle for the rhombohedral lattice in degrees.

        Returns:
            Rhombohedral lattice of dimensions a x a x a.
        """
        return Lattice.from_parameters(a, a, a, alpha, alpha, alpha)

    @staticmethod
    def from_lengths_and_angles(abc, ang):
        """
        Create a Lattice using unit cell lengths and angles (in degrees).

        Args:
            abc (3x1 array): Lattice parameters, e.g. (4, 4, 5).
            ang (3x1 array): Lattice angles in degrees, e.g., (90,90,120).

        Returns:
            A Lattice with the specified lattice parameters.
        """
        return Lattice.from_parameters(abc[0], abc[1], abc[2],
                                       ang[0], ang[1], ang[2])

    @staticmethod
    def from_parameters(a, b, c, alpha, beta, gamma):
        """
        Create a Lattice using unit cell lengths and angles (in degrees).

        Args:
            a (float): *a* lattice parameter.
            b (float): *b* lattice parameter.
            c (float): *c* lattice parameter.
            alpha (float): *alpha* angle in degrees.
            beta (float): *beta* angle in degrees.
            gamma (float): *gamma* angle in degrees.

        Returns:
            Lattice with the specified lattice parameters.
        """

        alpha_r = radians(alpha)
        beta_r = radians(beta)
        gamma_r = radians(gamma)
        val = (np.cos(alpha_r) * np.cos(beta_r) - np.cos(gamma_r))\
            / (np.sin(alpha_r) * np.sin(beta_r))
        # Sometimes rounding errors result in values slightly > 1.
        val = abs_cap(val)
        gamma_star = np.arccos(val)
        vector_a = [a * np.sin(beta_r), 0.0, a * np.cos(beta_r)]
        vector_b = [-b * np.sin(alpha_r) * np.cos(gamma_star),
                    b * np.sin(alpha_r) * np.sin(gamma_star),
                    b * np.cos(alpha_r)]
        vector_c = [0.0, 0.0, float(c)]
        return Lattice([vector_a, vector_b, vector_c])

    @classmethod
    def from_dict(cls, d, fmt=None, **kwargs):
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
            from pymatgen.io.abinit.abiobjects import lattice_from_abivars
            kwargs.update(d)
            return lattice_from_abivars(cls=cls, **kwargs)

        if "matrix" in d:
            return cls(d["matrix"])
        else:
            return cls.from_parameters(d["a"], d["b"], d["c"],
                                       d["alpha"], d["beta"], d["gamma"])

    @property
    def angles(self):
        """
        Returns the angles (alpha, beta, gamma) of the lattice.
        """
        return tuple(self._angles)

    @property
    def a(self):
        """
        *a* lattice parameter.ATATClusterExpansion
        """
        return self._lengths[0]

    @property
    def b(self):
        """
        *b* lattice parameter.
        """
        return self._lengths[1]

    @property
    def c(self):
        """
        *c* lattice parameter.
        """
        return self._lengths[2]

    @property
    def abc(self):
        """
        Lengths of the lattice vectors, i.e. (a, b, c)
        """
        return tuple(self._lengths)

    @property
    def alpha(self):
        """
        Angle alpha of lattice in degrees.
        """
        return self._angles[0]

    @property
    def beta(self):
        """
        Angle beta of lattice in degrees.
        """
        return self._angles[1]

    @property
    def gamma(self):
        """
        Angle gamma of lattice in degrees.
        """
        return self._angles[2]

    @property
    def volume(self):
        """
        Volume of the unit cell.
        """
        m = self._matrix
        return abs(np.dot(np.cross(m[0], m[1]), m[2]))

    @property
    def lengths_and_angles(self):
        """
        Returns (lattice lengths, lattice angles).
        """
        return tuple(self._lengths), tuple(self._angles)

    @property
    def reciprocal_lattice(self):
        """
        Return the reciprocal lattice. Note that this is the standard
        reciprocal lattice used for solid state physics with a factor of 2 *
        pi. If you are looking for the crystallographic reciprocal lattice,
        use the reciprocal_lattice_crystallographic property.
        The property is lazily generated for efficiency.
        """
        try:
            return self._reciprocal_lattice
        except AttributeError:
            v = np.linalg.inv(self._matrix).T
            self._reciprocal_lattice = Lattice(v * 2 * np.pi)
            return self._reciprocal_lattice

    @property
    def reciprocal_lattice_crystallographic(self):
        """
        Returns the *crystallographic* reciprocal lattice, i.e., no factor of
        2 * pi.
        """
        return Lattice(self.reciprocal_lattice.matrix / (2 * np.pi))

    @property
    def lll_matrix(self):
        if 0.75 not in self._lll_matrix_mappings:
            self._lll_matrix_mappings[0.75] = self._calculate_lll()
        return self._lll_matrix_mappings[0.75][0]

    @property
    def lll_mapping(self):
        if 0.75 not in self._lll_matrix_mappings:
            self._lll_matrix_mappings[0.75] = self._calculate_lll()
        return self._lll_matrix_mappings[0.75][1]

    @property
    def lll_inverse(self):
        if self._lll_inverse is not None:
            return self._lll_inverse
        else:
            self._lll_inverse = np.linalg.inv(self.lll_mapping)
            return self._lll_inverse

    def __repr__(self):
        outs = ["Lattice", "    abc : " + " ".join(map(repr, self._lengths)),
                " angles : " + " ".join(map(repr, self._angles)),
                " volume : " + repr(self.volume),
                "      A : " + " ".join(map(repr, self._matrix[0])),
                "      B : " + " ".join(map(repr, self._matrix[1])),
                "      C : " + " ".join(map(repr, self._matrix[2]))]
        return "\n".join(outs)

    def __eq__(self, other):
        """
        A lattice is considered to be equal to another if the internal matrix
        representation satisfies np.allclose(matrix1, matrix2) to be True.
        """
        if other is None:
            return False
        # shortcut the np.allclose if the memory addresses are the same
        # (very common in Structure.from_sites)
        return self is other or np.allclose(self.matrix, other.matrix)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return 7

    def __str__(self):
        return "\n".join([" ".join(["%.6f" % i for i in row])
                          for row in self._matrix])

    def as_dict(self, verbosity=0):
        """""
        Json-serialization dict representation of the Lattice.

        Args:
            verbosity (int): Verbosity level. Default of 0 only includes the
                matrix representation. Set to 1 for more details.
        """

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "matrix": self._matrix.tolist()}
        if verbosity > 0:
            d.update({
                "a": float(self.a),
                "b": float(self.b),
                "c": float(self.c),
                "alpha": float(self.alpha),
                "beta": float(self.beta),
                "gamma": float(self.gamma),
                "volume": float(self.volume)
            })

        return d

    def find_all_mappings(self, other_lattice, ltol=1e-5, atol=1,
                          skip_rotation_matrix=False):
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
        (lengths, angles) = other_lattice.lengths_and_angles
        (alpha, beta, gamma) = angles

        frac, dist, _ = self.get_points_in_sphere([[0, 0, 0]], [0, 0, 0],
                                                  max(lengths) * (1 + ltol),
                                                  zip_results=False)
        cart = self.get_cartesian_coords(frac)
        # this can't be broadcast because they're different lengths
        inds = [np.logical_and(dist / l < 1 + ltol,
                               dist / l > 1 / (1 + ltol)) for l in lengths]
        c_a, c_b, c_c = (cart[i] for i in inds)
        f_a, f_b, f_c = (frac[i] for i in inds)
        l_a, l_b, l_c = (np.sum(c ** 2, axis=-1) ** 0.5 for c in (c_a, c_b, c_c))

        def get_angles(v1, v2, l1, l2):
            x = np.inner(v1, v2) / l1[:, None] / l2
            x[x > 1] = 1
            x[x < -1] = -1
            angles = np.arccos(x) * 180. / pi
            return angles

        alphab = np.abs(get_angles(c_b, c_c, l_b, l_c) - alpha) < atol
        betab = np.abs(get_angles(c_a, c_c, l_a, l_c) - beta) < atol
        gammab = np.abs(get_angles(c_a, c_b, l_a, l_b) - gamma) < atol

        for i, all_j in enumerate(gammab):
            inds = np.logical_and(all_j[:, None],
                                  np.logical_and(alphab,
                                                 betab[i][None, :]))
            for j, k in np.argwhere(inds):
                scale_m = np.array((f_a[i], f_b[j], f_c[k]), dtype=np.int)
                if abs(np.linalg.det(scale_m)) < 1e-8:
                    continue

                aligned_m = np.array((c_a[i], c_b[j], c_c[k]))

                if skip_rotation_matrix:
                    rotation_m = None
                else:
                    rotation_m = np.linalg.solve(aligned_m,
                                                 other_lattice.matrix)

                yield Lattice(aligned_m), rotation_m, scale_m

    def find_mapping(self, other_lattice, ltol=1e-5, atol=1,
                     skip_rotation_matrix=False):
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
        for x in self.find_all_mappings(
                other_lattice, ltol, atol,
                skip_rotation_matrix=skip_rotation_matrix):
            return x

    def get_lll_reduced_lattice(self, delta=0.75):
        if delta not in self._lll_matrix_mappings:
            self._lll_matrix_mappings[delta] = self._calculate_lll()
        return Lattice(self._lll_matrix_mappings[delta][0])

    def _calculate_lll(self, delta=0.75):
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
        a = self._matrix.copy().T

        b = np.zeros((3, 3))  # Vectors after the Gram-Schmidt process
        u = np.zeros((3, 3))  # Gram-Schmidt coeffieicnts
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
                    mapping[:, k - 1] = mapping[:, k - 1] - q * \
                        mapping[:, i - 1]
                    uu = list(u[i - 1, 0:(i - 1)])
                    uu.append(1)
                    # Update the GS coefficients.
                    u[k - 1, 0:i] = u[k - 1, 0:i] - q * np.array(uu)

            # Check the Lovasz condition.
            if dot(b[:, k - 1], b[:, k - 1]) >=\
                    (delta - abs(u[k - 1, k - 2]) ** 2) *\
                    dot(b[:, (k - 2)], b[:, (k - 2)]):
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
                    u[s - 1, 0:(s - 1)] = dot(a[:, s - 1].T,
                                              b[:, 0:(s - 1)]) / m[0:(s - 1)]
                    b[:, s - 1] = a[:, s - 1] - dot(b[:, 0:(s - 1)],
                                                    u[s - 1, 0:(s - 1)].T)
                    m[s - 1] = dot(b[:, s - 1], b[:, s - 1])

                if k > 2:
                    k -= 1
                else:
                    # We have to do p/q, so do lstsq(q.T, p.T).T instead.
                    p = dot(a[:, k:3].T, b[:, (k - 2):k])
                    q = np.diag(m[(k - 2):k])
                    result = np.linalg.lstsq(q.T, p.T)[0].T
                    u[k:3, (k - 2):k] = result

        return a.T, mapping.T

    def get_lll_frac_coords(self, frac_coords):
        """
        Given fractional coordinates in the lattice basis, returns corresponding
        fractional coordinates in the lll basis.
        """
        return np.dot(frac_coords, self.lll_inverse)

    def get_frac_coords_from_lll(self, lll_frac_coords):
        """
        Given fractional coordinates in the lll basis, returns corresponding
        fractional coordinates in the lattice basis.
        """
        return np.dot(lll_frac_coords, self.lll_mapping)

    def get_niggli_reduced_lattice(self, tol=1e-5):
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
        a = matrix[0]
        b = matrix[1]
        c = matrix[2]
        e = tol * self.volume ** (1 / 3)

        # Define metric tensor
        G = [[dot(a, a), dot(a, b), dot(a, c)],
             [dot(a, b), dot(b, b), dot(b, c)],
             [dot(a, c), dot(b, c), dot(c, c)]]
        G = np.array(G)

        # This sets an upper limit on the number of iterations.
        for count in range(100):
            # The steps are labelled as Ax as per the labelling scheme in the
            # paper.
            (A, B, C, E, N, Y) = (G[0, 0], G[1, 1], G[2, 2],
                                  2 * G[1, 2], 2 * G[0, 2], 2 * G[0, 1])

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

            (A, B, C, E, N, Y) = (G[0, 0], G[1, 1], G[2, 2],
                                  2 * G[1, 2], 2 * G[0, 2], 2 * G[0, 1])

            # A5
            if abs(E) > B + e or (abs(E - B) < e and 2 * N < Y - e) or\
                    (abs(E + B) < e and Y < -e):
                M = [[1, 0, 0], [0, 1, -E / abs(E)], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
                continue

            # A6
            if abs(N) > A + e or (abs(A - N) < e and 2 * E < Y - e) or\
                    (abs(A + N) < e and Y < -e):
                M = [[1, 0, -N / abs(N)], [0, 1, 0], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
                continue

            # A7
            if abs(Y) > A + e or (abs(A - Y) < e and 2 * E < N - e) or\
                    (abs(A + Y) < e and N < -e):
                M = [[1, -Y / abs(Y), 0], [0, 1, 0], [0, 0, 1]]
                G = dot(transpose(M), dot(G, M))
                continue

            # A8
            if E + N + Y + A + B < -e or\
                    (abs(E + N + Y + A + B) < e < Y + (A + N) * 2):
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
            else:
                return Lattice(-mapped[0].matrix)

        raise ValueError("can't find niggli")

    def scale(self, new_volume):
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

        geo_factor = abs(np.dot(np.cross(versors[0], versors[1]), versors[2]))

        ratios = self.abc / self.c

        new_c = (new_volume / (geo_factor * np.prod(ratios))) ** (1/3.)

        return Lattice(versors * (new_c * ratios))

    def get_wigner_seitz_cell(self):
        """
        Returns the Wigner-Seitz cell for the given lattice.

        Returns:
            A list of list of coordinates.
            Each element in the list is a "facet" of the boundary of the
            Wigner Seitz cell. For instance, a list of four coordinates will
            represent a square facet.
        """
        vec1 = self.matrix[0]
        vec2 = self.matrix[1]
        vec3 = self.matrix[2]

        list_k_points = []
        for i, j, k in itertools.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
            list_k_points.append(i * vec1 + j * vec2 + k * vec3)
        from scipy.spatial import Voronoi
        tess = Voronoi(list_k_points)
        to_return = []
        for r in tess.ridge_dict:
            if r[0] == 13 or r[1] == 13:
                to_return.append([tess.vertices[i] for i in tess.ridge_dict[r]])

        return to_return

    def get_brillouin_zone(self):
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

    def dot(self, coords_a, coords_b, frac_coords=False):
        """
        Compute the scalar product of vector(s).

        Args:
            coords_a, coords_b: Array-like objects with the coordinates.
            frac_coords (bool): Boolean stating whether the vector
                corresponds to fractional or cartesian coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        coords_a, coords_b = np.reshape(coords_a, (-1, 3)), \
            np.reshape(coords_b, (-1, 3))

        if len(coords_a) != len(coords_b):
            raise ValueError("")

        if np.iscomplexobj(coords_a) or np.iscomplexobj(coords_b):
            raise TypeError("Complex array!")

        if not frac_coords:
            cart_a, cart_b = coords_a, coords_b
        else:
            cart_a = np.reshape([self.get_cartesian_coords(vec)
                                 for vec in coords_a], (-1, 3))
            cart_b = np.reshape([self.get_cartesian_coords(vec)
                                 for vec in coords_b], (-1, 3))

        return np.array([np.dot(a, b) for a, b in zip(cart_a, cart_b)])

    def norm(self, coords, frac_coords=True):
        """
        Compute the norm of vector(s).

        Args:
            coords:
                Array-like object with the coordinates.
            frac_coords:
                Boolean stating whether the vector corresponds to fractional or
                cartesian coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        return np.sqrt(self.dot(coords, coords, frac_coords=frac_coords))

    def get_points_in_sphere(self, frac_points, center, r, zip_results=True):
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
                [(fcoord, dist, index) ...] since most of the time, subsequent
                processing requires the distance.
            else:
                fcoords, dists, inds
        """
        # TODO: refactor to use lll matrix (nmax will be smaller)
        recp_len = np.array(self.reciprocal_lattice.abc) / (2 * pi)
        nmax = float(r) * recp_len + 0.01

        pcoords = self.get_fractional_coords(center)
        center = np.array(center)

        n = len(frac_points)
        fcoords = np.array(frac_points) % 1
        indices = np.arange(n)

        mins = np.floor(pcoords - nmax)
        maxes = np.ceil(pcoords + nmax)
        arange = np.arange(start=mins[0], stop=maxes[0])
        brange = np.arange(start=mins[1], stop=maxes[1])
        crange = np.arange(start=mins[2], stop=maxes[2])
        arange = arange[:, None] * np.array([1, 0, 0])[None, :]
        brange = brange[:, None] * np.array([0, 1, 0])[None, :]
        crange = crange[:, None] * np.array([0, 0, 1])[None, :]
        images = arange[:, None, None] + brange[None, :, None] +\
            crange[None, None, :]

        shifted_coords = fcoords[:, None, None, None, :] + \
            images[None, :, :, :, :]

        cart_coords = self.get_cartesian_coords(fcoords)
        cart_images = self.get_cartesian_coords(images)
        coords = cart_coords[:, None, None, None, :] + \
            cart_images[None, :, :, :, :]
        coords -= center[None, None, None, None, :]
        coords **= 2
        d_2 = np.sum(coords, axis=4)

        within_r = np.where(d_2 <= r ** 2)
        if zip_results:
            return list(zip(shifted_coords[within_r], np.sqrt(d_2[within_r]),
                            indices[within_r[0]]))
        else:
            return shifted_coords[within_r], np.sqrt(d_2[within_r]), \
                indices[within_r[0]]

    def get_all_distances(self, fcoords1, fcoords2):
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
            2d array of cartesian distances. E.g the distance between
            fcoords1[i] and fcoords2[j] is distances[i,j]
        """
        v, d2 = pbc_shortest_vectors(self, fcoords1, fcoords2, return_d2=True)
        return np.sqrt(d2)

    def is_hexagonal(self, hex_angle_tol=5, hex_length_tol=0.01):
        lengths, angles = self.lengths_and_angles
        right_angles = [i for i in range(3)
                        if abs(angles[i] - 90) < hex_angle_tol]
        hex_angles = [i for i in range(3)
                      if abs(angles[i] - 60) < hex_angle_tol or
                      abs(angles[i] - 120) < hex_angle_tol]

        return (len(right_angles) == 2 and len(hex_angles) == 1
                and abs(lengths[right_angles[0]] -
                        lengths[right_angles[1]]) < hex_length_tol)

    @deprecated(get_points_in_sphere)
    def get_all_distance_and_image(self, frac_coords1, frac_coords2):
        """
        Gets distance between two frac_coords and nearest periodic images.

        Args:
            frac_coords1 (3x1 array): Reference fcoords to get distance from.
            frac_coords2 (3x1 array): fcoords to get distance from.

        Returns:
            [(distance, jimage)] List of distance and periodic lattice
            translations of the other site for which the distance applies.
            This means that the distance between frac_coords1 and (jimage +
            frac_coords2) is equal to distance.
        """
        # The following code is heavily vectorized to maximize speed.
        # Get the image adjustment necessary to bring coords to unit_cell.
        adj1 = np.floor(frac_coords1)
        adj2 = np.floor(frac_coords2)
        # Shift coords to unitcell
        coord1 = frac_coords1 - adj1
        coord2 = frac_coords2 - adj2

        n = self._get_mic_range(coord1, coord2)
        ranges = [list(range(-i, i + 1)) for i in n]
        images = np.array(list(itertools.product(*ranges)))

        # Create tiled cartesian coords for computing distances.
        vec = np.tile(coord2 - coord1, (len(images), 1)) + images
        vec = self.get_cartesian_coords(vec)
        # Compute distances manually.
        dist = np.sqrt(np.sum(vec ** 2, 1)).tolist()
        return list(zip(dist, adj1 - adj2 + images))

    def get_distance_and_image(self, frac_coords1, frac_coords2, jimage=None):
        """
        Gets distance between two frac_coords assuming periodic boundary
        conditions. If the index jimage is not specified it selects the j
        image nearest to the i atom and returns the distance and jimage
        indices in terms of lattice vector translations. If the index jimage
        is specified it returns the distance between the frac_coords1 and
        the specified jimage of frac_coords2, and the given jimage is also
        returned.

        Args:
            fcoords1 (3x1 array): Reference fcoords to get distance from.
            fcoords2 (3x1 array): fcoords to get distance from.
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
            v, d2 = pbc_shortest_vectors(self, frac_coords1, frac_coords2,
                                         return_d2=True)
            fc = self.get_fractional_coords(v[0][0]) + frac_coords1 - \
                 frac_coords2
            fc = np.array(np.round(fc), dtype=np.int)
            return (np.sqrt(d2[0, 0]), fc)

        mapped_vec = self.get_cartesian_coords(jimage + frac_coords2
                                               - frac_coords1)
        return np.linalg.norm(mapped_vec), jimage

    @deprecated(get_points_in_sphere)
    def _get_mic_range(self, fcoords1, fcoords2):

        if self.is_orthogonal:
            return 1, 1, 1

        if self._diags is None:
            self._diags = np.sqrt((np.dot([[1, 1, 1], [-1, 1, 1], [1, -1, 1],
                                     [-1, -1, 1],
                                     ], self._matrix) ** 2).sum(1))

        vecs = np.vstack(fcoords1[:, None] - fcoords2[None, :])
        vecs = vecs - np.round(vecs)

        d = dot(vecs, self._matrix)
        d = np.sqrt(np.sum(d ** 2, axis=1))
        d = np.max(d)

        cutoff = min(d, max(self._diags) / 2)
        n = np.array(np.ceil(cutoff * np.prod(self._lengths) /
                             (self.volume * self._lengths)))

        if any(n > 50):
            n_new = np.minimum(n, 50)
            warnings.warn("Cell is highly skewed and requires a search "
                          "through image range of +- %s. For efficiency, "
                          "we will limit the search to %s images. Recommend "
                          "you do a niggli or LLL reduction of the cell "
                          "before computing distances" % (n, n_new))
            n = n_new.astype(dtype=np.int)
        return n.astype(dtype=int)