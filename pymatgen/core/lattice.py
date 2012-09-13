#!/usr/bin/env python

"""
This module defines the classes relating to 3D lattices.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import numpy as np
import numpy.linalg as npl
from numpy import pi

from pymatgen.serializers.json_coders import MSONable


class Lattice(MSONable):
    """
    A lattice object.  Essentially a matrix with conversion matrices. In
    general, it is assumed that length units are in Angstroms and angles are in
    degrees unless otherwise stated.
    """

    def __init__(self, matrix):
        """
        Create a lattice from any sequence of 9 numbers. Note that the sequence
        is assumed to be read one row at a time. Each row represents one
        lattice vector.

        Args:
            matrix:
                Sequence of numbers in any form. Examples of acceptable
                input.
                i) An actual numpy array.
                ii) [[1, 0, 0],[0, 1, 0], [0, 0, 1]]
                iii) [1, 0, 0 , 0, 1, 0, 0, 0, 1]
                iv) (1, 0, 0, 0, 1, 0, 0, 0, 1)

                Each row should correspond to a lattice vector.
                E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a lattice with
                lattice vectors [10,0,0], [20,10,0] and [0,0,30].
        """

        self._matrix = np.array(matrix).reshape((3, 3))
        #Store these matrices for faster access
        self._md2c = np.transpose(self._matrix)
        self._mc2d = npl.inv(np.transpose(self._matrix))

        lengths = np.sum(self._matrix ** 2, axis=1) ** 0.5
        angles = np.zeros((3), float)
        for i in xrange(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            angles[i] = np.dot(self._matrix[j], self._matrix[k]) / \
                (lengths[j] * lengths[k])
        angles = np.arccos(angles) * 180. / pi
        angles = np.around(angles, 9)
        self._angles = tuple(angles)
        self._lengths = tuple(lengths)

    @property
    def md2c(self):
        """Matrix for converting direct to cartesian coordinates"""
        return np.copy(self._md2c)

    @property
    def mc2d(self):
        """Matrix for converting cartesian to direct coordinates"""
        return np.copy(self._mc2d)

    @property
    def matrix(self):
        """Copy of matrix representing the Lattice"""
        return np.copy(self._matrix)

    def get_cartesian_coords(self, fractional_coords):
        """
        Returns the cartesian coordinates given fractional coordinates.

        Args:
            Fractional_coords:
                Fractional coords.

        Returns:
            Cartesian coordinates
        """
        return np.transpose(np.dot(self._md2c,
                                   np.transpose(fractional_coords)))

    def get_fractional_coords(self, cart_coords):
        """
        Returns the fractional coordinates given cartesian coordinates.

        Args:
            cartesian_coords:
                Cartesian coords.

        Returns:
            Fractional coordinates.
        """
        return np.transpose(np.dot(self._mc2d, np.transpose(cart_coords)))

    @staticmethod
    def cubic(a):
        """
        Convenience constructor for a cubic lattice.

        Args:
            a:
                The *a* lattice parameter of the cubic cell.

        Returns:
            Cubic lattice of dimensions a x a x a.
        """
        return Lattice(np.array([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]]))

    @staticmethod
    def tetragonal(a, c):
        """
        Convenience constructor for a tetragonal lattice.

        Args:
            a:
                The *a* lattice parameter of the tetragonal cell.
            c:
                The *c* lattice parameter of the tetragonal cell.

        Returns:
            Tetragonal lattice of dimensions a x a x c.
        """
        return Lattice.from_parameters(a, a, c, 90, 90, 90)

    @staticmethod
    def orthorhombic(a, b, c):
        """
        Convenience constructor for an orthorhombic lattice.

        Args:
            a:
                The *a* lattice parameter of the orthorhombic cell.
            b:
                The *b* lattice parameter of the orthorhombic cell.
            c:
                The *c* lattice parameter of the orthorhombic cell.

        Returns:
            Orthorhombic lattice of dimensions a x b x c.
        """
        return Lattice.from_parameters(a, b, c, 90, 90, 90)

    @staticmethod
    def monoclinic(a, b, c, alpha):
        """
        Convenience constructor for a monoclinic lattice.

        Args:
            a:
                The *a* lattice parameter of the monoclinc cell.
            b:
                The *b* lattice parameter of the monoclinc cell.
            c:
                The *c* lattice parameter of the monoclinc cell.
            alpha:
                The *alpha* angle between lattice vectors b and c.

        Returns:
            Monoclinic lattice of dimensions a x b x c with angle alpha between
            lattice vectors b and c.
        """
        return Lattice.from_parameters(a, b, c, alpha, 90, 90)

    @staticmethod
    def hexagonal(a, c):
        """
        Convenience constructor for a hexagonal lattice.

        Args:
            a:
                The *a* lattice parameter of the hexagonal cell.
            c:
                The *c* lattice parameter of the hexagonal cell.

        Returns:
            Hexagonal lattice of dimensions a x a x c.
        """
        return Lattice.from_parameters(a, a, c, 90.0, 90.0, 120.0)

    @staticmethod
    def rhombohedral(a, alpha):
        """
        Convenience constructor for a rhombohedral lattice.

        Args:
            a:
                The *a* lattice parameter of the rhombohedral cell.
            alpha:
                Angle for the rhombohedral lattice.

        Returns:
            Rhombohedral lattice of dimensions a x a x a.
        """
        return Lattice.from_parameters(a, a, a, alpha, alpha, alpha)

    @staticmethod
    def from_lengths_and_angles(abc, ang):
        """
        Create a Lattice using unit cell lengths and angles (in degrees).

        Args:
            abc:
                lattice parameters, e.g. (4, 4, 5).
            ang:
                lattice angles in degrees, e.g., (90,90,120).

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
            a:
                The *a* lattice parameter of the monoclinc cell.
            b:
                The *b* lattice parameter of the monoclinc cell.
            c:
                The *c* lattice parameter of the monoclinc cell.
            alpha:
                The *alpha* angle.
            beta:
                The *beta* angle.
            gamma:
                The *gamma* angle.

        Returns:
            A Lattice with the specified lattice parameters.
        """
        to_r = lambda degrees: np.radians(degrees)

        alpha_r = to_r(alpha)
        beta_r = to_r(beta)
        gamma_r = to_r(gamma)

        gamma_star = np.arccos((np.cos(alpha_r) * np.cos(beta_r)
                               - np.cos(gamma_r)) / (np.sin(alpha_r)
                               * np.sin(beta_r)))
        vector_a = [a * np.sin(to_r(beta)), 0.0, a * np.cos(to_r(beta))]
        vector_b = [-b * np.sin(to_r(alpha)) * np.cos(gamma_star), b
                    * np.sin(to_r(alpha)) * np.sin(gamma_star), b
                    * np.cos(to_r(alpha))]
        vector_c = np.array([0.0, 0.0, float(c)])
        return Lattice([vector_a, vector_b, vector_c])

    @staticmethod
    def from_dict(d):
        """
        Create a Lattice from a dictionary containing the a, b, c, alpha, beta,
        and gamma parameters.
        """
        a = d["a"]
        b = d["b"]
        c = d["c"]
        alpha = d["alpha"]
        beta = d["beta"]
        gamma = d["gamma"]
        return Lattice.from_parameters(a, b, c, alpha, beta, gamma)

    @property
    def angles(self):
        """
        Returns the angles (alpha, beta, gamma) of the lattice.
        """
        return self._angles

    @property
    def a(self):
        """
        *a* lattice parameter.
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
        return self._lengths

    @property
    def alpha(self):
        """
        Angle alpha of lattice.
        """
        return self._angles[0]

    @property
    def beta(self):
        """
        Angle beta of lattice.
        """
        return self._angles[1]

    @property
    def gamma(self):
        """
        Angle gamma of lattice.
        """
        return self._angles[2]

    @property
    def volume(self):
        """
        Volume of the unit cell.
        """
        return abs(npl.det(self._matrix))

    @property
    def lengths_and_angles(self):
        """
        Returns (lattice lengths, lattice angles).
        """
        return self._lengths, self._angles

    @property
    def reciprocal_lattice(self):
        """
        Return the reciprocal lattice.
        """
        v = 2 * np.pi / self.volume
        k1 = np.cross(self._matrix[1], self._matrix[2]) * v
        k2 = np.cross(self._matrix[2], self._matrix[0]) * v
        k3 = np.cross(self._matrix[0], self._matrix[1]) * v
        return Lattice([k1, k2, k3])

    def __repr__(self):
        f = lambda x: "%0.6f" % x
        outs = []
        outs.append("Lattice")
        outs.append("    abc : " + " ".join(map(f, self.abc)))
        outs.append(" angles : " + " ".join(map(f, self.angles)))
        outs.append(" volume : %0.4f" % self.volume)
        outs.append("      A : " + " ".join(map(f, self._matrix[0])))
        outs.append("      B : " + " ".join(map(f, self._matrix[1])))
        outs.append("      C : " + " ".join(map(f, self._matrix[2])))
        return "\n".join(outs)

    def __eq__(self, other):
        """
        A lattice is considered to be equal to another if the internal matrix
        representation satisfies np.allclose(matrix1, matrix2) to be True.
        """
        if other is None:
            return False
        return np.allclose(self._matrix, other._matrix)

    def __hash__(self):
        return 7

    def __str__(self):
        return "\n".join([" ".join(["%.6f" % i for i in row])
                          for row in self._matrix])

    @property
    def to_dict(self):
        """""
        Json-serialization dict representation of the Lattice.
        """
        d = {"matrix": self._matrix.tolist(),
             "a": float(self.a),
             "b": float(self.b),
             "c": float(self.c),
             "alpha": float(self.alpha),
             "beta": float(self.beta),
             "gamma": float(self.gamma),
             "volume": float(self.volume)
             }
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    def get_primitive_lattice(self, lattice_type):
        """
        Returns the primitive lattice of the lattice given the type.

        Args:
            lattice_type:
                An alphabet indicating whether the lattice type, P - primitive,
                R - rhombohedral, A - A-centered, B - B-centered,
                C - C-centered, I - body-centered, F - F-centered.
        """
        if lattice_type == "P":
            return Lattice(self._matrix)
        conv_to_prim = {
            "R": np.array([[2, 1, 1], [-1, 1, 1], [-1, -2, 1]]) / 3,
            "A": np.array([[2, 0, 0], [0, 1, 1], [0, -1, 1]]) / 2,
            "B": np.array([[1, 0, 1], [0, 2, 0], [-1, 0, 1]]) / 2,
            "C": np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 2]]) / 2,
            "I": np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]) / 2,
            "F": np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]]) / 2
        }
        return Lattice(np.dot(conv_to_prim[lattice_type], self._matrix))

    def get_most_compact_basis_on_lattice(self):
        """
        This method returns the alternative basis corresponding to the shortest
        3 linearly independent translational operations permitted.
        This tends to create larger angles for every elongated cells and is
        beneficial for viewing crystal structure (especially when they are
        Niggli cells).
        """
        matrix = self.matrix
        a = matrix[0]
        b = matrix[1]
        c = matrix[2]
        while True:
            anychange = False
            """
            take care of c
            """
            if np.dot(a, b) > 0:
                diffvector = a - b
            else:
                diffvector = a + b
            diffnorm = np.linalg.norm(diffvector)
            if diffnorm < np.linalg.norm(a) or \
                    diffnorm < np.linalg.norm(b):
                if np.linalg.norm(a) < np.linalg.norm(b):
                    b = diffvector
                else:
                    a = diffvector
            anychange = True
            """
            take care of b
            """
            if np.dot(a, c) > 0:
                diffvector = a - c
            else:
                diffvector = a + c
            diffnorm = np.linalg.norm(diffvector)
            if diffnorm < np.linalg.norm(a) or \
                    diffnorm < np.linalg.norm(c):
                if np.linalg.norm(a) < np.linalg.norm(c):
                    c = diffvector
                else:
                    a = diffvector
            anychange = True
            """
            take care of a
            """
            if np.dot(c, b) > 0:
                diffvector = c - b
            else:
                diffvector = c + b
            diffnorm = np.linalg.norm(diffvector)
            if diffnorm < np.linalg.norm(c) or \
                    diffnorm < np.linalg.norm(b):
                if np.linalg.norm(c) < np.linalg.norm(b):
                    b = diffvector
                else:
                    c = diffvector
            anychange = True
            if anychange:
                break
        return Lattice([a, b, c])

    def get_lll_reduced_lattice(self, delta=0.75):
        """
        Performs a Lenstra-Lenstra-Lovasz lattice basis reduction to obtain a
        c-reduced basis. This method returns a basis which is as "good" as
        possible, with "good" defined by orthongonality of the lattice vectors.

        Args:
            delta:
                Reduction parameter. Default of 0.75 is usually fine.

        Returns:
            Reduced lattice.
        """
        # Transpose the lattice matrix first so that basis vectors are columns.
        # Makes life easier.
        a = np.transpose(self._matrix.copy())

        b = np.zeros((3, 3))  # Vectors after the Gram-Schmidt process
        u = np.zeros((3, 3))  # Gram-Schmidt coeffieicnts
        m = np.zeros(3)  # These are the norm squared of each vec.

        # We are going to use dot a lot. Let's make an alias.
        dot = np.dot

        b[:, 0] = a[:, 0]
        m[0] = dot(b[:, 0], b[:, 0])
        for i in xrange(1, 3):
            u[i, 0:i] = dot(a[:, i].T, b[:, 0:i]) / m[0:i]
            b[:, i] = a[:, i] - dot(b[:, 0:i], u[i, 0:i].T)
            m[i] = dot(b[:, i], b[:, i])

        k = 2

        while k <= 3:
            # Size reduction.
            for i in xrange(k - 1, 0, -1):
                q = round(u[k - 1, i - 1])
                if q != 0:
                    # Reduce the k-th basis vector.
                    a[:, k - 1] = a[:, k - 1] - q * a[:, i - 1]
                    uu = list(u[i - 1, 0:(i - 1)])
                    uu.append(1)
                    # Update the GS coefficients.
                    u[k - 1, 0:i] = u[k - 1, 0:i] - q * np.array(uu)

            # Check the Lovasz condition.
            if dot(b[:, k - 1], b[:, k - 1]) >= \
                    (delta - abs(u[k - 1, k - 2]) ** 2) * \
                    dot(b[:, (k - 2)], b[:, (k - 2)]):
                # Increment k if the Lovasz condition holds.
                k += 1
            else:
                #If the Lovasz condition fails,
                #swap the k-th and (k-1)-th basis vector
                v = a[:, k - 1].copy()
                a[:, k - 1] = a[:, k - 2].copy()
                a[:, k - 2] = v
                #Update the Gram-Schmidt coefficients
                for s in xrange(k - 1, k + 1):
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

        return Lattice(a.T)
