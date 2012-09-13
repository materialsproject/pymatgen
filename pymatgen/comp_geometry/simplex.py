#!/usr/bin/env python

'''
This module defines a class representing an arbitrary Simplex in arbitrary
dimensional space.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "May 15, 2012"

import itertools

import numpy as np


class Simplex(object):
    """
    A generalized simplex object. See http://en.wikipedia.org/wiki/Simplex.
    """

    def __init__(self, coords):
        """
        Args:
            coords:
                Coords of the vertices of the simplex. E.g.,
                [[1,2,3],[2,4,5],[6,7,8]].
        """
        self.simplex_dim = len(coords) - 1
        self.space_dim = len(coords[0])
        for c in coords:
            if len(c) != self.space_dim:
                raise ValueError("All coords must have the same space "
                                 "dimension.")
        self._coords = np.array(coords)

    def in_simplex(self, point, tolerance=1e-8):
        """
        Checks if a point is in the simplex using the standard barycentric
        coordinate system algorithm.

        Taking an arbitrary vertex as an origin, we compute the basis for the
        simplex from this origin by subtracting all other vertices from the
        origin. We then project the point into this coordinate system and
        determine the linear decomposition coefficients in this coordinate
        system.  If the coeffs satisfy that all coeffs >= 0 and
        sum(coeffs) <= 1, the composition is in the facet.

        For example, take a tetrahedron. For a tetrahedron, let's label
        the vertices as O, A, B anc C.  Let's call our point coordinate as X.
        We form the composition matrix M with vectors OA, OB and OB, transponse
        it, and solve for M'.a = OX where a are the coefficients.

        If (a >= 0).all() and sum(a) <= 1, X is in the tetrahedron.
        Note that in reality, the test needs to provide a tolerance (set to
        1e-8 by default) for numerical errors.

        Args:
            point:
                Point to test
            tolerance:
                Tolerance to test if point is in simplex.
        """
        origin = self._coords[0]
        bary_coords = np.array([self._coords[i] - origin
                                for i in xrange(1, self.simplex_dim + 1)])
        shifted_point = np.array(point) - origin
        coeffs = np.linalg.solve(bary_coords.transpose(), shifted_point)
        return (coeffs >= -tolerance).all() and sum(coeffs) <= (1 + tolerance)

    def __eq__(self, other):
        for p in itertools.permutations(self._coords):
            if np.allclose(p, other.coords):
                return True
        return False

    def __hash__(self):
        return 7

    def __repr__(self):
        output = ["{}-simplex in {}D space".format(self.simplex_dim,
                                                   self.space_dim)]
        output.append("Vertices:")
        for coord in self._coords:
            output.append("\t({})".format(", ".join([str(i) for i in coord])))
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()

    @property
    def coords(self):
        return self._coords.copy()
