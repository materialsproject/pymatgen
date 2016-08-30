# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Utilities for manipulating coordinates or list of coordinates, under periodic
boundary conditions or otherwise. Many of these are heavily vectorized in
numpy for performance.
"""

from six.moves import zip

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Nov 27, 2011"

import itertools
import numpy as np
import math
from . import coord_utils_cython as cuc


# array size threshold for looping instead of broadcasting
LOOP_THRESHOLD = 1e6


def find_in_coord_list(coord_list, coord, atol=1e-8):
    """
    Find the indices of matches of a particular coord in a coord_list.

    Args:
        coord_list: List of coords to test
        coord: Specific coordinates
        atol: Absolute tolerance. Defaults to 1e-8. Accepts both scalar and
            array.

    Returns:
        Indices of matches, e.g., [0, 1, 2, 3]. Empty list if not found.
    """
    if len(coord_list) == 0:
        return []
    diff = np.array(coord_list) - np.array(coord)[None, :]
    return np.where(np.all(np.abs(diff) < atol, axis=1))[0]


def in_coord_list(coord_list, coord, atol=1e-8):
    """
    Tests if a particular coord is within a coord_list.

    Args:
        coord_list: List of coords to test
        coord: Specific coordinates
        atol: Absolute tolerance. Defaults to 1e-8. Accepts both scalar and
            array.

    Returns:
        True if coord is in the coord list.
    """
    return len(find_in_coord_list(coord_list, coord, atol=atol)) > 0


def is_coord_subset(subset, superset, atol=1e-8):
    """
    Tests if all coords in subset are contained in superset.
    Doesn't use periodic boundary conditions

    Args:
        subset, superset: List of coords

    Returns:
        True if all of subset is in superset.
    """
    c1 = np.array(subset)
    c2 = np.array(superset)
    is_close = np.all(np.abs(c1[:, None, :] - c2[None, :, :]) < atol, axis=-1)
    any_close = np.any(is_close, axis=-1)
    return np.all(any_close)


def coord_list_mapping(subset, superset, atol=1e-8):
    """
    Gives the index mapping from a subset to a superset.
    Subset and superset cannot contain duplicate rows

    Args:
        subset, superset: List of coords

    Returns:
        list of indices such that superset[indices] = subset
    """
    c1 = np.array(subset)
    c2 = np.array(superset)
    inds = np.where(np.all(np.isclose(c1[:, None, :], c2[None, :, :], atol=atol),
                           axis=2))[1]
    result = c2[inds]
    if not np.allclose(c1, result, atol=atol):
        if not is_coord_subset(subset, superset):
            raise ValueError("subset is not a subset of superset")
    if not result.shape == c1.shape:
        raise ValueError("Something wrong with the inputs, likely duplicates "
                         "in superset")
    return inds


def coord_list_mapping_pbc(subset, superset, atol=1e-8):
    """
    Gives the index mapping from a subset to a superset.
    Superset cannot contain duplicate matching rows

    Args:
        subset, superset: List of frac_coords

    Returns:
        list of indices such that superset[indices] = subset
    """
    atol = np.array([1., 1. ,1.]) * atol
    return cuc.coord_list_mapping_pbc(subset, superset, atol)


def get_linear_interpolated_value(x_values, y_values, x):
    """
    Returns an interpolated value by linear interpolation between two values.
    This method is written to avoid dependency on scipy, which causes issues on
    threading servers.

    Args:
        x_values: Sequence of x values.
        y_values: Corresponding sequence of y values
        x: Get value at particular x

    Returns:
        Value at x.
    """
    a = np.array(sorted(zip(x_values, y_values), key=lambda d: d[0]))

    ind = np.where(a[:, 0] >= x)[0]

    if len(ind) == 0 or ind[0] == 0:
        raise ValueError("x is out of range of provided x_values")

    i = ind[0]
    x1, x2 = a[i - 1][0], a[i][0]
    y1, y2 = a[i - 1][1], a[i][1]

    return y1 + (y2 - y1) / (x2 - x1) * (x - x1)


def all_distances(coords1, coords2):
    """
    Returns the distances between two lists of coordinates

    Args:
        coords1: First set of cartesian coordinates.
        coords2: Second set of cartesian coordinates.

    Returns:
        2d array of cartesian distances. E.g the distance between
        coords1[i] and coords2[j] is distances[i,j]
    """
    c1 = np.array(coords1)
    c2 = np.array(coords2)
    z = (c1[:, None, :] - c2[None, :, :]) ** 2
    return np.sum(z, axis=-1) ** 0.5


def pbc_diff(fcoords1, fcoords2):
    """
    Returns the 'fractional distance' between two coordinates taking into
    account periodic boundary conditions.

    Args:
        fcoords1: First set of fractional coordinates. e.g., [0.5, 0.6,
            0.7] or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
            coord or any array of coords.
        fcoords2: Second set of fractional coordinates.

    Returns:
        Fractional distance. Each coordinate must have the property that
        abs(a) <= 0.5. Examples:
        pbc_diff([0.1, 0.1, 0.1], [0.3, 0.5, 0.9]) = [-0.2, -0.4, 0.2]
        pbc_diff([0.9, 0.1, 1.01], [0.3, 0.5, 0.9]) = [-0.4, -0.4, 0.11]
    """
    fdist = np.subtract(fcoords1, fcoords2)
    return fdist - np.round(fdist)

def pbc_shortest_vectors(lattice, fcoords1, fcoords2, mask=None, return_d2=False):
    """
    Returns the shortest vectors between two lists of coordinates taking into
    account periodic boundary conditions and the lattice.

    Args:
        lattice: lattice to use
        fcoords1: First set of fractional coordinates. e.g., [0.5, 0.6, 0.7]
            or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
            coord or any array of coords.
        fcoords2: Second set of fractional coordinates.
        mask (boolean array): Mask of matches that are not allowed.
            i.e. if mask[1,2] == True, then subset[1] cannot be matched
            to superset[2]
        return_d2 (boolean): whether to also return the squared distances

    Returns:
        array of displacement vectors from fcoords1 to fcoords2
        first index is fcoords1 index, second is fcoords2 index
    """
    return cuc.pbc_shortest_vectors(lattice, fcoords1, fcoords2, mask, return_d2)


def find_in_coord_list_pbc(fcoord_list, fcoord, atol=1e-8):
    """
    Get the indices of all points in a fractional coord list that are
    equal to a fractional coord (with a tolerance), taking into account
    periodic boundary conditions.

    Args:
        fcoord_list: List of fractional coords
        fcoord: A specific fractional coord to test.
        atol: Absolute tolerance. Defaults to 1e-8.

    Returns:
        Indices of matches, e.g., [0, 1, 2, 3]. Empty list if not found.
    """
    if len(fcoord_list) == 0:
        return []
    fcoords = np.tile(fcoord, (len(fcoord_list), 1))
    fdist = fcoord_list - fcoords
    fdist -= np.round(fdist)
    return np.where(np.all(np.abs(fdist) < atol, axis=1))[0]


def in_coord_list_pbc(fcoord_list, fcoord, atol=1e-8):
    """
    Tests if a particular fractional coord is within a fractional coord_list.

    Args:
        fcoord_list: List of fractional coords to test
        fcoord: A specific fractional coord to test.
        atol: Absolute tolerance. Defaults to 1e-8.

    Returns:
        True if coord is in the coord list.
    """
    return len(find_in_coord_list_pbc(fcoord_list, fcoord, atol=atol)) > 0


def is_coord_subset_pbc(subset, superset, atol=1e-8, mask=None):
    """
    Tests if all fractional coords in subset are contained in superset.

    Args:
        subset, superset: List of fractional coords
        atol (float or size 3 array): Tolerance for matching
        mask (boolean array): Mask of matches that are not allowed.
            i.e. if mask[1,2] == True, then subset[1] cannot be matched
            to superset[2]

    Returns:
        True if all of subset is in superset.
    """
    c1 = np.array(subset, dtype=np.float64)
    c2 = np.array(superset, dtype=np.float64)
    if mask is not None:
        m = np.array(mask, dtype=np.int)
    else:
        m = np.zeros((len(subset), len(superset)), dtype=np.int)
    atol = np.zeros(3, dtype=np.float64) + atol
    return cuc.is_coord_subset_pbc(c1, c2, atol, m)


def lattice_points_in_supercell(supercell_matrix):
    """
    Returns the list of points on the original lattice contained in the
    supercell in fractional coordinates (with the supercell basis).
    e.g. [[2,0,0],[0,1,0],[0,0,1]] returns [[0,0,0],[0.5,0,0]]

    Args:
        supercell_matrix: 3x3 matrix describing the supercell

    Returns:
        numpy array of the fractional coordinates
    """
    diagonals = np.array(
        [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1],
         [1, 1, 0], [1, 1, 1]])
    d_points = np.dot(diagonals, supercell_matrix)

    mins = np.min(d_points, axis=0)
    maxes = np.max(d_points, axis=0) + 1

    ar = np.arange(mins[0], maxes[0])[:, None] * \
         np.array([1, 0, 0])[None, :]
    br = np.arange(mins[1], maxes[1])[:, None] * \
         np.array([0, 1, 0])[None, :]
    cr = np.arange(mins[2], maxes[2])[:, None] * \
         np.array([0, 0, 1])[None, :]

    all_points = ar[:, None, None] + br[None, :, None] + cr[None, None, :]
    all_points = all_points.reshape((-1, 3))

    frac_points = np.dot(all_points, np.linalg.inv(supercell_matrix))

    tvects = frac_points[np.all(frac_points < 1 - 1e-10, axis=1)
                         & np.all(frac_points >= -1e-10, axis=1)]
    assert len(tvects) == round(abs(np.linalg.det(supercell_matrix)))
    return tvects


def barycentric_coords(coords, simplex):
    """
    Converts a list of coordinates to barycentric coordinates, given a
    simplex with d+1 points. Only works for d >= 2.

    Args:
        coords: list of n coords to transform, shape should be (n,d)
        simplex: list of coordinates that form the simplex, shape should be
            (d+1, d)

    Returns:
        a LIST of barycentric coordinates (even if the original input was 1d)
    """
    coords = np.atleast_2d(coords)

    t = np.transpose(simplex[:-1, :]) - np.transpose(simplex[-1, :])[:, None]
    all_but_one = np.transpose(
        np.linalg.solve(t, np.transpose(coords - simplex[-1])))
    last_coord = 1 - np.sum(all_but_one, axis=-1)[:, None]
    return np.append(all_but_one, last_coord, axis=-1)


def get_angle(v1, v2, units="degrees"):
    """
    Calculates the angle between two vectors.

    Args:
        v1: Vector 1
        v2: Vector 2
        units: "degrees" or "radians". Defaults to "degrees".

    Returns:
        Angle between them in degrees.
    """
    d = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
    d = min(d, 1)
    d = max(d, -1)
    angle = math.acos(d)
    if units == "degrees":
        return math.degrees(angle)
    elif units == "radians":
        return angle
    else:
        raise ValueError("Invalid units {}".format(units))


class Simplex(object):
    """
    A generalized simplex object. See http://en.wikipedia.org/wiki/Simplex.

    .. attribute: space_dim

        Dimension of the space. Usually, this is 1 more than the simplex_dim.

    .. attribute: simplex_dim

        Dimension of the simplex coordinate space.
    """

    def __init__(self, coords):
        """
        Initializes a Simplex from vertex coordinates.

        Args:
            coords ([[float]]): Coords of the vertices of the simplex. E.g.,
                [[1, 2, 3], [2, 4, 5], [6, 7, 8], [8, 9, 10].
        """
        self._coords = np.array(coords)
        self.space_dim, self.simplex_dim = self._coords.shape
        self.origin = self._coords[-1]
        if self.space_dim == self.simplex_dim + 1:
            # precompute attributes for calculating bary_coords
            self.T = self._coords[:-1] - self.origin
            self.T_inv = np.linalg.inv(self.T)

    @property
    def volume(self):
        """
        Volume of the simplex.
        """
        return abs(np.linalg.det(self.T)) / math.factorial(self.simplex_dim)

    def bary_coords(self, point):
        try:
            c = np.dot((point - self.origin), self.T_inv)
            return np.concatenate([c, [1 - np.sum(c)]])
        except AttributeError:
            raise ValueError('Simplex is not full-dimensional')

    def in_simplex(self, point, tolerance=1e-8):
        """
        Checks if a point is in the simplex using the standard barycentric
        coordinate system algorithm.

        Taking an arbitrary vertex as an origin, we compute the basis for the
        simplex from this origin by subtracting all other vertices from the
        origin. We then project the point into this coordinate system and
        determine the linear decomposition coefficients in this coordinate
        system.  If the coeffs satisfy that all coeffs >= 0, the composition
        is in the facet.

        Args:
            point ([float]): Point to test
            tolerance (float): Tolerance to test if point is in simplex.
        """
        return (self.bary_coords(point) >= -tolerance).all()

    def __eq__(self, other):
        for p in itertools.permutations(self._coords):
            if np.allclose(p, other.coords):
                return True
        return False

    def __hash__(self):
        return len(self._coords)

    def __repr__(self):
        output = ["{}-simplex in {}D space".format(self.simplex_dim,
                                                   self.space_dim),
                  "Vertices:"]
        for coord in self._coords:
            output.append("\t({})".format(", ".join(map(str, coord))))
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()

    @property
    def coords(self):
        """
        Returns a copy of the vertex coordinates in the simplex.
        """
        return self._coords.copy()
