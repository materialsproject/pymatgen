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

import numpy as np
import math

#array size threshold for looping instead of broadcasting
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


def coord_list_mapping(subset, superset):
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
    inds = np.where(np.all(np.isclose(c1[:, None, :], c2[None, :, :]),
                           axis=2))[1]
    result = c2[inds]
    if not np.allclose(c1, result):
        if not is_coord_subset(subset, superset):
            raise ValueError("subset is not a subset of superset")
    if not result.shape == c1.shape:
        raise ValueError("Something wrong with the inputs, likely duplicates "
                         "in superset")
    return inds


def coord_list_mapping_pbc(subset, superset, atol=1e-8):
    """
    Gives the index mapping from a subset to a superset.
    Subset and superset cannot contain duplicate rows

    Args:
        subset, superset: List of frac_coords

    Returns:
        list of indices such that superset[indices] = subset
    """
    c1 = np.array(subset)
    c2 = np.array(superset)

    diff = c1[:, None, :] - c2[None, :, :]
    diff -= np.round(diff)
    inds = np.where(np.all(np.abs(diff) < atol, axis = 2))[1]

    #verify result (its easier to check validity of the result than
    #the validity of inputs)
    test = c2[inds] - c1
    test -= np.round(test)
    if not np.allclose(test, 0):
        if not is_coord_subset_pbc(subset, superset):
            raise ValueError("subset is not a subset of superset")
    if not test.shape == c1.shape:
        raise ValueError("Something wrong with the inputs, likely duplicates "
                         "in superset")
    return inds


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

#create images, 2d array of all length 3 combinations of [-1,0,1]
r = np.arange(-1, 2)
arange = r[:, None] * np.array([1, 0, 0])[None, :]
brange = r[:, None] * np.array([0, 1, 0])[None, :]
crange = r[:, None] * np.array([0, 0, 1])[None, :]
images = arange[:, None, None] + brange[None, :, None] + \
    crange[None, None, :]
images = images.reshape((27, 3))

def pbc_shortest_vectors(lattice, fcoords1, fcoords2):
    """
    Returns the shortest vectors between two lists of coordinates taking into
    account periodic boundary conditions and the lattice.

    Args:
        lattice: lattice to use
        fcoords1: First set of fractional coordinates. e.g., [0.5, 0.6, 0.7]
            or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
            coord or any array of coords.
        fcoords2: Second set of fractional coordinates.

    Returns:
        array of displacement vectors from fcoords1 to fcoords2
        first index is fcoords1 index, second is fcoords2 index
    """
    #ensure correct shape
    fcoords1, fcoords2 = np.atleast_2d(fcoords1, fcoords2)

    #ensure that all points are in the unit cell
    fcoords1 = np.mod(fcoords1, 1)
    fcoords2 = np.mod(fcoords2, 1)

    #create images of f2
    shifted_f2 = fcoords2[:, None, :] + images[None, :, :]

    cart_f1 = lattice.get_cartesian_coords(fcoords1)
    cart_f2 = lattice.get_cartesian_coords(shifted_f2)

    if cart_f1.size * cart_f2.size < LOOP_THRESHOLD:
        #all vectors from f1 to f2
        vectors = cart_f2[None, :, :, :] - cart_f1[:, None, None, :]
        d_2 = np.sum(vectors ** 2, axis=-1)
        a, b = np.indices(vectors.shape[:2])
        return vectors[a, b, np.argmin(d_2, axis=-1)]
    else:
        shortest = np.zeros((len(cart_f1), len(cart_f2), 3))
        for i, c1 in enumerate(cart_f1):
            vectors = cart_f2[:, :, :] - c1[None, None, :]
            d_2 = np.sum(vectors ** 2, axis=-1)
            shortest[i] = vectors[np.arange(len(vectors)),
                                  np.argmin(d_2, axis=-1)]
        return shortest


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
    c1 = np.array(subset)
    c2 = np.array(superset)
    dist = c1[:, None, :] - c2[None, :, :]
    dist -= np.round(dist)
    if mask is not None:
        dist[np.array(mask)] = np.inf
    is_close = np.all(np.abs(dist) < atol, axis=-1)
    any_close = np.any(is_close, axis=-1)
    return np.all(any_close)


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
