#!/usr/bin/env python

"""
Utilities for manipulating coordinates or list of coordinates, under periodic
boundary conditions or otherwise. Many of these are heavily vectorized in
numpy for performance.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 27, 2011"

import numpy as np
import math


def find_in_coord_list(coord_list, coord, atol=1e-8):
    """
    Find the indices of matches of a particular coord in a coord_list.

    Args:
        coord_list:
            List of coords to test
        coord:
            Specific coordinates
        atol:
            Absolute tolerance. Defaults to 1e-8.
            Accepts both scalar and array

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
        coord_list:
            List of coords to test
        coord:
            Specific coordinates
        atol:
            Absolute tolerance. Defaults to 1e-8.
            Accepts both scalar and array
    Returns:
        True if coord is in the coord list.
    """
    return len(find_in_coord_list(coord_list, coord, atol=atol)) > 0


def get_linear_interpolated_value(x_values, y_values, x):
    """
    Returns an interpolated value by linear interpolation between two values.
    This method is written to avoid dependency on scipy, which causes issues on
    threading servers.

    Args:
        x_values:
            Sequence of x values.
        y_values:
            Corresponding sequence of y values
        x:
            Get value at particular x
    """
    val_dict = dict(zip(x_values, y_values))
    if x < min(x_values) or x > max(x_values):
        raise ValueError("x is out of range of provided x_values")

    for val in sorted(x_values):
        if val < x:
            x1 = val
        else:
            x2 = val
            break

    y1 = val_dict[x1]
    y2 = val_dict[x2]

    return y1 + (y2 - y1) / (x2 - x1) * (x - x1)


def pbc_diff(fcoords1, fcoords2):
    """
    Returns the 'fractional distance' between two coordinates taking into
    account periodic boundary conditions.

    Args:
        fcoords1:
            First set of fractional coordinates. e.g., [0.5, 0.6,
            0.7] or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
            coord or any array of coords.
        fcoords2:
            Second set of fractional coordinates.

    Returns:
        Fractional distance. Each coordinate must have the property that
        abs(a) <= 0.5. Examples:
        pbc_diff([0.1, 0.1, 0.1], [0.3, 0.5, 0.9]) = [-0.2, -0.4, 0.2]
        pbc_diff([0.9, 0.1, 1.01], [0.3, 0.5, 0.9]) = [-0.4, -0.4, 0.11]
    """
    fdist = np.subtract(fcoords1, fcoords2)
    return fdist - np.round(fdist)


def find_in_coord_list_pbc(fcoord_list, fcoord, atol=1e-8):
    """
    Get the indices of all points in a fracitonal coord list that are
    equal to a fractional coord (with a tolerance), taking into account
    periodic boundary conditions.

    Args:
        fcoord_list:
            List of fractional coords
        fcoord:
            A specific fractional coord to test.
        atol:
            Absolute tolerance. Defaults to 1e-8.

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
        fcoord_list:
            List of fractional coords to test
        fcoord:
            A specific fractional coord to test.
        atol:
            Absolute tolerance. Defaults to 1e-8.

    Returns:
        True if coord is in the coord list.
    """
    return len(find_in_coord_list_pbc(fcoord_list, fcoord, atol=atol)) > 0


def get_points_in_sphere_pbc(lattice, frac_points, center, r):
    """
    Find all points within a sphere from the point taking into account
    periodic boundary conditions. This includes sites in other periodic images.

    Algorithm:

    1. place sphere of radius r in crystal and determine minimum supercell
       (parallelpiped) which would contain a sphere of radius r. for this
       we need the projection of a_1 on a unit vector perpendicular
       to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
       determine how many a_1"s it will take to contain the sphere.

       Nxmax = r * length_of_b_1 / (2 Pi)

    2. keep points falling within r.

    Args:
        lattice:
            The lattice/basis for the periodic boundary conditions.
        frac_points:
            All points in the lattice in fractional coordinates.
        center:
            cartesian coordinates of center of sphere.
        r:
            radius of sphere.

    Returns:
        [(fcoord, dist) ...] since most of the time, subsequent processing
        requires the distance.
    """
    recp_len = np.array(lattice.reciprocal_lattice.abc)
    sr = r + 0.15
    nmax = sr * recp_len / (2 * math.pi)
    pcoords = lattice.get_fractional_coords(center)
    floor = math.floor

    n = len(frac_points)
    fcoords = np.array(frac_points)
    frac_2_cart = lattice.get_cartesian_coords
    pts = np.tile(center, (n, 1))
    indices = np.array(range(n))

    arange = np.arange(start=int(floor(pcoords[0] - nmax[0])),
                       stop=int(floor(pcoords[0] + nmax[0])) + 1)
    brange = np.arange(start=int(floor(pcoords[1] - nmax[1])),
                       stop=int(floor(pcoords[1] + nmax[1])) + 1)
    crange = np.arange(start=int(floor(pcoords[2] - nmax[2])),
                       stop=int(floor(pcoords[2] + nmax[2])) + 1)

    arange = arange[:, None] * np.array([1, 0, 0])[None, :]
    brange = brange[:, None] * np.array([0, 1, 0])[None, :]
    crange = crange[:, None] * np.array([0, 0, 1])[None, :]

    images = arange[:, None, None] + brange[None, :, None] + \
        crange[None, None, :]

    shifted_coords = fcoords[:, None, None, None, :] + images[None, :, :, :, :]
    coords = frac_2_cart(shifted_coords)
    dists = np.sqrt(np.sum((coords - pts[:, None, None, None, :]) ** 2, axis=4))
    within_r = np.where(dists <= r)

    d = [shifted_coords[within_r], dists[within_r], indices[within_r[0]]]

    return np.transpose(d)
    
