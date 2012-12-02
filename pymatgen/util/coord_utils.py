#!/usr/bin/env python

"""
Utilities for manipulating coordinates or list of coordinates,
in periodic boundary conditions or otherwise.
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
import itertools


def in_coord_list(coord_list, coord, **kwargs):
    """
    Tests if a particular coord is within a coord_list using numpy.allclose.

    Args:
        coord_list:
            List of coords to test
        coord:
            Specific coordinates
        kwargs:
            keyword arguments supported by numpy.allclose. Please refer to
            documentation for numpy.allclose.
    """
    for test in coord_list:
        if np.allclose(test, coord, **kwargs):
            return True
    return False


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


def pbc_diff(fcoord1, fcoord2):
    """
    Returns the 'fractional distance' between two coordinates taking into
    account periodic boundary conditions.

    Args:
        fcoord1:
            First fractional coordinate.
        fcoord2:
            Second fractional coordinate.

    Returns:
        Fractional distance. Each coordinate must have the property that
        abs(a) <= 0.5. Examples:
        pbc_diff([0.1, 0.1, 0.1], [0.3, 0.5, 0.9]) = [-0.2, -0.4, 0.2]
        pbc_diff([0.9, 0.1, 1.01], [0.3, 0.5, 0.9]) = [-0.4, -0.4, 0.11]
    """
    fdist = np.mod(fcoord1, 1) - np.mod(fcoord2, 1)
    return [a if abs(a) <= 0.5 else a - a/abs(a) for a in fdist]


def in_coord_list_pbc(coord_list, coord, **kwargs):
    """
    Tests if a particular coord is within a coord_list using numpy.allclose.

    Args:
        coord_list:
            List of coords to test
        coord:
            Specific coordinates
        kwargs:
            keyword arguments supported by numpy.allclose. Please refer to
            documentation for numpy.allclose.
    """
    for test in coord_list:
        fdiff = pbc_diff(test, coord)
        if np.allclose(fdiff, [0,0,0], **kwargs):
            return True
    return False


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
        [(site, dist) ...] since most of the time, subsequent processing
        requires the distance.
    """
    recp_len = lattice.reciprocal_lattice.abc
    sr = r + 0.15
    nmax = [sr * l / (2 * math.pi) for l in recp_len]
    pcoords = lattice.get_fractional_coords(center)
    axis_ranges = []
    floor = math.floor
    for i in range(3):
        rangemax = int(floor(pcoords[i] + nmax[i]))
        rangemin = int(floor(pcoords[i] - nmax[i]))
        axis_ranges.append(range(rangemin, rangemax + 1))
    neighbors = []
    n = len(frac_points)
    fcoords = np.array(frac_points)
    frac_2_cart = lattice.get_cartesian_coords
    pts = np.tile(center, (n, 1))
    indices = np.array(range(n))
    for image in itertools.product(*axis_ranges):
        shift = np.tile(image, (n, 1))
        shifted_coords = fcoords + shift
        coords = frac_2_cart(shifted_coords)
        dists = np.sqrt(np.sum((coords - pts) ** 2, axis=1))
        within_r = dists <= r
        d = [shifted_coords[within_r], dists[within_r], indices[within_r]]
        neighbors.extend(np.transpose(d))
    return neighbors
