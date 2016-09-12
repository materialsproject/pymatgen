# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module contains some utility functions and classes that are used in the chemenv package.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"


import math

import numpy as np
from numpy.linalg import norm
from scipy.spatial import ConvexHull
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline

from pymatgen.analysis.chemenv.utils.chemenv_errors import SolidAngleError


def get_lower_and_upper_f(surface_calculation_options):
    mindist = surface_calculation_options['distance_bounds']['lower']
    maxdist = surface_calculation_options['distance_bounds']['upper']
    minang = surface_calculation_options['angle_bounds']['lower']
    maxang = surface_calculation_options['angle_bounds']['upper']
    if surface_calculation_options['type'] == 'standard_elliptic':
        lower_and_upper_functions = quarter_ellipsis_functions(xx=(mindist, maxang), yy=(maxdist, minang))
    elif surface_calculation_options['type'] == 'standard_diamond':
        deltadist = surface_calculation_options['distance_bounds']['delta']
        deltaang = surface_calculation_options['angle_bounds']['delta']
        lower_and_upper_functions = diamond_functions(xx=(mindist, maxang), yy=(maxdist, minang),
                                                      x_y0=deltadist, y_x0=deltaang)
    elif surface_calculation_options['type'] == 'standard_spline':
        lower_points = surface_calculation_options['lower_points']
        upper_points = surface_calculation_options['upper_points']
        degree = surface_calculation_options['degree']
        lower_and_upper_functions = spline_functions(lower_points=lower_points,
                                                     upper_points=upper_points,
                                                     degree=degree)
    else:
        raise ValueError('Surface calculation of type "{}" '
                         'is not implemented'.format(surface_calculation_options['type']))
    return lower_and_upper_functions


def function_comparison(f1, f2, x1, x2, numpoints_check=500):
    """
        Method that compares two functions

        Args:
            f1: First function to compare
            f2: Second function to compare
            x1: Lower bound of the interval to compare
            x2: Upper bound of the interval to compare
            numpoints_check: Number of points used to compare the functions

        Returns:
            Whether the function are equal ("="), f1 is always lower than f2 ("<"), f1 is always larger than f2 (">"),
             f1 is always lower than or equal to f2 ("<"), f1 is always larger than or equal to f2 (">") on the
             interval [x1, x2]. If the two functions cross, a RuntimeError is thrown (i.e. we expect to compare
             functions that do not cross...)
        """
    xx = np.linspace(x1, x2, num=numpoints_check)
    y1 = f1(xx)
    y2 = f2(xx)
    if np.all(y1 < y2):
        return '<'
    elif np.all(y1 > y2):
        return '>'
    elif np.all(y1 == y2):
        return '='
    elif np.all(y1 <= y2):
        return '<='
    elif np.all(y1 >= y2):
        return '>='
    else:
        raise RuntimeError('Error in comparing functions f1 and f2 ...')


def quarter_ellipsis_functions(xx, yy):
    """
    Method that creates two quarter-ellipse functions based on points xx and yy. The ellipsis is supposed to
    be aligned with the axes. The two ellipsis pass through the two points xx and yy.

    Args:
        xx:
            First point
        yy:
            Second point

    Returns:
        A dictionary with the lower and upper quarter ellipsis functions.
    """
    npxx = np.array(xx)
    npyy = np.array(yy)
    if np.any(npxx == npyy):
        raise RuntimeError('Invalid points for quarter_ellipsis_functions')
    if np.all(npxx < npyy) or np.all(npxx > npyy):
        if npxx[0] < npyy[0]:
            p1 = npxx
            p2 = npyy
        else:
            p1 = npyy
            p2 = npxx
        c_lower = np.array([p1[0], p2[1]])
        c_upper = np.array([p2[0], p1[1]])
        b2 = (p2[1] - p1[1]) ** 2
    else:
        if npxx[0] < npyy[0]:
            p1 = npxx
            p2 = npyy
        else:
            p1 = npyy
            p2 = npxx
        c_lower = np.array([p2[0], p1[1]])
        c_upper = np.array([p1[0], p2[1]])
        b2 = (p1[1] - p2[1]) ** 2
    b2overa2 = b2 / (p2[0] - p1[0]) ** 2

    def lower(x):
        return c_lower[1] - np.sqrt(b2 - b2overa2 * (x - c_lower[0]) ** 2)

    def upper(x):
        return c_upper[1] + np.sqrt(b2 - b2overa2 * (x - c_upper[0]) ** 2)

    return {'lower': lower, 'upper': upper}


def spline_functions(lower_points, upper_points, degree=3):
    """
    Method that creates two (upper and lower) spline functions based on points lower_points and upper_points.

    Args:
        lower_points:
            Points defining the lower function.
        upper_points:
            Points defining the upper function.
        degree:
            Degree for the spline function

    Returns:
        A dictionary with the lower and upper spline functions.
    """
    lower_xx = np.array([pp[0] for pp in lower_points])
    lower_yy = np.array([pp[1] for pp in lower_points])
    upper_xx = np.array([pp[0] for pp in upper_points])
    upper_yy = np.array([pp[1] for pp in upper_points])

    lower_spline = UnivariateSpline(lower_xx, lower_yy, k=degree, s=0)
    upper_spline = UnivariateSpline(upper_xx, upper_yy, k=degree, s=0)

    def lower(x):
        return lower_spline(x)

    def upper(x):
        return upper_spline(x)

    return {'lower': lower, 'upper': upper}

def diamond_functions(xx, yy, y_x0, x_y0):
    """
    Method that creates two upper and lower functions based on points xx and yy as well as intercepts defined by
     y_x0 and x_y0. The resulting functions form kind of a distorted diamond-like structure aligned from
     point xx to point yy.

    Schematically :

    xx is symbolized by x, yy is symbolized by y, y_x0 is equal to the distance from x to a, x_y0 is equal to the
     distance from x to b, the lines a-p and b-q are parallel to the line x-y such that points p and q are
     obtained automatically.
    In case of an increasing diamond the lower function is x-b-q and the upper function is a-p-y while in case of a
     decreasing diamond, the lower function is a-p-y and the upper function is x-b-q.

           Increasing diamond      |     Decreasing diamond
                     p--y                    x----b
                    /  /|                    |\    \
                   /  / |                    | \    q
                  /  /  |                    a  \   |
                 a  /   |                     \  \  |
                 | /    q                      \  \ |
                 |/    /                        \  \|
                 x----b                          p--y

    Args:
        xx:
            First point
        yy:
            Second point

    Returns:
        A dictionary with the lower and upper diamond functions.
    """
    npxx = np.array(xx)
    npyy = np.array(yy)
    if np.any(npxx == npyy):
        raise RuntimeError('Invalid points for diamond_functions')
    if np.all(npxx < npyy) or np.all(npxx > npyy):
        if npxx[0] < npyy[0]:
            p1 = npxx
            p2 = npyy
        else:
            p1 = npyy
            p2 = npxx
    else:
        if npxx[0] < npyy[0]:
            p1 = npxx
            p2 = npyy
        else:
            p1 = npyy
            p2 = npxx
    slope = (p2[1]-p1[1]) / (p2[0]- p1[0])
    if slope > 0.0:
        x_bpoint = p1[0] + x_y0
        myy = p1[1]
        bq_intercept = myy - slope*x_bpoint
        myx = p1[0]
        myy = p1[1] + y_x0
        ap_intercept = myy - slope*myx
        x_ppoint = (p2[1] - ap_intercept) / slope
        def lower(x):
            return np.where(x <= x_bpoint, p1[1] * np.ones_like(x), slope * x + bq_intercept)

        def upper(x):
            return np.where(x >= x_ppoint, p2[1] * np.ones_like(x), slope * x + ap_intercept)
    else:
        x_bpoint = p1[0] + x_y0
        myy = p1[1]
        bq_intercept = myy - slope * x_bpoint
        myx = p1[0]
        myy = p1[1] - y_x0
        ap_intercept = myy - slope * myx
        x_ppoint = (p2[1] - ap_intercept) / slope
        def lower(x):
            return np.where(x >= x_ppoint, p2[1] * np.ones_like(x), slope * x + ap_intercept)

        def upper(x):
            return np.where(x <= x_bpoint, p1[1] * np.ones_like(x), slope * x + bq_intercept)

    return {'lower': lower, 'upper': upper}


def rectangle_surface_intersection(rectangle, f_lower, f_upper,
                                   bounds_lower=None, bounds_upper=None,
                                   check=True, numpoints_check=500):
    """
    Method to calculate the surface of the intersection of a rectangle (aligned with axes) and another surface
    defined by two functions f_lower and f_upper.

    Args:
        rectangle:
            Rectangle defined as : ((x1, x2), (y1, y2)).
        f_lower:
            Function defining the lower bound of the surface.
        f_upper:
            Function defining the upper bound of the surface.
        bounds_lower:
            Interval in which the f_lower function is defined.
        bounds_upper:
            Interval in which the f_upper function is defined.
        check:
            Whether to check if f_lower is always lower than f_upper.
        numpoints_check:
            Number of points used to check whether f_lower is always lower than f_upper

    Returns:
        The surface of the intersection of the rectangle and the surface defined by f_lower and f_upper.
    """
    x1 = np.min(rectangle[0])
    x2 = np.max(rectangle[0])
    y1 = np.min(rectangle[1])
    y2 = np.max(rectangle[1])
    # Check that f_lower is allways lower than f_upper between x1 and x2 if no bounds are given or between the bounds
    #  of the f_lower and f_upper functions if they are given.
    if check:
        if bounds_lower is not None:
            if bounds_upper is not None:
                if not all(np.array(bounds_lower) == np.array(bounds_upper)):
                    raise ValueError('Bounds should be identical for both f_lower and f_upper')
                if not '<' in function_comparison(f1=f_lower, f2=f_upper,
                                                  x1=bounds_lower[0], x2=bounds_lower[1],
                                                  numpoints_check=numpoints_check):
                    raise RuntimeError('Function f_lower is not allways lower or equal to function f_upper within '
                                       'the domain defined by the functions bounds.')
            else:
                raise ValueError('Bounds are given for f_lower but not for f_upper')
        elif bounds_upper is not None:
            if bounds_lower is None:
                raise ValueError('Bounds are given for f_upper but not for f_lower')
            else:
                if not '<' in function_comparison(f1=f_lower, f2=f_upper,
                                                  x1=bounds_lower[0], x2=bounds_lower[1],
                                                  numpoints_check=numpoints_check):
                    raise RuntimeError('Function f_lower is not allways lower or equal to function f_upper within '
                                       'the domain defined by the functions bounds.')
        else:
            if not '<' in function_comparison(f1=f_lower, f2=f_upper,
                                              x1=x1, x2=x2, numpoints_check=numpoints_check):
                raise RuntimeError('Function f_lower is not allways lower or equal to function f_upper within '
                                   'the domain defined by x1 and x2.')
    if bounds_lower is None:
        raise NotImplementedError('Bounds should be given right now ...')
    else:
        if x2 < bounds_lower[0] or x1 > bounds_lower[1]:
            return (0.0, 0.0)
        if x1 < bounds_lower[0]:
            xmin = bounds_lower[0]
        else:
            xmin = x1
        if x2 > bounds_lower[1]:
            xmax = bounds_lower[1]
        else:
            xmax = x2
        def diff(x):
            flwx = f_lower(x)
            fupx = f_upper(x)
            minup = np.min([fupx, y2 * np.ones_like(fupx)], axis=0)
            maxlw = np.max([flwx, y1 * np.ones_like(flwx)], axis=0)
            zeros = np.zeros_like(fupx)
            upper = np.where(y2 >= flwx, np.where(y1 <= fupx, minup, zeros), zeros)
            lower = np.where(y1 <= fupx, np.where(y2 >= flwx, maxlw, zeros), zeros)
            return upper-lower
        return quad(diff, xmin, xmax)

def my_solid_angle(center, coords):
    """
    Helper method to calculate the solid angle of a set of coords from the
    center.

    Args:
        center:
            Center to measure solid angle from.
        coords:
            List of coords to determine solid angle.

    Returns:
        The solid angle.
    """
    o = np.array(center)
    r = [np.array(c) - o for c in coords]
    r.append(r[0])
    n = [np.cross(r[i + 1], r[i]) for i in range(len(r) - 1)]
    n.append(np.cross(r[1], r[0]))
    phi = 0.0
    for i in range(len(n) - 1):
        try:
            value = math.acos(-np.dot(n[i], n[i + 1]) / (np.linalg.norm(n[i]) * np.linalg.norm(n[i + 1])))
        except ValueError:
            mycos = -np.dot(n[i], n[i + 1]) / (np.linalg.norm(n[i]) * np.linalg.norm(n[i + 1]))
            if 0.999999999999 < mycos < 1.000000000001:
                value = math.acos(1.0)
            elif -0.999999999999 > mycos > -1.000000000001:
                value = math.acos(-1.0)
            else:
                raise SolidAngleError(mycos)
        phi += value
    return phi + (3 - len(r)) * math.pi


def vectorsToMatrix(aa, bb):
    """
    Performs the vector multiplication of the elements of two vectors, constructing the 3x3 matrix.
    :param aa: One vector of size 3
    :param bb: Another vector of size 3
    :return: A 3x3 matrix M composed of the products of the elements of aa and bb :
     M_ij = aa_i * bb_j
    """
    MM = np.zeros([3, 3], np.float)
    for ii in range(3):
        for jj in range(3):
            MM[ii, jj] = aa[ii] * bb[jj]
    return MM


def matrixTimesVector(MM, aa):
    """

    :param MM: A matrix of size 3x3
    :param aa: A vector of size 3
    :return: A vector of size 3 which is the product of the matrix by the vector
    """
    bb = np.zeros(3, np.float)
    for ii in range(3):
        bb[ii] = np.sum(MM[ii, :] * aa)
    return bb


def rotateCoords(coords, R):
    """
    Rotate the list of points using rotation matrix R
    :param coords: List of points to be rotated
    :param R: Rotation matrix
    :return: List of rotated points
    """
    newlist = list()
    for pp in coords:
        rpp = matrixTimesVector(R, pp)
        newlist.append(rpp)
    return newlist


def matrixMultiplication(AA, BB):
    """
    Performs the multiplication of two matrix of size 3x3
    :param AA: One matrix of size 3x3
    :param BB: Another matrix of size 3x3
    :return: A matrix of size 3x3
    """
    MM = np.zeros([3, 3], np.float)
    for ii in range(3):
        for jj in range(3):
            MM[ii, jj] = np.sum(AA[ii, :] * BB[:, jj])
    return MM


def changebasis(uu, vv, nn, pps):
    """
    For a list of points given in standard coordinates (in terms of e1, e2 and e3), returns the same list
    expressed in the basis (uu, vv, nn), which is supposed to be orthonormal.
    :param uu: First vector of the basis
    :param vv: Second vector of the basis
    :param nn: Third vector of the bais
    :param pps: List of points in basis (e1, e2, e3)
    :return: List of points in basis (uu, vv, nn)
    """
    MM = np.zeros([3, 3], np.float)
    for ii in range(3):
        MM[ii, 0] = uu[ii]
        MM[ii, 1] = vv[ii]
        MM[ii, 2] = nn[ii]
    PP = np.linalg.inv(MM)
    newpps = list()
    for pp in pps:
        newpps.append(matrixTimesVector(PP, pp))
    return newpps


def collinear(p1, p2, p3=None, tolerance=0.25):
    """
    Checks if the three points p1, p2 and p3 are collinear or not within a given tolerance. The collinearity is
    checked by computing the area of the triangle defined by the three points p1, p2 and p3. If the area of this
    triangle is less than (tolerance x largest_triangle), then the three points are considered collinear. The
    largest_triangle is defined as the right triangle whose legs are the two smallest distances between the three
     points ie, its area is : 0.5 x (min(|p2-p1|,|p3-p1|,|p3-p2|) x secondmin(|p2-p1|,|p3-p1|,|p3-p2|))
    :param p1: First point
    :param p2: Second point
    :param p3: Third point (origin [0.0, 0.0, 0.0 if not given])
    :param tolerance: Area tolerance for the collinearity test (0.25 gives about 0.125 deviation from the line)
    :return: True if the three points are considered as collinear within the given tolerance, False otherwise
    """
    if p3 is None:
        triangle_area = 0.5 * np.linalg.norm(np.cross(p1, p2))
        dist = np.sort([np.linalg.norm(p2 - p1), np.linalg.norm(p1), np.linalg.norm(p2)])
    else:
        triangle_area = 0.5 * np.linalg.norm(np.cross(p1 - p3, p2 - p3))
        dist = np.sort([np.linalg.norm(p2 - p1), np.linalg.norm(p3 - p1), np.linalg.norm(p3 - p2)])
    largest_triangle_area = 0.5 * dist[0] * dist[1]
    return triangle_area < tolerance * largest_triangle_area


def anticlockwise_sort(pps):
    """
    Sort a list of 2D points in anticlockwise order
    :param pps: List of points to be sorted
    :return: Sorted list of points
    """
    newpps = list()
    angles = np.zeros(len(pps), np.float)
    for ipp, pp in enumerate(pps):
        angles[ipp] = np.arctan2(pp[1], pp[0])
    iisorted = np.argsort(angles)
    for ii in range(len(pps)):
        newpps.append(pps[iisorted[ii]])
    return newpps


def anticlockwise_sort_indices(pps):
    """
    Returns the indices that would sort a list of 2D points in anticlockwise order
    :param pps: List of points to be sorted
    :return: Indices of the sorted list of points
    """
    angles = np.zeros(len(pps), np.float)
    for ipp, pp in enumerate(pps):
        angles[ipp] = np.arctan2(pp[1], pp[0])
    return np.argsort(angles)


def sort_separation(separation):
    if len(separation[0]) > len(separation[2]):
        return [sorted(separation[2]), sorted(separation[1]), sorted(separation[0])]
    return [sorted(separation[0]), sorted(separation[1]), sorted(separation[2])]


def separation_in_list(separation_indices, separation_indices_list):
    """
    Checks if the separation indices of a plane are already in the list
    :param separation_indices: list of separation indices (three arrays of integers)
    :param separation_indices_list: list of the list of separation indices to be compared to
    :return: True if the separation indices are already in the list, False otherwise
    """
    sorted_separation = sort_separation(separation_indices)
    for sep in separation_indices_list:
        if len(sep[1]) == len(sorted_separation[1]) and np.allclose(sorted_separation[1], sep[1]):
            return True
    return False


def is_anion_cation_bond(valences, ii, jj):
    """
    Checks if two given sites are an anion and a cation.
    :param valences: list of site valences
    :param ii: index of a site
    :param jj: index of another site
    :return: True if one site is an anion and the other is a cation (from the valences)
    """
    if valences == 'undefined':
        return True
    if valences[ii] == 0 or valences[jj] == 0:
        return True
    return (valences[ii] > 0 > valences[jj]) or (valences[jj] > 0 > valences[ii])


class Plane(object):
    """
    Class used to describe a plane
    """

    TEST_2D_POINTS = [np.array([0, 0], np.float),
                      np.array([1, 0], np.float),
                      np.array([0, 1], np.float),
                      np.array([-1, 0], np.float),
                      np.array([0, -1], np.float),
                      np.array([0, 2], np.float),
                      np.array([2, 0], np.float),
                      np.array([0, -2], np.float),
                      np.array([-2, 0], np.float),
                      np.array([1, 1], np.float),
                      np.array([2, 2], np.float),
                      np.array([-1, -1], np.float),
                      np.array([-2, -2], np.float),
                      np.array([1, 2], np.float),
                      np.array([1, -2], np.float),
                      np.array([-1, 2], np.float),
                      np.array([-1, -2], np.float),
                      np.array([2, 1], np.float),
                      np.array([2, -1], np.float),
                      np.array([-2, 1], np.float),
                      np.array([-2, -1], np.float)]

    def __init__(self, coefficients, p1=None, p2=None, p3=None):
        """
        Initializes a plane from the 4 coefficients a, b, c and d of ax + by + cz + d = 0
        :param coefficients: abcd coefficients of the plane
        """
        #Initializes the normal vector
        self.normal_vector = np.array([coefficients[0], coefficients[1], coefficients[2]], np.float)
        normv = np.linalg.norm(self.normal_vector)
        self.normal_vector /= normv
        nonzeros = np.argwhere(self.normal_vector != 0.0).flatten()
        zeros = list(set(range(3))-set(nonzeros))
        if len(nonzeros) == 0:
            raise ValueError("Normal vector is equal to 0.0")
        if self.normal_vector[nonzeros[0]] < 0.0:
            self.normal_vector = -self.normal_vector
            dd = -np.float(coefficients[3]) / normv
        else:
            dd = np.float(coefficients[3]) / normv
        self._coefficients = np.array([self.normal_vector[0],
                                      self.normal_vector[1],
                                      self.normal_vector[2],
                                      dd], np.float)
        self._crosses_origin = np.isclose(dd, 0.0, atol=1e-7, rtol=0.0)
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        #Initializes 3 points belonging to the plane (useful for some methods)
        if self.p1 is None:
            self.init_3points(nonzeros, zeros)
        self.vector_to_origin = dd * self.normal_vector

    def init_3points(self, nonzeros, zeros):
        if len(nonzeros) == 3:
            self.p1 = np.array([-self.d / self.a, 0.0, 0.0], np.float)
            self.p2 = np.array([0.0, -self.d / self.b, 0.0], np.float)
            self.p3 = np.array([0.0, 0.0, -self.d / self.c], np.float)
        elif len(nonzeros) == 2:
            self.p1 = np.zeros(3, np.float)
            self.p1[nonzeros[1]] = -self.d / self.coefficients[nonzeros[1]]
            self.p2 = np.array(self.p1)
            self.p2[zeros[0]] = 1.0
            self.p3 = np.zeros(3, np.float)
            self.p3[nonzeros[0]] = -self.d / self.coefficients[nonzeros[0]]
        elif len(nonzeros) == 1:
            self.p1 = np.zeros(3, np.float)
            self.p1[nonzeros[0]] = -self.d / self.coefficients[nonzeros[0]]
            self.p2 = np.array(self.p1)
            self.p2[zeros[0]] = 1.0
            self.p3 = np.array(self.p1)
            self.p3[zeros[1]] = 1.0

    def __str__(self):
        """
        String representation of the Plane object
        :return: String representation of the Plane object
        """
        outs = ['Plane object']
        outs.append('  => Normal vector : {nn}'.format(nn=self.normal_vector))
        outs.append('  => Equation of the plane ax + by + cz + d = 0')
        outs.append('     with a = {v}'.format(v=self._coefficients[0]))
        outs.append('          b = {v}'.format(v=self._coefficients[1]))
        outs.append('          c = {v}'.format(v=self._coefficients[2]))
        outs.append('          d = {v}'.format(v=self._coefficients[3]))
        return '\n'.join(outs)

    def is_in_plane(self, pp, dist_tolerance):
        """
        Determines if point pp is in the plane within the tolerance dist_tolerance
        :param pp: point to be tested
        :param dist_tolerance: tolerance on the distance to the plane within which point pp is considered in the plane
        :return: True if pp is in the plane, False otherwise
        """
        return np.abs(np.dot(self.normal_vector, pp) + self._coefficients[3]) <= dist_tolerance

    def is_same_plane_as(self, plane):
        """
        Checks whether the plane is identical to another Plane "plane"
        :param plane: Plane to be compared to
        :return: True if the two facets are identical, False otherwise
        """
        return np.allclose(self._coefficients, plane.coefficients)

    def is_in_list(self, plane_list):
        """
        Checks whether the plane is identical to one of the Planes in the plane_list list of Planes
        :param plane_list: List of Planes to be compared to
        :return: True if the plane is in the list, False otherwise
        """
        for plane in plane_list:
            if self.is_same_plane_as(plane):
                return True
        return False

    def indices_separate(self, points, dist_tolerance):
        """
        Returns three lists containing the indices of the points lying on one side of the plane, on the plane
        and on the other side of the plane. The dist_tolerance parameter controls the tolerance to which a point
        is considered to lie on the plane or not (distance to the plane)
        :param points: list of points
        :param dist_tolerance: tolerance to which a point is considered to lie on the plane
        or not (distance to the plane)
        :return: The lists of indices of the points on one side of the plane, on the plane and
        on the other side of the plane
        """
        side1 = list()
        inplane = list()
        side2 = list()
        for ip, pp in enumerate(points):
            if self.is_in_plane(pp, dist_tolerance):
                inplane.append(ip)
            else:
                if np.dot(pp + self.vector_to_origin, self.normal_vector) < 0.0:
                    side1.append(ip)
                else:
                    side2.append(ip)
        return [side1, inplane, side2]

    def distance_to_point(self, point):
        """
        Computes the absolute distance from the plane to the point
        :param point:
        :return:
        """
        return np.abs(np.dot(self.normal_vector, point) + self.d)

    def projectionpoints(self, pps):
        """
        Projects each points in the point list pps on plane and returns the list of projected points
        :param pps: List of points to project on plane
        :return: List of projected point on plane
        """
        return [pp - np.dot(pp - self.p1, self.normal_vector) * self.normal_vector for pp in pps]

    def orthonormal_vectors(self):
        """
        Returns a list of three orthogonal vectors, the two first being parallel to the plane and the
        third one is the normal vector of the plane
        :return: List of orthogonal vectors
        :raise: ValueError if all the coefficients are zero or if there is some other strange error
        """
        e1 = np.array((self.p2 - self.p1) / norm(self.p2 - self.p1))
        e2 = np.cross(self.normal_vector, e1)
        return [e1, e2, np.array(self.normal_vector)]

    def project_and_to2dim_ordered_indices(self, pps, plane_center='mean'):
        """
        Projects each points in the point list pps on plane and returns the indices that would sort the
        list of projected points in anticlockwise order
        :param pps: List of points to project on plane
        :return: List of indices that would sort the list of projected points
        """
        pp2d = self.project_and_to2dim(pps, plane_center)
        return anticlockwise_sort_indices(pp2d)

    def project_and_to2dim(self, pps, plane_center):
        """
        Projects the list of points pps to the plane and changes the basis from 3D to the 2D basis of the plane
        :param pps: List of points to be projected
        :return: :raise:
        """
        proj = self.projectionpoints(pps)
        [u1, u2, u3] = self.orthonormal_vectors()
        PP = np.array([[u1[0], u2[0], u3[0]],
                       [u1[1], u2[1], u3[1]],
                       [u1[2], u2[2], u3[2]]])
        xypps = list()
        for pp in proj:
            xyzpp = np.dot(pp, PP)
            xypps.append(xyzpp[0:2])
        if str(plane_center) == str('mean'):
            mean = np.zeros(2, np.float)
            for pp in xypps:
                mean += pp
            mean /= len(xypps)
            xypps = [pp - mean for pp in xypps]
        elif plane_center is not None:
            projected_plane_center = self.projectionpoints([plane_center])[0]
            xy_projected_plane_center = np.dot(projected_plane_center, PP)[0:2]
            xypps = [pp - xy_projected_plane_center for pp in xypps]
        return xypps

    def fit_error(self, points, fit='least_square_distance'):
        if fit == 'least_square_distance':
            return self.fit_least_square_distance_error(points)
        if fit == 'maximum_distance':
            return self.fit_maximum_distance_error(points)

    def fit_least_square_distance_error(self, points):
        return np.sum([self.distance_to_point(pp)**2.0 for pp in points])

    def fit_maximum_distance_error(self, points):
        return np.max([self.distance_to_point(pp) for pp in points])

    @property
    def coefficients(self):
        return np.copy(self._coefficients)

    @property
    def abcd(self):
        return self._coefficients[0], self._coefficients[1], self._coefficients[2], self._coefficients[3]

    @property
    def a(self):
        return self._coefficients[0]

    @property
    def b(self):
        return self._coefficients[1]

    @property
    def c(self):
        return self._coefficients[2]

    @property
    def d(self):
        return self._coefficients[3]

    @property
    def distance_to_origin(self):
        return self._coefficients[3]

    @property
    def crosses_origin(self):
        return self._crosses_origin

    @classmethod
    def from_2points_and_origin(cls, p1, p2):
        return cls.from_3points(p1, p2, np.zeros(3))

    @classmethod
    def from_3points(cls, p1, p2, p3):
        nn = np.cross(p1 - p3, p2 - p3)
        normal_vector = nn / norm(nn)
        nonzeros = np.argwhere(normal_vector != 0.0)
        if normal_vector[nonzeros[0, 0]] < 0.0:
            normal_vector = -normal_vector
        dd = - np.dot(normal_vector, p1)
        coefficients = np.array([normal_vector[0],
                                 normal_vector[1],
                                 normal_vector[2],
                                 dd], np.float)
        return cls(coefficients, p1=p1, p2=p2, p3=p3)

    @classmethod
    def from_npoints(cls, points, best_fit='least_square_distance'):
        if len(points) == 2:
            return cls.from_2points_and_origin(points[0], points[1])
        if len(points) == 3:
            return cls.from_3points(points[0], points[1], points[2])
        if best_fit == 'least_square_distance':
            return cls.from_npoints_least_square_distance(points)
        elif best_fit == 'maximum_distance':
            return cls.from_npoints_maximum_distance(points)

    @classmethod
    def from_npoints_least_square_distance(cls, points):
        mean_point = np.array([sum([pp[ii] for pp in points]) for ii in range(3)], np.float)
        mean_point /= len(points)
        AA = np.zeros((len(points), 3), np.float)
        for ii, pp in enumerate(points):
            for jj in range(3):
                AA[ii, jj] = pp[jj] - mean_point[jj]
        [UU, SS, Vt] = np.linalg.svd(AA)
        imin = np.argmin(SS)
        normal_vector = Vt[imin]
        nonzeros = np.argwhere(normal_vector != 0.0)
        if normal_vector[nonzeros[0, 0]] < 0.0:
            normal_vector = -normal_vector
        dd = - np.dot(normal_vector, mean_point)
        coefficients = np.array([normal_vector[0],
                                 normal_vector[1],
                                 normal_vector[2],
                                 dd], np.float)
        return cls(coefficients)

    @classmethod
    def perpendicular_bisector(cls, p1, p2):
        middle_point = 0.5*(p1+p2)
        normal_vector = p2-p1
        dd = -np.dot(normal_vector, middle_point)
        return cls(np.array([normal_vector[0],
                             normal_vector[1],
                             normal_vector[2],
                             dd], np.float))

    @classmethod
    def from_npoints_maximum_distance(cls, points):
        convex_hull = ConvexHull(points)
        heights = []
        ipoints_heights = []
        for isimplex, simplex in enumerate(convex_hull.simplices):
            cc = convex_hull.equations[isimplex]
            plane = Plane.from_coefficients(cc[0], cc[1], cc[2], cc[3])
            distances = [plane.distance_to_point(pp) for pp in points]
            ipoint_height = np.argmax(distances)
            heights.append(distances[ipoint_height])
            ipoints_heights.append(ipoint_height)
        imin_height = np.argmin(heights)
        normal_vector = convex_hull.equations[imin_height, 0:3]
        cc = convex_hull.equations[imin_height]
        highest_point = points[ipoints_heights[imin_height]]
        middle_point = (Plane.from_coefficients(cc[0], cc[1], cc[2], cc[3]).projectionpoints([highest_point])[0] +
                        highest_point) / 2
        dd = - np.dot(normal_vector, middle_point)
        return cls(np.array([normal_vector[0],
                             normal_vector[1],
                             normal_vector[2],
                             dd], np.float))

    @classmethod
    def from_coefficients(cls, a, b, c, d):
        return cls(np.array([a, b, c, d], np.float))