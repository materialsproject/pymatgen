# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module contains some math utils that are used in the chemenv package.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

from functools import reduce
from math import sqrt

import numpy as np
from scipy.special import erf

##############################################################
# cartesian product of lists ##################################
##############################################################


def _append_es2sequences(sequences, es):
    result = []
    if not sequences:
        for e in es:
            result.append([e])
    else:
        for e in es:
            result += [seq + [e] for seq in sequences]
    return result


def _cartesian_product(lists):
    """
    given a list of lists,
    returns all the possible combinations taking one element from each list
    The list does not have to be of equal length
    """
    return reduce(_append_es2sequences, lists, [])


def prime_factors(n):
    """Lists prime factors of a given natural integer, from greatest to smallest
    :param n: Natural integer
    :rtype : list of all prime factors of the given natural n
    """
    i = 2
    while i <= sqrt(n):
        if n % i == 0:
            l = prime_factors(n / i)
            l.append(i)
            return l
        i += 1
    return [n]  # n is prime


def _factor_generator(n):
    """
    From a given natural integer, returns the prime factors and their multiplicity
    :param n: Natural integer
    :return:
    """
    p = prime_factors(n)
    factors = {}
    for p1 in p:
        try:
            factors[p1] += 1
        except KeyError:
            factors[p1] = 1
    return factors


def divisors(n):
    """
    From a given natural integer, returns the list of divisors in ascending order
    :param n: Natural integer
    :return: List of divisors of n in ascending order
    """
    factors = _factor_generator(n)
    _divisors = []
    listexponents = [[k ** x for x in range(0, factors[k] + 1)] for k in list(factors.keys())]
    listfactors = _cartesian_product(listexponents)
    for f in listfactors:
        _divisors.append(reduce(lambda x, y: x * y, f, 1))
    _divisors.sort()
    return _divisors


def get_center_of_arc(p1, p2, radius):
    """
    :param p1:
    :param p2:
    :param radius:
    :return:
    """
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]
    dd = np.sqrt(dx * dx + dy * dy)
    radical = np.power((radius / dd), 2) - 0.25
    if radical < 0:
        raise ValueError("Impossible to find center of arc because the arc is ill-defined")
    tt = np.sqrt(radical)
    if radius > 0:
        tt = -tt
    return (p1[0] + p2[0]) / 2 - tt * dy, (p1[1] + p2[1]) / 2 + tt * dx


def get_linearly_independent_vectors(vectors_list):
    """
    :param vectors_list:
    :return:
    """
    independent_vectors_list = []
    for vector in vectors_list:
        if np.any(vector != 0):
            if len(independent_vectors_list) == 0:
                independent_vectors_list.append(np.array(vector))
            elif len(independent_vectors_list) == 1:
                rank = np.linalg.matrix_rank(np.array([independent_vectors_list[0], vector, [0, 0, 0]]))
                if rank == 2:
                    independent_vectors_list.append(np.array(vector))
            elif len(independent_vectors_list) == 2:
                mm = np.array([independent_vectors_list[0], independent_vectors_list[1], vector])
                if np.linalg.det(mm) != 0:
                    independent_vectors_list.append(np.array(vector))
        if len(independent_vectors_list) == 3:
            break
    return independent_vectors_list


def scale_and_clamp(xx, edge0, edge1, clamp0, clamp1):
    """
    :param xx:
    :param edge0:
    :param edge1:
    :param clamp0:
    :param clamp1:
    :return:
    """
    return np.clip((xx - edge0) / (edge1 - edge0), clamp0, clamp1)


# Step function based on the cumulative distribution function of the normal law
def normal_cdf_step(xx, mean, scale):
    """
    :param xx:
    :param mean:
    :param scale:
    :return:
    """
    return 0.5 * (1.0 + erf((xx - mean) / (np.sqrt(2.0) * scale)))


# SMOOTH STEP FUNCTIONS
# Set of smooth step functions that allow to smoothly go from y = 0.0 (1.0) to y = 1.0 (0.0) by changing x
# from 0.0 to 1.0 respectively when inverse is False (True).
# (except if edges is given in which case a the values are first scaled and clamped to the interval given by edges)
# The derivative at x = 0.0 and x = 1.0 have to be 0.0


def smoothstep(xx, edges=None, inverse=False):
    """
    :param xx:
    :param edges:
    :param inverse:
    :return:
    """
    if edges is None:
        xx_clipped = np.clip(xx, 0.0, 1.0)
        if inverse:
            return 1.0 - xx_clipped * xx_clipped * (3.0 - 2.0 * xx_clipped)
        return xx_clipped * xx_clipped * (3.0 - 2.0 * xx_clipped)
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return smoothstep(xx_scaled_and_clamped, inverse=inverse)


def smootherstep(xx, edges=None, inverse=False):
    """
    :param xx:
    :param edges:
    :param inverse:
    :return:
    """
    if edges is None:
        xx_clipped = np.clip(xx, 0.0, 1.0)
        if inverse:
            return 1.0 - xx_clipped * xx_clipped * xx_clipped * (xx_clipped * (xx_clipped * 6 - 15) + 10)
        return xx_clipped * xx_clipped * xx_clipped * (xx_clipped * (xx_clipped * 6 - 15) + 10)
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return smootherstep(xx_scaled_and_clamped, inverse=inverse)


def cosinus_step(xx, edges=None, inverse=False):
    """
    :param xx:
    :param edges:
    :param inverse:
    :return:
    """
    if edges is None:
        xx_clipped = np.clip(xx, 0.0, 1.0)
        if inverse:
            return (np.cos(xx_clipped * np.pi) + 1.0) / 2.0
        return 1.0 - (np.cos(xx_clipped * np.pi) + 1.0) / 2.0
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return cosinus_step(xx_scaled_and_clamped, inverse=inverse)


def power3_step(xx, edges=None, inverse=False):
    """
    :param xx:
    :param edges:
    :param inverse:
    :return:
    """
    return smoothstep(xx, edges=edges, inverse=inverse)


def powern_parts_step(xx, edges=None, inverse=False, nn=2):
    """
    :param xx:
    :param edges:
    :param inverse:
    :param nn:
    :return:
    """
    if edges is None:
        aa = np.power(0.5, 1.0 - nn)
        xx_clipped = np.clip(xx, 0.0, 1.0)
        if np.mod(nn, 2) == 0:
            if inverse:
                return 1.0 - np.where(
                    xx_clipped < 0.5,
                    aa * np.power(xx_clipped, nn),
                    1.0 - aa * np.power(xx_clipped - 1.0, nn),
                )
            return np.where(
                xx_clipped < 0.5,
                aa * np.power(xx_clipped, nn),
                1.0 - aa * np.power(xx_clipped - 1.0, nn),
            )
        if inverse:
            return 1.0 - np.where(
                xx_clipped < 0.5,
                aa * np.power(xx_clipped, nn),
                1.0 + aa * np.power(xx_clipped - 1.0, nn),
            )
        return np.where(
            xx_clipped < 0.5,
            aa * np.power(xx_clipped, nn),
            1.0 + aa * np.power(xx_clipped - 1.0, nn),
        )
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return powern_parts_step(xx_scaled_and_clamped, inverse=inverse, nn=nn)


# FINITE DECREASING FUNCTIONS
# Set of decreasing functions that allow to smoothly go from y = 1.0 to y = 0.0 by changing x from 0.0 to 1.0
# The derivative at x = 1.0 has to be 0.0


def powern_decreasing(xx, edges=None, nn=2):
    """
    :param xx:
    :param edges:
    :param nn:
    :return:
    """
    if edges is None:
        aa = 1.0 / np.power(-1.0, nn)
        return aa * np.power(xx - 1.0, nn)
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return powern_decreasing(xx_scaled_and_clamped, nn=nn)


def power2_decreasing_exp(xx, edges=None, alpha=1.0):
    """
    :param xx:
    :param edges:
    :param alpha:
    :return:
    """
    if edges is None:
        aa = 1.0 / np.power(-1.0, 2)
        return aa * np.power(xx - 1.0, 2) * np.exp(-alpha * xx)
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return power2_decreasing_exp(xx_scaled_and_clamped, alpha=alpha)


# INFINITE TO FINITE DECREASING FUNCTIONS
# Set of decreasing functions that allow to smoothly go from y = + Inf to y = 0.0 by changing x from 0.0 to 1.0
# The derivative at x = 1.0 has to be 0.0


def power2_tangent_decreasing(xx, edges=None, prefactor=None):
    """
    :param xx:
    :param edges:
    :param prefactor:
    :return:
    """
    if edges is None:
        if prefactor is None:
            aa = 1.0 / np.power(-1.0, 2)
        else:
            aa = prefactor
        return -aa * np.power(xx - 1.0, 2) * np.tan((xx - 1.0) * np.pi / 2.0)  # pylint: disable=E1130

    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return power2_tangent_decreasing(xx_scaled_and_clamped, prefactor=prefactor)


def power2_inverse_decreasing(xx, edges=None, prefactor=None):
    """
    :param xx:
    :param edges:
    :param prefactor:
    :return:
    """
    if edges is None:
        if prefactor is None:
            aa = 1.0 / np.power(-1.0, 2)
        else:
            aa = prefactor
        return np.where(np.isclose(xx, 0.0), aa * float("inf"), aa * np.power(xx - 1.0, 2) / xx)
        # return aa * np.power(xx-1.0, 2) / xx if xx != 0 else aa * float("inf")
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return power2_inverse_decreasing(xx_scaled_and_clamped, prefactor=prefactor)


def power2_inverse_power2_decreasing(xx, edges=None, prefactor=None):
    """
    :param xx:
    :param edges:
    :param prefactor:
    :return:
    """
    if edges is None:
        if prefactor is None:
            aa = 1.0 / np.power(-1.0, 2)
        else:
            aa = prefactor
        return np.where(
            np.isclose(xx, 0.0),
            aa * float("inf"),
            aa * np.power(xx - 1.0, 2) / xx ** 2.0,
        )
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return power2_inverse_power2_decreasing(xx_scaled_and_clamped, prefactor=prefactor)


def power2_inverse_powern_decreasing(xx, edges=None, prefactor=None, powern=2.0):
    """
    :param xx:
    :param edges:
    :param prefactor:
    :param powern:
    :return:
    """
    if edges is None:
        if prefactor is None:
            aa = 1.0 / np.power(-1.0, 2)
        else:
            aa = prefactor
        return aa * np.power(xx - 1.0, 2) / xx ** powern

    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return power2_inverse_powern_decreasing(xx_scaled_and_clamped, prefactor=prefactor, powern=powern)
