"""This module contains some math utils that are used in the chemenv package."""

from __future__ import annotations

import math
import operator
from functools import reduce
from typing import TYPE_CHECKING

import numpy as np
from scipy.special import erf

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

##############################################################
# Cartesian product of lists #################################
##############################################################


def _append_es2sequences(sequences, es):
    result = []
    if not sequences:
        for e in es:
            result.append([e])
    else:
        for e in es:
            result += [[*seq, e] for seq in sequences]
    return result


def _cartesian_product(lists):
    """
    given a list of lists,
    returns all the possible combinations taking one element from each list
    The list does not have to be of equal length.
    """
    return reduce(_append_es2sequences, lists, [])


def prime_factors(n: int) -> list[int]:
    """Lists prime factors of a given natural integer, from greatest to smallest.

    Args:
        n: Natural integer

    Returns:
        list of all prime factors of the given natural n.
    """
    idx = 2
    while idx <= math.sqrt(n):
        if n % idx == 0:
            lst = prime_factors(n // idx)
            lst.append(idx)
            return lst
        idx += 1
    return [n]  # n is prime


def _factor_generator(n: int) -> dict[int, int]:
    """From a given natural integer, returns the prime factors and their multiplicity.

    Args:
        n: Natural integer
    """
    p = prime_factors(n)
    factors: dict[int, int] = {}
    for p1 in p:
        try:
            factors[p1] += 1
        except KeyError:
            factors[p1] = 1
    return factors


def divisors(n):
    """From a given natural integer, returns the list of divisors in ascending order.

    Args:
        n: Natural integer

    Returns:
        List of divisors of n in ascending order.
    """
    factors = _factor_generator(n)
    _divisors = []
    exponents = [[k**x for x in range(factors[k] + 1)] for k in list(factors)]
    factors = _cartesian_product(exponents)
    for factor in factors:
        _divisors.append(reduce(operator.mul, factor, 1))
    _divisors.sort()
    return _divisors


def get_center_of_arc(p1, p2, radius):
    """
    Args:
        p1:
        p2:
        radius:
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


def get_linearly_independent_vectors(vectors: list[ArrayLike]) -> list[np.ndarray]:
    """
    Args:
        vectors (list[ArrayLike]): List of vectors.
    """
    independent_vectors: list[np.ndarray] = []
    for vector in vectors:
        if np.any(vector != 0):
            if len(independent_vectors) == 0:
                independent_vectors.append(np.array(vector))
            elif len(independent_vectors) == 1:
                rank = np.linalg.matrix_rank(np.array([independent_vectors[0], vector, [0, 0, 0]]))
                if rank == 2:
                    independent_vectors.append(np.array(vector))
            elif len(independent_vectors) == 2:
                mm = np.array([independent_vectors[0], independent_vectors[1], vector])
                if np.linalg.det(mm) != 0:
                    independent_vectors.append(np.array(vector))
        if len(independent_vectors) == 3:
            break
    return independent_vectors


def scale_and_clamp(xx, edge0, edge1, clamp0, clamp1):
    """
    Args:
        xx:
        edge0:
        edge1:
        clamp0:
        clamp1:
    """
    return np.clip((xx - edge0) / (edge1 - edge0), clamp0, clamp1)


# Step function based on the cumulative distribution function of the normal law
def normal_cdf_step(xx, mean, scale):
    """
    Args:
        xx:
        mean:
        scale:
    """
    return 0.5 * (1.0 + erf((xx - mean) / (np.sqrt(2.0) * scale)))


# SMOOTH STEP FUNCTIONS
# Set of smooth step functions that allow to smoothly go from y = 0.0 (1.0) to y = 1.0 (0.0) by changing x
# from 0.0 to 1.0 respectively when inverse is False (True).
# (except if edges is given in which case a the values are first scaled and clamped to the interval given by edges)
# The derivative at x = 0.0 and x = 1.0 have to be 0.0


def smoothstep(xx, edges=None, inverse=False):
    """
    Args:
        xx:
        edges:
        inverse:
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
    Args:
        xx:
        edges:
        inverse:
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
    Args:
        xx:
        edges:
        inverse:
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
    Args:
        xx:
        edges:
        inverse:
    """
    return smoothstep(xx, edges=edges, inverse=inverse)


def powern_parts_step(xx, edges=None, inverse=False, nn=2):
    """
    Args:
        xx:
        edges:
        inverse:
        nn:
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
    Args:
        xx:
        edges:
        nn:
    """
    if edges is None:
        aa = 1.0 / np.power(-1.0, nn)
        return aa * np.power(xx - 1.0, nn)
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return powern_decreasing(xx_scaled_and_clamped, nn=nn)


def power2_decreasing_exp(xx, edges=None, alpha=1.0):
    """
    Args:
        xx:
        edges:
        alpha:
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
    Args:
        xx:
        edges:
        prefactor:
    """
    if edges is None:
        aa = 1.0 / np.power(-1.0, 2) if prefactor is None else prefactor
        return -aa * np.power(xx - 1.0, 2) * np.tan((xx - 1.0) * np.pi / 2.0)

    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return power2_tangent_decreasing(xx_scaled_and_clamped, prefactor=prefactor)


def power2_inverse_decreasing(xx, edges=None, prefactor=None):
    """
    Args:
        xx:
        edges:
        prefactor:
    """
    if edges is None:
        aa = 1.0 / np.power(-1.0, 2) if prefactor is None else prefactor
        return np.where(np.isclose(xx, 0.0), aa * float("inf"), aa * np.power(xx - 1.0, 2) / xx)
        # return aa * np.power(xx-1.0, 2) / xx if xx != 0 else aa * float("inf")
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return power2_inverse_decreasing(xx_scaled_and_clamped, prefactor=prefactor)


def power2_inverse_power2_decreasing(xx, edges=None, prefactor=None):
    """
    Args:
        xx:
        edges:
        prefactor:
    """
    if edges is None:
        aa = 1.0 / np.power(-1.0, 2) if prefactor is None else prefactor
        return np.where(
            np.isclose(xx, 0.0),
            aa * float("inf"),
            aa * np.power(xx - 1.0, 2) / xx**2.0,
        )
    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return power2_inverse_power2_decreasing(xx_scaled_and_clamped, prefactor=prefactor)


def power2_inverse_powern_decreasing(xx, edges=None, prefactor=None, powern=2.0):
    """
    Args:
        xx:
        edges:
        prefactor:
        powern:
    """
    if edges is None:
        aa = 1.0 / np.power(-1.0, 2) if prefactor is None else prefactor
        return aa * np.power(xx - 1.0, 2) / xx**powern

    xx_scaled_and_clamped = scale_and_clamp(xx, edges[0], edges[1], 0.0, 1.0)
    return power2_inverse_powern_decreasing(xx_scaled_and_clamped, prefactor=prefactor, powern=powern)
