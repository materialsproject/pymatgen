"""This module provides utilities for basic math operations."""

from __future__ import annotations

import numpy as np


def abs_cap(val, max_abs_val=1):
    """
    Returns the value with its absolute value capped at max_abs_val.
    Particularly useful in passing values to trigonometric functions where
    numerical errors may result in an argument > 1 being passed in.

    Args:
        val (float): Input value.
        max_abs_val (float): The maximum absolute value for val. Defaults to 1.

    Returns:
        val if abs(val) < 1 else sign of val * max_abs_val.
    """
    return max(min(val, max_abs_val), -max_abs_val)


def minloc(seq):
    """
    Return the index of the (first) minimum in seq.

    >>> assert minloc(range(3)) == 0
    """
    return min(enumerate(seq), key=lambda s: s[1])[0]


def maxloc(seq):
    """
    Return the index of the (first) maximum in seq.

    >>> assert maxloc([1,3,2,3]) == 1
    """
    return max(enumerate(seq), key=lambda s: s[1])[0]


def min_max_indexes(seq):
    """
    Uses enumerate, max, and min to return the indices of the values
    in a list with the maximum and minimum value.

    Args:
        seq: A sequence of numbers.
    """
    lst = sorted(enumerate(seq), key=lambda tup: tup[1])
    return lst[0][0], lst[-1][0]


def strictly_increasing(values):
    """True if values are strictly increasing."""
    return all(x < y for x, y in zip(values, values[1:]))


def strictly_decreasing(values):
    """True if values are strictly decreasing."""
    return all(x > y for x, y in zip(values, values[1:]))


def round_to_sigfigs(num, sig_figs):
    """
    Rounds a number rounded to a specific number of significant
    figures instead of to a specific precision.
    """
    if not isinstance(sig_figs, int):
        raise TypeError("Number of significant figures must be integer")

    if sig_figs < 1:
        raise ValueError("Number of significant figures must be positive")

    if num == 0:
        return num

    prec = int(sig_figs - np.ceil(np.log10(np.absolute(num))))
    return round(num, prec)


def make_symmetric_matrix_from_upper_tri(val):
    """
    Given a symmetric matrix in upper triangular matrix form as flat array indexes as:
    [A_xx,A_yy,A_zz,A_xy,A_xz,A_yz]
    This will generate the full matrix:
    [[A_xx,A_xy,A_xz],[A_xy,A_yy,A_yz],[A_xz,A_yz,A_zz].
    """
    idx = [0, 3, 4, 1, 5, 2]
    val = np.array(val)[idx]
    mask = ~np.tri(3, k=-1, dtype=bool)
    out = np.zeros((3, 3), dtype=val.dtype)
    out[mask] = val
    # pylint: disable=E1137
    out.T[mask] = val
    return out
