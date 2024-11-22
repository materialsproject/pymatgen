"""This module provides utilities for basic math operations."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray


def round_to_sigfigs(num: float, sig_figs: int) -> float:
    """Rounds a number to a specific number of significant
    figures instead of to a specific precision.

    Returns:
        float: rounded number.
    """
    if not isinstance(sig_figs, int):
        raise TypeError("Number of significant figures must be integer")

    if sig_figs < 1:
        raise ValueError("Number of significant figures must be positive")

    if num == 0:
        return num

    prec = int(sig_figs - np.ceil(np.log10(np.absolute(num))))
    return round(num, prec)


def make_symmetric_matrix_from_upper_tri(val: ArrayLike) -> NDArray:
    """Construct a 3x3 symmetric matrix from its upper triangular
        elements in flat array form:
        [A_xx, A_yy, A_zz, A_xy, A_xz, A_yz].

    To the full symmetric matrix:
        [[A_xx, A_xy, A_xz],
         [A_xy, A_yy, A_yz],
         [A_xz, A_yz, A_zz]]

    Args:
        val (ArrayLike): Flattened upper triangular elements.

    Returns:
        NDArray: The symmetric matrix.
    """
    # Check shape of input array, this function is designed for 3x3 array only
    if (array := np.asarray(val)).shape != (6,):
        raise ValueError(f"Expect val of length 6, got {array.shape}")

    return np.array(
        [
            [array[0], array[3], array[4]],
            [array[3], array[1], array[5]],
            [array[4], array[5], array[2]],
        ]
    )
