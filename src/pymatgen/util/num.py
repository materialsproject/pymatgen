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
    """Given a symmetric matrix in upper triangular matrix form as a flat array:
        [A_xx, A_yy, A_zz, A_xy, A_xz, A_yz].

    This will generate the full matrix:
        [[A_xx, A_xy, A_xz], [A_xy, A_yy, A_yz], [A_xz, A_yz, A_zz]].

    Args:
        val (ArrayLike): Flattened upper triangular elements.

    Returns:
        NDArray: The symmetric matrix.
    """
    idx = [0, 3, 4, 1, 5, 2]  # TODO: can this be less hard-coded (magic)?
    val = np.asarray(val)[idx]
    mask = ~np.tri(3, k=-1, dtype=bool)
    out = np.zeros((3, 3), dtype=val.dtype)
    out[mask] = val

    out.T[mask] = val
    return out
