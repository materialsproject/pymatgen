from __future__ import annotations

import math

import numpy as np
import pytest

from pymatgen.util.num import make_symmetric_matrix_from_upper_tri, round_to_sigfigs


class TestRoundToSigfigs:
    def test_round_to_sigfigs(self):
        vals = [424.2425, 2.3425356, 0.000042535636653, 0.23, 2.468e6, 0, -1.392156]
        rounded_vals = [
            [400.0, 420.0, 424.0, 424.2, 424.24],
            [2.0, 2.3, 2.34, 2.343, 2.3425],
            [4e-5, 4.3e-5, 4.25e-5, 4.254e-5, 4.2536e-5],
            [0.2, 0.23, 0.23, 0.23, 0.23],
            [2e6, 2.5e6, 2.47e6, 2.468e6, 2.468e6],
            [0, 0, 0, 0, 0],
            [-1, -1.4, -1.39, -1.392, -1.3922],
        ]

        for idx_val, val in enumerate(vals):
            for idx_sig, sig in enumerate(range(1, 6)):
                assert math.isclose(round_to_sigfigs(val, sig), rounded_vals[idx_val][idx_sig])

    def test_exceptions(self):
        with pytest.raises(ValueError, match="Number of significant figures must be positive"):
            round_to_sigfigs(3.5, -2)

        with pytest.raises(TypeError, match="Number of significant figures must be integer"):
            round_to_sigfigs(3.5, 3.5)


def test_make_symmetric_matrix_from_upper_tri():
    # Test regular array
    val = [1, 2, 3, 4, 5, 6]  # A_xx, A_yy, A_zz, A_xy, A_xz, A_yz
    expected = np.array(
        [
            [1, 4, 5],  # A_xx, A_xy, A_xz
            [4, 2, 6],  # A_xy, A_yy, A_yz
            [5, 6, 3],  # A_xz, A_yz, A_zz
        ]
    )
    result = make_symmetric_matrix_from_upper_tri(val)
    assert isinstance(result, np.ndarray)
    np.testing.assert_array_equal(result, expected)

    # Test all zeros
    val = [0, 0, 0, 0, 0, 0]
    expected = np.zeros((3, 3))
    result = make_symmetric_matrix_from_upper_tri(val)
    np.testing.assert_array_equal(result, expected)

    # Test same elements
    val = [1, 1, 1, 1, 1, 1]
    expected = np.array(
        [
            [1, 1, 1],
            [1, 1, 1],
            [1, 1, 1],
        ]
    )
    result = make_symmetric_matrix_from_upper_tri(val)
    np.testing.assert_array_equal(result, expected)

    # Test too few elements (should be >= 6)
    with pytest.raises(IndexError, match="index 5 is out of bounds for axis 0 with size 5"):
        make_symmetric_matrix_from_upper_tri([1, 2, 3, 4, 5])
