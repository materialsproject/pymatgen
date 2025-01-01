from __future__ import annotations

import math
import re

import numpy as np
import pytest
from numpy.testing import assert_array_equal

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


class TestMakeSymmetricMatrixFromUpperTri:
    def test_convert(self):
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
        assert_array_equal(result, expected)

    def test_invalid_val_length(self):
        for length in [x for x in range(10) if x != 6]:
            with pytest.raises(ValueError, match=re.escape(f"Expect val of length 6, got ({length},)")):
                make_symmetric_matrix_from_upper_tri(range(length))
