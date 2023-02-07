# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import random
import unittest

import pytest

from pymatgen.util.num import abs_cap, min_max_indexes, round_to_sigfigs

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "9/25/14"


class FuncTestCase(unittest.TestCase):
    def test_abs_cap(self):
        assert abs_cap(1.000000001) == 1.0
        assert abs_cap(-1.000000001) == -1.0

        v = random.uniform(-1, 1)
        assert abs_cap(v) == v

        assert abs_cap(1.000000001, 2) == 1.000000001
        assert abs_cap(-2.000000001, 2) == -2.0

    def test_min_max_indexes(self):
        val = ["b", "a", "m", "z", "y"]
        min_ind, max_ind = min_max_indexes(val)
        assert min_ind == 1
        assert max_ind == 3

    def test_round(self):
        vals = [424.2425, 2.3425356, 0.000042535636653, 0.23, 2.468e6, 0, -1.392156]
        sigfigs = range(1, 6)
        rounded_vals = [
            [400.0, 420.0, 424.0, 424.2, 424.24],
            [2.0, 2.3, 2.34, 2.343, 2.3425],
            [4e-5, 4.3e-5, 4.25e-5, 4.254e-5, 4.2536e-5],
            [0.2, 0.23, 0.23, 0.23, 0.23],
            [2e6, 2.5e6, 2.47e6, 2.468e6, 2.468e6],
            [0, 0, 0, 0, 0],
            [-1, -1.4, -1.39, -1.392, -1.3922],
        ]

        for v, val in enumerate(vals):
            for s, sig in enumerate(sigfigs):
                assert round_to_sigfigs(val, sig) == rounded_vals[v][s]
        with pytest.raises(ValueError):
            round_to_sigfigs(3.5, -2)
        with pytest.raises(TypeError):
            round_to_sigfigs(3.5, 3.5)


if __name__ == "__main__":
    unittest.main()
