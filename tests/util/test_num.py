from __future__ import annotations

import pytest

from pymatgen.util.num import maxloc, min_max_indexes, minloc, round_to_sigfigs


def test_minloc():
    assert minloc(range(3)) == 0


def test_maxloc():
    assert maxloc([1, 3, 2, 3]) == 1


def test_min_max_indexes():
    val = ["b", "a", "m", "z", "y"]
    min_ind, max_ind = min_max_indexes(val)
    assert min_ind == 1
    assert max_ind == 3


def test_round():
    vals = [424.2425, 2.3425356, 0.000042535636653, 0.23, 2.468e6, 0, -1.392156]
    sig_figs = range(1, 6)
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
        for s, sig in enumerate(sig_figs):
            assert round_to_sigfigs(val, sig) == rounded_vals[v][s]
    with pytest.raises(ValueError, match="Number of significant figures must be positive"):
        round_to_sigfigs(3.5, -2)
    with pytest.raises(TypeError, match="Number of significant figures must be integer"):
        round_to_sigfigs(3.5, 3.5)
