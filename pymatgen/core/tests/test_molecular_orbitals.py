from __future__ import annotations

import pytest

from pymatgen.core.molecular_orbitals import MolecularOrbitals
from pymatgen.util.testing import PymatgenTest

test_case = MolecularOrbitals("NaCl")


class MolecularOrbitalTestCase(PymatgenTest):
    def test_max_electronegativity(self):
        test_elec_neg = 2.23
        assert test_elec_neg == test_case.max_electronegativity()

    def test_aos_as_list(self):
        test_list = [
            ["Cl", "1s", -100.369229],
            ["Na", "1s", -37.719975],
            ["Cl", "2s", -9.187993],
            ["Cl", "2p", -7.039982],
            ["Na", "2s", -2.063401],
            ["Na", "2p", -1.060636],
            ["Cl", "3s", -0.754458],
            ["Cl", "3p", -0.32038],
            ["Na", "3s", -0.103415],
        ]
        assert test_list == test_case.aos_as_list()

    def test_obtain_band_edges(self):
        test_edges = {
            "HOMO": ["Cl", "3p", -0.32038],
            "LUMO": ["Na", "3s", -0.103415],
            "metal": False,
        }
        for key, val in test_edges.items():
            assert val == test_case.obtain_band_edges()[key]

    # test for raising ValueError for fractional composition
    def test_fractional_compositions(self):
        with pytest.raises(ValueError, match="composition subscripts must be integers"):
            MolecularOrbitals("Na0.5Cl0.5")
