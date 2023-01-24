# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import unittest

from pytest import approx

from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.structure_prediction.dopant_predictor import (
    get_dopants_from_shannon_radii,
    get_dopants_from_substitution_probabilities,
)
from pymatgen.core.periodic_table import Species
from pymatgen.core.structure import Structure


class DopantPredictionTest(unittest.TestCase):
    def setUp(self):
        self.tin_dioxide = Structure(
            [3.24, 0, 0, 0, 4.83, 0, 0, 0, 4.84],
            ["O", "O", "O", "O", "Sn", "Sn"],
            [
                [0.5, 0.19, 0.80],
                [0.5, 0.80, 0.19],
                [0, 0.30, 0.30],
                [0, 0.69, 0.69],
                [0.5, 0.50, 0.50],
                [0, 0, 0],
            ],
        )
        self.tin_dioxide.add_oxidation_state_by_element({"Sn": 4, "O": -2})

    def test_dopants_from_substitution_probabilities(self):
        dopants = get_dopants_from_substitution_probabilities(self.tin_dioxide, num_dopants=5)

        assert "n_type" in dopants
        assert "p_type" in dopants

        assert len(dopants["n_type"]) <= 5

        assert len(dopants["p_type"]) <= 5

        assert dopants["n_type"][0]["probability"] == approx(0.06692682583342474)
        assert dopants["n_type"][0]["dopant_species"] == Species("F", -1)
        assert dopants["n_type"][0]["original_species"] == Species("O", -2)

        assert dopants["p_type"][0]["probability"] == approx(0.023398867249112935)
        assert dopants["p_type"][0]["dopant_species"] == Species("Co", 2)
        assert dopants["p_type"][0]["original_species"] == Species("Sn", 4)

        # test oxidation sign matching
        dopants = get_dopants_from_substitution_probabilities(self.tin_dioxide, num_dopants=15, match_oxi_sign=False)
        assert dopants["n_type"][14]["dopant_species"] == Species("Li", 1)
        assert dopants["n_type"][14]["original_species"] == Species("O", -2)

        dopants = get_dopants_from_substitution_probabilities(self.tin_dioxide, num_dopants=15, match_oxi_sign=True)
        assert dopants["n_type"][14]["dopant_species"] != Species("Li", 1)

    def test_dopants_from_shannon_radii(self):
        bonded_structure = CrystalNN().get_bonded_structure(self.tin_dioxide)

        dopants = get_dopants_from_shannon_radii(bonded_structure, num_dopants=5)

        assert "n_type" in dopants
        assert "p_type" in dopants

        assert len(dopants["n_type"]) <= 5
        assert len(dopants["p_type"]) <= 5

        assert dopants["n_type"][0]["radii_diff"] == approx(0.04)
        assert dopants["n_type"][0]["dopant_species"] == Species("U", 6)
        assert dopants["n_type"][0]["original_species"] == Species("Sn", 4)

        assert dopants["p_type"][0]["radii_diff"] == approx(0.0)
        assert dopants["p_type"][0]["dopant_species"] == Species("Ni", 2)
        assert dopants["p_type"][0]["original_species"] == Species("Sn", 4)
