# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest

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

        self.assertTrue("n_type" in dopants)
        self.assertTrue("p_type" in dopants)

        self.assertTrue(len(dopants["n_type"]) <= 5)

        self.assertTrue(len(dopants["p_type"]) <= 5)

        self.assertAlmostEqual(dopants["n_type"][0]["probability"], 0.06692682583342474)
        self.assertEqual(dopants["n_type"][0]["dopant_species"], Species("F", -1))
        self.assertEqual(dopants["n_type"][0]["original_species"], Species("O", -2))

        self.assertAlmostEqual(dopants["p_type"][0]["probability"], 0.023398867249112935)
        self.assertEqual(dopants["p_type"][0]["dopant_species"], Species("Co", 2))
        self.assertEqual(dopants["p_type"][0]["original_species"], Species("Sn", 4))

        # test oxidation sign matching
        dopants = get_dopants_from_substitution_probabilities(self.tin_dioxide, num_dopants=15, match_oxi_sign=False)
        self.assertEqual(dopants["n_type"][14]["dopant_species"], Species("Li", 1))
        self.assertEqual(dopants["n_type"][14]["original_species"], Species("O", -2))

        dopants = get_dopants_from_substitution_probabilities(self.tin_dioxide, num_dopants=15, match_oxi_sign=True)
        self.assertNotEqual(dopants["n_type"][14]["dopant_species"], Species("Li", 1))

    def test_dopants_from_shannon_radii(self):
        bonded_structure = CrystalNN().get_bonded_structure(self.tin_dioxide)

        dopants = get_dopants_from_shannon_radii(bonded_structure, num_dopants=5)

        self.assertTrue("n_type" in dopants)
        self.assertTrue("p_type" in dopants)

        self.assertTrue(len(dopants["n_type"]) <= 5)
        self.assertTrue(len(dopants["p_type"]) <= 5)

        self.assertAlmostEqual(dopants["n_type"][0]["radii_diff"], 0.04)
        self.assertEqual(dopants["n_type"][0]["dopant_species"], Species("U", 6))
        self.assertEqual(dopants["n_type"][0]["original_species"], Species("Sn", 4))

        self.assertAlmostEqual(dopants["p_type"][0]["radii_diff"], 0.0)
        self.assertEqual(dopants["p_type"][0]["dopant_species"], Species("Ni", 2))
        self.assertEqual(dopants["p_type"][0]["original_species"], Species("Sn", 4))
