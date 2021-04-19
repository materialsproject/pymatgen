# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest
import warnings

import pandas as pd

from pymatgen.core.structure import Structure
from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "magnetic_orderings")


class HeisenbergMapperTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.df = pd.read_json(os.path.join(test_dir, "mag_orderings_test_cases.json"))

        # Good tests
        cls.Mn3Al = pd.read_json(os.path.join(test_dir, "Mn3Al.json"))

        cls.compounds = [cls.Mn3Al]

        cls.hms = []
        for c in cls.compounds:
            ordered_structures = list(c["structure"])
            ordered_structures = [Structure.from_dict(d) for d in ordered_structures]
            epa = list(c["energy_per_atom"])
            energies = [e * len(s) for (e, s) in zip(epa, ordered_structures)]

            hm = HeisenbergMapper(ordered_structures, energies, cutoff=5.0, tol=0.02)
            cls.hms.append(hm)

    def setUp(self):
        pass

    def tearDown(self):
        warnings.simplefilter("default")

    def test_graphs(self):
        for hm in self.hms:
            sgraphs = hm.sgraphs
            self.assertEqual(len(sgraphs), 7)

    def test_sites(self):
        for hm in self.hms:
            unique_site_ids = hm.unique_site_ids
            self.assertEqual(unique_site_ids[(0, 1)], 0)

    def test_nn_interactions(self):
        for hm in self.hms:
            num_interacts = len(hm.nn_interactions)
            self.assertEqual(num_interacts, 3)

            dists = hm.dists
            self.assertEqual(dists["nn"], 2.51)

    def test_exchange_params(self):
        for hm in self.hms:
            ex_params = hm.get_exchange()
            J_nn = round(18.052116895702852, 3)
            self.assertEqual(round(ex_params["0-1-nn"], 3), J_nn)

    def test_mean_field(self):
        for hm in self.hms:
            j_avg = hm.estimate_exchange()
            value = round(52.54997149705518, 3)
            self.assertEqual(round(j_avg, 3), value)

            mft_t = hm.get_mft_temperature(j_avg)
            value = round(292.90252668100584)
            self.assertEqual(round(mft_t), value)

    def test_get_igraph(self):
        for hm in self.hms:
            igraph = hm.get_interaction_graph()
            self.assertEqual(len(igraph), 6)

    def test_heisenberg_model(self):
        for hm in self.hms:
            hmodel = hm.get_heisenberg_model()
            self.assertEqual(hmodel.formula, "Mn3Al")


if __name__ == "__main__":
    unittest.main()
