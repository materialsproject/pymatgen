from __future__ import annotations

import unittest
import warnings

import pandas as pd

from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper
from pymatgen.core.structure import Structure
from pymatgen.util.testing import TEST_FILES_DIR

test_dir = f"{TEST_FILES_DIR}/magnetic_orderings"


class TestHeisenbergMapper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.df = pd.read_json(f"{test_dir}/mag_orderings_test_cases.json")

        # Good tests
        cls.Mn3Al = pd.read_json(f"{test_dir}/Mn3Al.json")

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
            assert len(sgraphs) == 7

    def test_sites(self):
        for hm in self.hms:
            unique_site_ids = hm.unique_site_ids
            assert unique_site_ids[(0, 1)] == 0

    def test_nn_interactions(self):
        for hm in self.hms:
            num_interacts = len(hm.nn_interactions)
            assert num_interacts == 3

            dists = hm.dists
            assert dists["nn"] == 2.51

    def test_exchange_params(self):
        for hm in self.hms:
            ex_params = hm.get_exchange()
            J_nn = round(18.052116895702852, 3)
            assert round(ex_params["0-1-nn"], 3) == J_nn

    def test_mean_field(self):
        for hm in self.hms:
            j_avg = hm.estimate_exchange()
            value = round(52.54997149705518, 3)
            assert round(j_avg, 3) == value

            mft_t = hm.get_mft_temperature(j_avg)
            value = round(292.90252668100584)
            assert round(mft_t) == value

    def test_get_igraph(self):
        for hm in self.hms:
            igraph = hm.get_interaction_graph()
            assert len(igraph) == 6

    def test_heisenberg_model(self):
        for hm in self.hms:
            hmodel = hm.get_heisenberg_model()
            assert hmodel.formula == "Mn3Al"
