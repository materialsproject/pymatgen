from __future__ import annotations

import pandas as pd
from pytest import approx

from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper, HeisenbergModel
from pymatgen.core.structure import Structure
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/analysis/magnetic_orderings"


class TestHeisenbergMapper:
    @classmethod
    def setup_class(cls):
        cls.df = pd.read_json(f"{TEST_DIR}/mag_orderings_test_cases.json")

        # Good tests
        cls.Mn3Al = pd.read_json(f"{TEST_DIR}/Mn3Al.json")

        cls.compounds = [cls.Mn3Al]

        cls.hms = []
        for c in cls.compounds:
            ordered_structures = list(c["structure"])
            ordered_structures = [Structure.from_dict(d) for d in ordered_structures]
            epa = list(c["energy_per_atom"])
            energies = [e * len(s) for (e, s) in zip(epa, ordered_structures, strict=True)]

            hm = HeisenbergMapper(ordered_structures, energies, cutoff=5.0, tol=0.02)
            cls.hms.append(hm)

    def test_graphs(self):
        for hm in self.hms:
            struct_graphs = hm.sgraphs
            assert len(struct_graphs) == 7

    def test_sites(self):
        for hm in self.hms:
            unique_site_ids = hm.unique_site_ids
            assert unique_site_ids[0, 1] == 0

    def test_nn_interactions(self):
        for hm in self.hms:
            n_interacts = len(hm.nn_interactions)
            assert n_interacts == 3

            dists = hm.dists
            assert dists["nn"] == approx(2.51)

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

    def test_as_from_dict_round_trip(self):
        # HeisenbergModel must survive MSON round-trips. as_dict() serializes
        # ex_mat with jsanitize (a DataFrame becomes a dict), so from_dict()
        # must accept a dict and not only a string.
        # https://github.com/materialsproject/pymatgen/issues/4664
        for hm in self.hms:
            model = hm.get_heisenberg_model()

            # First round-trip: ex_mat starts as a JSON string, becomes a DataFrame.
            model_rt = HeisenbergModel.from_dict(model.as_dict())
            assert isinstance(model_rt.ex_mat, pd.DataFrame)
            assert model_rt.formula == model.formula

            # Second round-trip: ex_mat is now a DataFrame, which jsanitize
            # serializes to a dict. This raised
            # "ValueError: malformed node or string" before the fix.
            model_rt2 = HeisenbergModel.from_dict(model_rt.as_dict())
            assert isinstance(model_rt2.ex_mat, pd.DataFrame)
            pd.testing.assert_frame_equal(model_rt.ex_mat, model_rt2.ex_mat)
