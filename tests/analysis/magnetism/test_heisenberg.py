from __future__ import annotations

import pandas as pd
import pytest
from pytest import approx

from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper, HeisenbergModel
from pymatgen.core.structure import Structure
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/analysis/magnetic_orderings"


@pytest.fixture(scope="module")
def mn3al_data() -> pd.DataFrame:
    return pd.read_json(f"{TEST_DIR}/Mn3Al.json")


@pytest.fixture
def mn3al_inputs(mn3al_data) -> tuple[list[Structure], list[float]]:
    """Fresh Mn3Al orderings and total energies for each test."""
    structures = [Structure.from_dict(dct) for dct in mn3al_data["structure"]]
    energies = [epa * len(struct) for epa, struct in zip(mn3al_data["energy_per_atom"], structures, strict=True)]
    return structures, energies


@pytest.fixture
def mn3al_mapper(mn3al_inputs) -> HeisenbergMapper:
    """A fresh HeisenbergMapper per test, so no test depends on state left by another."""
    structures, energies = mn3al_inputs
    return HeisenbergMapper(structures, energies, cutoff=5.0, tol=0.02)


class TestHeisenbergMapper:
    def test_graphs(self, mn3al_mapper):
        assert len(mn3al_mapper.sgraphs) == 7

    def test_sites(self, mn3al_mapper):
        assert mn3al_mapper.unique_site_ids[0, 1] == 0

    def test_nn_interactions(self, mn3al_mapper):
        assert len(mn3al_mapper.nn_interactions) == 3
        assert mn3al_mapper.dists["nn"] == approx(2.51)

    def test_exchange_params(self, mn3al_mapper):
        ex_params = mn3al_mapper.get_exchange()
        assert ex_params["0-1-nn"] == approx(18.052116895702852, abs=1e-3)

    def test_mean_field(self, mn3al_mapper):
        # get_exchange() is deliberately not called first: get_mft_temperature
        # must compute the exchange parameters it needs on demand
        assert mn3al_mapper.ex_params is None
        j_avg = mn3al_mapper.estimate_exchange()
        assert j_avg == approx(52.54997149705518, abs=1e-3)

        mft_t = mn3al_mapper.get_mft_temperature(j_avg)
        assert round(mft_t) == 293

    def test_mean_field_only_average_exchange(self, mn3al_inputs):
        # with no cutoff only <J> can be fit, although Mn3Al has two magnetic
        # sublattices; this used to raise ValueError parsing "<J>" as a site index
        mapper = HeisenbergMapper(*mn3al_inputs)
        assert len(mapper.unique_site_ids) == 2
        j_avg = mapper.estimate_exchange()
        assert mapper.get_mft_temperature(j_avg) == approx(2 * j_avg / 3 / 0.0861733)
        assert set(mapper.ex_params) == {"<J>"}

    def test_get_igraph(self, mn3al_mapper):
        # also must not depend on get_exchange() having been called first
        igraph = mn3al_mapper.get_interaction_graph()
        assert len(igraph) == 6

    def test_heisenberg_model(self, mn3al_mapper):
        hmodel = mn3al_mapper.get_heisenberg_model()
        assert hmodel.formula == "Mn3Al"

    def test_as_from_dict_round_trip(self, mn3al_mapper):
        # HeisenbergModel must survive MSON round-trips. as_dict() serializes
        # ex_mat with jsanitize (a DataFrame becomes a dict), so from_dict()
        # must accept a dict and not only a string.
        # https://github.com/materialsproject/pymatgen/issues/4664
        model = mn3al_mapper.get_heisenberg_model()

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
