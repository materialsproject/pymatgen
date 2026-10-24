from __future__ import annotations

import pandas as pd
import pytest
from monty.serialization import loadfn
from pytest import approx

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.magnetism.heisenberg import HeisenbergMapper, HeisenbergModel, HeisenbergScreener
from pymatgen.core import Lattice
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

    def test_exchange_params_all_neighbors(self, mn3al_mapper):
        ex_params = mn3al_mapper.get_exchange()
        assert ex_params == approx(
            {
                "E0": -10.414980531535573,
                "0-1-nn": 18.052116895887593,
                "1-1-nnn": 7.113737574059087,
                "0-0-nnnn": 11.856870050536372,
                "1-1-nnnn": 0.29880538841098314,
            },
            abs=1e-6,
        )
        assert mn3al_mapper.wyckoff_ids == {0: "2b", 1: "4c"}
        assert mn3al_mapper.dists == approx({"nn": 2.51, "nnn": 2.9, "nnnn": 4.1})

    def test_get_low_energy_orderings(self, mn3al_mapper):
        # all Mn3Al orderings are ferrimagnetic, so the FM-like and AFM-like
        # states are picked by largest and smallest net moment
        fm_struct, afm_struct, fm_e, afm_e = mn3al_mapper.get_low_energy_orderings()
        assert fm_e == approx(-10.4357914284)
        assert afm_e == approx(-10.414978463333332)
        assert abs(sum(fm_struct.site_properties["magmom"])) == approx(2.37)
        assert abs(sum(afm_struct.site_properties["magmom"])) == approx(0.025)

    def test_small_moments_warning(self, mn3al_mapper, caplog):
        mn3al_mapper.estimate_exchange()
        assert "Local magnetic moments are small" in caplog.text

    def test_interaction_graph_to_file(self, mn3al_mapper, tmp_path):
        filename = tmp_path / "Mn3Al_igraph.json"
        igraph = mn3al_mapper.get_interaction_graph(filename=str(filename))
        reloaded = loadfn(filename)
        assert isinstance(reloaded, StructureGraph)
        assert reloaded.graph.number_of_edges() == igraph.graph.number_of_edges()
        weights = {round(data.get("weight", 0), 3) for *_, data in igraph.graph.edges(data=True)}
        # bonds between site pairs with no fitted J carry zero weight
        assert weights == {0, 18.052, 7.114, 11.857, 0.299}

    @pytest.mark.parametrize(
        ("i", "j", "dist", "expected"),
        [
            (0, 2, 2.51, 18.052116895887593),  # 2b-4c nn
            (2, 3, 2.90, 7.113737574059087),  # 4c-4c nnn
            (0, 1, 4.10, 11.856870050536372),  # 2b-2b nnnn
            (0, 2, 3.50, 0),  # not a mapped neighbor shell
        ],
    )
    def test_get_j_exc(self, mn3al_mapper, i, j, dist, expected):
        assert mn3al_mapper._get_j_exc(i, j, dist) == approx(expected, abs=1e-6)
        model = mn3al_mapper.get_heisenberg_model()
        assert model._get_j_exc(i, j, dist) == approx(expected, abs=1e-6)

    def test_get_unique_sites_ignores_wyckoff_property(self, mn3al_mapper):
        struct = mn3al_mapper.ordered_structures[0].copy()
        struct.add_site_property("wyckoff", ["x"] * len(struct))
        assert HeisenbergMapper._get_unique_sites(struct) == (
            mn3al_mapper.unique_site_ids,
            mn3al_mapper.wyckoff_ids,
        )

    def test_too_few_orderings(self, mn3al_inputs):
        structures, energies = mn3al_inputs
        with pytest.raises(SystemExit, match="We need at least 2 unique orderings"):
            HeisenbergMapper(structures[:1], energies[:1])


class TestHeisenbergScreener:
    def test_cleanup(self, mn3al_inputs):
        structures, energies = mn3al_inputs
        screener = HeisenbergScreener(structures, energies)
        # two of the nine orderings are degenerate in energy and are dropped
        assert len(structures) == 9
        assert len(screener.screened_structures) == len(screener.screened_energies) == 7
        # energies are per magnetic ion and sorted
        assert screener.screened_energies == sorted(screener.screened_energies)
        assert screener.screened_energies[0] == approx(-10.575664673333334)
        # nonmagnetic Al is stripped
        assert all(struct.composition.reduced_formula == "Mn" for struct in screener.screened_structures)

    def test_screen(self, mn3al_inputs):
        cleaned = HeisenbergScreener(*mn3al_inputs, screen=False)
        screened = HeisenbergScreener(*mn3al_inputs, screen=True)

        # the ground state and first excited state are kept in place...
        assert screened.screened_energies[:2] == approx(cleaned.screened_energies[:2])
        assert sorted(screened.screened_energies) == approx(cleaned.screened_energies)

        # ...and the rest prioritized by fewest moments below 1 muB
        n_below_1ub = [
            sum(abs(magmom) < 1 for magmom in struct.site_properties["magmom"])
            for struct in screened.screened_structures[2:]
        ]
        assert n_below_1ub == sorted(n_below_1ub)


@pytest.fixture
def bcc_mapper() -> HeisenbergMapper:
    """Single-sublattice bcc magnet with FM and AFM orderings, 0.1 eV apart."""
    lattice = Lattice.cubic(2.87)
    coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
    fm_struct = Structure(lattice, ["Fe", "Fe"], coords, site_properties={"magmom": [3.0, 3.0]})
    afm_struct = Structure(lattice, ["Fe", "Fe"], coords, site_properties={"magmom": [3.0, -3.0]})
    return HeisenbergMapper([afm_struct, fm_struct], [-0.9, -1.0])


class TestHeisenbergMapperSingleSublattice:
    # J = (E_AFM - E_FM) / |m|^2 = 0.05 eV / 9 per magnetic ion
    j_avg = 50 / 9  # meV
    t_mft = 2 * j_avg / 3 / 0.0861733  # K

    def test_mapping(self, bcc_mapper):
        assert bcc_mapper.unique_site_ids == {(0, 1): 0}
        assert bcc_mapper.dists == approx({"nn": 2.49, "nnn": 0, "nnnn": 0})
        # a single NN interaction cannot be fit, so ex_mat is left empty
        assert bcc_mapper.ex_mat.empty
        assert list(bcc_mapper.ex_mat.columns) == ["E", "E0", "0-0-nn"]

    def test_get_low_energy_orderings(self, bcc_mapper):
        fm_struct, afm_struct, fm_e, afm_e = bcc_mapper.get_low_energy_orderings()
        assert fm_struct.site_properties["magmom"] == [3.0, 3.0]
        assert sorted(afm_struct.site_properties["magmom"]) == [-3.0, 3.0]
        assert (fm_e, afm_e) == approx((-0.5, -0.45))

    def test_exchange_and_mean_field(self, bcc_mapper, caplog):
        assert bcc_mapper.estimate_exchange() == approx(self.j_avg)
        assert bcc_mapper.get_exchange() == approx({"<J>": self.j_avg})
        assert bcc_mapper.get_mft_temperature(self.j_avg) == approx(self.t_mft)
        assert self.t_mft == approx(42.9797, abs=1e-4)
        assert "too high" not in caplog.text

        assert bcc_mapper.get_mft_temperature(300) == approx(2 * 300 / 3 / 0.0861733)
        assert "This mean field estimate is too high" in caplog.text

    def test_interaction_graph(self, bcc_mapper, tmp_path, caplog):
        filename = tmp_path / "igraph"
        igraph = bcc_mapper.get_interaction_graph(filename=str(filename))
        assert "Only <J> is available" in caplog.text
        # 8 nearest neighbors per site, all carrying <J>
        assert igraph.graph.number_of_edges() == 8
        weights = [data["weight"] for *_, data in igraph.graph.edges(data=True)]
        assert weights == approx([self.j_avg] * 8)
        # .json is appended to the filename
        assert loadfn(tmp_path / "igraph.json").graph.number_of_edges() == 8

    def test_heisenberg_model(self, bcc_mapper):
        model = bcc_mapper.get_heisenberg_model()
        assert model.formula == "Fe"
        assert model.javg == approx(self.j_avg)
        assert model.ex_params == approx({"<J>": self.j_avg})
        assert model._get_j_exc(0, 1, 2.49) == approx(self.j_avg)
        assert model._get_j_exc(0, 1, 4.0) == 0

        model_rt = HeisenbergModel.from_dict(model.as_dict())
        assert model_rt.unique_site_ids == {(0, 1): 0}
        assert model_rt.nn_interactions == {"nn": {0: 0}, "nnn": {}, "nnnn": {}}
        assert model_rt.javg == approx(self.j_avg)
        assert model_rt.igraph.graph.number_of_edges() == 8

    def test_from_dict_legacy_fields(self, bcc_mapper):
        dct = bcc_mapper.get_heisenberg_model().as_dict()
        # older serializations: int site keys and an empty ex_mat string
        dct["unique_site_ids"] = {"0": 0}
        dct["ex_mat"] = ""
        model = HeisenbergModel.from_dict(dct)
        assert model.unique_site_ids == {(0,): 0}
        assert model.ex_mat.empty
        assert list(model.ex_mat.columns) == ["E", "E0"]
