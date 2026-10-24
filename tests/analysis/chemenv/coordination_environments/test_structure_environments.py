from __future__ import annotations

import json
import os

import matplotlib.pyplot as plt
import numpy as np
import orjson
import pytest
from matplotlib.figure import Figure
from monty.json import MontyEncoder
from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
    MultiWeightsChemenvStrategy,
    SimplestChemenvStrategy,
)
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import (
    ChemicalEnvironments,
    LightStructureEnvironments,
    StructureEnvironments,
)
from pymatgen.analysis.chemenv.utils.chemenv_errors import ChemenvError
from pymatgen.core import Species, Structure
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

__author__ = "waroquiers"

TEST_DIR = f"{TEST_FILES_DIR}/analysis/chemenv/structure_environments"


class TestStructureEnvironments(MatSciTest):
    def test_structure_environments(self):
        with open(f"{TEST_DIR}/se_mp-7000.json", "rb") as file:
            dct = orjson.loads(file.read())

        struct_envs = StructureEnvironments.from_dict(dct)
        isite = 6
        _csm_and_maps_fig, csm_and_maps_ax = struct_envs.get_csm_and_maps(isite=isite)
        assert_allclose(csm_and_maps_ax.lines[0].get_xydata().flatten(), [0, 0.53499332])
        assert_allclose(csm_and_maps_ax.lines[1].get_xydata().flatten(), [1, 0.47026441])
        assert_allclose(csm_and_maps_ax.lines[2].get_xydata().flatten(), [2, 0.00988778], atol=1e-8)

        _envs_fig, envs_ax = struct_envs.get_environments_figure(isite=isite)
        assert_allclose(
            envs_ax.patches[0].get_xy(),
            [
                [1, 1],
                [1, 0.99301365],
                [1.00179228, 0.99301365],
                [1.00179228, 1],
                [1, 1],
            ],
        )
        assert_allclose(
            envs_ax.patches[1].get_xy(),
            [
                [1, 0.99301365],
                [1, 0],
                [1.00179228, 0],
                [1.00179228, 0.99301365],
                [1, 0.99301365],
            ],
        )
        assert_allclose(
            envs_ax.patches[2].get_xy(),
            [
                [1.00179228, 1],
                [1.00179228, 0.99301365],
                [2.25, 0.99301365],
                [2.25, 1],
                [1.00179228, 1],
            ],
        )
        assert_allclose(
            envs_ax.patches[3].get_xy(),
            [
                [1.00179228, 0.99301365],
                [1.00179228, 0],
                [2.22376156, 0],
                [2.22376156, 0.0060837],
                [2.25, 0.0060837],
                [2.25, 0.99301365],
                [1.00179228, 0.99301365],
            ],
            atol=1e-8,
        )
        assert_allclose(
            envs_ax.patches[4].get_xy(),
            [
                [2.22376156, 0.0060837],
                [2.22376156, 0],
                [2.25, 0],
                [2.25, 0.0060837],
                [2.22376156, 0.0060837],
            ],
            atol=1e-8,
        )

        struct_envs.save_environments_figure(isite=isite, imagename="image.png")
        assert os.path.isfile("image.png")

        assert len(struct_envs.differences_wrt(struct_envs)) == 0

        assert struct_envs == struct_envs  # noqa: PLR0124

        ce = struct_envs.ce_list[isite][4][0]

        assert len(ce) == 4

        symbol, min_geometry = ce.minimum_geometry(symmetry_measure_type="csm_wocs_ctwocc")
        assert symbol == "T:4"
        assert min_geometry["symmetry_measure"] == approx(0.00988778424054)

        assert_allclose(
            min_geometry["other_symmetry_measures"]["rotation_matrix_wcs_csc"],
            [
                [-0.8433079817973094, -0.19705747216466898, 0.5000000005010193],
                [0.4868840909509757, 0.11377118475194581, 0.8660254034951744],
                [-0.22754236927612112, 0.9737681809261427, 1.3979531202869064e-13],
            ],
        )
        assert min_geometry["detailed_voronoi_index"] == {"index": 0, "cn": 4}
        assert min_geometry["other_symmetry_measures"]["scaling_factor_wocs_ctwocc"] == approx(1.627060)

        assert "csm1 (with central site) : 0.00988" in str(ce)
        assert "csm2 (without central site) : 0.00981" in str(ce)
        assert "csm1 (with central site) : 12.987" in str(ce)
        assert "csm2 (without central site) : 11.827" in str(ce)
        assert "csm1 (with central site) : 32.466" in str(ce)
        assert "csm2 (without central site) : 32.466" in str(ce)
        assert "csm1 (with central site) : 34.644" in str(ce)
        assert "csm2 (without central site) : 32.466" in str(ce)

        min_geoms = ce.minimum_geometries(symmetry_measure_type="csm_wocs_ctwocc", max_csm=12)
        assert len(min_geoms) == 2
        min_geoms = ce.minimum_geometries(symmetry_measure_type="csm_wocs_ctwcc", max_csm=12)
        assert len(min_geoms) == 1
        min_geoms = ce.minimum_geometries(n=3)
        assert len(min_geoms) == 3

        ce2 = struct_envs.ce_list[7][4][0]

        assert ce.is_close_to(ce2, rtol=0.01, atol=1e-4)
        assert not ce.is_close_to(ce2, rtol=0, atol=1e-8)

        assert ce != ce2
        assert ce != ce2

    def test_light_structure_environments(self):
        with open(f"{TEST_DIR}/se_mp-7000.json", "rb") as file:
            dct = orjson.loads(file.read())

        struct_envs = StructureEnvironments.from_dict(dct)

        strategy = SimplestChemenvStrategy()
        lse = LightStructureEnvironments.from_structure_environments(
            structure_environments=struct_envs, strategy=strategy, valences="undefined"
        )
        isite = 6
        nb_set = lse.neighbors_sets[isite][0]
        neighb_coords = [
            np.array([0.2443798, 1.80409653, -1.13218359]),
            np.array([1.44020353, 1.11368738, 1.13218359]),
            np.array([2.75513098, 2.54465207, -0.70467298]),
            np.array([0.82616785, 3.65833945, 0.70467298]),
        ]
        neighb_indices = [0, 3, 5, 1]
        neighb_images = [[0, 0, -1], [0, 0, 0], [0, 0, -1], [0, 0, 0]]

        assert_allclose(neighb_coords, nb_set.neighb_coords)
        assert_allclose(neighb_coords, [s.coords for s in nb_set.neighb_sites])
        nb_sai = nb_set.neighb_sites_and_indices
        assert_allclose(neighb_coords, [sai["site"].coords for sai in nb_sai])
        assert_allclose(neighb_indices, [sai["index"] for sai in nb_sai])
        nb_iai = nb_set.neighb_indices_and_images
        assert_allclose(neighb_indices, [iai["index"] for iai in nb_iai])
        assert_array_equal(neighb_images, [iai["image_cell"] for iai in nb_iai])

        assert len(nb_set) == 4
        assert hash(nb_set) == 4

        assert nb_set == nb_set  # noqa: PLR0124

        assert (
            str(nb_set) == "Neighbors Set for site #6 :\n"
            " - Coordination number : 4\n"
            " - Neighbors sites indices : 0, 1, 2, 3\n"
        )

        stats = lse.get_statistics()

        neighbors = lse.strategy.get_site_neighbors(site=lse.structure[isite])
        assert_allclose(neighbors[0].coords, [0.2443798, 1.80409653, -1.13218359])
        assert_allclose(neighbors[1].coords, [1.44020353, 1.11368738, 1.13218359])
        assert_allclose(neighbors[2].coords, [2.75513098, 2.54465207, -0.70467298])
        assert_allclose(neighbors[3].coords, [0.82616785, 3.65833945, 0.70467298])

        site_idx, d_equiv_site, d_this_site, sym = lse.strategy.equivalent_site_index_and_transform(neighbors[0])
        assert site_idx == 0
        assert_allclose(d_equiv_site, [0, 0, 0])
        assert_allclose(d_this_site, [0, 0, -1])
        assert_allclose(sym.affine_matrix, np.eye(4))

        site_idx, d_equiv_site, d_this_site, sym = lse.strategy.equivalent_site_index_and_transform(neighbors[1])
        assert site_idx == 3
        assert_allclose(d_equiv_site, [0, 0, 0])
        assert_allclose(d_this_site, [0, 0, 0], atol=1e-9)
        # assert_allclose(sym.affine_matrix, np.eye(4))

        assert stats["atom_coordination_environments_present"] == {"Si": {"T:4": 3}}
        assert stats["coordination_environments_atom_present"] == {"T:4": {"Si": 3}}
        assert stats["fraction_atom_coordination_environments_present"] == {"Si": {"T:4": 1}}

        site_info_ce = lse.get_site_info_for_specie_ce(specie=Species("Si4+"), ce_symbol="T:4")
        assert_allclose(site_info_ce["fractions"], [1, 1, 1])
        assert_allclose(
            site_info_ce["csms"],
            [0.009887784240541068, 0.009887786546730826, 0.009887787384385317],
        )
        assert site_info_ce["isites"] == [6, 7, 8]

        site_info_allces = lse.get_site_info_for_specie_allces(specie=Species("Si4+"))

        assert site_info_allces["T:4"] == site_info_ce

        assert not lse.contains_only_one_anion("I-")
        assert not lse.contains_only_one_anion_atom("I")
        assert lse.site_contains_environment(isite=isite, ce_symbol="T:4")
        assert not lse.site_contains_environment(isite=isite, ce_symbol="S:4")
        assert not lse.structure_contains_atom_environment(atom_symbol="Si", ce_symbol="S:4")
        assert lse.structure_contains_atom_environment(atom_symbol="Si", ce_symbol="T:4")
        assert not lse.structure_contains_atom_environment(atom_symbol="O", ce_symbol="T:4")
        assert lse.uniquely_determines_coordination_environments
        assert lse == lse  # noqa: PLR0124

        envs = lse.strategy.get_site_coordination_environments(lse.structure[6])
        assert len(envs) == 1
        assert envs[0][0] == "T:4"

        multi_strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()

        lse_multi = LightStructureEnvironments.from_structure_environments(
            strategy=multi_strategy,
            structure_environments=struct_envs,
            valences="undefined",
        )
        assert lse_multi.coordination_environments[isite][0]["csm"] == approx(0.009887784240541068)
        assert lse_multi.coordination_environments[isite][0]["ce_fraction"] == approx(1)
        assert lse_multi.coordination_environments[isite][0]["ce_symbol"] == "T:4"

    def test_from_structure_environments(self):
        # https://github.com/materialsproject/pymatgen/issues/2756
        struct = Structure.from_file(f"{TEST_DIR}/mp-554015.json.gz")
        strategy = SimplestChemenvStrategy(distance_cutoff=1.4, angle_cutoff=0.3)
        local_geom_finder = LocalGeometryFinder()
        local_geom_finder.setup_structure(structure=struct)
        envs = local_geom_finder.compute_structure_environments()

        lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=envs)
        assert len(lse.coordination_environments) == len(struct)


def _load_se(name: str) -> StructureEnvironments:
    with open(f"{TEST_DIR}/se_{name}.json", "rb") as file:
        return StructureEnvironments.from_dict(orjson.loads(file.read()))


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def se_7000() -> StructureEnvironments:
    """SiO2 (mp-7000): only the 3 Si sites (6, 7, 8) have computed environments."""
    return _load_se("mp-7000")


@pytest.fixture
def se_5020() -> StructureEnvironments:
    """BaTiO3 (mp-5020): Ba (0) and Ti (1) computed, info contains sites_info."""
    return _load_se("mp-5020")


@pytest.fixture
def se_743972() -> StructureEnvironments:
    """K2NaNb2Fe7Si8H4O31 (mp-743972), 55 sites, cations 0-19 computed."""
    return _load_se("mp-743972")


@pytest.fixture
def lse_5020(se_5020) -> LightStructureEnvironments:
    return LightStructureEnvironments.from_structure_environments(
        strategy=SimplestChemenvStrategy(),
        structure_environments=se_5020,
        valences=[2, 4, -2, -2, -2],
    )


@pytest.fixture
def lse_743972_multi(se_743972) -> LightStructureEnvironments:
    return LightStructureEnvironments.from_structure_environments(
        strategy=MultiWeightsChemenvStrategy.stats_article_weights_parameters(),
        structure_environments=se_743972,
    )


CSM_TYPES = ("csm_wcs_ctwcc", "csm_wcs_ctwocc", "csm_wcs_csc", "csm_wocs_ctwcc", "csm_wocs_ctwocc", "csm_wocs_csc")


def _make_ce(mp_symbol: str = "T:4", csm: float = 0.5, **kwargs) -> ChemicalEnvironments:
    ce = ChemicalEnvironments()
    ce.add_coord_geom(
        mp_symbol,
        csm,
        algo=kwargs.get("algo", "EXPLICIT_PERMUTATIONS"),
        permutation=kwargs.get("permutation", [0, 1, 2, 3]),
        detailed_voronoi_index=kwargs.get("detailed_voronoi_index", {"index": 0, "cn": 4}),
        other_symmetry_measures=dict.fromkeys(CSM_TYPES, csm),
    )
    return ce


class TestSENeighborsSet:
    def test_init_validation(self, se_7000):
        nb_cls = StructureEnvironments.NeighborsSet
        kwargs = {"structure": se_7000.structure, "isite": 6, "detailed_voronoi": se_7000.voronoi}
        with pytest.raises(TypeError, match="must be list/tuple of ints"):
            nb_cls(site_voronoi_indices={1, 4}, **kwargs)
        with pytest.raises(ValueError, match="contains 2 duplicates"):
            nb_cls(site_voronoi_indices=[1, 1, 4, 4], **kwargs)

        nb_set = nb_cls(site_voronoi_indices=(6, 1, 5, 4), **kwargs)
        assert nb_set.site_voronoi_indices == (1, 4, 5, 6)
        assert nb_set.sources == [{"origin": "UNKNOWN"}]
        assert nb_set.source == {"origin": "UNKNOWN"}
        assert nb_set == se_7000.neighbors_sets[6][4][0]
        assert nb_set != "not a neighbors set"

        nb_set = nb_cls(site_voronoi_indices=[1, 4], sources={"origin": "from_hints"}, **kwargs)
        assert nb_set.sources == [{"origin": "from_hints"}]
        nb_set.add_source({"origin": "from_hints"})
        assert len(nb_set.sources) == 1
        nb_set.add_source({"origin": "other"})
        assert len(nb_set.sources) == 2
        with pytest.raises(RuntimeError, match="Number of sources different from 1"):
            _ = nb_set.source
        with pytest.raises(NotImplementedError, match="Nothing implemented for other sources"):
            nb_set.voronoi_grid_surface_points(other_origins="RAISE")
        assert nb_set.voronoi_grid_surface_points() is None

    def test_geometry_properties(self, se_7000):
        nb_set = se_7000.neighbors_sets[6][4][0]
        assert nb_set.site_voronoi_indices == (1, 4, 5, 6)
        assert nb_set.get_neighb_voronoi_indices([3, 2, 1, 0]) == [6, 5, 4, 1]
        coords = nb_set.coords
        assert len(coords) == 5
        assert_allclose(coords[0], se_7000.structure[6].coords)
        assert_allclose(nb_set.distances, [1.62843998, 1.62843998, 1.62552659, 1.62552659])
        assert_allclose(nb_set.angles, [3.11738988, 3.11738988, 3.09561071, 3.09561071])

        info = nb_set.info
        assert info["distances_mean"] == approx(1.626983283472464)
        assert info["normalized_angles_sum"] == approx(3.986027306187217)
        assert info["normalized_distances_min"] == approx(1.0)
        assert info["normalized_angles_max"] == approx(1.0)

        assert nb_set.distance_plateau() == approx(1.2219692770485542)
        assert nb_set.angle_plateau() == approx(0.9765834191971131)
        assert se_7000.neighbors_sets[6][6][0].distance_plateau() == approx(0.10464832558746995)
        assert se_7000.neighbors_sets[6][6][0].angle_plateau() == approx(0.006083703662383561)

        assert (
            str(nb_set) == "Neighbors Set for site #6 :\n - Coordination number : 4\n - Voronoi indices : 1, 4, 5, 6\n"
        )

    def test_voronoi_grid_surface_points(self, se_7000):
        nb_set = se_7000.neighbors_sets[6][4][0]
        assert_allclose(
            nb_set.voronoi_grid_surface_points(),
            [
                (1.00179228, 0.99301365),
                (1.00179228, 0.0),
                (2.22376156, 0.0),
                (2.22376156, 0.0060837),
                (2.25, 0.0060837),
                (2.25, 0.99301365),
            ],
            atol=1e-8,
        )
        # No source for an additional condition that was not computed
        assert nb_set.voronoi_grid_surface_points(additional_condition=99) is None

    def test_dict_round_trip(self, se_7000):
        nb_set = se_7000.neighbors_sets[6][4][0]
        dct = orjson.loads(orjson.dumps(nb_set.as_dict()))
        assert dct["isite"] == 6
        assert dct["site_voronoi_indices"] == [1, 4, 5, 6]
        new_nb_set = StructureEnvironments.NeighborsSet.from_dict(
            dct, structure=se_7000.structure, detailed_voronoi=se_7000.voronoi
        )
        assert new_nb_set == nb_set
        assert new_nb_set.sources == nb_set.sources
        assert hash(new_nb_set) == 4


class TestStructureEnvironmentsMethods:
    def test_init_neighbors_sets(self, se_7000):
        # Oxygen sites have no Voronoi information: nothing is done
        assert se_7000.init_neighbors_sets(isite=0) is None
        assert se_7000.neighbors_sets[0] is None

        for ac in (StructureEnvironments.AC.ONLY_ACB, StructureEnvironments.AC.ONLY_ACB_AND_NO_E2SEB):
            with pytest.raises(ChemenvError, match="Valences are not given"):
                se_7000.init_neighbors_sets(isite=6, additional_conditions=[ac])

        # Conditions that do not need valences work without them
        se_7000.neighbors_sets[6] = None
        se_7000.ce_list[6] = None
        se_7000.init_neighbors_sets(isite=6, additional_conditions=[StructureEnvironments.AC.NONE])
        assert sorted(se_7000.neighbors_sets[6]) == [0, 2, 4, 6]
        nb_set_t4 = se_7000.neighbors_sets[6][4][0]
        assert nb_set_t4.site_voronoi_indices == (1, 4, 5, 6)
        assert {src["ac"] for src in nb_set_t4.sources} == {StructureEnvironments.AC.NONE}
        assert se_7000.ce_list[6][4] == [None]

    def test_update_coordination_environments(self, se_7000):
        nb_set = se_7000.neighbors_sets[6][4][0]
        ce = _make_ce()
        se_7000.update_coordination_environments(isite=6, cn=4, nb_set=nb_set, ce=ce)
        assert se_7000.ce_list[6][4][0] is ce
        assert se_7000.get_coordination_environments(isite=6, cn=4, nb_set=nb_set) is ce

        # ce_list for this site not yet initialized: it is created and the CE appended
        se_7000.ce_list[6] = None
        se_7000.update_coordination_environments(isite=6, cn=4, nb_set=nb_set, ce=ce)
        assert se_7000.ce_list[6] == {4: [ce]}

        other_nb_set = StructureEnvironments.NeighborsSet(
            structure=se_7000.structure,
            isite=6,
            detailed_voronoi=se_7000.voronoi,
            site_voronoi_indices=[0, 1, 2, 3],
        )
        with pytest.raises(ValueError, match="Neighbors set not found"):
            se_7000.update_coordination_environments(isite=6, cn=4, nb_set=other_nb_set, ce=ce)

        se_7000.ce_list[6] = {4: []}
        se_7000.neighbors_sets[6][4].insert(0, other_nb_set)
        with pytest.raises(ValueError, match="Neighbors set not yet in ce_list"):
            se_7000.update_coordination_environments(isite=6, cn=4, nb_set=nb_set, ce=ce)

    def test_add_neighbors_set(self, se_7000):
        new_source = {"origin": "from_hints"}
        nb_set = StructureEnvironments.NeighborsSet(
            structure=se_7000.structure,
            isite=6,
            detailed_voronoi=se_7000.voronoi,
            site_voronoi_indices=[1, 4, 5, 6],
            sources=new_source,
        )
        n_sources = len(se_7000.neighbors_sets[6][4][0].sources)
        se_7000.add_neighbors_set(isite=6, nb_set=nb_set)
        # Same neighbors set: only the source is added
        assert len(se_7000.neighbors_sets[6][4]) == 1
        assert len(se_7000.neighbors_sets[6][4][0].sources) == n_sources + 1
        assert se_7000.neighbors_sets[6][4][0].sources[-1] == new_source

        nb_set_3 = StructureEnvironments.NeighborsSet(
            structure=se_7000.structure,
            isite=6,
            detailed_voronoi=se_7000.voronoi,
            site_voronoi_indices=[1, 4, 5],
            sources=new_source,
        )
        se_7000.add_neighbors_set(isite=6, nb_set=nb_set_3)
        assert se_7000.neighbors_sets[6][3] == [nb_set_3]
        assert se_7000.ce_list[6][3] == [None]

    def test_get_coordination_environments_missing(self, se_7000):
        nb_set = se_7000.neighbors_sets[6][4][0]
        assert se_7000.get_coordination_environments(isite=0, cn=4, nb_set=nb_set) is None
        assert se_7000.get_coordination_environments(isite=6, cn=5, nb_set=nb_set) is None
        nb_set_6 = se_7000.neighbors_sets[6][6][0]
        assert se_7000.get_coordination_environments(isite=6, cn=4, nb_set=nb_set_6) is None
        ce = se_7000.get_coordination_environments(isite=6, cn=4, nb_set=nb_set)
        assert ce.minimum_geometry()[0] == "T:4"

    def test_get_csm(self, se_7000):
        csm = se_7000.get_csm(6, "T:4")
        assert csm["symmetry_measure"] == approx(0.009887784240541068)
        assert csm["algo"] == "EXPLICIT_PERMUTATIONS"
        assert se_7000.get_csms(6, "C:12") == []
        # Two neighbors sets with 2 neighbors -> two csms for a linear environment
        assert len(se_7000.get_csms(6, "L:2")) == 2
        with pytest.raises(ChemenvError, match="Number of csms for site #6 with mp_symbol 'L:2' = 2"):
            se_7000.get_csm(6, "L:2")

    def test_update_site_info(self, se_7000, se_5020):
        assert "sites_info" not in se_7000.info
        se_7000.update_site_info(isite=6, info_dict={"time": 1.5})
        assert len(se_7000.info["sites_info"]) == 9
        assert se_7000.info["sites_info"][6] == {"time": 1.5}
        assert se_7000.info["sites_info"][0] == {}

        se_5020.update_site_info(isite=0, info_dict={"extra": 1})
        assert se_5020.info["sites_info"][0]["extra"] == 1
        assert se_5020.info["sites_info"][0]["time"] == approx(0.44231447700000004)

    def test_plot_csm_and_maps(self, se_7000, monkeypatch):
        shown = []
        monkeypatch.setattr(plt, "show", lambda: shown.append(True))
        assert se_7000.plot_csm_and_maps(isite=6) is None
        assert shown == [True]

        fig, ax = se_7000.get_csm_and_maps(
            isite=6, max_csm=1.0, figsize=(4, 3), symmetry_measure_type="csm_wocs_ctwocc"
        )
        assert_allclose(fig.get_size_inches(), [4, 3])
        # Only the two linear CN 2 sets and the T:4 set have a csm below 1.0 (CN 6 is dropped)
        assert [txt.get_text() for txt in ax.texts] == ["A:2", "A:2", "T:4"]
        assert_allclose(ax.lines[2].get_xydata().flatten(), [2, 0.00981014], atol=1e-8)
        assert [lbl.get_text() for lbl in ax.get_xticklabels()] == ["(2, 0)", "(2, 1)", "(4, 0)"]

    @pytest.mark.parametrize(
        ("max_wd", "expected_top_label"),
        [(1.1, 1.25), (1.3, 1.4), (1.4, 1.5), (1.5, 1.6), (1.7, 1.8), (1.9, 2.0), (2.2, 2.5), (3.0, 4.0)],
    )
    def test_csm_and_maps_distance_axis(self, se_7000, monkeypatch, max_wd, expected_top_label):
        monkeypatch.setattr(
            StructureEnvironments.NeighborsSet,
            "normalized_distances",
            property(lambda self: [1.0] * (len(self) - 1) + [max_wd]),
        )
        fig, ax = se_7000.get_csm_and_maps(isite=6)
        ax_distang = fig.axes[1]
        assert ax_distang is not ax
        labels = [float(lbl.get_text()) for lbl in ax_distang.get_yticklabels()]
        assert labels[0] == approx(0.0)
        assert labels[-1] == approx(expected_top_label)

    @pytest.mark.parametrize(
        ("plot_type", "xlabel", "ylabel", "first_patch_x"),
        [
            (
                {"distance_parameter": ("initial_normalized", None), "angle_parameter": ("initial_normalized", None)},
                "Distance parameter : $\\alpha$",
                "Angle parameter : $\\gamma$",
                1.0,
            ),
            (
                {
                    "distance_parameter": ("one_minus_inverse_alpha_power_n", None),
                    "angle_parameter": ("one_minus_gamma", None),
                },
                "Distance parameter : $1.0-\\frac{1.0}{\\alpha^{3}}$",
                "Angle parameter : $1.0-\\gamma$",
                0.0,
            ),
            (
                {
                    "distance_parameter": ("one_minus_inverse_alpha_power_n", {"exponent": 2}),
                    "angle_parameter": ("initial_normalized_inverted", None),
                },
                "Distance parameter : $1.0-\\frac{1.0}{\\alpha^{2}}$",
                "Angle parameter : $\\gamma$",
                0.0,
            ),
        ],
    )
    def test_environments_figure_plot_types(self, se_7000, plot_type, xlabel, ylabel, first_patch_x):
        _fig, ax = se_7000.get_environments_figure(isite=6, plot_type=plot_type, title="SiO2", figsize=(5, 4))
        assert ax.get_title() == "SiO2"
        assert ax.get_xlabel() == xlabel
        assert ax.get_ylabel() == ylabel
        assert len(ax.patches) == 5
        assert ax.patches[0].get_xy()[0][0] == approx(first_patch_x)

    @pytest.mark.parametrize("param", ["distance_parameter", "angle_parameter"])
    def test_environments_figure_wrong_plot_type(self, se_7000, param):
        plot_type = {
            "distance_parameter": ("initial_normalized", None),
            "angle_parameter": ("initial_normalized", None),
        }
        plot_type[param] = ("unknown", None)
        with pytest.raises(ValueError, match=f'Wrong value for {param.replace("_", " ")} plot type "unknown"'):
            se_7000.get_environments_figure(isite=6, plot_type=plot_type)

    def test_environments_figure_strategy_and_missing_ce(self, se_7000):
        # Missing CE (None) and CE without any geometry are both drawn in white with the CN as label
        se_7000.ce_list[6][2][0] = None
        se_7000.ce_list[6][2][1] = ChemicalEnvironments()
        strategy = SimplestChemenvStrategy(distance_cutoff=1.4, angle_cutoff=0.3)
        _fig, ax = se_7000.get_environments_figure(isite=6, strategy=strategy)
        assert len(ax.patches) == 5
        white_patches = [patch for patch in ax.patches if tuple(patch.get_facecolor()) == (1.0, 1.0, 1.0, 1.0)]
        assert len(white_patches) == 3  # CN 0 (no CE computed) and the two CN 2 sets
        # The strategy adds its (distance, angle) cutoff point to the plot
        assert_allclose(ax.lines[0].get_xydata(), [[1.4, 0.3]])

        # A strategy without visualization support is silently ignored
        _fig, ax = se_7000.get_environments_figure(isite=6, strategy=object())
        assert len(ax.lines) == 0

        # Neighbors sets that do not come from the Voronoi grid (e.g. from hints) are not drawn
        se_7000.add_neighbors_set(
            isite=6,
            nb_set=StructureEnvironments.NeighborsSet(
                structure=se_7000.structure,
                isite=6,
                detailed_voronoi=se_7000.voronoi,
                site_voronoi_indices=[1, 4, 5],
                sources={"origin": "from_hints"},
            ),
        )
        _fig, ax = se_7000.get_environments_figure(isite=6)
        assert len(ax.patches) == 5

    def test_environments_figure_uneven_points(self, se_7000, monkeypatch):
        monkeypatch.setattr(
            StructureEnvironments.NeighborsSet,
            "voronoi_grid_surface_points",
            lambda self, *args, **kwargs: [(1.0, 0.0), (1.5, 0.0), (1.5, 1.0)],
        )
        with pytest.raises(RuntimeError, match="Uneven number of surface points"):
            se_7000.get_environments_figure(isite=6)

    def test_plot_and_save_environments(self, se_7000, tmp_path, monkeypatch):
        shown = []
        monkeypatch.setattr(Figure, "show", lambda self, *args, **kwargs: shown.append(self))
        se_7000.plot_environments(isite=6, title="Si env")
        assert len(shown) == 1
        assert shown[0].axes[0].get_title() == "Si env"

        img_path = tmp_path / "env.png"
        se_7000.save_environments_figure(isite=6, imagename=img_path, figsize=(4, 3))
        assert img_path.is_file()
        assert img_path.stat().st_size > 0

    def test_eq(self, se_7000, se_5020):
        assert se_7000 == _load_se("mp-7000")
        assert se_7000 != se_5020
        assert se_7000 != "not a StructureEnvironments"

        for attr, value in [
            ("sites_map", [0] * 9),
            ("equivalent_sites", []),
            ("info", {}),
            ("valences", [0]),
            ("ce_list", [None]),
        ]:
            modified = _load_se("mp-7000")
            setattr(modified, attr, value)
            assert se_7000 != modified, attr

        modified = _load_se("mp-7000")
        modified.structure = modified.structure.copy()
        modified.structure.perturb(0.01, seed=0)
        assert se_7000 != modified

        modified = _load_se("mp-7000")
        modified.voronoi.normalized_angle_tolerance += 0.1
        assert se_7000 != modified

        modified = _load_se("mp-7000")
        del modified.neighbors_sets[6][6]
        assert se_7000 != modified

        modified = _load_se("mp-7000")
        modified.ce_list[6][4][0] = _make_ce()
        assert se_7000 != modified

    def test_differences_wrt_structure(self, se_7000, se_5020):
        diffs = se_7000.differences_wrt(se_5020)
        assert [diff["difference"] for diff in diffs] == ["structure", "PREVIOUS DIFFERENCE IS DISMISSIVE"]

    @pytest.mark.parametrize(("delta", "comparison"), [(1e-12, "__eq__"), (0.1, "is_close_to")])
    def test_differences_wrt_voronoi(self, se_7000, delta, comparison):
        other = _load_se("mp-7000")
        other.voronoi.normalized_angle_tolerance += delta
        other.valences = [0] * 9
        other.info = {}
        diffs = se_7000.differences_wrt(other)
        assert [(diff["difference"], diff["comparison"]) for diff in diffs] == [
            ("valences", "__eq__"),
            ("info", "__eq__"),
            ("voronoi", comparison),
            ("PREVIOUS DIFFERENCE IS DISMISSIVE", "differences_wrt"),
        ]

    def test_differences_wrt_neighbors_and_ces(self, se_7000):
        other = _load_se("mp-7000")
        # site 6: self has neighbors, other does not
        other.neighbors_sets[6] = None
        # site 7: different coordination numbers and different nb sets for CN 4
        del other.neighbors_sets[7][6]
        other.neighbors_sets[7][4][0] = StructureEnvironments.NeighborsSet(
            structure=other.structure,
            isite=7,
            detailed_voronoi=other.voronoi,
            site_voronoi_indices=[0, 1, 2, 3],
        )
        # site 8: CE close to the original one (CN 4) and CE far from the original one (CN 2)
        ce_close = ChemicalEnvironments.from_dict(other.ce_list[8][4][0].as_dict())
        ce_close.coord_geoms["T:4"]["symmetry_measure"] += 1e-12
        other.ce_list[8][4][0] = ce_close
        other.ce_list[8][2][0] = _make_ce("L:2", csm=5.0, permutation=[0, 1])

        diffs = {(diff["difference"], diff["comparison"]) for diff in se_7000.differences_wrt(other)}
        assert diffs == {
            ("neighbors_sets[site_idx=6]", "has_neighbors"),
            ("neighbors_sets[site_idx=7]", "coordination_numbers"),
            ("neighbors_sets[site_idx=7][cn=4]", "neighbors_sets"),
            ("ce_list[site_idx=8][cn=4][inb_set=0]", "__eq__"),
            ("ce_list[site_idx=8][cn=2][inb_set=0]", "is_close_to"),
        }

        # Reverse direction: self has no neighbors for site 6, other has
        diffs = other.differences_wrt(se_7000)
        assert diffs[0]["difference"] == "neighbors_sets[site_idx=6]"
        assert diffs[0]["self"] == "None"
        assert diffs[0]["other"] == {0, 2, 4, 6}

    def test_dict_json_round_trip(self, se_5020):
        dct = se_5020.as_dict()
        assert dct["@class"] == "StructureEnvironments"
        assert set(dct["neighbors_sets"][0]) == {"3", "9", "12"}
        assert dct["neighbors_sets"][2] is None
        assert set(dct["info"]["sites_info"][0]["nb_sets_info"]) == {"3", "9", "12"}
        assert dct["info"]["sites_info"][2] == {}

        new_se = StructureEnvironments.from_dict(orjson.loads(orjson.dumps(dct)))
        assert new_se == se_5020
        assert new_se.differences_wrt(se_5020) == []
        assert new_se.info["sites_info"][1]["nb_sets_info"][6][0]["time"] == approx(0.03131982499999997)
        assert new_se.get_csm(0, "C:12")["symmetry_measure"] == approx(se_5020.get_csm(0, "C:12")["symmetry_measure"])

        # "None" strings and missing sites_info are accepted when deserializing
        dct = orjson.loads(orjson.dumps(dct))
        dct["ce_list"][2] = "None"
        dct["ce_list"][0]["3"] = ["None"]
        del dct["info"]["sites_info"]
        new_se = StructureEnvironments.from_dict(dct)
        assert new_se.ce_list[2] is None
        assert new_se.ce_list[0][3] == [None]
        assert "sites_info" not in new_se.info


class TestLightStructureEnvironmentsMethods:
    def test_neighbors_set(self, lse_5020):
        nb_set = lse_5020.neighbors_sets[1][0]
        assert len(nb_set) == 6
        assert nb_set != "not a neighbors set"
        with pytest.raises(ValueError, match="Set of neighbors contains duplicates"):
            LightStructureEnvironments.NeighborsSet(
                structure=lse_5020.structure,
                isite=1,
                all_nbs_sites=lse_5020._all_nbs_sites,
                all_nbs_sites_indices=[0, 0],
            )

    def test_valences_origin(self, se_5020, lse_5020):
        assert lse_5020.valences == [2, 4, -2, -2, -2]
        assert lse_5020.valences_origin == "user-specified"
        lse = LightStructureEnvironments.from_structure_environments(
            strategy=SimplestChemenvStrategy(),
            structure_environments=se_5020,
            valences=[2, 4, -2, -2, -2],
            valences_origin="bond-valence-analysis",
        )
        assert lse.valences_origin == "bond-valence-analysis"

    def test_statistics_with_valences(self, lse_5020):
        stats = lse_5020.get_statistics()
        assert stats["anion_list"] == {"O2-": 3}
        assert stats["anion_atom_list"] == {"O": 3}
        assert stats["cation_list"] == {"Ba2+": 1, "Ti4+": 1}
        assert stats["cation_atom_list"] == {"Ba": 1, "Ti": 1}
        assert stats["neutral_list"] == {}
        assert stats["atom_coordination_environments_present"] == {"Ba": {"C:12": 1}, "Ti": {"O:6": 1}}
        assert stats["ion_coordination_environments_present"] == {"Ba": {2: {"C:12": 1}}, "Ti": {4: {"O:6": 1}}}
        assert stats["fraction_ion_coordination_environments_present"] == {
            "Ba": {2: {"C:12": 1}},
            "Ti": {4: {"O:6": 1}},
        }
        assert stats["coordination_environments_ion_present"] == {"C:12": {"Ba": {2: 1}}, "O:6": {"Ti": {4: 1}}}

        all_stats = lse_5020.get_statistics(statistics_fields="ALL", bson_compatible=True)
        assert all_stats["valences_origin"] == "user-specified"
        assert all_stats["anion_number"] == 1
        assert all_stats["cation_number"] == 2
        assert all_stats["cation_atom_number"] == 2
        assert all_stats["count_ion_present"] == {"Ba": {"2": 1}, "Ti": {"4": 1}}
        assert all_stats["fraction_coordination_environments_ion_present"] == {
            "C:12": {"Ba": {"2": 1}},
            "O:6": {"Ti": {"4": 1}},
        }

        assert lse_5020.contains_only_one_anion("O2-")
        assert not lse_5020.contains_only_one_anion("F-")
        assert lse_5020.contains_only_one_anion_atom("O")
        assert not lse_5020.contains_only_one_anion_atom("F")

    def test_statistics_neutral(self, se_7000):
        lse = LightStructureEnvironments.from_structure_environments(
            strategy=SimplestChemenvStrategy(),
            structure_environments=se_7000,
            valences=[0] * 9,
        )
        # Statistics are set up lazily by any of the queries
        assert lse.statistics_dict is None
        assert not lse.contains_only_one_anion("O2-")
        assert lse.statistics_dict is not None
        lse.statistics_dict = None
        assert not lse.contains_only_one_anion_atom("O")
        stats = lse.get_statistics(statistics_fields=["neutral_list", "neutral_atom_list"])
        assert stats == {"neutral_list": {"O0+": 6, "Si0+": 3}, "neutral_atom_list": {"O": 6, "Si": 3}}

    def test_site_info_for_specie(self, lse_5020, lse_743972_multi):
        assert lse_5020.get_site_info_for_specie_ce(specie=Species("Ti", 4), ce_symbol="O:6") == {
            "isites": [1],
            "fractions": [1],
            "csms": [approx(0.23337176618721564)],
        }
        assert lse_5020.get_site_info_for_specie_ce(specie=Species("Ti", 3), ce_symbol="O:6")["isites"] == []
        allces = lse_5020.get_site_info_for_specie_allces(specie=Species("Ba", 2))
        assert list(allces) == ["C:12"]
        assert allces["C:12"]["isites"] == [0]
        # A different element with the same oxidation state must not match
        assert lse_5020.get_site_info_for_specie_allces(specie=Species("Sr", 2)) == {}

        # Structure with more sites than len("undefined") and several CEs per site
        si_ces = lse_743972_multi.get_site_info_for_specie_allces(specie=Species("Si", 4))
        assert sorted(si_ces) == ["T:4", "T:5"]
        assert si_ces["T:4"]["isites"] == list(range(12, 20))
        assert si_ces["T:5"]["isites"] == [16, 17, 18, 18, 19, 19]
        si_ces = lse_743972_multi.get_site_info_for_specie_allces(specie=Species("Si", 4), min_fraction=0.1)
        assert list(si_ces) == ["T:4"]
        assert si_ces["T:4"]["fractions"] == approx([1, 1, 1, 1, 0.985, 0.985, 0.979, 0.979], abs=1e-3)
        # O environments were not computed
        assert lse_743972_multi.get_site_info_for_specie_allces(specie=Species("O", -2)) == {}

    def test_clear_environments(self, lse_5020, lse_743972_multi):
        lse = lse_743972_multi
        with pytest.raises(ValueError, match="Coordination environments have not been determined for site 20"):
            lse.site_has_clear_environment(isite=20)
        assert lse.site_has_clear_environment(isite=3)
        assert not lse.site_has_clear_environment(isite=2)

        conds_frac = [{"target": "ce_fraction", "minvalue": 0.99}]
        conds_csm = [{"target": "csm", "maxvalue": 1.0}]
        conds_nces = [{"target": "number_of_ces", "maxnumber": 2}]
        assert lse.site_has_clear_environment(isite=2, conditions=conds_frac)
        assert not lse.site_has_clear_environment(isite=16, conditions=conds_frac)
        assert lse.site_has_clear_environment(isite=3, conditions=conds_csm)
        assert not lse.site_has_clear_environment(isite=2, conditions=conds_csm)
        assert lse.site_has_clear_environment(isite=16, conditions=conds_nces)
        assert not lse.site_has_clear_environment(isite=18, conditions=conds_nces)
        with pytest.raises(ValueError, match="Target 'wrong' for condition of clear environment is not allowed"):
            lse.site_has_clear_environment(isite=3, conditions=[{"target": "wrong"}])

        # Sites 0 and 1 have empty CE lists, sites >= 20 have none
        assert not lse.structure_has_clear_environments()
        assert not lse.structure_has_clear_environments(skip_none=False, skip_empty=True)
        assert not lse.structure_has_clear_environments(skip_empty=True)
        conds_loose = [{"target": "ce_fraction", "minvalue": 0.4}]
        assert lse.structure_has_clear_environments(skip_empty=True, conditions=conds_loose)
        assert not lse.structure_has_clear_environments(skip_none=False, skip_empty=True, conditions=conds_loose)
        assert lse_5020.structure_has_clear_environments()

        assert sorted(lse.clear_environments()) == ["O:6", "T:4"]
        assert sorted(lse.clear_environments(conditions=conds_frac)) == ["BS_1:10", "O:6", "T:4"]
        assert sorted(lse_5020.clear_environments()) == ["C:12", "O:6"]

    def test_environments_identified(self, lse_5020, lse_743972_multi):
        assert lse_5020.environments_identified() == {"C:12", "O:6"}
        assert lse_743972_multi.environments_identified() == {
            "BO_1:8",
            "BO_2:8",
            "BS_1:10",
            "FO:7",
            "O:6",
            "PA:10",
            "PB:7",
            "ST:7",
            "T:4",
            "T:5",
        }
        assert lse_5020.structure_contains_atom_environment("Ba", "C:12")
        assert not lse_5020.structure_contains_atom_environment("Ba", "O:6")
        assert lse_5020.uniquely_determines_coordination_environments
        assert not lse_743972_multi.uniquely_determines_coordination_environments

    def test_eq(self, se_5020, lse_5020):
        assert lse_5020 != "not a LightStructureEnvironments"
        other = LightStructureEnvironments.from_structure_environments(
            strategy=SimplestChemenvStrategy(),
            structure_environments=se_5020,
            valences=[2, 4, -2, -2, -2],
        )
        assert lse_5020 == other
        other.valences = "undefined"
        assert lse_5020 != other
        other.valences = lse_5020.valences
        other._all_nbs_sites = other._all_nbs_sites[:-1]
        assert lse_5020 != other

    @pytest.mark.parametrize("lse_name", ["lse_5020", "lse_743972_multi"])
    def test_dict_json_round_trip(self, request, lse_name):
        lse = request.getfixturevalue(lse_name)
        dct = lse.as_dict()
        assert dct["@class"] == "LightStructureEnvironments"
        new_lse = LightStructureEnvironments.from_dict(json.loads(json.dumps(dct, cls=MontyEncoder)))
        assert new_lse == lse
        assert new_lse.strategy == lse.strategy
        for site_idx, nb_sets in enumerate(lse.neighbors_sets):
            if nb_sets is None:
                assert new_lse.neighbors_sets[site_idx] is None
                continue
            for nb_set, new_nb_set in zip(nb_sets, new_lse.neighbors_sets[site_idx], strict=True):
                assert_allclose(new_nb_set.neighb_coords, nb_set.neighb_coords)
                assert [nb["index"] for nb in new_nb_set.neighb_indices_and_images] == [
                    nb["index"] for nb in nb_set.neighb_indices_and_images
                ]

    def test_from_dict_without_image_cell(self, lse_5020):
        dct = orjson.loads(orjson.dumps(lse_5020.as_dict()))
        for nb_site in dct["all_nbs_sites"]:
            del nb_site["image_cell"]
        new_lse = LightStructureEnvironments.from_dict(dct)
        for new_nb, nb in zip(new_lse._all_nbs_sites, lse_5020._all_nbs_sites, strict=True):
            assert_array_equal(new_nb["image_cell"], nb["image_cell"])

        dct["all_nbs_sites"][0]["site"]["abc"][0] += 0.25
        with pytest.raises(ValueError, match="differences between one site in a periodic image cell is not integer"):
            LightStructureEnvironments.from_dict(dct)


class TestChemicalEnvironments:
    def test_init_and_container(self):
        with pytest.raises(NotImplementedError, match="coord_geoms argument is not yet implemented"):
            ChemicalEnvironments(coord_geoms={"T:4": {}})
        ce = ChemicalEnvironments()
        assert len(ce) == 0
        assert ce.minimum_geometry() is None
        assert ce.minimum_geometries() == []
        assert str(ce) == "Chemical environments object :\n => No coordination in it <=\n"

        ce = _make_ce("T:4", csm=0.5)
        ce.add_coord_geom("S:4", 2.0, permutation=[0, 1, 2, 3], other_symmetry_measures=dict.fromkeys(CSM_TYPES, 2.0))
        assert len(ce) == 2
        assert ce["S:4"]["symmetry_measure"] == approx(2.0)
        assert ce["S:4"]["algo"] == "UNKNOWN"
        assert [symbol for symbol, _ in ce] == ["T:4", "S:4"]
        assert ce.minimum_geometry(max_csm=0.1) is None
        assert ce.minimum_geometry(max_csm=1.0)[0] == "T:4"
        assert [symbol for symbol, _ in ce.minimum_geometries(n=1, max_csm=5.0)] == ["T:4"]
        assert [symbol for symbol, _ in ce.minimum_geometries(max_csm=5.0)] == ["T:4", "S:4"]
        assert " => Coordination 4 <=" in str(ce)

    def test_add_coord_geom_errors(self):
        ce = _make_ce("T:4")
        with pytest.raises(ChemenvError, match="Coordination geometry with mp_symbol 'XX:4' is not valid"):
            ce.add_coord_geom("XX:4", 1.0, permutation=[0, 1, 2, 3])
        with pytest.raises(ChemenvError, match="already present and override is set to False"):
            ce.add_coord_geom("T:4", 1.0, permutation=[0, 1, 2, 3])
        ce.add_coord_geom("T:4", 1.0, permutation=[3, 2, 1, 0], override=True)
        assert ce["T:4"]["symmetry_measure"] == approx(1.0)
        assert ce["T:4"]["permutation"] == [3, 2, 1, 0]

    def test_is_close_to_and_eq(self):
        ce = _make_ce("T:4", csm=0.5)
        assert ce == _make_ce("T:4", csm=0.5)
        assert ce != "not a ChemicalEnvironments"
        assert ce != _make_ce("S:4", csm=0.5)
        assert not ce.is_close_to(_make_ce("S:4", csm=0.5))
        assert ce.is_close_to(_make_ce("T:4", csm=0.5 + 1e-10))
        assert ce != _make_ce("T:4", csm=0.5 + 1e-10)
        assert ce != _make_ce("T:4", algo="SEPARATION_PLANE")
        assert ce != _make_ce("T:4", permutation=[1, 0, 2, 3])
        assert ce != _make_ce("T:4", detailed_voronoi_index={"index": 1, "cn": 4})
        other = _make_ce("T:4", csm=0.5)
        other.coord_geoms["T:4"]["other_symmetry_measures"]["csm_wocs_csc"] = 0.6
        assert ce != other
        assert not ce.is_close_to(other, atol=0.01)
        assert ce.is_close_to(other, atol=0.2)

    def test_dict_round_trip(self, se_7000):
        ce = se_7000.ce_list[6][4][0]
        new_ce = ChemicalEnvironments.from_dict(orjson.loads(orjson.dumps(ce.as_dict())))
        assert new_ce == ce
        assert new_ce["T:4"]["local2perfect_map"] == {0: 0, 1: 1, 2: 2, 3: 3}

        # Missing maps and other symmetry measures
        dct = _make_ce("T:4").as_dict()
        del dct["coord_geoms"]["T:4"]["other_symmetry_measures"]
        new_ce = ChemicalEnvironments.from_dict(dct)
        assert new_ce["T:4"]["local2perfect_map"] is None
        assert new_ce["T:4"]["perfect2local_map"] is None
        assert new_ce["T:4"]["other_symmetry_measures"] is None
