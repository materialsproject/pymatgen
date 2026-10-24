from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import orjson
import pytest
from matplotlib.figure import Figure
from pytest import approx

from pymatgen.analysis.chemenv.coordination_environments import chemenv_strategies
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
    AdditionalConditionInt,
    AngleCutoffFloat,
    AngleNbSetWeight,
    CNBiasNbSetWeight,
    CSMFloat,
    DeltaCSMNbSetWeight,
    DistanceAngleAreaNbSetWeight,
    DistanceCutoffFloat,
    MultiWeightsChemenvStrategy,
    NormalizedAngleDistanceNbSetWeight,
    SelfCSMNbSetWeight,
    SimpleAbundanceChemenvStrategy,
    SimplestChemenvStrategy,
    TargetedPenaltiedAbundanceChemenvStrategy,
    WeightedNbSetChemenvStrategy,
)
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import StructureEnvironments
from pymatgen.analysis.chemenv.utils.chemenv_errors import EquivalentSiteSearchError
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

__author__ = "waroquiers"

SE_DIR = TEST_FILES_DIR / "analysis" / "chemenv" / "structure_environments"
SI_SITES = (6, 7, 8)  # Tetrahedral Si sites in SiO2 (mp-7000); sites 0-5 are O


def _load_se(mp_id: str) -> StructureEnvironments:
    with open(SE_DIR / f"se_{mp_id}.json", "rb") as file:
        return StructureEnvironments.from_dict(orjson.loads(file.read()))


@pytest.fixture(scope="module")
def se_sio2() -> StructureEnvironments:
    return _load_se("mp-7000")


@pytest.fixture(scope="module")
def se_batio3() -> StructureEnvironments:
    return _load_se("mp-5020")


class TestStrategyOptions(MatSciTest):
    def test_options(self):
        # DistanceCutoffFloat
        with pytest.raises(ValueError, match=r"Distance cutoff should be between 1 and \+infinity"):
            DistanceCutoffFloat(0.5)
        dc1 = DistanceCutoffFloat(1.2)
        dc1_dict = dc1.as_dict()
        dc2 = DistanceCutoffFloat.from_dict(dc1_dict)
        assert dc1 == dc2

        # AngleCutoffFloat
        with pytest.raises(ValueError, match="Angle cutoff should be between 0 and 1, got 1.2"):
            AngleCutoffFloat(1.2)
        ac1 = AngleCutoffFloat(0.3)
        ac1_dict = ac1.as_dict()
        ac2 = AngleCutoffFloat.from_dict(ac1_dict)
        assert ac1 == ac2

        # CSMFloat
        with pytest.raises(
            ValueError,
            match="Continuous symmetry measure limits should be between 0 and 100, got 100.1",
        ):
            CSMFloat(100.1)
        csm1 = CSMFloat(0.458)
        csm1_dict = csm1.as_dict()
        csm2 = CSMFloat.from_dict(csm1_dict)
        assert csm1 == csm2

        # AdditionalConditions
        with pytest.raises(ValueError, match="Additional condition 5 is not allowed"):
            AdditionalConditionInt(5)
        with pytest.raises(ValueError, match="Additional condition 0.458 is not an integer"):
            AdditionalConditionInt(0.458)
        acd1 = AdditionalConditionInt(3)
        acd1_dict = acd1.as_dict()
        acd2 = AdditionalConditionInt.from_dict(acd1_dict)
        assert acd1 == acd2

    def test_strategies(self):
        simplest_strategy = SimplestChemenvStrategy()
        assert simplest_strategy.uniquely_determines_coordination_environments
        assert simplest_strategy.continuous_symmetry_measure_cutoff == 10
        assert simplest_strategy.distance_cutoff == approx(1.4)
        assert simplest_strategy.angle_cutoff == approx(0.3)

        simplest_strategy = SimplestChemenvStrategy(
            distance_cutoff=1.3,
            angle_cutoff=0.45,
            continuous_symmetry_measure_cutoff=8.26,
        )
        assert simplest_strategy.continuous_symmetry_measure_cutoff == approx(8.26)
        assert simplest_strategy.distance_cutoff == approx(1.3)
        assert simplest_strategy.angle_cutoff == approx(0.45)

        simplest_strategy.set_option("distance_cutoff", 1.5)
        assert simplest_strategy.distance_cutoff == approx(1.5)

        with pytest.raises(ValueError, match=r"Distance cutoff should be between 1 and \+infinity"):
            simplest_strategy.set_option("distance_cutoff", 0.5)

        simplest_strategy.set_option("angle_cutoff", 0.2)
        assert simplest_strategy.angle_cutoff == approx(0.2)

        with pytest.raises(ValueError, match="Angle cutoff should be between 0 and 1, got 1.5"):
            simplest_strategy.set_option("angle_cutoff", 1.5)

        simplest_strategy.setup_options(
            {
                "distance_cutoff": 1.4,
                "additional_condition": 3,
                "continuous_symmetry_measure_cutoff": 8.5,
            }
        )
        assert simplest_strategy.distance_cutoff == approx(1.4)
        assert simplest_strategy.continuous_symmetry_measure_cutoff == approx(8.5)
        assert simplest_strategy.additional_condition == 3

        with pytest.raises(
            ValueError,
            match="Continuous symmetry measure limits should be between 0 and 100, got -0.1",
        ):
            simplest_strategy.setup_options({"continuous_symmetry_measure_cutoff": -0.1})

        with pytest.raises(
            ValueError,
            match="Continuous symmetry measure limits should be between 0 and 100, got 100.1",
        ):
            simplest_strategy.setup_options({"continuous_symmetry_measure_cutoff": 100.1})


class TestSimplestChemenvStrategy:
    @pytest.mark.parametrize("isite", SI_SITES)
    def test_tetrahedral_si(self, se_sio2, isite):
        strategy = SimplestChemenvStrategy(structure_environments=se_sio2)
        site = se_sio2.structure[isite]

        ce_symbol, ce_dict = strategy.get_site_coordination_environment(site)
        assert ce_symbol == "T:4"
        assert ce_dict["symmetry_measure"] == approx(0.0098878, abs=1e-6)

        neighbors = strategy.get_site_neighbors(site)
        assert len(neighbors) == 4
        assert {nb.specie.symbol for nb in neighbors} == {"O"}
        assert sorted(site.distance_from_point(nb.coords) for nb in neighbors) == approx(
            [1.6255, 1.6255, 1.6284, 1.6284], abs=1e-4
        )

        assert strategy.get_site_coordination_environments(site) == [(ce_symbol, ce_dict)]
        assert strategy.get_site_coordination_environments(site, return_maps=True) == [((ce_symbol, ce_dict), (4, 0))]

        fractions = strategy.get_site_coordination_environments_fractions(site, return_strategy_dict_info=True)
        assert fractions == [
            {"ce_symbol": "T:4", "ce_dict": ce_dict, "ce_fraction": 1, "ce_map": (4, 0), "strategy_info": {}}
        ]
        assert strategy.get_site_coordination_environments_fractions(site, return_maps=False) == [
            {"ce_symbol": "T:4", "ce_dict": ce_dict, "ce_fraction": 1}
        ]

    def test_ce_fractions_and_neighbors(self, se_sio2):
        strategy = SimplestChemenvStrategy(structure_environments=se_sio2)
        site = se_sio2.structure[6]
        (ce_and_nbs,) = strategy.get_site_ce_fractions_and_neighbors(site)
        assert ce_and_nbs["ce_symbol"] == "T:4"
        assert ce_and_nbs["ce_fraction"] == 1
        assert ce_and_nbs["ce_map"] == (4, 0)
        assert len(ce_and_nbs["neighbors"]) == 4
        for nb in ce_and_nbs["neighbors"]:
            assert set(nb) == {"site", "index"}
            assert se_sio2.structure[nb["index"]].specie.symbol == "O"
            assert site.distance_from_point(nb["site"].coords) == approx(1.627, abs=2e-3)

    @pytest.mark.parametrize(
        ("distance_cutoff", "angle_cutoff", "ce_symbol", "cn_map"),
        [
            (1.0, 0.3, "A:2", (2, 0)),
            (1.4, 1.0, "A:2", (2, 1)),
            (2.0, 0.3, "T:4", (4, 0)),
            (10.0, 0.0, "PP:6", (6, 0)),
        ],
    )
    def test_cutoffs(self, se_sio2, distance_cutoff, angle_cutoff, ce_symbol, cn_map):
        strategy = SimplestChemenvStrategy(
            structure_environments=se_sio2, distance_cutoff=distance_cutoff, angle_cutoff=angle_cutoff
        )
        site = se_sio2.structure[6]
        ce, ce_map = strategy.get_site_coordination_environment(site, return_map=True)
        assert ce[0] == ce_symbol
        assert ce_map == cn_map
        assert len(strategy.get_site_neighbors(site)) == cn_map[0]

    def test_no_neighbors_set_for_cutoffs(self, se_sio2):
        # The smallest distance cutoff together with the largest angle cutoff only selects the
        # empty neighbors set, for which no coordination environment is available
        strategy = SimplestChemenvStrategy(structure_environments=se_sio2, distance_cutoff=1.0, angle_cutoff=1.0)
        site = se_sio2.structure[6]
        assert strategy.get_site_coordination_environment(site) is None
        assert strategy.get_site_coordination_environments_fractions(site) is None
        assert strategy.get_site_ce_fractions_and_neighbors(site) is None

    def test_site_without_environments(self, se_sio2):
        # Environments were only computed for the Si cations
        strategy = SimplestChemenvStrategy(structure_environments=se_sio2)
        assert strategy.get_site_coordination_environments_fractions(se_sio2.structure[0]) is None

    @pytest.mark.parametrize(("isite", "ce_symbol"), [(0, "C:12"), (1, "O:6")])
    def test_batio3(self, se_batio3, isite, ce_symbol):
        strategy = SimplestChemenvStrategy(structure_environments=se_batio3)
        site = se_batio3.structure[isite]
        assert strategy.get_site_coordination_environment(site)[0] == ce_symbol
        neighbors = strategy.get_site_neighbors(site)
        assert len(neighbors) == int(ce_symbol.split(":")[1])
        assert {nb.specie.symbol for nb in neighbors} == {"O"}

    def test_periodic_image(self, se_sio2):
        strategy = SimplestChemenvStrategy(structure_environments=se_sio2)
        site = se_sio2.structure[6]
        translation = np.array([1, 0, -1])
        image = site.__class__(site.species, site.frac_coords + translation, site.lattice)
        isite, dequivsite, dthissite, sym = strategy.equivalent_site_index_and_transform(image)
        assert isite == 6
        assert dequivsite == approx([0, 0, 0])
        assert dthissite == approx(translation)
        assert sym.affine_matrix == approx(np.eye(4))
        # The neighbors of the image are the translated neighbors of the site
        nbs = sorted(tuple(np.round(nb.frac_coords, 6)) for nb in strategy.get_site_neighbors(site))
        image_nbs = sorted(
            tuple(np.round(nb.frac_coords - translation, 6)) for nb in strategy.get_site_neighbors(image)
        )
        np.testing.assert_allclose(image_nbs, nbs, atol=1e-5)

    def test_symmetry_equivalent_site(self):
        """Neighbors of a site are obtained from those of a symmetry-equivalent site when only the latter
        has been computed (as done when only symmetry-independent sites are analyzed).
        """
        se = _load_se("mp-7000")
        reference = SimplestChemenvStrategy(structure_environments=_load_se("mp-7000"))
        se.sites_map[7] = 6
        strategy = SimplestChemenvStrategy(structure_environments=se)
        site = se.structure[7]

        isite, _dequivsite, _dthissite, sym = strategy.equivalent_site_index_and_transform(site)
        assert isite == 6
        assert sym.affine_matrix != approx(np.eye(4))
        assert strategy.get_site_coordination_environment(site)[0] == "T:4"

        def cart_coords(sites):
            return sorted(tuple(np.round(nb.coords, 4)) for nb in sites)

        neighbors = strategy.get_site_neighbors(site)
        assert sorted(site.distance_from_point(nb.coords) for nb in neighbors) == approx(
            [1.6255, 1.6255, 1.6284, 1.6284], abs=1e-4
        )
        np.testing.assert_allclose(cart_coords(neighbors), cart_coords(reference.get_site_neighbors(site)), atol=1e-3)

    def test_equivalent_site_not_found(self):
        se = _load_se("mp-7000")
        se.sites_map[0] = 6  # No symmetry operation maps a Si site onto an O site
        strategy = SimplestChemenvStrategy(structure_environments=se)
        with pytest.raises(EquivalentSiteSearchError, match="Equivalent site could not be found"):
            strategy.equivalent_site_index_and_transform(se.structure[0])

    def test_prepare_symmetries(self, se_sio2, monkeypatch):
        strategy = SimplestChemenvStrategy(structure_environments=se_sio2)
        assert len(strategy.symops) > 1
        spg_analyzer = strategy.spg_analyzer
        strategy.prepare_symmetries()  # Cached: nothing is recomputed
        assert strategy.spg_analyzer is spg_analyzer

        def failing_analyzer(*args, **kwargs):
            raise RuntimeError("symmetry analysis failed")

        monkeypatch.setattr(chemenv_strategies, "SpacegroupAnalyzer", failing_analyzer)
        strategy = SimplestChemenvStrategy(structure_environments=se_sio2)
        assert strategy.symops == []
        assert strategy._symops_cache_fp is None

    def test_wrong_voronoi_container(self):
        with pytest.raises(TypeError, match='Voronoi Container not of type "DetailedVoronoiContainer"'):
            SimplestChemenvStrategy(structure_environments=SimpleNamespace(voronoi=None))

    def test_str(self):
        strategy = SimplestChemenvStrategy(distance_cutoff=1.5, angle_cutoff=0.2)
        out = str(strategy)
        assert "Chemenv Strategy 'SimplestChemenvStrategy'" in out
        assert SimplestChemenvStrategy.STRATEGY_DESCRIPTION in out
        assert "   - distance_cutoff : 1.5\n" in out
        assert "   - angle_cutoff : 0.2\n" in out
        assert "   - additional_condition : 1\n" in out
        assert "   - continuous_symmetry_measure_cutoff : 10.0\n" in out

    def test_as_dict_from_dict(self):
        strategy = SimplestChemenvStrategy(
            distance_cutoff=1.6,
            angle_cutoff=0.25,
            additional_condition=3,
            continuous_symmetry_measure_cutoff=5.5,
            symmetry_measure_type="csm_wocs_ctwocc",
        )
        dct = json.loads(json.dumps(strategy.as_dict()))
        assert dct["distance_cutoff"] == approx(1.6)
        assert dct["symmetry_measure_type"] == "csm_wocs_ctwocc"
        new_strategy = SimplestChemenvStrategy.from_dict(dct)
        assert new_strategy == strategy
        assert new_strategy.additional_condition == 3

    def test_eq(self):
        strategy = SimplestChemenvStrategy()
        assert strategy == SimplestChemenvStrategy()
        assert strategy != SimplestChemenvStrategy(symmetry_measure_type="csm_wocs_ctwocc")
        assert strategy != SimplestChemenvStrategy(continuous_symmetry_measure_cutoff=5)
        assert strategy.__eq__(MultiWeightsChemenvStrategy()) is NotImplemented

    def test_add_strategy_visualization_to_subplot(self):
        strategy = SimplestChemenvStrategy(distance_cutoff=1.6, angle_cutoff=0.25)
        ax = Figure().add_subplot()
        strategy.add_strategy_visualization_to_subplot(ax)
        assert len(ax.lines) == 2
        for line in ax.lines:
            assert list(line.get_xdata()) == approx([1.6])
            assert list(line.get_ydata()) == approx([0.25])
        assert [line.get_marker() for line in ax.lines] == ["o", "x"]


@pytest.mark.parametrize(
    "cls",
    [SimpleAbundanceChemenvStrategy, TargetedPenaltiedAbundanceChemenvStrategy],
)
def test_abundance_strategies_not_instantiable(cls):
    # These strategies do not implement get_site_coordination_environments_fractions (and their
    # constructors raise NotImplementedError anyway)
    with pytest.raises(TypeError, match="Can't instantiate abstract class"):
        cls()


def _angle_weight(nb_set) -> float:
    """Hand-computed AngleNbSetWeight(aa=1): solid angles summed and normalized by the full sphere."""
    return sum(nb_set.angles) / (4 * np.pi)


class TestWeightedNbSetChemenvStrategy:
    def test_missing_weights(self):
        with pytest.raises(ValueError, match="nb_set_weights=None must be provided"):
            WeightedNbSetChemenvStrategy()

    def test_not_unique(self):
        strategy = WeightedNbSetChemenvStrategy(nb_set_weights=[AngleNbSetWeight()])
        assert strategy.uniquely_determines_coordination_environments is False
        assert strategy.get_site_coordination_environment(site=None) is None
        assert strategy.get_site_neighbors(site=None) is None

    def test_fractions(self, se_sio2):
        weights = [AngleNbSetWeight(), SelfCSMNbSetWeight()]
        strategy = WeightedNbSetChemenvStrategy(structure_environments=se_sio2, nb_set_weights=weights)
        site = se_sio2.structure[6]
        fractions = strategy.get_site_coordination_environments_fractions(site, return_strategy_dict_info=True)

        # Hand-computed fractions of the neighbors sets from the product of the weights
        products = {}
        for cn, nb_sets in se_sio2.neighbors_sets[6].items():
            for inb, nb_set in enumerate(nb_sets):
                self_csm = weights[1].weight(nb_set, se_sio2, cn_map=(cn, inb), additional_info={})
                products[cn, inb] = _angle_weight(nb_set) * self_csm
        total = sum(products.values())

        assert [(frac["ce_symbol"], frac["ce_map"]) for frac in fractions] == [
            ("T:4", (4, 0)),
            ("A:2", (2, 1)),
            ("A:2", (2, 0)),
        ]
        assert sum(frac["ce_fraction"] for frac in fractions) == approx(1)
        for frac in fractions:
            info = frac["strategy_info"]
            cn, inb = frac["ce_map"]
            assert info["AngleWeight"] == approx(_angle_weight(se_sio2.neighbors_sets[6][cn][inb]))
            assert info["NbSetFraction"] == approx(products[frac["ce_map"]] / total)
            assert info["Product"] == approx(info["AngleWeight"] * info["SelfCSMWeight"])
            assert info["Fraction"] == approx(info["NbSetFraction"] * info["CEFraction"])
            assert frac["ce_fraction"] == approx(info["Fraction"])
        assert fractions[0]["ce_fraction"] == approx(0.5470581045602808)

        # Unordered, without maps, and with a minimum fraction
        unordered = strategy.get_site_coordination_environments_fractions(site, ordered=False, return_maps=False)
        assert sorted(frac["ce_fraction"] for frac in unordered) == approx(
            sorted(frac["ce_fraction"] for frac in fractions)
        )
        assert all("ce_map" not in frac for frac in unordered)
        major = strategy.get_site_coordination_environments_fractions(
            site, min_fraction=0.3, return_strategy_dict_info=True
        )
        assert [(frac["ce_symbol"], frac["ce_map"]) for frac in major] == [("T:4", (4, 0))]
        assert major[0]["strategy_info"] == fractions[0]["strategy_info"]

    def test_fractions_zero_weights(self, se_sio2):
        cn_bias = CNBiasNbSetWeight(cn_weights=dict.fromkeys(range(14), 0), initialization_options={"type": "explicit"})
        strategy = WeightedNbSetChemenvStrategy(structure_environments=se_sio2, nb_set_weights=[cn_bias])
        assert strategy.get_site_coordination_environments_fractions(se_sio2.structure[6]) == []

    def test_site_without_environments(self, se_sio2):
        strategy = WeightedNbSetChemenvStrategy(structure_environments=se_sio2, nb_set_weights=[AngleNbSetWeight()])
        assert strategy.get_site_coordination_environments_fractions(se_sio2.structure[0]) is None

    @pytest.mark.xfail(
        reason="get_site_coordination_environment() of WeightedNbSetChemenvStrategy only accepts `site`",
        raises=TypeError,
        strict=True,
    )
    def test_get_site_coordination_environments(self, se_sio2):
        strategy = WeightedNbSetChemenvStrategy(structure_environments=se_sio2, nb_set_weights=[AngleNbSetWeight()])
        strategy.get_site_coordination_environments(se_sio2.structure[6])

    @pytest.mark.xfail(
        reason="Crashes when no geometry of a weighted neighbors set is below the max_csm of the ce_estimator",
        raises=TypeError,
        strict=True,
    )
    def test_fractions_unclear_environment(self, se_sio2):
        # The 6-neighbors set of Si only has geometries with csm > 8 (PP:6 at 22.3) and gets a non-zero angle weight
        strategy = WeightedNbSetChemenvStrategy(structure_environments=se_sio2, nb_set_weights=[AngleNbSetWeight()])
        strategy.get_site_coordination_environments_fractions(se_sio2.structure[6])

    def test_as_dict_from_dict(self):
        strategy = WeightedNbSetChemenvStrategy(
            additional_condition=2,
            symmetry_measure_type="csm_wocs_ctwocc",
            nb_set_weights=[AngleNbSetWeight(aa=2), SelfCSMNbSetWeight()],
            ce_estimator={"function": "power2_inverse_power2_decreasing", "options": {"max_csm": 5.0}},
        )
        dct = json.loads(json.dumps(strategy.as_dict()))
        new_strategy = WeightedNbSetChemenvStrategy.from_dict(dct)
        assert new_strategy == strategy
        assert [weight["name"] for weight in new_strategy.ordered_weights] == ["AngleWeight", "SelfCSMWeight"]

    def test_eq(self):
        strategy = WeightedNbSetChemenvStrategy(nb_set_weights=[AngleNbSetWeight()])
        assert strategy == WeightedNbSetChemenvStrategy(nb_set_weights=[AngleNbSetWeight()])
        assert strategy != WeightedNbSetChemenvStrategy(nb_set_weights=[AngleNbSetWeight(aa=2)])
        assert strategy != WeightedNbSetChemenvStrategy(nb_set_weights=[AngleNbSetWeight()], additional_condition=0)
        assert strategy.__eq__(SimplestChemenvStrategy()) is NotImplemented


class TestMultiWeightsChemenvStrategy:
    @pytest.fixture(scope="class")
    def all_weights_strategy(self):
        return MultiWeightsChemenvStrategy(
            additional_condition=2,
            dist_ang_area_weight=DistanceAngleAreaNbSetWeight(),
            self_csm_weight=SelfCSMNbSetWeight(),
            delta_csm_weight=DeltaCSMNbSetWeight.delta_cn_specifics(),
            cn_bias_weight=CNBiasNbSetWeight.linearly_equidistant(weight_cn1=1, weight_cn13=2),
            angle_weight=AngleNbSetWeight(aa=2),
            normalized_angle_distance_weight=NormalizedAngleDistanceNbSetWeight(average_type="geometric", aa=1, bb=1),
        )

    def test_ordered_weights(self, all_weights_strategy):
        assert [weight["name"] for weight in all_weights_strategy.ordered_weights] == [
            "DistAngleAreaWeight",
            "SelfCSMWeight",
            "DeltaCSMWeight",
            "CNBiasWeight",
            "AngleWeight",
            "NormAngleDistWeight",
        ]
        assert all_weights_strategy.uniquely_determines_coordination_environments is False

    @pytest.mark.parametrize("use_all_weights", [True, False])
    def test_as_dict_from_dict(self, all_weights_strategy, use_all_weights):
        strategy = all_weights_strategy if use_all_weights else MultiWeightsChemenvStrategy()
        dct = json.loads(json.dumps(strategy.as_dict()))
        new_strategy = MultiWeightsChemenvStrategy.from_dict(dct)
        assert new_strategy == strategy
        assert new_strategy.as_dict() == strategy.as_dict()

    def test_eq(self, all_weights_strategy):
        stats_strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
        assert stats_strategy == MultiWeightsChemenvStrategy.stats_article_weights_parameters()
        assert stats_strategy != all_weights_strategy
        assert stats_strategy.__eq__(SimplestChemenvStrategy()) is NotImplemented

    @pytest.mark.parametrize("isite", SI_SITES)
    def test_stats_article_fractions(self, se_sio2, isite):
        strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
        strategy.set_structure_environments(se_sio2)
        site = se_sio2.structure[isite]

        (fraction,) = strategy.get_site_coordination_environments_fractions(site, return_strategy_dict_info=True)
        assert fraction["ce_symbol"] == "T:4"
        assert fraction["ce_fraction"] == approx(1)
        assert fraction["ce_map"] == (4, 0)
        info = fraction["strategy_info"]
        assert info["DistAngleAreaWeight"] == 1
        assert info["DeltaCSMWeight"] == approx(1)
        assert info["SelfCSMWeight"] == approx(0.9962974, abs=1e-6)
        assert info["Product"] == approx(info["SelfCSMWeight"])

        (ce_and_nbs,) = strategy.get_site_ce_fractions_and_neighbors(site)
        assert ce_and_nbs["ce_symbol"] == "T:4"
        assert len(ce_and_nbs["neighbors"]) == 4

    def test_stats_article_return_all(self, se_sio2):
        strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
        strategy.set_structure_environments(se_sio2)
        fractions = strategy.get_site_coordination_environments_fractions(
            se_sio2.structure[6], return_all=True, return_strategy_dict_info=True
        )
        # All neighbors sets (except the empty one) are returned, even those with zero weight
        assert {frac["ce_map"] for frac in fractions} == {(2, 0), (2, 1), (4, 0), (6, 0)}
        assert sum(frac["ce_fraction"] for frac in fractions) == approx(1)
        assert fractions[0]["ce_symbol"] == "T:4"
        assert fractions[0]["ce_fraction"] == approx(1)
        # All geometries of the 6-neighbors set are above max_csm: environment is unclear
        (unclear,) = (frac for frac in fractions if frac["ce_map"] == (6, 0))
        assert unclear["ce_symbol"] == "UNCLEAR:6"
        assert unclear["ce_dict"] is None
        assert unclear["ce_fraction"] == 0
        assert unclear["strategy_info"]["CEFraction"] is None
        symbols_4 = {frac["ce_symbol"] for frac in fractions if frac["ce_map"] == (4, 0)}
        assert symbols_4 == {"T:4", "SS:4", "S:4", "SY:4"}
