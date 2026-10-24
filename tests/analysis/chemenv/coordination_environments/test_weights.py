from __future__ import annotations

import json

import orjson
import pytest
from pytest import approx

from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
    AngleNbSetWeight,
    AnglePlateauNbSetWeight,
    CNBiasNbSetWeight,
    DeltaCSMNbSetWeight,
    DeltaDistanceNbSetWeight,
    DistanceAngleAreaNbSetWeight,
    DistanceNbSetWeight,
    DistancePlateauNbSetWeight,
    NormalizedAngleDistanceNbSetWeight,
    SelfCSMNbSetWeight,
)
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import StructureEnvironments
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

__author__ = "waroquiers"

struct_env_dir = f"{TEST_FILES_DIR}/analysis/chemenv/structure_environments"


class FakeNbSet:
    def __init__(self, cn=None):
        self.cn = cn

    def __len__(self):
        return self.cn


class DummyStructureEnvironments:
    pass


class DummyVoronoiContainer:
    pass


class TestStrategyWeights(MatSciTest):
    def test_angle_weight(self):
        fake_nb_set = FakeNbSet()
        dummy_se = DummyStructureEnvironments()

        # Angles for a given fake nb_set with 5 neighbors
        fake_nb_set.angles = [
            1.8595833644514066,
            2.622518848090717,
            3.08570351705799,
            2.2695472184920042,
            2.2695338778592387,
        ]
        angle_weight = AngleNbSetWeight(aa=1)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.9634354419021528) < 1e-8
        angle_weight = AngleNbSetWeight(aa=2)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.92820785071319645) < 1e-8
        angle_weight = AngleNbSetWeight(aa=0.5)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.98154747307613843) < 1e-8

        assert AngleNbSetWeight(1) != AngleNbSetWeight(2)

        # nb_set with no neighbor
        fake_nb_set.angles = []
        angle_weight = AngleNbSetWeight(aa=1)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw) < 1e-8
        angle_weight = AngleNbSetWeight(aa=2)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw) < 1e-8

        # nb_set with one neighbor
        fake_nb_set.angles = [3.08570351705799]
        angle_weight = AngleNbSetWeight(aa=1)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.24555248382791284) < 1e-8
        angle_weight = AngleNbSetWeight(aa=2)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.060296022314057396) < 1e-8
        angle_weight = AngleNbSetWeight(aa=0.5)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.49553252549950022) < 1e-8

        # nb_set with 6 neighbors (sum of the angles is 4*pi, i.e. the full sphere)
        fake_nb_set.angles = [
            1.8595833644514066,
            0.459483788407816,
            2.622518848090717,
            3.08570351705799,
            2.2695472184920042,
            2.2695338778592387,
        ]
        angle_weight = AngleNbSetWeight(aa=1)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 1) < 1e-8
        angle_weight = AngleNbSetWeight(aa=2)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 1) < 1e-8
        angle_weight = AngleNbSetWeight(aa=0.5)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 1) < 1e-8

    def test_normalized_angle_distance_weight(self):
        fake_nb_set = FakeNbSet()
        dummy_se = DummyStructureEnvironments()

        nadw1 = NormalizedAngleDistanceNbSetWeight(average_type="geometric", aa=1, bb=1)
        nadw2 = NormalizedAngleDistanceNbSetWeight(average_type="arithmetic", aa=1, bb=1)
        nadw3 = NormalizedAngleDistanceNbSetWeight(average_type="geometric", aa=0, bb=1)
        nadw4 = NormalizedAngleDistanceNbSetWeight(average_type="arithmetic", aa=1, bb=0)
        nadw5 = NormalizedAngleDistanceNbSetWeight(average_type="arithmetic", aa=0.1, bb=0.1)
        nadw6 = NormalizedAngleDistanceNbSetWeight(average_type="arithmetic", aa=0, bb=0.1)
        nadw7 = NormalizedAngleDistanceNbSetWeight(average_type="arithmetic", aa=0.1, bb=0)
        nadw8 = NormalizedAngleDistanceNbSetWeight(average_type="arithmetic", aa=2, bb=0)
        nadw9 = NormalizedAngleDistanceNbSetWeight(average_type="arithmetic", aa=0, bb=2)
        nadw10 = NormalizedAngleDistanceNbSetWeight(average_type="arithmetic", aa=2, bb=2)
        nadw11 = NormalizedAngleDistanceNbSetWeight(average_type="geometric", aa=1, bb=2)
        nadw12 = NormalizedAngleDistanceNbSetWeight(average_type="geometric", aa=2, bb=1)
        assert nadw11 != nadw12
        with pytest.raises(ValueError, match="Both exponents are 0."):
            NormalizedAngleDistanceNbSetWeight(average_type="arithmetic", aa=0, bb=0)
        with pytest.raises(
            ValueError,
            match="Average type is 'arithmetix' while it should be 'geometric' or 'arithmetic'",
        ):
            NormalizedAngleDistanceNbSetWeight(average_type="arithmetix", aa=1, bb=1)

        fake_nb_set.normalized_distances = [
            1.2632574171572457,
            1.1231971151388764,
            1,
            1.1887986376446249,
            1.188805134890625,
        ]
        fake_nb_set.normalized_angles = [
            0.6026448601336767,
            0.8498933334305273,
            1,
            0.7355039801931018,
            0.7354996568248028,
        ]
        w1 = nadw1.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w1 - 0.67310887189488189) < 1e-8
        w2 = nadw2.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w2 - 0.69422258996523023) < 1e-8
        w3 = nadw3.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w3 - 0.8700949310182079) < 1e-8
        w4 = nadw4.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w4 - 0.7847083661164217) < 1e-8
        w5 = nadw5.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w5 - 0.96148050989126843) < 1e-8
        w6 = nadw6.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w6 - 0.98621181678741754) < 1e-8
        w7 = nadw7.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w7 - 0.97479580875402994) < 1e-8
        w8 = nadw8.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w8 - 0.63348507114489783) < 1e-8
        w9 = nadw9.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w9 - 0.7668954450583646) < 1e-8
        w10 = nadw10.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w10 - 0.51313920014833292) < 1e-8
        w11 = nadw11.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w11 - 0.585668617459) < 1e-8
        w12 = nadw12.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w12 - 0.520719679281) < 1e-8

    def test_CN_bias_weight(self):
        fake_nb_set = FakeNbSet()
        dummy_se = DummyStructureEnvironments()
        bias_weight1 = CNBiasNbSetWeight.linearly_equidistant(weight_cn1=1, weight_cn13=13)
        bias_weight2 = CNBiasNbSetWeight.geometrically_equidistant(weight_cn1=1, weight_cn13=1.1**12)
        bias_weight3 = CNBiasNbSetWeight.explicit(
            cn_weights={
                1: 1,
                2: 3,
                3: 3.2,
                4: 4,
                5: 4.1,
                6: 4.2,
                7: 4.3,
                8: 4.4,
                9: 4.5,
                10: 4.6,
                11: 4.6,
                12: 4.7,
                13: 4.8,
            }
        )
        with pytest.raises(ValueError, match="Weights should be provided for CN 1 to 13"):
            CNBiasNbSetWeight.explicit(cn_weights={1: 1, 13: 2})

        fake_nb_set.cn = 1
        w1 = bias_weight1.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w1 - 1) < 1e-8
        w2 = bias_weight2.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w2 - 1) < 1e-8
        w3 = bias_weight3.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w3 - 1) < 1e-8
        fake_nb_set.cn = 7
        w1 = bias_weight1.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w1 - 7) < 1e-8
        w2 = bias_weight2.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w2 - 1.1**6) < 1e-8
        w3 = bias_weight3.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w3 - 4.3) < 1e-8
        fake_nb_set.cn = 13
        w1 = bias_weight1.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w1 - 13) < 1e-8
        w2 = bias_weight2.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w2 - 1.1**12) < 1e-8
        w3 = bias_weight3.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w3 - 4.8) < 1e-8

        bias_weight4 = CNBiasNbSetWeight.from_description(
            {"type": "linearly_equidistant", "weight_cn1": 2, "weight_cn13": 26}
        )
        for idx in range(1, 14):
            assert bias_weight4.cn_weights[idx] == approx(2 * idx)

        bias_weight5 = CNBiasNbSetWeight.from_description(
            {"type": "geometrically_equidistant", "weight_cn1": 1, "weight_cn13": 13}
        )
        assert bias_weight5.cn_weights[1] == approx(1)
        assert bias_weight5.cn_weights[3] == approx(1.5334062370163877)
        assert bias_weight5.cn_weights[9] == approx(5.5287748136788739)
        assert bias_weight5.cn_weights[12] == approx(10.498197520079623)

        cn_weights = dict.fromkeys(range(1, 14), 0)
        cn_weights[6] = 2
        cn_weights[4] = 1
        bias_weight6 = CNBiasNbSetWeight.from_description({"type": "explicit", "cn_weights": cn_weights})

        assert bias_weight6.cn_weights[1] == approx(0)
        assert bias_weight6.cn_weights[4] == approx(1)
        assert bias_weight6.cn_weights[6] == approx(2)

    def test_self_csms_weight(self):
        # Get the StructureEnvironments for K2NaNb2Fe7Si8H4O31 (mp-743972)
        with open(f"{struct_env_dir}/se_mp-743972.json", "rb") as file:
            dct = orjson.loads(file.read())
        struct_envs = StructureEnvironments.from_dict(dct)

        # Get neighbors sets for which we get the weights
        cn_maps = [(12, 3), (12, 2), (13, 2), (12, 0), (12, 1)]
        nb_sets = {cn_map: struct_envs.neighbors_sets[0][cn_map[0]][cn_map[1]] for cn_map in cn_maps}

        effective_csm_estimator = {
            "function": "power2_inverse_decreasing",
            "options": {"max_csm": 8},
        }
        weight_estimator = {
            "function": "power2_decreasing_exp",
            "options": {"max_csm": 8, "alpha": 1},
        }
        weight_estimator2 = {
            "function": "power2_decreasing_exp",
            "options": {"max_csm": 8.1, "alpha": 1},
        }
        symmetry_measure_type = "csm_wcs_ctwcc"
        self_weight = SelfCSMNbSetWeight(
            effective_csm_estimator=effective_csm_estimator,
            weight_estimator=weight_estimator,
            symmetry_measure_type=symmetry_measure_type,
        )
        self_weight2 = SelfCSMNbSetWeight(
            effective_csm_estimator=effective_csm_estimator,
            weight_estimator=weight_estimator2,
            symmetry_measure_type=symmetry_measure_type,
        )
        assert self_weight != self_weight2

        additional_info = {}
        cn_map = (12, 3)
        self_w = self_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.11671945916431022) < 1e-8
        cn_map = (12, 2)
        self_w = self_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w) < 1e-8
        cn_map = (12, 0)
        self_w = self_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w) < 1e-8
        cn_map = (12, 1)
        self_w = self_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w) < 1e-8
        cn_map = (13, 2)
        self_w = self_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.14204073172729198) < 1e-8

        # Get the StructureEnvironments for SiO2 (mp-7000)
        with open(f"{struct_env_dir}/se_mp-7000.json", "rb") as file:
            dct = orjson.loads(file.read())
        struct_envs = StructureEnvironments.from_dict(dct)

        # Get neighbors sets for which we get the weights
        cn_maps = [(2, 0), (4, 0)]
        nb_sets = {cn_map: struct_envs.neighbors_sets[6][cn_map[0]][cn_map[1]] for cn_map in cn_maps}

        effective_csm_estimator = {
            "function": "power2_inverse_decreasing",
            "options": {"max_csm": 8},
        }

        weight_estimator = {
            "function": "power2_decreasing_exp",
            "options": {"max_csm": 8, "alpha": 1},
        }
        symmetry_measure_type = "csm_wcs_ctwcc"
        self_weight = SelfCSMNbSetWeight(
            effective_csm_estimator=effective_csm_estimator,
            weight_estimator=weight_estimator,
            symmetry_measure_type=symmetry_measure_type,
        )

        additional_info = {}
        cn_map = (2, 0)
        self_w = self_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.8143992162836029) < 1e-8
        cn_map = (4, 0)
        self_w = self_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.99629742352359496) < 1e-8

    def test_delta_csms_weight(self):
        # Get the StructureEnvironments for K2NaNb2Fe7Si8H4O31 (mp-743972)
        with open(f"{struct_env_dir}/se_mp-743972.json", "rb") as file:
            dct = orjson.loads(file.read())
        struct_envs = StructureEnvironments.from_dict(dct)

        # Get neighbors sets for which we get the weights
        cn_maps = [(12, 3), (12, 2), (13, 2), (12, 0), (12, 1), (13, 0), (13, 1)]
        nb_sets = {cn_map: struct_envs.neighbors_sets[0][cn_map[0]][cn_map[1]] for cn_map in cn_maps}

        effective_csm_estimator = {
            "function": "power2_inverse_decreasing",
            "options": {"max_csm": 8},
        }
        weight_estimator = {
            "function": "smootherstep",
            "options": {"delta_csm_min": 0.5, "delta_csm_max": 3},
        }
        symmetry_measure_type = "csm_wcs_ctwcc"
        delta_weight = DeltaCSMNbSetWeight(
            effective_csm_estimator=effective_csm_estimator,
            weight_estimator=weight_estimator,
            symmetry_measure_type=symmetry_measure_type,
        )

        additional_info = {}
        cn_map = (12, 3)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w) < 1e-8
        cn_map = (12, 2)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w) < 1e-8
        cn_map = (12, 0)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w) < 1e-8
        cn_map = (12, 1)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w) < 1e-8
        cn_map = (13, 2)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 1) < 1e-8
        cn_map = (13, 0)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w) < 1e-8
        cn_map = (13, 1)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w) < 1e-8

        effective_csm_estimator = {
            "function": "power2_inverse_decreasing",
            "options": {"max_csm": 8},
        }

        weight_estimator = {
            "function": "smootherstep",
            "options": {"delta_csm_min": -1, "delta_csm_max": 3},
        }
        symmetry_measure_type = "csm_wcs_ctwcc"
        delta_weight = DeltaCSMNbSetWeight(
            effective_csm_estimator=effective_csm_estimator,
            weight_estimator=weight_estimator,
            symmetry_measure_type=symmetry_measure_type,
        )

        additional_info = {}
        cn_map = (12, 3)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.040830741048481355) < 1e-8
        cn_map = (13, 2)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 1) < 1e-8
        cn_map = (13, 0)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.103515625) < 1e-8
        cn_map = (13, 1)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.103515625) < 1e-8

        # Get the StructureEnvironments for SiO2 (mp-7000)
        with open(f"{struct_env_dir}/se_mp-7000.json", "rb") as file:
            dct = orjson.loads(file.read())
        struct_envs = StructureEnvironments.from_dict(dct)

        # Get neighbors sets for which we get the weights
        cn_maps = [(2, 0), (4, 0)]
        nb_sets = {cn_map: struct_envs.neighbors_sets[6][cn_map[0]][cn_map[1]] for cn_map in cn_maps}

        effective_csm_estimator = {
            "function": "power2_inverse_decreasing",
            "options": {"max_csm": 8},
        }

        weight_estimator = {
            "function": "smootherstep",
            "options": {"delta_csm_min": 0.5, "delta_csm_max": 3},
        }
        symmetry_measure_type = "csm_wcs_ctwcc"
        delta_weight = DeltaCSMNbSetWeight(
            effective_csm_estimator=effective_csm_estimator,
            weight_estimator=weight_estimator,
            symmetry_measure_type=symmetry_measure_type,
        )

        additional_info = {}
        cn_map = (2, 0)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w) < 1e-8
        cn_map = (4, 0)
        delta_w = delta_weight.weight(
            nb_set=nb_sets[cn_map],
            structure_environments=struct_envs,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 1) < 1e-8

    def test_dist_angle_area_weight(self):
        surface_definition = {
            "type": "standard_elliptic",
            "distance_bounds": {"lower": 1.2, "upper": 1.8},
            "angle_bounds": {"lower": 0.2, "upper": 0.8},
        }
        da_area_weight = DistanceAngleAreaNbSetWeight(
            weight_type="has_intersection",
            surface_definition=surface_definition,
            nb_sets_from_hints="fallback_to_source",
            other_nb_sets="0_weight",
            additional_condition=DistanceAngleAreaNbSetWeight.AC.ONLY_ACB,
        )

        d1, d2, a1, a2 = 1.05, 1.15, 0.05, 0.08
        assert not da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.05, 1.15, 0.1, 0.2
        assert not da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.9, 1.95, 0.1, 0.2
        assert not da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.05, 1.95, 0.05, 0.25
        assert da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.05, 1.95, 0.75, 0.9
        assert da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.1, 1.9, 0.1, 0.9
        assert da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.23, 1.77, 0.48, 0.52
        assert da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.23, 1.24, 0.48, 0.52
        assert not da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.4, 1.6, 0.4, 0.6
        assert da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.6, 1.9, 0.7, 0.9
        assert not da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.5, 1.6, 0.75, 0.78
        assert not da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.5, 1.6, 0.75, 0.95
        assert not da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.4, 1.6, 0.1, 0.9
        assert da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)
        d1, d2, a1, a2 = 1.4, 1.6, 0.3, 0.7
        assert da_area_weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2)

    def test_dist_nb_set_weight(self):
        dnbset_weight = DistanceNbSetWeight()
        dnbset_weight2 = DistanceNbSetWeight(
            weight_function={
                "function": "smoothstep",
                "options": {"lower": 1.2, "upper": 1.3},
            }
        )

        fake_nb_set1 = FakeNbSet(cn=1)
        fake_nb_set1.site_voronoi_indices = {0}
        fake_nb_set2 = FakeNbSet(cn=2)
        fake_nb_set2.site_voronoi_indices = {0, 1}
        fake_nb_set3 = FakeNbSet(cn=3)
        fake_nb_set3.site_voronoi_indices = {0, 1, 2}
        fake_nb_set4 = FakeNbSet(cn=4)
        fake_nb_set4.site_voronoi_indices = {0, 1, 2, 3}
        fake_nb_set5 = FakeNbSet(cn=5)
        fake_nb_set5.site_voronoi_indices = {0, 1, 2, 3, 4}
        fake_nb_set5_m2 = FakeNbSet(cn=4)
        fake_nb_set5_m2.site_voronoi_indices = {0, 1, 3, 4}
        fake_nb_set6 = FakeNbSet(cn=6)
        fake_nb_set6.site_voronoi_indices = {0, 1, 2, 3, 4, 5}
        fake_nb_set7 = FakeNbSet(cn=7)
        fake_nb_set7.site_voronoi_indices = {0, 1, 2, 3, 4, 5, 6}
        fake_nb_set1.isite = 0
        fake_nb_set2.isite = 0
        fake_nb_set3.isite = 0
        fake_nb_set4.isite = 0
        fake_nb_set5.isite = 0
        fake_nb_set5_m2.isite = 0
        fake_nb_set6.isite = 0
        fake_nb_set7.isite = 0
        dummy_se = DummyStructureEnvironments()
        dummy_se.neighbors_sets = []
        dummy_se.neighbors_sets.append({})
        dummy_se.neighbors_sets[0][1] = [fake_nb_set1]
        dummy_se.neighbors_sets[0][2] = [fake_nb_set2]
        dummy_se.neighbors_sets[0][3] = [fake_nb_set3]
        dummy_se.neighbors_sets[0][4] = [fake_nb_set4, fake_nb_set5_m2]
        dummy_se.neighbors_sets[0][5] = [fake_nb_set5]
        dummy_se.neighbors_sets[0][6] = [fake_nb_set6]
        dummy_se.neighbors_sets[0][7] = [fake_nb_set7]

        dummy_voronoi = DummyVoronoiContainer()
        dummy_voronoi.voronoi_list2 = []
        dummy_voronoi.voronoi_list2.append([])
        dummy_voronoi.voronoi_list2[0].append({"normalized_distance": 1})  # 0
        dummy_voronoi.voronoi_list2[0].append({"normalized_distance": 1.2})  # 1
        dummy_voronoi.voronoi_list2[0].append({"normalized_distance": 1.225})  # 2
        dummy_voronoi.voronoi_list2[0].append({"normalized_distance": 1.25})  # 3
        dummy_voronoi.voronoi_list2[0].append({"normalized_distance": 1.275})  # 4
        dummy_voronoi.voronoi_list2[0].append({"normalized_distance": 1.3})  # 5
        dummy_voronoi.voronoi_list2[0].append({"normalized_distance": 1.8})  # 6
        # Following fake neighbor dict is not in the neighbors sets
        dummy_voronoi.voronoi_list2[0].append({"normalized_distance": 1.55})  # 7

        for fake_nb_set in [
            fake_nb_set1,
            fake_nb_set2,
            fake_nb_set3,
            fake_nb_set4,
            fake_nb_set5,
            fake_nb_set5_m2,
            fake_nb_set6,
            fake_nb_set7,
        ]:
            fake_nb_set.normalized_distances = [
                dummy_voronoi.voronoi_list2[0][ivoro_nb]["normalized_distance"]
                for ivoro_nb in fake_nb_set.site_voronoi_indices
            ]

        dummy_se.voronoi = dummy_voronoi

        cn_map1 = (1, 0)
        cn_map2 = (2, 0)
        cn_map3 = (3, 0)
        cn_map4 = (4, 0)
        cn_map5 = (5, 0)
        cn_map5_m2 = (4, 1)
        cn_map6 = (6, 0)
        cn_map7 = (7, 0)

        weight1 = dnbset_weight.weight(fake_nb_set1, dummy_se, cn_map=cn_map1, additional_info=None)
        assert abs(weight1) < 1e-8
        weight2 = dnbset_weight.weight(fake_nb_set2, dummy_se, cn_map=cn_map2, additional_info=None)
        assert abs(weight2 - 0.103515625) < 1e-8
        weight3 = dnbset_weight.weight(fake_nb_set3, dummy_se, cn_map=cn_map3, additional_info=None)
        assert abs(weight3 - 0.5) < 1e-8
        weight4 = dnbset_weight.weight(fake_nb_set4, dummy_se, cn_map=cn_map4, additional_info=None)
        assert abs(weight4 - 0.896484375) < 1e-8
        weight5 = dnbset_weight.weight(fake_nb_set5, dummy_se, cn_map=cn_map5, additional_info=None)
        assert abs(weight5 - 1) < 1e-8
        weight5_m2 = dnbset_weight.weight(fake_nb_set5_m2, dummy_se, cn_map=cn_map5_m2, additional_info=None)
        assert abs(weight5_m2 - 0.103515625) < 1e-8
        weight7 = dnbset_weight.weight(fake_nb_set7, dummy_se, cn_map=cn_map7, additional_info=None)
        assert abs(weight7 - 1) < 1e-8

        weight_2_3 = dnbset_weight2.weight(fake_nb_set3, dummy_se, cn_map=cn_map3, additional_info=None)
        assert abs(weight_2_3 - 0.5) < 1e-8
        weight_2_4 = dnbset_weight2.weight(fake_nb_set4, dummy_se, cn_map=cn_map4, additional_info=None)
        assert abs(weight_2_4 - 0.84375) < 1e-8
        weight_2_2 = dnbset_weight2.weight(fake_nb_set2, dummy_se, cn_map=cn_map2, additional_info=None)
        assert abs(weight_2_2 - 0.15625) < 1e-8

        dnbset_weight3 = DistanceNbSetWeight(
            weight_function={
                "function": "smoothstep",
                "options": {"lower": 1.5, "upper": 1.7},
            },
            nbs_source="nb_sets",
        )
        dnbset_weight4 = DistanceNbSetWeight(
            weight_function={
                "function": "smoothstep",
                "options": {"lower": 1.5, "upper": 1.7},
            },
            nbs_source="voronoi",
        )

        weight_3_6 = dnbset_weight3.weight(fake_nb_set6, dummy_se, cn_map=cn_map6, additional_info=None)
        assert abs(weight_3_6 - 1) < 1e-8
        weight_4_6 = dnbset_weight4.weight(fake_nb_set6, dummy_se, cn_map=cn_map6, additional_info=None)
        assert abs(weight_4_6 - 0.15625) < 1e-8

        delta_dnb_set_weight = DeltaDistanceNbSetWeight(
            weight_function={
                "function": "smootherstep",
                "options": {"lower": 0.05, "upper": 0.15},
            }
        )

        weight_delta1 = delta_dnb_set_weight.weight(fake_nb_set1, dummy_se, cn_map=cn_map1, additional_info=None)
        assert abs(weight_delta1 - 1) < 1e-8
        weight_delta2 = delta_dnb_set_weight.weight(fake_nb_set2, dummy_se, cn_map=cn_map2, additional_info=None)
        assert abs(weight_delta2) < 1e-8
        weight_delta3 = delta_dnb_set_weight.weight(fake_nb_set3, dummy_se, cn_map=cn_map3, additional_info=None)
        assert abs(weight_delta3) < 1e-8

        delta_dnb_set_weight2 = DeltaDistanceNbSetWeight(
            weight_function={
                "function": "smootherstep",
                "options": {"lower": 0.1, "upper": 0.3},
            }
        )

        weight_delta1 = delta_dnb_set_weight2.weight(fake_nb_set1, dummy_se, cn_map=cn_map1, additional_info=None)
        assert abs(weight_delta1 - 0.5) < 1e-8
        weight_delta2 = delta_dnb_set_weight2.weight(fake_nb_set2, dummy_se, cn_map=cn_map2, additional_info=None)
        assert abs(weight_delta2) < 1e-8
        weight_delta3 = delta_dnb_set_weight2.weight(fake_nb_set3, dummy_se, cn_map=cn_map3, additional_info=None)
        assert abs(weight_delta3) < 1e-8

        delta_dnb_set_weight3 = DeltaDistanceNbSetWeight(
            weight_function={
                "function": "smoothstep",
                "options": {"lower": 0.1, "upper": 0.5},
            }
        )

        weight_delta1 = delta_dnb_set_weight3.weight(fake_nb_set1, dummy_se, cn_map=cn_map1, additional_info=None)
        assert abs(weight_delta1 - 0.15625) < 1e-8
        weight_delta6 = delta_dnb_set_weight3.weight(fake_nb_set6, dummy_se, cn_map=cn_map6, additional_info=None)
        assert abs(weight_delta6 - 0.31640625) < 1e-8

        delta_dnb_set_weight4 = DeltaDistanceNbSetWeight(
            weight_function={
                "function": "smoothstep",
                "options": {"lower": 0.1, "upper": 0.5},
            },
            nbs_source="nb_sets",
        )

        weight_delta1 = delta_dnb_set_weight4.weight(fake_nb_set1, dummy_se, cn_map=cn_map1, additional_info=None)
        assert abs(weight_delta1 - 0.15625) < 1e-8
        weight_delta6 = delta_dnb_set_weight4.weight(fake_nb_set6, dummy_se, cn_map=cn_map6, additional_info=None)
        assert abs(weight_delta6 - 1) < 1e-8


@pytest.fixture(scope="module")
def se_sio2():
    """StructureEnvironments of SiO2 (mp-7000); sites 6, 7 and 8 are tetrahedral Si."""
    with open(f"{struct_env_dir}/se_mp-7000.json", "rb") as file:
        return StructureEnvironments.from_dict(orjson.loads(file.read()))


def _nb_set(se, cn_map, isite=6):
    return se.neighbors_sets[isite][cn_map[0]][cn_map[1]]


class TestWeightSerialization:
    @pytest.mark.parametrize(
        "weight",
        [
            AngleNbSetWeight(aa=2),
            NormalizedAngleDistanceNbSetWeight(average_type="arithmetic", aa=2, bb=0.5),
            SelfCSMNbSetWeight(
                weight_estimator={"function": "power2_decreasing_exp", "options": {"max_csm": 5.0, "alpha": 2}},
                symmetry_measure_type="csm_wocs_ctwocc",
            ),
            DeltaCSMNbSetWeight.delta_cn_specifics(),
            DeltaCSMNbSetWeight(weight_estimator=None),
            CNBiasNbSetWeight.linearly_equidistant(weight_cn1=1, weight_cn13=13),
            DistanceAngleAreaNbSetWeight(additional_condition=0),
            DistancePlateauNbSetWeight(),
            AnglePlateauNbSetWeight(),
            DeltaDistanceNbSetWeight(nbs_source="nb_sets"),
        ],
        ids=lambda weight: type(weight).__name__,
    )
    def test_json_round_trip(self, weight):
        # Go through an actual JSON string so that e.g. int keys become str keys
        dct = json.loads(json.dumps(weight.as_dict()))
        assert dct["@class"] == type(weight).__name__
        new_weight = type(weight).from_dict(dct)
        assert new_weight == weight
        assert new_weight.as_dict() == weight.as_dict()

    def test_distance_nb_set_weight_round_trip(self):
        weight = DistanceNbSetWeight(
            weight_function={"function": "smoothstep", "options": {"lower": 1.1, "upper": 1.5}},
            nbs_source="nb_sets",
        )
        new_weight = DistanceNbSetWeight.from_dict(weight.as_dict())
        assert new_weight.weight_function == weight.weight_function
        assert new_weight.nbs_source == "nb_sets"
        assert new_weight.weight_rf.eval(1.3) == approx(0.5)

    @pytest.mark.parametrize(
        "weight",
        [
            AngleNbSetWeight(),
            NormalizedAngleDistanceNbSetWeight(average_type="geometric", aa=1, bb=1),
            SelfCSMNbSetWeight(),
            DeltaCSMNbSetWeight(),
            DistanceAngleAreaNbSetWeight(),
        ],
        ids=lambda weight: type(weight).__name__,
    )
    def test_eq_with_foreign_object(self, weight):
        assert weight.__eq__(object()) is NotImplemented
        assert weight != "not a weight"

    def test_cn_bias_eq(self):
        lin = CNBiasNbSetWeight.linearly_equidistant(weight_cn1=1, weight_cn13=13)
        assert lin.__eq__(object()) is NotImplemented
        assert lin == CNBiasNbSetWeight.from_description(
            {"type": "linearly_equidistant", "weight_cn1": 1, "weight_cn13": 13}
        )
        assert lin != CNBiasNbSetWeight.linearly_equidistant(weight_cn1=1, weight_cn13=14)

    def test_plateau_and_delta_distance_eq(self):
        # These weights compare equal whenever they are of the same type
        assert DistancePlateauNbSetWeight() == DistancePlateauNbSetWeight(distance_function={"type": "other"})
        assert AnglePlateauNbSetWeight() != DistancePlateauNbSetWeight()
        assert DeltaDistanceNbSetWeight() == DeltaDistanceNbSetWeight(nbs_source="nb_sets")


class TestWeightValidation:
    def test_normalized_angle_distance_both_exponents_zero(self):
        with pytest.raises(ValueError, match="Both exponents are 0"):
            NormalizedAngleDistanceNbSetWeight(average_type="geometric", aa=0, bb=0)

    @pytest.mark.parametrize(
        ("kwargs", "error", "match"),
        [
            ({"weight_type": "has_intersection_smoothstep"}, NotImplementedError, ""),
            ({"weight_type": "area"}, ValueError, "Weight type is 'area'"),
            ({"other_nb_sets": "1_weight"}, ValueError, 'Other nb_sets should be "0_weight"'),
            ({"nb_sets_from_hints": "ignore"}, ValueError, "Nb_sets from hints should fallback to source"),
        ],
    )
    def test_dist_angle_area_invalid_options(self, kwargs, error, match):
        with pytest.raises(error, match=match):
            DistanceAngleAreaNbSetWeight(**kwargs)

    @pytest.mark.parametrize("cls", [DistanceNbSetWeight, DeltaDistanceNbSetWeight])
    def test_invalid_nbs_source(self, cls):
        with pytest.raises(ValueError, match='"nbs_source" should be one of'):
            cls(nbs_source="neighbors")

    def test_cn_bias_from_description_unknown_type(self):
        with pytest.raises(RuntimeError, match="Cannot initialize Weights"):
            CNBiasNbSetWeight.from_description({"type": "random"})

    def test_cn_bias_geometrically_equidistant(self):
        weight = CNBiasNbSetWeight.geometrically_equidistant(weight_cn1=1, weight_cn13=4096)
        assert weight.initialization_options == {
            "type": "geometrically_equidistant",
            "weight_cn1": 1,
            "weight_cn13": 4096,
        }
        # 4096 = 2**12, hence the weight doubles with each additional neighbor
        for cn in range(1, 14):
            assert weight.cn_weights[cn] == approx(2 ** (cn - 1))
        assert weight.weight(FakeNbSet(cn=5), structure_environments=None) == approx(16)


class TestNormalizedAngleDistanceFunctions:
    normalized_angles = (1.0, 0.5, 0.25)
    normalized_distances = (1.0, 1.25, 2.0)

    @pytest.mark.parametrize(
        ("aa", "bb", "fda_name"),
        [
            (0, 1, "invdist"),
            (0, 2, "invndist"),
            (1, 0, "ang"),
            (2, 0, "angn"),
            (1, 1, "anginvdist"),
            (1, 2, "anginvndist"),
            (2, 1, "angninvdist"),
            (2, 2, "angninvndist"),
        ],
    )
    @pytest.mark.parametrize("average_type", ["arithmetic", "geometric"])
    def test_fda_selection_and_weight(self, aa, bb, fda_name, average_type):
        weight = NormalizedAngleDistanceNbSetWeight(average_type=average_type, aa=aa, bb=bb)
        assert weight.fda.__name__ == fda_name

        nb_set = FakeNbSet(cn=3)
        nb_set.normalized_angles = list(self.normalized_angles)
        nb_set.normalized_distances = list(self.normalized_distances)
        per_neighbor = [
            ang**aa / dist**bb for ang, dist in zip(self.normalized_angles, self.normalized_distances, strict=True)
        ]
        if average_type == "arithmetic":
            expected = sum(per_neighbor) / 3
        else:
            expected = (per_neighbor[0] * per_neighbor[1] * per_neighbor[2]) ** (1 / 3)
        assert weight.weight(nb_set=nb_set, structure_environments=None) == approx(expected)


class TestDistanceAngleAreaRectangles:
    # Default surface: distance in [1.2, 1.8], angle in [0.1, 0.8]
    weight = DistanceAngleAreaNbSetWeight()

    @pytest.mark.parametrize(
        ("d1", "d2", "a1", "a2", "expected"),
        [
            (1.0, 1.1, 0.0, 1.0, False),  # case 1: rectangle left of the area
            (1.9, 2.0, 0.0, 1.0, False),  # case 6: rectangle right of the area
            (1.0, 1.5, 0.3, 0.6, True),  # case 2: overlaps the lower-left part
            (1.0, 1.5, 0.85, 0.9, False),  # case 2: above amax
            (1.0, 2.0, 0.3, 0.6, True),  # case 3: spans the whole distance range
            (1.0, 2.0, 0.0, 0.05, False),  # case 3: below amin
            (1.4, 1.6, 0.3, 0.5, True),  # case 4: inside the area
            (1.4, 1.6, 0.0, 0.01, False),  # case 4: below the lower curve
            (1.4, 1.6, 0.79, 0.8, False),  # case 4: above the upper curve
            (1.5, 2.0, 0.2, 0.4, True),  # case 5: overlaps the right part
            (1.5, 2.0, 0.0, 0.05, False),  # case 5: below amin
        ],
    )
    def test_rectangle_crosses_area(self, d1, d2, a1, a2, expected):
        assert self.weight.rectangle_crosses_area(d1=d1, d2=d2, a1=a1, a2=a2) is expected

    def test_rectangle_inconsistent(self):
        # d1 inside the area but d2 below dmin is not a valid rectangle
        with pytest.raises(ValueError, match="Should not reach this point"):
            self.weight.rectangle_crosses_area(d1=1.5, d2=1.0, a1=0.2, a2=0.4)


def _dist_ang_source(d_min, d_next, a_next, a_max, ac=1):
    return {
        "origin": "dist_ang_ac_voronoi",
        "ac": ac,
        "dp_dict": {"min": d_min, "next": d_next},
        "ap_dict": {"next": a_next, "max": a_max},
    }


class TestDistanceAngleAreaFromHints:
    weight = DistanceAngleAreaNbSetWeight()

    @pytest.fixture
    def dummy_se(self):
        src_in = FakeNbSet(cn=4)
        src_in.sources = [_dist_ang_source(1.0, 1.5, 0.3, 0.6)]
        src_out = FakeNbSet(cn=4)
        src_out.sources = [_dist_ang_source(1.9, 2.0, 0.3, 0.6)]
        src_other_ac = FakeNbSet(cn=4)
        src_other_ac.sources = [_dist_ang_source(1.0, 1.5, 0.3, 0.6, ac=0)]
        dummy_se = DummyStructureEnvironments()
        dummy_se.neighbors_sets = [{4: [src_in, src_out, src_other_ac]}]
        return dummy_se

    @staticmethod
    def _hinted_nb_set(*cn_map_sources):
        nb_set = FakeNbSet(cn=5)
        nb_set.isite = 0
        nb_set.sources = [{"origin": "nb_set_hints", "cn_map_source": cn_map} for cn_map in cn_map_sources]
        return nb_set

    @pytest.mark.parametrize(("cn_map_source", "expected"), [((4, 0), 1), ((4, 1), 0), ((4, 2), 0)])
    def test_fallback_to_source(self, dummy_se, cn_map_source, expected):
        nb_set = self._hinted_nb_set(cn_map_source)
        assert self.weight.weight(nb_set, dummy_se, cn_map=(5, 0)) == expected

    def test_no_source(self, dummy_se):
        nb_set = self._hinted_nb_set()
        assert self.weight.weight(nb_set, dummy_se, cn_map=(5, 0)) == 0

    def test_multiple_hints(self, dummy_se):
        nb_set = self._hinted_nb_set((4, 0), (4, 1))
        with pytest.raises(ValueError, match="Found multiple hints sources for nb_set"):
            self.weight.weight(nb_set, dummy_se, cn_map=(5, 0))


class TestDeltaDistanceNoCompetitors:
    @pytest.mark.parametrize("nbs_source", ["voronoi", "nb_sets"])
    @pytest.mark.parametrize("cn", [0, 2])
    def test_no_other_neighbors(self, nbs_source, cn):
        # All the voronoi neighbors belong to the nb_set: no competitor -> weight of 1,
        # even for an inverse weight function (for which f(+inf) would be 0)
        nb_set = FakeNbSet(cn=cn)
        nb_set.isite = 0
        nb_set.site_voronoi_indices = {0, 1}
        nb_set.normalized_distances = [1.0, 1.1][:cn]
        dummy_voronoi = DummyVoronoiContainer()
        dummy_voronoi.voronoi_list2 = [[{"normalized_distance": 1.0}, {"normalized_distance": 1.1}]]
        dummy_se = DummyStructureEnvironments()
        dummy_se.voronoi = dummy_voronoi
        dummy_se.neighbors_sets = [{cn: [nb_set]}]
        weight = DeltaDistanceNbSetWeight(
            weight_function={"function": "inverse_smootherstep", "options": {"lower": 0.1, "upper": 0.2}},
            nbs_source=nbs_source,
        )
        assert weight.weight(nb_set, dummy_se, cn_map=(cn, 0)) == 1


class TestWeightsOnRealEnvironments:
    """Weights of the neighbors sets of the tetrahedral Si site 6 in SiO2 (mp-7000)."""

    @pytest.mark.parametrize(
        ("cn_map", "expected"),
        [
            # plateau = 0.0018 < lower = 0.2 -> inverse_smootherstep gives 1
            ((2, 0), 1.0),
            # plateau = 1.22 > upper = 0.4 -> inverse_smootherstep gives 0
            ((2, 1), 0.0),
            ((4, 0), 0.0),
        ],
    )
    def test_distance_plateau(self, se_sio2, cn_map, expected):
        nb_set = _nb_set(se_sio2, cn_map)
        weight = DistancePlateauNbSetWeight()
        assert weight.weight(nb_set, se_sio2, cn_map=cn_map) == approx(expected)

    @pytest.mark.parametrize(("cn_map", "expected"), [((2, 0), 0.0), ((2, 1), 1.0), ((6, 0), 1.0)])
    def test_angle_plateau(self, se_sio2, cn_map, expected):
        nb_set = _nb_set(se_sio2, cn_map)
        weight = AnglePlateauNbSetWeight()
        assert weight.weight(nb_set, se_sio2, cn_map=cn_map) == approx(expected)

    @pytest.mark.parametrize("nbs_source", ["voronoi", "nb_sets"])
    def test_distance_nb_set_weight(self, se_sio2, nbs_source):
        # The 4 O neighbors of Si are at normalized distances 1.0 and 1.0018, the next O at 2.2238.
        weight = DistanceNbSetWeight(
            weight_function={"function": "smootherstep", "options": {"lower": 1.2, "upper": 1.3}},
            nbs_source=nbs_source,
        )
        assert weight.weight(_nb_set(se_sio2, (4, 0)), se_sio2, cn_map=(4, 0)) == approx(1.0)
        # For the 2-neighbors set, the closest neighbor not in the set is at 1.0018 < 1.2
        assert weight.weight(_nb_set(se_sio2, (2, 0)), se_sio2, cn_map=(2, 0)) == approx(0.0)

    def test_distance_nb_set_weight_no_candidate(self, se_sio2):
        # The 6-neighbors set contains all the neighbors of the other neighbors sets
        weight = DistanceNbSetWeight(nbs_source="nb_sets")
        assert weight.weight(_nb_set(se_sio2, (6, 0)), se_sio2, cn_map=(6, 0)) == 1

    @pytest.mark.parametrize("nbs_source", ["voronoi", "nb_sets"])
    def test_delta_distance_weight(self, se_sio2, nbs_source):
        weight = DeltaDistanceNbSetWeight(nbs_source=nbs_source)  # smootherstep between 0.1 and 0.2
        # Gap between the 4 O first neighbors (1.0018) and the next ones (2.2238)
        assert weight.weight(_nb_set(se_sio2, (4, 0)), se_sio2, cn_map=(4, 0)) == approx(1)
        # No gap between the 2-neighbors set (1.0) and the other first neighbors (1.0018)
        assert weight.weight(_nb_set(se_sio2, (2, 0)), se_sio2, cn_map=(2, 0)) == approx(0)
        # Empty neighbors set
        assert weight.weight(_nb_set(se_sio2, (0, 0)), se_sio2, cn_map=(0, 0)) == 0

    def test_delta_distance_weight_6_neighbors(self, se_sio2):
        nb_set = _nb_set(se_sio2, (6, 0))
        voronoi = se_sio2.voronoi.voronoi_list2[6]
        delta = min(
            nb["normalized_distance"] for inb, nb in enumerate(voronoi) if inb not in nb_set.site_voronoi_indices
        ) - max(nb_set.normalized_distances)
        assert delta == approx(0.1046483, abs=1e-6)
        scaled = (delta - 0.1) / 0.1
        expected = 6 * scaled**5 - 15 * scaled**4 + 10 * scaled**3
        assert DeltaDistanceNbSetWeight().weight(nb_set, se_sio2, cn_map=(6, 0)) == approx(expected)
        # No other neighbors set contains neighbors outside of the 6-neighbors set
        assert DeltaDistanceNbSetWeight(nbs_source="nb_sets").weight(nb_set, se_sio2, cn_map=(6, 0)) == 1

    def test_self_csm_weight(self, se_sio2):
        weight = SelfCSMNbSetWeight()
        additional_info = {}
        weight_t4 = weight.weight(_nb_set(se_sio2, (4, 0)), se_sio2, cn_map=(4, 0), additional_info=additional_info)
        assert weight_t4 == approx(0.996297423523595)
        # Empty neighbors set: no coordination environment -> effective CSM of 100 and zero weight
        assert weight.weight(_nb_set(se_sio2, (0, 0)), se_sio2, cn_map=(0, 0), additional_info=additional_info) == 0
        assert additional_info["effective_csms"][6][(0, 0)] == 100
        assert set(additional_info["self_csms_weights"][6]) == {(4, 0), (0, 0)}
        # Cached effective csm is reused
        additional_info["effective_csms"][6][(4, 0)] = 100
        assert weight.weight(_nb_set(se_sio2, (4, 0)), se_sio2, cn_map=(4, 0), additional_info=additional_info) == 0

    def test_delta_csm_weight(self, se_sio2):
        weight = DeltaCSMNbSetWeight()
        additional_info = {}
        # T:4 is much more symmetric than the 6-neighbors set -> no penalty
        assert weight.weight(_nb_set(se_sio2, (4, 0)), se_sio2, cn_map=(4, 0), additional_info=additional_info) == 1
        assert additional_info["delta_csms_weights"][6][(4, 0)] == 1
        assert additional_info["delta_csms_cn_map2"][6][(4, 0)] is None
        # The 2-neighbors set is penalized by the much more symmetric T:4 environment
        assert weight.weight(_nb_set(se_sio2, (2, 0)), se_sio2, cn_map=(2, 0), additional_info=additional_info) == 0
        assert additional_info["delta_csms_weights"][6][(2, 0)] == 0
        assert additional_info["delta_csms_cn_map2"][6][(2, 0)][0] > 2
        assert additional_info["delta_csms"][6][(2, 0)] < 0

    def test_delta_csm_weight_without_estimators(self, se_sio2):
        weight = DeltaCSMNbSetWeight(weight_estimator=None)
        assert weight.weight(_nb_set(se_sio2, (2, 0)), se_sio2, cn_map=(2, 0), additional_info={}) == 1

    def test_delta_cn_specifics(self):
        weight = DeltaCSMNbSetWeight.delta_cn_specifics(
            delta_csm_mins=[float(dcn) for dcn in range(1, 13)],
            delta_csm_maxs=[float(dcn + 2) for dcn in range(1, 13)],
        )
        assert set(weight.delta_cn_weight_estimators) == set(range(1, 13))
        assert weight.delta_cn_weight_estimators[3]["options"] == {"delta_csm_min": 3.0, "delta_csm_max": 5.0}
        assert weight.weight_estimator["options"] == {"delta_csm_min": 12.0, "delta_csm_max": 14.0}
        default = DeltaCSMNbSetWeight.delta_cn_specifics()
        assert default.delta_cn_weight_estimators[3]["options"] == {"delta_csm_min": 1.0, "delta_csm_max": 5.75}
        assert default != weight

    def test_dist_angle_area_weight(self, se_sio2):
        weight = DistanceAngleAreaNbSetWeight()
        assert weight.weight(_nb_set(se_sio2, (4, 0)), se_sio2, cn_map=(4, 0)) == 1
        assert weight.weight(_nb_set(se_sio2, (2, 0)), se_sio2, cn_map=(2, 0)) == 0
        assert weight.weight(_nb_set(se_sio2, (6, 0)), se_sio2, cn_map=(6, 0)) == 0
