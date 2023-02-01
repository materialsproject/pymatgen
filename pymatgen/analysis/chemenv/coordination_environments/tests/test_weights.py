from __future__ import annotations

import json
import os
import unittest

import pytest

from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
    AngleNbSetWeight,
    CNBiasNbSetWeight,
    DeltaCSMNbSetWeight,
    DeltaDistanceNbSetWeight,
    DistanceAngleAreaNbSetWeight,
    DistanceNbSetWeight,
    NormalizedAngleDistanceNbSetWeight,
    SelfCSMNbSetWeight,
)
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import (
    StructureEnvironments,
)
from pymatgen.util.testing import PymatgenTest

__author__ = "waroquiers"

se_files_dir = os.path.join(
    PymatgenTest.TEST_FILES_DIR,
    "chemenv",
    "structure_environments_files",
)


class FakeNbSet:
    def __init__(self, cn=None):
        self.cn = cn

    def __len__(self):
        return self.cn


class DummyStructureEnvironments:
    pass


class DummyVoronoiContainer:
    pass


class StrategyWeightsTest(PymatgenTest):
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
        angle_weight = AngleNbSetWeight(aa=1.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.9634354419021528) < 1e-8
        angle_weight = AngleNbSetWeight(aa=2.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.92820785071319645) < 1e-8
        angle_weight = AngleNbSetWeight(aa=0.5)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.98154747307613843) < 1e-8

        assert AngleNbSetWeight(1.0) != AngleNbSetWeight(2.0)

        # nb_set with no neighbor
        fake_nb_set.angles = []
        angle_weight = AngleNbSetWeight(aa=1.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.0) < 1e-8
        angle_weight = AngleNbSetWeight(aa=2.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.0) < 1e-8

        # nb_set with one neighbor
        fake_nb_set.angles = [3.08570351705799]
        angle_weight = AngleNbSetWeight(aa=1.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 0.24555248382791284) < 1e-8
        angle_weight = AngleNbSetWeight(aa=2.0)
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
        angle_weight = AngleNbSetWeight(aa=1.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 1.0) < 1e-8
        angle_weight = AngleNbSetWeight(aa=2.0)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 1.0) < 1e-8
        angle_weight = AngleNbSetWeight(aa=0.5)
        aw = angle_weight.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(aw - 1.0) < 1e-8

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
            ValueError, match="Average type is 'arithmetix' while it should be 'geometric' or 'arithmetic'"
        ):
            NormalizedAngleDistanceNbSetWeight(average_type="arithmetix", aa=1, bb=1)

        fake_nb_set.normalized_distances = [
            1.2632574171572457,
            1.1231971151388764,
            1.0,
            1.1887986376446249,
            1.188805134890625,
        ]
        fake_nb_set.normalized_angles = [
            0.6026448601336767,
            0.8498933334305273,
            1.0,
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
        bias_weight1 = CNBiasNbSetWeight.linearly_equidistant(weight_cn1=1.0, weight_cn13=13.0)
        bias_weight2 = CNBiasNbSetWeight.geometrically_equidistant(weight_cn1=1.0, weight_cn13=1.1**12)
        bias_weight3 = CNBiasNbSetWeight.explicit(
            cn_weights={
                1: 1.0,
                2: 3.0,
                3: 3.2,
                4: 4.0,
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
            CNBiasNbSetWeight.explicit(cn_weights={1: 1.0, 13: 2.0})

        fake_nb_set.cn = 1
        w1 = bias_weight1.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w1 - 1.0) < 1e-8
        w2 = bias_weight2.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w2 - 1.0) < 1e-8
        w3 = bias_weight3.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w3 - 1.0) < 1e-8
        fake_nb_set.cn = 7
        w1 = bias_weight1.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w1 - 7.0) < 1e-8
        w2 = bias_weight2.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w2 - 1.1**6) < 1e-8
        w3 = bias_weight3.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w3 - 4.3) < 1e-8
        fake_nb_set.cn = 13
        w1 = bias_weight1.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w1 - 13.0) < 1e-8
        w2 = bias_weight2.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w2 - 1.1**12) < 1e-8
        w3 = bias_weight3.weight(nb_set=fake_nb_set, structure_environments=dummy_se)
        assert abs(w3 - 4.8) < 1e-8

        bias_weight4 = CNBiasNbSetWeight.from_description(
            {"type": "linearly_equidistant", "weight_cn1": 2.0, "weight_cn13": 26.0}
        )
        for cn in range(1, 14):
            assert round(abs(bias_weight4.cn_weights[cn] - 2.0 * cn), 7) == 0

        bias_weight5 = CNBiasNbSetWeight.from_description(
            {
                "type": "geometrically_equidistant",
                "weight_cn1": 1.0,
                "weight_cn13": 13.0,
            }
        )
        assert round(abs(bias_weight5.cn_weights[1] - 1.0), 7) == 0
        assert round(abs(bias_weight5.cn_weights[3] - 1.5334062370163877), 7) == 0
        assert round(abs(bias_weight5.cn_weights[9] - 5.5287748136788739), 7) == 0
        assert round(abs(bias_weight5.cn_weights[12] - 10.498197520079623), 7) == 0

        cn_weights = {cn: 0.0 for cn in range(1, 14)}
        cn_weights[6] = 2.0
        cn_weights[4] = 1.0
        bias_weight6 = CNBiasNbSetWeight.from_description({"type": "explicit", "cn_weights": cn_weights})

        assert round(abs(bias_weight6.cn_weights[1] - 0.0), 7) == 0
        assert round(abs(bias_weight6.cn_weights[4] - 1.0), 7) == 0
        assert round(abs(bias_weight6.cn_weights[6] - 2.0), 7) == 0

    def test_self_csms_weight(self):
        # Get the StructureEnvironments for K2NaNb2Fe7Si8H4O31 (mp-743972)
        with open(os.path.join(se_files_dir, "se_mp-743972.json")) as f:
            dd = json.load(f)
        se = StructureEnvironments.from_dict(dd)

        # Get neighbors sets for which we get the weights
        cn_maps = [(12, 3), (12, 2), (13, 2), (12, 0), (12, 1)]
        nbsets = {cn_map: se.neighbors_sets[0][cn_map[0]][cn_map[1]] for cn_map in cn_maps}

        effective_csm_estimator = {
            "function": "power2_inverse_decreasing",
            "options": {"max_csm": 8.0},
        }
        weight_estimator = {
            "function": "power2_decreasing_exp",
            "options": {"max_csm": 8.0, "alpha": 1.0},
        }
        weight_estimator2 = {
            "function": "power2_decreasing_exp",
            "options": {"max_csm": 8.1, "alpha": 1.0},
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
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.11671945916431022) < 1e-8
        cn_map = (12, 2)
        self_w = self_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.0) < 1e-8
        cn_map = (12, 0)
        self_w = self_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.0) < 1e-8
        cn_map = (12, 1)
        self_w = self_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.0) < 1e-8
        cn_map = (13, 2)
        self_w = self_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.14204073172729198) < 1e-8

        # Get the StructureEnvironments for SiO2 (mp-7000)
        with open(os.path.join(se_files_dir, "se_mp-7000.json")) as f:
            dd = json.load(f)
        se = StructureEnvironments.from_dict(dd)

        # Get neighbors sets for which we get the weights
        cn_maps = [(2, 0), (4, 0)]
        nbsets = {cn_map: se.neighbors_sets[6][cn_map[0]][cn_map[1]] for cn_map in cn_maps}

        effective_csm_estimator = {
            "function": "power2_inverse_decreasing",
            "options": {"max_csm": 8.0},
        }

        weight_estimator = {
            "function": "power2_decreasing_exp",
            "options": {"max_csm": 8.0, "alpha": 1.0},
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
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.8143992162836029) < 1e-8
        cn_map = (4, 0)
        self_w = self_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(self_w - 0.99629742352359496) < 1e-8

    def test_delta_csms_weight(self):
        # Get the StructureEnvironments for K2NaNb2Fe7Si8H4O31 (mp-743972)
        with open(os.path.join(se_files_dir, "se_mp-743972.json")) as f:
            dd = json.load(f)
        se = StructureEnvironments.from_dict(dd)

        # Get neighbors sets for which we get the weights
        cn_maps = [(12, 3), (12, 2), (13, 2), (12, 0), (12, 1), (13, 0), (13, 1)]
        nbsets = {cn_map: se.neighbors_sets[0][cn_map[0]][cn_map[1]] for cn_map in cn_maps}

        effective_csm_estimator = {
            "function": "power2_inverse_decreasing",
            "options": {"max_csm": 8.0},
        }
        weight_estimator = {
            "function": "smootherstep",
            "options": {"delta_csm_min": 0.5, "delta_csm_max": 3.0},
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
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.0) < 1e-8
        cn_map = (12, 2)
        delta_w = delta_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.0) < 1e-8
        cn_map = (12, 0)
        delta_w = delta_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.0) < 1e-8
        cn_map = (12, 1)
        delta_w = delta_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.0) < 1e-8
        cn_map = (13, 2)
        delta_w = delta_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 1.0) < 1e-8
        cn_map = (13, 0)
        delta_w = delta_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.0) < 1e-8
        cn_map = (13, 1)
        delta_w = delta_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.0) < 1e-8

        effective_csm_estimator = {
            "function": "power2_inverse_decreasing",
            "options": {"max_csm": 8.0},
        }

        weight_estimator = {
            "function": "smootherstep",
            "options": {"delta_csm_min": -1.0, "delta_csm_max": 3.0},
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
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.040830741048481355) < 1e-8
        cn_map = (13, 2)
        delta_w = delta_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 1.0) < 1e-8
        cn_map = (13, 0)
        delta_w = delta_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.103515625) < 1e-8
        cn_map = (13, 1)
        delta_w = delta_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.103515625) < 1e-8

        # Get the StructureEnvironments for SiO2 (mp-7000)
        with open(os.path.join(se_files_dir, "se_mp-7000.json")) as f:
            dd = json.load(f)
        se = StructureEnvironments.from_dict(dd)

        # Get neighbors sets for which we get the weights
        cn_maps = [(2, 0), (4, 0)]
        nbsets = {cn_map: se.neighbors_sets[6][cn_map[0]][cn_map[1]] for cn_map in cn_maps}

        effective_csm_estimator = {
            "function": "power2_inverse_decreasing",
            "options": {"max_csm": 8.0},
        }

        weight_estimator = {
            "function": "smootherstep",
            "options": {"delta_csm_min": 0.5, "delta_csm_max": 3.0},
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
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 0.0) < 1e-8
        cn_map = (4, 0)
        delta_w = delta_weight.weight(
            nb_set=nbsets[cn_map],
            structure_environments=se,
            cn_map=cn_map,
            additional_info=additional_info,
        )
        assert abs(delta_w - 1.0) < 1e-8

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
        dummy_voronoi.voronoi_list2[0].append({"normalized_distance": 1.0})  # 0
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

        myweight1 = dnbset_weight.weight(fake_nb_set1, dummy_se, cn_map=cn_map1, additional_info=None)
        assert abs(myweight1 - 0.0) < 1e-8
        myweight2 = dnbset_weight.weight(fake_nb_set2, dummy_se, cn_map=cn_map2, additional_info=None)
        assert abs(myweight2 - 0.103515625) < 1e-8
        myweight3 = dnbset_weight.weight(fake_nb_set3, dummy_se, cn_map=cn_map3, additional_info=None)
        assert abs(myweight3 - 0.5) < 1e-8
        myweight4 = dnbset_weight.weight(fake_nb_set4, dummy_se, cn_map=cn_map4, additional_info=None)
        assert abs(myweight4 - 0.896484375) < 1e-8
        myweight5 = dnbset_weight.weight(fake_nb_set5, dummy_se, cn_map=cn_map5, additional_info=None)
        assert abs(myweight5 - 1.0) < 1e-8
        myweight5_m2 = dnbset_weight.weight(fake_nb_set5_m2, dummy_se, cn_map=cn_map5_m2, additional_info=None)
        assert abs(myweight5_m2 - 0.103515625) < 1e-8
        myweight7 = dnbset_weight.weight(fake_nb_set7, dummy_se, cn_map=cn_map7, additional_info=None)
        assert abs(myweight7 - 1.0) < 1e-8

        myweight_2_3 = dnbset_weight2.weight(fake_nb_set3, dummy_se, cn_map=cn_map3, additional_info=None)
        assert abs(myweight_2_3 - 0.5) < 1e-8
        myweight_2_4 = dnbset_weight2.weight(fake_nb_set4, dummy_se, cn_map=cn_map4, additional_info=None)
        assert abs(myweight_2_4 - 0.84375) < 1e-8
        myweight_2_2 = dnbset_weight2.weight(fake_nb_set2, dummy_se, cn_map=cn_map2, additional_info=None)
        assert abs(myweight_2_2 - 0.15625) < 1e-8

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

        myweight_3_6 = dnbset_weight3.weight(fake_nb_set6, dummy_se, cn_map=cn_map6, additional_info=None)
        assert abs(myweight_3_6 - 1.0) < 1e-8
        myweight_4_6 = dnbset_weight4.weight(fake_nb_set6, dummy_se, cn_map=cn_map6, additional_info=None)
        assert abs(myweight_4_6 - 0.15625) < 1e-8

        deltadnbset_weight = DeltaDistanceNbSetWeight(
            weight_function={
                "function": "smootherstep",
                "options": {"lower": 0.05, "upper": 0.15},
            }
        )

        myweightdelta1 = deltadnbset_weight.weight(fake_nb_set1, dummy_se, cn_map=cn_map1, additional_info=None)
        assert abs(myweightdelta1 - 1.0) < 1e-8
        myweightdelta2 = deltadnbset_weight.weight(fake_nb_set2, dummy_se, cn_map=cn_map2, additional_info=None)
        assert abs(myweightdelta2 - 0.0) < 1e-8
        myweightdelta3 = deltadnbset_weight.weight(fake_nb_set3, dummy_se, cn_map=cn_map3, additional_info=None)
        assert abs(myweightdelta3 - 0.0) < 1e-8

        deltadnbset_weight2 = DeltaDistanceNbSetWeight(
            weight_function={
                "function": "smootherstep",
                "options": {"lower": 0.1, "upper": 0.3},
            }
        )

        myweightdelta1 = deltadnbset_weight2.weight(fake_nb_set1, dummy_se, cn_map=cn_map1, additional_info=None)
        assert abs(myweightdelta1 - 0.5) < 1e-8
        myweightdelta2 = deltadnbset_weight2.weight(fake_nb_set2, dummy_se, cn_map=cn_map2, additional_info=None)
        assert abs(myweightdelta2 - 0.0) < 1e-8
        myweightdelta3 = deltadnbset_weight2.weight(fake_nb_set3, dummy_se, cn_map=cn_map3, additional_info=None)
        assert abs(myweightdelta3 - 0.0) < 1e-8

        deltadnbset_weight3 = DeltaDistanceNbSetWeight(
            weight_function={
                "function": "smoothstep",
                "options": {"lower": 0.1, "upper": 0.5},
            }
        )

        myweightdelta1 = deltadnbset_weight3.weight(fake_nb_set1, dummy_se, cn_map=cn_map1, additional_info=None)
        assert abs(myweightdelta1 - 0.15625) < 1e-8
        myweightdelta6 = deltadnbset_weight3.weight(fake_nb_set6, dummy_se, cn_map=cn_map6, additional_info=None)
        assert abs(myweightdelta6 - 0.31640625) < 1e-8

        deltadnbset_weight4 = DeltaDistanceNbSetWeight(
            weight_function={
                "function": "smoothstep",
                "options": {"lower": 0.1, "upper": 0.5},
            },
            nbs_source="nb_sets",
        )

        myweightdelta1 = deltadnbset_weight4.weight(fake_nb_set1, dummy_se, cn_map=cn_map1, additional_info=None)
        assert abs(myweightdelta1 - 0.15625) < 1e-8
        myweightdelta6 = deltadnbset_weight4.weight(fake_nb_set6, dummy_se, cn_map=cn_map6, additional_info=None)
        assert abs(myweightdelta6 - 1.0) < 1e-8


if __name__ == "__main__":
    unittest.main()
