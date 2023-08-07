from __future__ import annotations

import pytest
from pytest import approx

from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
    AdditionalConditionInt,
    AngleCutoffFloat,
    CSMFloat,
    DistanceCutoffFloat,
    SimplestChemenvStrategy,
)
from pymatgen.util.testing import PymatgenTest

__author__ = "waroquiers"


class TestStrategyOptions(PymatgenTest):
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
            ValueError, match="Continuous symmetry measure limits should be between 0 and 100, got 100.1"
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
            ValueError, match="Continuous symmetry measure limits should be between 0 and 100, got -0.1"
        ):
            simplest_strategy.setup_options({"continuous_symmetry_measure_cutoff": -0.1})

        with pytest.raises(
            ValueError, match="Continuous symmetry measure limits should be between 0 and 100, got 100.1"
        ):
            simplest_strategy.setup_options({"continuous_symmetry_measure_cutoff": 100.1})
