from __future__ import annotations

import unittest

import numpy as np
import pytest

from pymatgen.analysis.chemenv.utils.func_utils import (
    CSMFiniteRatioFunction,
    CSMInfiniteRatioFunction,
    DeltaCSMRatioFunction,
)

__author__ = "waroquiers"


class FuncUtilsTest(unittest.TestCase):
    def test_CSMFiniteRatioFunction(self):
        max_csm = 8.0
        alpha = 1.0
        csm_finite_ratio = CSMFiniteRatioFunction(
            function="power2_decreasing_exp",
            options_dict={"max_csm": max_csm, "alpha": alpha},
        )
        assert csm_finite_ratio.evaluate(0.0) == 1.0
        assert round(abs(csm_finite_ratio.evaluate(2.0) - 0.43807544047766522), 7) == 0
        assert round(abs(csm_finite_ratio.evaluate(4.0) - 0.15163266492815836), 7) == 0
        assert csm_finite_ratio.evaluate(8.0) == 0.0
        assert csm_finite_ratio.evaluate(9.0) == 0.0

        max_csm = 8.0
        alpha = 2.0
        csm_finite_ratio = CSMFiniteRatioFunction(
            function="power2_decreasing_exp",
            options_dict={"max_csm": max_csm, "alpha": alpha},
        )
        assert csm_finite_ratio.evaluate(0.0) == 1.0
        assert round(abs(csm_finite_ratio.evaluate(4.0) - 0.091969860292860584), 7) == 0
        assert csm_finite_ratio.evaluate(8.0) == 0.0
        assert csm_finite_ratio.evaluate(9.0) == 0.0

        max_csm = 4.0
        alpha = 1.0
        csm_finite_ratio = CSMFiniteRatioFunction(
            function="power2_decreasing_exp",
            options_dict={"max_csm": max_csm, "alpha": alpha},
        )
        assert csm_finite_ratio.evaluate(0.0) == 1.0
        assert round(abs(csm_finite_ratio.evaluate(1.0) - 0.43807544047766522), 7) == 0
        assert round(abs(csm_finite_ratio.evaluate(2.0) - 0.15163266492815836), 7) == 0
        assert csm_finite_ratio.evaluate(4.0) == 0.0
        assert csm_finite_ratio.evaluate(4.5) == 0.0

        with pytest.raises(ValueError):
            CSMFiniteRatioFunction(
                function="powern_decreasing",
                options_dict={"max_csm": max_csm, "nn": 2},
            )

    def test_CSMInfiniteRatioFunction(self):
        max_csm = 8.0
        with pytest.raises(ValueError):
            CSMInfiniteRatioFunction(
                function="power2_inverse_decreasing",
                options_dict={"max_csm": max_csm, "nn": 2},
            )

        with pytest.raises(ValueError):
            CSMInfiniteRatioFunction(
                function="power2_tangent_decreasing",
                options_dict={"max_csm": max_csm},
            )

        csm_infinite_ratio = CSMInfiniteRatioFunction(
            function="power2_inverse_decreasing", options_dict={"max_csm": max_csm}
        )
        # csm_infinite_ratio = CSMInfiniteRatioFunction(function='power2_inverse_decreasing')
        assert csm_infinite_ratio.evaluate(0.0) == np.inf
        assert round(abs(csm_infinite_ratio.evaluate(2.0) - 2.25), 7) == 0
        assert csm_infinite_ratio.evaluate(4.0) == 0.5
        assert csm_infinite_ratio.evaluate(8.0) == 0.0
        assert csm_infinite_ratio.evaluate(9.0) == 0.0

        csm_infinite_ratio = CSMInfiniteRatioFunction(
            function="power2_inverse_power2_decreasing",
            options_dict={"max_csm": max_csm},
        )
        assert csm_infinite_ratio.evaluate(0.0) == np.inf
        assert csm_infinite_ratio.evaluate(2.0) == 9.0
        assert csm_infinite_ratio.evaluate(4.0) == 1.0
        assert csm_infinite_ratio.evaluate(8.0) == 0.0
        assert csm_infinite_ratio.evaluate(9.0) == 0.0

        max_csm = 12.0
        csm_infinite_ratio = CSMInfiniteRatioFunction(
            function="power2_inverse_power2_decreasing",
            options_dict={"max_csm": max_csm},
        )

        assert csm_infinite_ratio.evaluate(0.0) == np.inf
        assert csm_infinite_ratio.evaluate(3.0) == 9.0
        assert csm_infinite_ratio.evaluate(6.0) == 1.0
        assert csm_infinite_ratio.evaluate(12.0) == 0.0
        assert csm_infinite_ratio.evaluate(20.0) == 0.0

    def test_DeltaCSMRatioFunction(self):
        with pytest.raises(ValueError):
            DeltaCSMRatioFunction(function="smoothstep", options_dict={})
        with pytest.raises(ValueError):
            DeltaCSMRatioFunction(function="smootherstep", options_dict={})
        with pytest.raises(ValueError):
            DeltaCSMRatioFunction(
                function="smootherstep",
                options_dict={"delta_csm_min": 1.0},
            )
        with pytest.raises(ValueError):
            DeltaCSMRatioFunction(
                function="smootherstep",
                options_dict={"delta_csm_max": 1.0},
            )

        delta_csm_ratio_function = DeltaCSMRatioFunction(
            function="smootherstep",
            options_dict={"delta_csm_min": 1.0, "delta_csm_max": 4.0},
        )
        assert delta_csm_ratio_function.evaluate(0.0) == 0.0
        assert delta_csm_ratio_function.evaluate(1.0) == 0.0
        assert delta_csm_ratio_function.evaluate(2.5) == 0.5
        assert delta_csm_ratio_function.evaluate(4.0) == 1.0
        assert delta_csm_ratio_function.evaluate(5.0) == 1.0

        delta_csm_ratio_function = DeltaCSMRatioFunction(
            function="smootherstep",
            options_dict={"delta_csm_min": 2.0, "delta_csm_max": 8.0},
        )
        assert delta_csm_ratio_function.evaluate(0.0) == 0.0
        assert delta_csm_ratio_function.evaluate(2.0) == 0.0
        assert delta_csm_ratio_function.evaluate(5.0) == 0.5
        assert delta_csm_ratio_function.evaluate(8.0) == 1.0
        assert delta_csm_ratio_function.evaluate(12.0) == 1.0


if __name__ == "__main__":
    unittest.main()
