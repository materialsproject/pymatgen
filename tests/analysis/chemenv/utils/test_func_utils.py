from __future__ import annotations

import unittest

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.chemenv.utils.func_utils import (
    CSMFiniteRatioFunction,
    CSMInfiniteRatioFunction,
    DeltaCSMRatioFunction,
)

__author__ = "waroquiers"


class TestFuncUtils(unittest.TestCase):
    def test_CSMFiniteRatioFunction(self):
        max_csm = 8
        alpha = 1
        csm_finite_ratio = CSMFiniteRatioFunction(
            function="power2_decreasing_exp",
            options_dict={"max_csm": max_csm, "alpha": alpha},
        )
        assert csm_finite_ratio.evaluate(0) == 1
        assert csm_finite_ratio.evaluate(2) == approx(0.4380754)
        assert csm_finite_ratio.evaluate(4) == approx(0.1516326)
        assert csm_finite_ratio.evaluate(8) == 0
        assert csm_finite_ratio.evaluate(9) == 0

        max_csm = 8
        alpha = 2
        csm_finite_ratio = CSMFiniteRatioFunction(
            function="power2_decreasing_exp",
            options_dict={"max_csm": max_csm, "alpha": alpha},
        )
        assert csm_finite_ratio.evaluate(0) == 1
        assert csm_finite_ratio.evaluate(4) == approx(0.09196986)
        assert csm_finite_ratio.evaluate(8) == 0
        assert csm_finite_ratio.evaluate(9) == 0

        max_csm = 4
        alpha = 1
        csm_finite_ratio = CSMFiniteRatioFunction(
            function="power2_decreasing_exp",
            options_dict={"max_csm": max_csm, "alpha": alpha},
        )
        assert csm_finite_ratio.evaluate(0) == 1
        assert csm_finite_ratio.evaluate(1) == approx(0.4380754)
        assert csm_finite_ratio.evaluate(2) == approx(0.1516326)
        assert csm_finite_ratio.evaluate(4) == 0
        assert csm_finite_ratio.evaluate(4.5) == 0

        with pytest.raises(
            ValueError,
            match="function='powern_decreasing' is not allowed in RatioFunction of type CSMFiniteRatioFunction",
        ):
            CSMFiniteRatioFunction(
                function="powern_decreasing",
                options_dict={"max_csm": max_csm, "nn": 2},
            )

    def test_CSMInfiniteRatioFunction(self):
        max_csm = 8
        with pytest.raises(ValueError, match="Option 'nn' not allowed for function 'power2_inverse_decreasing' in "):
            CSMInfiniteRatioFunction(
                function="power2_inverse_decreasing",
                options_dict={"max_csm": max_csm, "nn": 2},
            )

        with pytest.raises(
            ValueError, match="function='power2_tangent_decreasing' is not allowed in RatioFunction of type "
        ):
            CSMInfiniteRatioFunction(
                function="power2_tangent_decreasing",
                options_dict={"max_csm": max_csm},
            )

        csm_infinite_ratio = CSMInfiniteRatioFunction(
            function="power2_inverse_decreasing", options_dict={"max_csm": max_csm}
        )
        # csm_infinite_ratio = CSMInfiniteRatioFunction(function='power2_inverse_decreasing')
        assert csm_infinite_ratio.evaluate(0) == np.inf
        assert csm_infinite_ratio.evaluate(2) == approx(2.25)
        assert csm_infinite_ratio.evaluate(4) == 0.5
        assert csm_infinite_ratio.evaluate(8) == 0
        assert csm_infinite_ratio.evaluate(9) == 0

        csm_infinite_ratio = CSMInfiniteRatioFunction(
            function="power2_inverse_power2_decreasing",
            options_dict={"max_csm": max_csm},
        )
        assert csm_infinite_ratio.evaluate(0) == np.inf
        assert csm_infinite_ratio.evaluate(2) == 9
        assert csm_infinite_ratio.evaluate(4) == 1
        assert csm_infinite_ratio.evaluate(8) == 0
        assert csm_infinite_ratio.evaluate(9) == 0

        max_csm = 12
        csm_infinite_ratio = CSMInfiniteRatioFunction(
            function="power2_inverse_power2_decreasing",
            options_dict={"max_csm": max_csm},
        )

        assert csm_infinite_ratio.evaluate(0) == np.inf
        assert csm_infinite_ratio.evaluate(3) == 9
        assert csm_infinite_ratio.evaluate(6) == 1
        assert csm_infinite_ratio.evaluate(12) == 0
        assert csm_infinite_ratio.evaluate(20) == 0

    def test_DeltaCSMRatioFunction(self):
        with pytest.raises(ValueError, match="function='smoothstep' is not allowed in RatioFunction of typ"):
            DeltaCSMRatioFunction(function="smoothstep", options_dict={})
        with pytest.raises(
            ValueError, match="Options 'delta_csm_min' and \"delta_csm_max\" should be provided for function"
        ):
            DeltaCSMRatioFunction(function="smootherstep", options_dict={})
        with pytest.raises(
            ValueError,
            match="Options 'delta_csm_min' and \"delta_csm_max\" should be provided for function 'smootherstep'",
        ):
            DeltaCSMRatioFunction(
                function="smootherstep",
                options_dict={"delta_csm_min": 1},
            )
        with pytest.raises(
            ValueError,
            match="Options 'delta_csm_min' and \"delta_csm_max\" should be provided for function 'smootherstep'",
        ):
            DeltaCSMRatioFunction(
                function="smootherstep",
                options_dict={"delta_csm_max": 1},
            )

        delta_csm_ratio_function = DeltaCSMRatioFunction(
            function="smootherstep",
            options_dict={"delta_csm_min": 1, "delta_csm_max": 4},
        )
        assert delta_csm_ratio_function.evaluate(0) == 0
        assert delta_csm_ratio_function.evaluate(1) == 0
        assert delta_csm_ratio_function.evaluate(2.5) == 0.5
        assert delta_csm_ratio_function.evaluate(4) == 1
        assert delta_csm_ratio_function.evaluate(5) == 1

        delta_csm_ratio_function = DeltaCSMRatioFunction(
            function="smootherstep",
            options_dict={"delta_csm_min": 2, "delta_csm_max": 8},
        )
        assert delta_csm_ratio_function.evaluate(0) == 0
        assert delta_csm_ratio_function.evaluate(2) == 0
        assert delta_csm_ratio_function.evaluate(5) == 0.5
        assert delta_csm_ratio_function.evaluate(8) == 1
        assert delta_csm_ratio_function.evaluate(12) == 1
