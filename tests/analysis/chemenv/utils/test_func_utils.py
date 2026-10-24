from __future__ import annotations

import re

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.chemenv.utils.func_utils import (
    CSMFiniteRatioFunction,
    CSMInfiniteRatioFunction,
    DeltaCSMRatioFunction,
    RatioFunction,
)

__author__ = "waroquiers"


class TestFuncUtils:
    def test_csm_finite_ratio_function(self):
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

    def test_csm_infinite_ratio_function(self):
        max_csm = 8
        with pytest.raises(
            ValueError,
            match="Option 'nn' not allowed for function 'power2_inverse_decreasing' in ",
        ):
            CSMInfiniteRatioFunction(
                function="power2_inverse_decreasing",
                options_dict={"max_csm": max_csm, "nn": 2},
            )

        with pytest.raises(
            ValueError,
            match="function='power2_tangent_decreasing' is not allowed in RatioFunction of type ",
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
        assert csm_infinite_ratio.evaluate(4) == approx(0.5)
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

    def test_delta_csm_ratio_function(self):
        with pytest.raises(
            ValueError,
            match="function='smoothstep' is not allowed in RatioFunction of typ",
        ):
            DeltaCSMRatioFunction(function="smoothstep", options_dict={})
        with pytest.raises(
            ValueError,
            match="Options 'delta_csm_min' and \"delta_csm_max\" should be provided for function",
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
        assert delta_csm_ratio_function.evaluate(2.5) == approx(0.5)
        assert delta_csm_ratio_function.evaluate(4) == 1
        assert delta_csm_ratio_function.evaluate(5) == 1

        delta_csm_ratio_function = DeltaCSMRatioFunction(
            function="smootherstep",
            options_dict={"delta_csm_min": 2, "delta_csm_max": 8},
        )
        assert delta_csm_ratio_function.evaluate(0) == 0
        assert delta_csm_ratio_function.evaluate(2) == 0
        assert delta_csm_ratio_function.evaluate(5) == approx(0.5)
        assert delta_csm_ratio_function.evaluate(8) == 1
        assert delta_csm_ratio_function.evaluate(12) == 1


@pytest.mark.parametrize(
    ("function", "options", "value", "expected"),
    [
        # (0.5 - 1)^2 * exp(-0.5)
        ("power2_decreasing_exp", {"max": 2, "alpha": 1}, 1.0, 0.25 * np.exp(-0.5)),
        ("smoothstep", {"lower": 0, "upper": 2}, 0.5, 0.15625),
        ("inverse_smoothstep", {"lower": 0, "upper": 2}, 0.5, 0.84375),
        ("smootherstep", {"lower": 0, "upper": 2}, 0.5, 0.103515625),
        ("inverse_smootherstep", {"lower": 0, "upper": 2}, 0.5, 0.896484375),
        # (0.5 - 1)^2 / 0.5
        ("power2_inverse_decreasing", {"max": 2}, 1.0, 0.5),
        # (0.5 - 1)^2 / 0.5^2
        ("power2_inverse_power2_decreasing", {"max": 2}, 1.0, 1.0),
    ],
)
def test_ratio_function(function, options, value, expected):
    ratio_function = RatioFunction(function=function, options_dict=options)
    assert ratio_function.evaluate(value) == approx(expected)
    from_dict = RatioFunction.from_dict({"function": function, "options": options})
    assert from_dict.function == function
    assert from_dict.evaluate(value) == approx(expected)


@pytest.mark.parametrize(
    ("function", "options", "match"),
    [
        (
            "power2_inverse_decreasing",
            None,
            (
                "Option \"max\" should be provided for function 'power2_inverse_decreasing' in RatioFunction of type "
                "'RatioFunction' while no option was provided."
            ),
        ),
        ("smoothstep", {}, "while no option was provided."),
        ("smoothstep", {"lower": 0}, "while only 'lower' was provided."),
        ("smoothstep", {"low": 0, "up": 1}, "while only 'low' and 'up' were provided."),
        ("smoothstep", {"lower": 0, "upper": 1, "extra": 2}, "Option 'extra' not allowed for function 'smoothstep'"),
    ],
)
def test_ratio_function_invalid_options(function, options, match):
    with pytest.raises(ValueError, match=re.escape(match)):
        RatioFunction(function=function, options_dict=options)


def test_csm_finite_ratio_function_steps():
    smootherstep = CSMFiniteRatioFunction("smootherstep", {"lower_csm": 1, "upper_csm": 3})
    assert smootherstep.evaluate(0.5) == 1
    assert smootherstep.evaluate(2) == approx(0.5)
    assert smootherstep.evaluate(2.5) == approx(0.103515625)
    assert smootherstep.evaluate(4) == 0
    smoothstep = CSMFiniteRatioFunction("smoothstep", {"lower_csm": 1, "upper_csm": 3})
    assert smoothstep.evaluate(1) == 1
    assert smoothstep.evaluate(2) == approx(0.5)
    assert smoothstep.evaluate(3) == 0


def test_csm_finite_ratio_function_fractions_and_mean():
    rf = CSMFiniteRatioFunction("power2_decreasing_exp", {"max_csm": 8, "alpha": 1})
    w2 = rf.evaluate(2)
    assert w2 == approx(0.4380754)
    assert rf.fractions([0, 2]) == approx([1 / (1 + w2), w2 / (1 + w2)])
    assert rf.ratios([0, 2]) == approx([1 / (1 + w2), w2 / (1 + w2)])
    assert rf.fractions([]) is None
    assert rf.fractions([8, 9]) is None  # all weights vanish
    assert rf.mean_estimator([]) is None
    assert rf.mean_estimator([3.0]) == 3.0
    assert rf.mean_estimator([8, 9]) is None
    assert rf.mean_estimator([0, 2]) == approx(2 * w2 / (1 + w2))


def test_csm_infinite_ratio_function_fractions_and_mean():
    rf = CSMInfiniteRatioFunction("power2_inverse_decreasing", {"max_csm": 8})
    # f(2) = (0.25 - 1)^2 / 0.25 = 2.25 ; f(4) = (0.5 - 1)^2 / 0.5 = 0.5
    assert rf.fractions([2, 4]) == approx([2.25 / 2.75, 0.5 / 2.75])
    assert rf.ratios([2, 4]) == approx([2.25 / 2.75, 0.5 / 2.75])
    assert rf.fractions([2, 0.0, 4]) == [0.0, 1.0, 0.0]
    with pytest.raises(RuntimeError, match="Should not have more than one continuous symmetry measure"):
        rf.fractions([0.0, 1e-12, 3])
    assert rf.fractions([]) is None
    assert rf.fractions([8, 9]) is None
    assert rf.mean_estimator([]) is None
    assert rf.mean_estimator([5.0]) == 5.0
    assert rf.mean_estimator([8, 9]) is None
    assert rf.mean_estimator([2, 4]) == approx((2 * 2.25 + 4 * 0.5) / 2.75)
