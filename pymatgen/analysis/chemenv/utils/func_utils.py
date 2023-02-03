# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module contains some utility functions and classes that are used in the chemenv package.
"""

from __future__ import annotations

import numpy as np

from pymatgen.analysis.chemenv.utils.math_utils import (
    power2_decreasing_exp,
    power2_inverse_decreasing,
    power2_inverse_power2_decreasing,
    smootherstep,
    smoothstep,
)

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"


class AbstractRatioFunction:
    """
    Abstract class for all ratio functions
    """

    ALLOWED_FUNCTIONS: dict[str, list] = {}

    def __init__(self, function, options_dict=None):
        """Constructor for AbstractRatioFunction

        :param function: Ration function name.
        :param options_dict: Dictionary containing the parameters for the ratio function.
        """
        if function not in self.ALLOWED_FUNCTIONS:
            raise ValueError(f'Function {function!r} is not allowed in RatioFunction of type "{type(self).__name__}"')
        self.eval = object.__getattribute__(self, function)
        self.function = function
        self.setup_parameters(options_dict=options_dict)

    def setup_parameters(self, options_dict):
        """Set up the parameters for this ratio function.

        :param options_dict: Dictionary containing the parameters for the ratio function.
        :return: None.
        """
        function_options = self.ALLOWED_FUNCTIONS[self.function]
        if len(function_options) > 0:
            # Check if there are missing options
            if options_dict is None:
                missing_options = True
            else:
                missing_options = False
                for op in function_options:
                    if op not in options_dict:
                        missing_options = True
                        break
            # If there are missing options, raise an error
            if missing_options:
                if len(function_options) == 1:
                    opts = f'Option "{function_options[0]}"'
                else:
                    opts1 = ", ".join(f"{op!r}" for op in function_options[:-1])
                    opts = f'Options {opts1} and "{function_options[-1]}"'
                if options_dict is None or len(options_dict) == 0:
                    missing = "no option was provided."
                else:
                    optgiven = list(options_dict)
                    if len(options_dict) == 1:
                        missing = f"only {optgiven[0]!r} was provided."
                    else:
                        missing1 = ", ".join(f"{miss!r}" for miss in optgiven[:-1])
                        missing = f"only {missing1} and {optgiven[-1]!r} were provided."
                raise ValueError(
                    f"{opts} should be provided for function {self.function!r} in RatioFunction of "
                    f"type {type(self).__name__!r} while {missing}"
                )
            # Setup the options and raise an error if a wrong option is provided
            for key, val in options_dict.items():
                if key not in function_options:
                    raise ValueError(
                        f"Option {key!r} not allowed for function {self.function!r} in RatioFunction of "
                        f'type "{type(self).__name__}"'
                    )
                setattr(self, key, val)

    def evaluate(self, value):
        """Evaluate the ratio function for the given value.

        :param value: Value for which ratio function has to be evaluated.
        :return: Ratio function corresponding to the value.
        """
        return self.eval(value)

    @classmethod
    def from_dict(cls, dd):
        """Construct ratio function from dict.

        :param dd: Dict representation of the ratio function
        :return: Ratio function object.
        """
        return cls(function=dd["function"], options_dict=dd["options"])


class RatioFunction(AbstractRatioFunction):
    """Concrete implementation of a series of ratio functions."""

    ALLOWED_FUNCTIONS = {
        "power2_decreasing_exp": ["max", "alpha"],
        "smoothstep": ["lower", "upper"],
        "smootherstep": ["lower", "upper"],
        "inverse_smoothstep": ["lower", "upper"],
        "inverse_smootherstep": ["lower", "upper"],
        "power2_inverse_decreasing": ["max"],
        "power2_inverse_power2_decreasing": ["max"],
    }

    def power2_decreasing_exp(self, vals):
        """Get the evaluation of the ratio function f(x)=exp(-a*x)*(x-1)^2.

        The values (i.e. "x"), are scaled to the "max" parameter. The "a" constant
        correspond to the "alpha" parameter.

        :param vals: Values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the values.
        """
        return power2_decreasing_exp(vals, edges=[0.0, self.__dict__["max"]], alpha=self.__dict__["alpha"])

    def smootherstep(self, vals):
        """Get the evaluation of the smootherstep ratio function: f(x)=6*x^5-15*x^4+10*x^3.

        The values (i.e. "x"), are scaled between the "lower" and "upper" parameters.

        :param vals: Values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the values.
        """
        return smootherstep(vals, edges=[self.__dict__["lower"], self.__dict__["upper"]])

    def smoothstep(self, vals):
        """Get the evaluation of the smoothstep ratio function: f(x)=3*x^2-2*x^3.

        The values (i.e. "x"), are scaled between the "lower" and "upper" parameters.

        :param vals: Values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the values.
        """
        return smoothstep(vals, edges=[self.__dict__["lower"], self.__dict__["upper"]])

    def inverse_smootherstep(self, vals):
        """Get the evaluation of the "inverse" smootherstep ratio function: f(x)=1-(6*x^5-15*x^4+10*x^3).

        The values (i.e. "x"), are scaled between the "lower" and "upper" parameters.

        :param vals: Values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the values.
        """
        return smootherstep(vals, edges=[self.__dict__["lower"], self.__dict__["upper"]], inverse=True)

    def inverse_smoothstep(self, vals):
        """Get the evaluation of the "inverse" smoothstep ratio function: f(x)=1-(3*x^2-2*x^3).

        The values (i.e. "x"), are scaled between the "lower" and "upper" parameters.

        :param vals: Values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the values.
        """
        return smoothstep(vals, edges=[self.__dict__["lower"], self.__dict__["upper"]], inverse=True)

    def power2_inverse_decreasing(self, vals):
        """Get the evaluation of the ratio function f(x)=(x-1)^2 / x.

        The values (i.e. "x"), are scaled to the "max" parameter.

        :param vals: Values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the values.
        """
        return power2_inverse_decreasing(vals, edges=[0.0, self.__dict__["max"]])

    def power2_inverse_power2_decreasing(self, vals):
        """Get the evaluation of the ratio function f(x)=(x-1)^2 / x^2.

        The values (i.e. "x"), are scaled to the "max" parameter.

        :param vals: Values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the values.
        """
        return power2_inverse_power2_decreasing(vals, edges=[0.0, self.__dict__["max"]])


class CSMFiniteRatioFunction(AbstractRatioFunction):
    """Concrete implementation of a series of ratio functions applied to the continuous symmetry measure (CSM).

    Uses "finite" ratio functions.

    See the following reference for details:
    ChemEnv: a fast and robust coordination environment identification tool,
    D. Waroquiers et al., Acta Cryst. B 76, 683 (2020).
    """

    ALLOWED_FUNCTIONS = {
        "power2_decreasing_exp": ["max_csm", "alpha"],
        "smoothstep": ["lower_csm", "upper_csm"],
        "smootherstep": ["lower_csm", "upper_csm"],
    }

    def power2_decreasing_exp(self, vals):
        """Get the evaluation of the ratio function f(x)=exp(-a*x)*(x-1)^2.

        The CSM values (i.e. "x"), are scaled to the "max_csm" parameter. The "a" constant
        correspond to the "alpha" parameter.

        :param vals: CSM values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the CSM values.
        """
        return power2_decreasing_exp(vals, edges=[0.0, self.__dict__["max_csm"]], alpha=self.__dict__["alpha"])

    def smootherstep(self, vals):
        """Get the evaluation of the smootherstep ratio function: f(x)=6*x^5-15*x^4+10*x^3.

        The CSM values (i.e. "x"), are scaled between the "lower_csm" and "upper_csm" parameters.

        :param vals: CSM values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the CSM values.
        """
        return smootherstep(
            vals,
            edges=[self.__dict__["lower_csm"], self.__dict__["upper_csm"]],
            inverse=True,
        )

    def smoothstep(self, vals):
        """Get the evaluation of the smoothstep ratio function: f(x)=3*x^2-2*x^3.

        The CSM values (i.e. "x"), are scaled between the "lower_csm" and "upper_csm" parameters.

        :param vals: CSM values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the CSM values.
        """
        return smootherstep(
            vals,
            edges=[self.__dict__["lower_csm"], self.__dict__["upper_csm"]],
            inverse=True,
        )

    def fractions(self, data):
        """Get the fractions from the CSM ratio function applied to the data.

        :param data: List of CSM values to estimate fractions.
        :return: Corresponding fractions for each CSM.
        """
        if len(data) == 0:
            return None
        total = np.sum([self.eval(dd) for dd in data])
        if total > 0.0:
            return [self.eval(dd) / total for dd in data]
        return None

    def mean_estimator(self, data):
        """Get the weighted CSM using this CSM ratio function applied to the data.

        :param data: List of CSM values to estimate the weighted CSM.
        :return: Weighted CSM from this ratio function.
        """
        if len(data) == 0:
            return None
        if len(data) == 1:
            return data[0]
        fractions = self.fractions(data)
        if fractions is None:
            return None
        return np.sum(np.array(fractions) * np.array(data))

    ratios = fractions


class CSMInfiniteRatioFunction(AbstractRatioFunction):
    """Concrete implementation of a series of ratio functions applied to the continuous symmetry measure (CSM).

    Uses "infinite" ratio functions.

    See the following reference for details:
    ChemEnv: a fast and robust coordination environment identification tool,
    D. Waroquiers et al., Acta Cryst. B 76, 683 (2020).
    """

    ALLOWED_FUNCTIONS = {
        "power2_inverse_decreasing": ["max_csm"],
        "power2_inverse_power2_decreasing": ["max_csm"],
    }

    def power2_inverse_decreasing(self, vals):
        """Get the evaluation of the ratio function f(x)=(x-1)^2 / x.

        The CSM values (i.e. "x"), are scaled to the "max_csm" parameter. The "a" constant
        correspond to the "alpha" parameter.

        :param vals: CSM values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the CSM values.
        """
        return power2_inverse_decreasing(vals, edges=[0.0, self.__dict__["max_csm"]])

    def power2_inverse_power2_decreasing(self, vals):
        """Get the evaluation of the ratio function f(x)=(x-1)^2 / x^2.

        The CSM values (i.e. "x"), are scaled to the "max_csm" parameter. The "a" constant
        correspond to the "alpha" parameter.

        :param vals: CSM values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the CSM values.
        """
        return power2_inverse_power2_decreasing(vals, edges=[0.0, self.__dict__["max_csm"]])

    def fractions(self, data):
        """Get the fractions from the CSM ratio function applied to the data.

        :param data: List of CSM values to estimate fractions.
        :return: Corresponding fractions for each CSM.
        """
        if len(data) == 0:
            return None
        close_to_zero = np.isclose(data, 0.0, atol=1e-10).tolist()
        nzeros = close_to_zero.count(True)
        if nzeros == 1:
            fractions = [0.0] * len(data)
            fractions[close_to_zero.index(True)] = 1.0
            return fractions
        if nzeros > 1:
            raise RuntimeError("Should not have more than one continuous symmetry measure with value equal to 0.0")
        fractions = self.eval(np.array(data))
        total = np.sum(fractions)
        if total > 0.0:
            return fractions / total
        return None

    def mean_estimator(self, data):
        """Get the weighted CSM using this CSM ratio function applied to the data.

        :param data: List of CSM values to estimate the weighted CSM.
        :return: Weighted CSM from this ratio function.
        """
        if len(data) == 0:
            return None
        if len(data) == 1:
            return data[0]
        fractions = self.fractions(data)
        if fractions is None:
            return None
        return np.sum(np.array(fractions) * np.array(data))

    ratios = fractions


class DeltaCSMRatioFunction(AbstractRatioFunction):
    """
    Concrete implementation of a series of ratio functions applied to differences of
    continuous symmetry measures (DeltaCSM).

    Uses "finite" ratio functions.

    See the following reference for details:
    ChemEnv: a fast and robust coordination environment identification tool,
    D. Waroquiers et al., Acta Cryst. B 76, 683 (2020).
    """

    ALLOWED_FUNCTIONS = {"smootherstep": ["delta_csm_min", "delta_csm_max"]}

    def smootherstep(self, vals):
        """Get the evaluation of the smootherstep ratio function: f(x)=6*x^5-15*x^4+10*x^3.

        The DeltaCSM values (i.e. "x"), are scaled between the "delta_csm_min" and "delta_csm_max" parameters.

        :param vals: DeltaCSM values for which the ratio function has to be evaluated.
        :return: Result of the ratio function applied to the DeltaCSM values.
        """
        return smootherstep(vals, edges=[self.__dict__["delta_csm_min"], self.__dict__["delta_csm_max"]])
