# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module contains some utility functions and classes that are used in the chemenv package.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

from typing import Dict

import numpy as np
from pymatgen.analysis.chemenv.utils.math_utils import power2_inverse_decreasing, power2_decreasing_exp
from pymatgen.analysis.chemenv.utils.math_utils import smoothstep, smootherstep
from pymatgen.analysis.chemenv.utils.math_utils import power2_inverse_power2_decreasing


class AbstractRatioFunction:
    ALLOWED_FUNCTIONS = {}  # type: Dict[str, list]

    def __init__(self, function, options_dict=None):
        if function not in self.ALLOWED_FUNCTIONS:
            raise ValueError('Function "{}" is not allowed in RatioFunction of '
                             'type "{}"'.format(function, self.__class__.__name__))
        self.eval = object.__getattribute__(self, function)
        self.function = function
        self.setup_parameters(options_dict=options_dict)

    def setup_parameters(self, options_dict):
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
                    opts = 'Option "{}"'.format(function_options[0])
                else:
                    opts1 = ', '.join(['"{}"'.format(op) for op in function_options[:-1]])
                    opts = 'Options {}'.format(' and '.join([opts1,
                                                             '"{}"'.format(function_options[-1])]))
                if options_dict is None or len(options_dict) == 0:
                    missing = 'no option was provided.'
                else:
                    optgiven = list(options_dict.keys())
                    if len(options_dict) == 1:
                        missing = 'only "{}" was provided.'.format(optgiven[0])
                    else:
                        missing1 = ', '.join(['"{}"'.format(miss) for miss in optgiven[:-1]])
                        missing = 'only {} were provided.'.format(' and '.join([missing1,
                                                                                '"{}"'.format(optgiven[-1])]))
                raise ValueError('{} should be provided for function "{}" in RatioFunction of '
                                 'type "{}" while {}'.format(opts,
                                                             self.function,
                                                             self.__class__.__name__,
                                                             missing))
            # Setup the options and raise an error if a wrong option is provided
            for key, val in options_dict.items():
                if key not in function_options:
                    raise ValueError('Option "{}" not allowed for function "{}" in RatioFunction of '
                                     'type "{}"'.format(key, self.function, self.__class__.__name__))
                self.__setattr__(key, val)

    def evaluate(self, value):
        return self.eval(value)

    @classmethod
    def from_dict(cls, dd):
        return cls(function=dd['function'], options_dict=dd['options'])


class RatioFunction(AbstractRatioFunction):
    ALLOWED_FUNCTIONS = {'power2_decreasing_exp': ['max', 'alpha'],
                         'smoothstep': ['lower', 'upper'],
                         'smootherstep': ['lower', 'upper'],
                         'inverse_smoothstep': ['lower', 'upper'],
                         'inverse_smootherstep': ['lower', 'upper'],
                         'power2_inverse_decreasing': ['max'],
                         'power2_inverse_power2_decreasing': ['max']
                         }

    def power2_decreasing_exp(self, vals):
        return power2_decreasing_exp(vals, edges=[0.0, self.__dict__['max']], alpha=self.__dict__['alpha'])

    def smootherstep(self, vals):
        return smootherstep(vals, edges=[self.__dict__['lower'], self.__dict__['upper']])

    def smoothstep(self, vals):
        return smoothstep(vals, edges=[self.__dict__['lower'], self.__dict__['upper']])

    def inverse_smootherstep(self, vals):
        return smootherstep(vals, edges=[self.__dict__['lower'], self.__dict__['upper']], inverse=True)

    def inverse_smoothstep(self, vals):
        return smoothstep(vals, edges=[self.__dict__['lower'], self.__dict__['upper']], inverse=True)

    def power2_inverse_decreasing(self, vals):
        return power2_inverse_decreasing(vals, edges=[0.0, self.__dict__['max']])

    def power2_inverse_power2_decreasing(self, vals):
        return power2_inverse_power2_decreasing(vals, edges=[0.0, self.__dict__['max']])


class CSMFiniteRatioFunction(AbstractRatioFunction):
    ALLOWED_FUNCTIONS = {'power2_decreasing_exp': ['max_csm', 'alpha'],
                         'smoothstep': ['lower_csm', 'upper_csm'],
                         'smootherstep': ['lower_csm', 'upper_csm']
                         }

    def power2_decreasing_exp(self, vals):
        return power2_decreasing_exp(vals, edges=[0.0, self.__dict__['max_csm']], alpha=self.__dict__['alpha'])

    def smootherstep(self, vals):
        return smootherstep(vals, edges=[self.__dict__['lower_csm'], self.__dict__['upper_csm']], inverse=True)

    def smoothstep(self, vals):
        return smootherstep(vals, edges=[self.__dict__['lower_csm'], self.__dict__['upper_csm']], inverse=True)

    def fractions(self, data):
        if len(data) == 0:
            return None
        total = np.sum([self.eval(dd) for dd in data])
        if total > 0.0:
            return [self.eval(dd) / total for dd in data]
        else:
            return None

    def mean_estimator(self, data):
        if len(data) == 0:
            return None
        elif len(data) == 1:
            return data[0]
        else:
            fractions = self.fractions(data)
            if fractions is None:
                return None
            return np.sum(np.array(fractions) * np.array(data))

    ratios = fractions


class CSMInfiniteRatioFunction(AbstractRatioFunction):
    ALLOWED_FUNCTIONS = {'power2_inverse_decreasing': ['max_csm'],
                         'power2_inverse_power2_decreasing': ['max_csm']}

    def power2_inverse_decreasing(self, vals):
        return power2_inverse_decreasing(vals, edges=[0.0, self.__dict__['max_csm']])

    def power2_inverse_power2_decreasing(self, vals):
        return power2_inverse_power2_decreasing(vals, edges=[0.0, self.__dict__['max_csm']])

    def fractions(self, data):
        if len(data) == 0:
            return None
        close_to_zero = np.isclose(data, 0.0, atol=1e-10).tolist()
        nzeros = close_to_zero.count(True)
        if nzeros == 1:
            fractions = [0.0] * len(data)
            fractions[close_to_zero.index(True)] = 1.0
            return fractions
        elif nzeros > 1:
            raise RuntimeError('Should not have more than one continuous symmetry measure with value equal to 0.0')
        else:
            fractions = self.eval(np.array(data))
            total = np.sum(fractions)
            if total > 0.0:
                return fractions / total
            else:
                return None

    def mean_estimator(self, data):
        if len(data) == 0:
            return None
        elif len(data) == 1:
            return data[0]
        else:
            fractions = self.fractions(data)
            if fractions is None:
                return None
            return np.sum(np.array(fractions) * np.array(data))

    ratios = fractions


class DeltaCSMRatioFunction(AbstractRatioFunction):
    ALLOWED_FUNCTIONS = {'smootherstep': ['delta_csm_min', 'delta_csm_max']}

    def smootherstep(self, vals):
        return smootherstep(vals, edges=[self.__dict__['delta_csm_min'], self.__dict__['delta_csm_max']])
