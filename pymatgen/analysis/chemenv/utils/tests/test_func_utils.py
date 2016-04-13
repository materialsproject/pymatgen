#!/usr/bin/env python


__author__ = 'waroquiers'

import unittest2 as unittest
from pymatgen.analysis.chemenv.utils.func_utils import CSMFiniteRatioFunction
from pymatgen.analysis.chemenv.utils.func_utils import CSMInfiniteRatioFunction
from pymatgen.analysis.chemenv.utils.func_utils import DeltaCSMRatioFunction
import numpy as np


class FuncUtilsTest(unittest.TestCase):

    def test_CSMFiniteRatioFunction(self):
        max_csm = 8.0
        alpha = 1.0
        csm_finite_ratio = CSMFiniteRatioFunction(function='power2_decreasing_exp',
                                                  options_dict={'max_csm': max_csm,
                                                                'alpha': alpha})
        self.assertEqual(csm_finite_ratio.evaluate(0.0), 1.0)
        self.assertEqual(csm_finite_ratio.evaluate(2.0), 0.43807544047766522)
        self.assertEqual(csm_finite_ratio.evaluate(4.0), 0.15163266492815836)
        self.assertEqual(csm_finite_ratio.evaluate(8.0), 0.0)
        self.assertEqual(csm_finite_ratio.evaluate(9.0), 0.0)

        max_csm = 8.0
        alpha = 2.0
        csm_finite_ratio = CSMFiniteRatioFunction(function='power2_decreasing_exp',
                                                  options_dict={'max_csm': max_csm,
                                                                'alpha': alpha})
        self.assertEqual(csm_finite_ratio.evaluate(0.0), 1.0)
        self.assertEqual(csm_finite_ratio.evaluate(4.0), 0.091969860292860584)
        self.assertEqual(csm_finite_ratio.evaluate(8.0), 0.0)
        self.assertEqual(csm_finite_ratio.evaluate(9.0), 0.0)

        max_csm = 4.0
        alpha = 1.0
        csm_finite_ratio = CSMFiniteRatioFunction(function='power2_decreasing_exp',
                                                  options_dict={'max_csm': max_csm,
                                                                'alpha': alpha})
        self.assertEqual(csm_finite_ratio.evaluate(0.0), 1.0)
        self.assertEqual(csm_finite_ratio.evaluate(1.0), 0.43807544047766522)
        self.assertEqual(csm_finite_ratio.evaluate(2.0), 0.15163266492815836)
        self.assertEqual(csm_finite_ratio.evaluate(4.0), 0.0)
        self.assertEqual(csm_finite_ratio.evaluate(4.5), 0.0)

        self.assertRaises(ValueError, CSMFiniteRatioFunction,
                          function='powern_decreasing',
                          options_dict={'max_csm': max_csm,
                                        'nn': 2})

    def test_CSMInfiniteRatioFunction(self):
        max_csm = 8.0
        self.assertRaises(ValueError, CSMInfiniteRatioFunction,
                          function='power2_inverse_decreasing',
                          options_dict={'max_csm': max_csm,
                                        'nn': 2})

        self.assertRaises(ValueError, CSMInfiniteRatioFunction,
                          function='power2_tangent_decreasing',
                          options_dict={'max_csm': max_csm})

        csm_infinite_ratio = CSMInfiniteRatioFunction(function='power2_inverse_decreasing',
                                                      options_dict={'max_csm': max_csm})
        # csm_infinite_ratio = CSMInfiniteRatioFunction(function='power2_inverse_decreasing')
        self.assertEqual(csm_infinite_ratio.evaluate(0.0), np.inf)
        self.assertEqual(csm_infinite_ratio.evaluate(2.0), 2.25)
        self.assertEqual(csm_infinite_ratio.evaluate(4.0), 0.5)
        self.assertEqual(csm_infinite_ratio.evaluate(8.0), 0.0)
        self.assertEqual(csm_infinite_ratio.evaluate(9.0), 0.0)

        csm_infinite_ratio = CSMInfiniteRatioFunction(function='power2_inverse_power2_decreasing',
                                                      options_dict={'max_csm': max_csm})
        self.assertEqual(csm_infinite_ratio.evaluate(0.0), np.inf)
        self.assertEqual(csm_infinite_ratio.evaluate(2.0), 9.0)
        self.assertEqual(csm_infinite_ratio.evaluate(4.0), 1.0)
        self.assertEqual(csm_infinite_ratio.evaluate(8.0), 0.0)
        self.assertEqual(csm_infinite_ratio.evaluate(9.0), 0.0)

        max_csm = 12.0
        csm_infinite_ratio = CSMInfiniteRatioFunction(function='power2_inverse_power2_decreasing',
                                                      options_dict={'max_csm': max_csm})

        self.assertEqual(csm_infinite_ratio.evaluate(0.0), np.inf)
        self.assertEqual(csm_infinite_ratio.evaluate(3.0), 9.0)
        self.assertEqual(csm_infinite_ratio.evaluate(6.0), 1.0)
        self.assertEqual(csm_infinite_ratio.evaluate(12.0), 0.0)
        self.assertEqual(csm_infinite_ratio.evaluate(20.0), 0.0)

    def test_DeltaCSMRatioFunction(self):
        self.assertRaises(ValueError, DeltaCSMRatioFunction,
                          function='smoothstep',
                          options_dict={})
        self.assertRaises(ValueError, DeltaCSMRatioFunction,
                          function='smootherstep',
                          options_dict={})
        self.assertRaises(ValueError, DeltaCSMRatioFunction,
                          function='smootherstep',
                          options_dict={'delta_csm_min': 1.0})
        self.assertRaises(ValueError, DeltaCSMRatioFunction,
                          function='smootherstep',
                          options_dict={'delta_csm_max': 1.0})

        delta_csm_ratio_function = DeltaCSMRatioFunction(function='smootherstep',
                                                         options_dict={'delta_csm_min': 1.0, 'delta_csm_max': 4.0})
        self.assertEqual(delta_csm_ratio_function.evaluate(0.0), 0.0)
        self.assertEqual(delta_csm_ratio_function.evaluate(1.0), 0.0)
        self.assertEqual(delta_csm_ratio_function.evaluate(2.5), 0.5)
        self.assertEqual(delta_csm_ratio_function.evaluate(4.0), 1.0)
        self.assertEqual(delta_csm_ratio_function.evaluate(5.0), 1.0)

        delta_csm_ratio_function = DeltaCSMRatioFunction(function='smootherstep',
                                                         options_dict={'delta_csm_min': 2.0, 'delta_csm_max': 8.0})
        self.assertEqual(delta_csm_ratio_function.evaluate(0.0), 0.0)
        self.assertEqual(delta_csm_ratio_function.evaluate(2.0), 0.0)
        self.assertEqual(delta_csm_ratio_function.evaluate(5.0), 0.5)
        self.assertEqual(delta_csm_ratio_function.evaluate(8.0), 1.0)
        self.assertEqual(delta_csm_ratio_function.evaluate(12.0), 1.0)

if __name__ == "__main__":
    unittest.main()