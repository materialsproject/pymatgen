# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
from pymatgen.analysis.hhi.hhi import HHIModel


class HHIModelTest(unittest.TestCase):

    def test_hhi(self):
        hhi = HHIModel()
        self.assertEqual(hhi.get_hhi("He"), (3200, 3900))
        self.assertEqual(hhi.get_hhi_production("He"), 3200)
        self.assertEqual(hhi.get_hhi_reserve("He"), 3900)

        self.assertAlmostEqual(hhi.get_hhi_production("Li2O"), 1614.96, 1)
        self.assertAlmostEqual(hhi.get_hhi_reserve("Li2O"), 2218.90, 1)

        self.assertEqual(hhi.get_hhi_designation(1400), "low")
        self.assertEqual(hhi.get_hhi_designation(1800), "medium")
        self.assertEqual(hhi.get_hhi_designation(3000), "high")
        self.assertEqual(hhi.get_hhi_designation(None), None)


if __name__ == "__main__":
    unittest.main()