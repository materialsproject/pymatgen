# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest

from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.cost import CostAnalyzer, CostDBCSV, CostDBElements


class CostAnalyzerTest(unittest.TestCase):
    def setUp(self):
        self.ca1 = CostAnalyzer(CostDBCSV(os.path.join(PymatgenTest.TEST_FILES_DIR, "costdb_1.csv")))
        self.ca2 = CostAnalyzer(CostDBCSV(os.path.join(PymatgenTest.TEST_FILES_DIR, "costdb_2.csv")))

    def test_cost_per_kg(self):
        self.assertAlmostEqual(self.ca1.get_cost_per_kg("Ag"), 3, 3)
        self.assertAlmostEqual(self.ca1.get_cost_per_kg("O"), 1, 3)
        self.assertAlmostEqual(self.ca1.get_cost_per_kg("AgO"), 2.7416, 3)
        self.assertAlmostEqual(self.ca2.get_cost_per_kg("AgO"), 1.5, 3)

    def test_cost_per_mol(self):
        self.assertAlmostEqual(self.ca1.get_cost_per_mol("Ag"), 0.3236, 3)
        self.assertAlmostEqual(self.ca1.get_cost_per_mol("O"), 0.0160, 3)
        self.assertAlmostEqual(self.ca1.get_cost_per_mol("AgO"), 0.3396, 3)
        self.assertAlmostEqual(self.ca2.get_cost_per_mol("AgO"), 0.1858, 3)

    def test_sanity(self):
        self.assertEqual(self.ca1.get_cost_per_kg("Ag"), self.ca2.get_cost_per_kg("Ag"))


class CostDBTest(unittest.TestCase):
    def test_sanity(self):
        ca = CostAnalyzer(CostDBElements())
        self.assertGreater(ca.get_cost_per_kg("PtO"), ca.get_cost_per_kg("MgO"))


if __name__ == "__main__":
    unittest.main()
