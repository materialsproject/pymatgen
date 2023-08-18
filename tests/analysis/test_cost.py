from __future__ import annotations

import unittest

from pytest import approx

from pymatgen.analysis.cost import CostAnalyzer, CostDBCSV, CostDBElements
from pymatgen.util.testing import TEST_FILES_DIR


class TestCostAnalyzer(unittest.TestCase):
    def setUp(self):
        self.ca1 = CostAnalyzer(CostDBCSV(f"{TEST_FILES_DIR}/costdb_1.csv"))
        self.ca2 = CostAnalyzer(CostDBCSV(f"{TEST_FILES_DIR}/costdb_2.csv"))

    def test_cost_per_kg(self):
        assert self.ca1.get_cost_per_kg("Ag") == approx(3, rel=1e-3)
        assert self.ca1.get_cost_per_kg("O") == approx(1, rel=1e-3)
        assert self.ca1.get_cost_per_kg("AgO") == approx(2.7416, rel=1e-3)
        assert self.ca2.get_cost_per_kg("AgO") == approx(1.5, rel=1e-3)

    def test_cost_per_mol(self):
        assert self.ca1.get_cost_per_mol("Ag") == approx(0.3236, rel=1e-3)
        assert self.ca1.get_cost_per_mol("O") == approx(0.0160, rel=1e-3)
        assert self.ca1.get_cost_per_mol("AgO") == approx(0.3396, rel=1e-3)
        assert self.ca2.get_cost_per_mol("AgO") == approx(0.1858, rel=1e-3)

    def test_sanity(self):
        assert self.ca1.get_cost_per_kg("Ag") == self.ca2.get_cost_per_kg("Ag")


class TestCostDB(unittest.TestCase):
    def test_sanity(self):
        ca = CostAnalyzer(CostDBElements())
        assert ca.get_cost_per_kg("PtO") > ca.get_cost_per_kg("MgO")
