# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Specie
from pymatgen.analysis.bond_valence import BVAnalyzer, calculate_bv_sum, \
    calculate_bv_sum_unordered
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class BVAnalyzerTest(PymatgenTest):

    def setUp(self):
        self.analyzer = BVAnalyzer()

    def test_get_valence(self):
        s = Structure.from_file(os.path.join(test_dir, "LiMn2O4.json"))
        ans = [1, 1, 3, 3, 4, 4, -2, -2, -2, -2, -2, -2, -2, -2]
        self.assertEqual(self.analyzer.get_valences(s), ans)
        s = self.get_structure("LiFePO4")
        ans = [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2, -2, -2, -2,
               - 2, -2, -2, -2, -2, -2, -2, -2, -2]
        self.assertEqual(self.analyzer.get_valences(s), ans)
        s = self.get_structure("Li3V2(PO4)3")
        ans = [1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, -2, -2, -2, -2,
               - 2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,
               - 2, -2, -2, -2]
        self.assertEqual(self.analyzer.get_valences(s), ans)
        s = Structure.from_file(os.path.join(test_dir, "Li4Fe3Mn1(PO4)4.json"))
        ans = [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2, -2, -2, -2,
               - 2, -2, -2, -2, -2, -2, -2, -2, -2]
        self.assertEqual(self.analyzer.get_valences(s), ans)
        s = self.get_structure("NaFePO4")
        ans = [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2, -2, -2, -2,
               - 2, -2, -2, -2, -2, -2, -2, -2, -2]
        self.assertEqual(self.analyzer.get_valences(s), ans)

    def test_get_oxi_state_structure(self):
        s = Structure.from_file(os.path.join(test_dir, "LiMn2O4.json"))
        news = self.analyzer.get_oxi_state_decorated_structure(s)
        self.assertIn(Specie("Mn", 3), news.composition.elements)
        self.assertIn(Specie("Mn", 4), news.composition.elements)


class BondValenceSumTest(PymatgenTest):

    def test_calculate_bv_sum(self):
        s = Structure.from_file(os.path.join(test_dir, "LiMn2O4.json"))
        neighbors = s.get_neighbors(s[0], 3.0)
        bv_sum = calculate_bv_sum(s[0], neighbors)
        self.assertAlmostEqual(bv_sum, 0.7723402182087497, places=5)

    def test_calculate_bv_sum_unordered(self):
        s = Structure.from_file(os.path.join(test_dir, "LiMn2O4.json"))
        s[0].species = Composition("Li0.5Na0.5")
        neighbors = s.get_neighbors(s[0], 3.0)
        bv_sum = calculate_bv_sum_unordered(s[0], neighbors)
        self.assertAlmostEqual(bv_sum, 1.5494662306918852, places=5)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
