#!/usr/bin/env python

'''
Created on Oct 24, 2012

@author: shyue
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Oct 24, 2012"

import unittest
import os

from pymatgen.io.cifio import CifParser
from pymatgen.core.periodic_table import Specie
from pymatgen.analysis.bond_valence import BVAnalyzer


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class BVAnalyzerTest(unittest.TestCase):

    def setUp(self):
        self.analyzer = BVAnalyzer()

    def test_get_valence(self):
        parser = CifParser(os.path.join(test_dir, "LiMn2O4.cif"))
        s = parser.get_structures()[0]
        ans = [1, 1, 3, 3, 4, 4, -2, -2, -2, -2, -2, -2, -2, -2]
        self.assertEqual(self.analyzer.get_valences(s), ans)
        parser = CifParser(os.path.join(test_dir, "LiFePO4.cif"))
        s = parser.get_structures()[0]
        ans = [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2, -2, -2, -2,
               - 2, -2, -2, -2, -2, -2, -2, -2, -2]
        self.assertEqual(self.analyzer.get_valences(s), ans)
        parser = CifParser(os.path.join(test_dir, "Li3V2(PO4)3.cif"))
        s = parser.get_structures()[0]
        ans = [1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, -2, -2, -2, -2,
               - 2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,
               - 2, -2, -2, -2]
        self.assertEqual(self.analyzer.get_valences(s), ans)
        parser = CifParser(os.path.join(test_dir, "Li4Fe3Mn1(PO4)4.cif"))
        s = parser.get_structures()[0]
        ans = [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2, -2, -2, -2,
               - 2, -2, -2, -2, -2, -2, -2, -2, -2]
        self.assertEqual(self.analyzer.get_valences(s), ans)
        parser = CifParser(os.path.join(test_dir, "NaFePO4.cif"))
        s = parser.get_structures()[0]
        ans = [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2, -2, -2, -2,
               - 2, -2, -2, -2, -2, -2, -2, -2, -2]
        self.assertEqual(self.analyzer.get_valences(s), ans)

    def test_get_oxi_state_structure(self):
        parser = CifParser(os.path.join(test_dir, "LiMn2O4.cif"))
        s = parser.get_structures()[0]
        news = self.analyzer.get_oxi_state_decorated_structure(s)
        self.assertIn(Specie("Mn", 3), news.composition.elements)
        self.assertIn(Specie("Mn", 4), news.composition.elements)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
