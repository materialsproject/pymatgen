# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import unittest

from pymatgen.core.bonds import CovalentBond, get_bond_length
from pymatgen.core.sites import Site


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 26, 2012"


class CovalentBondTest(unittest.TestCase):

    def test_length(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0.7, 0.6])
        self.assertAlmostEqual(CovalentBond(site1, site2).length,
                               0.92195444572928864)

    def test_is_bonded(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0, 1])
        self.assertTrue(CovalentBond.is_bonded(site1, site2))
        site2 = Site("H", [0, 0, 1.5])
        self.assertFalse(CovalentBond.is_bonded(site1, site2))

    def test_str(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0.7, 0.6])
        self.assertIsNotNone(CovalentBond(site1, site2))


class FuncTest(unittest.TestCase):

    def test_get_bond_length(self):
        self.assertAlmostEqual(get_bond_length("C", "C", 1), 1.54)
        self.assertAlmostEqual(get_bond_length("C", "C", 2), 1.34)
        self.assertAlmostEqual(get_bond_length("C", "H", 1), 1.08)
        self.assertEqual(get_bond_length("C", "H", 2), None)
        self.assertAlmostEqual(get_bond_length("C", "Br", 1), 1.85)

if __name__ == "__main__":
    unittest.main()
