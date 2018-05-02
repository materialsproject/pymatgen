# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import unittest
import warnings
from pymatgen.core.bonds import CovalentBond, get_bond_length, \
    obtain_all_bond_lengths, get_bond_order
from pymatgen.core.sites import Site
from pymatgen.core.periodic_table import Element


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 26, 2012"


class CovalentBondTest(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_length(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0.7, 0.6])
        self.assertAlmostEqual(CovalentBond(site1, site2).length,
                               0.92195444572928864)

    def test_get_bond_order(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0, 1.08])
        self.assertAlmostEqual(
            CovalentBond(site1, site2).get_bond_order(), 1)
        bond = CovalentBond(Site("C", [0, 0, 0]), Site("Br", [0, 0, 2]))
        self.assertAlmostEqual(
            bond.get_bond_order(0.5, 1.9), 0.894736842105263)

    def test_is_bonded(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0, 1])
        self.assertTrue(CovalentBond.is_bonded(site1, site2))
        site2 = Site("H", [0, 0, 1.5])
        self.assertFalse(CovalentBond.is_bonded(site1, site2))
        site1 = Site("U", [0, 0, 0])
        self.assertRaises(ValueError, CovalentBond.is_bonded, site1, site2)
        self.assertTrue(CovalentBond.is_bonded(site1, site2, default_bl=2))

    def test_str(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0.7, 0.6])
        self.assertIsNotNone(CovalentBond(site1, site2))


class FuncTest(unittest.TestCase):

    def test_get_bond_length(self):
        self.assertAlmostEqual(get_bond_length("C", "C", 1), 1.54)
        self.assertAlmostEqual(get_bond_length("C", "C", 2), 1.34)
        self.assertAlmostEqual(get_bond_length("C", "H", 1), 1.08)
        self.assertEqual(get_bond_length("C", "H", 2), 0.95)
        self.assertAlmostEqual(get_bond_length("C", "Br", 1), 1.85)

    def test_obtain_all_bond_lengths(self):
        self.assertDictEqual(obtain_all_bond_lengths('C', 'C'),
                             {1.0: 1.54, 2.0: 1.34, 3.0: 1.2})
        self.assertRaises(ValueError, obtain_all_bond_lengths,
                          'Br', Element('C'))
        self.assertDictEqual(obtain_all_bond_lengths(
            'C', Element('Br'), 1.76), {1: 1.76})
        bond_lengths_dict = obtain_all_bond_lengths('C', 'N')
        bond_lengths_dict[4] = 999
        self.assertDictEqual(obtain_all_bond_lengths('C', 'N'),
                             {1.0: 1.47, 2.0: 1.3, 3.0: 1.16})

    def test_get_bond_order(self):
        self.assertAlmostEqual(get_bond_order(
            'C', 'C', 1), 3)
        self.assertAlmostEqual(get_bond_order(
            'C', 'C', 1.2), 3)
        self.assertAlmostEqual(get_bond_order(
            'C', 'C', 1.25), 2.642857142857143)
        self.assertAlmostEqual(get_bond_order(
            'C', 'C', 1.34), 2)
        self.assertAlmostEqual(get_bond_order(
            'C', 'C', 1.4), 1.7)  # bond length in benzene
        self.assertAlmostEqual(get_bond_order(
            'C', 'C', 1.54), 1)
        self.assertAlmostEqual(get_bond_order(
            'C', 'C', 2.5), 0)
        self.assertAlmostEqual(get_bond_order(
            'C', 'C', 9999), 0)
        self.assertAlmostEqual(get_bond_order(
            'C', 'Br', 1.9, default_bl=1.9), 1)
        self.assertAlmostEqual(get_bond_order(
            'C', 'Br', 2, default_bl=1.9), 0.7368421052631575)
        self.assertAlmostEqual(get_bond_order(
            'C', 'Br', 1.9, tol=0.5, default_bl=1.9), 1)
        self.assertAlmostEqual(get_bond_order(
            'C', 'Br', 2, tol=0.5, default_bl=1.9), 0.894736842105263)
        self.assertRaises(ValueError, get_bond_order, 'C', 'Br', 1.9)


if __name__ == "__main__":
    unittest.main()
