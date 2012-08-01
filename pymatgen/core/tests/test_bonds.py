#!/usr/bin/env python

'''
Created on Jul 26, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 26, 2012"

import unittest

from pymatgen.core.bonds import CovalentBond
from pymatgen.core.sites import Site


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

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
