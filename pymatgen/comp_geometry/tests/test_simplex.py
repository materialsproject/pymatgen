#!/usr/bin/env python

'''
Created on May 15, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "May 15, 2012"

import unittest

from pymatgen.comp_geometry.simplex import Simplex


class SimplexTest(unittest.TestCase):

    def setUp(self):
        coords = []
        coords.append([0, 0, 0])
        coords.append([0, 1, 0])
        coords.append([0, 0, 1])
        coords.append([1, 0, 0])
        self.simplex = Simplex(coords)

    def test_in_simplex(self):
        self.assertTrue(self.simplex.in_simplex([0.1, 0.1, 0.1]))
        self.assertFalse(self.simplex.in_simplex([0.6, 0.6, 0.6]))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
