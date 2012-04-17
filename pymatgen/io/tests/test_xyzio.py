#!/usr/bin/env python

'''
Created on Apr 17, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 17, 2012"

import unittest

from pymatgen.core.structure import Molecule
from pymatgen.io.xyzio import XYZ

class XYZTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
              [0.000000, 0.000000, 1.089000],
              [1.026719, 0.000000, -0.363000],
              [-0.513360, -0.889165, -0.363000],
              [-0.513360 , 0.889165 , -0.363000]]
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)
        self.xyz = XYZ(self.mol)

    def test_str(self):
        ans = """5
H4 C1
C 0.000000 0.000000 0.000000
H 0.000000 0.000000 1.089000
H 1.026719 0.000000 -0.363000
H -0.513360 -0.889165 -0.363000
H -0.513360 0.889165 -0.363000"""
        self.assertEqual(str(self.xyz), ans)

    def test_from_string(self):
        ans = """5
H4 C1
C 0.000000 0.000000 0.000000
H 0.000000 0.000000 1.089000
H 1.026719 0.000000 -0.363000
H -0.513360 -0.889165 -0.363000
H -0.513360 0.889165 -0.363000"""
        xyz = XYZ.from_string(ans)
        mol = xyz.molecule
        sp = ["C", "H", "H", "H", "H"]
        for i, site in enumerate(mol):
            self.assertEqual(site.specie.symbol, sp[i])


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
