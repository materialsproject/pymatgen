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
import os

from pymatgen.core.structure import Molecule
from pymatgen.io.xyzio import XYZ
from pymatgen.io.vaspio.vasp_input import Poscar
import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)),
                        '..', 'test_files')


class XYZTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
              [0.000000, 0.000000, 1.089000],
              [1.026719, 0.000000, -0.363000],
              [-0.513360, -0.889165, -0.363000],
              [-0.513360, 0.889165, -0.363000]]
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

    def test_init_from_structure(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        struct = poscar.structure
        xyz = XYZ(struct)
        ans = """24
Fe4 P4 O16
Fe 2.277347 4.550379 2.260125
Fe 2.928536 1.516793 4.639870
Fe 7.483231 4.550379 0.119620
Fe 8.134420 1.516793 2.499364
P 0.985089 1.516793 1.990624
P 4.220794 4.550379 4.370369
P 6.190973 1.516793 0.389120
P 9.426677 4.550379 2.768865
O 0.451582 4.550379 3.365614
O 1.006219 1.516793 3.528306
O 1.725331 0.279529 1.358282
O 1.725331 2.754057 1.358282
O 3.480552 3.313115 3.738027
O 3.480552 5.787643 3.738027
O 4.199665 4.550379 1.148562
O 4.754301 1.516793 0.985870
O 5.657466 4.550379 3.773620
O 6.212102 1.516793 3.610928
O 6.931215 0.279529 1.021463
O 6.931215 2.754057 1.021463
O 8.686436 3.313115 3.401208
O 8.686436 5.787643 3.401208
O 9.405548 4.550379 1.231183
O 9.960184 1.516793 1.393875"""
        self.assertEqual(str(xyz), ans)


if __name__ == "__main__":
    unittest.main()
