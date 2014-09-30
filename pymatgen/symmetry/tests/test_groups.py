#!/usr/bin/env python

"""
TODO: Modify unittest doc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "4/10/14"

import unittest
import numpy as np

from pymatgen.symmetry.groups import PointGroup, SpaceGroup


class PointGroupTest(unittest.TestCase):

    def test_order(self):
        order = {"mmm": 8, "432": 24, "-6m2": 12}
        for k, v in order.items():
            pg = PointGroup(k)
            self.assertEqual(order[k], len(pg.symmetry_ops))


class SpaceGroupTest(unittest.TestCase):

    def test_order_symm_ops(self):
        for name in SpaceGroup.SG_SYMBOLS:
            sg = SpaceGroup(name)
            self.assertEqual(len(sg.symmetry_ops), sg.order)

    def test_crystal_system(self):
        sg = SpaceGroup("R-3c")
        self.assertEqual(sg.crystal_system, "trigonal")
        sg = SpaceGroup("R-3cH")
        self.assertEqual(sg.crystal_system, "trigonal")

    def test_get_orbit(self):
        sg = SpaceGroup("Fm-3m")
        p = np.random.random_integers(0, 100, size=(3,))
        p /= 100
        self.assertLessEqual(len(sg.get_orbit(p)), sg.order)

if __name__ == '__main__':
    unittest.main()
