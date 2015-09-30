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

from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.groups import PointGroup, SpaceGroup


class PointGroupTest(unittest.TestCase):

    def test_order(self):
        order = {"mmm": 8, "432": 24, "-6m2": 12}
        for k, v in order.items():
            pg = PointGroup(k)
            self.assertEqual(order[k], len(pg.symmetry_ops))

    def test_get_orbit(self):
        pg = PointGroup("mmm")
        self.assertEqual(len(pg.get_orbit([0.1, 0.1, 0.1])), 8)
        self.assertEqual(len(pg.get_orbit([0, 0, 0.1])), 2)
        self.assertEqual(len(pg.get_orbit([1.2, 1.2, 1])), 8)

    def test_is_sub_super_group(self):
        pgmmm = PointGroup("mmm")
        pgmm2 = PointGroup("mm2")
        pg222 = PointGroup("222")
        pg4 = PointGroup("4")
        self.assertTrue(pgmmm.is_supergroup(pgmm2))
        self.assertTrue(pgmm2.is_subgroup(pgmmm))
        self.assertTrue(pgmmm.is_supergroup(pg222))
        self.assertFalse(pgmmm.is_supergroup(pg4))
        pgm3m = PointGroup("m-3m")
        pg6mmm = PointGroup("6/mmm")
        pg3m = PointGroup("-3m")
        # TODO: Fix the test below.
        #self.assertTrue(pg3m.is_subgroup(pgm3m))
        self.assertTrue(pg3m.is_subgroup(pg6mmm))
        self.assertFalse(pgm3m.is_supergroup(pg6mmm))


class SpaceGroupTest(unittest.TestCase):

    def test_abbrev_symbols(self):
        sg = SpaceGroup("P2/c")
        self.assertEqual(sg.int_number, 13)
        sg = SpaceGroup("R-3mH")
        self.assertEqual(sg.int_number, 166)

    def test_attr(self):

        sg = SpaceGroup("Fm-3m")
        self.assertEqual(sg.full_symbol, "F4/m-32/m")
        self.assertEqual(sg.patterson_symmetry, "Fm-3m")
        self.assertEqual(sg.point_group, "m-3m")


    def test_full_symbols(self):
        sg = SpaceGroup("P2/m2/m2/m")
        self.assertEqual(sg.symbol, "Pmmm")

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

    def test_is_compatible(self):
        cubic = Lattice.cubic(1)
        hexagonal = Lattice.hexagonal(1, 2)
        rhom = Lattice.rhombohedral(3, 80)
        tet = Lattice.tetragonal(1, 2)
        ortho = Lattice.orthorhombic(1, 2, 3)
        sg = SpaceGroup("Fm-3m")
        self.assertTrue(sg.is_compatible(cubic))
        self.assertFalse(sg.is_compatible(hexagonal))
        sg = SpaceGroup("R-3mH")
        self.assertFalse(sg.is_compatible(cubic))
        self.assertTrue(sg.is_compatible(hexagonal))
        sg = SpaceGroup("R-3m")
        self.assertTrue(sg.is_compatible(cubic))
        self.assertTrue(sg.is_compatible(rhom))
        self.assertFalse(sg.is_compatible(hexagonal))
        sg = SpaceGroup("Pnma")
        self.assertTrue(sg.is_compatible(cubic))
        self.assertTrue(sg.is_compatible(tet))
        self.assertTrue(sg.is_compatible(ortho))
        self.assertFalse(sg.is_compatible(rhom))
        self.assertFalse(sg.is_compatible(hexagonal))
        sg = SpaceGroup("P12/c1")
        self.assertTrue(sg.is_compatible(cubic))
        self.assertTrue(sg.is_compatible(tet))
        self.assertTrue(sg.is_compatible(ortho))
        self.assertFalse(sg.is_compatible(rhom))
        self.assertFalse(sg.is_compatible(hexagonal))
        sg = SpaceGroup("P-1")
        self.assertTrue(sg.is_compatible(cubic))
        self.assertTrue(sg.is_compatible(tet))
        self.assertTrue(sg.is_compatible(ortho))
        self.assertTrue(sg.is_compatible(rhom))
        self.assertTrue(sg.is_compatible(hexagonal))

if __name__ == '__main__':
    unittest.main()
