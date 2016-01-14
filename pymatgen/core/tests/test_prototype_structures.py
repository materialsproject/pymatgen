# coding: utf-8

from __future__ import division, print_function, unicode_literals

import unittest
#import os
from math import sqrt
import numpy

from pymatgen.core.prototype_structures import PrototypeStructure
#from pymatgen import Element, Structure, Lattice
from pymatgen.util.testing import PymatgenTest


class PrototypeStructureTest(PymatgenTest):

    def test_init(self):
        # Each newly implemented prototype should be added here to assure
        # flawless construction.
        self.assertIsNotNone(
                PrototypeStructure("cubic", supercell_scaling=[2,2,2]),
                "Initialization of \"cubic\" prototype cell failed.")
        self.assertIsNotNone(
                PrototypeStructure("primitive cubic"),
                "Initialization of \"primitive cubic\" prototype cell failed.")
        self.assertIsNotNone(
                PrototypeStructure("pc"),
                "Initialization of \"pc\" prototype cell failed.")
        self.assertIsNotNone(
                PrototypeStructure("c"),
                "Initialization of \"c\" prototype cell failed.")
        self.assertIsNotNone(
                PrototypeStructure("body-centered cubic"),
                "Initialization of \"body-centered cubic\" prototype cell" \
                " failed.")
        self.assertIsNotNone(
                PrototypeStructure("body centered cubic"),
                "Initialization of \"body centered cubic\" prototype cell" \
                " failed.")
        self.assertIsNotNone(
                PrototypeStructure("bcc"),
                "Initialization of \"bcc\" prototype cell failed.")
        self.assertIsNotNone(
                PrototypeStructure("face-centered cubic"),
                "Initialization of \"face-centered cubic\" prototype cell"
                " failed.")
        self.assertIsNotNone(
                PrototypeStructure("face centered cubic"),
                "Initialization of \"face centered cubic\" prototype cell"
                " failed.")
        self.assertIsNotNone(
                PrototypeStructure("fcc"),
                "Initialization of \"fcc\" prototype cell failed.")
        self.assertIsNotNone(
                PrototypeStructure("hexagonal closed packed"),
                "Initialization of \"hexagonal closed packed\" prototype"
                " cell failed.")
        self.assertIsNotNone(
                PrototypeStructure("hcp"),
                "Initialization of \"hcp\" prototype cell failed.")
        self.assertIsNotNone(
                PrototypeStructure("diamond"),
                "Initialization of \"diamond\" prototype cell failed.")
        self.assertIsNotNone(
                PrototypeStructure("d"),
                "Initialization of \"d\" prototype cell failed.")
        self.assertIsNotNone(
                PrototypeStructure("c", supercell_scaling=[2, 2, 2]),
                "Initialization with explicitly specifying supercell"
                " scaling values failed.")
        self.assertIsNotNone(
                PrototypeStructure("c", enveloping_box_factor=2.0),
                "Initialization with enveloping box failed.")
        cubic_envel2 = PrototypeStructure("c", enveloping_box_factor=2.0)
        self.assertAlmostEqual(cubic_envel2.lattice.a, 2.0) 
        self.assertAlmostEqual(cubic_envel2.lattice.b, 2.0)
        self.assertAlmostEqual(cubic_envel2.lattice.c, 2.0)
        self.assertAlmostEqual(cubic_envel2.sites[0].coords[0], 0.5)
        self.assertAlmostEqual(cubic_envel2.sites[0].coords[1], 0.5)
        self.assertAlmostEqual(cubic_envel2.sites[0].coords[2], 0.5)
        ortho_envel2 = PrototypeStructure(
                "c", lengths=[3.0, 2.0, 1.0], enveloping_box_factor=2.0)
        self.assertAlmostEqual(ortho_envel2.lattice.a, 6.0)
        self.assertAlmostEqual(ortho_envel2.lattice.b, 6.0)
        self.assertAlmostEqual(ortho_envel2.lattice.c, 6.0)
        self.assertAlmostEqual(ortho_envel2.sites[0].coords[0], 1.5)
        self.assertAlmostEqual(ortho_envel2.sites[0].coords[1], 2.0)
        self.assertAlmostEqual(ortho_envel2.sites[0].coords[2], 2.5)
        hcp_envel2 = PrototypeStructure(
                "hcp", enveloping_box_factor=2.0)
        self.assertAlmostEqual(hcp_envel2.lattice.a, 1.633*2.0)
        self.assertAlmostEqual(hcp_envel2.lattice.b, 1.633*2.0)
        self.assertAlmostEqual(hcp_envel2.lattice.c, 1.633*2.0)
        hcp_envel2_c1 = PrototypeStructure(
                "hcp", lengths=[1.0, 1.0, 1.0], enveloping_box_factor=2.0)
        self.assertAlmostEqual(hcp_envel2_c1.lattice.a, 1.5*2.0)
        self.assertAlmostEqual(hcp_envel2_c1.lattice.b, 1.5*2.0)
        self.assertAlmostEqual(hcp_envel2_c1.lattice.c, 1.5*2.0)
        self.assertNotAlmostEqual(hcp_envel2_c1.lattice.c, 1.51*2.0)
        hex_envel2 = PrototypeStructure(
                "c", lengths=[1.0, 1.0, 1.633], angles=[90.0, 90.0, 120.0],
                enveloping_box_factor=2.0)
        self.assertAlmostEqual(
                hex_envel2.sites[0].coords[0], 0.5-0.75+1.633)
        self.assertAlmostEqual(
                round(hex_envel2.sites[0].coords[1] - \
                      (-0.8660331248136633/2.0+1.633), 4), 0.0)
        self.assertAlmostEqual(
                hex_envel2.sites[0].coords[2], 0.5*1.633)
        with self.assertRaises(ValueError):
            PrototypeStructure("nonsense")
        with self.assertRaises(ValueError):
            PrototypeStructure("c", supercell_scaling=[0, 1, 1])


    def test_rotate_sites_around_center_of_enveloping_cuboid(self):
        cubic_envel2 = PrototypeStructure("c", enveloping_box_factor=2.0)
        self.assertIsNotNone(
                cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                        [0.0, 0.0, 1.0], 90.0))
        rot_struct = cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                [0.0, 0.0, 1.0], 90.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[0], 1.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[1], 0.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[2], 0.5)
        rot_struct = cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                [0.0, 0.0, 1.0], 180.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[0], 1.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[1], 1.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[2], 0.5)
        rot_struct = cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                [0.0, 0.0, 1.0], 270.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[0], 0.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[1], 1.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[2], 0.5)
        rot_struct = cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                [0.0, 0.0, 1.0], -90.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[0], 0.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[1], 1.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[2], 0.5)
        rot_struct = cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                [0.0, 0.0, 1.0], 45.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[0], 1.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[1], 0.5-(sqrt(0.5)-0.5))
        self.assertAlmostEqual(rot_struct.sites[0].coords[2], 0.5)
        rot_struct = cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                [0.0, 1.0, 0.0], 90.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[0], 0.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[1], 0.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[2], 1.5)
        rot_struct = cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                [1.0, 0.0, 0.0], 90.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[0], 0.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[1], 1.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[2], 0.5)
        rot_struct = cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                [1.0, 1.0, 1.0], 90.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[0], 0.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[1], 0.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[2], 0.5)
        rot_struct = cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                [1.0, 1.0, 1.0], 90.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[0], 0.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[1], 0.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[2], 0.5)
        cubic_envel2 = PrototypeStructure("c", supercell_scaling=[3, 3, 3], \
                enveloping_box_factor=2.0)
        self.assertIsNotNone(
                cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                        [0.0, 0.0, 1.0], 90.0))
        rot_struct = cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                [0.0, 0.0, 1.0], 90.0)
        self.assertAlmostEqual(rot_struct.sites[0].coords[0], 4.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[1], 1.5)
        self.assertAlmostEqual(rot_struct.sites[0].coords[2], 1.5)

        with self.assertRaises(TypeError):
            self.assertIsNotNone(
                cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                        [0.0, 0.0], 90.0))
        with self.assertRaises(ValueError):
            self.assertIsNotNone(
                cubic_envel2.rotate_sites_around_center_of_enveloping_cuboid(
                        [0.0, 0.0, 0.0], 90.0))


if __name__ == '__main__':
    unittest.main()
