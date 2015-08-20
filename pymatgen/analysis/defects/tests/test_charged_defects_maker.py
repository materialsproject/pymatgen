# coding: utf-8

from __future__ import division, unicode_literals

__author__ = "Nils E. R. Zimmermann"
__copyright__ = "Copyright 2015, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Nils E. R. Zimmermann"
__email__ = "n.zimmermann@tuhh.de"
__date__ = "Aug 12, 2015"

import unittest

from pymatgen.analysis.defects.charged_defects_maker import get_optimized_sc_scale, \
        ChargedDefectsStructures
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite

class SupercellScalingTest(unittest.TestCase):

    def setUp(self):

        self.H_struct = Structure(
                Lattice.orthorhombic(10.0, 4.55, 1.5), ['H'],
                [[0.0, 0.0, 0.0]],
                validate_proximity=False, to_unit_cell=False,
                coords_are_cartesian=True, site_properties=None)
        self.H_cubic_struct = Structure(
                Lattice.cubic(10.0), ['H'],
                [[0.0, 0.0, 0.0]],
                validate_proximity=False, to_unit_cell=False,
                coords_are_cartesian=True, site_properties=None)

    def test_get_optimized_sc_scale(self):

        with self.assertRaises(ValueError):
            get_optimized_sc_scale(None, 1)
        with self.assertRaises(ValueError):
            get_optimized_sc_scale(self.H_struct, 0)
        with self.assertRaises(ValueError):
            get_optimized_sc_scale(self.H_struct, None)
        self.assertEqual(get_optimized_sc_scale(self.H_struct, 10), [1, 2, 5])
        self.assertEqual(get_optimized_sc_scale(
                self.H_cubic_struct, 8), [2, 2, 2])
        self.assertEqual(get_optimized_sc_scale(
                self.H_cubic_struct, 12), [2, 2, 2])

    def tearDown(self):

        del self.H_struct
        del self.H_cubic_struct


class ChargedDefectsStructuresTest(unittest.TestCase):


    def setUp(self):

        # A simple and lean primitive cubic structure.
        self.C_cubic_struct = Structure(
                Lattice.cubic(2.0),
                ['C'], [[0.0, 0.0, 0.0]],
                validate_proximity=False, to_unit_cell=False,
                coords_are_cartesian=False, site_properties=None)

        # Coordinates from mp-66.
        self.diamond_struct = Structure(
                Lattice.cubic(3.574593),
                ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'],
                [[0.0, 0.0, 0.5], [0.75, 0.75, 0.75],
                 [0.0, 0.5, 0.0], [0.75, 0.25, 0.25],
                 [0.5, 0.0, 0.0], [0.25, 0.75, 0.25],
                 [0.5, 0.5, 0.5], [0.25, 0.25, 0.75],],
                validate_proximity=False, to_unit_cell=False,
                coords_are_cartesian=False, site_properties=None)

        # Coordinates from mp-2534.
        self.GaAs_struct = Structure(
                Lattice.cubic(5.750183),
                ['Ga', 'Ga', 'Ga', 'Ga', 'As', 'As', 'As', 'As'],
                [[0.0, 0.0, 0.0],    [0.0, 0.5, 0.5],
                 [0.5, 0.0, 0.5],    [0.5, 0.5, 0.0],
                 [0.25, 0.75, 0.25], [0.25, 0.25, 0.75],
                 [0.75, 0.75, 0.75], [0.75, 0.25, 0.25],],
                validate_proximity=False, to_unit_cell=False,
                coords_are_cartesian=False, site_properties=None)


    def test_init(self):

        # Test constructor and its optional arguments for expected successes.
        self.assertIsNotNone(
                ChargedDefectsStructures(self.C_cubic_struct),
                "Initialization with primitive cubic structure failed")
        self.assertIsNotNone(
                ChargedDefectsStructures(self.C_cubic_struct, cellmax=2),
                "Initialization with input cellmax failed")
        self.assertIsNotNone(
                ChargedDefectsStructures(self.C_cubic_struct, cellmax=2,
                        min_max_oxi={"C": (-2, 2)}),
                "Initialization with input min_max_oxi failed (cubic C)")
        self.assertIsNotNone(
                ChargedDefectsStructures(self.C_cubic_struct, cellmax=2,
                        substitutions={"C": "N"}),
                "Initialization with subsitutions failed")
        self.assertIsNotNone(
                ChargedDefectsStructures(self.C_cubic_struct, cellmax=2,
                        oxi_states={"C": 0}),
                "Initialization with input oxi_states failed")
        self.assertIsNotNone(
                ChargedDefectsStructures(self.C_cubic_struct, cellmax=2,
                        antisites_flag=False),
                "Initialization with antisite-generation turned off failed")
        self.assertIsNotNone(
                ChargedDefectsStructures(self.C_cubic_struct, cellmax=2,
                        standardized=True),
                "Initialization with standardizing input structure failed")
        self.assertIsNotNone(
                ChargedDefectsStructures(self.C_cubic_struct, cellmax=2,
                        charge_states='conservative'),
                "Initialization with conservative charge state setup failed")
        self.assertIsNotNone(
                ChargedDefectsStructures(self.diamond_struct, cellmax=8),
                "Initialization with diamond structure failed")
        self.assertIsNotNone(
                ChargedDefectsStructures(self.GaAs_struct, cellmax=8),
                "Initialization with GaAs structure failed")
        self.assertIsNotNone(
                ChargedDefectsStructures(self.GaAs_struct, cellmax=8,
                        min_max_oxi={"Ga": (-2 ,2), "As": (-2 ,2)}),
                "Initialization with input min_max_oxi failed (GaAs)")

        # Test constructor and its optional arguments for expected failures.
        with self.assertRaises(ValueError):
            ChargedDefectsStructures(self.C_cubic_struct, cellmax=1)
        with self.assertRaises(ValueError):
            ChargedDefectsStructures(self.C_cubic_struct, cellmax=2, \
                    min_max_oxi={"C": (2, -2)})
        with self.assertRaises(ValueError):
            ChargedDefectsStructures(self.C_cubic_struct, cellmax=2,
                    charge_states='nonsense')
        with self.assertRaises(ValueError):
            ChargedDefectsStructures(self.C_cubic_struct, cellmax=2, \
                    substitutions={"N": "C"})
        with self.assertRaises(ValueError):
            ChargedDefectsStructures(self.C_cubic_struct, cellmax=2, \
                    oxi_states={"C": 0, "N": 0})
        with self.assertRaises(ValueError):
            ChargedDefectsStructures(self.C_cubic_struct, cellmax=2, \
                    oxi_states={"N": 0})
        with self.assertRaises(ValueError):
            ChargedDefectsStructures(self.C_cubic_struct, cellmax=2, \
                    oxi_states={"C": 1})



    def test_make_interstitial(self):

        inter_site = PeriodicSite(
                "C", [0.0, 0.0, 0.5], self.C_cubic_struct.lattice)
        cds_C_cubic = ChargedDefectsStructures(
                self.C_cubic_struct, cellmax=2)
        self.assertIsNotNone(
                cds_C_cubic.make_interstitial(inter_site, [1, 1, 2]),
                'Inserting an interstitial site into a primitive cubic' \
                ' structure failed')


    def tearDown(self):

        del self.C_cubic_struct
        del self.diamond_struct
        del self.GaAs_struct


if __name__ == '__main__':
    unittest.main()
