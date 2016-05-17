# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on Mar 18, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 18, 2012"

import unittest2 as unittest
import os

from pymatgen.apps.borg.hive import VaspToComputedEntryDrone, \
    SimpleVaspToComputedEntryDrone, GaussianToComputedEntryDrone
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.entries.compatibility import MITCompatibility


class VaspToComputedEntryDroneTest(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__), "..", "..",
                                     "..", "..", 'test_files')
        self.drone = VaspToComputedEntryDrone(data=["efermi"])
        self.structure_drone = VaspToComputedEntryDrone(True)

    def test_get_valid_paths(self):
        for path in os.walk(self.test_dir):
            if path[0] == self.test_dir:
                self.assertTrue(len(self.drone.get_valid_paths(path)) > 0)

    def test_assimilate(self):
        entry = self.drone.assimilate(self.test_dir)
        for p in ["hubbards", "is_hubbard", "potcar_spec", "run_type"]:
            self.assertIn(p, entry.parameters)
        self.assertAlmostEqual(entry.data["efermi"], 1.8301027)
        self.assertEqual(entry.composition.reduced_formula, "LiFe4(PO4)4")
        self.assertAlmostEqual(entry.energy, -269.38319884)
        entry = self.structure_drone.assimilate(self.test_dir)
        self.assertEqual(entry.composition.reduced_formula, "LiFe4(PO4)4")
        self.assertAlmostEqual(entry.energy, -269.38319884)
        self.assertIsInstance(entry, ComputedStructureEntry)
        self.assertIsNotNone(entry.structure)
        self.assertEqual(len(entry.parameters["history"]), 2)
        compat = MITCompatibility(check_potcar_hash=False)
        self.assertIsNone(compat.process_entry(entry))

    def test_to_from_dict(self):
        d = self.structure_drone.as_dict()
        drone = VaspToComputedEntryDrone.from_dict(d)
        self.assertEqual(type(drone), VaspToComputedEntryDrone)


class SimpleVaspToComputedEntryDroneTest(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__), "..", "..",
                                     "..", "..", 'test_files')
        self.drone = SimpleVaspToComputedEntryDrone()
        self.structure_drone = SimpleVaspToComputedEntryDrone(True)

    def test_get_valid_paths(self):
        for path in os.walk(self.test_dir):
            if path[0] == self.test_dir:
                self.assertTrue(len(self.drone.get_valid_paths(path)) > 0)

    def test_to_from_dict(self):
        d = self.structure_drone.as_dict()
        drone = SimpleVaspToComputedEntryDrone.from_dict(d)
        self.assertEqual(type(drone), SimpleVaspToComputedEntryDrone)


class GaussianToComputedEntryDroneTest(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__), "..", "..",
                                     "..", "..", 'test_files', "molecules")
        self.drone = GaussianToComputedEntryDrone(data=["corrections"])
        self.structure_drone = GaussianToComputedEntryDrone(True)

    def test_get_valid_paths(self):
        for path in os.walk(self.test_dir):
            if path[0] == self.test_dir:
                self.assertTrue(len(self.drone.get_valid_paths(path)) > 0)

    def test_assimilate(self):
        test_file = os.path.join(self.test_dir, "methane.log")
        entry = self.drone.assimilate(test_file)
        for p in ["functional", "basis_set", "charge", "spin_mult", 'route']:
            self.assertIn(p, entry.parameters)
        for p in ["corrections"]:
            self.assertIn(p, entry.data)

        self.assertEqual(entry.composition.reduced_formula, "H4C")
        self.assertAlmostEqual(entry.energy, -39.9768775602)
        entry = self.structure_drone.assimilate(test_file)
        self.assertEqual(entry.composition.reduced_formula, "H4C")
        self.assertAlmostEqual(entry.energy, -39.9768775602)
        self.assertIsInstance(entry, ComputedStructureEntry)
        self.assertIsNotNone(entry.structure)
        for p in ["properly_terminated", "stationary_type"]:
            self.assertIn(p, entry.data)

    def test_to_from_dict(self):
        d = self.structure_drone.as_dict()
        drone = GaussianToComputedEntryDrone.from_dict(d)
        self.assertEqual(type(drone), GaussianToComputedEntryDrone)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
