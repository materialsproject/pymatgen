#!/usr/bin/env python

'''
Created on Mar 18, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 18, 2012"

import unittest
import os

from pymatgen.borg.hive import VaspToComputedEntryDrone
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.entries.compatibility import MITCompatibility
import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')

class VaspToComputedEntryDroneTest(unittest.TestCase):

    def setUp(self):
        self.drone = VaspToComputedEntryDrone(data=["efermi"])
        self.structure_drone = VaspToComputedEntryDrone(True)

    def test_get_valid_paths(self):
        for path in os.walk(test_dir):
            self.assertTrue(len(self.drone.get_valid_paths(path)) > 0)

    def test_assimilate_and_convert(self):
        d = self.drone.assimilate(test_dir)
        for p in ["hubbards", "is_hubbard", "potcar_symbols", "run_type"]:
            self.assertIn(p, d["parameters"])
        self.assertAlmostEqual(d["data"]["efermi"], 1.8301027)
        entry = self.drone.convert(d)
        self.assertEqual(entry.composition.reduced_formula, "LiFe4(PO4)4")
        self.assertAlmostEqual(entry.energy, -269.38319884)
        d = self.structure_drone.assimilate(test_dir)
        entry = self.structure_drone.convert(d)
        self.assertEqual(entry.composition.reduced_formula, "LiFe4(PO4)4")
        self.assertAlmostEqual(entry.energy, -269.38319884)
        self.assertIsInstance(entry, ComputedStructureEntry)
        self.assertIsNotNone(entry.structure)
        compat = MITCompatibility()
        self.assertIsNone(compat.process_entry(entry))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
