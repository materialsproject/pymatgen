# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Created on Jan 25, 2012
"""

__author__ = "Jimmy Shen"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"
__date__ = "April 1, 2019"

import unittest
import os
import json

from monty.serialization import loadfn
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.apps.battery.insertion_battery import InsertionElectrode
from pymatgen import MontyEncoder, MontyDecoder

test_dir = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "..", 'test_files')


class HopTest(unittest.TestCase):
    def setUp(self):
        self.test_ents_MOF = loadfn(test_dir + '/Mn6O5F7_cat_migration.json')
        self.aeccar_MOF = Chgcar.from_file(test_dir + '/AECCAR_Mn6O5F7.vasp')
        self.mpa_MOF = MigrationPathAnalyzer(
            base_entry=self.test_ents_MOF['ent_base'],
            single_cat_entries=self.test_ents_MOF['one_cation'],
            base_aeccar=self.aeccar_MOF)
        self.assertEqual(True, False)

    def test_sum_tuple(self):
        self.assertEqual(sum((1, 2, 2)), 6, "Should be 6")

if __name__ == '__main__':
    unittest.main()
