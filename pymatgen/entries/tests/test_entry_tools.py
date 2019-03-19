# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
FIXME: Proper module docstring
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Nov 9, 2012"

import unittest
from pathlib import Path
from monty.serialization import loadfn, dumpfn
import os
from pymatgen.core.periodic_table import Element
from pymatgen.entries.entry_tools import group_entries_by_structure, EntrySet

test_dir = Path(__file__).absolute().parent / ".." / ".." / ".." / 'test_files'


class FuncTest(unittest.TestCase):

    def test_group_entries_by_structure(self):

        entries = loadfn(str(test_dir / "TiO2_entries.json"))
        groups = group_entries_by_structure(entries)
        self.assertEqual(sorted([len(g) for g in groups]),
                         [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4])
        self.assertLess(len(groups), len(entries))
        # Make sure no entries are left behind
        self.assertEqual(sum([len(g) for g in groups]), len(entries))


class EntrySetTest(unittest.TestCase):

    def setUp(self):
        entries = loadfn(str(test_dir / "Li-Fe-P-O_entries.json"))
        self.entry_set = EntrySet(entries)

    def test_chemsys(self):
        self.assertEqual(self.entry_set.chemsys, {'Fe', 'Li', 'O', 'P'})

    def test_get_subset(self):
        entries = self.entry_set.get_subset_in_chemsys(["Li", "O"])
        for e in entries:
            self.assertTrue(set([Element.Li, Element.O]).issuperset(e.composition.keys()))
        self.assertRaises(ValueError, self.entry_set.get_subset_in_chemsys, ["Fe", "F"])

    def test_remove_non_ground_states(self):
        l = len(self.entry_set)
        self.entry_set.remove_non_ground_states()
        self.assertLess(len(self.entry_set), l)

    def test_as_dict(self):
        dumpfn(self.entry_set, "temp_entry_set.json")
        entry_set = loadfn("temp_entry_set.json")
        self.assertEqual(len(entry_set), len(self.entry_set))
        os.remove("temp_entry_set.json")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
