# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest
from pathlib import Path

from monty.serialization import dumpfn, loadfn

from pymatgen.core.periodic_table import Element
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.entries.entry_tools import (
    EntrySet,
    group_entries_by_composition,
    group_entries_by_structure,
)
from pymatgen.util.testing import PymatgenTest

test_dir = Path(__file__).absolute().parent / ".." / ".." / ".." / "test_files"


class FuncTest(unittest.TestCase):
    def test_group_entries_by_structure(self):
        entries = loadfn(os.path.join(PymatgenTest.TEST_FILES_DIR, "TiO2_entries.json"))
        groups = group_entries_by_structure(entries)
        self.assertEqual(sorted(len(g) for g in groups), [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4])
        self.assertLess(len(groups), len(entries))
        # Make sure no entries are left behind
        self.assertEqual(sum(len(g) for g in groups), len(entries))

    def test_group_entries_by_composition(self):
        entries = [
            ComputedEntry("Na", -2),
            ComputedEntry("Na", -5),
            ComputedEntry("Cl", -1),
            ComputedEntry("Cl", -10),
            ComputedEntry("NaCl", -20),
            ComputedEntry("NaCl", -21),
            ComputedEntry("Na2Cl2", -50),
        ]

        groups = group_entries_by_composition(entries)
        self.assertEqual(sorted(len(g) for g in groups), [2, 2, 3])
        self.assertLess(len(groups), len(entries))
        # Make sure no entries are left behind
        self.assertEqual(sum(len(g) for g in groups), len(entries))
        # test sorting by energy
        for g in groups:
            assert g == sorted(g, key=lambda e: e.energy_per_atom)


class EntrySetTest(unittest.TestCase):
    def setUp(self):
        entries = loadfn(os.path.join(PymatgenTest.TEST_FILES_DIR, "Li-Fe-P-O_entries.json"))
        self.entry_set = EntrySet(entries)

    def test_chemsys(self):
        self.assertEqual(self.entry_set.chemsys, {"Fe", "Li", "O", "P"})

    def test_get_subset(self):
        entries = self.entry_set.get_subset_in_chemsys(["Li", "O"])
        for e in entries:
            self.assertTrue({Element.Li, Element.O}.issuperset(e.composition.keys()))
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
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
