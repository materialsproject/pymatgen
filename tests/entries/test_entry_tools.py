from __future__ import annotations

import os
import unittest

import pytest
from monty.serialization import dumpfn, loadfn

from pymatgen.core.periodic_table import Element
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.entries.entry_tools import EntrySet, group_entries_by_composition, group_entries_by_structure
from pymatgen.util.testing import TEST_FILES_DIR


class TestFunc(unittest.TestCase):
    def test_group_entries_by_structure(self):
        entries = loadfn(f"{TEST_FILES_DIR}/TiO2_entries.json")
        groups = group_entries_by_structure(entries)
        assert sorted(len(g) for g in groups) == [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4]
        assert len(groups) < len(entries)
        # Make sure no entries are left behind
        assert sum(len(g) for g in groups) == len(entries)

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
        assert sorted(len(g) for g in groups) == [2, 2, 3]
        assert len(groups) < len(entries)
        # Make sure no entries are left behind
        assert sum(len(g) for g in groups) == len(entries)
        # test sorting by energy
        for g in groups:
            assert g == sorted(g, key=lambda e: e.energy_per_atom)


class TestEntrySet(unittest.TestCase):
    def setUp(self):
        entries = loadfn(f"{TEST_FILES_DIR}/Li-Fe-P-O_entries.json")
        self.entry_set = EntrySet(entries)

    def test_chemsys(self):
        assert self.entry_set.chemsys == {"Fe", "Li", "O", "P"}

    def test_get_subset(self):
        entries = self.entry_set.get_subset_in_chemsys(["Li", "O"])
        for e in entries:
            assert {Element.Li, Element.O}.issuperset(e.composition)
        with pytest.raises(ValueError) as exc:  # noqa: PT011
            self.entry_set.get_subset_in_chemsys(["Fe", "F"])
        assert "['F', 'Fe'] is not a subset of ['Fe', 'Li', 'O', 'P'], extra: {'F'}" in str(exc.value)

    def test_remove_non_ground_states(self):
        length = len(self.entry_set)
        self.entry_set.remove_non_ground_states()
        assert len(self.entry_set) < length

    def test_as_dict(self):
        dumpfn(self.entry_set, "temp_entry_set.json")
        entry_set = loadfn("temp_entry_set.json")
        assert len(entry_set) == len(self.entry_set)
        os.remove("temp_entry_set.json")
