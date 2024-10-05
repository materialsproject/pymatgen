from __future__ import annotations

import re
from itertools import starmap

import pytest
from monty.serialization import dumpfn, loadfn

from pymatgen.core import Element
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.entries.entry_tools import EntrySet, group_entries_by_composition, group_entries_by_structure
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/entries"


class TestFunc(PymatgenTest):
    def test_group_entries_by_structure(self):
        entries = loadfn(f"{TEST_DIR}/TiO2_entries.json")
        groups = group_entries_by_structure(entries)
        assert sorted(len(g) for g in groups) == [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4]
        assert len(groups) < len(entries)
        # Make sure no entries are left behind
        assert sum(len(g) for g in groups) == len(entries)

    def test_group_entries_by_composition(self):
        entries = [
            *starmap(
                ComputedEntry,
                [("Na", -2), ("Na", -5), ("Cl", -1), ("Cl", -10), ("NaCl", -20), ("NaCl", -21), ("Na2Cl2", -50)],
            )
        ]

        groups = group_entries_by_composition(entries)
        assert sorted(len(g) for g in groups) == [2, 2, 3]
        assert len(groups) < len(entries)
        # Make sure no entries are left behind
        assert sum(len(g) for g in groups) == len(entries)
        # test sorting by energy
        for group in groups:
            assert group == sorted(group, key=lambda e: e.energy_per_atom)


class TestEntrySet(PymatgenTest):
    def setUp(self):
        entries = loadfn(f"{TEST_DIR}/Li-Fe-P-O_entries.json")
        self.entry_set = EntrySet(entries)

    def test_chemsys(self):
        assert self.entry_set.chemsys == {"Fe", "Li", "O", "P"}

    def test_get_subset(self):
        entries = self.entry_set.get_subset_in_chemsys(["Li", "O"])
        for ent in entries:
            assert {Element.Li, Element.O}.issuperset(ent.composition)
        with pytest.raises(
            ValueError, match=re.escape("['F', 'Fe'] is not a subset of ['Fe', 'Li', 'O', 'P'], extra: {'F'}")
        ):
            self.entry_set.get_subset_in_chemsys(["Fe", "F"])

    def test_remove_non_ground_states(self):
        length = len(self.entry_set)
        self.entry_set.remove_non_ground_states()
        assert len(self.entry_set) < length

    def test_as_dict(self):
        dumpfn(self.entry_set, f"{self.tmp_path}/temp_entry_set.json")
        entry_set = loadfn(f"{self.tmp_path}/temp_entry_set.json")
        assert len(entry_set) == len(self.entry_set)

    def test_ground_states(self):
        ground_states = self.entry_set.ground_states
        assert len(ground_states) < len(self.entry_set)

        # Check if ground states have the lowest energy per atom for each composition
        for gs in ground_states:
            same_comp_entries = [
                ent for ent in self.entry_set if ent.composition.reduced_formula == gs.composition.reduced_formula
            ]
            assert gs.energy_per_atom <= min(entry.energy_per_atom for entry in same_comp_entries)
