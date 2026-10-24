from __future__ import annotations

import os
import re
from itertools import starmap

import pytest
from monty.serialization import dumpfn, loadfn

from pymatgen.analysis.compatibility.entry_tools import (
    EntrySet,
    group_entries_by_composition,
    group_entries_by_structure,
)
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core import Composition, Element, Lattice, Structure
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/entries"


pytestmark = pytest.mark.skipif(not os.path.exists(TEST_DIR), reason="Requires test files")


class TestFunc(MatSciTest):
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
                [
                    ("Na", -2),
                    ("Na", -5),
                    ("Cl", -1),
                    ("Cl", -10),
                    ("NaCl", -20),
                    ("NaCl", -21),
                    ("Na2Cl2", -50),
                ],
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

    def test_group_entries_by_composition_unsorted(self):
        entries = [ComputedEntry("NaCl", -20), ComputedEntry("NaCl", -21), ComputedEntry("Na", -2)]
        groups = group_entries_by_composition(entries, sort_by_e_per_atom=False)
        assert [[e.energy for e in g] for g in groups] == [[-2], [-20, -21]]

    def test_group_entries_by_structure_species_to_remove(self):
        def rocksalt(cation, a):
            struct = Structure.from_spacegroup("Fm-3m", Lattice.cubic(a), [cation, "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
            return ComputedStructureEntry(struct, -10.0)

        entries = [rocksalt("Na", 5.64), rocksalt("Li", 5.13), rocksalt("K", 6.29)]
        assert len(group_entries_by_structure(entries)) == 3
        # with the cation removed, all three Cl fcc frameworks match
        groups = group_entries_by_structure(entries, species_to_remove=["Na", "Li", "K"])
        assert len(groups) == 1
        assert sorted(e.reduced_formula for e in groups[0]) == ["KCl", "LiCl", "NaCl"]
        # the original entries are untouched
        assert all(len(e.structure) == 8 for e in entries)

    def test_group_entries_by_structure_parallel(self):
        entries = loadfn(f"{TEST_DIR}/TiO2_entries.json")
        serial = group_entries_by_structure(entries)
        parallel = group_entries_by_structure(entries, ncpus=2)
        assert sorted(len(g) for g in parallel) == [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4]

        def as_ids(groups):
            return sorted(tuple(sorted(e.entry_id for e in g)) for g in groups)

        assert as_ids(parallel) == as_ids(serial)


class TestEntrySet(MatSciTest):
    def setup_method(self):
        entries = loadfn(f"{TEST_DIR}/Li-Fe-P-O_entries.json")
        self.entry_set = EntrySet(entries)

    def test_chemsys(self):
        assert self.entry_set.chemsys == {"Fe", "Li", "O", "P"}

    def test_get_subset(self):
        entries = self.entry_set.get_subset_in_chemsys(["Li", "O"])
        for ent in entries:
            assert {Element.Li, Element.O}.issuperset(ent.composition)
        with pytest.raises(
            ValueError,
            match=re.escape("['F', 'Fe'] is not a subset of ['Fe', 'Li', 'O', 'P'], extra: {'F'}"),
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

    def test_is_ground_state(self):
        ground_states = self.entry_set.ground_states
        non_ground_states = set(self.entry_set) - ground_states
        assert non_ground_states
        # is_ground_state recomputes ground_states on every call, so only sample a few entries
        for entry in sorted(ground_states, key=str)[:5]:
            assert self.entry_set.is_ground_state(entry)
        for entry in sorted(non_ground_states, key=str)[:5]:
            assert not self.entry_set.is_ground_state(entry)

    def test_mutable_set_ops(self):
        entry = ComputedEntry("Li2O", -14.0)
        length = len(self.entry_set)
        assert entry not in self.entry_set

        self.entry_set.add(entry)
        assert entry in self.entry_set
        assert len(self.entry_set) == length + 1
        self.entry_set.add(entry)
        assert len(self.entry_set) == length + 1

        self.entry_set.discard(entry)
        assert entry not in self.entry_set
        assert len(self.entry_set) == length
        # discarding an absent entry is a no-op
        self.entry_set.discard(entry)
        assert len(self.entry_set) == length

    @pytest.mark.parametrize(
        ("latexify_names", "expected_names"),
        [(False, {"Li2O", "LiFePO4", "Fe"}), (True, {"Li_{2}O", "LiFePO_{4}", "Fe"})],
    )
    def test_csv_round_trip(self, latexify_names, expected_names):
        entries = [
            PDEntry(Composition("Li2O"), -14.3),
            PDEntry(Composition("LiFePO4"), -53.5),
            PDEntry(Composition("Fe"), -8.4),
        ]
        path = f"{self.tmp_path}/entries.csv"
        EntrySet(entries).to_csv(path, latexify_names=latexify_names)

        with open(path, encoding="utf-8") as file:
            header = file.readline().strip()
        # elements sorted by electronegativity
        assert header == "Name,Li,Fe,P,O,Energy"

        loaded = EntrySet.from_csv(path)
        assert len(loaded) == 3
        assert {e.name for e in loaded} == expected_names
        energies = {e.composition.reduced_formula: e.energy for e in loaded}
        assert energies == {"Li2O": -14.3, "LiFePO4": -53.5, "Fe": -8.4}
        lfp = next(e for e in loaded if e.composition.reduced_formula == "LiFePO4")
        assert lfp.composition == Composition("LiFePO4")
        assert loaded.chemsys == {"Li", "Fe", "P", "O"}
