from __future__ import annotations

import json
import unittest

from monty.json import MontyDecoder

from pymatgen.alchemy.filters import (
    ContainsSpecieFilter,
    RemoveDuplicatesFilter,
    RemoveExistingFilter,
    SpecieProximityFilter,
)
from pymatgen.alchemy.transmuters import StandardTransmuter
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Species
from pymatgen.core.structure import Structure
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestContainsSpecieFilter(PymatgenTest):
    def test_filtering(self):
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = Lattice([[3.0, 0.0, 0.0], [1.0, 3.0, 0], [0, -2.0, 3.0]])
        struct = Structure(lattice, [{"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25}] * 4, coords)

        species1 = [Species("Si", 5), Species("Mg", 2)]
        f1 = ContainsSpecieFilter(species1, strict_compare=True, AND=False)
        assert not f1.test(struct), "Incorrect filter"
        f2 = ContainsSpecieFilter(species1, strict_compare=False, AND=False)
        assert f2.test(struct), "Incorrect filter"
        species2 = [Species("Si", 4), Species("Mg", 2)]
        f3 = ContainsSpecieFilter(species2, strict_compare=True, AND=False)
        assert f3.test(struct), "Incorrect filter"
        f4 = ContainsSpecieFilter(species2, strict_compare=False, AND=False)
        assert f4.test(struct), "Incorrect filter"

        species3 = [Species("Si", 5), Species("O", -2)]
        f5 = ContainsSpecieFilter(species3, strict_compare=True, AND=True)
        assert not f5.test(struct), "Incorrect filter"
        f6 = ContainsSpecieFilter(species3, strict_compare=False, AND=True)
        assert f6.test(struct), "Incorrect filter"
        species4 = [Species("Si", 4), Species("Mg", 2)]
        f7 = ContainsSpecieFilter(species4, strict_compare=True, AND=True)
        assert not f7.test(struct), "Incorrect filter"
        f8 = ContainsSpecieFilter(species4, strict_compare=False, AND=True)
        assert not f8.test(struct), "Incorrect filter"

    def test_to_from_dict(self):
        species1 = ["Si5+", "Mg2+"]
        f1 = ContainsSpecieFilter(species1, strict_compare=True, AND=False)
        d = f1.as_dict()
        assert isinstance(ContainsSpecieFilter.from_dict(d), ContainsSpecieFilter)


class TestSpecieProximityFilter(PymatgenTest):
    def test_filter(self):
        struct = self.get_structure("Li10GeP2S12")
        sf = SpecieProximityFilter({"Li": 1})
        assert sf.test(struct)
        sf = SpecieProximityFilter({"Li": 2})
        assert not sf.test(struct)
        sf = SpecieProximityFilter({"P": 1})
        assert sf.test(struct)
        sf = SpecieProximityFilter({"P": 5})
        assert not sf.test(struct)

    def test_to_from_dict(self):
        sf = SpecieProximityFilter({"Li": 1})
        d = sf.as_dict()
        assert isinstance(SpecieProximityFilter.from_dict(d), SpecieProximityFilter)


class TestRemoveDuplicatesFilter(unittest.TestCase):
    def setUp(self):
        with open(f"{TEST_FILES_DIR}/TiO2_entries.json") as file:
            entries = json.load(file, cls=MontyDecoder)
        self._struct_list = [e.structure for e in entries]
        self._sm = StructureMatcher()

    def test_filter(self):
        transmuter = StandardTransmuter.from_structures(self._struct_list)
        dup_filter = RemoveDuplicatesFilter()
        transmuter.apply_filter(dup_filter)
        assert len(transmuter.transformed_structures) == 11

    def test_to_from_dict(self):
        fil = RemoveDuplicatesFilter()
        d = fil.as_dict()
        assert isinstance(RemoveDuplicatesFilter().from_dict(d), RemoveDuplicatesFilter)


class TestRemoveExistingFilter(unittest.TestCase):
    def setUp(self):
        with open(f"{TEST_FILES_DIR}/TiO2_entries.json") as fp:
            entries = json.load(fp, cls=MontyDecoder)
        self._struct_list = [e.structure for e in entries]
        self._sm = StructureMatcher()
        self._exisiting_structures = self._struct_list[:-1]

    def test_filter(self):
        fil = RemoveExistingFilter(self._exisiting_structures)
        transmuter = StandardTransmuter.from_structures(self._struct_list)
        transmuter.apply_filter(fil)
        assert len(transmuter.transformed_structures) == 1
        assert self._sm.fit(
            self._struct_list[-1],
            transmuter.transformed_structures[-1].final_structure,
        )
