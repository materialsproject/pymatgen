# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import json
import os
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
from pymatgen.util.testing import PymatgenTest


class ContainsSpecieFilterTest(PymatgenTest):
    def test_filtering(self):
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = Lattice([[3.0, 0.0, 0.0], [1.0, 3.0, 0.00], [0.00, -2.0, 3.0]])
        s = Structure(
            lattice,
            [
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
            ],
            coords,
        )

        species1 = [Species("Si", 5), Species("Mg", 2)]
        f1 = ContainsSpecieFilter(species1, strict_compare=True, AND=False)
        self.assertFalse(f1.test(s), "Incorrect filter")
        f2 = ContainsSpecieFilter(species1, strict_compare=False, AND=False)
        self.assertTrue(f2.test(s), "Incorrect filter")
        species2 = [Species("Si", 4), Species("Mg", 2)]
        f3 = ContainsSpecieFilter(species2, strict_compare=True, AND=False)
        self.assertTrue(f3.test(s), "Incorrect filter")
        f4 = ContainsSpecieFilter(species2, strict_compare=False, AND=False)
        self.assertTrue(f4.test(s), "Incorrect filter")

        species3 = [Species("Si", 5), Species("O", -2)]
        f5 = ContainsSpecieFilter(species3, strict_compare=True, AND=True)
        self.assertFalse(f5.test(s), "Incorrect filter")
        f6 = ContainsSpecieFilter(species3, strict_compare=False, AND=True)
        self.assertTrue(f6.test(s), "Incorrect filter")
        species4 = [Species("Si", 4), Species("Mg", 2)]
        f7 = ContainsSpecieFilter(species4, strict_compare=True, AND=True)
        self.assertFalse(f7.test(s), "Incorrect filter")
        f8 = ContainsSpecieFilter(species4, strict_compare=False, AND=True)
        self.assertFalse(f8.test(s), "Incorrect filter")

    def test_to_from_dict(self):
        species1 = ["Si5+", "Mg2+"]
        f1 = ContainsSpecieFilter(species1, strict_compare=True, AND=False)
        d = f1.as_dict()
        self.assertIsInstance(ContainsSpecieFilter.from_dict(d), ContainsSpecieFilter)


class SpecieProximityFilterTest(PymatgenTest):
    def test_filter(self):
        s = self.get_structure("Li10GeP2S12")
        sf = SpecieProximityFilter({"Li": 1})
        self.assertTrue(sf.test(s))
        sf = SpecieProximityFilter({"Li": 2})
        self.assertFalse(sf.test(s))
        sf = SpecieProximityFilter({"P": 1})
        self.assertTrue(sf.test(s))
        sf = SpecieProximityFilter({"P": 5})
        self.assertFalse(sf.test(s))

    def test_to_from_dict(self):
        sf = SpecieProximityFilter({"Li": 1})
        d = sf.as_dict()
        self.assertIsInstance(SpecieProximityFilter.from_dict(d), SpecieProximityFilter)


class RemoveDuplicatesFilterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "TiO2_entries.json"), "r") as fp:
            entries = json.load(fp, cls=MontyDecoder)
        self._struct_list = [e.structure for e in entries]
        self._sm = StructureMatcher()

    def test_filter(self):
        transmuter = StandardTransmuter.from_structures(self._struct_list)
        fil = RemoveDuplicatesFilter()
        transmuter.apply_filter(fil)
        self.assertEqual(len(transmuter.transformed_structures), 11)

    def test_to_from_dict(self):
        fil = RemoveDuplicatesFilter()
        d = fil.as_dict()
        self.assertIsInstance(RemoveDuplicatesFilter().from_dict(d), RemoveDuplicatesFilter)


class RemoveExistingFilterTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "TiO2_entries.json"), "r") as fp:
            entries = json.load(fp, cls=MontyDecoder)
        self._struct_list = [e.structure for e in entries]
        self._sm = StructureMatcher()
        self._exisiting_structures = self._struct_list[:-1]

    def test_filter(self):
        fil = RemoveExistingFilter(self._exisiting_structures)
        transmuter = StandardTransmuter.from_structures(self._struct_list)
        transmuter.apply_filter(fil)
        self.assertEqual(len(transmuter.transformed_structures), 1)
        self.assertTrue(
            self._sm.fit(
                self._struct_list[-1],
                transmuter.transformed_structures[-1].final_structure,
            )
        )


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
