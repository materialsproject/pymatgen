# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.entries.computed_entries import ComputedEntry, \
    ComputedStructureEntry

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

filepath = os.path.join(test_dir, 'vasprun.xml')
vasprun = Vasprun(filepath)


class ComputedEntryTest(unittest.TestCase):

    def setUp(self):
        self.entry = ComputedEntry(vasprun.final_structure.composition,
                                   vasprun.final_energy,
                                   parameters=vasprun.incar)
        self.entry2 = ComputedEntry({"Fe": 2, "O": 3}, 2.3)
        self.entry3 = ComputedEntry("Fe2O3", 2.3)
        self.entry4 = ComputedEntry("Fe2O3", 2.3, entry_id=1)
        self.entry5 = ComputedEntry("Fe6O9", 6.9)

    def test_energy(self):
        self.assertAlmostEqual(self.entry.energy, -269.38319884)
        self.entry.correction = 1.0
        self.assertAlmostEqual(self.entry.energy, -268.38319884)
        self.assertAlmostEqual(self.entry3.energy_per_atom, 2.3 / 5)

    def test_composition(self):
        self.assertEqual(self.entry.composition.reduced_formula, "LiFe4(PO4)4")
        self.assertEqual(self.entry2.composition.reduced_formula, "Fe2O3")
        self.assertEqual(self.entry5.composition.reduced_formula, "Fe2O3")
        self.assertEqual(self.entry5.composition.get_reduced_formula_and_factor()[1], 3)

    def test_normalize(self):
        entry = ComputedEntry("Fe6O9", 6.9, correction=1)
        entry.normalize()
        self.assertEqual(entry.composition.formula, "Fe2 O3")
        self.assertAlmostEqual(entry.uncorrected_energy, 6.9 / 3)
        self.assertAlmostEqual(entry.correction, 1 / 3)
        self.assertAlmostEqual(entry.energy * 3, 6.9 + 1)
        entry.normalize("atom")
        self.assertEqual(entry.composition.formula, "Fe0.4 O0.6")
        self.assertAlmostEqual(entry.uncorrected_energy, 6.9 / 15)
        self.assertAlmostEqual(entry.correction, 1 / 15)
        self.assertAlmostEqual(entry.energy * 15, 6.9 + 1)

    def test_to_from_dict(self):
        d = self.entry.as_dict()
        e = ComputedEntry.from_dict(d)
        self.assertAlmostEqual(e.energy, -269.38319884)

    def test_entry_id(self):
        self.assertEqual(self.entry4.entry_id, 1)
        self.assertEqual(self.entry2.entry_id, None)

    def test_str(self):
        self.assertIsNotNone(str(self.entry))

    def test_sulfide_energy(self):
        self.entry = ComputedEntry("BaS", -10.21249155)
        self.assertAlmostEqual(self.entry.energy, -10.21249155)
        self.assertAlmostEqual(self.entry.energy_per_atom, -10.21249155 / 2)
        self.entry.correction = 1.0
        self.assertAlmostEqual(self.entry.energy, -9.21249155)

    def test_is_element(self):
        entry = ComputedEntry("Fe3", 2.3)
        self.assertTrue(entry.is_element)


class ComputedStructureEntryTest(unittest.TestCase):

    def setUp(self):
        self.entry = ComputedStructureEntry(vasprun.final_structure,
                                            vasprun.final_energy,
                                            parameters=vasprun.incar)

    def test_energy(self):
        self.assertAlmostEqual(self.entry.energy, -269.38319884)
        self.entry.correction = 1.0
        self.assertAlmostEqual(self.entry.energy, -268.38319884)

    def test_composition(self):
        self.assertEqual(self.entry.composition.reduced_formula, "LiFe4(PO4)4")

    def test_to_from_dict(self):
        d = self.entry.as_dict()
        e = ComputedStructureEntry.from_dict(d)
        self.assertAlmostEqual(e.energy, -269.38319884)

    def test_str(self):
        self.assertIsNotNone(str(self.entry))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
