# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import pytest
import os

from collections import defaultdict
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.entries.computed_entries import ComputedEntry, \
    ComputedStructureEntry, EnergyAdjustment, ConstantEnergyAdjustment, \
    CompositionEnergyAdjustment, TempEnergyAdjustment, ManualEnergyAdjustment

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

filepath = os.path.join(test_dir, 'vasprun.xml')
vasprun = Vasprun(filepath)


def testEnergyAdjustment():
    ea = EnergyAdjustment(10)
    assert ea.name == "Manual adjustment"
    assert ea.cls == "None"
    assert ea.description == ""
    ead = ea.as_dict()
    ea2 = EnergyAdjustment.from_dict(ead)
    assert str(ead) == str(ea2.as_dict())


def testManualEnergyAdjustment():
    ea = ManualEnergyAdjustment(10)
    assert ea.name == "Manual energy adjustment"
    assert ea.value == 10
    assert ea.description == "Manual energy adjustment (10.000 eV)"


def testConstantEnergyAdjustment():
    ea = ConstantEnergyAdjustment(8)
    assert ea.name == "Constant energy adjustment"
    assert ea.value == 8
    assert ea.description == "Constant energy adjustment (8.000 eV)"


def testCompositionEnergyAdjustment():
    ea = CompositionEnergyAdjustment(2, 2, "H")
    assert ea.name == "H composition energy adjustment"
    assert ea.value == 4
    assert ea.description == "H Composition-based energy adjustment (2.000 eV/atom x 2 atoms)"


def testTempEnergyAdjustment():
    ea = TempEnergyAdjustment(-0.1, 298, 5, "entropy")
    assert ea.name == "entropy"
    assert ea.value == -0.1 * 298 * 5
    assert ea.n_atoms == 5
    assert ea.temp == 298
    assert ea.description == "Temperature-based entropy adjustment (-0.1000 eV/K/atom x 298 K x 5 atoms)"


class ComputedEntryTest(unittest.TestCase):

    def setUp(self):
        self.entry = ComputedEntry(vasprun.final_structure.composition,
                                   vasprun.final_energy,
                                   parameters=vasprun.incar)
        self.entry2 = ComputedEntry({"Fe": 2, "O": 3}, 2.3)
        self.entry3 = ComputedEntry("Fe2O3", 2.3)
        self.entry4 = ComputedEntry("Fe2O3", 2.3, entry_id=1)
        self.entry5 = ComputedEntry("Fe6O9", 6.9)
        ea = ConstantEnergyAdjustment(-5, name="Dummy adjustment")
        self.entry6 = ComputedEntry("Fe6O9", 6.9, correction=-10)
        self.entry7 = ComputedEntry("Fe6O9", 6.9, energy_adjustments=[ea])

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
        self.assertAlmostEqual(entry.energy_adjustments[0].value, 1/3)
        entry.normalize("atom")
        self.assertEqual(entry.composition.formula, "Fe0.4 O0.6")
        self.assertAlmostEqual(entry.uncorrected_energy, 6.9 / 15)
        self.assertAlmostEqual(entry.correction, 1 / 15)
        self.assertAlmostEqual(entry.energy * 15, 6.9 + 1)
        self.assertAlmostEqual(entry.energy_adjustments[0].value, 1/15)

    def test_normalize_energy_adjustments(self):
        ealist = [ManualEnergyAdjustment(5),
                  ConstantEnergyAdjustment(5),
                  CompositionEnergyAdjustment(1, 5, "Na"),
                  TempEnergyAdjustment(0.005, 100, 10)
                  ]
        entry = ComputedEntry("Na5Cl5", 6.9, energy_adjustments=ealist)
        assert entry.correction == 20
        entry.normalize()
        assert entry.correction == 4
        for ea in entry.energy_adjustments:
            assert ea.value == 1

    def test_to_from_dict(self):
        d = self.entry.as_dict()
        e = ComputedEntry.from_dict(d)
        self.assertAlmostEqual(e.energy, -269.38319884)

    def test_to_from_dict_with_adjustment(self):
        """
        Legacy case where adjustment was provided manually
        """
        d = self.entry6.as_dict()
        e = ComputedEntry.from_dict(d)
        self.assertAlmostEqual(e.uncorrected_energy, 6.9)
        self.assertEqual(e.energy_adjustments[0].value, self.entry6.energy_adjustments[0].value)

    def test_to_from_dict_with_adjustment_2(self):
        """
        Modern case where correction was provided manually
        """
        d = self.entry7.as_dict()
        e = ComputedEntry.from_dict(d)
        self.assertAlmostEqual(e.uncorrected_energy, 6.9)
        self.assertEqual(e.energy_adjustments[0].value, self.entry7.energy_adjustments[0].value)

    def test_to_from_dict_with_adjustment_3(self):
        """
        Legacy case where the entry was serialized before the energy_adjustment
        attribute was part of ComputedEntry
        """
        # same as entry6
        d = {'@module': 'pymatgen.entries.computed_entries',
             '@class': 'ComputedEntry',
             'energy': 6.9,
             'composition': defaultdict(float, {'Fe': 6.0, 'O': 9.0}),
             'parameters': {},
             'data': {},
             'entry_id': None,
             'correction': -10}
        e = ComputedEntry.from_dict(d)
        self.assertAlmostEqual(e.uncorrected_energy, 6.9)
        self.assertAlmostEqual(e.correction, -10)
        assert len(e.energy_adjustments) == 1

    def test_conflicting_correction_adjustment(self):
        """
        Should raise a ValueError if a user tries to manually set both the correction
        and energy_adjustment, even if the values match.
        """
        ea = ConstantEnergyAdjustment(-10, name="Dummy adjustment")
        with pytest.raises(ValueError, match="Argument conflict!"):
            ComputedEntry("Fe6O9", 6.9, correction=-10, energy_adjustments=[ea])

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
