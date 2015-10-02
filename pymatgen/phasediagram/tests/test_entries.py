# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os

from pymatgen.phasediagram.entries import PDEntryIO, PDEntry, \
    GrandPotPDEntry, TransformedPDEntry
from pymatgen.core.periodic_table import Element, DummySpecie
from pymatgen.core.composition import Composition


class PDEntryTest(unittest.TestCase):
    '''
    Test all functions using a ficitious entry
    '''
    def setUp(self):
        comp = Composition("LiFeO2")
        self.entry = PDEntry(comp, 53)
        self.gpentry = GrandPotPDEntry(self.entry, {Element('O'): 1.5})

    def test_get_energy(self):
        self.assertEqual(self.entry.energy, 53, "Wrong energy!")
        self.assertEqual(self.gpentry.energy, 50, "Wrong energy!")

    def test_get_energy_per_atom(self):
        self.assertEqual(self.entry.energy_per_atom, 53.0 / 4,
                          "Wrong energy per atom!")
        self.assertEqual(self.gpentry.energy_per_atom, 50.0 / 2,
                          "Wrong energy per atom!")

    def test_get_name(self):
        self.assertEqual(self.entry.name, 'LiFeO2', "Wrong name!")
        self.assertEqual(self.gpentry.name, 'LiFeO2', "Wrong name!")

    def test_get_composition(self):
        comp = self.entry.composition
        expected_comp = Composition('LiFeO2')
        self.assertEqual(comp, expected_comp, "Wrong composition!")
        comp = self.gpentry.composition
        expected_comp = Composition("LiFe")
        self.assertEqual(comp, expected_comp, "Wrong composition!")

    def test_is_element(self):
        self.assertFalse(self.entry.is_element)
        self.assertFalse(self.gpentry.is_element)

    def test_to_from_dict(self):
        d = self.entry.as_dict()
        gpd = self.gpentry.as_dict()
        entry = PDEntry.from_dict(d)

        self.assertEqual(entry.name, 'LiFeO2', "Wrong name!")
        self.assertEqual(entry.energy_per_atom, 53.0 / 4)
        gpentry = GrandPotPDEntry.from_dict(gpd)
        self.assertEqual(gpentry.name, 'LiFeO2', "Wrong name!")
        self.assertEqual(gpentry.energy_per_atom, 50.0 / 2)

    def test_str(self):
        self.assertIsNotNone(str(self.entry))


class TransformedPDEntryTest(unittest.TestCase):
    '''
    Test all functions using a ficitious entry
    '''
    def setUp(self):
        comp = Composition("LiFeO2")
        entry = PDEntry(comp, 53)
        self.transformed_entry = TransformedPDEntry({DummySpecie('Xa'): 1,
                                                     DummySpecie("Xb"): 1},
                                                    entry)

    def test_get_energy(self):
        self.assertEqual(self.transformed_entry.energy, 53, "Wrong energy!")
        self.assertEqual(self.transformed_entry.original_entry.energy, 53.0)

    def test_get_energy_per_atom(self):
        self.assertEqual(self.transformed_entry.energy_per_atom, 53.0 / 2)

    def test_get_name(self):
        self.assertEqual(self.transformed_entry.name, 'LiFeO2', "Wrong name!")

    def test_get_composition(self):
        comp = self.transformed_entry.composition
        expected_comp = Composition({DummySpecie('Xa'): 1,
                                     DummySpecie('Xb'): 1})
        self.assertEqual(comp, expected_comp, "Wrong composition!")

    def test_is_element(self):
        self.assertFalse(self.transformed_entry.is_element)

    def test_to_from_dict(self):
        d = self.transformed_entry.as_dict()
        entry = TransformedPDEntry.from_dict(d)
        self.assertEqual(entry.name, 'LiFeO2', "Wrong name!")
        self.assertEqual(entry.energy_per_atom, 53.0 / 2)

    def test_str(self):
        self.assertIsNotNone(str(self.transformed_entry))


class PDEntryIOTestCase(unittest.TestCase):

    def test_read_csv(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements,
         entries) = PDEntryIO.from_csv(os.path.join(module_dir,
                                                    "pdentries_test.csv"))
        self.assertEqual(elements,
                         [Element('Li'), Element('Fe'), Element('O')],
                         "Wrong elements!")
        self.assertEqual(len(entries), 492, "Wrong number of entries!")


if __name__ == '__main__':
    unittest.main()
