# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest

from pymatgen.core.ion import Ion
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition, CompositionError
import random


class IonTest(unittest.TestCase):

    def setUp(self):
        self.comp = list()
        self.comp.append(Ion.from_formula("Li+"))
        self.comp.append(Ion.from_formula("MnO4-"))
        self.comp.append(Ion.from_formula("Mn++"))
        self.comp.append(Ion.from_formula("PO3-2"))
        self.comp.append(Ion.from_formula("Fe(CN)6-3"))
        self.comp.append(Ion.from_formula("Fe(CN)6----"))
        self.comp.append(Ion.from_formula("Fe2((PO4)3(CO3)5)2-3"))
        self.comp.append(Ion.from_formula("Ca[2+]"))
        self.comp.append(Ion.from_formula("NaOH(aq)"))

    def test_init_(self):
        c = Composition({'Fe': 4, 'O': 16, 'P': 4})
        charge = 4
        self.assertEqual("Fe4 P4 O16 +4", Ion(c, charge).formula)
        f = {1: 1, 8: 1}
        charge = -1
        self.assertEqual("H1 O1 -1", Ion(Composition(f), charge).formula)
        self.assertEqual("S2 O3 -2", Ion(Composition(S=2, O=3), -2).formula)

    def test_formula(self):
        correct_formulas = ['Li1 +1', 'Mn1 O4 -1', 'Mn1 +2', 'P1 O3 -2',\
                            'Fe1 C6 N6 -3', 'Fe1 C6 N6 -4', \
                            'Fe2 P6 C10 O54 -3', 'Ca1 +2', 'Na1 H1 O1']
        all_formulas = [c.formula for c in self.comp]
        self.assertEqual(all_formulas, correct_formulas)
        self.assertRaises(CompositionError, Ion.from_formula, "(co2)(po4)2")

    def test_mixed_valence(self):
        comp = Ion(Composition({"Fe2+": 2, "Fe3+": 4, "Li+": 8}))
        self.assertEqual(comp.reduced_formula, "Li4Fe3(aq)")
        self.assertEqual(comp.alphabetical_formula, "Fe6 Li8")
        self.assertEqual(comp.formula, "Li8 Fe6")

    def test_alphabetical_formula(self):
        correct_formulas = ['Li1 +1', 'Mn1 O4 -1', 'Mn1 +2', 'O3 P1 -2',\
                            'C6 Fe1 N6 -3', 'C6 Fe1 N6 -4', \
                            'C10 Fe2 O54 P6 -3', 'Ca1 +2', "H1 Na1 O1"]
        all_formulas = [c.alphabetical_formula for c in self.comp]
        self.assertEqual(all_formulas, correct_formulas)

    def test_num_atoms(self):
        correct_num_atoms = [1, 5, 1, 4, 13, 13, 72, 1, 3]
        all_natoms = [c.num_atoms for c in self.comp]
        self.assertEqual(all_natoms, correct_num_atoms)

    def test_anonymized_formula(self):
        expected_formulas = ['A+1', 'AB4-1', 'A+2', 'AB3-2',
                             'AB6C6-3', 'AB6C6-4', 'AB3C5D27-3',
                             'A+2', 'ABC']
        for i in range(len(self.comp)):
            self.assertEqual(self.comp[i].anonymized_formula,
                             expected_formulas[i])

    def test_from_dict(self):
        sym_dict = {"P": 1, "O": 4, 'charge': -2}
        self.assertEqual(Ion.from_dict(sym_dict).reduced_formula,
                         "PO4[2-]",
                         "Creation form sym_amount dictionary failed!")

    def test_as_dict(self):
        c = Ion.from_dict({'Mn': 1, 'O': 4, 'charge': -1})
        d = c.as_dict()
        correct_dict = {'Mn': 1.0, 'O': 4.0, 'charge': -1.0}
        self.assertEqual(d, correct_dict)
        self.assertEqual(d['charge'], correct_dict['charge'])
        correct_dict = {'Mn': 1.0, 'O': 4.0, 'charge': -1}
        d = c.to_reduced_dict
        self.assertEqual(d, correct_dict)
        self.assertEqual(d['charge'], correct_dict['charge'])

    def test_equals(self):
        random_z = random.randint(1, 92)
        fixed_el = Element.from_Z(random_z)
        other_z = random.randint(1, 92)
        while other_z == random_z:
            other_z = random.randint(1, 92)
        comp1 = Ion(Composition({fixed_el: 1, Element.from_Z(other_z): 0}), 1)
        other_z = random.randint(1, 92)
        while other_z == random_z:
            other_z = random.randint(1, 92)
        comp2 = Ion(Composition({fixed_el: 1, Element.from_Z(other_z): 0}), 1)
        self.assertEqual(comp1, comp2,
                         "Composition equality test failed. " +
                         "%s should be equal to %s" % (comp1.formula,
                                                       comp2.formula))
        self.assertEqual(comp1.__hash__(), comp2.__hash__(),
                         "Hashcode equality test failed!")

    def test_equality(self):
        self.assertTrue(self.comp[0] == (self.comp[0]))
        self.assertFalse(self.comp[0] == (self.comp[1]))
        self.assertFalse(self.comp[0] != (self.comp[0]))
        self.assertTrue(self.comp[0] != (self.comp[1]))

    def test_mul(self):
        self.assertEqual((self.comp[1] * 4).formula, "Mn4 O16 -4",
                         "Incorrect composition after addition!")

    def test_len(self):
        self.assertEqual(len(self.comp[1]), 2, "Lengths are not equal!")

if __name__ == '__main__':
    unittest.main()
