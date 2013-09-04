#!/usr/bin/env python

"""
Created on Nov 10, 2012

@author: shyue
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Nov 10, 2012"

import unittest


from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition, CompositionError
import random


class CompositionTest(unittest.TestCase):

    def setUp(self):
        self.comp = list()
        self.comp.append(Composition.from_formula("Li3Fe2(PO4)3"))
        self.comp.append(Composition.from_formula("Li3Fe(PO4)O"))
        self.comp.append(Composition.from_formula("LiMn2O4"))
        self.comp.append(Composition.from_formula("Li4O4"))
        self.comp.append(Composition.from_formula("Li3Fe2Mo3O12"))
        self.comp.append(Composition.from_formula("Li3Fe2((PO4)3(CO3)5)2"))
        self.comp.append(Composition.from_formula("Li1.5Si0.5"))
        self.comp.append(Composition.from_formula("ZnOH"))

        self.indeterminate_comp = []
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("Co1",
                                                                       True)
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("Co1",
                                                                       False)
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("co2o3")
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("ncalu")
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("calun")
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula(
                "liCoo2n (pO4)2")
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula(
                "(co)2 (PO)4")
        )
        self.indeterminate_comp.append(
            Composition.ranked_compositions_from_indeterminate_formula("Fee3"))

    def test_init_(self):
        self.assertRaises(CompositionError, Composition, {"H": -0.1})
        f = {'Fe': 4, 'Li': 4, 'O': 16, 'P': 4}
        self.assertEqual("Li4 Fe4 P4 O16", Composition(f).formula)
        f = {None: 4, 'Li': 4, 'O': 16, 'P': 4}
        self.assertRaises(ValueError, Composition, f)
        f = {1: 2, 8: 1}
        self.assertEqual("H2 O1", Composition(f).formula)
        self.assertEqual("Na2 O1", Composition(Na=2, O=1).formula)

    def test_average_electroneg(self):
        val = [2.7224999999999997, 2.4160000000000004, 2.5485714285714285,
               2.21, 2.718, 3.08, 1.21, 2.43]
        for i, c in enumerate(self.comp):
            self.assertAlmostEqual(c.average_electroneg,
                                   val[i])

    def test_formula(self):
        correct_formulas = ['Li3 Fe2 P3 O12', 'Li3 Fe1 P1 O5', 'Li1 Mn2 O4',
                            'Li4 O4', 'Li3 Fe2 Mo3 O12', 'Li3 Fe2 P6 C10 O54',
                            'Li1.5 Si0.5', 'Zn1 H1 O1']
        all_formulas = [c.formula for c in self.comp]
        self.assertEqual(all_formulas, correct_formulas)
        self.assertRaises(CompositionError, Composition.from_formula,
                          "(co2)(po4)2")

    def test_mixed_valence(self):
        comp = Composition({"Fe2+": 2, "Fe3+": 4, "Li+": 8})
        self.assertEqual(comp.reduced_formula, "Li4Fe3")
        self.assertEqual(comp.alphabetical_formula, "Fe6 Li8")
        self.assertEqual(comp.formula, "Li8 Fe6")

    def test_indeterminate_formula(self):
        correct_formulas = [["Co1"], ["Co1", "C1 O1"], ["Co2 O3", "C1 O5"],
                            ["N1 Ca1 Lu1", "U1 Al1 C1 N1"],
                            ["N1 Ca1 Lu1", "U1 Al1 C1 N1"],
                            ["Li1 Co1 P2 N1 O10", "Li1 P2 C1 N1 O11",
                             "Li1 Co1 Po8 N1 O2", "Li1 Po8 C1 N1 O3"],
                            ["Co2 P4 O4", "Co2 Po4", "P4 C2 O6",
                             "Po4 C2 O2"], []]
        for i, c in enumerate(correct_formulas):
            self.assertEqual([Composition.from_formula(comp) for comp in c],
                             self.indeterminate_comp[i])

    def test_alphabetical_formula(self):
        correct_formulas = ['Fe2 Li3 O12 P3', 'Fe1 Li3 O5 P1', 'Li1 Mn2 O4',
                            'Li4 O4', 'Fe2 Li3 Mo3 O12', 'C10 Fe2 Li3 O54 P6',
                            'Li1.5 Si0.5', 'H1 O1 Zn1']
        all_formulas = [c.alphabetical_formula for c in self.comp]
        self.assertEqual(all_formulas, correct_formulas)

    def test_reduced_composition(self):
        correct_reduced_formulas = ['Li3Fe2(PO4)3', 'Li3FePO5', 'LiMn2O4',
                                    'Li2O2', 'Li3Fe2(MoO4)3',
                                    'Li3Fe2P6(C5O27)2', 'Li1.5Si0.5', 'ZnHO']
        for i in xrange(len(self.comp)):
            self.assertEqual(self.comp[i]
                             .get_reduced_composition_and_factor()[0],
                             Composition
                             .from_formula(correct_reduced_formulas[i]))

    def test_reduced_formula(self):
        correct_reduced_formulas = ['Li3Fe2(PO4)3', 'Li3FePO5', 'LiMn2O4',
                                    'Li2O2', 'Li3Fe2(MoO4)3',
                                    'Li3Fe2P6(C5O27)2', 'Li1.5Si0.5', 'ZnHO']
        all_formulas = [c.reduced_formula for c in self.comp]
        self.assertEqual(all_formulas, correct_reduced_formulas)

    def test_num_atoms(self):
        correct_num_atoms = [20, 10, 7, 8, 20, 75, 2, 3]
        all_natoms = [c.num_atoms for c in self.comp]
        self.assertEqual(all_natoms, correct_num_atoms)

    def test_weight(self):
        correct_weights = [417.427086, 187.63876199999999, 180.81469, 91.7616,
                           612.3258, 1302.430172, 24.454250000000002, 82.41634]
        all_weights = [c.weight for c in self.comp]
        self.assertAlmostEqual(all_weights, correct_weights, 5)

    def test_get_atomic_fraction(self):
        correct_at_frac = {"Li": 0.15, "Fe": 0.1, "P": 0.15, "O": 0.6}
        for el in ["Li", "Fe", "P", "O"]:
            self.assertEqual(self.comp[0].get_atomic_fraction(el),
                             correct_at_frac[el],
                             "Wrong computed atomic fractions")
        self.assertEqual(self.comp[0].get_atomic_fraction("S"), 0,
                         "Wrong computed atomic fractions")

    def test_anonymized_formula(self):
        expected_formulas = ['A2B3C3D12', 'ABC3D5', 'AB2C4', 'A2B2',
                             'A2B3C3D12', 'A2B3C6D10E54', 'A0.5B1.5', 'ABC']
        for i in xrange(len(self.comp)):
            self.assertEqual(self.comp[i].anonymized_formula,
                             expected_formulas[i])

    def test_get_wt_fraction(self):
        correct_wt_frac = {"Li": 0.0498841610868, "Fe": 0.267567687258,
                           "P": 0.222604831158, "O": 0.459943320496}
        for el in ["Li", "Fe", "P", "O"]:
            self.assertAlmostEqual(correct_wt_frac[el],
                                   self.comp[0].get_wt_fraction(el),
                                   5, "Wrong computed weight fraction")
        self.assertEqual(self.comp[0].get_wt_fraction(Element("S")), 0,
                         "Wrong computed weight fractions")

    def test_from_dict(self):
        sym_dict = {"Fe": 6, "O": 8}
        self.assertEqual(Composition.from_dict(sym_dict).reduced_formula,
                         "Fe3O4",
                         "Creation form sym_amount dictionary failed!")
        comp = Composition({"Fe2+": 2, "Fe3+": 4, "O2-": 8})
        comp2 = Composition.from_dict(comp.to_dict)
        self.assertEqual(comp, comp2)

    def test_to_dict(self):
        c = Composition.from_dict({'Fe': 4, 'O': 6})
        d = c.to_dict
        correct_dict = {'Fe': 4.0, 'O': 6.0}
        self.assertEqual(d['Fe'], correct_dict['Fe'])
        self.assertEqual(d['O'], correct_dict['O'])
        correct_dict = {'Fe': 2.0, 'O': 3.0}
        d = c.to_reduced_dict
        self.assertEqual(d['Fe'], correct_dict['Fe'])
        self.assertEqual(d['O'], correct_dict['O'])

    def test_add(self):
        self.assertEqual((self.comp[0] + self.comp[2]).formula,
                         "Li4 Mn2 Fe2 P3 O16",
                         "Incorrect composition after addition!")
        self.assertEqual((self.comp[3] + {"Fe": 4, "O": 4}).formula,
                         "Li4 Fe4 O8", "Incorrect composition after addition!")

    def test_sub(self):
        self.assertEqual((self.comp[0]
                          - Composition.from_formula("Li2O")).formula,
                         "Li1 Fe2 P3 O11",
                         "Incorrect composition after addition!")
        self.assertEqual((self.comp[0] - {"Fe": 2, "O": 3}).formula,
                         "Li3 P3 O9")

    def test_mul(self):
        self.assertEqual((self.comp[0] * 4).formula, "Li12 Fe8 P12 O48")
        self.assertEqual((3 * self.comp[1]).formula, "Li9 Fe3 P3 O15")

    def test_equals(self):
        random_z = random.randint(1, 92)
        fixed_el = Element.from_Z(random_z)
        other_z = random.randint(1, 92)
        while other_z == random_z:
            other_z = random.randint(1, 92)
        comp1 = Composition({fixed_el: 1, Element.from_Z(other_z): 0})
        other_z = random.randint(1, 92)
        while other_z == random_z:
            other_z = random.randint(1, 92)
        comp2 = Composition({fixed_el: 1, Element.from_Z(other_z): 0})
        self.assertEqual(comp1, comp2,
                         "Composition equality test failed. " +
                         "%s should be equal to %s" % (comp1.formula,
                                                       comp2.formula))
        self.assertEqual(comp1.__hash__(), comp2.__hash__(),
                         "Hashcode equality test failed!")

    def test_almost_equals(self):
        c1 = Composition({'Fe': 2.0, 'O': 3.0, 'Mn': 0})
        c2 = Composition({'O': 3.2, 'Fe': 1.9, 'Zn': 0})
        c3 = Composition({'Ag': 2.0, 'O': 3.0})
        c4 = Composition({'Fe': 2.0, 'O': 3.0, 'Ag': 2.0})
        self.assertTrue(c1.almost_equals(c2, rtol=0.1))
        self.assertFalse(c1.almost_equals(c2, rtol=0.01))
        self.assertFalse(c1.almost_equals(c3, rtol=0.1))
        self.assertFalse(c1.almost_equals(c4, rtol=0.1))

    def test_equality(self):
        self.assertTrue(self.comp[0].__eq__(self.comp[0]))
        self.assertFalse(self.comp[0].__eq__(self.comp[1]))
        self.assertFalse(self.comp[0].__ne__(self.comp[0]))
        self.assertTrue(self.comp[0].__ne__(self.comp[1]))

    def test_get_fractional_composition(self):
        for c in self.comp:
            self.assertAlmostEqual(c.get_fractional_composition().num_atoms, 1)

    def test_init_numerical_tolerance(self):
        self.assertEqual(Composition({'B':1, 'C':-1e-12}), Composition('B'))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
