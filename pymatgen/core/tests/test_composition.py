# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Created on Nov 10, 2012

@author: shyue
"""

from pymatgen.util.testing import PymatgenTest


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Nov 10, 2012"

import unittest

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.composition import Composition, CompositionError, \
    ChemicalPotential

import random


class CompositionTest(PymatgenTest):

    def setUp(self):
        self.comp = list()
        self.comp.append(Composition("Li3Fe2(PO4)3"))
        self.comp.append(Composition("Li3Fe(PO4)O"))
        self.comp.append(Composition("LiMn2O4"))
        self.comp.append(Composition("Li4O4"))
        self.comp.append(Composition("Li3Fe2Mo3O12"))
        self.comp.append(Composition("Li3Fe2((PO4)3(CO3)5)2"))
        self.comp.append(Composition("Li1.5Si0.5"))
        self.comp.append(Composition("ZnOH"))

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

    def test_immutable(self):
        try:
            self.comp[0]["Fe"] = 1
        except Exception as ex:
            self.assertIsInstance(ex, TypeError)

        try:
            del self.comp[0]["Fe"]
        except Exception as ex:
            self.assertIsInstance(ex, TypeError)

    def test_in(self):
        self.assertIn("Fe", self.comp[0])
        self.assertNotIn("Fe", self.comp[2])
        self.assertIn(Element("Fe"), self.comp[0])
        self.assertEqual(self.comp[0]["Fe"], 2)
        self.assertEqual(self.comp[0]["Mn"], 0)
        self.assertRaises(TypeError, self.comp[0].__getitem__, "Hello")
        self.assertRaises(TypeError, self.comp[0].__getitem__, "Vac")

    def test_hill_formula(self):
        c = Composition("CaCO3")
        self.assertEqual(c.hill_formula, "C Ca O3")
        c = Composition("C2H5OH")
        self.assertEqual(c.hill_formula, "C2 H6 O")

    def test_init_(self):
        self.assertRaises(CompositionError, Composition, {"H": -0.1})
        f = {'Fe': 4, 'Li': 4, 'O': 16, 'P': 4}
        self.assertEqual("Li4 Fe4 P4 O16", Composition(f).formula)
        f = {None: 4, 'Li': 4, 'O': 16, 'P': 4}
        self.assertRaises(TypeError, Composition, f)
        f = {1: 2, 8: 1}
        self.assertEqual("H2 O1", Composition(f).formula)
        self.assertEqual("Na2 O1", Composition(Na=2, O=1).formula)

        c = Composition({'S': Composition.amount_tolerance / 2})
        self.assertEqual(len(c.elements), 0)

    def test_average_electroneg(self):
        val = [2.7224999999999997, 2.4160000000000004, 2.5485714285714285,
               2.21, 2.718, 3.08, 1.21, 2.43]
        for i, c in enumerate(self.comp):
            self.assertAlmostEqual(c.average_electroneg,
                                   val[i])

    def test_total_electrons(self):
        test_cases = {'C': 6, 'SrTiO3': 84}
        for item in test_cases.keys():
            c = Composition(item)
            self.assertAlmostEqual(c.total_electrons, test_cases[item])

    def test_formula(self):
        correct_formulas = ['Li3 Fe2 P3 O12', 'Li3 Fe1 P1 O5', 'Li1 Mn2 O4',
                            'Li4 O4', 'Li3 Fe2 Mo3 O12', 'Li3 Fe2 P6 C10 O54',
                            'Li1.5 Si0.5', 'Zn1 H1 O1']
        all_formulas = [c.formula for c in self.comp]
        self.assertEqual(all_formulas, correct_formulas)
        self.assertRaises(CompositionError, Composition,
                          "(co2)(po4)2")

        self.assertEqual(Composition("K Na 2").reduced_formula, "KNa2")

        self.assertEqual(Composition("K3 Na 2").reduced_formula, "K3Na2")

        self.assertEqual(Composition("Na 3 Zr (PO 4) 3").reduced_formula,
                         "Na3Zr(PO4)3")

    def test_iupac_formula(self):
        correct_formulas = ['Li3 Fe2 P3 O12', 'Li3 Fe1 P1 O5', 'Li1 Mn2 O4',
                            'Li4 O4', 'Li3 Mo3 Fe2 O12', 'Li3 Fe2 C10 P6 O54',
                            'Li1.5 Si0.5', 'Zn1 H1 O1']
        all_formulas = [c.iupac_formula for c in self.comp]
        self.assertEqual(all_formulas, correct_formulas)

    def test_mixed_valence(self):
        comp = Composition({"Fe2+": 2, "Fe3+": 4, "Li+": 8})
        self.assertEqual(comp.reduced_formula, "Li4Fe3")
        self.assertEqual(comp.alphabetical_formula, "Fe6 Li8")
        self.assertEqual(comp.formula, "Li8 Fe6")

    def test_indeterminate_formula(self):
        correct_formulas = [["Co1"], ["Co1", "C1 O1"], ["Co2 O3", "C1 O5"],
                            ["N1 Ca1 Lu1", "U1 Al1 C1 N1"],
                            ["N1 Ca1 Lu1", "U1 Al1 C1 N1"],
                            ["Li1 Co1 P2 N1 O10", "Li1 Co1 Po8 N1 O2",
                             "Li1 P2 C1 N1 O11", "Li1 Po8 C1 N1 O3"],
                            ["Co2 P4 O4", "Co2 Po4", "P4 C2 O6",
                             "Po4 C2 O2"], []]
        for i, c in enumerate(correct_formulas):
            self.assertEqual([Composition(comp) for comp in c],
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
        for i in range(len(self.comp)):
            self.assertEqual(self.comp[i]
                             .get_reduced_composition_and_factor()[0],
                             Composition(correct_reduced_formulas[i]))

    def test_reduced_formula(self):
        correct_reduced_formulas = ['Li3Fe2(PO4)3', 'Li3FePO5', 'LiMn2O4',
                                    'Li2O2', 'Li3Fe2(MoO4)3',
                                    'Li3Fe2P6(C5O27)2', 'Li1.5Si0.5', 'ZnHO']
        all_formulas = [c.reduced_formula for c in self.comp]
        self.assertEqual(all_formulas, correct_reduced_formulas)

        # test iupac reduced formula (polyanions should still appear at the end)
        all_formulas = [c.get_reduced_formula_and_factor(iupac_ordering=True)[0]
                        for c in self.comp]
        self.assertEqual(all_formulas, correct_reduced_formulas)
        self.assertEqual(
            Composition('H6CN').get_integer_formula_and_factor(
                iupac_ordering=True)[0],
            'CNH6')

        # test rounding
        c = Composition({'Na': 2 - Composition.amount_tolerance / 2, 'Cl': 2})
        self.assertEqual('NaCl', c.reduced_formula)

    def test_integer_formula(self):
        correct_reduced_formulas = ['Li3Fe2(PO4)3', 'Li3FePO5', 'LiMn2O4',
                                    'Li2O2', 'Li3Fe2(MoO4)3',
                                    'Li3Fe2P6(C5O27)2', 'Li3Si', 'ZnHO']
        all_formulas = [c.get_integer_formula_and_factor()[0]
                        for c in self.comp]
        self.assertEqual(all_formulas, correct_reduced_formulas)
        self.assertEqual(Composition('Li0.5O0.25').get_integer_formula_and_factor(),
                         ('Li2O', 0.25))
        self.assertEqual(Composition('O0.25').get_integer_formula_and_factor(),
                         ('O2', 0.125))
        formula, factor = Composition(
            "Li0.16666667B1.0H1.0").get_integer_formula_and_factor()
        self.assertEqual(formula, 'Li(BH)6')
        self.assertAlmostEqual(factor, 1 / 6)

        # test iupac reduced formula (polyanions should still appear at the end)
        all_formulas = [c.get_integer_formula_and_factor(iupac_ordering=True)[0]
                        for c in self.comp]
        self.assertEqual(all_formulas, correct_reduced_formulas)
        self.assertEqual(
            Composition('H6CN0.5').get_integer_formula_and_factor(
                iupac_ordering=True),
            ('C2NH12', 0.5))

    def test_num_atoms(self):
        correct_num_atoms = [20, 10, 7, 8, 20, 75, 2, 3]

        all_natoms = [c.num_atoms for c in self.comp]
        self.assertEqual(all_natoms, correct_num_atoms)

    def test_weight(self):
        correct_weights = [417.427086, 187.63876199999999, 180.81469, 91.7616,
                           612.3258, 1302.430172, 24.454250000000002, 82.41634]
        all_weights = [c.weight for c in self.comp]
        self.assertArrayAlmostEqual(all_weights, correct_weights, 5)

    def test_get_atomic_fraction(self):
        correct_at_frac = {"Li": 0.15, "Fe": 0.1, "P": 0.15, "O": 0.6}
        for el in ["Li", "Fe", "P", "O"]:
            self.assertEqual(self.comp[0].get_atomic_fraction(el),
                             correct_at_frac[el],
                             "Wrong computed atomic fractions")
        self.assertEqual(self.comp[0].get_atomic_fraction("S"), 0,
                         "Wrong computed atomic fractions")

    def test_anonymized_formula(self):
        expected_formulas = ['A2B3C3D12', 'ABC3D5', 'AB2C4', 'AB',
                             'A2B3C3D12', 'A2B3C6D10E54', 'A0.5B1.5', 'ABC']
        for i in range(len(self.comp)):
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
        comp2 = Composition.from_dict(comp.as_dict())
        self.assertEqual(comp, comp2)

    def test_as_dict(self):
        c = Composition.from_dict({'Fe': 4, 'O': 6})
        d = c.as_dict()
        correct_dict = {'Fe': 4.0, 'O': 6.0}
        self.assertEqual(d['Fe'], correct_dict['Fe'])
        self.assertEqual(d['O'], correct_dict['O'])
        correct_dict = {'Fe': 2.0, 'O': 3.0}
        d = c.to_reduced_dict
        self.assertEqual(d['Fe'], correct_dict['Fe'])
        self.assertEqual(d['O'], correct_dict['O'])

    def test_pickle(self):
        for c in self.comp:
            self.serialize_with_pickle(c, test_eq=True)
            self.serialize_with_pickle(c.to_data_dict, test_eq=True)

    def test_to_data_dict(self):
        comp = Composition('Fe0.00009Ni0.99991')
        d = comp.to_data_dict
        self.assertAlmostEqual(d["reduced_cell_composition"]["Fe"], 9e-5)

    def test_add(self):
        self.assertEqual((self.comp[0] + self.comp[2]).formula,
                         "Li4 Mn2 Fe2 P3 O16",
                         "Incorrect composition after addition!")
        self.assertEqual((self.comp[3] + {"Fe": 4, "O": 4}).formula,
                         "Li4 Fe4 O8", "Incorrect composition after addition!")

    def test_sub(self):
        self.assertEqual((self.comp[0]
                          - Composition("Li2O")).formula,
                         "Li1 Fe2 P3 O11",
                         "Incorrect composition after addition!")
        self.assertEqual((self.comp[0] - {"Fe": 2, "O": 3}).formula,
                         "Li3 P3 O9")

        self.assertRaises(CompositionError, Composition('O').__sub__,
                          Composition('H'))

        # check that S is completely removed by subtraction
        c1 = Composition({'S': 1 + Composition.amount_tolerance / 2, 'O': 1})
        c2 = Composition({'S': 1})
        self.assertEqual(len((c1 - c2).elements), 1)

    def test_mul(self):
        self.assertEqual((self.comp[0] * 4).formula, "Li12 Fe8 P12 O48")
        self.assertEqual((3 * self.comp[1]).formula, "Li9 Fe3 P3 O15")

    def test_div(self):
        self.assertEqual((self.comp[0] / 4).formula, 'Li0.75 Fe0.5 P0.75 O3')

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

    def test_comparisons(self):
        c1 = Composition({'S': 1})
        c1_1 = Composition({'S': 1.00000000000001})
        c2 = Composition({'S': 2})
        c3 = Composition({'O': 1})
        c4 = Composition({'O': 1, 'S': 1})
        self.assertFalse(c1 > c2)
        self.assertFalse(c1_1 > c1)
        self.assertFalse(c1_1 < c1)
        self.assertTrue(c1 > c3)
        self.assertTrue(c3 < c1)
        self.assertTrue(c4 > c1)
        self.assertEqual(sorted([c1, c1_1, c2, c4, c3]),
                         [c3, c1, c1_1, c4, c2])

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

    def test_fractional_composition(self):
        for c in self.comp:
            self.assertAlmostEqual(c.fractional_composition.num_atoms, 1)

    def test_init_numerical_tolerance(self):
        self.assertEqual(Composition({'B': 1, 'C': -1e-12}), Composition('B'))

    def test_negative_compositions(self):
        self.assertEqual(Composition('Li-1(PO-1)4', allow_negative=True).formula,
                         'Li-1 P4 O-4')
        self.assertEqual(Composition('Li-1(PO-1)4', allow_negative=True).reduced_formula,
                         'Li-1(PO-1)4')
        self.assertEqual(Composition('Li-2Mg4', allow_negative=True).reduced_composition,
                         Composition('Li-1Mg2', allow_negative=True))
        self.assertEqual(Composition('Li-2.5Mg4', allow_negative=True).reduced_composition,
                         Composition('Li-2.5Mg4', allow_negative=True))

        # test math
        c1 = Composition('LiCl', allow_negative=True)
        c2 = Composition('Li')
        self.assertEqual(c1 - 2 * c2, Composition({'Li': -1, 'Cl': 1},
                                                  allow_negative=True))
        self.assertEqual((c1 + c2).allow_negative, True)
        self.assertEqual(c1 / -1, Composition('Li-1Cl-1', allow_negative=True))

        # test num_atoms
        c1 = Composition('Mg-1Li', allow_negative=True)
        self.assertEqual(c1.num_atoms, 2)
        self.assertEqual(c1.get_atomic_fraction('Mg'), 0.5)
        self.assertEqual(c1.get_atomic_fraction('Li'), 0.5)
        self.assertEqual(c1.fractional_composition,
                         Composition('Mg-0.5Li0.5', allow_negative=True))

        # test copy
        self.assertEqual(c1.copy(), c1)

        # test species
        c1 = Composition({'Mg': 1, 'Mg2+': -1}, allow_negative=True)
        self.assertEqual(c1.num_atoms, 2)
        self.assertEqual(c1.element_composition, Composition())
        self.assertEqual(c1.average_electroneg, 1.31)

    def test_special_formulas(self):
        special_formulas = {"LiO": "Li2O2", "NaO": "Na2O2", "KO": "K2O2",
                            "HO": "H2O2", "CsO": "Cs2O2", "RbO": "Rb2O2",
                            "O": "O2",  "N": "N2", "F": "F2", "Cl": "Cl2",
                            "H": "H2"}
        for k, v in special_formulas.items():
            self.assertEqual(Composition(k).reduced_formula, v)

    def test_oxi_state_guesses(self):
        self.assertEqual(Composition("LiFeO2").oxi_state_guesses(),
                         ({"Li": 1, "Fe": 3, "O": -2},))

        self.assertEqual(Composition("Fe4O5").oxi_state_guesses(),
                         ({"Fe": 2.5, "O": -2},))

        self.assertEqual(Composition("V2O3").oxi_state_guesses(),
                         ({"V": 3, "O": -2},))

        # all_oxidation_states produces *many* possible responses
        self.assertEqual(len(Composition("MnO").oxi_state_guesses(
            all_oxi_states=True)), 4)

        # can't balance b/c missing V4+
        self.assertEqual(Composition("VO2").oxi_state_guesses(
            oxi_states_override={"V": [2, 3, 5]}), [])

        # missing V4+, but can balance due to additional sites
        self.assertEqual(Composition("V2O4").oxi_state_guesses(
            oxi_states_override={"V": [2, 3, 5]}), ({"V": 4, "O": -2},))

        # multiple solutions - Mn/Fe = 2+/4+ or 3+/3+ or 4+/2+
        self.assertEqual(len(Composition("MnFeO3").oxi_state_guesses(
            oxi_states_override={"Mn": [2, 3, 4], "Fe": [2, 3, 4]})), 3)

        # multiple solutions prefers 3/3 over 2/4 or 4/2
        self.assertEqual(Composition("MnFeO3").oxi_state_guesses(
            oxi_states_override={"Mn": [2, 3, 4], "Fe": [2, 3, 4]})[0],
                         {"Mn": 3, "Fe": 3, "O": -2})

        # target charge of 1
        self.assertEqual(Composition("V2O6").oxi_state_guesses(
            oxi_states_override={"V": [2, 3, 4, 5]}, target_charge=-2),
            ({"V": 5, "O": -2},))

        # max_sites for very large composition - should timeout if incorrect
        self.assertEqual(Composition("Li10000Fe10000P10000O40000").
                         oxi_state_guesses(max_sites=7)[0],
                         {"Li": 1, "Fe": 2, "P": 5, "O": -2})

        # max_sites for very large composition - should timeout if incorrect
        self.assertEqual(Composition("Li10000Fe10000P10000O40000").
                         oxi_state_guesses(max_sites=-1)[0],
                         {"Li": 1, "Fe": 2, "P": 5, "O": -2})

        # negative max_sites less than -1 - should throw error if cannot reduce
        # to under the abs(max_sites) number of sites. Will also timeout if
        # incorrect.
        self.assertEqual(
            Composition("Sb10000O10000F10000").oxi_state_guesses(
                max_sites=-3)[0],
            {"Sb": 3, "O": -2, "F": -1})
        self.assertRaises(ValueError, Composition("LiOF").oxi_state_guesses,
                          max_sites=-2)

        self.assertRaises(ValueError, Composition("V2O3").
                          oxi_state_guesses, max_sites=1)

    def test_oxi_state_decoration(self):
        # Basic test: Get compositions where each element is in a single charge state
        decorated = Composition("H2O").add_charges_from_oxi_state_guesses()
        self.assertIn(Specie("H", 1), decorated)
        self.assertEqual(2, decorated.get(Specie("H", 1)))

        # Test: More than one charge state per element
        decorated = Composition("Fe3O4").add_charges_from_oxi_state_guesses()
        self.assertEqual(1, decorated.get(Specie("Fe", 2)))
        self.assertEqual(2, decorated.get(Specie("Fe", 3)))
        self.assertEqual(4, decorated.get(Specie("O", -2)))

        # Test: No possible charge states
        #   It should return an uncharged composition
        decorated = Composition("NiAl").add_charges_from_oxi_state_guesses()
        self.assertEqual(1, decorated.get(Specie("Ni", 0)))
        self.assertEqual(1, decorated.get(Specie("Al", 0)))

    def test_Metallofullerene(self):
        # Test: Parse Metallofullerene formula (e.g. Y3N@C80)
        formula = "Y3N@C80"
        sym_dict = {"Y": 3, "N": 1, "C": 80}
        cmp = Composition(formula)
        cmp2 = Composition.from_dict(sym_dict)
        self.assertEqual(cmp, cmp2)

    def test_contains_element_type(self):

        formula = "EuTiO3"
        cmp = Composition(formula)
        self.assertTrue(cmp.contains_element_type("lanthanoid"))
        self.assertFalse(cmp.contains_element_type("noble_gas"))
        self.assertTrue(cmp.contains_element_type("f-block"))
        self.assertFalse(cmp.contains_element_type("s-block"))

    def test_chemical_system(self):

        formula = "NaCl"
        cmp = Composition(formula)
        self.assertEqual(cmp.chemical_system, "Cl-Na")

    def test_is_valid(self):

        formula = "NaCl"
        cmp = Composition(formula)
        self.assertTrue(cmp.valid)

        formula = "NaClX"
        cmp = Composition(formula)
        self.assertFalse(cmp.valid)

        self.assertRaises(ValueError,
                          Composition,
                          "NaClX", strict=True)

    def test_remove_charges(self):
        cmp1 = Composition({'Al3+': 2.0, 'O2-': 3.0})

        cmp2 = Composition({'Al': 2.0, 'O': 3.0})
        self.assertNotEqual(str(cmp1), str(cmp2))

        cmp1 = cmp1.remove_charges()
        self.assertEqual(str(cmp1), str(cmp2))

        cmp1 = cmp1.remove_charges()
        self.assertEqual(str(cmp1), str(cmp2))

        cmp1 = Composition({'Fe3+': 2.0, 'Fe2+': 3.0, 'O2-': 6.0})
        cmp2 = Composition({'Fe': 5.0, 'O': 6.0})
        self.assertNotEqual(str(cmp1), str(cmp2))

        cmp1 = cmp1.remove_charges()
        self.assertEqual(str(cmp1), str(cmp2))


class ChemicalPotentialTest(unittest.TestCase):

    def test_init(self):
        d = {'Fe': 1, Element('Fe'): 1}
        self.assertRaises(ValueError, ChemicalPotential, d)
        for k in ChemicalPotential(Fe=1).keys():
            self.assertIsInstance(k, Element)

    def test_math(self):
        fepot = ChemicalPotential({'Fe': 1})
        opot = ChemicalPotential({'O': 2.1})
        pots = ChemicalPotential({'Fe': 1, 'O': 2.1})
        potsx2 = ChemicalPotential({'Fe': 2, 'O': 4.2})
        feo2 = Composition('FeO2')

        # test get_energy()
        self.assertAlmostEqual(pots.get_energy(feo2), 5.2)
        self.assertAlmostEqual(fepot.get_energy(feo2, False), 1)
        self.assertRaises(ValueError, fepot.get_energy, feo2)

        # test multiplication
        self.assertRaises(NotImplementedError, lambda: (pots * pots))
        self.assertDictEqual(pots * 2, potsx2)
        self.assertDictEqual(2 * pots, potsx2)

        # test division
        self.assertDictEqual(potsx2 / 2, pots)
        self.assertRaises(NotImplementedError, lambda: (pots / pots))
        self.assertRaises(NotImplementedError, lambda: (pots / feo2))

        # test add/subtract
        self.assertDictEqual(pots + pots, potsx2)
        self.assertDictEqual(potsx2 - pots, pots)
        self.assertDictEqual(fepot + opot, pots)
        self.assertDictEqual(fepot - opot, pots - opot - opot)


if __name__ == "__main__":
    unittest.main()
