# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest

from pymatgen import Composition
from pymatgen.analysis.reaction_calculator import Reaction, BalancedReaction, \
    ReactionError, ComputedReaction
from pymatgen.entries.computed_entries import ComputedEntry


class ReactionTest(unittest.TestCase):

    def test_init(self):
        reactants = [Composition("Fe"),
                     Composition("O2")]
        products = [Composition("Fe2O3")]
        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn), "2.000 Fe + 1.500 O2 -> 1.000 Fe2O3",
                         "Wrong reaction obtained!")
        self.assertEqual(rxn.normalized_repr, "4 Fe + 3 O2 -> 2 Fe2O3",
                         "Wrong normalized reaction obtained!")

        d = rxn.as_dict()
        rxn = Reaction.from_dict(d)
        self.assertEqual(rxn.normalized_repr, "4 Fe + 3 O2 -> 2 Fe2O3",
                         "Wrong normalized reaction obtained!")

        reactants = [Composition("Fe"), Composition("O"), Composition("Mn"),
                     Composition("P")]
        products = [Composition("FeP"), Composition("MnO")]
        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn),
                         "1.000 Fe + 0.500 O2 + 1.000 Mn + 1.000 P -> 1.000 FeP + 1.000 MnO",
                         "Wrong reaction obtained!")
        self.assertEqual(rxn.normalized_repr,
                         "2 Fe + O2 + 2 Mn + 2 P -> 2 FeP + 2 MnO",
                         "Wrong normalized reaction obtained!")
        reactants = [Composition("FePO4")]
        products = [Composition("FePO4")]

        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn), "1.000 FePO4 -> 1.000 FePO4",
                         "Wrong reaction obtained!")

        products = [Composition("Ti2 O4"), Composition("O1")]
        reactants = [Composition("Ti1 O2")]
        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn), "2.000 TiO2 -> 2.000 TiO2",
                         "Wrong reaction obtained!")

        reactants = [Composition("FePO4"), Composition("Li")]
        products = [Composition("LiFePO4")]
        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn), "1.000 FePO4 + 1.000 Li -> 1.000 LiFePO4",
                         "Wrong reaction obtained!")

        reactants = [Composition("MgO")]
        products = [Composition("MgO")]

        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn), "1.000 MgO -> 1.000 MgO",
                         "Wrong reaction obtained!")

        reactants = [Composition("Mg")]
        products = [Composition("Mg")]

        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn), "1.000 Mg -> 1.000 Mg",
                         "Wrong reaction obtained!")

        reactants = [Composition("FePO4"), Composition("LiPO3")]
        products = [Composition("LiFeP2O7")]

        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn),
                         "1.000 FePO4 + 1.000 LiPO3 -> 1.000 LiFeP2O7",
                         "Wrong reaction obtained!")

        reactants = [Composition("Na"), Composition("K2O")]
        products = [Composition("Na2O"), Composition("K")]
        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn),
                         "1.000 Na + 0.500 K2O -> 0.500 Na2O + 1.000 K",
                         "Wrong reaction obtained!")

        # Test for an old bug which has a problem when excess product is
        # defined.
        products = [Composition("FePO4"), Composition("O")]
        reactants = [Composition("FePO4")]
        rxn = Reaction(reactants, products)

        self.assertEqual(str(rxn), "1.000 FePO4 -> 1.000 FePO4",
                         "Wrong reaction obtained!")

        products = list(map(Composition, ['La8Ti8O12', 'O2', 'LiCrO2']))
        reactants = [Composition('LiLa3Ti3CrO12')]
        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn),
                         "1.000 LiLa3Ti3CrO12 -> 1.500 La2Ti2O3 + 2.750 O2 + 1.000 LiCrO2",
                         "Wrong reaction obtained!")

    def test_equals(self):
        reactants = [Composition("Fe"),
                     Composition("O2")]
        products = [Composition("Fe2O3")]
        rxn = Reaction(reactants, products)
        reactants = [Composition("O2"), Composition("Fe")]
        products = [Composition("Fe2O3")]
        rxn2 = Reaction(reactants, products)
        self.assertTrue(rxn == rxn2)

    def test_normalize_to(self):
        products = [Composition("Fe"), Composition("O2")]
        reactants = [Composition("Fe2O3")]
        rxn = Reaction(reactants, products)
        rxn.normalize_to(Composition("Fe"), 3)
        self.assertEqual(str(rxn), "1.500 Fe2O3 -> 3.000 Fe + 2.250 O2",
                         "Wrong reaction obtained!")

    def test_calculate_energy(self):
        reactants = [Composition("MgO"), Composition("Al2O3")]
        products = [Composition("MgAl2O4")]
        energies = {Composition("MgO"): -0.1, Composition("Al2O3"): -0.2,
                    Composition("MgAl2O4"): -0.5}
        rxn = Reaction(reactants, products)
        self.assertEqual(str(rxn),
                         "1.000 MgO + 1.000 Al2O3 -> 1.000 MgAl2O4")
        self.assertEqual(rxn.normalized_repr, "MgO + Al2O3 -> MgAl2O4")
        self.assertAlmostEqual(rxn.calculate_energy(energies), -0.2, 5)

    def test_as_entry(self):
        reactants = [Composition("MgO"), Composition("Al2O3")]
        products = [Composition("MgAl2O4")]
        energies = {Composition("MgO"): -0.1, Composition("Al2O3"): -0.2,
                    Composition("MgAl2O4"): -0.5}
        rxn = Reaction(reactants, products)
        entry = rxn.as_entry(energies)
        self.assertEqual(entry.name,
                         "1.000 MgO + 1.000 Al2O3 -> 1.000 MgAl2O4")
        self.assertAlmostEqual(entry.energy, -0.2, 5)

        products = [Composition("Fe"), Composition("O2")]
        reactants = [Composition("Fe2O3")]
        rxn = Reaction(reactants, products)
        energies = {Composition("Fe"): 0, Composition("O2"): 0,
                    Composition("Fe2O3"): 0.5}
        entry = rxn.as_entry(energies)
        self.assertEqual(entry.composition.formula, "Fe1.33333333 O2")
        self.assertAlmostEqual(entry.energy, -0.333333, 5)

    def test_products_reactants(self):
        reactants = [Composition("Li3Fe2(PO4)3"), Composition("Fe2O3"),
                     Composition("O2")]
        products = [Composition("LiFePO4")]
        energies = {Composition("Li3Fe2(PO4)3"): -0.1,
                    Composition("Fe2O3"): -0.2,
                    Composition("O2"): -0.2,
                    Composition("LiFePO4"): -0.5}
        rxn = Reaction(reactants, products)

        self.assertIn(Composition("O2"), rxn.products, "O not in products!")
        self.assertIn(Composition("Li3Fe2(PO4)3"), rxn.reactants,
                      "Li3Fe2(PO4)4 not in reactants!")
        self.assertEqual(str(rxn),
                         "0.333 Li3Fe2(PO4)3 + 0.167 Fe2O3 -> 0.250 O2 + 1.000 LiFePO4")
        self.assertEqual(rxn.normalized_repr,
                         "4 Li3Fe2(PO4)3 + 2 Fe2O3 -> 3 O2 + 12 LiFePO4")
        self.assertAlmostEqual(rxn.calculate_energy(energies), -0.48333333, 5)

    def test_to_from_dict(self):
        reactants = [Composition("Fe"), Composition("O2")]
        products = [Composition("Fe2O3")]
        rxn = Reaction(reactants, products)
        d = rxn.as_dict()
        rxn = Reaction.from_dict(d)
        self.assertEqual(rxn.normalized_repr, "4 Fe + 3 O2 -> 2 Fe2O3")


class BalancedReactionTest(unittest.TestCase):
    def test_init(self):
        rct = {Composition('K2SO4'): 3, Composition('Na2S'): 1,
               Composition('Li'): 24}
        prod = {Composition('KNaS'): 2, Composition('K2S'): 2,
                Composition('Li2O'): 12}
        rxn = BalancedReaction(rct, prod)
        self.assertIsNotNone(str(rxn))

        #Test unbalanced exception
        rct = {Composition('K2SO4'): 1,
               Composition('Na2S'): 1, Composition('Li'): 24}
        prod = {Composition('KNaS'): 2, Composition('K2S'): 2,
                Composition('Li2O'): 12}
        self.assertRaises(ReactionError, BalancedReaction, rct, prod)

    def test_to_from_dict(self):
        rct = {Composition('K2SO4'): 3,
               Composition('Na2S'): 1, Composition('Li'): 24}
        prod = {Composition('KNaS'): 2, Composition('K2S'): 2,
                Composition('Li2O'): 12}
        rxn = BalancedReaction(rct, prod)
        d = rxn.as_dict()
        new_rxn = BalancedReaction.from_dict(d)
        for comp in new_rxn.all_comp:
            self.assertEqual(new_rxn.get_coeff(comp), rxn.get_coeff(comp))

    def test_from_string(self):
        rxn = BalancedReaction({Composition("Li"): 4, Composition("O2"): 1},
                               {Composition("Li2O"): 2})
        self.assertEqual(rxn,
                         BalancedReaction.from_string("4 Li + O2 -> 2Li2O"))

        rxn = BalancedReaction({Composition("Li(NiO2)3"): 1}, {Composition("O2"): 0.5,
                               Composition("Li(NiO2)2"): 1, Composition("NiO"): 1})

        self.assertEqual(rxn,
                         BalancedReaction.from_string("1.000 Li(NiO2)3 -> 0.500 O2 + 1.000 Li(NiO2)2 + 1.000 NiO"))

    def test_remove_spectator_species(self):
        rxn = BalancedReaction({Composition("Li"): 4, Composition("O2"): 1, Composition('Na'): 1},
                               {Composition("Li2O"): 2, Composition('Na'): 1})

        self.assertTrue(Composition('Na') not in rxn.all_comp)


class ComputedReactionTest(unittest.TestCase):

    def setUp(self):
        d = [{"correction": 0.0, "data": {}, "energy": -108.56492362,
              "parameters": {}, "composition": {"Li": 54}},
             {"correction": 0.0, "data": {}, "energy": -577.94689128,
              "parameters": {}, "composition": {"O": 32, "Li": 64}},
             {"correction": 0.0, "data": {}, "energy": -17.02844794,
              "parameters": {}, "composition": {"O": 2}},
             {"correction": 0.0, "data": {}, "energy": -959.64693323,
              "parameters": {}, "composition": {"O": 72, "Li": 72}}]
        entries = []
        for e in d:
            entries.append(ComputedEntry.from_dict(e))
        rcts = list(filter(lambda e: e.composition.reduced_formula in ["Li", "O2"],
                      entries))
        prods = list(filter(lambda e: e.composition.reduced_formula == "Li2O2",
                       entries))

        self.rxn = ComputedReaction(rcts, prods)

    def test_calculated_reaction_energy(self):
        self.assertAlmostEqual(self.rxn.calculated_reaction_energy,
                               - 5.60748821935)

    def test_init(self):
        self.assertEqual(str(self.rxn), "1.000 O2 + 2.000 Li -> 1.000 Li2O2")

    def test_to_from_dict(self):
        d = self.rxn.as_dict()
        new_rxn = ComputedReaction.from_dict(d)
        self.assertEqual(str(new_rxn), "1.000 O2 + 2.000 Li -> 1.000 Li2O2")

    def test_all_entries(self):
        for c, e in zip(self.rxn.coeffs, self.rxn.all_entries):
            if c > 0:
                self.assertEqual(e.composition.reduced_formula, "Li2O2")
                self.assertAlmostEqual(e.energy, -959.64693323)


if __name__ == '__main__':
    unittest.main()
