# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os

from numbers import Number

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.phasediagram.maker import PhaseDiagram, GrandPotentialPhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.phasediagram.entries import PDEntryIO, PDEntry


class PDAnalyzerTest(unittest.TestCase):

    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements, entries) = PDEntryIO.from_csv(os.path.join(module_dir,
                                                              "pdentries_test.csv"))
        self.pd = PhaseDiagram(entries)
        self.analyzer = PDAnalyzer(self.pd)

    def test_get_e_above_hull(self):
        for entry in self.pd.stable_entries:
            self.assertLess(self.analyzer.get_e_above_hull(entry), 1e-11,
                            "Stable entries should have e above hull of zero!")
        for entry in self.pd.all_entries:
            if entry not in self.pd.stable_entries:
                e_ah = self.analyzer.get_e_above_hull(entry)
                self.assertGreaterEqual(e_ah, 0)
                self.assertTrue(isinstance(e_ah, Number))

    def test_get_equilibrium_reaction_energy(self):
        for entry in self.pd.stable_entries:
            self.assertLessEqual(
                self.analyzer.get_equilibrium_reaction_energy(entry), 0,
                "Stable entries should have negative equilibrium reaction energy!")

    def test_get_decomposition(self):
        for entry in self.pd.stable_entries:
            self.assertEqual(len(self.analyzer.get_decomposition(entry.composition)), 1,
                              "Stable composition should have only 1 decomposition!")
        dim = len(self.pd.elements)
        for entry in self.pd.all_entries:
            ndecomp = len(self.analyzer.get_decomposition(entry.composition))
            self.assertTrue(ndecomp > 0 and ndecomp <= dim,
                            "The number of decomposition phases can at most be equal to the number of components.")

        #Just to test decomp for a ficitious composition
        ansdict = {entry.composition.formula: amt
                   for entry, amt in
                   self.analyzer.get_decomposition(Composition("Li3Fe7O11")).items()}
        expected_ans = {"Fe2 O2": 0.0952380952380949,
                        "Li1 Fe1 O2": 0.5714285714285714,
                        "Fe6 O8": 0.33333333333333393}
        for k, v in expected_ans.items():
            self.assertAlmostEqual(ansdict[k], v)

    def test_get_transition_chempots(self):
        for el in self.pd.elements:
            self.assertLessEqual(len(self.analyzer.get_transition_chempots(el)),
                                 len(self.pd.facets))

    def test_get_element_profile(self):
        for el in self.pd.elements:
            for entry in self.pd.stable_entries:
                if not (entry.composition.is_element):
                    self.assertLessEqual(len(self.analyzer.get_element_profile(el, entry.composition)),
                                         len(self.pd.facets))

        expected = [{'evolution': 1.0,
                     'chempot': -4.2582781416666666,
                     'reaction': 'Li2O + 0.5 O2 -> Li2O2'},
                    {'evolution': 0,
                     'chempot': -5.0885906699999968,
                     'reaction': 'Li2O -> Li2O'},
                    {'evolution': -1.0,
                     'chempot': -10.487582010000001,
                     'reaction': 'Li2O -> 2 Li + 0.5 O2'}]
        result = self.analyzer.get_element_profile(Element('O'), Composition('Li2O'))
        for d1, d2 in zip(expected, result):
            self.assertAlmostEqual(d1['evolution'], d2['evolution'])
            self.assertAlmostEqual(d1['chempot'], d2['chempot'])
            self.assertEqual(d1['reaction'], str(d2['reaction']))

    def test_get_get_chempot_range_map(self):
        elements = [el for el in self.pd.elements if el.symbol != "Fe"]
        self.assertEqual(len(self.analyzer.get_chempot_range_map(elements)), 10)

    def test_getmu_vertices_stability_phase(self):
        results = self.analyzer.getmu_vertices_stability_phase(Composition("LiFeO2"), Element("O"))
        self.assertAlmostEqual(len(results), 6)
        test_equality = False
        for c in results:
            if abs(c[Element("O")]+7.115) < 1e-2 and abs(c[Element("Fe")]+6.596) < 1e-2 and \
                    abs(c[Element("Li")]+3.931) < 1e-2:
                test_equality = True
        self.assertTrue(test_equality,"there is an expected vertex missing in the list")


    def test_getmu_range_stability_phase(self):
        results = self.analyzer.get_chempot_range_stability_phase(
            Composition("LiFeO2"), Element("O"))
        self.assertAlmostEqual(results[Element("O")][1], -4.4501812249999997)
        self.assertAlmostEqual(results[Element("Fe")][0], -6.5961470999999996)
        self.assertAlmostEqual(results[Element("Li")][0], -3.6250022625000007)

    def test_get_hull_energy(self):
        for entry in self.pd.stable_entries:
            h_e = self.analyzer.get_hull_energy(entry.composition)
            self.assertAlmostEqual(h_e, entry.energy)
            n_h_e = self.analyzer.get_hull_energy(entry.composition.fractional_composition)
            self.assertAlmostEqual(n_h_e, entry.energy_per_atom)

    def test_1d_pd(self):
        entry = PDEntry('H', 0)
        pd = PhaseDiagram([entry])
        pda = PDAnalyzer(pd)
        decomp, e = pda.get_decomp_and_e_above_hull(PDEntry('H', 1))
        self.assertAlmostEqual(e, 1)
        self.assertAlmostEqual(decomp[entry], 1.0)

    def test_get_critical_compositions_fractional(self):
        c1 = Composition('Fe2O3').fractional_composition
        c2 = Composition('Li3FeO4').fractional_composition
        c3 = Composition('Li2O').fractional_composition

        comps = self.analyzer.get_critical_compositions(c1, c2)
        expected = [Composition('Fe2O3').fractional_composition,
                    Composition('Li0.3243244Fe0.1621621O0.51351349'),
                    Composition('Li3FeO4').fractional_composition]
        for crit, exp in zip(comps, expected):
            self.assertTrue(crit.almost_equals(exp, rtol=0, atol=1e-5))

        comps = self.analyzer.get_critical_compositions(c1, c3)
        expected = [Composition('Fe0.4O0.6'),
                    Composition('LiFeO2').fractional_composition,
                    Composition('Li5FeO4').fractional_composition,
                    Composition('Li2O').fractional_composition]
        for crit, exp in zip(comps, expected):
            self.assertTrue(crit.almost_equals(exp, rtol=0, atol=1e-5))

    def test_get_critical_compositions(self):
        c1 = Composition('Fe2O3')
        c2 = Composition('Li3FeO4')
        c3 = Composition('Li2O')

        comps = self.analyzer.get_critical_compositions(c1, c2)
        expected = [Composition('Fe2O3'),
                    Composition('Li0.3243244Fe0.1621621O0.51351349') * 7.4,
                    Composition('Li3FeO4')]
        for crit, exp in zip(comps, expected):
            self.assertTrue(crit.almost_equals(exp, rtol=0, atol=1e-5))

        comps = self.analyzer.get_critical_compositions(c1, c3)
        expected = [Composition('Fe2O3'),
                    Composition('LiFeO2'),
                    Composition('Li5FeO4') / 3,
                    Composition('Li2O')]
        for crit, exp in zip(comps, expected):
            self.assertTrue(crit.almost_equals(exp, rtol=0, atol=1e-5))

        # Don't fail silently if input compositions aren't in phase diagram
        # Can be very confusing if you're working with a GrandPotentialPD
        self.assertRaises(ValueError, self.analyzer.get_critical_compositions,
                          Composition('Xe'), Composition('Mn'))

        # For the moment, should also fail even if compositions are in the gppd
        # because it isn't handled properly
        gppd = GrandPotentialPhaseDiagram(self.pd.all_entries, {'Xe': 1},
                                          self.pd.elements + [Element('Xe')])
        pda = PDAnalyzer(gppd)
        self.assertRaises(ValueError, pda.get_critical_compositions,
                          Composition('Fe2O3'), Composition('Li3FeO4Xe'))

        # check that the function still works though
        comps = pda.get_critical_compositions(c1, c2)
        expected = [Composition('Fe2O3'),
                    Composition('Li0.3243244Fe0.1621621O0.51351349') * 7.4,
                    Composition('Li3FeO4')]
        for crit, exp in zip(comps, expected):
            self.assertTrue(crit.almost_equals(exp, rtol=0, atol=1e-5))

    def test_get_composition_chempots(self):
        c1 = Composition('Fe3.1O4')
        c2 = Composition('Fe3.2O4.1Li0.01')

        e1 = self.analyzer.get_hull_energy(c1)
        e2 = self.analyzer.get_hull_energy(c2)

        cp = self.analyzer.get_composition_chempots(c1)
        calc_e2 = e1 + sum(cp[k] * v for k, v in (c2 - c1).items())
        self.assertAlmostEqual(e2, calc_e2)


if __name__ == '__main__':
    unittest.main()
