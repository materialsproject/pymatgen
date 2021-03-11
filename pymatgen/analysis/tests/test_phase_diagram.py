# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import unittest
import warnings
from numbers import Number
from pathlib import Path
from collections import OrderedDict

import numpy as np

from pymatgen.analysis.phase_diagram import (
    CompoundPhaseDiagram,
    GrandPotentialPhaseDiagram,
    GrandPotPDEntry,
    PDEntry,
    PDPlotter,
    PhaseDiagram,
    ReactionDiagram,
    TransformedPDEntry,
    tet_coord,
    triangular_coord,
    uniquelines,
    BasePhaseDiagram,
)
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import DummySpecies, Element
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.entries.entry_tools import EntrySet

module_dir = Path(__file__).absolute().parent


class PDEntryTest(unittest.TestCase):
    def setUp(self):
        comp = Composition("LiFeO2")
        self.entry = PDEntry(comp, 53)
        self.gpentry = GrandPotPDEntry(self.entry, {Element("O"): 1.5})

    def test_get_energy(self):
        self.assertEqual(self.entry.energy, 53, "Wrong energy!")
        self.assertEqual(self.gpentry.energy, 50, "Wrong energy!")

    def test_get_chemical_energy(self):
        self.assertEqual(self.gpentry.chemical_energy, 3, "Wrong energy!")

    def test_get_energy_per_atom(self):
        self.assertEqual(self.entry.energy_per_atom, 53.0 / 4, "Wrong energy per atom!")
        self.assertEqual(self.gpentry.energy_per_atom, 50.0 / 2, "Wrong energy per atom!")

    def test_get_name(self):
        self.assertEqual(self.entry.name, "LiFeO2", "Wrong name!")
        self.assertEqual(self.gpentry.name, "LiFeO2", "Wrong name!")

    def test_get_composition(self):
        comp = self.entry.composition
        expected_comp = Composition("LiFeO2")
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

        self.assertEqual(entry.name, "LiFeO2", "Wrong name!")
        self.assertEqual(entry.energy_per_atom, 53.0 / 4)
        gpentry = GrandPotPDEntry.from_dict(gpd)
        self.assertEqual(gpentry.name, "LiFeO2", "Wrong name!")
        self.assertEqual(gpentry.energy_per_atom, 50.0 / 2)

        d_anon = d.copy()
        del d_anon["name"]
        try:
            entry = PDEntry.from_dict(d_anon)
        except KeyError:
            self.fail("Should not need to supply name!")

    def test_str(self):
        self.assertIsNotNone(str(self.entry))

    def test_read_csv(self):
        entries = EntrySet.from_csv(str(module_dir / "pdentries_test.csv"))
        self.assertEqual(entries.chemsys, {"Li", "Fe", "O"}, "Wrong elements!")
        self.assertEqual(len(entries), 490, "Wrong number of entries!")


class TransformedPDEntryTest(unittest.TestCase):
    def setUp(self):
        comp = Composition("LiFeO2")
        entry = PDEntry(comp, 53)

        terminal_compositions = ["Li2O", "FeO", "LiO8"]
        terminal_compositions = [Composition(c) for c in terminal_compositions]

        sp_mapping = OrderedDict()
        for i, comp in enumerate(terminal_compositions):
            sp_mapping[comp] = DummySpecies("X" + chr(102 + i))

        self.transformed_entry = TransformedPDEntry(entry, sp_mapping)

    def test_get_energy(self):
        self.assertEqual(self.transformed_entry.energy, 53, "Wrong energy!")
        self.assertAlmostEqual(self.transformed_entry.original_entry.energy, 53.0, 11)

    def test_get_energy_per_atom(self):
        self.assertAlmostEqual(self.transformed_entry.energy_per_atom, 53.0 / (23 / 15), 11)

    def test_get_name(self):
        self.assertEqual(self.transformed_entry.name, "LiFeO2", "Wrong name!")

    def test_get_composition(self):
        comp = self.transformed_entry.composition
        expected_comp = Composition({DummySpecies("Xf"): 14 / 30, DummySpecies("Xg"): 1.0, DummySpecies("Xh"): 2 / 30})
        self.assertEqual(comp, expected_comp, "Wrong composition!")

    def test_is_element(self):
        self.assertFalse(self.transformed_entry.is_element)

    def test_to_from_dict(self):
        d = self.transformed_entry.as_dict()
        entry = TransformedPDEntry.from_dict(d)
        self.assertEqual(entry.name, "LiFeO2", "Wrong name!")
        self.assertAlmostEqual(entry.energy_per_atom, 53.0 / (23 / 15), 11)

    def test_str(self):
        self.assertIsNotNone(str(self.transformed_entry))

    def test_normalize(self):
        norm_entry = self.transformed_entry.normalize(mode="atom")
        expected_comp = Composition(
            {DummySpecies("Xf"): 7 / 23, DummySpecies("Xg"): 15 / 23, DummySpecies("Xh"): 1 / 23}
        )
        self.assertEqual(norm_entry.composition, expected_comp, "Wrong composition!")


class PhaseDiagramTest(unittest.TestCase):
    def setUp(self):
        self.entries = EntrySet.from_csv(str(module_dir / "pdentries_test.csv"))
        self.pd = PhaseDiagram(self.entries)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_init(self):
        # Ensure that a bad set of entries raises a PD error. Remove all Li
        # from self.entries.
        entries = filter(
            lambda e: (not e.composition.is_element) or e.composition.elements[0] != Element("Li"),
            self.entries,
        )
        self.assertRaises(ValueError, PhaseDiagram, entries)

    def test_dim1(self):
        # Ensure that dim 1 PDs can eb generated.
        for el in ["Li", "Fe", "O2"]:
            entries = [e for e in self.entries if e.composition.reduced_formula == el]
            pd = PhaseDiagram(entries)
            self.assertEqual(len(pd.stable_entries), 1)

            for e in entries:
                decomp, ehull = pd.get_decomp_and_e_above_hull(e)
                self.assertGreaterEqual(ehull, 0)
            plotter = PDPlotter(pd)
            lines, stable_entries, unstable_entries = plotter.pd_plot_data
            self.assertEqual(lines[0][1], [0, 0])

    def test_ordering(self):
        # Test sorting of elements
        entries = [ComputedEntry(Composition(formula), 0) for formula in ["O", "N", "Fe"]]
        pd = PhaseDiagram(entries)
        sorted_elements = (Element("Fe"), Element("N"), Element("O"))
        self.assertEqual(tuple(pd.elements), sorted_elements)

        entries.reverse()
        pd = PhaseDiagram(entries)
        self.assertEqual(tuple(pd.elements), sorted_elements)

        # Test manual specification of order
        ordering = [Element(elt_string) for elt_string in ["O", "N", "Fe"]]
        pd = PhaseDiagram(entries, elements=ordering)
        self.assertEqual(tuple(pd.elements), tuple(ordering))

    def test_stable_entries(self):
        stable_formulas = [ent.composition.reduced_formula for ent in self.pd.stable_entries]
        expected_stable = [
            "Fe2O3",
            "Li5FeO4",
            "LiFeO2",
            "Fe3O4",
            "Li",
            "Fe",
            "Li2O",
            "O2",
            "FeO",
        ]
        for formula in expected_stable:
            self.assertTrue(formula in stable_formulas, formula + " not in stable entries!")

    def test_get_formation_energy(self):
        stable_formation_energies = {
            ent.composition.reduced_formula: self.pd.get_form_energy(ent) for ent in self.pd.stable_entries
        }
        expected_formation_energies = {
            "Li5FeO4": -164.8117344866667,
            "Li2O2": -14.119232793333332,
            "Fe2O3": -16.574164339999996,
            "FeO": -5.7141519966666685,
            "Li": 0.0,
            "LiFeO2": -7.732752316666666,
            "Li2O": -6.229303868333332,
            "Fe": 0.0,
            "Fe3O4": -22.565714456666683,
            "Li2FeO3": -45.67166036000002,
            "O2": 0.0,
        }
        for formula, energy in expected_formation_energies.items():
            self.assertAlmostEqual(energy, stable_formation_energies[formula], 7)

    def test_all_entries_hulldata(self):
        self.assertEqual(len(self.pd.all_entries_hulldata), 490)

    def test_planar_inputs(self):
        e1 = PDEntry("H", 0)
        e2 = PDEntry("He", 0)
        e3 = PDEntry("Li", 0)
        e4 = PDEntry("Be", 0)
        e5 = PDEntry("B", 0)
        e6 = PDEntry("Rb", 0)

        pd = PhaseDiagram([e1, e2, e3, e4, e5, e6], map(Element, ["Rb", "He", "B", "Be", "Li", "H"]))

        self.assertEqual(len(pd.facets), 1)

    def test_str(self):
        self.assertIsNotNone(str(self.pd))

    def test_get_e_above_hull(self):
        for entry in self.pd.stable_entries:
            self.assertLess(
                self.pd.get_e_above_hull(entry),
                1e-11,
                "Stable entries should have e above hull of zero!",
            )

        for entry in self.pd.all_entries:
            if entry not in self.pd.stable_entries:
                e_ah = self.pd.get_e_above_hull(entry)
                self.assertTrue(isinstance(e_ah, Number))
                self.assertGreaterEqual(e_ah, 0)

    def test_get_equilibrium_reaction_energy(self):
        for entry in self.pd.stable_entries:
            self.assertLessEqual(
                self.pd.get_equilibrium_reaction_energy(entry),
                0,
                "Stable entries should have negative equilibrium reaction energy!",
            )

    def test_get_phase_separation_energy(self):
        for entry in self.pd.unstable_entries:
            if entry.composition.fractional_composition not in [
                e.composition.fractional_composition for e in self.pd.stable_entries
            ]:
                self.assertGreaterEqual(
                    self.pd.get_phase_separation_energy(entry),
                    0,
                    "Unstable entries should have positive decomposition energy!",
                )
            else:
                if entry.is_element:
                    el_ref = self.pd.el_refs[entry.composition.elements[0]]
                    e_d = entry.energy_per_atom - el_ref.energy_per_atom
                    self.assertAlmostEqual(self.pd.get_phase_separation_energy(entry), e_d, 7)
                # NOTE the remaining materials would require explicit tests as they
                # could be either positive or negative
                pass

        for entry in self.pd.stable_entries:
            if entry.composition.is_element:
                self.assertEqual(
                    self.pd.get_phase_separation_energy(entry),
                    0,
                    "Stable elemental entries should have decomposition energy of zero!",
                )
            else:
                self.assertLessEqual(
                    self.pd.get_phase_separation_energy(entry),
                    0,
                    "Stable entries should have negative decomposition energy!",
                )
                self.assertAlmostEqual(
                    self.pd.get_phase_separation_energy(entry, stable_only=True),
                    self.pd.get_equilibrium_reaction_energy(entry),
                    7,
                    (
                        "Using `stable_only=True` should give decomposition energy equal to "
                        "equilibrium reaction energy!"
                    ),
                )

        # Test that we get correct behaviour with a polymorph
        toy_entries = {
            "Li": 0.0,
            "Li2O": -5,
            "LiO2": -4,
            "O2": 0.0,
        }

        toy_pd = PhaseDiagram([PDEntry(c, e) for c, e in toy_entries.items()])

        # stable entry
        self.assertAlmostEqual(
            toy_pd.get_phase_separation_energy(PDEntry("Li2O", -5)),
            -1.0,
            7,
        )
        # polymorph
        self.assertAlmostEqual(
            toy_pd.get_phase_separation_energy(PDEntry("Li2O", -4)),
            -2.0 / 3.0,
            7,
        )

        # Test that the method works for novel entries
        novel_stable_entry = PDEntry("Li5FeO4", -999)
        self.assertLess(
            self.pd.get_phase_separation_energy(novel_stable_entry),
            0,
            "Novel stable entries should have negative decomposition energy!",
        )

        novel_unstable_entry = PDEntry("Li5FeO4", 999)
        self.assertGreater(
            self.pd.get_phase_separation_energy(novel_unstable_entry),
            0,
            "Novel unstable entries should have positive decomposition energy!",
        )

        duplicate_entry = PDEntry("Li2O", -14.31361175)
        scaled_dup_entry = PDEntry("Li4O2", -14.31361175 * 2)
        stable_entry = [e for e in self.pd.stable_entries if e.name == "Li2O"][0]

        self.assertEqual(
            self.pd.get_phase_separation_energy(duplicate_entry),
            self.pd.get_phase_separation_energy(stable_entry),
            "Novel duplicates of stable entries should have same decomposition energy!",
        )

        self.assertEqual(
            self.pd.get_phase_separation_energy(scaled_dup_entry),
            self.pd.get_phase_separation_energy(stable_entry),
            "Novel scaled duplicates of stable entries should have same decomposition energy!",
        )

    def test_get_decomposition(self):
        for entry in self.pd.stable_entries:
            self.assertEqual(
                len(self.pd.get_decomposition(entry.composition)),
                1,
                "Stable composition should have only 1 decomposition!",
            )
        dim = len(self.pd.elements)
        for entry in self.pd.all_entries:
            ndecomp = len(self.pd.get_decomposition(entry.composition))
            self.assertTrue(
                ndecomp > 0 and ndecomp <= dim,
                "The number of decomposition phases can at most be equal to the number of components.",
            )

        # Just to test decomp for a ficitious composition
        ansdict = {
            entry.composition.formula: amt for entry, amt in self.pd.get_decomposition(Composition("Li3Fe7O11")).items()
        }
        expected_ans = {
            "Fe2 O2": 0.0952380952380949,
            "Li1 Fe1 O2": 0.5714285714285714,
            "Fe6 O8": 0.33333333333333393,
        }
        for k, v in expected_ans.items():
            self.assertAlmostEqual(ansdict[k], v, 7)

    def test_get_transition_chempots(self):
        for el in self.pd.elements:
            self.assertLessEqual(len(self.pd.get_transition_chempots(el)), len(self.pd.facets))

    def test_get_element_profile(self):
        for el in self.pd.elements:
            for entry in self.pd.stable_entries:
                if not (entry.composition.is_element):
                    self.assertLessEqual(
                        len(self.pd.get_element_profile(el, entry.composition)),
                        len(self.pd.facets),
                    )

        expected = [
            {
                "evolution": 1.0,
                "chempot": -4.2582781416666666,
                "reaction": "Li2O + 0.5 O2 -> Li2O2",
            },
            {
                "evolution": 0,
                "chempot": -5.0885906699999968,
                "reaction": "Li2O -> Li2O",
            },
            {
                "evolution": -1.0,
                "chempot": -10.487582010000001,
                "reaction": "Li2O -> 2 Li + 0.5 O2",
            },
        ]
        result = self.pd.get_element_profile(Element("O"), Composition("Li2O"))
        for d1, d2 in zip(expected, result):
            self.assertAlmostEqual(d1["evolution"], d2["evolution"])
            self.assertAlmostEqual(d1["chempot"], d2["chempot"])
            self.assertEqual(d1["reaction"], str(d2["reaction"]))

    def test_get_get_chempot_range_map(self):
        elements = [el for el in self.pd.elements if el.symbol != "Fe"]
        self.assertEqual(len(self.pd.get_chempot_range_map(elements)), 10)

    def test_getmu_vertices_stability_phase(self):
        results = self.pd.getmu_vertices_stability_phase(Composition("LiFeO2"), Element("O"))
        self.assertAlmostEqual(len(results), 6)
        test_equality = False
        for c in results:
            if (
                abs(c[Element("O")] + 7.115) < 1e-2
                and abs(c[Element("Fe")] + 6.596) < 1e-2
                and abs(c[Element("Li")] + 3.931) < 1e-2
            ):
                test_equality = True
        self.assertTrue(test_equality, "there is an expected vertex missing in the list")

    def test_getmu_range_stability_phase(self):
        results = self.pd.get_chempot_range_stability_phase(Composition("LiFeO2"), Element("O"))
        self.assertAlmostEqual(results[Element("O")][1], -4.4501812249999997)
        self.assertAlmostEqual(results[Element("Fe")][0], -6.5961470999999996)
        self.assertAlmostEqual(results[Element("Li")][0], -3.6250022625000007)

    def test_get_hull_energy(self):
        for entry in self.pd.stable_entries:
            h_e = self.pd.get_hull_energy(entry.composition)
            self.assertAlmostEqual(h_e, entry.energy)
            n_h_e = self.pd.get_hull_energy(entry.composition.fractional_composition)
            self.assertAlmostEqual(n_h_e, entry.energy_per_atom)

    def test_1d_pd(self):
        entry = PDEntry("H", 0)
        pd = PhaseDiagram([entry])
        decomp, e = pd.get_decomp_and_e_above_hull(PDEntry("H", 1))
        self.assertAlmostEqual(e, 1)
        self.assertAlmostEqual(decomp[entry], 1.0)

    def test_get_critical_compositions_fractional(self):
        c1 = Composition("Fe2O3").fractional_composition
        c2 = Composition("Li3FeO4").fractional_composition
        c3 = Composition("Li2O").fractional_composition

        comps = self.pd.get_critical_compositions(c1, c2)
        expected = [
            Composition("Fe2O3").fractional_composition,
            Composition("Li0.3243244Fe0.1621621O0.51351349"),
            Composition("Li3FeO4").fractional_composition,
        ]
        for crit, exp in zip(comps, expected):
            self.assertTrue(crit.almost_equals(exp, rtol=0, atol=1e-5))

        comps = self.pd.get_critical_compositions(c1, c3)
        expected = [
            Composition("Fe0.4O0.6"),
            Composition("LiFeO2").fractional_composition,
            Composition("Li5FeO4").fractional_composition,
            Composition("Li2O").fractional_composition,
        ]
        for crit, exp in zip(comps, expected):
            self.assertTrue(crit.almost_equals(exp, rtol=0, atol=1e-5))

    def test_get_critical_compositions(self):
        c1 = Composition("Fe2O3")
        c2 = Composition("Li3FeO4")
        c3 = Composition("Li2O")

        comps = self.pd.get_critical_compositions(c1, c2)
        expected = [
            Composition("Fe2O3"),
            Composition("Li0.3243244Fe0.1621621O0.51351349") * 7.4,
            Composition("Li3FeO4"),
        ]
        for crit, exp in zip(comps, expected):
            self.assertTrue(crit.almost_equals(exp, rtol=0, atol=1e-5))

        comps = self.pd.get_critical_compositions(c1, c3)
        expected = [
            Composition("Fe2O3"),
            Composition("LiFeO2"),
            Composition("Li5FeO4") / 3,
            Composition("Li2O"),
        ]
        for crit, exp in zip(comps, expected):
            self.assertTrue(crit.almost_equals(exp, rtol=0, atol=1e-5))

        # Don't fail silently if input compositions aren't in phase diagram
        # Can be very confusing if you're working with a GrandPotentialPD
        self.assertRaises(
            ValueError,
            self.pd.get_critical_compositions,
            Composition("Xe"),
            Composition("Mn"),
        )

        # For the moment, should also fail even if compositions are in the gppd
        # because it isn't handled properly
        gppd = GrandPotentialPhaseDiagram(self.pd.all_entries, {"Xe": 1}, self.pd.elements + [Element("Xe")])
        self.assertRaises(
            ValueError,
            gppd.get_critical_compositions,
            Composition("Fe2O3"),
            Composition("Li3FeO4Xe"),
        )

        # check that the function still works though
        comps = gppd.get_critical_compositions(c1, c2)
        expected = [
            Composition("Fe2O3"),
            Composition("Li0.3243244Fe0.1621621O0.51351349") * 7.4,
            Composition("Li3FeO4"),
        ]
        for crit, exp in zip(comps, expected):
            self.assertTrue(crit.almost_equals(exp, rtol=0, atol=1e-5))

        # case where the endpoints are identical
        self.assertEqual(self.pd.get_critical_compositions(c1, c1 * 2), [c1, c1 * 2])

    def test_get_composition_chempots(self):
        c1 = Composition("Fe3.1O4")
        c2 = Composition("Fe3.2O4.1Li0.01")

        e1 = self.pd.get_hull_energy(c1)
        e2 = self.pd.get_hull_energy(c2)

        cp = self.pd.get_composition_chempots(c1)
        calc_e2 = e1 + sum(cp[k] * v for k, v in (c2 - c1).items())
        self.assertAlmostEqual(e2, calc_e2)

    def test_get_all_chempots(self):
        c1 = Composition("Fe3.1O4")
        c2 = Composition("FeO")

        cp1 = self.pd.get_all_chempots(c1)
        cpresult = {
            Element("Li"): -4.077061954999998,
            Element("Fe"): -6.741593864999999,
            Element("O"): -6.969907375000003,
        }

        for elem, energy in cpresult.items():
            self.assertAlmostEqual(cp1["Fe3O4-FeO-LiFeO2"][elem], energy)

        cp2 = self.pd.get_all_chempots(c2)
        cpresult = {
            Element("O"): -7.115354140000001,
            Element("Fe"): -6.5961471,
            Element("Li"): -3.9316151899999987,
        }

        for elem, energy in cpresult.items():
            self.assertAlmostEqual(cp2["FeO-LiFeO2-Fe"][elem], energy)

    def test_to_from_dict(self):

        # test round-trip for other entry types such as ComputedEntry
        entry = ComputedEntry("H", 0.0, 0.0, entry_id="test")
        pd = PhaseDiagram([entry])
        d = pd.as_dict()
        pd_roundtrip = PhaseDiagram.from_dict(d)
        self.assertEqual(pd.all_entries[0].entry_id, pd_roundtrip.all_entries[0].entry_id)


class GrandPotentialPhaseDiagramTest(unittest.TestCase):
    def setUp(self):
        self.entries = EntrySet.from_csv(str(module_dir / "pdentries_test.csv"))
        self.pd = GrandPotentialPhaseDiagram(self.entries, {Element("O"): -5})
        self.pd6 = GrandPotentialPhaseDiagram(self.entries, {Element("O"): -6})

    def test_stable_entries(self):
        stable_formulas = [ent.original_entry.composition.reduced_formula for ent in self.pd.stable_entries]
        expected_stable = ["Li5FeO4", "Li2FeO3", "LiFeO2", "Fe2O3", "Li2O2"]
        for formula in expected_stable:
            self.assertTrue(formula in stable_formulas, "{} not in stable entries!".format(formula))
        self.assertEqual(len(self.pd6.stable_entries), 4)

    def test_get_formation_energy(self):
        stable_formation_energies = {
            ent.original_entry.composition.reduced_formula: self.pd.get_form_energy(ent)
            for ent in self.pd.stable_entries
        }
        expected_formation_energies = {
            "Fe2O3": 0.0,
            "Li5FeO4": -5.305515040000046,
            "Li2FeO3": -2.3424741500000152,
            "LiFeO2": -0.43026396250000154,
            "Li2O2": 0.0,
        }
        for formula, energy in expected_formation_energies.items():
            self.assertAlmostEqual(
                energy,
                stable_formation_energies[formula],
                7,
                "Calculated formation for {} is not correct!".format(formula),
            )

    def test_str(self):
        self.assertIsNotNone(str(self.pd))


class BasePhaseDiagramTest(PhaseDiagramTest):
    def setUp(self):
        self.entries = EntrySet.from_csv(str(module_dir / "pdentries_test.csv"))
        self.pd = BasePhaseDiagram.from_entries(self.entries)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_init(self):
        pass

    def test_as_dict_from_dict(self):
        dd = self.pd.as_dict()
        new_pd = BasePhaseDiagram.from_dict(dd)
        new_dd = new_pd.as_dict()
        self.assertEqual(new_dd, dd)


class CompoundPhaseDiagramTest(unittest.TestCase):
    def setUp(self):
        self.entries = EntrySet.from_csv(str(module_dir / "pdentries_test.csv"))
        self.pd = CompoundPhaseDiagram(self.entries, [Composition("Li2O"), Composition("Fe2O3")])

    def test_stable_entries(self):
        stable_formulas = [ent.name for ent in self.pd.stable_entries]
        expected_stable = ["Fe2O3", "Li5FeO4", "LiFeO2", "Li2O"]
        for formula in expected_stable:
            self.assertTrue(formula in stable_formulas)

    def test_get_formation_energy(self):
        stable_formation_energies = {ent.name: self.pd.get_form_energy(ent) for ent in self.pd.stable_entries}
        expected_formation_energies = {
            "Li5FeO4": -7.0773284399999739,
            "Fe2O3": 0,
            "LiFeO2": -0.47455929750000081,
            "Li2O": 0,
        }
        for formula, energy in expected_formation_energies.items():
            self.assertAlmostEqual(energy, stable_formation_energies[formula], 7)

    def test_str(self):
        self.assertIsNotNone(str(self.pd))


class ReactionDiagramTest(unittest.TestCase):
    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        self.entries = list(EntrySet.from_csv(os.path.join(module_dir, "reaction_entries_test.csv")).entries)
        for e in self.entries:
            if e.composition.reduced_formula == "VPO5":
                entry1 = e
            elif e.composition.reduced_formula == "H4(CO)3":
                entry2 = e
        self.rd = ReactionDiagram(entry1=entry1, entry2=entry2, all_entries=self.entries[2:])

    def test_get_compound_pd(self):
        self.rd.get_compound_pd()

    def test_formula(self):
        for e in self.rd.rxn_entries:
            self.assertIn(Element.V, e.composition)
            self.assertIn(Element.O, e.composition)
            self.assertIn(Element.C, e.composition)
            self.assertIn(Element.P, e.composition)
            self.assertIn(Element.H, e.composition)
        # formed_formula = [e.composition.reduced_formula for e in
        #                   self.rd.rxn_entries]
        # expected_formula = [
        #     'V0.12707182P0.12707182H0.0441989C0.03314917O0.66850829',
        #     'V0.125P0.125H0.05C0.0375O0.6625',
        #     'V0.12230216P0.12230216H0.05755396C0.04316547O0.65467626',
        #     'V0.11340206P0.11340206H0.08247423C0.06185567O0.62886598',
        #     'V0.11267606P0.11267606H0.08450704C0.06338028O0.62676056',
        #     'V0.11229947P0.11229947H0.0855615C0.06417112O0.62566845',
        #     'V0.09677419P0.09677419H0.12903226C0.09677419O0.58064516',
        #     'V0.05882353P0.05882353H0.23529412C0.17647059O0.47058824',
        #     'V0.04225352P0.04225352H0.28169014C0.21126761O0.42253521']
        #
        # for formula in expected_formula:
        #     self.assertTrue(formula in formed_formula, "%s not in %s" % (formed_formula, expected_formula))


class PDPlotterTest(unittest.TestCase):
    def setUp(self):
        entries = list(EntrySet.from_csv(os.path.join(module_dir, "pdentries_test.csv")))

        self.pd_ternary = PhaseDiagram(entries)
        self.plotter_ternary_mpl = PDPlotter(self.pd_ternary, backend="matplotlib")
        self.plotter_ternary_plotly = PDPlotter(self.pd_ternary, backend="plotly")

        entrieslio = [e for e in entries if "Fe" not in e.composition]
        self.pd_binary = PhaseDiagram(entrieslio)
        self.plotter_binary_mpl = PDPlotter(self.pd_binary, backend="matplotlib")
        self.plotter_binary_plotly = PDPlotter(self.pd_binary, backend="plotly")

        entries.append(PDEntry("C", 0))
        self.pd_quaternary = PhaseDiagram(entries)
        self.plotter_quaternary_mpl = PDPlotter(self.pd_quaternary, backend="matplotlib")
        self.plotter_quaternary_plotly = PDPlotter(self.pd_quaternary, backend="plotly")

    def test_pd_plot_data(self):
        (lines, labels, unstable_entries) = self.plotter_ternary_mpl.pd_plot_data
        self.assertEqual(len(lines), 22)
        self.assertEqual(
            len(labels),
            len(self.pd_ternary.stable_entries),
            "Incorrect number of lines generated!",
        )
        self.assertEqual(
            len(unstable_entries),
            len(self.pd_ternary.all_entries) - len(self.pd_ternary.stable_entries),
            "Incorrect number of lines generated!",
        )
        (lines, labels, unstable_entries) = self.plotter_quaternary_mpl.pd_plot_data
        self.assertEqual(len(lines), 33)
        self.assertEqual(len(labels), len(self.pd_quaternary.stable_entries))
        self.assertEqual(
            len(unstable_entries),
            len(self.pd_quaternary.all_entries) - len(self.pd_quaternary.stable_entries),
        )
        (lines, labels, unstable_entries) = self.plotter_binary_mpl.pd_plot_data
        self.assertEqual(len(lines), 3)
        self.assertEqual(len(labels), len(self.pd_binary.stable_entries))

    def test_mpl_plots(self):
        # Some very basic ("non")-tests. Just to make sure the methods are callable.
        self.plotter_binary_mpl.get_plot().close()
        self.plotter_ternary_mpl.get_plot().close()
        self.plotter_quaternary_mpl.get_plot().close()
        self.plotter_ternary_mpl.get_contour_pd_plot().close()
        self.plotter_ternary_mpl.get_chempot_range_map_plot([Element("Li"), Element("O")]).close()
        self.plotter_ternary_mpl.plot_element_profile(Element("O"), Composition("Li2O")).close()

    def test_plotly_plots(self):
        # Also very basic tests. Ensures callability and 2D vs 3D properties.
        self.plotter_binary_plotly.get_plot()
        self.plotter_ternary_plotly.get_plot()
        self.plotter_quaternary_plotly.get_plot()


class UtilityFunctionTest(unittest.TestCase):
    def test_unique_lines(self):
        testdata = [
            [5, 53, 353],
            [399, 20, 52],
            [399, 400, 20],
            [13, 399, 52],
            [21, 400, 353],
            [393, 5, 353],
            [400, 393, 353],
            [393, 400, 399],
            [393, 13, 5],
            [13, 393, 399],
            [400, 17, 20],
            [21, 17, 400],
        ]
        expected_ans = {
            (5, 393),
            (21, 353),
            (353, 400),
            (5, 13),
            (17, 20),
            (21, 400),
            (17, 400),
            (52, 399),
            (393, 399),
            (20, 52),
            (353, 393),
            (5, 353),
            (5, 53),
            (13, 399),
            (393, 400),
            (13, 52),
            (53, 353),
            (17, 21),
            (13, 393),
            (20, 399),
            (399, 400),
            (20, 400),
        }
        self.assertEqual(uniquelines(testdata), expected_ans)

    def test_triangular_coord(self):
        coord = [0.5, 0.5]
        coord = triangular_coord(coord)
        self.assertTrue(np.allclose(coord, [0.75, 0.4330127]))

    def test_tet_coord(self):
        coord = [0.5, 0.5, 0.5]
        coord = tet_coord(coord)
        self.assertTrue(np.allclose(coord, [1.0, 0.57735027, 0.40824829]))


if __name__ == "__main__":
    unittest.main()
