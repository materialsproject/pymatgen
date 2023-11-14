from __future__ import annotations

import collections
import os
import unittest
from numbers import Number

import numpy as np
import pytest
from monty.serialization import dumpfn, loadfn
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.phase_diagram import (
    CompoundPhaseDiagram,
    GrandPotentialPhaseDiagram,
    GrandPotPDEntry,
    PatchedPhaseDiagram,
    PDEntry,
    PDPlotter,
    PhaseDiagram,
    ReactionDiagram,
    TransformedPDEntry,
    tet_coord,
    triangular_coord,
    uniquelines,
)
from pymatgen.core import Composition, DummySpecies, Element
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.entries.entry_tools import EntrySet
from pymatgen.util.testing import PymatgenTest

module_dir = os.path.dirname(os.path.abspath(__file__))


class TestPDEntry(unittest.TestCase):
    def setUp(self):
        comp = Composition("LiFeO2")
        self.entry = PDEntry(comp, 53, name="mp-757614")
        self.gp_entry = GrandPotPDEntry(self.entry, {Element("O"): 1.5})

    def test_get_energy(self):
        assert self.entry.energy == 53, "Wrong energy!"
        assert self.gp_entry.energy == 50, "Wrong energy!"

    def test_get_chemical_energy(self):
        assert self.gp_entry.chemical_energy == 3, "Wrong energy!"

    def test_get_energy_per_atom(self):
        assert self.entry.energy_per_atom == 53.0 / 4, "Wrong energy per atom!"
        assert self.gp_entry.energy_per_atom == 50.0 / 2, "Wrong energy per atom!"

    def test_get_name(self):
        assert self.entry.name == "mp-757614"
        assert self.gp_entry.name == "mp-757614"

    def test_composition(self):
        comp = self.entry.composition
        expected_comp = Composition("LiFeO2")
        assert comp == expected_comp
        comp = self.gp_entry.composition
        expected_comp = Composition("LiFe")
        assert comp == expected_comp

    def test_elements(self):
        expected_elements = list(map(Element, ["Li", "Fe", "O"]))
        assert self.entry.elements == expected_elements

    def test_is_element(self):
        assert not self.entry.is_element
        assert not self.gp_entry.is_element

    def test_as_from_dict(self):
        d = self.entry.as_dict()
        gpd = self.gp_entry.as_dict()
        entry = PDEntry.from_dict(d)

        assert entry.name == "mp-757614"
        assert entry.energy_per_atom == 53.0 / 4
        gpentry = GrandPotPDEntry.from_dict(gpd)
        assert gpentry.name == "mp-757614"
        assert gpentry.energy_per_atom == 50.0 / 2

        d_anon = d.copy()
        del d_anon["name"]
        try:
            entry = PDEntry.from_dict(d_anon)
        except KeyError:
            self.fail("Should not need to supply name!")

    def test_str(self):
        assert str(self.entry) == "PDEntry : Li1 Fe1 O2 (mp-757614) with energy = 53.0000"
        pde = self.entry.as_dict()
        del pde["name"]
        pde = PDEntry.from_dict(pde)
        assert str(pde) == "PDEntry : Li1 Fe1 O2 with energy = 53.0000"

    def test_read_csv(self):
        entries = EntrySet.from_csv(f"{module_dir}/pd_entries_test.csv")
        assert entries.chemsys == {"Li", "Fe", "O"}, "Wrong elements!"
        assert len(entries) == 490, "Wrong number of entries!"


class TestTransformedPDEntry(unittest.TestCase):
    def setUp(self):
        comp = Composition("LiFeO2")
        entry = PDEntry(comp, 53)

        terminal_compositions = ["Li2O", "FeO", "LiO8"]
        terminal_compositions = [Composition(c) for c in terminal_compositions]

        sp_mapping = {}
        for idx, comp in enumerate(terminal_compositions):
            sp_mapping[comp] = DummySpecies("X" + chr(102 + idx))

        self.transformed_entry = TransformedPDEntry(entry, sp_mapping)

    def test_get_energy(self):
        assert self.transformed_entry.energy == 53, "Wrong energy!"
        assert self.transformed_entry.original_entry.energy == approx(53.0)

    def test_get_energy_per_atom(self):
        assert self.transformed_entry.energy_per_atom == approx(53.0 / (23 / 15))

    def test_get_name(self):
        assert self.transformed_entry.name == "LiFeO2"

    def test_composition(self):
        comp = self.transformed_entry.composition
        expected_comp = Composition({DummySpecies("Xf"): 14 / 30, DummySpecies("Xg"): 1.0, DummySpecies("Xh"): 2 / 30})
        assert comp == expected_comp

    def test_elements(self):
        expected_elements = list(map(Element, ["Li", "Fe", "O"]))
        assert self.transformed_entry.elements == expected_elements

    def test_is_element(self):
        assert self.transformed_entry.is_element is False
        assert self.transformed_entry.original_entry.is_element is False
        iron = Composition("Fe")
        assert TransformedPDEntry(PDEntry(iron, 0), {iron: iron}).is_element is True

    def test_as_from_dict(self):
        dct = self.transformed_entry.as_dict()
        entry = TransformedPDEntry.from_dict(dct)
        assert entry.name == "LiFeO2" == self.transformed_entry.name
        assert entry.energy_per_atom == approx(53.0 / (23 / 15)) == self.transformed_entry.energy_per_atom

    def test_str(self):
        assert str(self.transformed_entry).startswith("TransformedPDEntry Xf0+0.46666667 Xg")
        assert str(self.transformed_entry).endswith("with original composition Li1 Fe1 O2, energy = 53.0000")

    def test_normalize(self):
        norm_entry = self.transformed_entry.normalize(mode="atom")
        expected_comp = Composition(
            {DummySpecies("Xf"): 7 / 23, DummySpecies("Xg"): 15 / 23, DummySpecies("Xh"): 1 / 23}
        )
        assert norm_entry.composition == expected_comp


class TestPhaseDiagram(PymatgenTest):
    def setUp(self):
        self.entries = EntrySet.from_csv(f"{module_dir}/pd_entries_test.csv")
        self.pd = PhaseDiagram(self.entries)

    def test_init(self):
        # Ensure that a bad set of entries raises a PD error. Remove all Li
        # from self.entries.
        entries = filter(
            lambda e: (not e.composition.is_element) or e.elements[0] != Element("Li"),
            self.entries,
        )
        with pytest.raises(ValueError, match=r"Missing terminal entries for elements \['Fe', 'Li', 'O'\]"):
            PhaseDiagram(entries)

    def test_repr(self):
        assert (
            str(self.pd) == "Li-Fe-O phase diagram\n11 stable phases: \nFe, FeO, "
            "Fe2O3, Fe3O4, LiFeO2, Li, Li2O, LiO, Li5FeO4, Li2FeO3, O"
        )

    def test_dim1(self):
        # Ensure that dim 1 PDs can be generated.
        for el in ("Li", "Fe", "O2"):
            entries = [entry for entry in self.entries if entry.composition.reduced_formula == el]
            pd = PhaseDiagram(entries)
            assert len(pd.stable_entries) == 1

            for entry in entries:
                e_hull = pd.get_e_above_hull(entry)
                assert e_hull >= 0

            plotter = PDPlotter(pd)
            lines, *_ = plotter.pd_plot_data
            assert lines[0][1] == [0, 0]

    def test_ordering(self):
        # Test sorting of elements
        entries = [ComputedEntry(Composition(formula), 0) for formula in ["O", "N", "Fe"]]
        pd = PhaseDiagram(entries)
        sorted_elements = (Element("Fe"), Element("N"), Element("O"))
        assert tuple(pd.elements) == sorted_elements

        entries.reverse()
        pd = PhaseDiagram(entries)
        assert tuple(pd.elements) == sorted_elements

        # Test manual specification of order
        ordering = [Element(elt_string) for elt_string in ["O", "N", "Fe"]]
        pd = PhaseDiagram(entries, elements=ordering)
        assert tuple(pd.elements) == tuple(ordering)

    def test_stable_entries(self):
        stable_formulas = [ent.composition.reduced_formula for ent in self.pd.stable_entries]
        expected_stable = "Fe2O3 Li5FeO4 LiFeO2 Fe3O4 Li Fe Li2O O2 FeO".split()
        for formula in expected_stable:
            assert formula in stable_formulas, f"{formula} not in stable entries!"

    def test_get_form_energy(self):
        expected_formation_energies = {
            "Li5FeO4": -164.8117344,
            "Li2O2": -14.119232793,
            "Fe2O3": -16.57416433,
            "FeO": -5.71415199,
            "Li": 0.0,
            "LiFeO2": -7.73275231,
            "Li2O": -6.229303868,
            "Fe": 0.0,
            "Fe3O4": -22.56571445,
            "Li2FeO3": -45.67166036,
            "O2": 0.0,
        }

        for entry in self.pd.stable_entries:
            formula = entry.composition.reduced_formula
            expected = expected_formation_energies[formula]
            n_atoms = entry.composition.num_atoms
            # test get_form_energy
            assert self.pd.get_form_energy(entry) == approx(expected), formula
            # test get_form_energy_per_atom
            assert self.pd.get_form_energy_per_atom(entry) == approx(expected / n_atoms), formula

    def test_get_reference_energy(self):
        expected_ref_energies = {
            "Li5FeO4": -265.5546721,
            "Li2O2": -24.685172046,
            "Fe2O3": -51.93425725,
            "FeO": -21.70885048,
            "Li": -1.91301487,
            "LiFeO2": -17.02571825,
            "Li2O": -8.084307881,
            "Fe": -6.5961471,
            "Fe3O4": -73.6431077,
            "Li2FeO3": -92.78804506,
            "O2": -25.54966885,
        }
        for entry in self.pd.stable_entries:
            formula = entry.composition.reduced_formula
            actual = self.pd.get_reference_energy(entry.composition)
            expected = expected_ref_energies[formula]
            assert actual == approx(expected), formula

    def test_all_entries_hulldata(self):
        assert len(self.pd.all_entries_hulldata) == 490

    def test_planar_inputs(self):
        elems = ["H", "He", "Li", "Be", "B", "Rb"]
        e1, e2, e3, e4, e5, e6 = (PDEntry(elem, 0) for elem in elems)

        pd = PhaseDiagram([e1, e2, e3, e4, e5, e6], map(Element, elems))

        assert len(pd.facets) == 1

    def test_str(self):
        assert (
            str(self.pd) == "Li-Fe-O phase diagram\n11 stable phases: \nFe, FeO, Fe2O3, Fe3O4,"
            " LiFeO2, Li, Li2O, LiO, Li5FeO4, Li2FeO3, O"
        )

    def test_get_e_above_hull(self):
        for entry in self.pd.all_entries:
            for entry in self.pd.stable_entries:
                decomp, e_hull = self.pd.get_decomp_and_e_above_hull(entry)
                assert e_hull < 1e-11, "Stable entries should have e_above_hull of zero!"
                assert decomp[entry] == 1, "Decomposition of stable entry should be itself."

            e_ah = self.pd.get_e_above_hull(entry)
            assert isinstance(e_ah, Number)
            assert e_ah >= 0

    def test_get_decomp_and_e_above_hull_on_error(self):
        for method, expected in (
            (self.pd.get_e_above_hull, None),
            (self.pd.get_decomp_and_e_above_hull, (None, None)),
        ):
            # test raises ValueError on entry with element not in the phase diagram
            U_entry = PDEntry("U", 0)
            with pytest.raises(ValueError, match="Unable to get decomposition for PDEntry : U1 with energy"):
                method(U_entry)

            # test raises ValueError on entry with very negative energy
            too_neg_entry = PDEntry("Li", -1e6)
            match_msg = "No valid decomposition found for PDEntry : Li1 with energy"
            with pytest.raises(ValueError, match=match_msg):
                method(too_neg_entry)

            with pytest.warns(UserWarning, match=match_msg):
                out = method(too_neg_entry, on_error="warn")
                assert out == expected

            out = method(too_neg_entry, on_error="ignore")
            assert out == expected

    def test_downstream_methods_can_also_ignore_errors(self):
        # test that downstream methods get_e_above_hull() and get_phase_separation_energy()
        # can also ignore errors
        too_neg_entry = PDEntry("Li", -1e6)
        exotic_entry = PDEntry("U W", -1e6)

        # get_e_above_hull
        with pytest.raises(ValueError, match="No valid decomposition found for PDEntry "):
            self.pd.get_e_above_hull(too_neg_entry, on_error="raise")
        assert self.pd.get_e_above_hull(too_neg_entry, on_error="ignore") is None

        with pytest.raises(ValueError, match="Unable to get decomposition for PDEntry"):
            self.pd.get_e_above_hull(exotic_entry, on_error="raise")
        assert self.pd.get_e_above_hull(exotic_entry, on_error="ignore") is None

        # get_phase_separation_energy
        with pytest.raises(ValueError, match="Unable to get decomposition for PDEntry"):
            self.pd.get_phase_separation_energy(exotic_entry, on_error="raise")
        assert self.pd.get_phase_separation_energy(exotic_entry, on_error="ignore") is None

    def test_get_equilibrium_reaction_energy(self):
        for entry in self.pd.stable_entries:
            assert (
                self.pd.get_equilibrium_reaction_energy(entry) <= 0
            ), "Stable entries should have negative equilibrium reaction energy!"

    def test_get_phase_separation_energy(self):
        for entry in self.pd.unstable_entries:
            if entry.composition.fractional_composition not in [
                e.composition.fractional_composition for e in self.pd.stable_entries
            ]:
                assert (
                    self.pd.get_phase_separation_energy(entry) >= 0
                ), "Unstable entries should have positive decomposition energy!"
            elif entry.is_element:
                el_ref = self.pd.el_refs[entry.elements[0]]
                e_d = entry.energy_per_atom - el_ref.energy_per_atom
                assert self.pd.get_phase_separation_energy(entry) == approx(e_d)
            # NOTE the remaining materials would require explicit tests as they
            # could be either positive or negative

        for entry in self.pd.stable_entries:
            if entry.composition.is_element:
                assert (
                    self.pd.get_phase_separation_energy(entry) == 0
                ), "Stable elemental entries should have decomposition energy of zero!"
            else:
                assert (
                    self.pd.get_phase_separation_energy(entry) <= 0
                ), "Stable entries should have negative decomposition energy!"
                assert self.pd.get_phase_separation_energy(entry, stable_only=True) == approx(
                    self.pd.get_equilibrium_reaction_energy(entry)
                ), "Using `stable_only=True` should give decomposition energy equal to equilibrium reaction energy!"

        # Test that we get correct behavior with a polymorph
        toy_entries = {"Li": 0.0, "Li2O": -5, "LiO2": -4, "O2": 0.0}

        toy_pd = PhaseDiagram([PDEntry(c, e) for c, e in toy_entries.items()])

        # stable entry
        assert toy_pd.get_phase_separation_energy(PDEntry("Li2O", -5)) == approx(-1.0)
        # polymorph
        assert toy_pd.get_phase_separation_energy(PDEntry("Li2O", -4)) == approx(-2.0 / 3.0)

        # Test that the method works for novel entries
        novel_stable_entry = PDEntry("Li5FeO4", -999)
        assert (
            self.pd.get_phase_separation_energy(novel_stable_entry) < 0
        ), "Novel stable entries should have negative decomposition energy!"

        novel_unstable_entry = PDEntry("Li5FeO4", 999)
        assert (
            self.pd.get_phase_separation_energy(novel_unstable_entry) > 0
        ), "Novel unstable entries should have positive decomposition energy!"

        duplicate_entry = PDEntry("Li2O", -14.31361175)
        scaled_dup_entry = PDEntry("Li4O2", -14.31361175 * 2)
        stable_entry = next(e for e in self.pd.stable_entries if e.name == "Li2O")

        assert self.pd.get_phase_separation_energy(duplicate_entry) == self.pd.get_phase_separation_energy(
            stable_entry
        ), "Novel duplicates of stable entries should have same decomposition energy!"

        assert self.pd.get_phase_separation_energy(scaled_dup_entry) == self.pd.get_phase_separation_energy(
            stable_entry
        ), "Novel scaled duplicates of stable entries should have same decomposition energy!"

    def test_get_decomposition(self):
        for entry in self.pd.stable_entries:
            assert (
                len(self.pd.get_decomposition(entry.composition)) == 1
            ), "Stable composition should have only 1 decomposition!"
        dim = len(self.pd.elements)
        for entry in self.pd.all_entries:
            n_decomp = len(self.pd.get_decomposition(entry.composition))
            assert n_decomp > 0
            assert (
                n_decomp <= dim
            ), "The number of decomposition phases can at most be equal to the number of components."

        # Just to test decomposition for a fictitious composition
        actual = {
            entry.composition.formula: amt for entry, amt in self.pd.get_decomposition(Composition("Li3Fe7O11")).items()
        }
        expected = {
            "Fe2 O2": 0.0952380952380949,
            "Li1 Fe1 O2": 0.5714285714285714,
            "Fe6 O8": 0.33333333333333393,
        }
        assert actual == approx(expected)

    def test_get_transition_chempots(self):
        for el in self.pd.elements:
            assert len(self.pd.get_transition_chempots(el)) <= len(self.pd.facets)

    def test_get_element_profile(self):
        for el in self.pd.elements:
            for entry in self.pd.stable_entries:
                if not entry.composition.is_element:
                    assert len(self.pd.get_element_profile(el, entry.composition)) <= len(self.pd.facets)

        expected = [
            {"evolution": 1.0, "chempot": -4.2582781416666666, "reaction": "Li2O + 0.5 O2 -> Li2O2"},
            {"evolution": 0, "chempot": -5.08859066, "reaction": "Li2O -> Li2O"},
            {"evolution": -1.0, "chempot": -10.48758201, "reaction": "Li2O -> 2 Li + 0.5 O2"},
        ]
        result = self.pd.get_element_profile(Element("O"), Composition("Li2O"))
        for d1, d2 in zip(expected, result):
            assert d1["evolution"] == approx(d2["evolution"])
            assert d1["chempot"] == approx(d2["chempot"])
            assert d1["reaction"] == str(d2["reaction"])

    def test_get_get_chempot_range_map(self):
        elements = [el for el in self.pd.elements if el.symbol != "Fe"]
        assert len(self.pd.get_chempot_range_map(elements)) == 10

    def test_getmu_vertices_stability_phase(self):
        results = self.pd.getmu_vertices_stability_phase(Composition("LiFeO2"), Element("O"))
        assert len(results) == approx(6)
        test_equality = False
        for c in results:
            if (
                abs(c[Element("O")] + 7.115) < 1e-2
                and abs(c[Element("Fe")] + 6.596) < 1e-2
                and abs(c[Element("Li")] + 3.931) < 1e-2
            ):
                test_equality = True
        assert test_equality, "there is an expected vertex missing in the list"

    def test_getmu_range_stability_phase(self):
        results = self.pd.get_chempot_range_stability_phase(Composition("LiFeO2"), Element("O"))
        assert results[Element("O")][1] == approx(-4.450181224)
        assert results[Element("Fe")][0] == approx(-6.5961470)
        assert results[Element("Li")][0] == approx(-3.6250022625)

    def test_get_hull_energy(self):
        for entry in self.pd.stable_entries:
            h_e = self.pd.get_hull_energy(entry.composition)
            assert h_e == approx(entry.energy)
            n_h_e = self.pd.get_hull_energy(entry.composition.fractional_composition)
            assert n_h_e == approx(entry.energy_per_atom)

    def test_get_hull_energy_per_atom(self):
        for entry in self.pd.stable_entries:
            h_e = self.pd.get_hull_energy_per_atom(entry.composition)
            assert h_e == approx(entry.energy_per_atom)

    def test_1d_pd(self):
        entry = PDEntry("H", 0)
        pd = PhaseDiagram([entry])
        decomp, e = pd.get_decomp_and_e_above_hull(PDEntry("H", 1))
        assert e == 1
        assert decomp[entry] == approx(1.0)

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
            assert crit.almost_equals(exp, rtol=0, atol=1e-5)

        comps = self.pd.get_critical_compositions(c1, c3)
        expected = [
            Composition("Fe0.4O0.6"),
            Composition("LiFeO2").fractional_composition,
            Composition("Li5FeO4").fractional_composition,
            Composition("Li2O").fractional_composition,
        ]
        for crit, exp in zip(comps, expected):
            assert crit.almost_equals(exp, rtol=0, atol=1e-5)

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
            assert crit.almost_equals(exp, rtol=0, atol=1e-5)

        comps = self.pd.get_critical_compositions(c1, c3)
        expected = [
            Composition("Fe2O3"),
            Composition("LiFeO2"),
            Composition("Li5FeO4") / 3,
            Composition("Li2O"),
        ]
        for crit, exp in zip(comps, expected):
            assert crit.almost_equals(exp, rtol=0, atol=1e-5)

        # Don't fail silently if input compositions aren't in phase diagram
        # Can be very confusing if you're working with a GrandPotentialPD
        with pytest.raises(ValueError, match="Xe1 has elements not in the phase diagram Li, Fe, O"):
            self.pd.get_critical_compositions(Composition("Xe"), Composition("Mn"))

        # For the moment, should also fail even if compositions are in the gppd
        # because it isn't handled properly
        gppd = GrandPotentialPhaseDiagram(self.pd.all_entries, {"Xe": 1}, [*self.pd.elements, Element("Xe")])
        with pytest.raises(ValueError, match="Li3 Fe1 O4 Xe1 has elements not in the phase diagram O, Fe, Li"):
            gppd.get_critical_compositions(Composition("Fe2O3"), Composition("Li3FeO4Xe"))

        # check that the function still works though
        comps = gppd.get_critical_compositions(c1, c2)
        expected = [
            Composition("Fe2O3"),
            Composition("Li0.3243244Fe0.1621621O0.51351349") * 7.4,
            Composition("Li3FeO4"),
        ]
        for crit, exp in zip(comps, expected):
            assert crit.almost_equals(exp, rtol=0, atol=1e-5)

        # case where the endpoints are identical
        assert self.pd.get_critical_compositions(c1, c1 * 2) == [c1, c1 * 2]

    def test_get_composition_chempots(self):
        c1 = Composition("Fe3.1O4")
        c2 = Composition("Fe3.2O4.1Li0.01")

        e1 = self.pd.get_hull_energy(c1)
        e2 = self.pd.get_hull_energy(c2)

        cp = self.pd.get_composition_chempots(c1)
        calc_e2 = e1 + sum(cp[k] * v for k, v in (c2 - c1).items())
        assert e2 == approx(calc_e2)

    def test_get_all_chempots(self):
        c1 = Composition("Fe3.1O4")
        c2 = Composition("FeO")

        cp1 = self.pd.get_all_chempots(c1)
        cp_result = {
            Element("Li"): -4.077061954,
            Element("Fe"): -6.741593864,
            Element("O"): -6.969907375,
        }

        for elem, energy in cp_result.items():
            assert cp1["Fe3O4-FeO-LiFeO2"][elem] == approx(energy)

        cp2 = self.pd.get_all_chempots(c2)
        cp_result = {
            Element("O"): -7.11535414,
            Element("Fe"): -6.5961471,
            Element("Li"): -3.93161518,
        }

        for elem, energy in cp_result.items():
            assert cp2["FeO-LiFeO2-Fe"][elem] == approx(energy)

    def test_get_plot(self):
        self.pd.get_plot()  # PDPlotter functionality is tested separately

    def test_as_from_dict(self):
        # test round-trip for other entry types such as ComputedEntry
        entry = ComputedEntry("H", 0.0, 0.0, entry_id="test")
        pd = PhaseDiagram([entry])
        pd_dict = pd.as_dict()
        pd_roundtrip = PhaseDiagram.from_dict(pd_dict)
        assert pd.all_entries[0].entry_id == pd_roundtrip.all_entries[0].entry_id
        dd = self.pd.as_dict()
        new_pd = PhaseDiagram.from_dict(dd)
        new_pd_dict = new_pd.as_dict()
        assert new_pd_dict == dd
        assert isinstance(pd.to_json(), str)

    def test_read_json(self):
        dumpfn(self.pd, f"{self.tmp_path}/pd.json")
        pd = loadfn(f"{self.tmp_path}/pd.json")
        assert isinstance(pd, PhaseDiagram)
        assert {*pd.as_dict()} == {*self.pd.as_dict()}

    def test_el_refs(self):
        # Create an imitation of pre_computed phase diagram with el_refs keys being
        # tuple[str, PDEntry] instead of tuple[Element, PDEntry].
        mock_el_refs = [(str(el), entry) for el, entry in self.pd.el_refs.items()]
        mock_computed_data = {**self.pd.computed_data, "el_refs": mock_el_refs}
        pd = PhaseDiagram(self.entries, computed_data=mock_computed_data)
        # Check the keys in el_refs dict have been updated to Element object via PhaseDiagram class.
        assert all(isinstance(el, Element) for el in pd.el_refs)

    def test_val_err_on_no_entries(self):
        # check that PhaseDiagram raises ValueError when building phase diagram with no entries
        for entries in [None, [], set(), ()]:
            with pytest.raises(ValueError, match="Unable to build phase diagram without entries."):
                PhaseDiagram(entries=entries)


class TestGrandPotentialPhaseDiagram(unittest.TestCase):
    def setUp(self):
        self.entries = EntrySet.from_csv(f"{module_dir}/pd_entries_test.csv")
        self.pd = GrandPotentialPhaseDiagram(self.entries, {Element("O"): -5})
        self.pd6 = GrandPotentialPhaseDiagram(self.entries, {Element("O"): -6})

    def test_stable_entries(self):
        stable_formulas = [ent.original_entry.composition.reduced_formula for ent in self.pd.stable_entries]
        expected_stable = ["Li5FeO4", "Li2FeO3", "LiFeO2", "Fe2O3", "Li2O2"]
        for formula in expected_stable:
            assert formula in stable_formulas, f"{formula} not in stable entries!"
        assert len(self.pd6.stable_entries) == 4

    def test_get_formation_energy(self):
        stable_formation_energies = {
            ent.original_entry.composition.reduced_formula: self.pd.get_form_energy(ent)
            for ent in self.pd.stable_entries
        }
        expected_formation_energies = {
            "Fe2O3": 0.0,
            "Li5FeO4": -5.30551504,
            "Li2FeO3": -2.34247415,
            "LiFeO2": -0.4302639625,
            "Li2O2": 0.0,
        }
        for formula, energy in expected_formation_energies.items():
            assert energy == approx(
                stable_formation_energies[formula]
            ), f"Calculated formation for {formula} is not correct!"

    def test_str(self):
        # using startswith since order of stable phases is random
        assert str(self.pd).startswith(
            "Fe-Li GrandPotentialPhaseDiagram with chempots = 'mu_O = -5.0000'5 stable phases: "
        )


class TestCompoundPhaseDiagram(unittest.TestCase):
    def setUp(self):
        self.entries = EntrySet.from_csv(f"{module_dir}/pd_entries_test.csv")
        self.pd = CompoundPhaseDiagram(self.entries, [Composition("Li2O"), Composition("Fe2O3")])

    def test_stable_entries(self):
        stable_formulas = [ent.name for ent in self.pd.stable_entries]
        expected_stable = ["Fe2O3", "Li5FeO4", "LiFeO2", "Li2O"]
        for formula in expected_stable:
            assert formula in stable_formulas

    def test_get_formation_energy(self):
        stable_formation_energies = {ent.name: self.pd.get_form_energy(ent) for ent in self.pd.stable_entries}
        expected_formation_energies = {
            "Li5FeO4": -7.0773284399999739,
            "Fe2O3": 0,
            "LiFeO2": -0.4745592975,
            "Li2O": 0,
        }
        for formula, energy in expected_formation_energies.items():
            assert energy == approx(stable_formation_energies[formula])

    def test_str(self):
        assert str(self.pd) == "Xf-Xg phase diagram\n4 stable phases: \nLiFeO2, Li2O, Li5FeO4, Fe2O3"


class TestPatchedPhaseDiagram(unittest.TestCase):
    def setUp(self):
        self.entries = EntrySet.from_csv(f"{module_dir}/reaction_entries_test.csv")
        # NOTE add He to test for correct behavior despite no patches involving He
        self.no_patch_entry = he_entry = PDEntry("He", -1.23)
        self.entries.add(he_entry)

        self.pd = PhaseDiagram(entries=self.entries)
        self.ppd = PatchedPhaseDiagram(entries=self.entries)

        # novel entries not in any of the patches
        self.novel_comps = [Composition("H5C2OP"), Composition("V2PH4C")]
        for c in self.novel_comps:
            assert c.chemical_system not in self.ppd.spaces

        self.novel_entries = [PDEntry(c, -39.8) for c in self.novel_comps]

    def test_get_stable_entries(self):
        assert self.pd.stable_entries == self.ppd.stable_entries

    def test_get_qhull_entries(self):
        # NOTE qhull_entry is an specially sorted list due to it's construction, we
        # can't mimic this in ppd therefore just test if sorted versions are equal.
        assert sorted(self.pd.qhull_entries, key=lambda e: e.composition) == sorted(
            self.ppd.qhull_entries, key=lambda e: e.composition
        )

    def test_get_decomposition(self):
        for comp in self.novel_comps:
            decomp_pd = self.pd.get_decomposition(comp)
            decomp_ppd = self.ppd.get_decomposition(comp)
            assert decomp_pd == approx(decomp_ppd)

    def test_get_phase_separation_energy(self):
        for entry in self.novel_entries:
            e_phase_sep_pd = self.pd.get_phase_separation_energy(entry)
            e_phase_sep_ppd = self.ppd.get_phase_separation_energy(entry)
            assert np.isclose(e_phase_sep_pd, e_phase_sep_ppd)

    def test_get_equilibrium_reaction_energy(self):
        for entry in self.pd.stable_entries:
            e_equi_rxn_pd = self.pd.get_equilibrium_reaction_energy(entry)
            e_equi_rxn_pdd = self.ppd.get_equilibrium_reaction_energy(entry)
            assert np.isclose(e_equi_rxn_pd, e_equi_rxn_pdd)

    def test_get_form_energy(self):
        for entry in self.pd.stable_entries:
            e_form_pd = self.pd.get_form_energy(entry)
            e_form_ppd = self.ppd.get_form_energy(entry)
            assert np.isclose(e_form_pd, e_form_ppd)

    def test_dimensionality(self):
        assert self.pd.dim == self.ppd.dim

        # test dims of sub PDs
        dim_counts = collections.Counter(pd.dim for pd in self.ppd.pds.values())
        assert dim_counts == {3: 7, 2: 6, 4: 2}

    def test_get_hull_energy(self):
        for comp in self.novel_comps:
            e_hull_pd = self.pd.get_hull_energy(comp)
            e_hull_ppd = self.ppd.get_hull_energy(comp)
            assert np.isclose(e_hull_pd, e_hull_ppd)

    def test_get_decomp_and_e_above_hull(self):
        for entry in self.pd.stable_entries:
            decomp_pd, e_above_hull_pd = self.pd.get_decomp_and_e_above_hull(entry)
            decomp_ppd, e_above_hull_ppd = self.ppd.get_decomp_and_e_above_hull(entry, check_stable=True)
            assert decomp_pd == decomp_ppd
            assert np.isclose(e_above_hull_pd, e_above_hull_ppd)

    def test_repr(self):
        assert repr(self.ppd) == str(self.ppd) == "PatchedPhaseDiagram covering 15 sub-spaces"

    def test_as_from_dict(self):
        ppd_dict = self.ppd.as_dict()
        assert ppd_dict["@module"] == self.ppd.__class__.__module__
        assert ppd_dict["@class"] == self.ppd.__class__.__name__
        assert ppd_dict["all_entries"] == [entry.as_dict() for entry in self.ppd.all_entries]
        assert ppd_dict["elements"] == [elem.as_dict() for elem in self.ppd.elements]
        # test round-trip dict serialization
        assert PatchedPhaseDiagram.from_dict(ppd_dict).as_dict() == ppd_dict

    def test_get_pd_for_entry(self):
        for entry in self.ppd.all_entries:
            if entry == self.no_patch_entry:
                continue
            pd = self.ppd.get_pd_for_entry(entry)
            # test that entry is in pd and pd can return valid decomp
            assert isinstance(pd.get_decomposition(entry.composition), dict)

        with pytest.raises(ValueError, match="No suitable PhaseDiagrams found for PDEntry"):
            self.ppd.get_pd_for_entry(self.no_patch_entry)

    def test_raises_on_missing_terminal_entries(self):
        entry = PDEntry("FeO", -1.23)
        with pytest.raises(ValueError, match=r"Missing terminal entries for elements \['Fe', 'O'\]"):
            PatchedPhaseDiagram(entries=[entry])

    def test_contains(self):
        for space in self.ppd.spaces:
            assert space in self.ppd
        unlikely_chem_space = frozenset(map(Element, "HBCNOFPS"))
        assert unlikely_chem_space not in self.ppd

    def test_getitem(self):
        chem_space = self.ppd.spaces[0]
        pd = self.ppd[chem_space]
        assert isinstance(pd, PhaseDiagram)
        assert chem_space in pd._qhull_spaces
        assert str(pd) == "V-C phase diagram\n4 stable phases: \nC, V, V6C5, V2C"

        with pytest.raises(KeyError, match="frozenset"):
            self.ppd[frozenset(map(Element, "HBCNOFPS"))]

    def test_iter(self):
        for pd in self.ppd:
            assert isinstance(pd, PhaseDiagram)
        assert len(self.ppd) == len(self.ppd.pds)

    def test_len(self):
        assert len(self.ppd) == len(self.ppd.pds)

    def test_setitem_and_delitem(self):
        unlikely_chem_space = frozenset(map(Element, "HBCNOFPS"))
        self.ppd[unlikely_chem_space] = self.pd
        assert unlikely_chem_space in self.ppd
        assert self.ppd[unlikely_chem_space] == self.pd
        del self.ppd[unlikely_chem_space]  # test __delitem__() and restore original state


class TestReactionDiagram(unittest.TestCase):
    def setUp(self):
        self.entries = list(EntrySet.from_csv(f"{module_dir}/reaction_entries_test.csv").entries)
        for e in self.entries:
            if e.composition.reduced_formula == "VPO5":
                entry1 = e
            elif e.composition.reduced_formula == "H4(CO)3":
                entry2 = e
        self.rd = ReactionDiagram(entry1=entry1, entry2=entry2, all_entries=self.entries[2:])

    def test_get_compound_pd(self):
        self.rd.get_compound_pd()

    def test_formula(self):
        for entry in self.rd.rxn_entries:
            assert Element.V in entry.composition
            assert Element.O in entry.composition
            assert Element.C in entry.composition
            assert Element.P in entry.composition
            assert Element.H in entry.composition
        # formed_formula = [e.composition.reduced_formula for e in self.rd.rxn_entries]
        # expected_formula = [
        #     "V0.12707182P0.12707182H0.0441989C0.03314917O0.66850829",
        #     "V0.125P0.125H0.05C0.0375O0.6625",
        #     "V0.12230216P0.12230216H0.05755396C0.04316547O0.65467626",
        #     "V0.11340206P0.11340206H0.08247423C0.06185567O0.62886598",
        #     "V0.11267606P0.11267606H0.08450704C0.06338028O0.62676056",
        #     "V0.11229947P0.11229947H0.0855615C0.06417112O0.62566845",
        #     "V0.09677419P0.09677419H0.12903226C0.09677419O0.58064516",
        #     "V0.05882353P0.05882353H0.23529412C0.17647059O0.47058824",
        #     "V0.04225352P0.04225352H0.28169014C0.21126761O0.42253521",
        # ]

        # # Please do not uncomment this test. This test is far too fragile because of numerical precision errors.
        # # Unless someone wants to make an effort to write a PROPER test which do not fail with changes in
        # # OS or numpy versions, DO NOT UNCOMMENT!
        # for formula in expected_formula:
        #     assert formula in formed_formula, f"{formed_formula=} not in {expected_formula=}"


class TestPDPlotter(unittest.TestCase):
    def setUp(self):
        entries = list(EntrySet.from_csv(f"{module_dir}/pd_entries_test.csv"))

        elemental_entries = [e for e in entries if e.elements == [Element("Li")]]
        self.pd_unary = PhaseDiagram(elemental_entries)
        self.plotter_unary_plotly = PDPlotter(self.pd_unary, backend="plotly")

        entries_LiO = [e for e in entries if "Fe" not in e.composition]
        self.pd_binary = PhaseDiagram(entries_LiO)
        self.plotter_binary_mpl = PDPlotter(self.pd_binary, backend="matplotlib")
        self.plotter_binary_plotly = PDPlotter(self.pd_binary, backend="plotly")

        self.pd_ternary = PhaseDiagram(entries)
        self.plotter_ternary_mpl = PDPlotter(self.pd_ternary, backend="matplotlib")
        self.plotter_ternary_plotly_2d = PDPlotter(self.pd_ternary, backend="plotly", ternary_style="2d")
        self.plotter_ternary_plotly_3d = PDPlotter(self.pd_ternary, backend="plotly", ternary_style="3d")

        entries.append(PDEntry("C", 0))
        self.pd_quaternary = PhaseDiagram(entries)
        self.plotter_quaternary_mpl = PDPlotter(self.pd_quaternary, backend="matplotlib")
        self.plotter_quaternary_plotly = PDPlotter(self.pd_quaternary, backend="plotly")

    def test_plot_pd_with_no_unstable(self):
        # https://github.com/materialsproject/pymatgen/issues/2885
        pd_entries = [PDEntry(comp, 0) for comp in ["Li", "Co", "O"]]
        pd = PhaseDiagram(pd_entries)
        plotter = PDPlotter(pd, backend="plotly", show_unstable=False)
        plotter.get_plot()

    def test_pd_plot_data(self):
        lines, labels, unstable_entries = self.plotter_ternary_mpl.pd_plot_data
        assert len(lines) == 22
        assert len(labels) == len(self.pd_ternary.stable_entries), "Incorrect number of lines generated!"
        assert len(unstable_entries) == len(self.pd_ternary.all_entries) - len(
            self.pd_ternary.stable_entries
        ), "Incorrect number of lines generated!"
        lines, labels, unstable_entries = self.plotter_quaternary_mpl.pd_plot_data
        assert len(lines) == 33
        assert len(labels) == len(self.pd_quaternary.stable_entries)
        assert len(unstable_entries) == len(self.pd_quaternary.all_entries) - len(self.pd_quaternary.stable_entries)
        lines, labels, unstable_entries = self.plotter_binary_mpl.pd_plot_data
        assert len(lines) == 3
        assert len(labels) == len(self.pd_binary.stable_entries)

    def test_mpl_plots(self):
        # Some very basic ("non")-tests. Just to make sure the methods are callable.
        self.plotter_binary_mpl.get_plot()
        self.plotter_ternary_mpl.get_plot()
        self.plotter_quaternary_mpl.get_plot()
        self.plotter_ternary_mpl.get_contour_pd_plot()
        self.plotter_ternary_mpl.get_chempot_range_map_plot([Element("Li"), Element("O")])
        self.plotter_ternary_mpl.plot_element_profile(Element("O"), Composition("Li2O"))

    def test_plotly_plots(self):
        # Also very basic tests. Ensures callability and 2D vs 3D properties.
        self.plotter_unary_plotly.get_plot()
        self.plotter_binary_plotly.get_plot()
        self.plotter_ternary_plotly_2d.get_plot()
        self.plotter_ternary_plotly_3d.get_plot()
        self.plotter_quaternary_plotly.get_plot()


class TestUtilityFunction(unittest.TestCase):
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
        expected = {
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
        assert uniquelines(testdata) == expected

    def test_triangular_coord(self):
        coord = [0.5, 0.5]
        coord = triangular_coord(coord)
        assert_allclose(coord, [0.75, 0.4330127])

    def test_tet_coord(self):
        coord = [0.5, 0.5, 0.5]
        coord = tet_coord(coord)
        assert_allclose(coord, [1.0, 0.57735027, 0.40824829])
