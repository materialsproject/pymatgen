# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Created on Nov 10, 2012

@author: Shyue Ping Ong
"""

from __future__ import annotations

import random
import unittest

import pytest

from pymatgen.core.composition import ChemicalPotential, Composition
from pymatgen.core.periodic_table import Element, Species
from pymatgen.util.testing import PymatgenTest


class CompositionTest(PymatgenTest):
    def setUp(self):
        self.comp = [
            Composition("Li3Fe2(PO4)3"),
            Composition("Li3Fe(PO4)O"),
            Composition("LiMn2O4"),
            Composition("Li4O4"),
            Composition("Li3Fe2Mo3O12"),
            Composition("Li3Fe2((PO4)3(CO3)5)2"),
            Composition("Li1.5Si0.5"),
            Composition("ZnOH"),
        ]

        self.indeterminate_comp = [
            Composition.ranked_compositions_from_indeterminate_formula("Co1", True),
            Composition.ranked_compositions_from_indeterminate_formula("Co1", False),
            Composition.ranked_compositions_from_indeterminate_formula("co2o3"),
            Composition.ranked_compositions_from_indeterminate_formula("ncalu"),
            Composition.ranked_compositions_from_indeterminate_formula("calun"),
            Composition.ranked_compositions_from_indeterminate_formula("liCoo2n (pO4)2"),
            Composition.ranked_compositions_from_indeterminate_formula("(co)2 (PO)4"),
            Composition.ranked_compositions_from_indeterminate_formula("Fee3"),
        ]

    def test_immutable(self):
        try:
            self.comp[0]["Fe"] = 1
        except Exception as ex:
            assert isinstance(ex, TypeError)

        try:
            del self.comp[0]["Fe"]
        except Exception as ex:
            assert isinstance(ex, TypeError)

    def test_in(self):
        assert "Fe" in self.comp[0]
        assert "Fe" not in self.comp[2]
        assert Element("Fe") in self.comp[0]
        assert self.comp[0]["Fe"] == 2
        assert self.comp[0]["Mn"] == 0
        with pytest.raises(TypeError):
            self.comp[0]["Hello"]
        with pytest.raises(TypeError):
            self.comp[0]["Vac"]

    def test_hill_formula(self):
        c = Composition("CaCO3")
        assert c.hill_formula == "C Ca O3"
        c = Composition("C2H5OH")
        assert c.hill_formula == "C2 H6 O"

    def test_init_(self):
        with pytest.raises(ValueError):
            Composition({"H": -0.1})
        f = {"Fe": 4, "Li": 4, "O": 16, "P": 4}
        assert Composition(f).formula == "Li4 Fe4 P4 O16"
        f = {None: 4, "Li": 4, "O": 16, "P": 4}
        with pytest.raises(TypeError):
            Composition(f)
        f = {1: 2, 8: 1}
        assert Composition(f).formula == "H2 O1"
        assert Composition(Na=2, O=1).formula == "Na2 O1"

        c = Composition({"S": Composition.amount_tolerance / 2})
        assert len(c.elements) == 0

    def test_average_electroneg(self):
        val = [
            2.7224999999999997,
            2.4160000000000004,
            2.5485714285714285,
            2.21,
            2.718,
            3.08,
            1.21,
            2.43,
        ]
        for i, c in enumerate(self.comp):
            assert round(abs(c.average_electroneg - val[i]), 7) == 0

    def test_total_electrons(self):
        test_cases = {"C": 6, "SrTiO3": 84}
        for key, val in test_cases.items():
            c = Composition(key)
            assert round(abs(c.total_electrons - val), 7) == 0

    def test_formula(self):
        correct_formulas = [
            "Li3 Fe2 P3 O12",
            "Li3 Fe1 P1 O5",
            "Li1 Mn2 O4",
            "Li4 O4",
            "Li3 Fe2 Mo3 O12",
            "Li3 Fe2 P6 C10 O54",
            "Li1.5 Si0.5",
            "Zn1 H1 O1",
        ]
        all_formulas = [c.formula for c in self.comp]
        assert all_formulas == correct_formulas
        with pytest.raises(ValueError):
            Composition("(co2)(po4)2")

        assert Composition("K Na 2").reduced_formula == "KNa2"

        assert Composition("K3 Na 2").reduced_formula == "K3Na2"

        assert Composition("Na 3 Zr (PO 4) 3").reduced_formula == "Na3Zr(PO4)3"

    def test_to_latex_html_unicode(self):
        assert self.comp[0].to_latex_string() == "Li$_{3}$Fe$_{2}$P$_{3}$O$_{12}$"
        assert self.comp[0].to_html_string() == "Li<sub>3</sub>Fe<sub>2</sub>P<sub>3</sub>O<sub>12</sub>"
        assert self.comp[0].to_unicode_string() == "Li₃Fe₂P₃O₁₂"

    def test_iupac_formula(self):
        correct_formulas = [
            "Li3 Fe2 P3 O12",
            "Li3 Fe1 P1 O5",
            "Li1 Mn2 O4",
            "Li4 O4",
            "Li3 Mo3 Fe2 O12",
            "Li3 Fe2 C10 P6 O54",
            "Li1.5 Si0.5",
            "Zn1 H1 O1",
        ]
        all_formulas = [c.iupac_formula for c in self.comp]
        assert all_formulas == correct_formulas

    def test_mixed_valence(self):
        comp = Composition({"Fe2+": 2, "Fe3+": 4, "Li+": 8})
        assert comp.reduced_formula == "Li4Fe3"
        assert comp.alphabetical_formula == "Fe6 Li8"
        assert comp.formula == "Li8 Fe6"

    def test_indeterminate_formula(self):
        correct_formulas = [
            ["Co1"],
            ["Co1", "C1 O1"],
            ["Co2 O3", "C1 O5"],
            ["N1 Ca1 Lu1", "U1 Al1 C1 N1"],
            ["N1 Ca1 Lu1", "U1 Al1 C1 N1"],
            [
                "Li1 Co1 P2 N1 O10",
                "Li1 Co1 Po8 N1 O2",
                "Li1 P2 C1 N1 O11",
                "Li1 Po8 C1 N1 O3",
            ],
            ["Co2 P4 O4", "Co2 Po4", "P4 C2 O6", "Po4 C2 O2"],
            [],
        ]
        for i, c in enumerate(correct_formulas):
            assert [Composition(comp) for comp in c] == self.indeterminate_comp[i]

    def test_alphabetical_formula(self):
        correct_formulas = [
            "Fe2 Li3 O12 P3",
            "Fe1 Li3 O5 P1",
            "Li1 Mn2 O4",
            "Li4 O4",
            "Fe2 Li3 Mo3 O12",
            "C10 Fe2 Li3 O54 P6",
            "Li1.5 Si0.5",
            "H1 O1 Zn1",
        ]
        all_formulas = [c.alphabetical_formula for c in self.comp]
        assert all_formulas == correct_formulas

    def test_reduced_composition(self):
        correct_reduced_formulas = [
            "Li3Fe2(PO4)3",
            "Li3FePO5",
            "LiMn2O4",
            "Li2O2",
            "Li3Fe2(MoO4)3",
            "Li3Fe2P6(C5O27)2",
            "Li1.5Si0.5",
            "ZnHO",
        ]
        for idx, comp in enumerate(self.comp):
            assert comp.reduced_composition == Composition(correct_reduced_formulas[idx])

    def test_reduced_formula(self):
        correct_reduced_formulas = [
            "Li3Fe2(PO4)3",
            "Li3FePO5",
            "LiMn2O4",
            "Li2O2",
            "Li3Fe2(MoO4)3",
            "Li3Fe2P6(C5O27)2",
            "Li1.5Si0.5",
            "ZnHO",
        ]
        all_formulas = [c.reduced_formula for c in self.comp]
        assert all_formulas == correct_reduced_formulas

        # test iupac reduced formula (polyanions should still appear at the end)
        all_formulas = [c.get_reduced_formula_and_factor(iupac_ordering=True)[0] for c in self.comp]
        assert all_formulas == correct_reduced_formulas
        assert Composition("H6CN").get_integer_formula_and_factor(iupac_ordering=True)[0] == "CNH6"

        # test rounding
        c = Composition({"Na": 2 - Composition.amount_tolerance / 2, "Cl": 2})
        assert c.reduced_formula == "NaCl"

    def test_integer_formula(self):
        correct_reduced_formulas = [
            "Li3Fe2(PO4)3",
            "Li3FePO5",
            "LiMn2O4",
            "Li2O2",
            "Li3Fe2(MoO4)3",
            "Li3Fe2P6(C5O27)2",
            "Li3Si",
            "ZnHO",
        ]
        all_formulas = [c.get_integer_formula_and_factor()[0] for c in self.comp]
        assert all_formulas == correct_reduced_formulas
        assert Composition("Li0.5O0.25").get_integer_formula_and_factor() == ("Li2O", 0.25)
        assert Composition("O0.25").get_integer_formula_and_factor() == ("O2", 0.125)
        formula, factor = Composition("Li0.16666667B1.0H1.0").get_integer_formula_and_factor()
        assert formula == "Li(BH)6"
        assert round(abs(factor - 1 / 6), 7) == 0

        # test iupac reduced formula (polyanions should still appear at the end)
        all_formulas = [c.get_integer_formula_and_factor(iupac_ordering=True)[0] for c in self.comp]
        assert all_formulas == correct_reduced_formulas
        assert Composition("H6CN0.5").get_integer_formula_and_factor(iupac_ordering=True) == ("C2NH12", 0.5)

    def test_num_atoms(self):
        correct_num_atoms = [20, 10, 7, 8, 20, 75, 2, 3]

        all_natoms = [c.num_atoms for c in self.comp]
        assert all_natoms == correct_num_atoms

    def test_weight(self):
        correct_weights = [
            417.427086,
            187.63876199999999,
            180.81469,
            91.7616,
            612.3258,
            1302.430172,
            24.454250000000002,
            82.41634,
        ]
        all_weights = [c.weight for c in self.comp]
        self.assertArrayAlmostEqual(all_weights, correct_weights, 5)

    def test_get_atomic_fraction(self):
        correct_at_frac = {"Li": 0.15, "Fe": 0.1, "P": 0.15, "O": 0.6}
        for el in ["Li", "Fe", "P", "O"]:
            assert self.comp[0].get_atomic_fraction(el) == correct_at_frac[el], "Wrong computed atomic fractions"
        assert self.comp[0].get_atomic_fraction("S") == 0, "Wrong computed atomic fractions"

    def test_anonymized_formula(self):
        expected_formulas = [
            "A2B3C3D12",
            "ABC3D5",
            "AB2C4",
            "AB",
            "A2B3C3D12",
            "A2B3C6D10E54",
            "A0.5B1.5",
            "ABC",
        ]
        for idx, comp in enumerate(self.comp):
            assert comp.anonymized_formula == expected_formulas[idx]

    def test_get_wt_fraction(self):
        correct_wt_frac = {
            "Li": 0.0498841610868,
            "Fe": 0.267567687258,
            "P": 0.222604831158,
            "O": 0.459943320496,
        }
        for el in ["Li", "Fe", "P", "O"]:
            assert (
                round(abs(correct_wt_frac[el] - self.comp[0].get_wt_fraction(el)), 5) == 0
            ), "Wrong computed weight fraction"
        assert self.comp[0].get_wt_fraction(Element("S")) == 0, "Wrong computed weight fractions"

    def test_from_dict(self):
        sym_dict = {"Fe": 6, "O": 8}
        assert Composition.from_dict(sym_dict).reduced_formula == "Fe3O4", "Creation form sym_amount dictionary failed!"
        comp = Composition({"Fe2+": 2, "Fe3+": 4, "O2-": 8})
        comp2 = Composition.from_dict(comp.as_dict())
        assert comp == comp2

    def test_from_weight_dict(self):
        weight_dict_list = [{"Ti": 90, "V": 6, "Al": 4}, {"Ni": 60, "Ti": 40}, {"H": 0.1119, "O": 0.8881}]
        formula_list = ["Ti87.6 V5.5 Al6.9", "Ti44.98 Ni55.02", "H2O"]

        for weight_dict, formula in zip(weight_dict_list, formula_list):
            c1 = Composition(formula).fractional_composition
            c2 = Composition.from_weight_dict(weight_dict).fractional_composition
            assert set(c1.elements) == set(c2.elements)
            for el in c1.elements:
                assert c1[el] == pytest.approx(c2[el], abs=1e-3)

    def test_tofrom_weight_dict(self):
        for c in self.comp:
            c2 = Composition().from_weight_dict(c.to_weight_dict)
            c.almost_equals(c2)

    def test_as_dict(self):
        c = Composition.from_dict({"Fe": 4, "O": 6})
        d = c.as_dict()
        correct_dict = {"Fe": 4.0, "O": 6.0}
        assert d["Fe"] == correct_dict["Fe"]
        assert d["O"] == correct_dict["O"]
        correct_dict = {"Fe": 2.0, "O": 3.0}
        d = c.to_reduced_dict
        assert isinstance(d, dict)
        assert d["Fe"] == correct_dict["Fe"]
        assert d["O"] == correct_dict["O"]

    def test_pickle(self):
        for c in self.comp:
            self.serialize_with_pickle(c, test_eq=True)
            self.serialize_with_pickle(c.to_data_dict, test_eq=True)

    def test_to_data_dict(self):
        comp = Composition("Fe0.00009Ni0.99991")
        d = comp.to_data_dict
        assert round(abs(d["reduced_cell_composition"]["Fe"] - 9e-5), 7) == 0

    def test_add(self):
        assert (self.comp[0] + self.comp[2]).formula == "Li4 Mn2 Fe2 P3 O16", "Incorrect composition after addition!"
        assert (self.comp[3] + {"Fe": 4, "O": 4}).formula == "Li4 Fe4 O8", "Incorrect composition after addition!"

        Fe = Element("Fe")
        assert self.comp[0].__add__(Fe) == NotImplemented  # pylint: disable=C2801

    def test_sub(self):
        assert (self.comp[0] - Composition("Li2O")).formula == "Li1 Fe2 P3 O11", "Incorrect composition after addition!"
        assert (self.comp[0] - {"Fe": 2, "O": 3}).formula == "Li3 P3 O9"

        with pytest.raises(ValueError):
            Composition("O") - Composition("H")

        # check that S is completely removed by subtraction
        c1 = Composition({"S": 1 + Composition.amount_tolerance / 2, "O": 1})
        c2 = Composition({"S": 1})
        assert len((c1 - c2).elements) == 1

        Fe = Element("Fe")
        assert self.comp[0].__add__(Fe) == NotImplemented  # pylint: disable=C2801

    def test_mul(self):
        assert (self.comp[0] * 4).formula == "Li12 Fe8 P12 O48"
        assert (3 * self.comp[1]).formula == "Li9 Fe3 P3 O15"

    def test_div(self):
        assert (self.comp[0] / 4).formula == "Li0.75 Fe0.5 P0.75 O3"

    def test_equals(self):
        # generate randomized compositions for robustness (tests might pass for specific elements
        # but fail for others)
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
        assert comp1 == comp2, f"Composition equality test failed. {comp1.formula} should be equal to {comp2.formula}"
        assert hash(comp1) == hash(comp2), "Hashcode equality test failed!"

        c1, c2 = self.comp[:2]
        assert c1 == c1
        assert c1 != c2

    def test_hash_robustness(self):
        c1 = Composition(f"O{0.2}Fe{0.8}Na{Composition.amount_tolerance*0.99}")
        c2 = Composition(f"O{0.2}Fe{0.8}Na{Composition.amount_tolerance*1.01}")
        c3 = Composition(f"O{0.2}Fe{0.8+Composition.amount_tolerance*0.99}")

        assert c1 == c3, "__eq__ not robust"
        assert (c1 == c3) == (hash(c1) == hash(c3)), "Hash doesn't match eq when true"
        assert hash(c1) != hash(c2), "Hash equal for different chemical systems"

    def test_comparisons(self):
        c1 = Composition({"S": 1})
        c1_1 = Composition({"S": 1.00000000000001})
        c2 = Composition({"S": 2})
        c3 = Composition({"O": 1})
        c4 = Composition({"O": 1, "S": 1})
        assert not c1 > c2
        assert not c1_1 > c1
        assert not c1_1 < c1
        assert c1 > c3
        assert c3 < c1
        assert c4 > c1
        assert sorted([c1, c1_1, c2, c4, c3]) == [c3, c1, c1_1, c4, c2]

        Fe = Element("Fe")
        assert c1 != Fe, NotImplemented
        assert c1 != Fe
        with pytest.raises(TypeError):
            c1 < Fe  # noqa: B015

    def test_almost_equals(self):
        c1 = Composition({"Fe": 2.0, "O": 3.0, "Mn": 0})
        c2 = Composition({"O": 3.2, "Fe": 1.9, "Zn": 0})
        c3 = Composition({"Ag": 2.0, "O": 3.0})
        c4 = Composition({"Fe": 2.0, "O": 3.0, "Ag": 2.0})
        assert c1.almost_equals(c2, rtol=0.1)
        assert not c1.almost_equals(c2, rtol=0.01)
        assert not c1.almost_equals(c3, rtol=0.1)
        assert not c1.almost_equals(c4, rtol=0.1)

    def test_equality(self):
        assert self.comp[0] == self.comp[0]
        assert self.comp[0] != self.comp[1]
        assert self.comp[0] == self.comp[0]
        assert self.comp[0] != self.comp[1]

    def test_fractional_composition(self):
        for c in self.comp:
            assert round(abs(c.fractional_composition.num_atoms - 1), 7) == 0

    def test_init_numerical_tolerance(self):
        assert Composition({"B": 1, "C": -1e-12}) == Composition("B")

    def test_negative_compositions(self):
        assert Composition("Li-1(PO-1)4", allow_negative=True).formula == "Li-1 P4 O-4"
        assert Composition("Li-1(PO-1)4", allow_negative=True).reduced_formula == "Li-1(PO-1)4"
        assert Composition("Li-2Mg4", allow_negative=True).reduced_composition == Composition(
            "Li-1Mg2", allow_negative=True
        )
        assert Composition("Li-2.5Mg4", allow_negative=True).reduced_composition == Composition(
            "Li-2.5Mg4", allow_negative=True
        )

        # test math
        c1 = Composition("LiCl", allow_negative=True)
        c2 = Composition("Li")
        assert c1 - 2 * c2 == Composition({"Li": -1, "Cl": 1}, allow_negative=True)
        assert (c1 + c2).allow_negative is True
        assert c1 / -1 == Composition("Li-1Cl-1", allow_negative=True)

        # test num_atoms
        c1 = Composition("Mg-1Li", allow_negative=True)
        assert c1.num_atoms == 2
        assert c1.get_atomic_fraction("Mg") == 0.5
        assert c1.get_atomic_fraction("Li") == 0.5
        assert c1.fractional_composition == Composition("Mg-0.5Li0.5", allow_negative=True)

        # test copy
        assert c1.copy() == c1

        # test species
        c1 = Composition({"Mg": 1, "Mg2+": -1}, allow_negative=True)
        assert c1.num_atoms == 2
        assert c1.element_composition == Composition()
        assert c1.average_electroneg == 1.31

    def test_special_formulas(self):
        special_formulas = {
            "LiO": "Li2O2",
            "NaO": "Na2O2",
            "KO": "K2O2",
            "HO": "H2O2",
            "CsO": "Cs2O2",
            "RbO": "Rb2O2",
            "O": "O2",
            "N": "N2",
            "F": "F2",
            "Cl": "Cl2",
            "H": "H2",
        }
        for k, v in special_formulas.items():
            assert Composition(k).reduced_formula == v

    def test_oxi_state_guesses(self):
        assert Composition("LiFeO2").oxi_state_guesses() == ({"Li": 1, "Fe": 3, "O": -2},)

        assert Composition("Fe4O5").oxi_state_guesses() == ({"Fe": 2.5, "O": -2},)

        assert Composition("V2O3").oxi_state_guesses() == ({"V": 3, "O": -2},)

        # all_oxidation_states produces *many* possible responses
        assert len(Composition("MnO").oxi_state_guesses(all_oxi_states=True)) == 4

        # can't balance b/c missing V4+
        assert Composition("VO2").oxi_state_guesses(oxi_states_override={"V": [2, 3, 5]}) == []

        # missing V4+, but can balance due to additional sites
        assert Composition("V2O4").oxi_state_guesses(oxi_states_override={"V": [2, 3, 5]}) == ({"V": 4, "O": -2},)

        # multiple solutions - Mn/Fe = 2+/4+ or 3+/3+ or 4+/2+
        assert len(Composition("MnFeO3").oxi_state_guesses(oxi_states_override={"Mn": [2, 3, 4], "Fe": [2, 3, 4]})) == 3

        # multiple solutions prefers 3/3 over 2/4 or 4/2
        assert Composition("MnFeO3").oxi_state_guesses(oxi_states_override={"Mn": [2, 3, 4], "Fe": [2, 3, 4]})[0] == {
            "Mn": 3,
            "Fe": 3,
            "O": -2,
        }

        # target charge of 1
        assert Composition("V2O6").oxi_state_guesses(oxi_states_override={"V": [2, 3, 4, 5]}, target_charge=-2) == (
            {"V": 5, "O": -2},
        )

        # max_sites for very large composition - should timeout if incorrect
        assert Composition("Li10000Fe10000P10000O40000").oxi_state_guesses(max_sites=7)[0] == {
            "Li": 1,
            "Fe": 2,
            "P": 5,
            "O": -2,
        }

        # max_sites for very large composition - should timeout if incorrect
        assert Composition("Li10000Fe10000P10000O40000").oxi_state_guesses(max_sites=-1)[0] == {
            "Li": 1,
            "Fe": 2,
            "P": 5,
            "O": -2,
        }

        # negative max_sites less than -1 - should throw error if cannot reduce
        # to under the abs(max_sites) number of sites. Will also timeout if
        # incorrect.
        assert Composition("Sb10000O10000F10000").oxi_state_guesses(max_sites=-3)[0] == {"Sb": 3, "O": -2, "F": -1}
        with pytest.raises(ValueError):
            Composition("LiOF").oxi_state_guesses(max_sites=-2)

        with pytest.raises(ValueError):
            Composition("V2O3").oxi_state_guesses(max_sites=1)

    def test_oxi_state_decoration(self):
        # Basic test: Get compositions where each element is in a single charge state
        decorated = Composition("H2O").add_charges_from_oxi_state_guesses()
        assert Species("H", 1) in decorated
        assert decorated.get(Species("H", 1)) == 2

        # Test: More than one charge state per element
        decorated = Composition("Fe3O4").add_charges_from_oxi_state_guesses()
        assert decorated.get(Species("Fe", 2)) == 1
        assert decorated.get(Species("Fe", 3)) == 2
        assert decorated.get(Species("O", -2)) == 4

        # Test: No possible charge states
        #   It should return an uncharged composition
        decorated = Composition("NiAl").add_charges_from_oxi_state_guesses()
        assert decorated.get(Species("Ni", 0)) == 1
        assert decorated.get(Species("Al", 0)) == 1

    def test_Metallofullerene(self):
        # Test: Parse Metallofullerene formula (e.g. Y3N@C80)
        formula = "Y3N@C80"
        sym_dict = {"Y": 3, "N": 1, "C": 80}
        cmp = Composition(formula)
        cmp2 = Composition.from_dict(sym_dict)
        assert cmp == cmp2

    def test_contains_element_type(self):
        formula = "EuTiO3"
        cmp = Composition(formula)
        assert cmp.contains_element_type("lanthanoid")
        assert not cmp.contains_element_type("noble_gas")
        assert cmp.contains_element_type("f-block")
        assert not cmp.contains_element_type("s-block")

    def test_chemical_system(self):
        assert Composition({"Na": 1, "Cl": 1}).chemical_system == "Cl-Na"
        assert Composition({"Na+": 1, "Cl-": 1}).chemical_system == "Cl-Na"

    def test_is_valid(self):
        formula = "NaCl"
        cmp = Composition(formula)
        assert cmp.valid

        formula = "NaClX"
        cmp = Composition(formula)
        assert not cmp.valid

        with pytest.raises(ValueError):
            Composition("NaClX", strict=True)

    def test_remove_charges(self):
        cmp1 = Composition({"Al3+": 2.0, "O2-": 3.0})

        cmp2 = Composition({"Al": 2.0, "O": 3.0})
        assert str(cmp1) != str(cmp2)

        cmp1 = cmp1.remove_charges()
        assert str(cmp1) == str(cmp2)

        cmp1 = cmp1.remove_charges()
        assert str(cmp1) == str(cmp2)

        cmp1 = Composition({"Fe3+": 2.0, "Fe2+": 3.0, "O2-": 6.0})
        cmp2 = Composition({"Fe": 5.0, "O": 6.0})
        assert str(cmp1) != str(cmp2)

        cmp1 = cmp1.remove_charges()
        assert str(cmp1) == str(cmp2)

    def test_replace(self):
        Fe2O3 = Composition("Fe2O3")
        Cu2O3 = Composition("Cu2O3")
        MgCuO3 = Composition("MgCuO3")
        Mg2Cu2O3 = Composition("Mg2Cu2O3")

        Cu2O3_repl = Fe2O3.replace({"Fe": "Cu"})
        assert Cu2O3_repl == Cu2O3

        # handles one-to-many substitutions
        MgCuO3_repl = Fe2O3.replace({"Fe": {"Cu": 0.5, "Mg": 0.5}})
        assert MgCuO3_repl == MgCuO3

        # handles unnormalized one-to-many substitutions
        Mg2Cu2O3_repl = Fe2O3.replace({"Fe": {"Cu": 1, "Mg": 1}})
        assert Mg2Cu2O3_repl == Mg2Cu2O3

        # leaves the composition unchanged when replacing non-existent species
        assert Fe2O3 == Fe2O3.replace({"Li": "Cu"})

        # check for complex substitutions where element is involved at
        # multiple places
        Ca2NF = Composition("Ca2NF")
        example_sub_1 = {"Ca": "Sr", "N": "O", "F": "O"}
        c_new_1 = Ca2NF.replace(example_sub_1)
        assert c_new_1 == Composition("Sr2O2")

        example_sub_2 = {"Ca": "Sr", "N": "F", "F": "Cl"}
        c_new_2 = Ca2NF.replace(example_sub_2)
        assert c_new_2 == Composition("Sr2ClF")

        example_sub_3 = {"Ca": "Sr", "N": "F", "F": "N"}
        c_new_3 = Ca2NF.replace(example_sub_3)
        assert c_new_3 == Composition("Sr2NF")

        # Check with oxidation-state decorated compositions
        Ca2NF_oxi = Ca2NF.add_charges_from_oxi_state_guesses()
        example_sub_4 = {"Ca2+": "Mg2+", "N3-": "O2-", "F-": "O2-"}
        c_new_4 = Ca2NF_oxi.replace(example_sub_4)
        assert c_new_4 == Composition("Mg2O2").add_charges_from_oxi_state_guesses()


class ChemicalPotentialTest(unittest.TestCase):
    def test_init(self):
        d = {"Fe": 1, Element("Fe"): 1}
        with pytest.raises(ValueError):
            ChemicalPotential(d)
        for k in ChemicalPotential(Fe=1):
            assert isinstance(k, Element)

    def test_math(self):
        fepot = ChemicalPotential({"Fe": 1})
        opot = ChemicalPotential({"O": 2.1})
        pots = ChemicalPotential({"Fe": 1, "O": 2.1})
        potsx2 = ChemicalPotential({"Fe": 2, "O": 4.2})
        feo2 = Composition("FeO2")

        # test get_energy()
        assert round(abs(pots.get_energy(feo2) - 5.2), 7) == 0
        assert round(abs(fepot.get_energy(feo2, False) - 1), 7) == 0
        with pytest.raises(ValueError):
            fepot.get_energy(feo2)

        # test multiplication
        with pytest.raises(TypeError):
            pots * pots
        assert pots * 2 == potsx2
        assert 2 * pots == potsx2

        # test division
        assert potsx2 / 2 == pots
        assert pots.__div__(pots) == NotImplemented
        assert pots.__div__(feo2) == NotImplemented

        # test add/subtract
        assert pots + pots == potsx2
        assert potsx2 - pots == pots
        assert fepot + opot == pots
        assert fepot - opot == pots - opot - opot


if __name__ == "__main__":
    unittest.main()
