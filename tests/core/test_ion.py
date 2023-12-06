from __future__ import annotations

import random
import unittest

import pytest

from pymatgen.core import Composition, Element
from pymatgen.core.ion import Ion


class TestIon(unittest.TestCase):
    def setUp(self):
        self.comp = []
        self.comp.append(Ion.from_formula("Li+"))
        self.comp.append(Ion.from_formula("MnO4-"))
        self.comp.append(Ion.from_formula("Mn++"))
        self.comp.append(Ion.from_formula("PO3-2"))
        self.comp.append(Ion.from_formula("Fe(CN)6-3"))
        self.comp.append(Ion.from_formula("Fe(CN)6----"))
        self.comp.append(Ion.from_formula("Fe2((PO4)3(CO3)5)2-3"))
        self.comp.append(Ion.from_formula("Ca[2+]"))
        self.comp.append(Ion.from_formula("NaOH(aq)"))

    def test_init(self):
        c = Composition({"Fe": 4, "O": 16, "P": 4})
        charge = 4
        assert Ion(c, charge).formula == "Fe4 P4 O16 +4"
        f = {1: 1, 8: 1}
        charge = -1
        assert Ion(Composition(f), charge).formula == "H1 O1 -1"
        assert Ion(Composition(S=2, O=3), -2).formula == "S2 O3 -2"

    def test_charge_from_formula(self):
        assert Ion.from_formula("Li+").charge == 1
        assert Ion.from_formula("Li[+]").charge == 1
        assert Ion.from_formula("Ca[2+]").charge == 2
        assert Ion.from_formula("Ca[+2]").charge == 2
        assert Ion.from_formula("Ca++").charge == 2
        assert Ion.from_formula("Ca[++]").charge == 2
        assert Ion.from_formula("Ca2+").charge == 1

        assert Ion.from_formula("Cl-").charge == -1
        assert Ion.from_formula("Cl[-]").charge == -1
        assert Ion.from_formula("SO4[-2]").charge == -2
        assert Ion.from_formula("SO4-2").charge == -2
        assert Ion.from_formula("SO42-").charge == -1
        assert Ion.from_formula("SO4--").charge == -2
        assert Ion.from_formula("SO4[--]").charge == -2

        assert Ion.from_formula("Na[+-+]").charge == 1

    def test_special_formulas(self):
        special_formulas = [
            ("Cl-", "Cl[-1]"),
            ("H+", "H[+1]"),
            ("F-", "F[-1]"),
            ("F2", "F2(aq)"),
            ("H2", "H2(aq)"),
            ("O3", "O3(aq)"),
            ("O2", "O2(aq)"),
            ("N2", "N2(aq)"),
            ("H4O4", "H2O2(aq)"),
            ("OH-", "OH[-1]"),
            ("CH3COO-", "CH3COO[-1]"),
            ("CH3COOH", "CH3COOH(aq)"),
            ("CH3OH", "CH3OH(aq)"),
            ("H4CO", "CH3OH(aq)"),
            ("C2H6O", "C2H5OH(aq)"),
            ("C3H8O", "C3H7OH(aq)"),
            ("C4H10O", "C4H9OH(aq)"),
            ("Fe(OH)4+", "Fe(OH)4[+1]"),
            ("Zr(OH)4", "Zr(OH)4(aq)"),
        ]

        for tup in special_formulas:
            assert Ion.from_formula(tup[0]).reduced_formula == tup[1]

        assert Ion.from_formula("Fe(OH)4+").get_reduced_formula_and_factor(hydrates=True) == ("FeO2.2H2O", 1)
        assert Ion.from_formula("Zr(OH)4").get_reduced_formula_and_factor(hydrates=True) == ("ZrO2.2H2O", 1)
        assert Ion.from_formula("O").get_reduced_formula_and_factor(hydrates=False) == ("O", 1)
        assert Ion.from_formula("O2").get_reduced_formula_and_factor(hydrates=False) == ("O2", 1)
        assert Ion.from_formula("O3").get_reduced_formula_and_factor(hydrates=False) == ("O3", 1)
        assert Ion.from_formula("O6").get_reduced_formula_and_factor(hydrates=False) == ("O3", 2)
        assert Ion.from_formula("N8").get_reduced_formula_and_factor(hydrates=False) == ("N2", 4)

    def test_formula(self):
        correct_formulas = [
            "Li1 +1",
            "Mn1 O4 -1",
            "Mn1 +2",
            "P1 O3 -2",
            "Fe1 C6 N6 -3",
            "Fe1 C6 N6 -4",
            "Fe2 P6 C10 O54 -3",
            "Ca1 +2",
            "Na1 H1 O1 (aq)",
        ]
        all_formulas = [c.formula for c in self.comp]
        assert all_formulas == correct_formulas
        with pytest.raises(ValueError, match="co2 is an invalid formula"):
            Ion.from_formula("(co2)(po4)2")

    def test_mixed_valence(self):
        comp = Ion(Composition({"Fe2+": 2, "Fe3+": 4, "Li+": 8}))
        assert comp.reduced_formula == "Li4Fe3(aq)"
        assert comp.alphabetical_formula == "Fe6 Li8 (aq)"
        assert comp.formula == "Li8 Fe6 (aq)"

    def test_oxi_state_guesses(self):
        i = Ion.from_formula("SO4-2")
        assert i.oxi_state_guesses()[0].get("S") == 6
        assert i.oxi_state_guesses()[0].get("O") == -2

    def test_alphabetical_formula(self):
        correct_formulas = [
            "Li1 +1",
            "Mn1 O4 -1",
            "Mn1 +2",
            "O3 P1 -2",
            "C6 Fe1 N6 -3",
            "C6 Fe1 N6 -4",
            "C10 Fe2 O54 P6 -3",
            "Ca1 +2",
            "H1 Na1 O1 (aq)",
        ]
        all_formulas = [c.alphabetical_formula for c in self.comp]
        assert all_formulas == correct_formulas

    def test_num_atoms(self):
        correct_num_atoms = [1, 5, 1, 4, 13, 13, 72, 1, 3]
        all_natoms = [c.num_atoms for c in self.comp]
        assert all_natoms == correct_num_atoms

    def test_anonymized_formula(self):
        expected_formulas = [
            "A+1",
            "AB4-1",
            "A+2",
            "AB3-2",
            "AB6C6-3",
            "AB6C6-4",
            "AB3C5D27-3",
            "A+2",
            "ABC(aq)",
        ]
        for i, _ in enumerate(self.comp):
            assert self.comp[i].anonymized_formula == expected_formulas[i]

    def test_from_dict(self):
        sym_dict = {"P": 1, "O": 4, "charge": -2}
        assert Ion.from_dict(sym_dict).reduced_formula == "PO4[-2]", "Creation form sym_amount dictionary failed!"

    def test_as_dict(self):
        c = Ion.from_dict({"Mn": 1, "O": 4, "charge": -1})
        d = c.as_dict()
        correct_dict = {"Mn": 1.0, "O": 4.0, "charge": -1.0}
        assert d == correct_dict
        assert d["charge"] == correct_dict["charge"]
        correct_dict = {"Mn": 1.0, "O": 4.0, "charge": -1}
        d = c.to_reduced_dict
        assert d == correct_dict
        assert d["charge"] == correct_dict["charge"]

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
        assert comp1 == comp2, f"Composition equality test failed. {comp1.formula} should be equal to {comp2.formula}"
        assert hash(comp1) == hash(comp2), "Hash equality test failed!"

    def test_equality(self):
        assert self.comp[0] == (self.comp[0])
        assert self.comp[0] != self.comp[1]
        assert self.comp[0] == self.comp[0]
        assert self.comp[0] != (self.comp[1])

    def test_mul(self):
        assert (self.comp[1] * 4).formula == "Mn4 O16 -4", "Incorrect composition after addition!"

    def test_len(self):
        assert len(self.comp[1]) == 2, "Lengths are not equal!"

    def test_to_latex_string(self):
        correct_latex = [
            "Li$^{+1}$",
            "MnO$_{4}$$^{-1}$",
            "Mn$^{+2}$",
            "PO$_{3}$$^{-2}$",
            "Fe(CN)$_{6}$$^{-3}$",
            "Fe(CN)$_{6}$$^{-4}$",
            "FeP$_{3}$C$_{5}$O$_{27}$$^{-3}$",
            "Ca$^{+2}$",
            "NaOH",
        ]
        all_latex = [c.to_latex_string() for c in self.comp]
        assert all_latex == correct_latex
