from __future__ import annotations

import numpy as np
import pytest

from pymatgen.core import Composition, Element
from pymatgen.core.ion import Ion


class TestIon:
    def setup_method(self):
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
        comp = Composition({"Fe": 4, "O": 16, "P": 4})
        charge = 4
        assert Ion(comp, charge).formula == "Fe4 P4 O16 +4"
        formula_dict = {1: 1, 8: 1}
        charge = -1
        assert Ion(Composition(formula_dict), charge).formula == "H1 O1 -1"
        assert Ion(Composition(S=2, O=3), -2).formula == "S2 O3 -2"

    def test_charge_from_formula(self):
        assert Ion.from_formula("Li+").charge == 1
        assert Ion.from_formula("Li[+]").charge == 1
        assert Ion.from_formula("Ca[2+]").charge == 2
        assert Ion.from_formula("Ca[+2]").charge == 2
        assert Ion.from_formula("Ca++").charge == 2
        assert Ion.from_formula("Ca[++]").charge == 2
        assert Ion.from_formula("Ca2+").charge == 1
        assert Ion.from_formula("C2O4-2").charge == -2
        assert Ion.from_formula("CO2").charge == 0

        assert Ion.from_formula("Cl-").charge == -1
        assert Ion.from_formula("Cl[-]").charge == -1
        assert Ion.from_formula("SO4[-2]").charge == -2
        assert Ion.from_formula("SO4-2").charge == -2
        assert Ion.from_formula("SO42-").charge == -1
        assert Ion.from_formula("SO4--").charge == -2
        assert Ion.from_formula("SO4[--]").charge == -2
        assert Ion.from_formula("N3-").charge == -1

        assert Ion.from_formula("Na[+-+]").charge == 1

    def test_get_reduced_formula_and_factor_special_formulas(self):
        special_formulas = [
            ("Cl-", "Cl[-1]"),
            ("H+", "H[+1]"),
            ("F-", "F[-1]"),
            ("I-", "I[-1]"),
            ("F2", "F2(aq)"),
            ("H2", "H2(aq)"),
            ("O3", "O3(aq)"),
            ("O2", "O2(aq)"),
            ("N2", "N2(aq)"),
            ("NaOH", "NaOH(aq)"),
            ("H4O4", "H2O2(aq)"),
            ("OH-", "OH[-1]"),
            ("H2PO4-", "H2PO4[-1]"),
            ("CH3COO-", "CH3COO[-1]"),
            ("CH3COOH", "CH3COOH(aq)"),
            ("CH3OH", "CH3OH(aq)"),
            ("H4CO", "CH3OH(aq)"),
            ("C2O4--", "C2O4[-2]"),
            ("CO2", "CO2(aq)"),
            ("CO3--", "CO3[-2]"),
            ("CH4", "CH4(aq)"),
            ("NH4+", "NH4[+1]"),
            ("NH3", "NH3(aq)"),
            ("N3-", "N3[-1]"),
            ("HCOO-", "HCO2[-1]"),
            ("C2H6O", "C2H5OH(aq)"),
            ("C3H8O", "C3H7OH(aq)"),
            ("C4H10O", "C4H9OH(aq)"),
            ("Fe(OH)4+", "Fe(OH)4[+1]"),
            ("Zr(OH)4", "Zr(OH)4(aq)"),
        ]
        for tup in special_formulas:
            assert Ion.from_formula(tup[0]).reduced_formula == tup[1], (
                f"Expected {tup[1]} but got {Ion.from_formula(tup[0]).reduced_formula}"
            )

        assert Ion.from_formula("Fe(OH)4+").get_reduced_formula_and_factor(hydrates=True) == ("FeO2.2H2O", 1)
        assert Ion.from_formula("Zr(OH)4").get_reduced_formula_and_factor(hydrates=True) == ("ZrO2.2H2O", 1)
        assert Ion.from_formula("O").get_reduced_formula_and_factor(hydrates=False) == (
            "O",
            1,
        )
        assert Ion.from_formula("O2").get_reduced_formula_and_factor(hydrates=False) == ("O2", 1)
        assert Ion.from_formula("O3").get_reduced_formula_and_factor(hydrates=False) == ("O3", 1)
        assert Ion.from_formula("O6").get_reduced_formula_and_factor(hydrates=False) == ("O3", 2)
        assert Ion.from_formula("N8").get_reduced_formula_and_factor(hydrates=False) == ("N2", 4)

    def test_get_reduced_formula_and_factor_non_int_amounts(self):
        """Ensure get_reduced_formula_and_factor returns early for non-integer amounts."""
        ion = Ion({"H": 0.5, "O": 1}, charge=-1)
        formula, factor = ion.get_reduced_formula_and_factor()

        assert formula == "H0.5O1-1"
        assert factor == 1

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

        with pytest.raises(ValueError, match="Invalid formula"):
            Ion.from_formula("Mn[2+3]")

        with pytest.raises(ValueError, match="co2 is an invalid formula"):
            Ion.from_formula("(co2)(po4)2")  # raised by Composition._parse_formula

    def test_mixed_valence(self):
        comp = Ion(Composition({"Fe2+": 2, "Fe3+": 4, "Li+": 8}))
        assert comp.reduced_formula == "Li4Fe3(aq)"
        assert comp.alphabetical_formula == "Fe6 Li8 (aq)"
        assert comp.formula == "Li8 Fe6 (aq)"

    def test_oxi_state_guesses(self):
        ion = Ion.from_formula("SO4-2")
        assert ion.oxi_state_guesses()[0].get("S") == 6
        assert ion.oxi_state_guesses()[0].get("O") == -2

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
        all_n_atoms = [c.num_atoms for c in self.comp]
        assert all_n_atoms == correct_num_atoms

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
        for idx, expected in enumerate(expected_formulas):
            assert self.comp[idx].anonymized_formula == expected

    def test_from_dict(self):
        sym_dict = {"P": 1, "O": 4, "charge": -2}
        ion = Ion.from_dict(sym_dict)
        assert ion.reduced_formula == "PO4[-2]"

        # Ensure input dict was not modified
        assert sym_dict == {"P": 1, "O": 4, "charge": -2}

    def test_as_dict(self):
        ion = Ion.from_dict({"Mn": 1, "O": 4, "charge": -1})
        assert ion.as_dict() == {"Mn": 1.0, "O": 4.0, "charge": -1.0}

    def test_as_reduced_dict(self):
        ion = Ion.from_dict({"Mn": 2, "O": 8, "charge": -2})
        assert ion.as_reduced_dict() == {"Mn": 1, "O": 4, "charge": -1}

        ion = Ion.from_dict({"Mn": 1, "O": 4, "charge": -1})
        assert ion.as_reduced_dict() == {"Mn": 1, "O": 4, "charge": -1}

    def test_equals(self):
        rng = np.random.default_rng()
        random_z = rng.integers(1, 93)
        fixed_el = Element.from_Z(random_z)
        other_z = rng.integers(1, 93)
        while other_z == random_z:
            other_z = rng.integers(1, 93)
        comp1 = Ion(Composition({fixed_el: 1, Element.from_Z(other_z): 0}), 1)
        other_z = rng.integers(1, 93)
        while other_z == random_z:
            other_z = rng.integers(1, 93)
        comp2 = Ion(Composition({fixed_el: 1, Element.from_Z(other_z): 0}), 1)
        assert comp1 == comp2, f"Composition equality test failed. {comp1.formula} should be equal to {comp2.formula}"
        assert hash(comp1) == hash(comp2), "Hash equality test failed!"

    def test_equality(self):
        assert self.comp[0] == (self.comp[0])
        assert self.comp[0] != self.comp[1]
        assert self.comp[0] == self.comp[0]
        assert self.comp[0] != (self.comp[1])

        assert Ion.from_dict({"Na": 1, "charge": 1}).__eq__("Na+") is NotImplemented

    def test_add(self):
        ion1 = Ion.from_dict({"Na": 1, "charge": 1})
        ion2 = Ion.from_dict({"Cl": 1, "charge": -1})
        result = ion1 + ion2

        assert isinstance(result, Ion)
        assert result.composition.as_dict() == {"Na": 1.0, "Cl": 1.0}
        assert result.charge == 0

    def test_sub(self):
        ion1 = Ion.from_dict({"Ca": 2, "charge": 2})
        ion2 = Ion.from_dict({"Ca": 1, "charge": 1})
        result = ion1 - ion2

        assert isinstance(result, Ion)
        assert result.composition.as_dict() == {"Ca": 1.0}
        assert result.charge == 1

    def test_mul(self):
        assert (self.comp[1] * 4).formula == "Mn4 O16 -4", "Incorrect composition after addition!"

    def test_str(self):
        assert str(Ion.from_dict({"Mn": 1, "O": 4, "charge": -1})) == "Mn1 O4 -1"

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
