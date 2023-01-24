# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import unittest

import pytest

from pymatgen.core import Structure
from pymatgen.util.string import (
    Stringify,
    charge_string,
    disordered_formula,
    formula_double_format,
    htmlify,
    latexify,
    latexify_spacegroup,
    transformation_to_string,
    unicodeify,
    unicodeify_spacegroup,
    unicodeify_species,
)


class SubStr(Stringify):
    def __str__(self):
        return "Fe8O12"


class SupStr(Stringify):
    STRING_MODE = "SUPERSCRIPT"

    def to_pretty_string(self) -> str:
        return "Fe2+"

    def __str__(self):
        return "Fe**2+"


class StringifyTest(unittest.TestCase):
    def test_to_latex_string(self):
        assert SubStr().to_latex_string() == "Fe$_{8}$O$_{12}$"
        assert SupStr().to_latex_string() == "Fe$^{2+}$"

    def test_to_html_string(self):
        assert SubStr().to_html_string() == "Fe<sub>8</sub>O<sub>12</sub>"
        assert SupStr().to_html_string() == "Fe<sup>2+</sup>"

    def test_to_unicode_string(self):
        assert SubStr().to_unicode_string() == "Fe₈O₁₂"
        assert SupStr().to_unicode_string() == "Fe²⁺"


class FuncTest(unittest.TestCase):
    def test_latexify(self):
        assert latexify("Li3Fe2(PO4)3") == "Li$_{3}$Fe$_{2}$(PO$_{4}$)$_{3}$"
        assert latexify("Li0.2Na0.8Cl") == "Li$_{0.2}$Na$_{0.8}$Cl"

    def test_latexify_spacegroup(self):
        assert latexify_spacegroup("Fd-3m") == "Fd$\\overline{3}$m"
        assert latexify_spacegroup("P2_1/c") == "P2$_{1}$/c"

    def test_htmlify(self):
        assert htmlify("Li3Fe2(PO4)3") == "Li<sub>3</sub>Fe<sub>2</sub>(PO<sub>4</sub>)<sub>3</sub>"
        assert htmlify("Li0.2Na0.8Cl") == "Li<sub>0.2</sub>Na<sub>0.8</sub>Cl"

    def test_unicodeify(self):
        assert unicodeify("Li3Fe2(PO4)3") == "Li₃Fe₂(PO₄)₃"
        with pytest.raises(ValueError):
            unicodeify("Li0.2Na0.8Cl")
        assert unicodeify_species("O2+") == "O²⁺"
        assert unicodeify_spacegroup("F-3m") == "F3̅m"

    def test_formula_double_format(self):
        assert formula_double_format(1.00) == ""
        assert formula_double_format(2.00) == "2"
        assert formula_double_format(2.10) == "2.1"
        assert formula_double_format(2.10000000002) == "2.1"

    def test_charge_string(self):
        assert charge_string(1) == "[+1]"
        assert charge_string(1, brackets=False) == "+1"
        assert charge_string(1, explicit_one=False) == "[+]"

        assert charge_string(-1) == "[-1]"
        assert charge_string(-1, brackets=False) == "-1"
        assert charge_string(-1, explicit_one=False) == "[-]"

        assert charge_string(2) == "[+2]"
        assert charge_string(-4) == "[-4]"
        assert charge_string(3.5, brackets=False) == "+3.5"

        assert charge_string(0) == "(aq)"

    def test_transformation_to_string(self):
        m = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        t = [0, 0, 0]
        s = "x,y,z"
        ms = "mx,my,mz"
        abc = "a,b,c"
        assert s == transformation_to_string(m, t)
        assert ms == transformation_to_string(m, t, c="m")
        assert abc == transformation_to_string(m, t, components=("a", "b", "c"))

        m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        t = [11, 12, 13]
        s = "x+2y+3z+11,4x+5y+6z+12,7x+8y+9z+13"
        assert s == transformation_to_string(m, t)

        m = [
            [-1 / 2, -2 / 3, -3 / 4],
            [-5 / 6, -6 / 7, -7 / 8],
            [-8 / 9, -9 / 10, -10 / 11],
        ]
        t = [-11 / 12, -12 / 13, -13 / 14]
        s = "-x/2-2y/3-3z/4-11/12,-5x/6-6y/7-7z/8-12/13,-8x/9-9y/10-10z/11-13/14"
        assert s == transformation_to_string(m, t)

    def test_disordered_formula(self):
        disordered_struct = Structure(
            [[10, 0, 0], [0, 10, 0], [0, 0, 10]],
            [{"Cu": 0.25, "Au": 0.75}],
            [[0, 0, 0]],
        )

        formula_plain = disordered_formula(disordered_struct, fmt="plain")
        formula_latex = disordered_formula(disordered_struct, fmt="LaTeX")
        formula_html = disordered_formula(disordered_struct, fmt="HTML")

        assert formula_plain == "CuxAu1-x x=0.25"
        assert formula_latex == "Cu_{x}Au_{1-x} x=0.25"
        assert formula_html == "Cu<sub>x</sub>Au<sub>1-x</sub> x=0.25"


if __name__ == "__main__":
    unittest.main()
