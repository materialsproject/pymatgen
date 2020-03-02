# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest

from pymatgen.util.string import formula_double_format, latexify, \
    latexify_spacegroup, transformation_to_string, htmlify, unicodeify, \
    disordered_formula, unicodeify_spacegroup, unicodeify_species

from pymatgen.core import Structure


class FuncTest(unittest.TestCase):

    def test_latexify(self):
        self.assertEqual(latexify("Li3Fe2(PO4)3"),
                         "Li$_{3}$Fe$_{2}$(PO$_{4}$)$_{3}$")
        self.assertEqual(latexify("Li0.2Na0.8Cl"),
                         "Li$_{0.2}$Na$_{0.8}$Cl")

    def test_latexify_spacegroup(self):
        self.assertEqual(latexify_spacegroup("Fd-3m"), "Fd$\\overline{3}$m")
        self.assertEqual(latexify_spacegroup("P2_1/c"), "P2$_{1}$/c")

    def test_htmlify(self):
        self.assertEqual(htmlify("Li3Fe2(PO4)3"),
                         "Li<sub>3</sub>Fe<sub>2</sub>(PO<sub>4</sub>)<sub>3</sub>")
        self.assertEqual(htmlify("Li0.2Na0.8Cl"),
                         "Li<sub>0.2</sub>Na<sub>0.8</sub>Cl")

    def test_unicodeify(self):
        self.assertEqual(unicodeify("Li3Fe2(PO4)3"),
                         "Li₃Fe₂(PO₄)₃")
        self.assertRaises(ValueError, unicodeify,
                          "Li0.2Na0.8Cl")
        self.assertEqual(unicodeify_species("O2+"), "O²⁺")
        self.assertEqual(unicodeify_spacegroup("F-3m"), "F̅3m")

    def test_formula_double_format(self):
        self.assertEqual(formula_double_format(1.00), "")
        self.assertEqual(formula_double_format(2.00), "2")
        self.assertEqual(formula_double_format(2.10), "2.1")
        self.assertEqual(formula_double_format(2.10000000002), "2.1")

    def test_transformation_to_string(self):
        m = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        t = [0, 0, 0]
        s = 'x,y,z'
        ms = 'mx,my,mz'
        abc = 'a,b,c'
        self.assertEqual(s, transformation_to_string(m, t))
        self.assertEqual(ms, transformation_to_string(m, t, c='m'))
        self.assertEqual(abc, transformation_to_string(m, t, components=('a', 'b', 'c')))

        m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        t = [11, 12, 13]
        s = 'x+2y+3z+11,4x+5y+6z+12,7x+8y+9z+13'
        self.assertEqual(s, transformation_to_string(m, t))

        m = [[-1 / 2, -2 / 3, -3 / 4], [-5 / 6, -6 / 7, -7 / 8], [-8 / 9, -9 / 10, -10 / 11]]
        t = [-11 / 12, -12 / 13, -13 / 14]
        s = '-x/2-2y/3-3z/4-11/12,-5x/6-6y/7-7z/8-12/13,-8x/9-9y/10-10z/11-13/14'
        self.assertEqual(s, transformation_to_string(m, t))

    def test_disordered_formula(self):
        disordered_struct = Structure([[10, 0, 0], [0, 10, 0], [0, 0, 10]],
                                      [{'Cu': 0.25, 'Au': 0.75}],
                                      [[0, 0, 0]])

        formula_plain = disordered_formula(disordered_struct, fmt='plain')
        formula_latex = disordered_formula(disordered_struct, fmt='LaTeX')
        formula_html = disordered_formula(disordered_struct, fmt='HTML')

        self.assertEqual(formula_plain, 'CuxAu1-x x=0.25')
        self.assertEqual(formula_latex, 'Cu_{x}Au_{1-x} x=0.25')
        self.assertEqual(formula_html, 'Cu<sub>x</sub>Au<sub>1-x</sub> x=0.25')


if __name__ == "__main__":
    unittest.main()
