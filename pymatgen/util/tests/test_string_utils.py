# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Aug 26, 2012"

import unittest

from pymatgen.util.string import formula_double_format, latexify, \
    latexify_spacegroup, transformation_to_string


class FuncTest(unittest.TestCase):

    def test_latexify(self):
        self.assertEqual(latexify("Li3Fe2(PO4)3"),
                         "Li$_{3}$Fe$_{2}$(PO$_{4}$)$_{3}$")
        self.assertEqual(latexify("Li0.2Na0.8Cl"),
                         "Li$_{0.2}$Na$_{0.8}$Cl")

    def test_latexify_spacegroup(self):
        self.assertEqual(latexify_spacegroup("Fd-3m"), "Fd$\\overline{3}$m")
        self.assertEqual(latexify_spacegroup("P2_1/c"), "P2$_{1}$/c")

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


if __name__ == "__main__":
    unittest.main()
