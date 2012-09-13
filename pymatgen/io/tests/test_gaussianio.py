#!/usr/bin/env python

'''
Created on Apr 17, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 17, 2012"

import unittest
import os

from pymatgen import Molecule, __file__
from pymatgen.io.gaussianio import GaussianInput, GaussianOutput

test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                        'test_files')


class GaussianInputTest(unittest.TestCase):

    def setUp(self):

        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]

        mol = Molecule(["C", "H", "H", "H", "H"], coords)
        self.gau = GaussianInput(mol,
                                 route_parameters={'SP': "", "SCF": "Tight"},
                                 input_parameters={"EPS": 12})

    def test_str_and_from_string(self):
        ans = """#P HF/6-31G(d) SP SCF=Tight Test

H4 C1

0 1
C
H 1 B1
H 1 B2 2 A2
H 1 B3 2 A3 3 D3
H 1 B4 2 A4 4 D4

B1=1.089000
B2=1.089000
A2=109.471221
B3=1.089000
A3=109.471213
D3=120.000017
B4=1.089000
A4=109.471213
D4=119.999966

EPS=12
"""
        self.assertEqual(str(self.gau), ans)
        gau = GaussianInput.from_string(ans)
        self.assertEqual(gau.functional, 'HF')
        self.assertEqual(gau.input_parameters['EPS'], 12)

    def test_from_file(self):
        filepath = os.path.join(test_dir, 'MethylPyrrolidine_drawn.gjf')
        gau = GaussianInput.from_file(filepath)
        self.assertEqual(gau.molecule.composition.formula, "H11 C5 N1")
        self.assertIn("opt", gau.route_parameters)
        self.assertEqual(gau.route_parameters["geom"], "connectivity")
        self.assertEqual(gau.functional, "b3lyp")
        self.assertEqual(gau.basis_set, "6-311+g(d,p)")


class GaussianOutputTest(unittest.TestCase):

    def setUp(self):
        self.gauout = GaussianOutput(os.path.join(test_dir, "methane.log"))

    def test_props(self):
        gau = self.gauout
        self.assertEqual(len(gau.energies), 3)
        self.assertAlmostEqual(gau.energies[-1], -39.9768775602)
        self.assertEqual(len(gau.structures), 4)
        for mol in gau.structures:
            self.assertEqual(mol.formula, 'H4 C1')
        self.assertIn("OPT", gau.route)
        self.assertEqual("Minimum", gau.stationary_type)
        self.assertEqual("HF", gau.functional)
        self.assertEqual("3-21G", gau.basis_set)
        self.assertEqual(17, gau.num_basis_func)


if __name__ == "__main__":
    unittest.main()
