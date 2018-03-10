# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest
import os

from pymatgen.io.pwscf import PWInput, PWInputError, PWOutput
from pymatgen.util.testing import PymatgenTest

__author__ = 'Shyue Ping Ong'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'ongsp@ucsd.edu'
__date__ = '3/28/15'


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

class PWInputTest(PymatgenTest):

    def test_init(self):
        s = self.get_structure("Li2O")
        self.assertRaises(
            PWInputError, PWInput,
            s,
            control={"calculation": "scf", "pseudo_dir": './'},
            pseudo={"Li": "Li.pbe-n-kjpaw_psl.0.1.UPF"}
        )

    def test_str(self):
        s = self.get_structure("Li2O")

        pw = PWInput(s,
                     control={"calculation": "scf", "pseudo_dir": './'},
                     pseudo={"Li": "Li.pbe-n-kjpaw_psl.0.1.UPF",
                             "O": "O.pbe-n-kjpaw_psl.0.1.UPF"},
                     system={"ecutwfc": 50})
        ans = """&CONTROL
  calculation = 'scf',
  pseudo_dir = './',
/
&SYSTEM
  ecutwfc = 50,
  ibrav = 0,
  nat = 3,
  ntyp = 2,
/
&ELECTRONS
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  Li  6.9410 Li.pbe-n-kjpaw_psl.0.1.UPF
  O  15.9994 O.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS crystal
  O 0.000000 0.000000 0.000000
  Li 0.750178 0.750178 0.750178
  Li 0.249822 0.249822 0.249822
K_POINTS automatic
  1 1 1 0 0 0
CELL_PARAMETERS angstrom
  2.917389 0.097894 1.520005
  0.964634 2.755036 1.520005
  0.133206 0.097894 3.286918
"""
        self.assertEqual(pw.__str__().strip(), ans.strip())


class PWOuputTest(PymatgenTest):

    def setUp(self):
        self.pwout = PWOutput(os.path.join(test_dir, "Si.pwscf.out"))

    def test_properties(self):
        self.assertAlmostEqual(self.pwout.final_energy, -93.45259708)

    def test_get_celldm(self):
        self.assertAlmostEqual(self.pwout.get_celldm(1), 10.323)
        for i in range(2, 7):
            self.assertAlmostEqual(self.pwout.get_celldm(i), 0)

if __name__ == '__main__':
    unittest.main()
