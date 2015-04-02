# coding: utf-8

from __future__ import division, unicode_literals

__author__ = 'Shyue Ping Ong'
__copyright__ = 'Copyright 2013, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'ongsp@ucsd.edu'
__date__ = '3/28/15'

import unittest

from pymatgen.io.pwscfio import PWInput, PWInputError
from pymatgen.util.testing import PymatgenTest


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
  Li 6.9410 Li.pbe-n-kjpaw_psl.0.1.UPF
  O 15.9994 O.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS crystal
  Li 0.250000 0.250000 0.250000
  Li 0.750000 0.750000 0.750000
  O 0.000000 0.000000 0.000000
K_POINTS automatic
  1 1 1 0 0 0
CELL_PARAMETERS angstrom
  -2.305000 -2.305000 0.000000
  -2.305000 0.000000 -2.305000
  0.000000 -2.305000 -2.305000
"""
        self.assertEqual(pw.__str__().strip(), ans.strip())


if __name__ == '__main__':
    unittest.main()
