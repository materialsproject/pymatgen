# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on Jul 30, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 30, 2012"

import unittest
import os

from pymatgen.io.smart import read_structure, read_mol, write_structure, \
    write_mol
from pymatgen.core.structure import Structure, Molecule
from pymatgen.analysis.structure_matcher import StructureMatcher

try:
    import openbabel as ob
except ImportError:
    ob = None


class MethodsTest(unittest.TestCase):

    def test_read_structure(self):
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files')
        for fname in ("Li2O.cif", "vasprun.xml",
                      "vasprun_Si_bands.xml", "Si.cssr"):
            filename = os.path.join(test_dir, fname)
            struct = read_structure(filename)
            self.assertIsInstance(struct, Structure)
            prim = read_structure(filename, primitive=True)
            self.assertLessEqual(len(prim), len(struct))
            sorted_s = read_structure(filename, sort=True)
            self.assertEqual(sorted_s, sorted_s.get_sorted_structure())

        m = StructureMatcher()
        for ext in [".cif", ".json", ".cssr"]:
            fn = "smartio_structure_test" + ext
            write_structure(struct, fn)
            back = read_structure(fn)
            self.assertTrue(m.fit(back, struct))
            os.remove(fn)

    def test_read_mol(self):
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files', "molecules")
        for fname in ("methane.log", "c60.xyz", "ethane.gjf"):
            filename = os.path.join(test_dir, fname)
            mol = read_mol(filename)
            self.assertIsInstance(mol, Molecule)

        for ext in [".xyz", ".json", ".gjf"]:
            fn = "smartio_mol_test" + ext
            write_mol(mol, fn)
            back = read_mol(fn)
            self.assertEqual(back, mol)
            os.remove(fn)

    @unittest.skipIf(ob is None, "No openbabel")
    def test_read_mol_babel(self):
        test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                'test_files', "molecules")
        for fname in ("ethane.mol", ):
            filename = os.path.join(test_dir, fname)
            mol = read_mol(filename)
            self.assertIsInstance(mol, Molecule)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
