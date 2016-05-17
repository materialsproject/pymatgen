# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on Apr 28, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 28, 2012"

import unittest2 as unittest
import os

from pymatgen.core.structure import Molecule
from pymatgen.io.xyz import XYZ
from pymatgen.io.babel import BabelMolAdaptor

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "molecules")

try:
    import openbabel as ob
    import pybel as pb
except ImportError:
    pb = None
    ob = None


@unittest.skipIf(not (pb and ob), "OpenBabel not present. Skipping...")
class BabelMolAdaptorTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)

    def test_init(self):
        adaptor = BabelMolAdaptor(self.mol)
        obmol = adaptor.openbabel_mol
        self.assertEqual(obmol.NumAtoms(), 5)

        adaptor = BabelMolAdaptor(adaptor.openbabel_mol)
        self.assertEqual(adaptor.pymatgen_mol.formula, "H4 C1")

    def test_from_file(self):
        adaptor = BabelMolAdaptor.from_file(
            os.path.join(test_dir, "Ethane_e.pdb"), "pdb")
        mol = adaptor.pymatgen_mol
        self.assertEqual(mol.formula, "H6 C2")

    def test_from_string(self):
        xyz = XYZ(self.mol)
        adaptor = BabelMolAdaptor.from_string(str(xyz), "xyz")
        mol = adaptor.pymatgen_mol
        self.assertEqual(mol.formula, "H4 C1")

    def test_localopt(self):
        self.mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(self.mol)
        adaptor.localopt()
        optmol = adaptor.pymatgen_mol
        for site in optmol[1:]:
            self.assertAlmostEqual(site.distance(optmol[0]), 1.09216, 2)

if __name__ == "__main__":
    unittest.main()
