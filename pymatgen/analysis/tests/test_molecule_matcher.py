#!/usr/bin/env python

"""
TODO: Modify module doc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "6/9/13"

import unittest
import os

from pymatgen.analysis.molecule_matcher import MoleculeMatcher
from pymatgen.analysis.molecule_matcher import IsomorphismMolAtomMapper
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Molecule
from pymatgen.io.babelio import BabelMolAdaptor

try:
    import openbabel as ob
except ImportError:
    ob = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "molecules", "molecule_matcher")


@unittest.skipIf(ob is None, "OpenBabel not present. Skipping...")
class MoleculeMatcherTest(unittest.TestCase):

    def test_fit(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        mol1 = Molecule(["C", "H", "H", "H", "H"], coords)
        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)
        rotcoords = [op.operate(c) for c in coords]
        mol2 = Molecule(["C", "H", "H", "H", "H"], rotcoords)
        mm = MoleculeMatcher(mapper=IsomorphismMolAtomMapper())
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "benzene1.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "benzene2.xyz")).pymatgen_mol
        self.assertTrue(mm.fit(mol1, mol2))
        
        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "benzene1.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "t2.xyz")).pymatgen_mol
        self.assertFalse(mm.fit(mol1, mol2))
        
        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "c1.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "c2.xyz")).pymatgen_mol
        self.assertTrue(mm.fit(mol1, mol2))
        
        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "t3.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "t4.xyz")).pymatgen_mol
        self.assertTrue(mm.fit(mol1, mol2))
        
        mm = MoleculeMatcher(tolerance=0.001, mapper=IsomorphismMolAtomMapper())
        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "t3.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "t4.xyz")).pymatgen_mol
        self.assertFalse(mm.fit(mol1, mol2))



if __name__ == '__main__':
    unittest.main()
