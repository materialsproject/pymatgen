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
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Molecule
try:
    import openbabel as ob
except ImportError:
    ob = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "molecules")


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
        mm = MoleculeMatcher()
        print mm.fit(mol1, mol2)
        self.assertTrue(mm.fit(mol1, mol2))


if __name__ == '__main__':
    unittest.main()