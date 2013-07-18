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
from pymatgen.analysis.molecule_matcher import InchiMolAtomMapper
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Molecule
from pymatgen.io.babelio import BabelMolAdaptor

try:
    import openbabel as ob
    ob.OBAlign
except (ImportError, AttributeError):
    ob = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "molecules", "molecule_matcher")


@unittest.skipIf(ob is None, "OpenBabel not present. Skipping...")
class MoleculeMatcherTest(unittest.TestCase):

    def test_fit(self):
        self.fit_with_mapper(IsomorphismMolAtomMapper())
        self.fit_with_mapper(InchiMolAtomMapper())

    def test_get_rmsd(self):
        mm = MoleculeMatcher()
        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "t3.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "t4.xyz")).pymatgen_mol
        self.assertEqual('{0:7.3}'.format(mm.get_rmsd(mol1, mol2)), "0.00488")

    def test_group_molecules(self):
        mm = MoleculeMatcher(tolerance=0.001)
        filename_list = None
        with open(os.path.join(test_dir, "mol_list.txt")) as f:
            filename_list = [line.strip() for line in f.readlines()]
        mol_list = [BabelMolAdaptor.from_file(os.path.join(test_dir, f)).pymatgen_mol\
                    for f in filename_list]
        mol_groups = mm.group_molecules(mol_list)
        filename_groups = [[filename_list[mol_list.index(m)] for m in g] for g \
                           in mol_groups]
        grouped_text = None
        with open(os.path.join(test_dir, "grouped_mol_list.txt")) as f:
            grouped_text = f.read().strip()
        self.assertEqual(str(filename_groups), grouped_text)

    def test_to_and_from_dict(self):
        mm = MoleculeMatcher(tolerance=0.5, mapper=InchiMolAtomMapper(angle_tolerance=50.0))
        d = mm.to_dict
        mm2 = MoleculeMatcher.from_dict(d)
        self.assertEqual(d, mm2.to_dict)

        mm = MoleculeMatcher(tolerance=0.5, mapper=IsomorphismMolAtomMapper())
        d = mm.to_dict
        mm2 = MoleculeMatcher.from_dict(d)
        self.assertEqual(d, mm2.to_dict)

    
    def fit_with_mapper(self, mapper):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        mol1 = Molecule(["C", "H", "H", "H", "H"], coords)
        op = SymmOp.from_origin_axis_angle([0, 0, 0], [0.1, 0.2, 0.3], 60)
        rotcoords = [op.operate(c) for c in coords]
        mol2 = Molecule(["C", "H", "H", "H", "H"], rotcoords)
        mm = MoleculeMatcher(mapper=mapper)
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

        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "j1.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "j2.xyz")).pymatgen_mol
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "ethene1.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "ethene2.xyz")).pymatgen_mol
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "toluene1.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "toluene2.xyz")).pymatgen_mol
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "cyclohexane1.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "cyclohexane2.xyz")).pymatgen_mol
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "oxygen1.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "oxygen2.xyz")).pymatgen_mol
        self.assertTrue(mm.fit(mol1, mol2))

        mm = MoleculeMatcher(tolerance=0.001, mapper=mapper)
        mol1 = BabelMolAdaptor.from_file(os.path.join(test_dir, "t3.xyz")).pymatgen_mol
        mol2 = BabelMolAdaptor.from_file(os.path.join(test_dir, "t4.xyz")).pymatgen_mol
        self.assertFalse(mm.fit(mol1, mol2))


if __name__ == '__main__':
    unittest.main()
