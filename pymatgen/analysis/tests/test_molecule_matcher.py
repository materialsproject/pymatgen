# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os

from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Molecule

try:
    import openbabel as ob
    from pymatgen.analysis.molecule_matcher import MoleculeMatcher
    from pymatgen.analysis.molecule_matcher import IsomorphismMolAtomMapper
    from pymatgen.analysis.molecule_matcher import InchiMolAtomMapper
except (ImportError, RuntimeError):
    ob = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "molecules", "molecule_matcher")

obalign_missing = (ob is None) or ('OBAlign' not in dir(ob))


@unittest.skipIf(obalign_missing, "OBAlign is missing, Skipping")
class MoleculeMatcherTest(unittest.TestCase):

    def test_fit(self):
        self.fit_with_mapper(IsomorphismMolAtomMapper())
        self.fit_with_mapper(InchiMolAtomMapper())

    def test_get_rmsd(self):
        mm = MoleculeMatcher()
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        self.assertEqual('{0:7.3}'.format(mm.get_rmsd(mol1, mol2)), "0.00488")

    def test_group_molecules(self):
        mm = MoleculeMatcher(tolerance=0.001)
        with open(os.path.join(test_dir, "mol_list.txt")) as f:
            filename_list = [line.strip() for line in f.readlines()]
        mol_list = [Molecule.from_file(os.path.join(test_dir, f))
                    for f in filename_list]
        mol_groups = mm.group_molecules(mol_list)
        filename_groups = [[filename_list[mol_list.index(m)] for m in g]
                           for g in mol_groups]
        with open(os.path.join(test_dir, "grouped_mol_list.txt")) as f:
            grouped_text = f.read().strip()
        self.assertEqual(str(filename_groups), grouped_text)

    def test_to_and_from_dict(self):
        mm = MoleculeMatcher(tolerance=0.5,
                             mapper=InchiMolAtomMapper(angle_tolerance=50.0))
        d = mm.as_dict()
        mm2 = MoleculeMatcher.from_dict(d)
        self.assertEqual(d, mm2.as_dict())

        mm = MoleculeMatcher(tolerance=0.5, mapper=IsomorphismMolAtomMapper())
        d = mm.as_dict()
        mm2 = MoleculeMatcher.from_dict(d)
        self.assertEqual(d, mm2.as_dict())

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

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "benzene2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "benzene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))
        self.assertFalse(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "c1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "c2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "j1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "j2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "ethene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "ethene2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "toluene1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "toluene2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "cyclohexane1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cyclohexane2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mol1 = Molecule.from_file(os.path.join(test_dir, "oxygen1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "oxygen2.xyz"))
        self.assertTrue(mm.fit(mol1, mol2))

        mm = MoleculeMatcher(tolerance=0.001, mapper=mapper)
        mol1 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t4.xyz"))
        self.assertFalse(mm.fit(mol1, mol2))

    def test_strange_inchi(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "k1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "k2.sdf"))
        self.assertTrue(mm.fit(mol1, mol2))

    def test_thiane(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "thiane1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "thiane2.sdf"))
        self.assertFalse(mm.fit(mol1, mol2))

    def test_thiane_ethynyl(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "thiane_ethynyl1.sdf"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "thiane_ethynyl2.sdf"))
        self.assertFalse(mm.fit(mol1, mol2))

    def test_cdi_23(self):
        mm = MoleculeMatcher(tolerance=0.05, mapper=InchiMolAtomMapper())
        mol1 = Molecule.from_file(os.path.join(test_dir, "cdi_23_1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "cdi_23_2.xyz"))
        self.assertFalse(mm.fit(mol1, mol2))


if __name__ == '__main__':
    unittest.main()
