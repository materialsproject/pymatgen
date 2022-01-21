# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest
from unittest import TestCase

from pymatgen.analysis.molecule_structure_comparator import MoleculeStructureComparator
from pymatgen.core.structure import Molecule
from pymatgen.util.testing import PymatgenTest

__author__ = "xiaohuiqu"


test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "structural_change")


class TestMoleculeStructureComparator(TestCase):
    def test_are_equal(self):
        msc1 = MoleculeStructureComparator()
        mol1 = Molecule.from_file(os.path.join(test_dir, "t1.xyz"))
        mol2 = Molecule.from_file(os.path.join(test_dir, "t2.xyz"))
        mol3 = Molecule.from_file(os.path.join(test_dir, "t3.xyz"))
        self.assertFalse(msc1.are_equal(mol1, mol2))
        self.assertTrue(msc1.are_equal(mol2, mol3))
        thio1 = Molecule.from_file(os.path.join(test_dir, "thiophene1.xyz"))
        thio2 = Molecule.from_file(os.path.join(test_dir, "thiophene2.xyz"))
        # noinspection PyProtectedMember
        msc2 = MoleculeStructureComparator(priority_bonds=msc1._get_bonds(thio1))
        self.assertTrue(msc2.are_equal(thio1, thio2))
        hal1 = Molecule.from_file(os.path.join(test_dir, "molecule_with_halogen_bonds_1.xyz"))
        hal2 = Molecule.from_file(os.path.join(test_dir, "molecule_with_halogen_bonds_2.xyz"))
        msc3 = MoleculeStructureComparator(priority_bonds=msc1._get_bonds(hal1))
        self.assertTrue(msc3.are_equal(hal1, hal2))

    def test_get_bonds(self):
        mol1 = Molecule.from_file(os.path.join(test_dir, "t1.xyz"))
        msc = MoleculeStructureComparator()
        # noinspection PyProtectedMember
        bonds = msc._get_bonds(mol1)
        bonds_ref = [
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 23),
            (3, 4),
            (3, 5),
            (5, 6),
            (5, 7),
            (7, 8),
            (7, 9),
            (7, 21),
            (9, 10),
            (9, 11),
            (9, 12),
            (12, 13),
            (12, 14),
            (12, 15),
            (15, 16),
            (15, 17),
            (15, 18),
            (18, 19),
            (18, 20),
            (18, 21),
            (21, 22),
            (21, 23),
            (23, 24),
            (23, 25),
        ]
        self.assertEqual(bonds, bonds_ref)
        mol2 = Molecule.from_file(os.path.join(test_dir, "MgBH42.xyz"))
        bonds = msc._get_bonds(mol2)
        self.assertEqual(bonds, [(1, 3), (2, 3), (3, 4), (3, 5), (6, 8), (7, 8), (8, 9), (8, 10)])
        msc = MoleculeStructureComparator(ignore_ionic_bond=False)
        bonds = msc._get_bonds(mol2)
        self.assertEqual(
            bonds,
            [
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 5),
                (0, 6),
                (0, 7),
                (0, 8),
                (0, 9),
                (1, 3),
                (2, 3),
                (3, 4),
                (3, 5),
                (6, 8),
                (7, 8),
                (8, 9),
                (8, 10),
            ],
        )

        mol1 = Molecule.from_file(os.path.join(test_dir, "molecule_with_halogen_bonds_1.xyz"))
        msc = MoleculeStructureComparator()
        # noinspection PyProtectedMember
        bonds = msc._get_bonds(mol1)
        self.assertEqual(
            bonds,
            [
                (0, 12),
                (0, 13),
                (0, 14),
                (0, 15),
                (1, 12),
                (1, 16),
                (1, 17),
                (1, 18),
                (2, 4),
                (2, 11),
                (2, 19),
                (3, 5),
                (3, 10),
                (3, 20),
                (4, 6),
                (4, 10),
                (5, 11),
                (5, 12),
                (6, 7),
                (6, 8),
                (6, 9),
            ],
        )

    def test_to_and_from_dict(self):
        msc1 = MoleculeStructureComparator()
        d1 = msc1.as_dict()
        d2 = MoleculeStructureComparator.from_dict(d1).as_dict()
        self.assertEqual(d1, d2)
        thio1 = Molecule.from_file(os.path.join(test_dir, "thiophene1.xyz"))
        # noinspection PyProtectedMember
        msc2 = MoleculeStructureComparator(bond_length_cap=0.2, priority_bonds=msc1._get_bonds(thio1), priority_cap=0.5)
        d1 = msc2.as_dict()
        d2 = MoleculeStructureComparator.from_dict(d1).as_dict()
        self.assertEqual(d1, d2)

    # def test_structural_change_in_geom_opt(self):
    #     qcout_path = os.path.join(test_dir, "mol_1_3_bond.qcout")
    #     qcout = QcOutput(qcout_path)
    #     mol1 = qcout.data[0]["molecules"][0]
    #     mol2 = qcout.data[0]["molecules"][-1]
    #     priority_bonds = [[0, 1], [0, 2], [1, 3], [1, 4], [1, 7], [2, 5], [2, 6], [2, 8], [4, 6], [4, 10], [6, 9]]
    #     msc = MoleculeStructureComparator(priority_bonds=priority_bonds)
    #     self.assertTrue(msc.are_equal(mol1, mol2))

    def test_get_13_bonds(self):
        priority_bonds = [
            [0, 1],
            [0, 2],
            [1, 3],
            [1, 4],
            [1, 7],
            [2, 5],
            [2, 6],
            [2, 8],
            [4, 6],
            [4, 10],
            [6, 9],
        ]
        bonds_13 = MoleculeStructureComparator.get_13_bonds(priority_bonds)
        ans = (
            (0, 3),
            (0, 4),
            (0, 5),
            (0, 6),
            (0, 7),
            (0, 8),
            (1, 2),
            (1, 6),
            (1, 10),
            (2, 4),
            (2, 9),
            (3, 4),
            (3, 7),
            (4, 7),
            (4, 9),
            (5, 6),
            (5, 8),
            (6, 8),
            (6, 10),
        )
        self.assertEqual(bonds_13, tuple(ans))


if __name__ == "__main__":
    unittest.main()
