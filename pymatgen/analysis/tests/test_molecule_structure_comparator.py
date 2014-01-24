import os
from unittest import TestCase
from pymatgen.analysis.molecule_structure_comparator import \
    MoleculeStructureComparator
from pymatgen.io.smartio import read_mol

__author__ = 'xiaohuiqu'


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "molecules", "structural_change")


class TestMoleculeStructureComparator(TestCase):
    def test_are_equal(self):
        msc1 = MoleculeStructureComparator()
        mol1 = read_mol(os.path.join(test_dir, "t1.xyz"))
        mol2 = read_mol(os.path.join(test_dir, "t2.xyz"))
        mol3 = read_mol(os.path.join(test_dir, "t3.xyz"))
        self.assertFalse(msc1.are_equal(mol1, mol2))
        self.assertTrue(msc1.are_equal(mol2, mol3))
        thio1 = read_mol(os.path.join(test_dir, "thiophene1.xyz"))
        thio2 = read_mol(os.path.join(test_dir, "thiophene2.xyz"))
        # noinspection PyProtectedMember
        msc2 = MoleculeStructureComparator(
            priority_bonds=msc1._get_bonds(thio1))
        self.assertTrue(msc2.are_equal(thio1, thio2))

    def test_get_bonds(self):
        mol = read_mol(os.path.join(test_dir, "t1.xyz"))
        msc = MoleculeStructureComparator()
        # noinspection PyProtectedMember
        bonds = msc._get_bonds(mol)
        bonds_ref = [(0, 1), (0, 2), (0, 3), (0, 23), (3, 4), (3, 5), (5, 6),
                     (5, 7), (7, 8), (7, 9), (7, 21), (9, 10), (9, 11),
                     (9, 12), (12, 13), (12, 14), (12, 15), (15, 16), (15, 17),
                     (15, 18), (18, 19), (18, 20), (18, 21), (21, 22),
                     (21, 23), (23, 24), (23, 25)]
        self.assertEqual(bonds, bonds_ref)

    def test_to_and_from_dict(self):
        msc1 = MoleculeStructureComparator()
        d1 = msc1.to_dict
        d2 = MoleculeStructureComparator.from_dict(d1).to_dict
        self.assertEqual(d1, d2)
        thio1 = read_mol(os.path.join(test_dir, "thiophene1.xyz"))
        # noinspection PyProtectedMember
        msc2 = MoleculeStructureComparator(
            bond_length_cap=0.2,
            priority_bonds=msc1._get_bonds(thio1),
            priority_cap=0.5)
        d1 = msc2.to_dict
        d2 = MoleculeStructureComparator.from_dict(d1).to_dict
        self.assertEqual(d1, d2)
