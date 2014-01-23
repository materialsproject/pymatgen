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
        pass

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
