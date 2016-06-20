# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

from pymatgen.core import Molecule

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class TestPackmolOutput(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ethanol_coords = [[0.00720, -0.56870, 0.00000],
                          [-1.28540, 0.24990, 0.00000],
                          [1.13040, 0.31470, 0.00000],
                          [0.03920, -1.19720, 0.89000],
                          [0.03920, -1.19720, -0.89000],
                          [-1.31750, 0.87840, 0.89000],
                          [-1.31750, 0.87840, -0.89000],
                          [-2.14220, -0.42390, -0.00000],
                          [1.98570, -0.13650, -0.00000]]
        water_coords = [[9.626, 6.787, 12.673],
                        [9.626, 8.420, 12.673],
                        [10.203, 7.604, 12.673]]
        cls.ethanol_atoms = ["C", "C", "O", "H", "H", "H", "H", "H", "H"]
        cls.water_atoms = ["H", "H", "O"]
        ethanol = Molecule(cls.ethanol_atoms, ethanol_coords)
        water = Molecule(cls.water_atoms, water_coords)
        cls.mols = [ethanol, water]
        cls.cocktail = Molecule.from_file(
            os.path.join(test_dir, "cocktail.xyz"))
        cls.packmol_config = [{"number": 1}, {"number": 15}]

    def test_packed_molecule(self):
        self.assertEqual(len(self.cocktail),
                         sum([len(mol) * self.packmol_config[i]["number"]
                              for i, mol in enumerate(self.mols)]))
        atoms = self.ethanol_atoms * self.packmol_config[0]["number"] + \
                self.water_atoms * self.packmol_config[1]["number"]
        atoms_ans = [str(site.specie) for site in self.cocktail]
        self.assertEqual(atoms, atoms_ans)


if __name__ == "__main__":
    unittest.main()
