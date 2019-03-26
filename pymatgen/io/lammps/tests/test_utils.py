# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
import os
import unittest

from pymatgen import Molecule
from pymatgen.io.lammps.utils import Polymer
from pymatgen.io.lammps.data import Topology


__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class TestPolymer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # head molecule
        cls.peo_head = Molecule.from_file(os.path.join(test_dir, "peo_head.xyz"))
        charges = [-0.1187, 0.0861, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861, 0.0861]
        cls.peo_head.add_site_property("charge", charges)
        s_head = 0
        s_tail = 5

        # chain molecule
        cls.peo_bulk = Molecule.from_file(os.path.join(test_dir, "peo_bulk.xyz"))
        charges = [-0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861, 0.0861]
        cls.peo_bulk.add_site_property("charge", charges)
        head = 0
        tail = 4

        # terminal molecule
        cls.peo_tail = Molecule.from_file(os.path.join(test_dir, "peo_tail.xyz"))
        charges = [-0.0326, 0.0861, 0.0861, -0.2792, -0.1187, 0.0861, 0.0861, 0.0861]
        cls.peo_tail.add_site_property("charge", charges)
        e_head = 0
        e_tail = 4

        cls.n_units = 25
        link_distance = 1.5075

        # create the polymer
        cls.peo_polymer = Polymer(cls.peo_head, s_head, s_tail,
                                  cls.peo_bulk, head, tail,
                                  cls.peo_tail, e_head, e_tail,
                                  cls.n_units, link_distance)

        # linear chain
        cls.peo_polymer_linear = Polymer(cls.peo_head, s_head, s_tail,
                                         cls.peo_bulk, head, tail,
                                         cls.peo_tail, e_head, e_tail,
                                         cls.n_units, link_distance, linear_chain=True)

    def test_polymer_chain_lengths(self):
        self.assertEqual(len(self.peo_polymer.molecule),
                         len(self.peo_head)+
                         (self.n_units-2)*len(self.peo_bulk)+
                         len(self.peo_tail))
        self.assertEqual(len(self.peo_polymer.molecule),
                         len(self.peo_polymer_linear.molecule))

    def test_polymer_chain_topologies(self):
        topology_random = Topology.from_bonding(self.peo_polymer.molecule)
        topology_linear = Topology.from_bonding(self.peo_polymer_linear.molecule)
        self.assertNotEqual(topology_linear.topologies["Bonds"],
                            topology_random.topologies["Bonds"])
        self.assertNotEqual(topology_linear.topologies["Angles"],
                            topology_random.topologies["Angles"])
        self.assertNotEqual(topology_linear.topologies["Dihedrals"],
                            topology_random.topologies["Dihedrals"])


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
