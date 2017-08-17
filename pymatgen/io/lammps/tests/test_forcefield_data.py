# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest
from collections import OrderedDict

import numpy as np

from pymatgen.io.lammps.data import LammpsForceFieldData
from pymatgen.io.lammps.force_field import ForceField
from pymatgen.io.lammps.topology import Topology
from pymatgen.core.structure import Molecule

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class TestLammpsForceFieldData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        polymer_chain = Molecule.from_file(os.path.join(test_dir, "polymer_chain.xyz"))
        polymer_linear = Molecule.from_file(os.path.join(test_dir, "polymer_linear.xyz"))
        cls.polymer_matrix = Molecule.from_file(os.path.join(test_dir, "polymer_matrix.xyz"))

        charges = [-0.1187, 0.0861, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.1187, 0.0861,
                   0.0861, 0.0861]
        polymer_linear.add_site_property("charge", charges)
        topology = Topology.from_molecule(polymer_linear)

        cls.polymer_linear_ff_decorated = Molecule.from_file(
            os.path.join(test_dir,"polymer_linear.xyz"))
        ff_map = ['C2', 'H3', 'H2', 'H3', 'O', 'C3', 'H2', 'H3', 'C2', 'H3',
                  'H2', 'O', 'C2', 'H3', 'H2', 'C3', 'H2', 'H3', 'O', 'C3',
                  'H2', 'H3', 'C2', 'H3', 'H2', 'O', 'C2', 'H3', 'H2', 'C3',
                  'H2', 'H3', 'O', 'C3', 'H2', 'H3', 'C2', 'H3', 'H2', 'O',
                  'C2', 'H3', 'H2', 'C3', 'H2', 'H3', 'O', 'C3', 'H2', 'H3', 'H2']
        cls.polymer_linear_ff_decorated.add_site_property("ff_map", ff_map)

        atoms = OrderedDict([("C", "C"), ("H", "H"), ("O", "O")])
        bonds = OrderedDict([((u'C', u'O'), [1000, 1.4115]),
                             ((u'C', u'H'), [1000, 1.1041]),
                             ((u'C', u'C'), [1000, 1.5075])])
        pairs = OrderedDict([((u'O', u'O'), [75844.8, 0.2461, 396.9]),
                             ((u'H', u'H'), [2649.6, 0.2674, 27.22]),
                             ((u'C', u'C'), [14976.0, 0.3236, 637.6])])
        angles = OrderedDict([((u'C', u'C', u'H'), [42.9, 110.1]),
                              ((u'H', u'C', u'H'), [38.5, 109.47]),
                              ((u'H', u'C', u'O'), [56.0, 109.48]),
                              ((u'C', u'C', u'O'), [86.0, 108.54]),
                              ((u'C', u'O', u'C'), [74.5, 108.05])])
        dihedrals = OrderedDict([((u'H', u'C', u'O', u'C'), [0.0, 0.0, -0.73, 0.0]),
                                 ((u'H', u'C', u'C', u'H'), [0.0, 0.0, 0.28, 0.0]),
                                 ((u'C', u'C', u'O', u'C'), [1.76, 0.67, 0.04, 0.0]),
                                 ((u'H', u'C', u'C', u'O'), [0.0, 0.0, 0.28, 0.0]),
                                 ((u'O', u'C', u'C', u'O'), [0.41, -2.1, -0.6, -0.82])])
        forcefield = ForceField(atoms, bonds, angles, dihedrals=dihedrals, pairs=pairs)

        cls.molecules = [polymer_chain] * 3
        cls.mols_number = [7, 3, 1]
        box_size = [[0.0, 50], [0.0, 50], [0.0, 50]]
        cls.topologies = [topology] * len(cls.molecules)

        cls.lammps_ff_data_1 = LammpsForceFieldData.from_forcefield_and_topology(
            cls.molecules, cls.mols_number, box_size, cls.polymer_matrix,
            forcefield, cls.topologies)

    def test_system_info(self):
        # check te molecule ids
        mol_ids = np.array(self.lammps_ff_data_1.atoms_data)[:, 1]
        mol_ids_ans = [i + 1 for i in range(sum(self.mols_number))]
        self.assertEqual(set(mol_ids.tolist()), set(mol_ids_ans))
        # check the size consistency of the polymer matrix
        self.assertEqual(len(self.polymer_matrix),
                         sum([len(mol) * self.mols_number[i]
                              for i, mol in enumerate(self.molecules)]))
        for top in self.topologies:
            self.assertEqual(len(self.lammps_ff_data_1.atoms_data),
                             sum([len(top.atoms)*mol_number
                                  for mol_number in self.mols_number]))
            self.assertEqual(len(self.lammps_ff_data_1.bonds_data),
                             sum([len(top.bonds) * mol_number
                                  for mol_number in self.mols_number]))
            self.assertEqual(len(self.lammps_ff_data_1.angles_data),
                             sum([len(top.angles) * mol_number
                                  for mol_number in self.mols_number]))
            self.assertEqual(len(self.lammps_ff_data_1.dihedrals_data),
                             sum([len(top.dihedrals) * mol_number
                                  for mol_number in self.mols_number]))

    def test_system_info_with_ff_map(self):
        natoms, natom_types, atomic_masses_dict = \
            LammpsForceFieldData.get_basic_system_info(self.polymer_linear_ff_decorated)

        ans_atom_masses_dict = OrderedDict([('H2', [1, 1.00794]),
                                            ('H3', [2, 1.00794]),
                                            ('C2', [3, 12.0107]),
                                            ('C3', [4, 12.0107]),
                                            ('O', [5, 15.9994])])

        self.assertEqual(natoms, 51)
        self.assertEqual(natom_types, 5)
        self.assertEqual(atomic_masses_dict.keys(), ans_atom_masses_dict.keys())
        for k, v in atomic_masses_dict.items():
            self.assertEqual(ans_atom_masses_dict[k][1], v[1])

    def test_atoms_data_with_ff_map(self):
        natoms, natom_types, atomic_masses_dict = \
            LammpsForceFieldData.get_basic_system_info(self.polymer_linear_ff_decorated)
        topology = Topology.from_molecule(self.polymer_linear_ff_decorated,
                                          ff_map="ff_map")
        atoms_data, molid_to_atomid = LammpsForceFieldData.get_atoms_data(
            [self.polymer_linear_ff_decorated], [1],
            self.polymer_linear_ff_decorated, atomic_masses_dict,
            [topology])

        for i, atom in enumerate(atoms_data):
            key = self.polymer_linear_ff_decorated[i].ff_map
            self.assertEqual(atom[2], atomic_masses_dict[key][0])

    def test_to_and_from_file(self):
        self.lammps_ff_data_1.write_data_file(
            os.path.join(test_dir,"lammps_ff_data.dat"))
        lammps_ff_data_2 = LammpsForceFieldData.from_file(
            os.path.join(test_dir,"lammps_ff_data.dat"))
        np.testing.assert_almost_equal(self.lammps_ff_data_1.bond_coeffs,
                         lammps_ff_data_2.bond_coeffs)
        np.testing.assert_almost_equal(self.lammps_ff_data_1.pair_coeffs,
                         lammps_ff_data_2.pair_coeffs)
        np.testing.assert_almost_equal(self.lammps_ff_data_1.angle_coeffs,
                         lammps_ff_data_2.angle_coeffs)
        np.testing.assert_almost_equal(self.lammps_ff_data_1.dihedral_coeffs,
                         lammps_ff_data_2.dihedral_coeffs)
        np.testing.assert_almost_equal(self.lammps_ff_data_1.atoms_data,
                                       lammps_ff_data_2.atoms_data, decimal=10)
        self.assertEqual(self.lammps_ff_data_1.bonds_data,
                         lammps_ff_data_2.bonds_data)
        self.assertEqual(self.lammps_ff_data_1.angles_data,
                         lammps_ff_data_2.angles_data)
        self.assertEqual(self.lammps_ff_data_1.dihedrals_data,
                         lammps_ff_data_2.dihedrals_data)

    def tearDown(self):
        for x in ["lammps_ff_data.dat"]:
            if os.path.exists(os.path.join(test_dir, x)):
                os.remove(os.path.join(test_dir, x))


if __name__ == "__main__":
    unittest.main()
