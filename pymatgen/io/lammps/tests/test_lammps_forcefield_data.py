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
        polymer_chain = Molecule.from_file(os.path.join(test_dir,"polymer_chain.xyz"))
        polymer_linear = Molecule.from_file(os.path.join(test_dir,"polymer_linear.xyz"))
        cls.polymer_matrix = Molecule.from_file(os.path.join(test_dir,"polymer_matrix.xyz"))
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

        atoms = OrderedDict([("C","C"), ("H","H"), ("O", "O")])
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
        forcefield =ForceField(atoms, bonds, angles, dihedrals=dihedrals, pairs=pairs)

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
