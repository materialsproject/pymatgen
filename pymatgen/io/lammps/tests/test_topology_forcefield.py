# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

from pymatgen.io.lammps.force_field import ForceField
from pymatgen.io.lammps.topology import Topology
from pymatgen.core.structure import Molecule

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'
__credits__ = 'Brandon Wood'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class TestLammpsTopology(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        polymer = Molecule.from_file(os.path.join(test_dir, "polymer.xyz"))
        cls.monomer = Molecule.from_file(os.path.join(test_dir, "monomer.xyz"))
        cls.topology = Topology.from_molecule(polymer, tol=0.1)
        cls.forcefield = ForceField.from_file(os.path.join(test_dir,
                                                       "ff_data.yaml"))

    def test_topology(self):
        atoms = [['H', 'H'], ['H', 'H'], ['H', 'H'], ['H', 'H'], ['H', 'H'],
                 ['C', 'C'], ['O', 'O'], ['C', 'C'], ['C', 'C'], ['O', 'O']]
        bonds = [[11, 19, ('C', 'H')], [12, 13, ('O', 'C')],
                 [13, 20, ('C', 'H')], [13, 21, ('C', 'H')],
                 [19, 22, ('H', 'C')], [22, 23, ('C', 'O')],
                 [22, 29, ('C', 'H')], [23, 24, ('O', 'C')],
                 [24, 25, ('C', 'C')], [24, 30, ('C', 'H')]]
        angles = [[8, 9, 10, ('C', 'O', 'C')], [9, 10, 11, ('O', 'C', 'C')],
                  [9, 10, 16, ('O', 'C', 'H')], [9, 10, 17, ('O', 'C', 'H')],
                  [11, 10, 16, ('C', 'C', 'H')], [11, 10, 17, ('C', 'C', 'H')],
                  [10, 11, 12, ('C', 'C', 'O')], [10, 11, 18, ('C', 'C', 'H')],
                  [10, 11, 19, ('C', 'C', 'H')], [16, 10, 17, ('H', 'C', 'H')]]
        dihedrals = [[8, 9, 10, 16, ('C', 'O', 'C', 'H')],
                     [8, 9, 10, 17, ('C', 'O', 'C', 'H')],
                     [15, 8, 9, 10, ('H', 'C', 'O', 'C')],
                     [9, 10, 11, 12, ('O', 'C', 'C', 'O')],
                     [9, 10, 11, 18, ('O', 'C', 'C', 'H')],
                     [9, 10, 11, 19, ('O', 'C', 'C', 'H')],
                     [10, 11, 12, 13, ('C', 'C', 'O', 'C')],
                     [10, 11, 19, 22, ('C', 'C', 'H', 'C')],
                     [16, 10, 11, 12, ('H', 'C', 'C', 'O')],
                     [16, 10, 11, 18, ('H', 'C', 'C', 'H')]]
        self.assertEqual(len(self.topology.atoms), 44)
        self.assertEqual(len(self.topology.bonds), 50)
        self.assertEqual(len(self.topology.angles), 123)
        self.assertEqual(len(self.topology.dihedrals), 175)
        self.assertEqual(self.topology.atoms[17:27], atoms)
        self.assertEqual(self.topology.bonds[17:27], bonds)
        self.assertEqual(self.topology.angles[17:27], angles)
        self.assertEqual(self.topology.dihedrals[17:27], dihedrals)
        monomer_charges = [-0.118700, -0.279200, -0.032600, -0.118700,
                          -0.279200, -0.032600, 0.086100, 0.086100,
                           0.086100, 0.086100,0.086100, 0.086100,
                           0.086100, 0.086100]
        monomer_atoms = [['C', 'C'], ['O', 'O'], ['C', 'C'], ['C', 'C'],
                         ['O', 'O'], ['C', 'C'], ['H', 'H'], ['H', 'H'],
                         ['H', 'H'], ['H', 'H'], ['H', 'H'], ['H', 'H'],
                         ['H', 'H'], ['H', 'H']]
        self.monomer.add_site_property("charge", monomer_charges)
        monomer_topology = Topology.from_molecule(self.monomer, tol=0.1)
        self.assertEqual(monomer_topology.atoms, monomer_atoms)
        self.assertEqual(monomer_topology.charges, monomer_charges)

    def test_forcefield(self):
        atoms = {'H': 'H', 'C': 'C', 'O': 'O'}
        bonds = {(u'C', u'O'): [1000, 1.4115], (u'C', u'H'): [1000, 1.1041],
                 (u'C', u'C'): [1000, 1.5075]}
        angles = {(u'C', u'C', u'H'): [42.9, 110.1],
                  (u'H', u'C', u'H'): [38.5, 109.47],
                  (u'H', u'C', u'O'): [56.0, 109.48],
                  (u'C', u'C', u'O'): [86.0, 108.54],
                  (u'C', u'O', u'C'): [74.5, 108.05]}
        dihedrals = {(u'H', u'C', u'O', u'C'): [0.0, 0.0, -0.73, 0.0],
                     (u'H', u'C', u'C', u'O'): [0.0, 0.0, 0.28, 0.0],
                     (u'C', u'C', u'O', u'C'): [1.76, 0.67, 0.04, 0.0],
                     (u'H', u'C', u'C', u'H'): [0.0, 0.0, 0.28, 0.0],
                     (u'O', u'C', u'C', u'O'): [0.41, -2.1, -0.6, -0.82]}
        pairs = {(u'C', u'C'): [14976.0, 0.3236, 637.6],
                 (u'H', u'H'): [2649.6, 0.2674, 27.22],
                 (u'O', u'O'): [75844.8, 0.2461, 396.9]}
        self.assertEqual(self.forcefield.atoms, atoms)
        self.assertEqual(self.forcefield.bonds, bonds)
        self.assertEqual(self.forcefield.angles, angles)
        self.assertEqual(self.forcefield.dihedrals, dihedrals)
        self.assertEqual(self.forcefield.pairs, pairs)


if __name__ == "__main__":
    unittest.main()
