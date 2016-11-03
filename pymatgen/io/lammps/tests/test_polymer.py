# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

from pymatgen import Molecule
from pymatgen.io.lammps.utils import Polymer
from pymatgen.io.lammps.topology import Topology


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
        topology_random = Topology.from_molecule(self.peo_polymer.molecule, tol=0.1)
        topology_linear = Topology.from_molecule(self.peo_polymer_linear.molecule, tol=0.1)
        self.assertNotEqual(topology_linear.bonds, topology_random.bonds)
        self.assertNotEqual(topology_linear.angles, topology_random.angles)
        self.assertNotEqual(topology_linear.dihedrals, topology_random.dihedrals)


if __name__ == "__main__":
    unittest.main()
