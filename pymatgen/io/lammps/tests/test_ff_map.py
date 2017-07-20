from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest
from collections import OrderedDict

import ruamel.yaml as yaml

from pymatgen.io.lammps.data import LammpsForceFieldData
from pymatgen.io.lammps.force_field import ForceField
from pymatgen.io.lammps.topology import Topology
from pymatgen.core.structure import Molecule

__author__ = 'Rishi Gurnani'
__email__ = 'rgurnani96@lbl.gov'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class TestLammpsForceFieldData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        polymer_chain = Molecule.from_file(os.path.join(test_dir, "polymer_chain.xyz"))
        polymer_linear = Molecule.from_file(os.path.join(test_dir, "polymer_linear.xyz"))
        with open(os.path.join(test_dir, "topology_data.yaml")) as f:
            ff_map_test_data = yaml.safe_load(f)
        cls.polymer_matrix = Molecule.from_file(os.path.join(test_dir,"polymer_matrix.xyz"))
        charges = [-0.1187, 0.0861, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.0326, 0.0861,
                   0.0861, -0.0326, 0.0861, 0.0861, -0.2792, -0.1187, 0.0861,
                   0.0861, 0.0861]

        ff_map = ["C3", "H3", "H3", "H3", "O", "C2", "H2",
                   "H2", "C2", "H2", "H2", "O", "C2", "H2",
                   "H2", "C2", "H2", "H2", "O", "C2", "H2",
                   "H2", "C2", "H2", "H2", "O", "C2", "H2",
                   "H2", "C2", "H2", "H2", "O", "C2", "H2",
                   "H2", "C2", "H2", "H2", "O", "C2", "H2",
                   "H2", "C2", "H2", "H2", "O", "C3", "H3", "H3", "H3"]
        polymer_linear.add_site_property("charge", charges)
        polymer_linear.add_site_property("ff_map", ff_map)
        cls.topology = Topology.from_molecule(polymer_linear)

        cls.tatoms = ff_map_test_data['atoms']
        cls.tbonds = ff_map_test_data['bonds']
        cls.tangles = ff_map_test_data['angles']
        cls.tdihedrals = ff_map_test_data['dihedrals']
        cls.atoms = OrderedDict([("C2","C"), ("C3","C"), ("H2", "H"), ("H3", "H"), ("O", "O")])
        cls.bonds = OrderedDict([((u'C2', u'C2'), [222.5, 1.53]),
                             ((u'C2', u'H2'), [309.0, 1.111]),
                             ((u'C2', u'O'), [360.0, 1.415]),
                             ((u'C3', u'H3'), [322.0, 1.111]),
                             ((u'C3', u'O'), [360.0, 1.415])])
        cls.pairs = OrderedDict([((u'C', u'C'), [-0.056, 2.01, -0.01, 1.9]),
                             ((u'H', u'H'), [-0.035, 1.34]),
                             ((u'O', u'O'), [-0.1, 1.65])])
        cls.angles = OrderedDict([((u'C2', u'C2', u'H2'), [26.5, 110.1]),
                              ((u'C2', u'C2', u'O'), [45.0, 111.5]),
                              ((u'C2', u'O', u'C2'), [95.0, 109.7]),
                              ((u'C3', u'O', u'C2'), [95.0, 109.7]),
                              ((u'H2', u'C2', u'H2'), [35.5, 109.0]),
                              ((u'H2', u'C2', u'O'), [60.0, 109.5]),
                              ((u'H3', u'C3', u'H3'), [35.5, 108.4]),
                              ((u'H3', u'C3', u'O'), [60.0, 109.5])])
        cls.dihedrals = OrderedDict([((u'C2', u'C2', u'O', u'C2'), [0.57, 1, 0, 0.0]),
                                 ((u'C2', u'O', u'C2', u'H2'), [0.284, 3, 0, 0.0]),
                                 ((u'C2', u'O', u'C3', u'H3'), [0.284, 3, 0, 0.0]),
                                 ((u'C3', u'O', u'C2', u'C2'), [0.57, 1, 0, 0.0]),
                                 ((u'H2', u'C2', u'C2', u'H2'), [0.19, 3, 0, 0.0]),
                                 ((u'H2', u'C2', u'O', u'C3'), [0.284, 3, 0, 0.0]),
                                 ((u'H2', u'O', u'C2', u'H2'), [0.284, 3, 0, 0.0]),
                                 ((u'H3', u'O', u'C3', u'H3'), [0.284, 3, 0, 0.0]),
                                 ((u'O', u'C2', u'C2', u'H2'), [0.19, 3, 0, 0.0]),
                                 ((u'O', u'C2', u'C2', u'O'), [1.16, 2, 0, 0.0])])
        cls.forcefield = ForceField.from_file(os.path.join(test_dir, "ffmap_data.yaml"))

        cls.molecules = [polymer_chain] * 3
        cls.mols_number = [7, 3, 1]
        box_size = [[0.0, 50], [0.0, 50], [0.0, 50]]
        cls.topologies = [cls.topology] * len(cls.molecules)

        cls.lammps_ff_data_1 = LammpsForceFieldData.from_forcefield_and_topology(
            cls.molecules, cls.mols_number, box_size, cls.polymer_matrix,
            cls.forcefield, cls.topologies)

    def test_ff_map(self):
        #check atoms
        self.assertEqual(self.tatoms, self.topology.atoms)
        #check pairs
        self.assertEqual(self.pairs, self.forcefield.pairs)
        #check bonds
        self.assertEqual(self.bonds, self.forcefield.bonds)
        self.assertEqual(self.tbonds, self.topology.bonds)
        #check angles
        self.assertEqual(self.angles, self.forcefield.angles)
        self.assertEqual(self.tangles, self.topology.angles)
        #check dihedrals
        self.assertEqual(self.dihedrals, self.forcefield.dihedrals)
        self.assertEqual(self.tdihedrals, self.topology.dihedrals)

    def tearDown(self):
        for x in ["lammps_ff_data.dat"]:
            if os.path.exists(os.path.join(test_dir, x)):
                os.remove(os.path.join(test_dir, x))


if __name__ == "__main__":
    unittest.main()
