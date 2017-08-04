from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest
from collections import OrderedDict

from pymatgen.io.lammps.force_field import ForceField
from pymatgen.io.lammps.topology import Topology
from pymatgen.core.structure import Molecule
from pymatgen.io.lammps.data import LammpsForceFieldData
from pymatgen.io.lammps.utils import PackmolRunner

__author__ = 'Rishi Gurnani'
__email__ = 'rgurnani96@lbl.gov'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


#DMOE stands for Dimethoxyethane
class TestLammpsForceFieldDataWithMap(unittest.TestCase):

    def setUp(self):
        self.dmoe = Molecule.from_file(os.path.join(test_dir,
                                                              "DMOE.xyz"))

        dmoe_charges = [-0.10, 0.09, 0.09, 0.09,
                        -0.34, -0.01, 0.09, 0.09,
                        -0.01, 0.09, 0.09, -0.34,
                        -0.10, 0.09, 0.09, 0.09]
        self.dmoe.add_site_property("charge", dmoe_charges)
        ff_map = ["C3", "H3", "H3", "H3", "O", "C2", "H2",
                    "H2", "C2", "H2", "H2", "O", "C3", "H3", "H3", "H3"]
        self.dmoe.add_site_property("ff_map", ff_map)
        self.topology = Topology.from_molecule(self.dmoe, ff_map = "ff_map")
        self.forcefield = ForceField.from_file(os.path.join(test_dir,
                                                            "ffmap_data.yaml"))

    def test_topology(self):
        tatoms = [['C', 'C3'],
                    ['H', 'H3'],
                    ['H', 'H3'],
                    ['H', 'H3'],
                    ['O', 'O'],
                    ['C', 'C2'],
                    ['H', 'H2'],
                    ['H', 'H2'],
                    ['C', 'C2'],
                    ['H', 'H2'],
                    ['H', 'H2'],
                    ['O', 'O'],
                    ['C', 'C3'],
                    ['H', 'H3'],
                    ['H', 'H3'],
                    ['H', 'H3']]
        tbonds = [[0, 1, ('C3', 'H3')], [0, 2, ('C3', 'H3')], [0, 3, ('C3', 'H3')],
                    [0, 4, ('C3', 'O')], [4, 5, ('O', 'C2')], [5, 6, ('C2', 'H2')],
                    [5, 7, ('C2', 'H2')], [5, 8, ('C2', 'C2')], [8, 9, ('C2', 'H2')],
                    [8, 10, ('C2', 'H2')], [8, 11, ('C2', 'O')], [11, 12, ('O', 'C3')],
                    [12, 13, ('C3', 'H3')], [12, 14, ('C3', 'H3')], [12, 15, ('C3', 'H3')]]

        self.assertEqual(tatoms, self.topology.atoms)
        self.assertEqual(tbonds, self.topology.bonds)
        # for i, angle in enumerate(self.topology.angles):
        #     for j in range(3):
        #         self.assertEqual(self.polymer_linear[angle[j]].ff_map,
        #                          angle[3][j])
        #         self.assertEqual(self.polymer_linear[angle[j]].specie.symbol,
        #                          angle[3][j][0])
        #
        # for i, dih in enumerate(self.topology.dihedrals):
        #     for j in range(4):
        #         self.assertEqual(self.polymer_linear[dih[j]].ff_map,
        #                          dih[4][j])
        #         self.assertEqual(self.polymer_linear[dih[j]].specie.symbol,
        #                          dih[4][j][0])

        # box_size = [[0, 50],
        #             [0, 50],
        #             [0, 50]]
        self.matrix_config = [{"number": 1, "inside box":[0,0,0,10,10,10]}]

        # constituent molecules
        self.molecules = [self.dmoe]
        #Polymer matrix settings
        output_file = os.path.expanduser("~/dmoe.pdb")
        self.pmr = PackmolRunner(self.molecules,
                            self.matrix_config,
                            tolerance=2.0,
                            filetype="pdb",
                            control_params={"nloop": 1000},
                            output_file=output_file)
        self.packed_mol_with_ff_map = self.pmr.restore_site_properties("ff_map")
        self.box_size = [[0.0, 10],
                    [0.0, 10],
                    [0.0, 10]]
        self.lammps_ff_data = LammpsForceFieldData.from_forcefield_and_topology(
        self.molecules, [1], self.box_size, self.packed_mol_with_ff_map,
        self.forcefield, [self.topology])


        #test pair coeffs
        natoms, natom_types, atomic_masses_dict = \
            LammpsForceFieldData.get_basic_system_info(self.dmoe)

        # atom_types = list(atomic_masses_dict.keys())
        ans_atomtype_ids = [x[0] for x in atomic_masses_dict.values()]

        pair_coeffs = self.lammps_ff_data.pair_coeffs
        atomtype_ids = [x[0] for x in pair_coeffs]

        self.assertEqual(ans_atomtype_ids, atomtype_ids)
    def test_ff_map(self):
        bonds = OrderedDict([((u'C2', u'C2'), [222.5, 1.53]),
                             ((u'C2', u'H2'), [309.0, 1.111]),
                             ((u'C2', u'O'), [360.0, 1.415]),
                             ((u'C3', u'H3'), [322.0, 1.111]),
                             ((u'C3', u'O'), [360.0, 1.415])])
        pairs = OrderedDict([((u'C2', u'C2'), [-0.056, 2.01, -0.01, 1.9]),
                             ((u'C3', u'C3'), [-0.078, 2.04, -0.01, 1.9]),
                             ((u'H2', u'H2'), [-0.035, 1.34]),
                             ((u'H3', u'H3'), [-0.024, 1.34]),
                             ((u'O', u'O'), [-0.1, 1.65])])
        angles = OrderedDict([((u'C2', u'C2', u'H2'), [26.5, 110.1]),
                              ((u'C2', u'C2', u'O'), [45.0, 111.5]),
                              ((u'C2', u'O', u'C2'), [95.0, 109.7]),
                              ((u'C3', u'O', u'C2'), [95.0, 109.7]),
                              ((u'H2', u'C2', u'H2'), [35.5, 109.0]),
                              ((u'H2', u'C2', u'O'), [60.0, 109.5]),
                              ((u'H3', u'C3', u'H3'), [35.5, 108.4]),
                              ((u'H3', u'C3', u'O'), [60.0, 109.5])])
        dihedrals = OrderedDict([((u'C2', u'C2', u'O', u'C2'), [0.57, 1, 0, 0.0]),
                                 ((u'C2', u'O', u'C2', u'H2'), [0.284, 3, 0, 0.0]),
                                 ((u'C2', u'O', u'C3', u'H3'), [0.284, 3, 0, 0.0]),
                                 ((u'C3', u'O', u'C2', u'C2'), [0.57, 1, 0, 0.0]),
                                 ((u'H2', u'C2', u'C2', u'H2'), [0.19, 3, 0, 0.0]),
                                 ((u'H2', u'C2', u'O', u'C3'), [0.284, 3, 0, 0.0]),
                                 ((u'H2', u'O', u'C2', u'H2'), [0.284, 3, 0, 0.0]),
                                 ((u'H3', u'O', u'C3', u'H3'), [0.284, 3, 0, 0.0]),
                                 ((u'O', u'C2', u'C2', u'H2'), [0.19, 3, 0, 0.0]),
                                 ((u'O', u'C2', u'C2', u'O'), [1.16, 2, 0, 0.0])])

        self.assertDictEqual(pairs, self.forcefield.pairs)
        self.assertDictEqual(bonds, self.forcefield.bonds)
        self.assertDictEqual(angles, self.forcefield.angles)
        self.assertDictEqual(dihedrals, self.forcefield.dihedrals)

    def tearDown(self):
        for x in ["lammps_ff_data.dat"]:
            if os.path.exists(os.path.join(test_dir, x)):
                os.remove(os.path.join(test_dir, x))


if __name__ == "__main__":
    unittest.main()
