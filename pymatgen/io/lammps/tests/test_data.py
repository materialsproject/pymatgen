# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import unittest
import os
import random

import numpy as np
from pymatgen import Molecule

from pymatgen.io.lammps.data import LammpsData, Topology


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class LammpsDataTest(unittest.TestCase):

    def test_from_file(self):
        peptide = LammpsData.from_file(filename=os.path.join(test_dir,
                                                             "data.peptide"),
                                       atom_style="full")

        # header stats
        self.assertEqual(len(peptide.atoms), 2004)
        topology = peptide.topology
        self.assertEqual(len(topology["Bonds"]), 1365)
        self.assertEqual(len(topology["Angles"]), 786)
        self.assertEqual(len(topology["Dihedrals"]), 207)
        self.assertEqual(len(topology["Impropers"]), 12)
        self.assertEqual(len(peptide.masses), 14)
        ff_coeffs = peptide.ff_coeffs
        self.assertEqual(len(ff_coeffs["Pair Coeffs"]), 14)
        self.assertEqual(len(ff_coeffs["Bond Coeffs"]), 18)
        self.assertEqual(len(ff_coeffs["Angle Coeffs"]), 31)
        self.assertEqual(len(ff_coeffs["Dihedral Coeffs"]), 21)
        self.assertEqual(len(ff_coeffs["Improper Coeffs"]), 2)
        # header box
        np.testing.assert_array_equal(peptide.box_bounds,
                                      [[36.840194, 64.211560],
                                       [41.013691, 68.385058],
                                       [29.768095, 57.139462]])
        # body  last line of each section
        self.assertDictEqual(peptide.masses[-1], {"id": 14, "mass": 1.0100})
        self.assertDictEqual(ff_coeffs["Pair Coeffs"][-1],
                             {"id": 14, "coeffs": [0.046000, 0.400014,
                                                   0.046000, 0.400014]})
        self.assertDictEqual(ff_coeffs["Bond Coeffs"][-1],
                             {"id": 18, "coeffs": [450.000000, 0.957200]})
        self.assertDictEqual(ff_coeffs["Angle Coeffs"][-1],
                             {"id": 31, "coeffs": [55.000000, 104.520000,
                                                   0.000000, 0.000000]})
        self.assertDictEqual(ff_coeffs["Dihedral Coeffs"][-1],
                             {"id": 21, "coeffs": [0.010000, 3, 0, 1.000000]})
        for c in ff_coeffs["Dihedral Coeffs"][-1]["coeffs"][1:2]:
            self.assertIsInstance(c, int)
        self.assertDictEqual(ff_coeffs["Improper Coeffs"][-1],
                             {"id": 2, "coeffs": [20.000000, 0.000000]})
        self.assertDictEqual(peptide.atoms[-1],
                             {"id": 2004, "molecule-ID": 641, "type": 14,
                              "q": 0.417, "x": 56.55074, "y": 49.75049,
                              "z": 48.61854, "nx": 1, "ny": 1, "nz": 1})
        self.assertDictEqual(peptide.velocities[-1],
                             {"id": 2004, "velocity": [-0.010076,
                                                       -0.005729,
                                                       -0.026032]})
        self.assertDictEqual(topology["Bonds"][-1],
                             {"id": 1365, "type": 18, "bond": [2002, 2003]})
        self.assertDictEqual(topology["Angles"][-1],
                             {"id": 786, "type": 31,
                              "angle": [2003, 2002, 2004]})
        self.assertDictEqual(topology["Dihedrals"][-1],
                             {"id": 207, "type": 3,
                              "dihedral": [64, 79, 80, 82]})
        self.assertDictEqual(topology["Impropers"][-1],
                             {"id": 12, "type": 2,
                              "improper": [79, 64, 80, 81]})
        # box_tilt and another atom_style
        quartz = LammpsData.from_file(filename=os.path.join(test_dir,
                                                            "data.quartz"),
                                      atom_style="atomic")
        np.testing.assert_array_equal(quartz.box_tilt, [-2.456700, 0.0, 0.0])
        self.assertDictEqual(quartz.atoms[-1],
                             {"id": 9, "type": 2, "x": 1.375998,
                              "y": -1.140800, "z": -2.443511})
        # sort_id
        nvt = LammpsData.from_file(filename=os.path.join(test_dir,
                                                         "nvt.data"),
                                   sort_id=True)
        atom_id = random.randint(1, 648)
        self.assertEqual(nvt.atoms[atom_id - 1]["id"], atom_id)
        # PairIJ Coeffs section
        virus = LammpsData.from_file(filename=os.path.join(test_dir,
                                                           "virus.data"),
                                     atom_style="angle")
        n = len(virus.masses)
        pairij = virus.ff_coeffs["PairIJ Coeffs"]
        self.assertEqual(len(pairij), n * (n + 1) / 2)
        self.assertDictEqual(pairij[-1],
                             {"id1": 4, "id2": 4, "coeffs": [1, 1, 1.1225]})


class TopologyTest(unittest.TestCase):

    def test_init(self):
        nsites = random.randint(1, 10)
        inner_velo = np.random.rand(nsites, 3)
        outer_velo = np.random.rand(nsites, 3)
        m = Molecule(["H"] * nsites, np.random.rand(nsites, 3) * 100,
                     site_properties={"velocities": inner_velo})
        # test set velocities from site properties
        topo0 = Topology(sites=m, velocities=None)
        np.testing.assert_array_equal(topo0.velocities, inner_velo)
        # test using a list of sites instead of SiteCollection
        topo1 = Topology(sites=m.sites, velocities=None)
        np.testing.assert_array_equal(topo1.velocities, inner_velo)
        # test overriding velocities
        topo2 = Topology(sites=m, velocities=outer_velo)
        np.testing.assert_array_equal(topo2.velocities, outer_velo)
        # test wrong velocity format
        wrong_nsites = np.random.rand(11, 3)
        self.assertRaises(Exception,
                          lambda: Topology(sites=m, velocities=wrong_nsites))
        wrong_dims = np.random.rand(nsites, 2)
        self.assertRaises(Exception,
                          lambda: Topology(sites=m, velocities=wrong_dims))

    def test_from_bonding(self):
        # He: no bonding topologies
        helium = Molecule(["He"], [[0, 0, 0]])
        topo_he = Topology.from_bonding(molecule=helium)
        self.assertIsNone(topo_he.topologies)
        # H2: 1 bond only
        hydrogen = Molecule(["H"] * 2, [[0, 0, 0], [0, 0, 0.7414]])
        topo_h = Topology.from_bonding(molecule=hydrogen)
        tp_h = topo_h.topologies
        self.assertListEqual(tp_h["Bonds"], [[0, 1]])
        self.assertIsNone(tp_h["Angles"])
        self.assertIsNone(tp_h["Dihedrals"])
        # water: 2 bonds and 1 angle only
        water = Molecule(["O", "H", "H"], [[0.0000, 0.0000, 0.1173],
                                           [0.0000, 0.7572, -0.4692],
                                           [0.0000, -0.7572, -0.4692]])
        topo_water = Topology.from_bonding(molecule=water)
        tp_water = topo_water.topologies
        self.assertListEqual(tp_water["Bonds"], [[0, 1], [0, 2]])
        self.assertListEqual(tp_water["Angles"], [[1, 0, 2]])
        self.assertIsNone(tp_water["Dihedrals"])
        # EtOH
        etoh = Molecule(["C", "C", "O", "H", "H", "H", "H", "H", "H"],
                        [[1.1879, -0.3829, 0.0000],
                         [0.0000, 0.5526, 0.0000],
                         [-1.1867, -0.2472, 0.0000],
                         [-1.9237, 0.3850, 0.0000],
                         [2.0985, 0.2306, 0.0000],
                         [1.1184, -1.0093, 0.8869],
                         [1.1184, -1.0093, -0.8869],
                         [-0.0227, 1.1812, 0.8852],
                         [-0.0227, 1.1812, -0.8852]])
        topo_etoh = Topology.from_bonding(molecule=etoh)
        tp_etoh = topo_etoh.topologies
        self.assertEqual(len(tp_etoh["Bonds"]), 8)
        etoh_bonds = [[0, 1], [0, 4], [0, 5], [0, 6],
                      [1, 2], [1, 7], [1, 8], [2, 3]]
        np.testing.assert_array_equal(tp_etoh["Bonds"], etoh_bonds)
        self.assertEqual(len(tp_etoh["Angles"]), 13)
        etoh_angles = [[1, 0, 4], [1, 0, 5], [1, 0, 6], [4, 0, 5], [4, 0, 6],
                       [5, 0, 6], [0, 1, 2], [0, 1, 7], [0, 1, 8], [2, 1, 7],
                       [2, 1, 8], [7, 1, 8], [1, 2, 3]]
        np.testing.assert_array_equal(tp_etoh["Angles"], etoh_angles)
        self.assertEqual(len(tp_etoh["Dihedrals"]), 12)
        etoh_dihedrals = [[4, 0, 1, 2], [4, 0, 1, 7], [4, 0, 1, 8],
                          [5, 0, 1, 2], [5, 0, 1, 7], [5, 0, 1, 8],
                          [6, 0, 1, 2], [6, 0, 1, 7], [6, 0, 1, 8],
                          [0, 1, 2, 3], [7, 1, 2, 3], [8, 1, 2, 3]]
        np.testing.assert_array_equal(tp_etoh["Dihedrals"], etoh_dihedrals)
        # bond flag to off
        topo_etoh0 = Topology.from_bonding(molecule=etoh, bond=False,
                                           angle=True, dihedral=True)
        self.assertIsNone(topo_etoh0.topologies)
        # angle or dihedral flag to off
        topo_etoh1 = Topology.from_bonding(molecule=etoh, angle=False,
                                           dihedral=True)
        self.assertIsNone(topo_etoh1.topologies["Angles"])
        topo_etoh2 = Topology.from_bonding(molecule=etoh, angle=True,
                                           dihedral=False)
        self.assertIsNone(topo_etoh2.topologies["Dihedrals"])


if __name__ == "__main__":
    unittest.main()
