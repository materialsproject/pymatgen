# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import unittest
import os
import random
import itertools

import numpy as np
from pymatgen import Molecule, Element

from pymatgen.io.lammps.data import LammpsData, Topology, ForceField


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

    @staticmethod
    def test_init():
        nsites = random.randint(1, 10)
        inner_charge = np.random.rand(nsites) - 0.5
        outer_charge = np.random.rand(nsites) - 0.5
        inner_velo = np.random.rand(nsites, 3) - 0.5
        outer_velo = np.random.rand(nsites, 3) - 0.5
        m = Molecule(["H"] * nsites, np.random.rand(nsites, 3) * 100,
                     site_properties={"charge": inner_charge,
                                      "velocities": inner_velo})
        # test set charges and velocities from site properties
        topo = Topology(sites=m)
        np.testing.assert_array_equal(topo.charges, inner_charge)
        np.testing.assert_array_equal(topo.velocities, inner_velo)
        # test using a list of sites instead of SiteCollection
        topo_from_list = Topology(sites=m.sites)
        np.testing.assert_array_equal(topo_from_list.charges, inner_charge)
        np.testing.assert_array_equal(topo_from_list.velocities, inner_velo)
        # test overriding charges and velocities
        topo_override = Topology(sites=m, charges=outer_charge,
                                 velocities=outer_velo)
        np.testing.assert_array_equal(topo_override.charges, outer_charge)
        np.testing.assert_array_equal(topo_override.velocities, outer_velo)

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
        topo_etoh1 = Topology.from_bonding(molecule=etoh, angle=False)
        self.assertIsNone(topo_etoh1.topologies["Angles"])
        topo_etoh2 = Topology.from_bonding(molecule=etoh, dihedral=False)
        self.assertIsNone(topo_etoh2.topologies["Dihedrals"])


class ForceFieldTest(unittest.TestCase):

    def test_init(self):
        masses = {"H0": 1.00794, "O0": 15.9994}
        masses_alt = {"H0": Element("H"), "O0": "O"}
        ff = ForceField(masses=masses_alt)
        self.assertDictEqual(ff.masses, masses)
        self.assertListEqual(ff.masses_data, [{"mass": 1.00794, "id": 1},
                                              {"mass": 15.9994, "id": 2}])
        self.assertDictEqual(ff._atom_map, {"H0": 1, "O0": 2})

    def test_get_coeffs_and_mapper(self):
        masses = {"C": 12, "H": 1}
        bonds = ["C-C", "C-H"]
        angles = ["H-C-H", "C-C-H", "C-C-C"]
        dihedrals = ["H-C-C-H", "C-C-C-H", "C-C-C-C"]
        impropers = ["C-H-H-C"]
        ff_coeffs = {
            "Bond Coeffs":
                {k: v for k, v in zip(bonds, np.random.rand(2, 2))},
            "Angle Coeffs":
                {k: v for k, v in zip(angles, np.random.rand(3, 4))},
            "Dihedral Coeffs":
                {k: v for k, v in zip(dihedrals, np.random.rand(3, 4))},
            "Improper Coeffs":
                {k: v for k, v in zip(impropers, np.random.rand(1, 2))}
        }
        ff = ForceField(masses=masses, ff_coeffs=ff_coeffs)
        masses_data, masses_map = ff.get_coeffs_and_mapper(section="Masses")
        self.assertListEqual(masses_data, [{"id": 1, "mass": 12},
                                           {"id": 2, "mass": 1}])
        self.assertDictEqual(masses_map, {"C": 1, "H": 2})
        bonds_data, \
            bonds_map = ff.get_coeffs_and_mapper(section="Bond Coeffs")
        self.assertDictEqual(bonds_map, {"C-C": 1, "C-H": 2})
        bond, i = random.sample(bonds_map.items(), 1)[0]
        np.testing.assert_array_equal(bonds_data[i - 1]["coeffs"],
                                      ff_coeffs["Bond Coeffs"][bond])
        angles_data, \
            angles_map = ff.get_coeffs_and_mapper(section="Angle Coeffs")
        self.assertDictEqual(angles_map, {"H-C-H": 1, "C-C-H": 2, "C-C-C": 3})
        angle, j = random.sample(angles_map.items(), 1)[0]
        np.testing.assert_array_equal(angles_data[j - 1]["coeffs"],
                                      ff_coeffs["Angle Coeffs"][angle])
        dihedrals_data, \
            dihedrals_map = ff.get_coeffs_and_mapper(section="Dihedral "
                                                             "Coeffs")
        self.assertDictEqual(dihedrals_map,
                             {"H-C-C-H": 1, "C-C-C-H": 2, "C-C-C-C": 3})
        dihedral, k = random.sample(dihedrals_map.items(), 1)[0]
        np.testing.assert_array_equal(dihedrals_data[k - 1]["coeffs"],
                                      ff_coeffs["Dihedral Coeffs"][dihedral])
        impropers_data, \
            impropers_map = ff.get_coeffs_and_mapper(section="Improper "
                                                             "Coeffs")
        self.assertDictEqual(impropers_map, {"C-H-H-C": 1})
        np.testing.assert_array_equal(impropers_data[0]["coeffs"],
                                      ff_coeffs["Improper Coeffs"]["C-H-H-C"])

    def test_get_pair_coeffs(self):
        masses = {"C": 12, "H": 1, "O": 16}
        pair_coeffs = np.random.rand(3, 2)
        pairij_coeffs = np.random.rand(6, 2)
        pair_keys = masses.keys()
        pairij_keys = ["-".join(k) for k in
                       itertools.combinations_with_replacement(masses.keys(),
                                                               2)]
        ff_coeffs = {
            "Pair Coeffs":
                {k: v for k, v in zip(pair_keys, pair_coeffs)},
            "PairIJ Coeffs":
                {k: v for k, v in zip(pairij_keys, pairij_coeffs)}
        }
        ff = ForceField(masses=masses, ff_coeffs=ff_coeffs)
        pair_data = ff.get_pair_coeffs(section="Pair Coeffs")
        p = random.randint(0, 2)
        self.assertEqual(pair_data[p]["id"], p + 1)
        np.testing.assert_array_equal(pair_data[p]["coeffs"], pair_coeffs[p])
        pairij_data = ff.get_pair_coeffs(section="PairIJ Coeffs")
        pij = random.randint(0, 5)
        np.testing.assert_array_equal(pairij_data[pij]["coeffs"],
                                      pairij_coeffs[pij])
        # sort id feature
        pair = list(zip(pair_keys, pair_coeffs))
        random.shuffle(pair)
        pairij = list(zip(pairij_keys, pairij_coeffs))
        random.shuffle(pairij)
        shuffled_coeffs = {"Pair Coeffs": {k: v for k, v in pair},
                           "PairIJ Coeffs": {k: v for k, v in pairij}}
        ff_sort = ForceField(masses=masses, ff_coeffs=shuffled_coeffs)
        pair_data_sorted = ff_sort.get_pair_coeffs(section="Pair Coeffs",
                                                   sort_id=True)
        self.assertListEqual([d["id"] for d in pair_data_sorted], [1, 2, 3])
        np.testing.assert_array_equal([d["coeffs"] for d in pair_data_sorted],
                                      pair_coeffs)
        pairij_data_sorted = ff.get_pair_coeffs(section="PairIJ Coeffs",
                                                sort_id=True)
        np.testing.assert_array_equal([d["coeffs"] for d in
                                       pairij_data_sorted], pairij_coeffs)


        


if __name__ == "__main__":
    unittest.main()
