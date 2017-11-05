# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import unittest
import os
import random

import numpy as np
from ruamel.yaml import YAML
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
        force_field = peptide.force_field
        self.assertEqual(len(force_field["Pair Coeffs"]), 14)
        self.assertEqual(len(force_field["Bond Coeffs"]), 18)
        self.assertEqual(len(force_field["Angle Coeffs"]), 31)
        self.assertEqual(len(force_field["Dihedral Coeffs"]), 21)
        self.assertEqual(len(force_field["Improper Coeffs"]), 2)
        # header box
        np.testing.assert_array_equal(peptide.box_bounds,
                                      [[36.840194, 64.211560],
                                       [41.013691, 68.385058],
                                       [29.768095, 57.139462]])
        # body  last line of each section
        self.assertDictEqual(peptide.masses[-1], {"id": 14, "mass": 1.0100})
        self.assertDictEqual(force_field["Pair Coeffs"][-1],
                             {"id": 14, "coeffs": [0.046000, 0.400014,
                                                   0.046000, 0.400014]})
        self.assertDictEqual(force_field["Bond Coeffs"][-1],
                             {"id": 18, "coeffs": [450.000000, 0.957200]})
        self.assertDictEqual(force_field["Angle Coeffs"][-1],
                             {"id": 31, "coeffs": [55.000000, 104.520000,
                                                   0.000000, 0.000000]})
        self.assertDictEqual(force_field["Dihedral Coeffs"][-1],
                             {"id": 21, "coeffs": [0.010000, 3, 0, 1.000000]})
        for c in force_field["Dihedral Coeffs"][-1]["coeffs"][1:2]:
            self.assertIsInstance(c, int)
        self.assertDictEqual(force_field["Improper Coeffs"][-1],
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
        pairij = virus.force_field["PairIJ Coeffs"]
        self.assertEqual(len(pairij), n * (n + 1) / 2)
        self.assertDictEqual(pairij[-1],
                             {"id1": 4, "id2": 4, "coeffs": [1, 1, 1.1225]})


class TopologyTest(unittest.TestCase):

    def test_init(self):
        inner_charge = np.random.rand(10) - 0.5
        outer_charge = np.random.rand(10) - 0.5
        inner_velo = np.random.rand(10, 3) - 0.5
        outer_velo = np.random.rand(10, 3) - 0.5
        m = Molecule(["H"] * 10, np.random.rand(10, 3) * 100,
                     site_properties={"ff_map": ["D"] * 10,
                                      "charge": inner_charge,
                                      "velocities": inner_velo})
        # q and v from site properties, while type from species_string
        topo = Topology(sites=m)
        self.assertListEqual(topo.types, ["H"] * 10)
        np.testing.assert_array_equal(topo.charges, inner_charge)
        np.testing.assert_array_equal(topo.velocities, inner_velo)
        # q and v from overriding, while type from site property
        topo_override = Topology(sites=m, atom_type="ff_map",
                                 charges=outer_charge,
                                 velocities=outer_velo)
        self.assertListEqual(topo_override.types, ["D"] * 10)
        np.testing.assert_array_equal(topo_override.charges, outer_charge)
        np.testing.assert_array_equal(topo_override.velocities, outer_velo)
        # test using a list of sites instead of SiteCollection
        topo_from_list = Topology(sites=m.sites)
        self.assertListEqual(topo_from_list.types, topo.types)
        np.testing.assert_array_equal(topo_from_list.charges, topo.charges)
        np.testing.assert_array_equal(topo_from_list.velocities,
                                      topo.velocities)

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

    @classmethod
    def setUpClass(cls):
        cls.benzene = ForceField.from_file(os.path.join(test_dir,
                                                        "ff_benzene.yaml"))
        cls.virus = ForceField.from_file(os.path.join(test_dir,
                                                      "ff_virus.yaml"))

    def test_init(self):
        mass_dict = {"H0": 1.00794, "O0": 15.9994}
        mass_dict_alt = {"H0": Element("H"), "O0": "O"}
        ff = ForceField(mass_dict=mass_dict_alt)
        self.assertDictEqual(ff.mass_dict, mass_dict)
        self.assertListEqual(ff.masses, [{"mass": 1.00794, "id": 1},
                                         {"mass": 15.9994, "id": 2}])
        self.assertDictEqual(ff.atom_map, {"H0": 1, "O0": 2})

    def test_get_coeffs_and_mapper(self):
        ff = self.benzene
        masses_data, masses_map = ff.get_coeffs_and_mapper(section="Masses")
        self.assertListEqual(masses_data, [{"id": 1, "mass": 12.01115},
                                           {"id": 2, "mass": 1.00797}])
        self.assertDictEqual(masses_map, {"C": 1, "H": 2})
        bonds_data, \
            bonds_map = ff.get_coeffs_and_mapper(section="Bond Coeffs")
        self.assertDictEqual(bonds_map, {"C-C": 1, "H-C": 2, "C-H": 2})
        self.assertListEqual(bonds_data[1]["coeffs"], [363.4164, 1.08])
        angles_data, \
            angles_map = ff.get_coeffs_and_mapper(section="Angle Coeffs")
        self.assertDictEqual(angles_map, {"C-C-C": 1, "H-C-C": 2, "C-C-H": 2})
        self.assertListEqual(angles_data[1]["coeffs"], [37.0, 120.0])
        dihedrals_data, \
            dihedrals_map = ff.get_coeffs_and_mapper(section="Dihedral "
                                                             "Coeffs")
        self.assertDictEqual(dihedrals_map,
                             {"C-C-C-C": 1, "H-C-C-C": 2, "H-C-C-H": 3,
                              "C-C-C-H": 2})
        self.assertListEqual(dihedrals_data[2]["coeffs"], [3.0, -1, 2])
        impropers_data, \
            impropers_map = ff.get_coeffs_and_mapper(section="Improper "
                                                             "Coeffs")
        self.assertDictEqual(impropers_map, {"H-C-C-C": 1, "C-C-C-H": 1})
        self.assertListEqual(impropers_data[0]["coeffs"], [0.37, -1, 2])

    def test_get_pair_coeffs(self):
        ff = self.virus
        pair = ff.get_pair_coeffs(section="Pair Coeffs")
        self.assertListEqual(pair, [{"id": 1, "coeffs": [1, 1, 1.1225]},
                                    {"id": 2, "coeffs": [1, 1.35, 4]},
                                    {"id": 3, "coeffs": [1, 2.1, 4]},
                                    {"id": 4, "coeffs": [1, 1, 1.1225]}])
        unsorted_pairij = ff.get_pair_coeffs(section="PairIJ Coeffs",
                                             sort_id=False)
        self.assertListEqual(
            unsorted_pairij,
            [{"id1": 1, "id2": 1, "coeffs": [1, 1, 1.1225]},
             {"id1": 1, "id2": 2, "coeffs": [1, 1.175, 1.31894]},
             {"id1": 1, "id2": 4, "coeffs": [1, 1, 1.1225]},
             {"id1": 2, "id2": 2, "coeffs": [1, 1.35, 4]},
             {"id1": 1, "id2": 3, "coeffs": [1, 1.55, 1.73988]},
             {"id1": 2, "id2": 3, "coeffs": [1, 1.725, 1.93631]},
             {"id1": 3, "id2": 3, "coeffs": [1, 2.1, 4]},
             {"id1": 3, "id2": 4, "coeffs": [1, 1.55, 1.73988]},
             {"id1": 2, "id2": 4, "coeffs": [1, 1.175, 1.31894]},
             {"id1": 4, "id2": 4, "coeffs": [1, 1, 1.1225]}]
        )
        sorted_pairij = ff.get_pair_coeffs(section="PairIJ Coeffs",
                                           sort_id=True)
        self.assertListEqual(sorted_pairij,
                             sorted(unsorted_pairij,
                                    key=lambda d: (d["id1"], d["id2"])))

    def test_to_file(self):
        filename = "ff_test.yaml"
        self.benzene.to_file(filename=filename)
        yaml = YAML(typ="safe")
        with open(filename, "r") as f:
            d = yaml.load(f)
        self.assertDictEqual(d["mass_dict"], self.benzene.mass_dict)
        self.assertDictEqual(d["ff_coeffs"], self.benzene.ff_coeffs)

    def test_from_file(self):
        ff = self.virus  # from_file already called in setUpClass
        self.assertDictEqual(ff.mass_dict,
                             {"A": 1.0, "B": 1.728, "C": 2.744, "D": 1.0})
        self.assertDictEqual(ff.ff_coeffs["Pair Coeffs"],
                             {"A": [1, 1, 1.1225], "B": [1, 1.35, 4],
                              "C": [1, 2.1, 4], "D": [1, 1, 1.1225]})
        self.assertDictEqual(ff.ff_coeffs["PairIJ Coeffs"],
                             {"A-A": [1, 1, 1.1225],
                              "A-B": [1, 1.175, 1.31894],
                              "A-D": [1, 1, 1.1225],
                              "B-B": [1, 1.35, 4],
                              "C-A": [1, 1.55, 1.73988],
                              "C-B": [1, 1.725, 1.93631],
                              "C-C": [1, 2.1, 4],
                              "C-D": [1, 1.55, 1.73988],
                              "D-B": [1, 1.175, 1.31894],
                              "D-D": [1, 1, 1.1225]})

    @classmethod
    def tearDownClass(cls):
        if os.path.exists("ff_test.yaml"):
            os.remove("ff_test.yaml")


if __name__ == "__main__":
    unittest.main()
