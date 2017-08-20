# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

import numpy as np

from pymatgen.io.lammps.data import LammpsData, parse_data_file
from pymatgen.core.structure import Molecule
from pymatgen.util.testing import PymatgenTest

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class TestLammpsDataMolecule(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        polymer_chain = Molecule.from_file(os.path.join(test_dir,"polymer_chain.xyz"))
        box_size = [[0.0, 20.0], [0.0, 20.0], [0.0, 20.0]]
        cls.lammps_data = LammpsData.from_structure(polymer_chain, box_size)

    def test_system_info(self):
        atomic_masses = [[1, 1.00794], [2, 12.0107], [3, 15.9994]]
        atoms_data = [[1, 1, 2, 0.0, 10.216511872506619, 11.338023345800135, 12.744427580409154],
                 [2, 1, 1, 0.0, 9.8598518725066189, 12.346833345800135, 12.744427580409154],
                 [3, 1, 1, 0.0, 9.844392872506619, 10.820737345800135, 11.842773580409153],
                 [4, 1, 1, 0.0, 9.844392872506619, 10.820737345800135, 13.646081580409154],
                 [5, 1, 3, 0.0, 11.724011872506619, 11.338023345800135, 12.744427580409154],
                 [6, 1, 2, 0.0, 12.161361872506619, 9.9959933458001355, 12.744427580409154],
                 [7, 1, 1, 0.0, 11.789241872506619, 9.4787033458001346, 11.842773580409153],
                 [8, 1, 1, 0.0, 11.789241872506619, 9.4787033458001346, 13.646081580409154],
                 [9, 1, 2, 0.0, 12.161361872506619, 9.9959933458001355, 11.236927580409153],
                 [10, 1, 1, 0.0, 11.057400872506619, 9.9837103458001355, 11.249412580409153],
                 [11, 1, 1, 0.0, 12.54163787250662, 8.959522345800135, 11.249412580409153],
                 [12, 1, 3, 0.0, 12.647630872506618, 10.700686345800134, 9.9961595804091541],
                 [13, 1, 2, 0.0, 12.161361872506619, 9.9959933458001355, 8.8739875804091533],
                 [14, 1, 1, 0.0, 11.057398872506619, 9.9837063458001349, 8.8864715804091539],
                 [15, 1, 1, 0.0, 12.541635872506619, 8.9595193458001354, 8.8864715804091539],
                 [16, 1, 2, 0.0, 12.161361872506619, 8.4884933458001353, 8.8739875804091533],
                 [17, 1, 1, 0.0, 11.524257872506618, 8.5009783458001351, 7.9723335804091535],
                 [18, 1, 1, 0.0, 11.524257872506618, 8.5009783458001351, 9.7756415804091539],
                 [19, 1, 3, 0.0, 13.017544872506619, 7.2477253458001352, 8.8739875804091533],
                 [20, 1, 2, 0.0, 12.161361872506619, 6.1255533458001352, 8.8739875804091533],
                 [21, 1, 1, 0.0, 11.524253872506618, 6.1380373458001349, 7.9723335804091535],
                 [22, 1, 1, 0.0, 11.524253872506618, 6.1380373458001349, 9.7756415804091539],
                 [23, 1, 2, 0.0, 10.653861872506619, 6.1255533458001352, 8.8739875804091533],
                 [24, 1, 1, 0.0, 10.666346872506619, 6.7626573458001351, 7.9723335804091535],
                 [25, 1, 1, 0.0, 10.666346872506619, 6.7626573458001351, 9.7756415804091539],
                 [26, 1, 3, 0.0, 9.413093872506618, 5.2693693458001345, 8.8739875804091533],
                 [27, 1, 2, 0.0, 8.2909218725066189, 6.1255533458001352, 8.8739875804091533],
                 [28, 1, 1, 0.0, 8.3034058725066195, 6.7626613458001348, 7.9723335804091535],
                 [29, 1, 1, 0.0, 8.3034058725066195, 6.7626613458001348, 9.7756415804091539],
                 [30, 1, 2, 0.0, 8.2909218725066189, 7.6330533458001355, 8.8739875804091533],
                 [31, 1, 1, 0.0, 8.9280258725066179, 7.6205673458001346, 7.9723335804091535],
                 [32, 1, 1, 0.0, 8.9280258725066179, 7.6205673458001346, 9.7756415804091539],
                 [33, 1, 3, 0.0, 7.4347378725066182, 8.8738213458001347, 8.8739875804091533],
                 [34, 1, 2, 0.0, 8.2909218725066189, 9.9959933458001355, 8.8739875804091533],
                 [35, 1, 1, 0.0, 8.9280298725066185, 9.9835093458001349, 7.9723335804091535],
                 [36, 1, 1, 0.0, 8.9280298725066185, 9.9835093458001349, 9.7756415804091539],
                 [37, 1, 2, 0.0, 8.2909218725066189, 11.503493345800134, 8.8739875804091533],
                 [38, 1, 1, 0.0, 8.9280258725066179, 11.491008345800134, 7.9723335804091535],
                 [39, 1, 1, 0.0, 8.9280258725066179, 11.491008345800134, 9.7756415804091539],
                 [40, 1, 3, 0.0, 7.4347378725066182, 12.744261345800135, 8.8739875804091533],
                 [41, 1, 2, 0.0, 8.2909218725066189, 13.866433345800136, 8.8739875804091533],
                 [42, 1, 1, 0.0, 8.9280298725066185, 13.853949345800135, 7.9723335804091535],
                 [43, 1, 1, 0.0, 8.9280298725066185, 13.853949345800135, 9.7756415804091539],
                 [44, 1, 2, 0.0, 8.2909218725066189, 13.866433345800136, 10.381487580409154],
                 [45, 1, 1, 0.0, 8.6711978725066192, 12.829962345800135, 10.369001580409153],
                 [46, 1, 1, 0.0, 7.186960872506619, 13.854150345800136, 10.369001580409153],
                 [47, 1, 3, 0.0, 8.777190872506619, 14.571127345800136, 11.622255580409155],
                 [48, 1, 2, 0.0, 8.2909218725066189, 13.866433345800136, 12.744427580409154],
                 [49, 1, 1, 0.0, 8.6711958725066189, 12.829959345800134, 12.731943580409153],
                 [50, 1, 1, 0.0, 7.1869578725066185, 13.854147345800135, 12.731943580409153],
                 [51, 1, 1, 0.0, 8.7837588725066187, 14.580646345800135, 13.427099580409154]]
        natom_types = 3
        natoms = 51
        np.testing.assert_almost_equal(self.lammps_data.atomic_masses,
                                       atomic_masses, decimal=10)
        np.testing.assert_almost_equal(self.lammps_data.atoms_data, atoms_data,
                                       decimal=6)
        self.assertEqual(self.lammps_data.natom_types, natom_types)
        self.assertEqual(self.lammps_data.natoms, natoms)

    def test_from_file(self):
        self.lammps_data.write_data_file(
            os.path.join(test_dir, "lammps_data.dat"))
        lammps_data = LammpsData.from_file(
            os.path.join(test_dir, "lammps_data.dat"))
        self.assertEqual(str(lammps_data), str(self.lammps_data))

    def tearDown(self):
        for x in ["lammps_data.dat"]:
            if os.path.exists(os.path.join(test_dir, x)):
                os.remove(os.path.join(test_dir, x))


class TestLammpsDataStructure(unittest.TestCase):

    def setUp(self):
        self.structure = PymatgenTest.get_structure("Li2O")
        self.lammps_data = LammpsData.from_structure(self.structure)
        neutral_structure = self.structure.copy()
        neutral_structure.remove_oxidation_states()
        self.lammps_data_neutral = LammpsData.from_structure(neutral_structure,
                                                             set_charge=False)

    def test_system_info(self):
        natom_types = 2
        natoms = 3
        atomic_masses = [[1, 6.941], [2, 15.9994]]
        box_size = [[0.0, 3.259762],
                    [0.0, 2.823037],
                    [0.0, 2.661585]]
        box_tilt = [1.629881, 1.629881, 0.941012]
        atoms_data = [[1, 1, 1.0, -1.1525, -1.1525, -1.1525],
                      [2, 1, 1.0, -3.4575, -3.4575, -3.4575],
                      [3, 2, -2.0, 0.0, 0.0, 0.0]]

        self.assertEqual(self.lammps_data.natom_types, natom_types)
        self.assertEqual(self.lammps_data.natoms, natoms)

        np.testing.assert_almost_equal(self.lammps_data.atomic_masses,
                                       atomic_masses, decimal=6)
        np.testing.assert_almost_equal(self.lammps_data.atoms_data, atoms_data,
                                       decimal=6)
        np.testing.assert_almost_equal(self.lammps_data.box_size, box_size,
                                       decimal=6)
        np.testing.assert_almost_equal(self.lammps_data.box_tilt, box_tilt,
                                       decimal=6)
        neutral_atoms_data = np.delete(atoms_data, 2, axis=1)
        np.testing.assert_almost_equal(self.lammps_data_neutral.atoms_data,
                                       neutral_atoms_data, decimal=6)

    def test_from_file(self):
        self.lammps_data.write_data_file(
            os.path.join(test_dir, "lammps_data.dat"))
        lammps_data = LammpsData.from_file(
            os.path.join(test_dir, "lammps_data.dat"), atom_style="atomic")
        self.assertEqual(str(lammps_data), str(self.lammps_data))

    def tearDown(self):
        for x in ["lammps_data.dat"]:
            if os.path.exists(os.path.join(test_dir, x)):
                os.remove(os.path.join(test_dir, x))


class TestLammpsDataParser(unittest.TestCase):

    def setUp(self):
        self.data = parse_data_file(os.path.join(test_dir, "const_pot.data"))

    def test_keys(self):
        ans = sorted(['natoms', 'nbonds', 'nangles', 'ndihedrals', 'nimpropers',
                      'atom-types', 'bond-types', 'angle-types', 'dihedral-types',
                      'improper-types', 'x', 'y', 'z', 'masses', 'bond-coeffs',
                      'angle-coeffs', 'atoms', 'bonds', 'angles'])
        self.assertEqual(ans, sorted(self.data.keys()))

    def test_values(self):
        self.assertEqual(self.data["natoms"], 432)
        self.assertEqual(self.data["nbonds"], 160)
        self.assertEqual(self.data["nangles"], 80)
        self.assertEqual(self.data["ndihedrals"], 0)
        self.assertEqual(self.data["nimpropers"], 0)
        self.assertEqual(self.data["atom-types"], 4)
        self.assertEqual(self.data["bond-types"], 2)
        self.assertEqual(self.data["angle-types"], 1)
        self.assertEqual(self.data["dihedral-types"], 0)
        self.assertEqual(self.data["improper-types"], 0)
        self.assertEqual(self.data["atom-types"], 4)

        np.testing.assert_almost_equal(self.data["x"], [0, 9.83800], decimal=6)
        np.testing.assert_almost_equal(self.data["y"], [0, 8.52000], decimal=6)
        np.testing.assert_almost_equal(self.data["z"], [-40.20000, 40.20000], decimal=6)
        np.testing.assert_almost_equal(self.data["masses"],
                                       [[1, 12.01070000],
                                        [2, 15.03450000],
                                        [3, 12.01000000],
                                        [4, 14.00670000]], decimal=6)
        np.testing.assert_almost_equal(self.data["bond-coeffs"],
                                       [[1, 380.0, 1.46],
                                        [2, 600.0, 1.157]], decimal=6)
        np.testing.assert_almost_equal(self.data["angle-coeffs"],
                                       [[1, 20.0, 180.0]], decimal=6)

        self.assertEqual(len(self.data["atoms"]), self.data["natoms"])
        self.assertEqual(len(self.data["bonds"]), self.data["nbonds"])
        self.assertEqual(len(self.data["angles"]), self.data["nangles"])

        # test random line from the data block
        np.testing.assert_almost_equal(
            self.data["atoms"][110],
            [111, 37, 1, 0.12900000, 4.869000, 5.574000, -13.992000],
            decimal=10)
        self.assertEqual(self.data["bonds"][76], [77, 2, 115, 117])
        self.assertEqual(self.data["angles"][68], [69, 1, 205, 207, 206])


if __name__ == "__main__":
    unittest.main()
