# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import json
import os
import tempfile
import unittest
import warnings

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.io.lobster import (
    Bandoverlaps,
    Charge,
    Cohpcar,
    Doscar,
    Fatband,
    Grosspop,
    Icohplist,
    Lobsterin,
    Lobsterout,
    Wavefunction,
)
from pymatgen.io.lobster.inputs import get_all_possible_basis_combinations
from pymatgen.io.vasp import Vasprun
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.util.testing import PymatgenTest

__author__ = "Janine George, Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.2"
__email__ = "janine.george@uclouvain.be, esters@uoregon.edu"
__date__ = "Dec 10, 2017"

test_dir_doscar = PymatgenTest.TEST_FILES_DIR

this_dir = os.path.dirname(os.path.abspath(__file__))


class CohpcarTest(PymatgenTest):
    def setUp(self):
        self.cohp_bise = Cohpcar(filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "COHPCAR.lobster.BiSe"))
        self.coop_bise = Cohpcar(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "COOPCAR.lobster.BiSe"),
            are_coops=True,
        )
        self.cohp_fe = Cohpcar(filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "COOPCAR.lobster"))
        self.coop_fe = Cohpcar(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "COOPCAR.lobster"),
            are_coops=True,
        )
        self.orb = Cohpcar(filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "COHPCAR.lobster.orbitalwise"))
        self.orb_notot = Cohpcar(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "COHPCAR.lobster.notot.orbitalwise")
        )

        # Lobster 3.1 (Test data is from prerelease of Lobster 3.1)
        self.cohp_KF = Cohpcar(filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "COHPCAR.lobster.KF"))
        self.coop_KF = Cohpcar(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "COHPCAR.lobster.KF"),
            are_coops=True,
        )

        # example with f electrons
        self.cohp_Na2UO4 = Cohpcar(filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "COHPCAR.lobster.Na2UO4"))
        self.coop_Na2UO4 = Cohpcar(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "COOPCAR.lobster.Na2UO4"),
            are_coops=True,
        )

    def test_attributes(self):
        self.assertFalse(self.cohp_bise.are_coops)
        self.assertTrue(self.coop_bise.are_coops)
        self.assertFalse(self.cohp_bise.is_spin_polarized)
        self.assertFalse(self.coop_bise.is_spin_polarized)
        self.assertFalse(self.cohp_fe.are_coops)
        self.assertTrue(self.coop_fe.are_coops)
        self.assertTrue(self.cohp_fe.is_spin_polarized)
        self.assertTrue(self.coop_fe.is_spin_polarized)
        self.assertEqual(len(self.cohp_bise.energies), 241)
        self.assertEqual(len(self.coop_bise.energies), 241)
        self.assertEqual(len(self.cohp_fe.energies), 301)
        self.assertEqual(len(self.coop_fe.energies), 301)
        self.assertEqual(len(self.cohp_bise.cohp_data), 12)
        self.assertEqual(len(self.coop_bise.cohp_data), 12)
        self.assertEqual(len(self.cohp_fe.cohp_data), 3)
        self.assertEqual(len(self.coop_fe.cohp_data), 3)

        # Lobster 3.1
        self.assertFalse(self.cohp_KF.are_coops)
        self.assertTrue(self.coop_KF.are_coops)
        self.assertFalse(self.cohp_KF.is_spin_polarized)
        self.assertFalse(self.coop_KF.is_spin_polarized)
        self.assertEqual(len(self.cohp_KF.energies), 6)
        self.assertEqual(len(self.coop_KF.energies), 6)
        self.assertEqual(len(self.cohp_KF.cohp_data), 7)
        self.assertEqual(len(self.coop_KF.cohp_data), 7)

    def test_energies(self):
        efermi_bise = 5.90043
        elim_bise = (-0.124679, 11.9255)
        efermi_fe = 9.75576
        elim_fe = (-0.277681, 14.7725)
        efermi_KF = -2.87475
        elim_KF = (-11.25000 + efermi_KF, 7.5000 + efermi_KF)

        self.assertEqual(self.cohp_bise.efermi, efermi_bise)
        self.assertEqual(self.coop_bise.efermi, efermi_bise)
        self.assertEqual(self.cohp_fe.efermi, efermi_fe)
        self.assertEqual(self.coop_fe.efermi, efermi_fe)
        # Lobster 3.1
        self.assertEqual(self.cohp_KF.efermi, efermi_KF)
        self.assertEqual(self.coop_KF.efermi, efermi_KF)

        self.assertAlmostEqual(self.cohp_bise.energies[0] + self.cohp_bise.efermi, elim_bise[0], places=4)
        self.assertAlmostEqual(self.cohp_bise.energies[-1] + self.cohp_bise.efermi, elim_bise[1], places=4)
        self.assertAlmostEqual(self.coop_bise.energies[0] + self.coop_bise.efermi, elim_bise[0], places=4)
        self.assertAlmostEqual(self.coop_bise.energies[-1] + self.coop_bise.efermi, elim_bise[1], places=4)

        self.assertAlmostEqual(self.cohp_fe.energies[0] + self.cohp_fe.efermi, elim_fe[0], places=4)
        self.assertAlmostEqual(self.cohp_fe.energies[-1] + self.cohp_fe.efermi, elim_fe[1], places=4)
        self.assertAlmostEqual(self.coop_fe.energies[0] + self.coop_fe.efermi, elim_fe[0], places=4)
        self.assertAlmostEqual(self.coop_fe.energies[-1] + self.coop_fe.efermi, elim_fe[1], places=4)

        # Lobster 3.1
        self.assertAlmostEqual(self.cohp_KF.energies[0] + self.cohp_KF.efermi, elim_KF[0], places=4)
        self.assertAlmostEqual(self.cohp_KF.energies[-1] + self.cohp_KF.efermi, elim_KF[1], places=4)
        self.assertAlmostEqual(self.coop_KF.energies[0] + self.coop_KF.efermi, elim_KF[0], places=4)
        self.assertAlmostEqual(self.coop_KF.energies[-1] + self.coop_KF.efermi, elim_KF[1], places=4)

    def test_cohp_data(self):
        lengths_sites_bise = {
            "1": (2.882308829886294, (0, 6)),
            "2": (3.1014396233274444, (0, 9)),
            "3": (2.8823088298862083, (1, 7)),
            "4": (3.1014396233275434, (1, 8)),
            "5": (3.0500070394403904, (2, 9)),
            "6": (2.9167594580335807, (2, 10)),
            "7": (3.05000703944039, (3, 8)),
            "8": (2.9167594580335803, (3, 11)),
            "9": (3.3752173204052101, (4, 11)),
            "10": (3.0729354518345948, (4, 5)),
            "11": (3.3752173204052101, (5, 10)),
        }
        lengths_sites_fe = {
            "1": (2.8318907764979082, (7, 6)),
            "2": (2.4524893531900283, (7, 8)),
        }
        # Lobster 3.1
        lengths_sites_KF = {
            "1": (2.7119923200622269, (0, 1)),
            "2": (2.7119923200622269, (0, 1)),
            "3": (2.7119923576010501, (0, 1)),
            "4": (2.7119923576010501, (0, 1)),
            "5": (2.7119923200622269, (0, 1)),
            "6": (2.7119923200622269, (0, 1)),
        }

        for data in [self.cohp_bise.cohp_data, self.coop_bise.cohp_data]:
            for bond in data:
                if bond != "average":
                    self.assertEqual(data[bond]["length"], lengths_sites_bise[bond][0])
                    self.assertEqual(data[bond]["sites"], lengths_sites_bise[bond][1])
                    self.assertEqual(len(data[bond]["COHP"][Spin.up]), 241)
                    self.assertEqual(len(data[bond]["ICOHP"][Spin.up]), 241)
        for data in [self.cohp_fe.cohp_data, self.coop_fe.cohp_data]:
            for bond in data:
                if bond != "average":
                    self.assertEqual(data[bond]["length"], lengths_sites_fe[bond][0])
                    self.assertEqual(data[bond]["sites"], lengths_sites_fe[bond][1])
                    self.assertEqual(len(data[bond]["COHP"][Spin.up]), 301)
                    self.assertEqual(len(data[bond]["ICOHP"][Spin.up]), 301)

        # Lobster 3.1
        for data in [self.cohp_KF.cohp_data, self.coop_KF.cohp_data]:
            for bond in data:
                if bond != "average":
                    self.assertEqual(data[bond]["length"], lengths_sites_KF[bond][0])
                    self.assertEqual(data[bond]["sites"], lengths_sites_KF[bond][1])
                    self.assertEqual(len(data[bond]["COHP"][Spin.up]), 6)
                    self.assertEqual(len(data[bond]["ICOHP"][Spin.up]), 6)

    def test_orbital_resolved_cohp(self):
        orbitals = [tuple((Orbital(i), Orbital(j))) for j in range(4) for i in range(4)]
        self.assertIsNone(self.cohp_bise.orb_res_cohp)
        self.assertIsNone(self.coop_bise.orb_res_cohp)
        self.assertIsNone(self.cohp_fe.orb_res_cohp)
        self.assertIsNone(self.coop_fe.orb_res_cohp)
        self.assertIsNone(self.orb_notot.cohp_data["1"]["COHP"])
        self.assertIsNone(self.orb_notot.cohp_data["1"]["ICOHP"])
        for orbs in self.orb.orb_res_cohp["1"]:
            orb_set = self.orb.orb_res_cohp["1"][orbs]["orbitals"]
            self.assertEqual(orb_set[0][0], 4)
            self.assertEqual(orb_set[1][0], 4)
            self.assertIn(tuple((orb_set[0][1], orb_set[1][1])), orbitals)

        # test d and f orbitals
        comparelist = [
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            5,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            7,
            7,
            7,
            7,
        ]
        comparelist2 = [
            "f0",
            "f0",
            "f0",
            "f0",
            "f1",
            "f1",
            "f1",
            "f1",
            "f2",
            "f2",
            "f2",
            "f2",
            "f3",
            "f3",
            "f3",
            "f3",
            "f_1",
            "f_1",
            "f_1",
            "f_1",
            "f_2",
            "f_2",
            "f_2",
            "f_2",
            "f_3",
            "f_3",
            "f_3",
            "f_3",
            "dx2",
            "dx2",
            "dx2",
            "dx2",
            "dxy",
            "dxy",
            "dxy",
            "dxy",
            "dxz",
            "dxz",
            "dxz",
            "dxz",
            "dyz",
            "dyz",
            "dyz",
            "dyz",
            "dz2",
            "dz2",
            "dz2",
            "dz2",
            "px",
            "px",
            "px",
            "px",
            "py",
            "py",
            "py",
            "py",
            "pz",
            "pz",
            "pz",
            "pz",
            "s",
            "s",
            "s",
            "s",
            "s",
            "s",
            "s",
            "s",
        ]
        for iorb, orbs in enumerate(sorted(self.cohp_Na2UO4.orb_res_cohp["49"])):
            orb_set = self.cohp_Na2UO4.orb_res_cohp["49"][orbs]["orbitals"]
            self.assertEqual(orb_set[0][0], comparelist[iorb])
            self.assertEqual(str(orb_set[0][1]), comparelist2[iorb])

        # The sum of the orbital-resolved COHPs should be approximately
        # the total COHP. Due to small deviations in the LOBSTER calculation,
        # the precision is not very high though.
        cohp = self.orb.cohp_data["1"]["COHP"][Spin.up]
        icohp = self.orb.cohp_data["1"]["ICOHP"][Spin.up]
        tot = np.sum(
            [self.orb.orb_res_cohp["1"][orbs]["COHP"][Spin.up] for orbs in self.orb.orb_res_cohp["1"]],
            axis=0,
        )
        self.assertArrayAlmostEqual(tot, cohp, decimal=3)
        tot = np.sum(
            [self.orb.orb_res_cohp["1"][orbs]["ICOHP"][Spin.up] for orbs in self.orb.orb_res_cohp["1"]],
            axis=0,
        )
        self.assertArrayAlmostEqual(tot, icohp, decimal=3)

        # Lobster 3.1
        cohp_KF = self.cohp_KF.cohp_data["1"]["COHP"][Spin.up]
        icohp_KF = self.cohp_KF.cohp_data["1"]["ICOHP"][Spin.up]
        tot_KF = np.sum(
            [self.cohp_KF.orb_res_cohp["1"][orbs]["COHP"][Spin.up] for orbs in self.cohp_KF.orb_res_cohp["1"]],
            axis=0,
        )
        self.assertArrayAlmostEqual(tot_KF, cohp_KF, decimal=3)
        tot_KF = np.sum(
            [self.cohp_KF.orb_res_cohp["1"][orbs]["ICOHP"][Spin.up] for orbs in self.cohp_KF.orb_res_cohp["1"]],
            axis=0,
        )
        self.assertArrayAlmostEqual(tot_KF, icohp_KF, decimal=3)

        # d and f orbitals
        cohp_Na2UO4 = self.cohp_Na2UO4.cohp_data["49"]["COHP"][Spin.up]
        icohp_Na2UO4 = self.cohp_Na2UO4.cohp_data["49"]["ICOHP"][Spin.up]
        tot_Na2UO4 = np.sum(
            [
                self.cohp_Na2UO4.orb_res_cohp["49"][orbs]["COHP"][Spin.up]
                for orbs in self.cohp_Na2UO4.orb_res_cohp["49"]
            ],
            axis=0,
        )
        self.assertArrayAlmostEqual(tot_Na2UO4, cohp_Na2UO4, decimal=3)
        tot_Na2UO4 = np.sum(
            [
                self.cohp_Na2UO4.orb_res_cohp["49"][orbs]["ICOHP"][Spin.up]
                for orbs in self.cohp_Na2UO4.orb_res_cohp["49"]
            ],
            axis=0,
        )
        self.assertArrayAlmostEqual(tot_Na2UO4, icohp_Na2UO4, decimal=3)


class IcohplistTest(unittest.TestCase):
    def setUp(self):
        self.icohp_bise = Icohplist(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "ICOHPLIST.lobster.BiSe")
        )
        self.icoop_bise = Icohplist(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "ICOOPLIST.lobster.BiSe"),
            are_coops=True,
        )
        self.icohp_fe = Icohplist(filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "ICOHPLIST.lobster"))
        # allow gzipped files
        self.icohp_gzipped = Icohplist(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "ICOHPLIST.lobster.gz")
        )
        self.icoop_fe = Icohplist(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "ICOHPLIST.lobster"),
            are_coops=True,
        )

    def test_attributes(self):
        self.assertFalse(self.icohp_bise.are_coops)
        self.assertTrue(self.icoop_bise.are_coops)
        self.assertFalse(self.icohp_bise.is_spin_polarized)
        self.assertFalse(self.icoop_bise.is_spin_polarized)
        self.assertEqual(len(self.icohp_bise.icohplist), 11)
        self.assertEqual(len(self.icoop_bise.icohplist), 11)
        self.assertFalse(self.icohp_fe.are_coops)
        self.assertTrue(self.icoop_fe.are_coops)
        self.assertTrue(self.icohp_fe.is_spin_polarized)
        self.assertTrue(self.icoop_fe.is_spin_polarized)
        self.assertEqual(len(self.icohp_fe.icohplist), 2)
        self.assertEqual(len(self.icoop_fe.icohplist), 2)

    def test_values(self):
        icohplist_bise = {
            "1": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -2.18042},
                "translation": [0, 0, 0],
            },
            "2": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.14347},
                "translation": [0, 0, 0],
            },
            "3": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -2.18042},
                "translation": [0, 0, 0],
            },
            "4": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.14348},
                "translation": [0, 0, 0],
            },
            "5": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.30006},
                "translation": [0, 0, 0],
            },
            "6": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.96843},
                "translation": [0, 0, 0],
            },
            "7": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.30006},
                "translation": [0, 0, 0],
            },
            "8": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -1.96843},
                "translation": [0, 0, 0],
            },
            "9": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.47531},
                "translation": [0, 0, 0],
            },
            "10": {
                "length": 3.07294,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -2.38796},
                "translation": [0, 0, 0],
            },
            "11": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.47531},
                "translation": [0, 0, 0],
            },
        }
        icooplist_bise = {
            "1": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.14245},
                "translation": [0, 0, 0],
            },
            "2": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.04118},
                "translation": [0, 0, 0],
            },
            "3": {
                "length": 2.88231,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.14245},
                "translation": [0, 0, 0],
            },
            "4": {
                "length": 3.10144,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.04118},
                "translation": [0, 0, 0],
            },
            "5": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.03516},
                "translation": [0, 0, 0],
            },
            "6": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.10745},
                "translation": [0, 0, 0],
            },
            "7": {
                "length": 3.05001,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.03516},
                "translation": [0, 0, 0],
            },
            "8": {
                "length": 2.91676,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.10745},
                "translation": [0, 0, 0],
            },
            "9": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.12395},
                "translation": [0, 0, 0],
            },
            "10": {
                "length": 3.07294,
                "number_of_bonds": 3,
                "icohp": {Spin.up: 0.24714},
                "translation": [0, 0, 0],
            },
            "11": {
                "length": 3.37522,
                "number_of_bonds": 3,
                "icohp": {Spin.up: -0.12395},
                "translation": [0, 0, 0],
            },
        }
        icooplist_fe = {
            "1": {
                "length": 2.83189,
                "number_of_bonds": 2,
                "icohp": {Spin.up: -0.10218, Spin.down: -0.19701},
                "translation": [0, 0, 0],
            },
            "2": {
                "length": 2.45249,
                "number_of_bonds": 1,
                "icohp": {Spin.up: -0.28485, Spin.down: -0.58279},
                "translation": [0, 0, 0],
            },
        }

        self.assertEqual(icohplist_bise, self.icohp_bise.icohplist)
        self.assertEqual(icooplist_fe, self.icoop_fe.icohplist)


class DoscarTest(unittest.TestCase):
    def setUp(self):
        # first for spin polarized version
        doscar = os.path.join(test_dir_doscar, "DOSCAR.lobster.spin")
        poscar = os.path.join(test_dir_doscar, "POSCAR.lobster.spin_DOS")
        # not spin polarized
        doscar2 = os.path.join(test_dir_doscar, "DOSCAR.lobster.nonspin")
        poscar2 = os.path.join(test_dir_doscar, "POSCAR.lobster.nonspin_DOS")
        self.DOSCAR_spin_pol = Doscar(doscar=doscar, structure_file=poscar)
        self.DOSCAR_nonspin_pol = Doscar(doscar=doscar2, structure_file=poscar2)

        with open(os.path.join(test_dir_doscar, "structure_KF.json"), "r") as f:
            data = json.load(f)

        self.structure = Structure.from_dict(data)

    def test_completedos(self):
        # first for spin polarized version
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_up = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02577]
        tdos_down = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02586]
        fermi = 0.0

        PDOS_F_2s_up = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        PDOS_F_2s_down = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        PDOS_F_2py_up = [0.00000, 0.00160, 0.00000, 0.25801, 0.00000, 0.00029]
        PDOS_F_2py_down = [0.00000, 0.00161, 0.00000, 0.25819, 0.00000, 0.00029]
        PDOS_F_2pz_up = [0.00000, 0.00161, 0.00000, 0.25823, 0.00000, 0.00029]
        PDOS_F_2pz_down = [0.00000, 0.00160, 0.00000, 0.25795, 0.00000, 0.00029]
        PDOS_F_2px_up = [0.00000, 0.00160, 0.00000, 0.25805, 0.00000, 0.00029]
        PDOS_F_2px_down = [0.00000, 0.00161, 0.00000, 0.25814, 0.00000, 0.00029]

        self.assertListEqual(energies_spin, self.DOSCAR_spin_pol.completedos.energies.tolist())
        self.assertListEqual(tdos_up, self.DOSCAR_spin_pol.completedos.densities[Spin.up].tolist())
        self.assertListEqual(tdos_down, self.DOSCAR_spin_pol.completedos.densities[Spin.down].tolist())
        self.assertAlmostEqual(fermi, self.DOSCAR_spin_pol.completedos.efermi)
        for coords, coords2 in zip(
            self.DOSCAR_spin_pol.completedos.structure.frac_coords,
            self.structure.frac_coords,
        ):
            for xyz, xyz2 in zip(coords, coords2):
                self.assertAlmostEqual(xyz, xyz2)
        self.assertListEqual(
            self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2s"][Spin.up].tolist(),
            PDOS_F_2s_up,
        )
        self.assertListEqual(
            self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2s"][Spin.down].tolist(),
            PDOS_F_2s_down,
        )
        self.assertListEqual(
            self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_y"][Spin.up].tolist(),
            PDOS_F_2py_up,
        )
        self.assertListEqual(
            self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_y"][Spin.down].tolist(),
            PDOS_F_2py_down,
        )
        self.assertListEqual(
            self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_z"][Spin.up].tolist(),
            PDOS_F_2pz_up,
        )
        self.assertListEqual(
            self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_z"][Spin.down].tolist(),
            PDOS_F_2pz_down,
        )
        self.assertListEqual(
            self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_x"][Spin.up].tolist(),
            PDOS_F_2px_up,
        )
        self.assertListEqual(
            self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]["2p_x"][Spin.down].tolist(),
            PDOS_F_2px_down,
        )

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]
        PDOS_F_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        PDOS_F_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        PDOS_F_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        PDOS_F_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        self.assertListEqual(energies_nonspin, self.DOSCAR_nonspin_pol.completedos.energies.tolist())

        self.assertListEqual(
            tdos_nonspin,
            self.DOSCAR_nonspin_pol.completedos.densities[Spin.up].tolist(),
        )

        self.assertAlmostEqual(fermi, self.DOSCAR_nonspin_pol.completedos.efermi)
        self.assertDictEqual(
            self.DOSCAR_nonspin_pol.completedos.structure.as_dict(),
            self.structure.as_dict(),
        )

        self.assertListEqual(
            self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2s"][Spin.up].tolist(),
            PDOS_F_2s,
        )
        self.assertListEqual(
            self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2p_y"][Spin.up].tolist(),
            PDOS_F_2py,
        )
        self.assertListEqual(
            self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2p_z"][Spin.up].tolist(),
            PDOS_F_2pz,
        )
        self.assertListEqual(
            self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]["2p_x"][Spin.up].tolist(),
            PDOS_F_2px,
        )

    def test_pdos(self):
        # first for spin polarized version

        PDOS_F_2s_up = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        PDOS_F_2s_down = [0.00000, 0.00159, 0.00000, 0.00011, 0.00000, 0.00069]
        PDOS_F_2py_up = [0.00000, 0.00160, 0.00000, 0.25801, 0.00000, 0.00029]
        PDOS_F_2py_down = [0.00000, 0.00161, 0.00000, 0.25819, 0.00000, 0.00029]
        PDOS_F_2pz_up = [0.00000, 0.00161, 0.00000, 0.25823, 0.00000, 0.00029]
        PDOS_F_2pz_down = [0.00000, 0.00160, 0.00000, 0.25795, 0.00000, 0.00029]
        PDOS_F_2px_up = [0.00000, 0.00160, 0.00000, 0.25805, 0.00000, 0.00029]
        PDOS_F_2px_down = [0.00000, 0.00161, 0.00000, 0.25814, 0.00000, 0.00029]

        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]["2s"][Spin.up].tolist(), PDOS_F_2s_up)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]["2s"][Spin.down].tolist(), PDOS_F_2s_down)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]["2p_y"][Spin.up].tolist(), PDOS_F_2py_up)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]["2p_y"][Spin.down].tolist(), PDOS_F_2py_down)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]["2p_z"][Spin.up].tolist(), PDOS_F_2pz_up)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]["2p_z"][Spin.down].tolist(), PDOS_F_2pz_down)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]["2p_x"][Spin.up].tolist(), PDOS_F_2px_up)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]["2p_x"][Spin.down].tolist(), PDOS_F_2px_down)

        # non spin
        PDOS_F_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        PDOS_F_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        PDOS_F_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        PDOS_F_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        self.assertListEqual(self.DOSCAR_nonspin_pol.pdos[0]["2s"][Spin.up].tolist(), PDOS_F_2s)
        self.assertListEqual(self.DOSCAR_nonspin_pol.pdos[0]["2p_y"][Spin.up].tolist(), PDOS_F_2py)
        self.assertListEqual(self.DOSCAR_nonspin_pol.pdos[0]["2p_z"][Spin.up].tolist(), PDOS_F_2pz)
        self.assertListEqual(self.DOSCAR_nonspin_pol.pdos[0]["2p_x"][Spin.up].tolist(), PDOS_F_2px)

    def test_tdos(self):
        # first for spin polarized version
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_up = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02577]
        tdos_down = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02586]
        fermi = 0.0

        self.assertListEqual(energies_spin, self.DOSCAR_spin_pol.tdos.energies.tolist())
        self.assertListEqual(tdos_up, self.DOSCAR_spin_pol.tdos.densities[Spin.up].tolist())
        self.assertListEqual(tdos_down, self.DOSCAR_spin_pol.tdos.densities[Spin.down].tolist())
        self.assertAlmostEqual(fermi, self.DOSCAR_spin_pol.tdos.efermi)

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]
        fermi = 0.0

        self.assertListEqual(energies_nonspin, self.DOSCAR_nonspin_pol.tdos.energies.tolist())
        self.assertListEqual(tdos_nonspin, self.DOSCAR_nonspin_pol.tdos.densities[Spin.up].tolist())
        self.assertAlmostEqual(fermi, self.DOSCAR_nonspin_pol.tdos.efermi)

    def test_energies(self):
        # first for spin polarized version
        energies_spin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]

        self.assertListEqual(energies_spin, self.DOSCAR_spin_pol.energies.tolist())

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        self.assertListEqual(energies_nonspin, self.DOSCAR_nonspin_pol.energies.tolist())

    def test_tdensities(self):
        # first for spin polarized version
        tdos_up = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02577]
        tdos_down = [0.00000, 0.79999, 0.00000, 0.79999, 0.00000, 0.02586]

        self.assertListEqual(tdos_up, self.DOSCAR_spin_pol.tdensities[Spin.up].tolist())
        self.assertListEqual(tdos_down, self.DOSCAR_spin_pol.tdensities[Spin.down].tolist())

        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]
        self.assertListEqual(tdos_nonspin, self.DOSCAR_nonspin_pol.tdensities[Spin.up].tolist())

    def test_itdensities(self):
        itdos_up = [1.99997, 4.99992, 4.99992, 7.99987, 7.99987, 8.09650]
        itdos_down = [1.99997, 4.99992, 4.99992, 7.99987, 7.99987, 8.09685]
        self.assertListEqual(itdos_up, self.DOSCAR_spin_pol.itdensities[Spin.up].tolist())
        self.assertListEqual(itdos_down, self.DOSCAR_spin_pol.itdensities[Spin.down].tolist())

        itdos_nonspin = [4.00000, 10.00000, 10.00000, 16.00000, 16.00000, 16.09067]
        self.assertListEqual(itdos_nonspin, self.DOSCAR_nonspin_pol.itdensities[Spin.up].tolist())

    def test_is_spin_polarized(self):
        # first for spin polarized version
        self.assertTrue(self.DOSCAR_spin_pol.is_spin_polarized)

        self.assertFalse(self.DOSCAR_nonspin_pol.is_spin_polarized)


class ChargeTest(PymatgenTest):
    def setUp(self):
        self.charge2 = Charge(filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "CHARGE.lobster.MnO"))
        # gzipped file
        self.charge = Charge(filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "CHARGE.lobster.MnO2.gz"))

    def testattributes(self):
        charge_Loewdin = [-1.25, 1.25]
        charge_Mulliken = [-1.30, 1.30]
        atomlist = ["O1", "Mn2"]
        types = ["O", "Mn"]
        num_atoms = 2
        self.assertArrayEqual(charge_Mulliken, self.charge2.Mulliken)
        self.assertArrayEqual(charge_Loewdin, self.charge2.Loewdin)
        self.assertArrayEqual(atomlist, self.charge2.atomlist)
        self.assertArrayEqual(types, self.charge2.types)
        self.assertArrayEqual(num_atoms, self.charge2.num_atoms)

    def test_get_structure_with_charges(self):
        structure_dict2 = {
            "lattice": {
                "c": 3.198244,
                "volume": 23.132361565928807,
                "b": 3.1982447183003364,
                "gamma": 60.00000011873414,
                "beta": 60.00000401737447,
                "alpha": 60.00000742944491,
                "matrix": [
                    [2.769761, 0.0, 1.599122],
                    [0.923254, 2.611356, 1.599122],
                    [0.0, 0.0, 3.198244],
                ],
                "a": 3.1982443884113985,
            },
            "@class": "Structure",
            "sites": [
                {
                    "xyz": [1.846502883732, 1.305680611356, 3.198248797366],
                    "properties": {"Loewdin Charges": -1.25, "Mulliken Charges": -1.3},
                    "abc": [0.499998, 0.500001, 0.500002],
                    "species": [{"occu": 1, "element": "O"}],
                    "label": "O",
                },
                {
                    "xyz": [0.0, 0.0, 0.0],
                    "properties": {"Loewdin Charges": 1.25, "Mulliken Charges": 1.3},
                    "abc": [0.0, 0.0, 0.0],
                    "species": [{"occu": 1, "element": "Mn"}],
                    "label": "Mn",
                },
            ],
            "charge": None,
            "@module": "pymatgen.core.structure",
        }
        s2 = Structure.from_dict(structure_dict2)
        self.assertEqual(
            s2,
            self.charge2.get_structure_with_charges(os.path.join(this_dir, "../../tests/POSCAR.MnO")),
        )


class LobsteroutTest(PymatgenTest):
    def setUp(self):
        warnings.simplefilter("ignore")
        self.lobsterout_normal = Lobsterout(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterout.normal")
        )
        # make sure .gz files are also read correctly
        self.lobsterout_normal = Lobsterout(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterout.normal2.gz")
        )
        self.lobsterout_fatband_grosspop_densityofenergies = Lobsterout(
            filename=os.path.join(
                PymatgenTest.TEST_FILES_DIR,
                "cohp",
                "lobsterout.fatband_grosspop_densityofenergy",
            )
        )
        self.lobsterout_saveprojection = Lobsterout(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterout.saveprojection")
        )
        self.lobsterout_skipping_all = Lobsterout(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterout.skipping_all")
        )
        self.lobsterout_twospins = Lobsterout(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterout.twospins")
        )
        self.lobsterout_GaAs = Lobsterout(filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterout.GaAs"))
        self.lobsterout_from_projection = Lobsterout(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterout_from_projection")
        )
        self.lobsterout_onethread = Lobsterout(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterout.onethread")
        )

    def tearDown(self):
        warnings.simplefilter("default")

    def testattributes(self):
        self.assertListEqual(
            self.lobsterout_normal.basis_functions,
            [
                [
                    "3s",
                    "4s",
                    "3p_y",
                    "3p_z",
                    "3p_x",
                    "3d_xy",
                    "3d_yz",
                    "3d_z^2",
                    "3d_xz",
                    "3d_x^2-y^2",
                ]
            ],
        )
        self.assertListEqual(self.lobsterout_normal.basis_type, ["pbeVaspFit2015"])
        self.assertListEqual(self.lobsterout_normal.chargespilling, [0.0268])
        self.assertEqual(self.lobsterout_normal.dftprogram, "VASP")
        self.assertListEqual(self.lobsterout_normal.elements, ["Ti"])
        self.assertTrue(self.lobsterout_normal.has_CHARGE)
        self.assertTrue(self.lobsterout_normal.has_COHPCAR)
        self.assertTrue(self.lobsterout_normal.has_COOPCAR)
        self.assertTrue(self.lobsterout_normal.has_DOSCAR)
        self.assertFalse(self.lobsterout_normal.has_Projection)
        self.assertTrue(self.lobsterout_normal.has_bandoverlaps)
        self.assertFalse(self.lobsterout_normal.has_density_of_energies)
        self.assertFalse(self.lobsterout_normal.has_fatbands)
        self.assertFalse(self.lobsterout_normal.has_grosspopulation)
        self.assertListEqual(
            self.lobsterout_normal.info_lines,
            [
                "There are more PAW bands than local basis functions available.",
                "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
                "the PAW bands from 21 and upwards will be ignored.",
            ],
        )
        self.assertListEqual(
            self.lobsterout_normal.info_orthonormalization,
            ["3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."],
        )
        self.assertFalse(self.lobsterout_normal.is_restart_from_projection)
        self.assertEqual(self.lobsterout_normal.lobster_version, "v3.1.0")
        self.assertEqual(self.lobsterout_normal.number_of_spins, 1)
        self.assertEqual(self.lobsterout_normal.number_of_threads, 8)
        self.assertDictEqual(
            self.lobsterout_normal.timing,
            {
                "walltime": {"h": "0", "min": "0", "s": "2", "ms": "702"},
                "usertime": {"h": "0", "min": "0", "s": "20", "ms": "330"},
                "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "310"},
            },
        )
        self.assertAlmostEqual(self.lobsterout_normal.totalspilling[0], [0.044000000000000004][0])
        self.assertListEqual(
            self.lobsterout_normal.warninglines,
            [
                "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
                "Generally, this is not a critical error. But to help you analyze it,",
                "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
                "Please check how much they deviate from the identity matrix and decide to",
                "use your results only, if you are sure that this is ok.",
            ],
        )

        self.assertListEqual(
            self.lobsterout_fatband_grosspop_densityofenergies.basis_functions,
            [
                [
                    "3s",
                    "4s",
                    "3p_y",
                    "3p_z",
                    "3p_x",
                    "3d_xy",
                    "3d_yz",
                    "3d_z^2",
                    "3d_xz",
                    "3d_x^2-y^2",
                ]
            ],
        )
        self.assertListEqual(
            self.lobsterout_fatband_grosspop_densityofenergies.basis_type,
            ["pbeVaspFit2015"],
        )
        self.assertListEqual(self.lobsterout_fatband_grosspop_densityofenergies.chargespilling, [0.0268])
        self.assertEqual(self.lobsterout_fatband_grosspop_densityofenergies.dftprogram, "VASP")
        self.assertListEqual(self.lobsterout_fatband_grosspop_densityofenergies.elements, ["Ti"])
        self.assertTrue(self.lobsterout_fatband_grosspop_densityofenergies.has_CHARGE)
        self.assertFalse(self.lobsterout_fatband_grosspop_densityofenergies.has_COHPCAR)
        self.assertFalse(self.lobsterout_fatband_grosspop_densityofenergies.has_COOPCAR)
        self.assertFalse(self.lobsterout_fatband_grosspop_densityofenergies.has_DOSCAR)
        self.assertFalse(self.lobsterout_fatband_grosspop_densityofenergies.has_Projection)
        self.assertTrue(self.lobsterout_fatband_grosspop_densityofenergies.has_bandoverlaps)
        self.assertTrue(self.lobsterout_fatband_grosspop_densityofenergies.has_density_of_energies)
        self.assertTrue(self.lobsterout_fatband_grosspop_densityofenergies.has_fatbands)
        self.assertTrue(self.lobsterout_fatband_grosspop_densityofenergies.has_grosspopulation)
        self.assertListEqual(
            self.lobsterout_fatband_grosspop_densityofenergies.info_lines,
            [
                "There are more PAW bands than local basis functions available.",
                "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
                "the PAW bands from 21 and upwards will be ignored.",
            ],
        )
        self.assertListEqual(
            self.lobsterout_fatband_grosspop_densityofenergies.info_orthonormalization,
            ["3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."],
        )
        self.assertFalse(self.lobsterout_fatband_grosspop_densityofenergies.is_restart_from_projection)
        self.assertEqual(self.lobsterout_fatband_grosspop_densityofenergies.lobster_version, "v3.1.0")
        self.assertEqual(self.lobsterout_fatband_grosspop_densityofenergies.number_of_spins, 1)
        self.assertEqual(self.lobsterout_fatband_grosspop_densityofenergies.number_of_threads, 8)
        self.assertDictEqual(
            self.lobsterout_fatband_grosspop_densityofenergies.timing,
            {
                "walltime": {"h": "0", "min": "0", "s": "4", "ms": "136"},
                "usertime": {"h": "0", "min": "0", "s": "18", "ms": "280"},
                "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "290"},
            },
        )
        self.assertAlmostEqual(
            self.lobsterout_fatband_grosspop_densityofenergies.totalspilling[0],
            [0.044000000000000004][0],
        )
        self.assertListEqual(
            self.lobsterout_fatband_grosspop_densityofenergies.warninglines,
            [
                "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
                "Generally, this is not a critical error. But to help you analyze it,",
                "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
                "Please check how much they deviate from the identity matrix and decide to",
                "use your results only, if you are sure that this is ok.",
            ],
        )

        self.assertListEqual(
            self.lobsterout_saveprojection.basis_functions,
            [
                [
                    "3s",
                    "4s",
                    "3p_y",
                    "3p_z",
                    "3p_x",
                    "3d_xy",
                    "3d_yz",
                    "3d_z^2",
                    "3d_xz",
                    "3d_x^2-y^2",
                ]
            ],
        )
        self.assertListEqual(self.lobsterout_saveprojection.basis_type, ["pbeVaspFit2015"])
        self.assertListEqual(self.lobsterout_saveprojection.chargespilling, [0.0268])
        self.assertEqual(self.lobsterout_saveprojection.dftprogram, "VASP")
        self.assertListEqual(self.lobsterout_saveprojection.elements, ["Ti"])
        self.assertTrue(self.lobsterout_saveprojection.has_CHARGE)
        self.assertFalse(self.lobsterout_saveprojection.has_COHPCAR)
        self.assertFalse(self.lobsterout_saveprojection.has_COOPCAR)
        self.assertFalse(self.lobsterout_saveprojection.has_DOSCAR)
        self.assertTrue(self.lobsterout_saveprojection.has_Projection)
        self.assertTrue(self.lobsterout_saveprojection.has_bandoverlaps)
        self.assertTrue(self.lobsterout_saveprojection.has_density_of_energies)
        self.assertFalse(self.lobsterout_saveprojection.has_fatbands)
        self.assertFalse(self.lobsterout_saveprojection.has_grosspopulation)
        self.assertListEqual(
            self.lobsterout_saveprojection.info_lines,
            [
                "There are more PAW bands than local basis functions available.",
                "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
                "the PAW bands from 21 and upwards will be ignored.",
            ],
        )
        self.assertListEqual(
            self.lobsterout_saveprojection.info_orthonormalization,
            ["3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."],
        )
        self.assertFalse(self.lobsterout_saveprojection.is_restart_from_projection)
        self.assertEqual(self.lobsterout_saveprojection.lobster_version, "v3.1.0")
        self.assertEqual(self.lobsterout_saveprojection.number_of_spins, 1)
        self.assertEqual(self.lobsterout_saveprojection.number_of_threads, 8)
        self.assertDictEqual(
            self.lobsterout_saveprojection.timing,
            {
                "walltime": {"h": "0", "min": "0", "s": "2", "ms": "574"},
                "usertime": {"h": "0", "min": "0", "s": "18", "ms": "250"},
                "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "320"},
            },
        )
        self.assertAlmostEqual(self.lobsterout_saveprojection.totalspilling[0], [0.044000000000000004][0])
        self.assertListEqual(
            self.lobsterout_saveprojection.warninglines,
            [
                "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
                "Generally, this is not a critical error. But to help you analyze it,",
                "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
                "Please check how much they deviate from the identity matrix and decide to",
                "use your results only, if you are sure that this is ok.",
            ],
        )

        self.assertListEqual(
            self.lobsterout_skipping_all.basis_functions,
            [
                [
                    "3s",
                    "4s",
                    "3p_y",
                    "3p_z",
                    "3p_x",
                    "3d_xy",
                    "3d_yz",
                    "3d_z^2",
                    "3d_xz",
                    "3d_x^2-y^2",
                ]
            ],
        )
        self.assertListEqual(self.lobsterout_skipping_all.basis_type, ["pbeVaspFit2015"])
        self.assertListEqual(self.lobsterout_skipping_all.chargespilling, [0.0268])
        self.assertEqual(self.lobsterout_skipping_all.dftprogram, "VASP")
        self.assertListEqual(self.lobsterout_skipping_all.elements, ["Ti"])
        self.assertFalse(self.lobsterout_skipping_all.has_CHARGE)
        self.assertFalse(self.lobsterout_skipping_all.has_COHPCAR)
        self.assertFalse(self.lobsterout_skipping_all.has_COOPCAR)
        self.assertFalse(self.lobsterout_skipping_all.has_DOSCAR)
        self.assertFalse(self.lobsterout_skipping_all.has_Projection)
        self.assertTrue(self.lobsterout_skipping_all.has_bandoverlaps)
        self.assertFalse(self.lobsterout_skipping_all.has_density_of_energies)
        self.assertFalse(self.lobsterout_skipping_all.has_fatbands)
        self.assertFalse(self.lobsterout_skipping_all.has_grosspopulation)
        self.assertListEqual(
            self.lobsterout_skipping_all.info_lines,
            [
                "There are more PAW bands than local basis functions available.",
                "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
                "the PAW bands from 21 and upwards will be ignored.",
            ],
        )
        self.assertListEqual(
            self.lobsterout_skipping_all.info_orthonormalization,
            ["3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."],
        )
        self.assertFalse(self.lobsterout_skipping_all.is_restart_from_projection)
        self.assertEqual(self.lobsterout_skipping_all.lobster_version, "v3.1.0")
        self.assertEqual(self.lobsterout_skipping_all.number_of_spins, 1)
        self.assertEqual(self.lobsterout_skipping_all.number_of_threads, 8)
        self.assertDictEqual(
            self.lobsterout_skipping_all.timing,
            {
                "walltime": {"h": "0", "min": "0", "s": "2", "ms": "117"},
                "usertime": {"h": "0", "min": "0", "s": "16", "ms": "79"},
                "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "320"},
            },
        )
        self.assertAlmostEqual(self.lobsterout_skipping_all.totalspilling[0], [0.044000000000000004][0])
        self.assertListEqual(
            self.lobsterout_skipping_all.warninglines,
            [
                "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
                "Generally, this is not a critical error. But to help you analyze it,",
                "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
                "Please check how much they deviate from the identity matrix and decide to",
                "use your results only, if you are sure that this is ok.",
            ],
        )

        self.assertListEqual(
            self.lobsterout_twospins.basis_functions,
            [
                [
                    "4s",
                    "4p_y",
                    "4p_z",
                    "4p_x",
                    "3d_xy",
                    "3d_yz",
                    "3d_z^2",
                    "3d_xz",
                    "3d_x^2-y^2",
                ]
            ],
        )
        self.assertListEqual(self.lobsterout_twospins.basis_type, ["pbeVaspFit2015"])
        self.assertAlmostEqual(self.lobsterout_twospins.chargespilling[0], 0.36619999999999997)
        self.assertAlmostEqual(self.lobsterout_twospins.chargespilling[1], 0.36619999999999997)
        self.assertEqual(self.lobsterout_twospins.dftprogram, "VASP")
        self.assertListEqual(self.lobsterout_twospins.elements, ["Ti"])
        self.assertTrue(self.lobsterout_twospins.has_CHARGE)
        self.assertTrue(self.lobsterout_twospins.has_COHPCAR)
        self.assertTrue(self.lobsterout_twospins.has_COOPCAR)
        self.assertTrue(self.lobsterout_twospins.has_DOSCAR)
        self.assertFalse(self.lobsterout_twospins.has_Projection)
        self.assertTrue(self.lobsterout_twospins.has_bandoverlaps)
        self.assertFalse(self.lobsterout_twospins.has_density_of_energies)
        self.assertFalse(self.lobsterout_twospins.has_fatbands)
        self.assertFalse(self.lobsterout_twospins.has_grosspopulation)
        self.assertListEqual(
            self.lobsterout_twospins.info_lines,
            [
                "There are more PAW bands than local basis functions available.",
                "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
                "the PAW bands from 19 and upwards will be ignored.",
            ],
        )
        self.assertListEqual(
            self.lobsterout_twospins.info_orthonormalization,
            ["60 of 294 k-points could not be orthonormalized with an accuracy of 1.0E-5."],
        )
        self.assertFalse(self.lobsterout_twospins.is_restart_from_projection)
        self.assertEqual(self.lobsterout_twospins.lobster_version, "v3.1.0")
        self.assertEqual(self.lobsterout_twospins.number_of_spins, 2)
        self.assertEqual(self.lobsterout_twospins.number_of_threads, 8)
        self.assertDictEqual(
            self.lobsterout_twospins.timing,
            {
                "walltime": {"h": "0", "min": "0", "s": "3", "ms": "71"},
                "usertime": {"h": "0", "min": "0", "s": "22", "ms": "660"},
                "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "310"},
            },
        )
        self.assertAlmostEqual(self.lobsterout_twospins.totalspilling[0], [0.2567][0])
        self.assertAlmostEqual(self.lobsterout_twospins.totalspilling[1], [0.2567][0])
        self.assertListEqual(
            self.lobsterout_twospins.warninglines,
            [
                "60 of 294 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
                "Generally, this is not a critical error. But to help you analyze it,",
                "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
                "Please check how much they deviate from the identity matrix and decide to",
                "use your results only, if you are sure that this is ok.",
            ],
        )

        self.assertListEqual(self.lobsterout_from_projection.basis_functions, [])
        self.assertListEqual(self.lobsterout_from_projection.basis_type, [])
        self.assertAlmostEqual(self.lobsterout_from_projection.chargespilling[0], 0.0177)
        self.assertEqual(self.lobsterout_from_projection.dftprogram, None)
        self.assertListEqual(self.lobsterout_from_projection.elements, [])
        self.assertTrue(self.lobsterout_from_projection.has_CHARGE)
        self.assertTrue(self.lobsterout_from_projection.has_COHPCAR)
        self.assertTrue(self.lobsterout_from_projection.has_COOPCAR)
        self.assertTrue(self.lobsterout_from_projection.has_DOSCAR)
        self.assertFalse(self.lobsterout_from_projection.has_Projection)
        self.assertFalse(self.lobsterout_from_projection.has_bandoverlaps)
        self.assertFalse(self.lobsterout_from_projection.has_density_of_energies)
        self.assertFalse(self.lobsterout_from_projection.has_fatbands)
        self.assertFalse(self.lobsterout_from_projection.has_grosspopulation)
        self.assertListEqual(self.lobsterout_from_projection.info_lines, [])
        self.assertListEqual(self.lobsterout_from_projection.info_orthonormalization, [])
        self.assertTrue(self.lobsterout_from_projection.is_restart_from_projection)
        self.assertEqual(self.lobsterout_from_projection.lobster_version, "v3.1.0")
        self.assertEqual(self.lobsterout_from_projection.number_of_spins, 1)
        self.assertEqual(self.lobsterout_from_projection.number_of_threads, 8)
        self.assertDictEqual(
            self.lobsterout_from_projection.timing,
            {
                "walltime": {"h": "0", "min": "2", "s": "1", "ms": "890"},
                "usertime": {"h": "0", "min": "15", "s": "10", "ms": "530"},
                "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "400"},
            },
        )
        self.assertAlmostEqual(self.lobsterout_from_projection.totalspilling[0], [0.1543][0])
        self.assertListEqual(self.lobsterout_from_projection.warninglines, [])

        self.assertListEqual(
            self.lobsterout_GaAs.basis_functions,
            [
                ["4s", "4p_y", "4p_z", "4p_x"],
                [
                    "4s",
                    "4p_y",
                    "4p_z",
                    "4p_x",
                    "3d_xy",
                    "3d_yz",
                    "3d_z^2",
                    "3d_xz",
                    "3d_x^2-y^2",
                ],
            ],
        )
        self.assertListEqual(self.lobsterout_GaAs.basis_type, ["Bunge", "Bunge"])
        self.assertAlmostEqual(self.lobsterout_GaAs.chargespilling[0], 0.0089)
        self.assertEqual(self.lobsterout_GaAs.dftprogram, "VASP")
        self.assertListEqual(self.lobsterout_GaAs.elements, ["As", "Ga"])
        self.assertTrue(self.lobsterout_GaAs.has_CHARGE)
        self.assertTrue(self.lobsterout_GaAs.has_COHPCAR)
        self.assertTrue(self.lobsterout_GaAs.has_COOPCAR)
        self.assertTrue(self.lobsterout_GaAs.has_DOSCAR)
        self.assertFalse(self.lobsterout_GaAs.has_Projection)
        self.assertFalse(self.lobsterout_GaAs.has_bandoverlaps)
        self.assertFalse(self.lobsterout_GaAs.has_density_of_energies)
        self.assertFalse(self.lobsterout_GaAs.has_fatbands)
        self.assertFalse(self.lobsterout_GaAs.has_grosspopulation)
        self.assertListEqual(
            self.lobsterout_GaAs.info_lines,
            [
                "There are more PAW bands than local basis functions available.",
                "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
                "the PAW bands from 14 and upwards will be ignored.",
            ],
        )
        self.assertListEqual(self.lobsterout_GaAs.info_orthonormalization, [])
        self.assertFalse(self.lobsterout_GaAs.is_restart_from_projection)
        self.assertEqual(self.lobsterout_GaAs.lobster_version, "v3.1.0")
        self.assertEqual(self.lobsterout_GaAs.number_of_spins, 1)
        self.assertEqual(self.lobsterout_GaAs.number_of_threads, 8)
        self.assertDictEqual(
            self.lobsterout_GaAs.timing,
            {
                "walltime": {"h": "0", "min": "0", "s": "2", "ms": "726"},
                "usertime": {"h": "0", "min": "0", "s": "12", "ms": "370"},
                "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "180"},
            },
        )
        self.assertAlmostEqual(self.lobsterout_GaAs.totalspilling[0], [0.0859][0])

        self.assertEqual(self.lobsterout_onethread.number_of_threads, 1)

    def test_get_doc(self):
        comparedict = {
            "restart_from_projection": False,
            "lobster_version": "v3.1.0",
            "threads": 8,
            "Dftprogram": "VASP",
            "chargespilling": [0.0268],
            "totalspilling": [0.044000000000000004],
            "elements": ["Ti"],
            "basistype": ["pbeVaspFit2015"],
            "basisfunctions": [
                [
                    "3s",
                    "4s",
                    "3p_y",
                    "3p_z",
                    "3p_x",
                    "3d_xy",
                    "3d_yz",
                    "3d_z^2",
                    "3d_xz",
                    "3d_x^2-y^2",
                ]
            ],
            "timing": {
                "walltime": {"h": "0", "min": "0", "s": "2", "ms": "702"},
                "usertime": {"h": "0", "min": "0", "s": "20", "ms": "330"},
                "sys_time": {"h": "0", "min": "0", "s": "0", "ms": "310"},
            },
            "warnings": [
                "3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5.",
                "Generally, this is not a critical error. But to help you analyze it,",
                "I dumped the band overlap matrices to the file bandOverlaps.lobster.",
                "Please check how much they deviate from the identity matrix and decide to",
                "use your results only, if you are sure that this is ok.",
            ],
            "orthonormalization": ["3 of 147 k-points could not be orthonormalized with an accuracy of 1.0E-5."],
            "infos": [
                "There are more PAW bands than local basis functions available.",
                "To prevent trouble in orthonormalization and Hamiltonian reconstruction",
                "the PAW bands from 21 and upwards will be ignored.",
            ],
            "hasDOSCAR": True,
            "hasCOHPCAR": True,
            "hasCOOPCAR": True,
            "hasCHARGE": True,
            "hasProjection": False,
            "hasbandoverlaps": True,
            "hasfatband": False,
            "hasGrossPopuliation": False,
            "hasDensityOfEnergies": False,
        }
        for key, item in self.lobsterout_normal.get_doc().items():
            if isinstance(item, str):
                self.assertTrue(comparedict[key], item)
            elif isinstance(item, int):
                self.assertEqual(comparedict[key], item)
            elif key in ("chargespilling", "totalspilling"):
                self.assertAlmostEqual(item[0], comparedict[key][0])
            elif isinstance(item, list):
                self.assertListEqual(item, comparedict[key])
            elif isinstance(item, dict):
                self.assertDictEqual(item, comparedict[key])


class FatbandTest(PymatgenTest):
    def setUp(self):
        warnings.simplefilter("ignore")
        self.fatband_SiO2_p_x = Fatband(
            filenames=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p_x"),
            Kpointsfile=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p_x/KPOINTS"),
            vasprun=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p_x/vasprun.xml"),
        )
        self.vasprun_SiO2_p_x = Vasprun(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p_x/vasprun.xml")
        )
        self.bs_symmline = self.vasprun_SiO2_p_x.get_band_structure(line_mode=True, force_hybrid_mode=True)
        self.fatband_SiO2_p = Fatband(
            filenames=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p"),
            Kpointsfile=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p/KPOINTS"),
            vasprun=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p/vasprun.xml"),
        )
        self.vasprun_SiO2_p = Vasprun(
            filename=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p/vasprun.xml")
        )
        self.bs_symmline2 = self.vasprun_SiO2_p.get_band_structure(line_mode=True, force_hybrid_mode=True)
        self.fatband_SiO2_spin = Fatband(
            filenames=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_Spin"),
            Kpointsfile=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_Spin/KPOINTS"),
            vasprun=os.path.join(
                PymatgenTest.TEST_FILES_DIR,
                "cohp",
                "Fatband_SiO2/Test_Spin/vasprun.xml",
            ),
        )
        self.vasprun_SiO2_spin = Vasprun(
            filename=os.path.join(
                PymatgenTest.TEST_FILES_DIR,
                "cohp",
                "Fatband_SiO2/Test_Spin/vasprun.xml",
            )
        )
        self.bs_symmline_spin = self.vasprun_SiO2_p.get_band_structure(line_mode=True, force_hybrid_mode=True)

    def tearDown(self):
        warnings.simplefilter("default")

    def test_attributes(self):
        self.assertAlmostEqual(list(self.fatband_SiO2_p_x.label_dict["M"])[0], 0.5)
        self.assertAlmostEqual(list(self.fatband_SiO2_p_x.label_dict["M"])[1], 0.0)
        self.assertAlmostEqual(list(self.fatband_SiO2_p_x.label_dict["M"])[2], 0.0)
        self.assertEqual(self.fatband_SiO2_p_x.efermi, self.vasprun_SiO2_p_x.efermi)
        lattice1 = self.bs_symmline.lattice_rec.as_dict()
        lattice2 = self.fatband_SiO2_p_x.lattice.as_dict()
        self.assertAlmostEqual(lattice1["matrix"][0][0], lattice2["matrix"][0][0])
        self.assertAlmostEqual(lattice1["matrix"][0][1], lattice2["matrix"][0][1])
        self.assertAlmostEqual(lattice1["matrix"][0][2], lattice2["matrix"][0][2])
        self.assertAlmostEqual(lattice1["matrix"][1][0], lattice2["matrix"][1][0])
        self.assertAlmostEqual(lattice1["matrix"][1][1], lattice2["matrix"][1][1])
        self.assertAlmostEqual(lattice1["matrix"][1][2], lattice2["matrix"][1][2])
        self.assertAlmostEqual(lattice1["matrix"][2][0], lattice2["matrix"][2][0])
        self.assertAlmostEqual(lattice1["matrix"][2][1], lattice2["matrix"][2][1])
        self.assertAlmostEqual(lattice1["matrix"][2][2], lattice2["matrix"][2][2])
        self.assertEqual(
            self.fatband_SiO2_p_x.eigenvals[Spin.up][1][1] - self.fatband_SiO2_p_x.efermi,
            -18.245,
        )
        self.assertEqual(self.fatband_SiO2_p_x.is_spinpolarized, False)
        self.assertAlmostEqual(self.fatband_SiO2_p_x.kpoints_array[3][0], 0.03409091)
        self.assertEqual(self.fatband_SiO2_p_x.kpoints_array[3][1], 0.0)
        self.assertEqual(self.fatband_SiO2_p_x.kpoints_array[3][2], 0.0)
        self.assertEqual(self.fatband_SiO2_p_x.nbands, 36)
        self.assertEqual(self.fatband_SiO2_p_x.p_eigenvals[Spin.up][2][1]["Si1"]["3p_x"], 0.002)
        self.assertAlmostEqual(self.fatband_SiO2_p_x.structure[0].frac_coords[0], 0.0)
        self.assertAlmostEqual(self.fatband_SiO2_p_x.structure[0].frac_coords[1], 0.47634315)
        self.assertAlmostEqual(self.fatband_SiO2_p_x.structure[0].frac_coords[2], 0.666667)
        self.assertEqual(self.fatband_SiO2_p_x.structure[0].species_string, "Si")
        self.assertAlmostEqual(self.fatband_SiO2_p_x.structure[0].coords[0], -1.19607309)
        self.assertAlmostEqual(self.fatband_SiO2_p_x.structure[0].coords[1], 2.0716597)
        self.assertAlmostEqual(self.fatband_SiO2_p_x.structure[0].coords[2], 3.67462144)

        self.assertAlmostEqual(list(self.fatband_SiO2_p.label_dict["M"])[0], 0.5)
        self.assertAlmostEqual(list(self.fatband_SiO2_p.label_dict["M"])[1], 0.0)
        self.assertAlmostEqual(list(self.fatband_SiO2_p.label_dict["M"])[2], 0.0)
        self.assertEqual(self.fatband_SiO2_p.efermi, self.vasprun_SiO2_p.efermi)
        lattice1 = self.bs_symmline2.lattice_rec.as_dict()
        lattice2 = self.fatband_SiO2_p.lattice.as_dict()
        self.assertAlmostEqual(lattice1["matrix"][0][0], lattice2["matrix"][0][0])
        self.assertAlmostEqual(lattice1["matrix"][0][1], lattice2["matrix"][0][1])
        self.assertAlmostEqual(lattice1["matrix"][0][2], lattice2["matrix"][0][2])
        self.assertAlmostEqual(lattice1["matrix"][1][0], lattice2["matrix"][1][0])
        self.assertAlmostEqual(lattice1["matrix"][1][1], lattice2["matrix"][1][1])
        self.assertAlmostEqual(lattice1["matrix"][1][2], lattice2["matrix"][1][2])
        self.assertAlmostEqual(lattice1["matrix"][2][0], lattice2["matrix"][2][0])
        self.assertAlmostEqual(lattice1["matrix"][2][1], lattice2["matrix"][2][1])
        self.assertAlmostEqual(lattice1["matrix"][2][2], lattice2["matrix"][2][2])
        self.assertEqual(
            self.fatband_SiO2_p.eigenvals[Spin.up][1][1] - self.fatband_SiO2_p.efermi,
            -18.245,
        )
        self.assertEqual(self.fatband_SiO2_p.is_spinpolarized, False)
        self.assertAlmostEqual(self.fatband_SiO2_p.kpoints_array[3][0], 0.03409091)
        self.assertEqual(self.fatband_SiO2_p.kpoints_array[3][1], 0.0)
        self.assertEqual(self.fatband_SiO2_p.kpoints_array[3][2], 0.0)
        self.assertEqual(self.fatband_SiO2_p.nbands, 36)
        self.assertEqual(self.fatband_SiO2_p.p_eigenvals[Spin.up][2][1]["Si1"]["3p"], 0.042)
        self.assertAlmostEqual(self.fatband_SiO2_p.structure[0].frac_coords[0], 0.0)
        self.assertAlmostEqual(self.fatband_SiO2_p.structure[0].frac_coords[1], 0.47634315)
        self.assertAlmostEqual(self.fatband_SiO2_p.structure[0].frac_coords[2], 0.666667)
        self.assertEqual(self.fatband_SiO2_p.structure[0].species_string, "Si")
        self.assertAlmostEqual(self.fatband_SiO2_p.structure[0].coords[0], -1.19607309)
        self.assertAlmostEqual(self.fatband_SiO2_p.structure[0].coords[1], 2.0716597)
        self.assertAlmostEqual(self.fatband_SiO2_p.structure[0].coords[2], 3.67462144)

        self.assertAlmostEqual(list(self.fatband_SiO2_spin.label_dict["M"])[0], 0.5)
        self.assertAlmostEqual(list(self.fatband_SiO2_spin.label_dict["M"])[1], 0.0)
        self.assertAlmostEqual(list(self.fatband_SiO2_spin.label_dict["M"])[2], 0.0)
        self.assertEqual(self.fatband_SiO2_spin.efermi, self.vasprun_SiO2_spin.efermi)
        lattice1 = self.bs_symmline_spin.lattice_rec.as_dict()
        lattice2 = self.fatband_SiO2_spin.lattice.as_dict()
        self.assertAlmostEqual(lattice1["matrix"][0][0], lattice2["matrix"][0][0])
        self.assertAlmostEqual(lattice1["matrix"][0][1], lattice2["matrix"][0][1])
        self.assertAlmostEqual(lattice1["matrix"][0][2], lattice2["matrix"][0][2])
        self.assertAlmostEqual(lattice1["matrix"][1][0], lattice2["matrix"][1][0])
        self.assertAlmostEqual(lattice1["matrix"][1][1], lattice2["matrix"][1][1])
        self.assertAlmostEqual(lattice1["matrix"][1][2], lattice2["matrix"][1][2])
        self.assertAlmostEqual(lattice1["matrix"][2][0], lattice2["matrix"][2][0])
        self.assertAlmostEqual(lattice1["matrix"][2][1], lattice2["matrix"][2][1])
        self.assertAlmostEqual(lattice1["matrix"][2][2], lattice2["matrix"][2][2])
        self.assertEqual(
            self.fatband_SiO2_spin.eigenvals[Spin.up][1][1] - self.fatband_SiO2_spin.efermi,
            -18.245,
        )
        self.assertEqual(
            self.fatband_SiO2_spin.eigenvals[Spin.down][1][1] - self.fatband_SiO2_spin.efermi,
            -18.245,
        )
        self.assertEqual(self.fatband_SiO2_spin.is_spinpolarized, True)
        self.assertAlmostEqual(self.fatband_SiO2_spin.kpoints_array[3][0], 0.03409091)
        self.assertEqual(self.fatband_SiO2_spin.kpoints_array[3][1], 0.0)
        self.assertEqual(self.fatband_SiO2_spin.kpoints_array[3][2], 0.0)
        self.assertEqual(self.fatband_SiO2_spin.nbands, 36)

        self.assertEqual(self.fatband_SiO2_spin.p_eigenvals[Spin.up][2][1]["Si1"]["3p"], 0.042)
        self.assertAlmostEqual(self.fatband_SiO2_spin.structure[0].frac_coords[0], 0.0)
        self.assertAlmostEqual(self.fatband_SiO2_spin.structure[0].frac_coords[1], 0.47634315)
        self.assertAlmostEqual(self.fatband_SiO2_spin.structure[0].frac_coords[2], 0.666667)
        self.assertEqual(self.fatband_SiO2_spin.structure[0].species_string, "Si")
        self.assertAlmostEqual(self.fatband_SiO2_spin.structure[0].coords[0], -1.19607309)
        self.assertAlmostEqual(self.fatband_SiO2_spin.structure[0].coords[1], 2.0716597)
        self.assertAlmostEqual(self.fatband_SiO2_spin.structure[0].coords[2], 3.67462144)

    def test_raises(self):
        with self.assertRaises(ValueError):
            self.fatband_SiO2_p_x = Fatband(
                filenames=[
                    os.path.join(
                        PymatgenTest.TEST_FILES_DIR,
                        "cohp",
                        "Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                    ),
                    os.path.join(
                        PymatgenTest.TEST_FILES_DIR,
                        "cohp",
                        "Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                    ),
                ],
                Kpointsfile=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p_x/KPOINTS"),
                vasprun=os.path.join(
                    PymatgenTest.TEST_FILES_DIR,
                    "cohp",
                    "Fatband_SiO2/Test_p_x/vasprun.xml",
                ),
            )

        with self.assertRaises(ValueError):
            self.fatband_SiO2_p_x = Fatband(
                filenames=[
                    os.path.join(
                        PymatgenTest.TEST_FILES_DIR,
                        "cohp",
                        "Fatband_SiO2/Test_p_x/FATBAND_si1_3p_x.lobster",
                    ),
                    os.path.join(
                        PymatgenTest.TEST_FILES_DIR,
                        "cohp",
                        "Fatband_SiO2/Test_p/FATBAND_si1_3p.lobster",
                    ),
                ],
                Kpointsfile=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p_x/KPOINTS"),
                vasprun=os.path.join(
                    PymatgenTest.TEST_FILES_DIR,
                    "cohp",
                    "Fatband_SiO2/Test_p_x/vasprun.xml",
                ),
            )

        with self.assertRaises(ValueError):
            self.fatband_SiO2_p_x = Fatband(
                filenames=".",
                Kpointsfile=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "Fatband_SiO2/Test_p_x/KPOINTS"),
                vasprun=os.path.join(
                    PymatgenTest.TEST_FILES_DIR,
                    "cohp",
                    "Fatband_SiO2/Test_p_x/vasprun.xml",
                ),
            )

    def test_get_bandstructure(self):
        bs_p = self.fatband_SiO2_p.get_bandstructure()
        atom1 = bs_p.structure[0]
        atom2 = self.bs_symmline2.structure[0]
        self.assertAlmostEqual(atom1.frac_coords[0], atom2.frac_coords[0])
        self.assertAlmostEqual(atom1.frac_coords[1], atom2.frac_coords[1])
        self.assertAlmostEqual(atom1.frac_coords[2], atom2.frac_coords[2])
        self.assertAlmostEqual(atom1.coords[0], atom2.coords[0])
        self.assertAlmostEqual(atom1.coords[1], atom2.coords[1])
        self.assertAlmostEqual(atom1.coords[2], atom2.coords[2])
        self.assertEqual(atom1.species_string, atom2.species_string)
        self.assertEqual(bs_p.efermi, self.bs_symmline2.efermi)
        branch1 = bs_p.branches[0]
        branch2 = self.bs_symmline2.branches[0]
        self.assertEqual(branch2["name"], branch1["name"])
        self.assertEqual(branch2["start_index"], branch1["start_index"])
        self.assertEqual(branch2["end_index"], branch1["end_index"])

        self.assertAlmostEqual(bs_p.distance[30], self.bs_symmline2.distance[30])
        lattice1 = bs_p.lattice_rec.as_dict()
        lattice2 = self.bs_symmline2.lattice_rec.as_dict()
        self.assertAlmostEqual(lattice1["matrix"][0][0], lattice2["matrix"][0][0])
        self.assertAlmostEqual(lattice1["matrix"][0][1], lattice2["matrix"][0][1])
        self.assertAlmostEqual(lattice1["matrix"][0][2], lattice2["matrix"][0][2])
        self.assertAlmostEqual(lattice1["matrix"][1][0], lattice2["matrix"][1][0])
        self.assertAlmostEqual(lattice1["matrix"][1][1], lattice2["matrix"][1][1])
        self.assertAlmostEqual(lattice1["matrix"][1][2], lattice2["matrix"][1][2])
        self.assertAlmostEqual(lattice1["matrix"][2][0], lattice2["matrix"][2][0])
        self.assertAlmostEqual(lattice1["matrix"][2][1], lattice2["matrix"][2][1])
        self.assertAlmostEqual(lattice1["matrix"][2][2], lattice2["matrix"][2][2])

        self.assertAlmostEqual(bs_p.kpoints[8].frac_coords[0], self.bs_symmline2.kpoints[8].frac_coords[0])
        self.assertAlmostEqual(bs_p.kpoints[8].frac_coords[1], self.bs_symmline2.kpoints[8].frac_coords[1])
        self.assertAlmostEqual(bs_p.kpoints[8].frac_coords[2], self.bs_symmline2.kpoints[8].frac_coords[2])
        self.assertAlmostEqual(bs_p.kpoints[8].cart_coords[0], self.bs_symmline2.kpoints[8].cart_coords[0])
        self.assertAlmostEqual(bs_p.kpoints[8].cart_coords[1], self.bs_symmline2.kpoints[8].cart_coords[1])
        self.assertAlmostEqual(bs_p.kpoints[8].cart_coords[2], self.bs_symmline2.kpoints[8].cart_coords[2])
        self.assertAlmostEqual(
            bs_p.kpoints[50].frac_coords[0],
            self.bs_symmline2.kpoints[50].frac_coords[0],
        )
        self.assertAlmostEqual(
            bs_p.kpoints[50].frac_coords[1],
            self.bs_symmline2.kpoints[50].frac_coords[1],
        )
        self.assertAlmostEqual(
            bs_p.kpoints[50].frac_coords[2],
            self.bs_symmline2.kpoints[50].frac_coords[2],
        )
        self.assertAlmostEqual(
            bs_p.kpoints[50].cart_coords[0],
            self.bs_symmline2.kpoints[50].cart_coords[0],
        )
        self.assertAlmostEqual(
            bs_p.kpoints[50].cart_coords[1],
            self.bs_symmline2.kpoints[50].cart_coords[1],
        )
        self.assertAlmostEqual(
            bs_p.kpoints[50].cart_coords[2],
            self.bs_symmline2.kpoints[50].cart_coords[2],
        )
        self.assertAlmostEqual(
            bs_p.get_band_gap()["energy"],
            self.bs_symmline2.get_band_gap()["energy"],
            places=2,
        )
        self.assertAlmostEqual(bs_p.get_projection_on_elements()[Spin.up][0][0]["Si"], 3 * (0.001 + 0.064))
        self.assertAlmostEqual(
            bs_p.get_projections_on_elements_and_orbitals({"Si": ["3p"]})[Spin.up][0][0]["Si"]["3p"],
            0.003,
        )
        self.assertAlmostEqual(
            bs_p.get_projections_on_elements_and_orbitals({"O": ["2p"]})[Spin.up][0][0]["O"]["2p"],
            0.002 * 3 + 0.003 * 3,
        )
        dict_here = bs_p.get_projections_on_elements_and_orbitals({"Si": ["3s", "3p"], "O": ["2s", "2p"]})[Spin.up][0][
            0
        ]
        self.assertAlmostEqual(dict_here["Si"]["3s"], 0.192)
        self.assertAlmostEqual(dict_here["Si"]["3p"], 0.003)
        self.assertAlmostEqual(dict_here["O"]["2s"], 0.792)
        self.assertAlmostEqual(dict_here["O"]["2p"], 0.015)

        bs_spin = self.fatband_SiO2_spin.get_bandstructure()
        self.assertAlmostEqual(
            bs_spin.get_projection_on_elements()[Spin.up][0][0]["Si"],
            3 * (0.001 + 0.064),
        )
        self.assertAlmostEqual(
            bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3p"]})[Spin.up][0][0]["Si"]["3p"],
            0.003,
        )
        self.assertAlmostEqual(
            bs_spin.get_projections_on_elements_and_orbitals({"O": ["2p"]})[Spin.up][0][0]["O"]["2p"],
            0.002 * 3 + 0.003 * 3,
        )
        dict_here = bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3s", "3p"], "O": ["2s", "2p"]})[Spin.up][
            0
        ][0]
        self.assertAlmostEqual(dict_here["Si"]["3s"], 0.192)
        self.assertAlmostEqual(dict_here["Si"]["3p"], 0.003)
        self.assertAlmostEqual(dict_here["O"]["2s"], 0.792)
        self.assertAlmostEqual(dict_here["O"]["2p"], 0.015)

        self.assertAlmostEqual(
            bs_spin.get_projection_on_elements()[Spin.up][0][0]["Si"],
            3 * (0.001 + 0.064),
        )
        self.assertAlmostEqual(
            bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3p"]})[Spin.down][0][0]["Si"]["3p"],
            0.003,
        )
        self.assertAlmostEqual(
            bs_spin.get_projections_on_elements_and_orbitals({"O": ["2p"]})[Spin.down][0][0]["O"]["2p"],
            0.002 * 3 + 0.003 * 3,
        )
        dict_here = bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3s", "3p"], "O": ["2s", "2p"]})[
            Spin.down
        ][0][0]
        self.assertAlmostEqual(dict_here["Si"]["3s"], 0.192)
        self.assertAlmostEqual(dict_here["Si"]["3p"], 0.003)
        self.assertAlmostEqual(dict_here["O"]["2s"], 0.792)
        self.assertAlmostEqual(dict_here["O"]["2p"], 0.015)
        bs_p_x = self.fatband_SiO2_p_x.get_bandstructure()
        self.assertAlmostEqual(
            bs_p_x.get_projection_on_elements()[Spin.up][0][0]["Si"],
            3 * (0.001 + 0.064),
            2,
        )


class LobsterinTest(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")
        self.Lobsterinfromfile = Lobsterin.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterin.1"))
        self.Lobsterinfromfile2 = Lobsterin.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterin.2"))
        self.Lobsterinfromfile3 = Lobsterin.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterin.3"))
        self.Lobsterinfromfile4 = Lobsterin.from_file(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "lobsterin.4.gz")
        )

    def test_from_file(self):
        # test read from file
        self.assertAlmostEqual(self.Lobsterinfromfile["cohpstartenergy"], -15.0)
        self.assertAlmostEqual(self.Lobsterinfromfile["cohpendenergy"], 5.0)
        self.assertAlmostEqual(self.Lobsterinfromfile["basisset"], "pbeVaspFit2015")
        self.assertAlmostEqual(self.Lobsterinfromfile["gaussiansmearingwidth"], 0.1)
        self.assertEqual(self.Lobsterinfromfile["basisfunctions"][0], "Fe 3d 4p 4s")
        self.assertEqual(self.Lobsterinfromfile["basisfunctions"][1], "Co 3d 4p 4s")
        self.assertEqual(self.Lobsterinfromfile["skipdos"], True)
        self.assertEqual(self.Lobsterinfromfile["skipcohp"], True)
        self.assertEqual(self.Lobsterinfromfile["skipcoop"], True)
        self.assertEqual(self.Lobsterinfromfile["skippopulationanalysis"], True)
        self.assertEqual(self.Lobsterinfromfile["skipgrosspopulation"], True)

        # test if comments are correctly removed
        self.assertDictEqual(self.Lobsterinfromfile, self.Lobsterinfromfile2)

    def test_getitem(self):
        # tests implementation of getitem, should be case independent
        self.assertAlmostEqual(self.Lobsterinfromfile["COHPSTARTENERGY"], -15.0)

    def test_setitem(self):
        # test implementation of setitem
        self.Lobsterinfromfile["skipCOHP"] = False
        self.assertEqual(self.Lobsterinfromfile["skipcohp"], False)

    def test_initialize_from_dict(self):
        # initialize from dict
        lobsterin1 = Lobsterin(
            {
                "cohpstartenergy": -15.0,
                "cohpendenergy": 5.0,
                "basisset": "pbeVaspFit2015",
                "gaussiansmearingwidth": 0.1,
                "basisfunctions": ["Fe 3d 4p 4s", "Co 3d 4p 4s"],
                "skipdos": True,
                "skipcohp": True,
                "skipcoop": True,
                "skippopulationanalysis": True,
                "skipgrosspopulation": True,
            }
        )
        self.assertAlmostEqual(lobsterin1["cohpstartenergy"], -15.0)
        self.assertAlmostEqual(lobsterin1["cohpendenergy"], 5.0)
        self.assertAlmostEqual(lobsterin1["basisset"], "pbeVaspFit2015")
        self.assertAlmostEqual(lobsterin1["gaussiansmearingwidth"], 0.1)
        self.assertEqual(lobsterin1["basisfunctions"][0], "Fe 3d 4p 4s")
        self.assertEqual(lobsterin1["basisfunctions"][1], "Co 3d 4p 4s")
        self.assertEqual(lobsterin1["skipdos"], True)
        self.assertEqual(lobsterin1["skipcohp"], True)
        self.assertEqual(lobsterin1["skipcoop"], True)
        self.assertEqual(lobsterin1["skippopulationanalysis"], True)
        self.assertEqual(lobsterin1["skipgrosspopulation"], True)
        with self.assertRaises(IOError):
            lobsterin2 = Lobsterin({"cohpstartenergy": -15.0, "cohpstartEnergy": -20.0})
        lobsterin2 = Lobsterin({"cohpstartenergy": -15.0})
        # can only calculate nbands if basis functions are provided
        with self.assertRaises(IOError):
            lobsterin2._get_nbands(structure=Structure.from_file(os.path.join(test_dir_doscar, "POSCAR.Fe3O4")))

    def test_standard_settings(self):
        # test standard settings
        for option in [
            "standard",
            "standard_from_projection",
            "standard_with_fatband",
            "onlyprojection",
            "onlydos",
            "onlycohp",
            "onlycoop",
            "onlycohpcoop",
        ]:

            lobsterin1 = Lobsterin.standard_calculations_from_vasp_files(
                os.path.join(test_dir_doscar, "POSCAR.Fe3O4"),
                os.path.join(test_dir_doscar, "INCAR.lobster"),
                os.path.join(test_dir_doscar, "POTCAR.Fe3O4"),
                option=option,
            )
            self.assertAlmostEqual(lobsterin1["cohpstartenergy"], -15.0)
            self.assertAlmostEqual(lobsterin1["cohpendenergy"], 5.0)
            self.assertAlmostEqual(lobsterin1["basisset"], "pbeVaspFit2015")
            self.assertAlmostEqual(lobsterin1["gaussiansmearingwidth"], 0.1)
            self.assertEqual(lobsterin1["basisfunctions"][0], "Fe 3d 4p 4s ")
            self.assertEqual(lobsterin1["basisfunctions"][1], "O 2p 2s ")

            if option in [
                "standard",
                "standard_with_fatband",
                "onlyprojection",
                "onlycohp",
                "onlycoop",
                "onlycohpcoop",
            ]:
                self.assertEqual(lobsterin1["saveProjectiontoFile"], True)
            if option in [
                "standard",
                "standard_with_fatband",
                "onlycohp",
                "onlycoop",
                "onlycohpcoop",
            ]:
                self.assertEqual(lobsterin1["cohpGenerator"], "from 0.1 to 6.0 orbitalwise")
            if option in ["standard"]:
                self.assertEqual("skipdos" not in lobsterin1, True)
                self.assertEqual("skipcohp" not in lobsterin1, True)
                self.assertEqual("skipcoop" not in lobsterin1, True)
            if option in ["standard_with_fatband"]:
                self.assertListEqual(lobsterin1["createFatband"], ["Fe 3d 4p 4s ", "O 2p 2s "])
                self.assertEqual("skipdos" not in lobsterin1, True)
                self.assertEqual("skipcohp" not in lobsterin1, True)
                self.assertEqual("skipcoop" not in lobsterin1, True)
            if option in ["standard_from_projection"]:
                self.assertTrue(lobsterin1["loadProjectionFromFile"], True)
            if option in ["onlyprojection", "onlycohp", "onlycoop", "onlycohpcoop"]:
                self.assertTrue(lobsterin1["skipdos"], True)
                self.assertTrue(lobsterin1["skipPopulationAnalysis"], True)
                self.assertTrue(lobsterin1["skipGrossPopulation"], True)
            if option in ["onlydos"]:
                self.assertTrue(lobsterin1["skipPopulationAnalysis"], True)
                self.assertTrue(lobsterin1["skipGrossPopulation"], True)
                self.assertTrue(lobsterin1["skipcohp"], True)
                self.assertTrue(lobsterin1["skipcoop"], True)
            if option in ["onlycohp"]:
                self.assertTrue(lobsterin1["skipcoop"], True)
            if option in ["onlycoop"]:
                self.assertTrue(lobsterin1["skipcohp"], True)
            if option in ["onlyprojection"]:
                self.assertTrue(lobsterin1["skipdos"], True)

        # test basis functions by dict
        lobsterin_new = Lobsterin.standard_calculations_from_vasp_files(
            os.path.join(test_dir_doscar, "POSCAR.Fe3O4"),
            os.path.join(test_dir_doscar, "INCAR.lobster"),
            dict_for_basis={"Fe": "3d 4p 4s", "O": "2s 2p"},
            option="standard",
        )
        self.assertListEqual(lobsterin_new["basisfunctions"], ["Fe 3d 4p 4s", "O 2s 2p"])

        # test gaussian smearing
        lobsterin_new = Lobsterin.standard_calculations_from_vasp_files(
            os.path.join(test_dir_doscar, "POSCAR.Fe3O4"),
            os.path.join(test_dir_doscar, "INCAR.lobster2"),
            dict_for_basis={"Fe": "3d 4p 4s", "O": "2s 2p"},
            option="standard",
        )
        self.assertTrue("gaussiansmearingwidth" not in lobsterin_new)

        # fatband and ISMEAR=-5 does not work together
        with self.assertRaises(ValueError):
            lobsterin_new = Lobsterin.standard_calculations_from_vasp_files(
                os.path.join(test_dir_doscar, "POSCAR.Fe3O4"),
                os.path.join(test_dir_doscar, "INCAR.lobster2"),
                dict_for_basis={"Fe": "3d 4p 4s", "O": "2s 2p"},
                option="standard_with_fatband",
            )

    def test_diff(self):
        # test diff
        self.assertDictEqual(self.Lobsterinfromfile.diff(self.Lobsterinfromfile2)["Different"], {})
        self.assertAlmostEqual(
            self.Lobsterinfromfile.diff(self.Lobsterinfromfile2)["Same"]["COHPSTARTENERGY"],
            -15.0,
        )

        # test diff in both directions
        for entry in self.Lobsterinfromfile.diff(self.Lobsterinfromfile3)["Same"].keys():
            self.assertTrue(entry in self.Lobsterinfromfile3.diff(self.Lobsterinfromfile)["Same"].keys())
        for entry in self.Lobsterinfromfile3.diff(self.Lobsterinfromfile)["Same"].keys():
            self.assertTrue(entry in self.Lobsterinfromfile.diff(self.Lobsterinfromfile3)["Same"].keys())
        for entry in self.Lobsterinfromfile.diff(self.Lobsterinfromfile3)["Different"].keys():
            self.assertTrue(entry in self.Lobsterinfromfile3.diff(self.Lobsterinfromfile)["Different"].keys())
        for entry in self.Lobsterinfromfile3.diff(self.Lobsterinfromfile)["Different"].keys():
            self.assertTrue(entry in self.Lobsterinfromfile.diff(self.Lobsterinfromfile3)["Different"].keys())

        self.assertEqual(
            self.Lobsterinfromfile.diff(self.Lobsterinfromfile3)["Different"]["SKIPCOHP"]["lobsterin1"],
            self.Lobsterinfromfile3.diff(self.Lobsterinfromfile)["Different"]["SKIPCOHP"]["lobsterin2"],
        )

    def test_get_basis(self):
        # get basis functions
        lobsterin1 = Lobsterin({})
        potcar = Potcar.from_file(os.path.join(test_dir_doscar, "POTCAR.Fe3O4"))
        Potcar_names = [name["symbol"] for name in potcar.spec]

        self.assertListEqual(
            lobsterin1.get_basis(
                Structure.from_file(os.path.join(test_dir_doscar, "Fe3O4.cif")),
                potcar_symbols=Potcar_names,
            ),
            ["Fe 3d 4p 4s ", "O 2p 2s "],
        )
        potcar = Potcar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "POTCAR.GaAs"))
        Potcar_names = [name["symbol"] for name in potcar.spec]
        self.assertListEqual(
            lobsterin1.get_basis(
                Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "POSCAR.GaAs")),
                potcar_symbols=Potcar_names,
            ),
            ["Ga 3d 4p 4s ", "As 4p 4s "],
        )

    def test_get_all_possible_basis_functions(self):
        potcar = Potcar.from_file(os.path.join(test_dir_doscar, "POTCAR.Fe3O4"))
        Potcar_names = [name["symbol"] for name in potcar.spec]
        result = Lobsterin.get_all_possible_basis_functions(
            Structure.from_file(os.path.join(test_dir_doscar, "Fe3O4.cif")),
            potcar_symbols=Potcar_names,
        )
        self.assertDictEqual(result[0], {"Fe": "3d 4s", "O": "2p 2s"})
        self.assertDictEqual(result[1], {"Fe": "3d 4s 4p", "O": "2p 2s"})

        potcar2 = Potcar.from_file(os.path.join(test_dir_doscar, "POT_GGA_PAW_PBE_54/POTCAR.Fe_pv.gz"))
        Potcar_names2 = [name["symbol"] for name in potcar2.spec]
        result2 = Lobsterin.get_all_possible_basis_functions(
            Structure.from_file(os.path.join(test_dir_doscar, "Fe.cif")),
            potcar_symbols=Potcar_names2,
        )
        self.assertDictEqual(result2[0], {"Fe": "3d 3p 4s"})

    def test_get_potcar_symbols(self):
        lobsterin1 = Lobsterin({})
        self.assertListEqual(
            lobsterin1._get_potcar_symbols(os.path.join(test_dir_doscar, "POTCAR.Fe3O4")),
            ["Fe", "O"],
        )
        self.assertListEqual(
            lobsterin1._get_potcar_symbols(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "POTCAR.GaAs")),
            ["Ga_d", "As"],
        )

    def test_write_lobsterin(self):
        # write lobsterin, read it and compare it
        outfile_path = tempfile.mkstemp()[1]
        lobsterin1 = Lobsterin.standard_calculations_from_vasp_files(
            os.path.join(test_dir_doscar, "POSCAR.Fe3O4"),
            os.path.join(test_dir_doscar, "INCAR.lobster"),
            os.path.join(test_dir_doscar, "POTCAR.Fe3O4"),
            option="standard",
        )
        lobsterin1.write_lobsterin(outfile_path)
        lobsterin2 = Lobsterin.from_file(outfile_path)
        self.assertDictEqual(lobsterin1.diff(lobsterin2)["Different"], {})

    def test_write_INCAR(self):
        # write INCAR and compare
        outfile_path = tempfile.mkstemp()[1]
        lobsterin1 = Lobsterin.standard_calculations_from_vasp_files(
            os.path.join(test_dir_doscar, "POSCAR.Fe3O4"),
            os.path.join(test_dir_doscar, "INCAR.lobster"),
            os.path.join(test_dir_doscar, "POTCAR.Fe3O4"),
            option="standard",
        )
        lobsterin1.write_INCAR(
            os.path.join(test_dir_doscar, "INCAR.lobster3"),
            outfile_path,
            os.path.join(test_dir_doscar, "POSCAR.Fe3O4"),
        )

        incar1 = Incar.from_file(os.path.join(test_dir_doscar, "INCAR.lobster3"))
        incar2 = Incar.from_file(outfile_path)

        self.assertDictEqual(
            incar1.diff(incar2)["Different"],
            {
                "ISYM": {"INCAR1": 2, "INCAR2": -1},
                "NBANDS": {"INCAR1": None, "INCAR2": 86},
                "NSW": {"INCAR1": 500, "INCAR2": 0},
                "LWAVE": {"INCAR1": False, "INCAR2": True},
            },
        )

    def test_write_KPOINTS(self):

        # line mode
        outfile_path = tempfile.mkstemp()[1]
        outfile_path2 = tempfile.mkstemp(prefix="POSCAR")[1]
        lobsterin1 = Lobsterin({})
        # test writing primitive cell
        lobsterin1.write_POSCAR_with_standard_primitive(
            POSCAR_input=os.path.join(test_dir_doscar, "POSCAR.Fe3O4"),
            POSCAR_output=outfile_path2,
        )

        lobsterin1.write_KPOINTS(
            POSCAR_input=outfile_path2,
            KPOINTS_output=outfile_path,
            kpoints_line_density=58,
        )
        kpoint = Kpoints.from_file(outfile_path)
        self.assertEqual(kpoint.num_kpts, 562)
        self.assertAlmostEqual(kpoint.kpts[-1][0], -0.5)
        self.assertAlmostEqual(kpoint.kpts[-1][1], 0.5)
        self.assertAlmostEqual(kpoint.kpts[-1][2], 0.5)
        self.assertEqual(kpoint.labels[-1], "T")
        kpoint2 = Kpoints.from_file(os.path.join(test_dir_doscar, "KPOINTS_band.lobster"))

        labels = []
        number = 0
        for label in kpoint.labels:
            if label is not None:
                if number != 0:
                    if label != labels[number - 1]:
                        labels.append(label)
                        number += 1
                else:
                    labels.append(label)
                    number += 1

        labels2 = []
        number2 = 0
        for label in kpoint2.labels:
            if label is not None:
                if number2 != 0:
                    if label != labels2[number2 - 1]:
                        labels2.append(label)
                        number2 += 1
                else:
                    labels2.append(label)
                    number2 += 1
        self.assertListEqual(labels, labels2)

        # without line mode
        lobsterin1.write_KPOINTS(POSCAR_input=outfile_path2, KPOINTS_output=outfile_path, line_mode=False)
        kpoint = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(os.path.join(test_dir_doscar, "IBZKPT.lobster"))

        for num_kpt, list_kpoint in enumerate(kpoint.kpts):
            self.assertAlmostEqual(list_kpoint[0], kpoint2.kpts[num_kpt][0])
            self.assertAlmostEqual(list_kpoint[1], kpoint2.kpts[num_kpt][1])
            self.assertAlmostEqual(list_kpoint[2], kpoint2.kpts[num_kpt][2])

        self.assertEqual(kpoint.num_kpts, 108)

        # without line mode, use grid instead of reciprocal density
        lobsterin1.write_KPOINTS(
            POSCAR_input=outfile_path2,
            KPOINTS_output=outfile_path,
            line_mode=False,
            from_grid=True,
            input_grid=[6, 6, 3],
        )
        kpoint = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(os.path.join(test_dir_doscar, "IBZKPT.lobster"))

        for num_kpt, list_kpoint in enumerate(kpoint.kpts):
            self.assertAlmostEqual(list_kpoint[0], kpoint2.kpts[num_kpt][0])
            self.assertAlmostEqual(list_kpoint[1], kpoint2.kpts[num_kpt][1])
            self.assertAlmostEqual(list_kpoint[2], kpoint2.kpts[num_kpt][2])

        self.assertEqual(kpoint.num_kpts, 108)

        #
        # #without line mode, using a certain grid, isym=0 instead of -1
        lobsterin1.write_KPOINTS(
            POSCAR_input=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "POSCAR.Li"),
            KPOINTS_output=outfile_path,
            line_mode=False,
            from_grid=True,
            input_grid=[3, 3, 3],
            isym=0,
        )

        kpoint1 = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "IBZKPT_3_3_3_Li"))
        for ikpoint, kpoint in enumerate(kpoint1.kpts):
            self.assertTrue(
                self.is_kpoint_in_list(
                    kpoint,
                    kpoint2.kpts,
                    kpoint1.kpts_weights[ikpoint],
                    kpoint2.kpts_weights,
                )
            )
        for ikpoint, kpoint in enumerate(kpoint2.kpts):
            self.assertTrue(
                self.is_kpoint_in_list(
                    kpoint,
                    kpoint1.kpts,
                    kpoint2.kpts_weights[ikpoint],
                    kpoint1.kpts_weights,
                )
            )

        lobsterin1.write_KPOINTS(
            POSCAR_input=os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "POSCAR.Li"),
            KPOINTS_output=outfile_path,
            line_mode=False,
            from_grid=True,
            input_grid=[2, 2, 2],
            isym=0,
        )

        kpoint1 = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "IBZKPT_2_2_2_Li"))
        for ikpoint, kpoint in enumerate(kpoint1.kpts):
            self.assertTrue(
                self.is_kpoint_in_list(
                    kpoint,
                    kpoint2.kpts,
                    kpoint1.kpts_weights[ikpoint],
                    kpoint2.kpts_weights,
                )
            )
        for ikpoint, kpoint in enumerate(kpoint2.kpts):
            self.assertTrue(
                self.is_kpoint_in_list(
                    kpoint,
                    kpoint1.kpts,
                    kpoint2.kpts_weights[ikpoint],
                    kpoint1.kpts_weights,
                )
            )

    def is_kpoint_in_list(self, kpoint, kpointlist, weight, weightlist):
        found = 0
        for ikpoint2, kpoint2 in enumerate(kpointlist):
            if (
                np.isclose(kpoint[0], kpoint2[0])
                and np.isclose(kpoint[1], kpoint2[1])
                and np.isclose(kpoint[2], kpoint2[2])
            ):
                if weight == weightlist[ikpoint2]:
                    found += 1
            elif (
                np.isclose(-kpoint[0], kpoint2[0])
                and np.isclose(-kpoint[1], kpoint2[1])
                and np.isclose(-kpoint[2], kpoint2[2])
            ):
                if weight == weightlist[ikpoint2]:
                    found += 1
        if found == 1:
            return True
        else:
            return False

    def test_MSONable_implementation(self):
        # tests as dict and from dict methods
        newLobsterin = Lobsterin.from_dict(self.Lobsterinfromfile.as_dict())
        self.assertDictEqual(newLobsterin, self.Lobsterinfromfile)
        newLobsterin.to_json()

    def tearDown(self):
        warnings.simplefilter("default")


class BandoverlapsTest(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")
        # test spin polarlized calc and non spinpolarized calc

        self.bandoverlaps1 = Bandoverlaps(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "bandOverlaps.lobster.1"))
        self.bandoverlaps2 = Bandoverlaps(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "bandOverlaps.lobster.2"))

    def test_attributes(self):
        # bandoverlapsdict
        self.assertAlmostEqual(
            self.bandoverlaps1.bandoverlapsdict[Spin.up]["0.5 0 0"]["maxDeviation"],
            0.000278953,
        )
        self.assertAlmostEqual(
            self.bandoverlaps1.bandoverlapsdict[Spin.up]["0.5 0 0"]["matrix"][-1][-1],
            0.0188058,
        )
        self.assertAlmostEqual(self.bandoverlaps1.bandoverlapsdict[Spin.up]["0.5 0 0"]["matrix"][0][0], 1)

        self.assertAlmostEqual(
            self.bandoverlaps1.bandoverlapsdict[Spin.down]["0.0261194 0.0261194 0.473881"]["maxDeviation"],
            4.31567e-05,
        )
        self.assertAlmostEqual(
            self.bandoverlaps1.bandoverlapsdict[Spin.down]["0.0261194 0.0261194 0.473881"]["matrix"][0][-1],
            4.0066e-07,
        )

        # maxDeviation
        self.assertAlmostEqual(self.bandoverlaps1.max_deviation[0], 0.000278953)
        self.assertAlmostEqual(self.bandoverlaps1.max_deviation[-1], 4.31567e-05)

        self.assertAlmostEqual(self.bandoverlaps2.max_deviation[0], 0.000473319)
        self.assertAlmostEqual(self.bandoverlaps2.max_deviation[-1], 1.48451e-05)

    def test_has_good_quality(self):
        self.assertFalse(self.bandoverlaps1.has_good_quality_maxDeviation(limit_maxDeviation=0.1))
        self.assertFalse(
            self.bandoverlaps1.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=9,
                number_occ_bands_spin_down=5,
                limit_deviation=0.1,
                spin_polarized=True,
            )
        )
        self.assertTrue(
            self.bandoverlaps1.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=3,
                number_occ_bands_spin_down=0,
                limit_deviation=0.001,
                spin_polarized=True,
            )
        )
        self.assertFalse(
            self.bandoverlaps1.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=1,
                number_occ_bands_spin_down=1,
                limit_deviation=0.000001,
                spin_polarized=True,
            )
        )
        self.assertFalse(
            self.bandoverlaps1.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=1,
                number_occ_bands_spin_down=0,
                limit_deviation=0.000001,
                spin_polarized=True,
            )
        )
        self.assertFalse(
            self.bandoverlaps1.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=0,
                number_occ_bands_spin_down=1,
                limit_deviation=0.000001,
                spin_polarized=True,
            )
        )
        self.assertFalse(
            self.bandoverlaps1.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=4,
                number_occ_bands_spin_down=4,
                limit_deviation=0.001,
                spin_polarized=True,
            )
        )

        self.assertTrue(self.bandoverlaps1.has_good_quality_maxDeviation(limit_maxDeviation=100))
        self.assertTrue(self.bandoverlaps2.has_good_quality_maxDeviation())
        self.assertFalse(self.bandoverlaps2.has_good_quality_maxDeviation(limit_maxDeviation=0.0000001))
        self.assertFalse(
            self.bandoverlaps2.has_good_quality_check_occupied_bands(
                number_occ_bands_spin_up=10, limit_deviation=0.0000001
            )
        )
        self.assertTrue(
            self.bandoverlaps2.has_good_quality_check_occupied_bands(number_occ_bands_spin_up=1, limit_deviation=0.1)
        )

        self.assertFalse(
            self.bandoverlaps2.has_good_quality_check_occupied_bands(number_occ_bands_spin_up=1, limit_deviation=1e-8)
        )
        self.assertTrue(
            self.bandoverlaps2.has_good_quality_check_occupied_bands(number_occ_bands_spin_up=10, limit_deviation=0.1)
        )

        self.assertTrue(
            self.bandoverlaps2.has_good_quality_check_occupied_bands(number_occ_bands_spin_up=1, limit_deviation=0.1)
        )


class GrosspopTest(unittest.TestCase):
    def setUp(self):
        self.grosspop1 = Grosspop(os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "GROSSPOP.lobster"))

    def testattributes(self):
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[0]["Mulliken GP"]["3s"], 0.52)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[0]["Mulliken GP"]["3p_y"], 0.38)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[0]["Mulliken GP"]["3p_z"], 0.37)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[0]["Mulliken GP"]["3p_x"], 0.37)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[0]["Mulliken GP"]["total"], 1.64)
        self.assertEqual(self.grosspop1.list_dict_grosspop[0]["element"], "Si")
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[0]["Loewdin GP"]["3s"], 0.61)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[0]["Loewdin GP"]["3p_y"], 0.52)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[0]["Loewdin GP"]["3p_z"], 0.52)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[0]["Loewdin GP"]["3p_x"], 0.52)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[0]["Loewdin GP"]["total"], 2.16)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[5]["Mulliken GP"]["2s"], 1.80)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[5]["Loewdin GP"]["2s"], 1.60)
        self.assertEqual(self.grosspop1.list_dict_grosspop[5]["element"], "O")
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[8]["Mulliken GP"]["2s"], 1.80)
        self.assertAlmostEqual(self.grosspop1.list_dict_grosspop[8]["Loewdin GP"]["2s"], 1.60)
        self.assertEqual(self.grosspop1.list_dict_grosspop[8]["element"], "O")

    def test_structure_with_grosspop(self):
        struct_dict = {
            "@module": "pymatgen.core.structure",
            "@class": "Structure",
            "charge": None,
            "lattice": {
                "matrix": [
                    [5.021897888834907, 4.53806e-11, 0.0],
                    [-2.5109484443388332, 4.349090983701526, 0.0],
                    [0.0, 0.0, 5.511929408565514],
                ],
                "a": 5.021897888834907,
                "b": 5.0218974974248045,
                "c": 5.511929408565514,
                "alpha": 90.0,
                "beta": 90.0,
                "gamma": 119.99999598960493,
                "volume": 120.38434608659402,
            },
            "sites": [
                {
                    "species": [{"element": "Si", "occu": 1}],
                    "abc": [-3e-16, 0.4763431475490085, 0.6666669999999968],
                    "xyz": [-1.1960730853096477, 2.0716596881533986, 3.674621443020128],
                    "label": "Si",
                    "properties": {"Total Mulliken GP": 1.64, "Total Loewdin GP": 2.16},
                },
                {
                    "species": [{"element": "Si", "occu": 1}],
                    "abc": [0.5236568524509936, 0.5236568524509926, 0.0],
                    "xyz": [1.3148758827683875, 2.277431295571896, 0.0],
                    "label": "Si",
                    "properties": {"Total Mulliken GP": 1.64, "Total Loewdin GP": 2.16},
                },
                {
                    "species": [{"element": "Si", "occu": 1}],
                    "abc": [0.4763431475490066, -1.2e-15, 0.3333330000000032],
                    "xyz": [
                        2.392146647037334,
                        2.1611518932482004e-11,
                        1.8373079655453863,
                    ],
                    "label": "Si",
                    "properties": {"Total Mulliken GP": 1.64, "Total Loewdin GP": 2.16},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.1589037798059321, 0.7440031622164922, 0.4613477252144715],
                    "xyz": [-1.0701550264153763, 3.235737444648381, 2.5429160941844473],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.2559968377835071, 0.4149006175894398, 0.7946807252144676],
                    "xyz": [0.2437959189219816, 1.8044405351020447, 4.380224059729795],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.5850993824105679, 0.8410962201940679, 0.1280147252144683],
                    "xyz": [0.8263601076506712, 3.6580039876980064, 0.7056081286390611],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.7440031622164928, 0.1589037798059326, 0.5386522747855285],
                    "xyz": [3.337308710918233, 0.6910869960638374, 2.969013314381067],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.4149006175894392, 0.2559968377835, 0.2053192747855324],
                    "xyz": [1.4407936739605638, 1.1133535390791505, 1.13170534883572],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
                {
                    "species": [{"element": "O", "occu": 1}],
                    "abc": [0.841096220194068, 0.5850993824105675, 0.8719852747855317],
                    "xyz": [2.754744948452184, 2.5446504486493, 4.806321279926453],
                    "label": "O",
                    "properties": {"Total Mulliken GP": 7.18, "Total Loewdin GP": 6.92},
                },
            ],
        }

        newstructure = self.grosspop1.get_structure_with_total_grosspop(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "POSCAR.SiO2")
        )
        for coords, coords2 in zip(newstructure.frac_coords, Structure.from_dict(struct_dict).frac_coords):
            for xyz, xyz2 in zip(coords, coords2):
                self.assertAlmostEqual(xyz, xyz2)


class TestUtils(PymatgenTest):
    def test_get_all_possible_basis_combinations(self):
        # this basis is just for testing (not correct)
        min_basis = ["Li 1s 2s ", "Na 1s 2s", "Si 1s 2s"]
        max_basis = ["Li 1s 2p 2s ", "Na 1s 2p 2s", "Si 1s 2s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        self.assertListEqual(
            combinations_basis,
            [
                ["Li 1s 2s", "Na 1s 2s", "Si 1s 2s"],
                ["Li 1s 2s", "Na 1s 2s 2p", "Si 1s 2s"],
                ["Li 1s 2s 2p", "Na 1s 2s", "Si 1s 2s"],
                ["Li 1s 2s 2p", "Na 1s 2s 2p", "Si 1s 2s"],
            ],
        )

        min_basis = ["Li 1s 2s"]
        max_basis = ["Li 1s 2s 2p 3s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        self.assertListEqual(
            combinations_basis,
            [["Li 1s 2s"], ["Li 1s 2s 2p"], ["Li 1s 2s 3s"], ["Li 1s 2s 2p 3s"]],
        )

        min_basis = ["Li 1s 2s", "Na 1s 2s"]
        max_basis = ["Li 1s 2s 2p 3s", "Na 1s 2s 2p 3s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        self.assertListEqual(
            combinations_basis,
            [
                ["Li 1s 2s", "Na 1s 2s"],
                ["Li 1s 2s", "Na 1s 2s 2p"],
                ["Li 1s 2s", "Na 1s 2s 3s"],
                ["Li 1s 2s", "Na 1s 2s 2p 3s"],
                ["Li 1s 2s 2p", "Na 1s 2s"],
                ["Li 1s 2s 2p", "Na 1s 2s 2p"],
                ["Li 1s 2s 2p", "Na 1s 2s 3s"],
                ["Li 1s 2s 2p", "Na 1s 2s 2p 3s"],
                ["Li 1s 2s 3s", "Na 1s 2s"],
                ["Li 1s 2s 3s", "Na 1s 2s 2p"],
                ["Li 1s 2s 3s", "Na 1s 2s 3s"],
                ["Li 1s 2s 3s", "Na 1s 2s 2p 3s"],
                ["Li 1s 2s 2p 3s", "Na 1s 2s"],
                ["Li 1s 2s 2p 3s", "Na 1s 2s 2p"],
                ["Li 1s 2s 2p 3s", "Na 1s 2s 3s"],
                ["Li 1s 2s 2p 3s", "Na 1s 2s 2p 3s"],
            ],
        )

        min_basis = ["Si 1s 2s 2p", "Na 1s 2s"]
        max_basis = ["Si 1s 2s 2p 3s", "Na 1s 2s 2p 3s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        self.assertListEqual(
            combinations_basis,
            [
                ["Si 1s 2s 2p", "Na 1s 2s"],
                ["Si 1s 2s 2p", "Na 1s 2s 2p"],
                ["Si 1s 2s 2p", "Na 1s 2s 3s"],
                ["Si 1s 2s 2p", "Na 1s 2s 2p 3s"],
                ["Si 1s 2s 2p 3s", "Na 1s 2s"],
                ["Si 1s 2s 2p 3s", "Na 1s 2s 2p"],
                ["Si 1s 2s 2p 3s", "Na 1s 2s 3s"],
                ["Si 1s 2s 2p 3s", "Na 1s 2s 2p 3s"],
            ],
        )


class WavefunctionTest(PymatgenTest):
    def test_parse_file(self):
        grid, points, real, imaginary, distance = Wavefunction._parse_file(
            os.path.join(
                test_dir_doscar,
                "cohp",
                "LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            )
        )
        self.assertArrayEqual([41, 41, 41], grid)
        self.assertAlmostEqual(points[4][0], 0.0000)
        self.assertAlmostEqual(points[4][1], 0.0000)
        self.assertAlmostEqual(points[4][2], 0.4000)
        self.assertAlmostEqual(real[8], 1.38863e-01)
        self.assertAlmostEqual(imaginary[8], 2.89645e-01)
        self.assertEqual(len(imaginary), 41 * 41 * 41)
        self.assertEqual(len(real), 41 * 41 * 41)
        self.assertEqual(len(points), 41 * 41 * 41)
        self.assertAlmostEqual(distance[0], 0.0000)

    def test_set_volumetric_data(self):
        wave1 = Wavefunction(
            filename=os.path.join(
                test_dir_doscar,
                "cohp",
                "LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1" ".gz",
            ),
            structure=Structure.from_file(os.path.join(test_dir_doscar, "cohp", "POSCAR_O.gz")),
        )

        wave1.set_volumetric_data(grid=wave1.grid, structure=wave1.structure)
        self.assertTrue(hasattr(wave1, "volumetricdata_real"))
        self.assertTrue(hasattr(wave1, "volumetricdata_imaginary"))

    def test_get_volumetricdata_real(self):
        wave1 = Wavefunction(
            filename=os.path.join(
                test_dir_doscar,
                "cohp",
                "LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            ),
            structure=Structure.from_file(os.path.join(test_dir_doscar, "cohp", "POSCAR_O.gz")),
        )
        volumetricdata_real = wave1.get_volumetricdata_real()
        self.assertAlmostEqual(volumetricdata_real.data["total"][0, 0, 0], -3.0966)

    def test_get_volumetricdata_imaginary(self):
        wave1 = Wavefunction(
            filename=os.path.join(
                test_dir_doscar,
                "cohp",
                "LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            ),
            structure=Structure.from_file(os.path.join(test_dir_doscar, "cohp", "POSCAR_O.gz")),
        )
        volumetricdata_imaginary = wave1.get_volumetricdata_imaginary()
        self.assertAlmostEqual(volumetricdata_imaginary.data["total"][0, 0, 0], -6.45895e00)

    def test_get_volumetricdata_density(self):
        wave1 = Wavefunction(
            filename=os.path.join(
                test_dir_doscar,
                "cohp",
                "LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            ),
            structure=Structure.from_file(os.path.join(test_dir_doscar, "cohp", "POSCAR_O.gz")),
        )
        volumetricdata_density = wave1.get_volumetricdata_density()
        self.assertAlmostEqual(
            volumetricdata_density.data["total"][0, 0, 0],
            (-3.0966 * -3.0966) + (-6.45895 * -6.45895),
        )

    def test_write_file(self):
        wave1 = Wavefunction(
            filename=os.path.join(
                test_dir_doscar,
                "cohp",
                "LCAOWaveFunctionAfterLSO1PlotOfSpin1Kpoint1band1.gz",
            ),
            structure=Structure.from_file(os.path.join(test_dir_doscar, "cohp", "POSCAR_O.gz")),
        )
        wave1.write_file(filename=os.path.join("wavecar_test.vasp"), part="real")
        self.assertTrue(os.path.isfile("wavecar_test.vasp"))

        wave1.write_file(filename=os.path.join("wavecar_test.vasp"), part="imaginary")
        self.assertTrue(os.path.isfile("wavecar_test.vasp"))
        os.remove("wavecar_test.vasp")
        wave1.write_file(filename=os.path.join("density.vasp"), part="density")
        self.assertTrue(os.path.isfile("density.vasp"))
        os.remove("density.vasp")

    def tearDown(self):
        warnings.simplefilter("default")


if __name__ == "__main__":
    unittest.main()
