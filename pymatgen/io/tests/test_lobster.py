# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest
import os
import json
import numpy as np
from pymatgen import Structure
from pymatgen.io.lobster import Cohpcar, Icohplist, Doscar, Charge
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.util.testing import PymatgenTest

__author__ = "Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.2"
__email__ = "esters@uoregon.edu"
__date__ = "Dec 10, 2017"

test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "..", "test_files", "cohp")
test_dir_doscar = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "..", "..", "..", "test_files")

this_dir = os.path.dirname(os.path.abspath(__file__))


class CohpcarTest(PymatgenTest):
    def setUp(self):
        self.cohp_bise = Cohpcar(filename=os.path.join(test_dir, "COHPCAR.lobster.BiSe"))
        self.coop_bise = Cohpcar(filename=os.path.join(test_dir, "COOPCAR.lobster.BiSe"),
                                 are_coops=True)
        self.cohp_fe = Cohpcar(filename=os.path.join(test_dir, "COOPCAR.lobster"))
        self.coop_fe = Cohpcar(filename=os.path.join(test_dir, "COOPCAR.lobster"), are_coops=True)
        self.orb = Cohpcar(filename=os.path.join(test_dir, "COHPCAR.lobster.orbitalwise"))
        self.orb_notot = Cohpcar(filename=os.path.join(test_dir, "COHPCAR.lobster.notot.orbitalwise"))

        # Lobster 3.1 (Test data is from prerelease of Lobster 3.1)
        self.cohp_KF = Cohpcar(filename=os.path.join(test_dir, "COHPCAR.lobster.KF"))
        self.coop_KF = Cohpcar(filename=os.path.join(test_dir, "COHPCAR.lobster.KF"), are_coops=True)

        # example with f electrons
        self.cohp_Na2UO4 = Cohpcar(filename=os.path.join(test_dir, "COHPCAR.lobster.Na2UO4"))
        self.coop_Na2UO4 = Cohpcar(filename=os.path.join(test_dir, "COOPCAR.lobster.Na2UO4"), are_coops=True)

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

        self.assertAlmostEqual(self.cohp_bise.energies[0]
                               + self.cohp_bise.efermi, elim_bise[0],
                               places=4)
        self.assertAlmostEqual(self.cohp_bise.energies[-1]
                               + self.cohp_bise.efermi, elim_bise[1],
                               places=4)
        self.assertAlmostEqual(self.coop_bise.energies[0]
                               + self.coop_bise.efermi, elim_bise[0],
                               places=4)
        self.assertAlmostEqual(self.coop_bise.energies[-1]
                               + self.coop_bise.efermi, elim_bise[1],
                               places=4)

        self.assertAlmostEqual(self.cohp_fe.energies[0] + self.cohp_fe.efermi,
                               elim_fe[0], places=4)
        self.assertAlmostEqual(self.cohp_fe.energies[-1] + self.cohp_fe.efermi,
                               elim_fe[1], places=4)
        self.assertAlmostEqual(self.coop_fe.energies[0] + self.coop_fe.efermi,
                               elim_fe[0], places=4)
        self.assertAlmostEqual(self.coop_fe.energies[-1] + self.coop_fe.efermi,
                               elim_fe[1], places=4)

        # Lobster 3.1
        self.assertAlmostEqual(self.cohp_KF.energies[0] + self.cohp_KF.efermi,
                               elim_KF[0], places=4)
        self.assertAlmostEqual(self.cohp_KF.energies[-1] + self.cohp_KF.efermi,
                               elim_KF[1], places=4)
        self.assertAlmostEqual(self.coop_KF.energies[0] + self.coop_KF.efermi,
                               elim_KF[0], places=4)
        self.assertAlmostEqual(self.coop_KF.energies[-1] + self.coop_KF.efermi,
                               elim_KF[1], places=4)

    def test_cohp_data(self):
        lengths_sites_bise = {"1": (2.882308829886294, (0, 6)),
                              "2": (3.1014396233274444, (0, 9)),
                              "3": (2.8823088298862083, (1, 7)),
                              "4": (3.1014396233275434, (1, 8)),
                              "5": (3.0500070394403904, (2, 9)),
                              "6": (2.9167594580335807, (2, 10)),
                              "7": (3.05000703944039, (3, 8)),
                              "8": (2.9167594580335803, (3, 11)),
                              "9": (3.3752173204052101, (4, 11)),
                              "10": (3.0729354518345948, (4, 5)),
                              "11": (3.3752173204052101, (5, 10))}
        lengths_sites_fe = {"1": (2.8318907764979082, (7, 6)),
                            "2": (2.4524893531900283, (7, 8))}
        # Lobster 3.1
        lengths_sites_KF = {"1": (2.7119923200622269, (0, 1)),
                            "2": (2.7119923200622269, (0, 1)),
                            "3": (2.7119923576010501, (0, 1)),
                            "4": (2.7119923576010501, (0, 1)),
                            "5": (2.7119923200622269, (0, 1)),
                            "6": (2.7119923200622269, (0, 1))}

        for data in [self.cohp_bise.cohp_data, self.coop_bise.cohp_data]:
            for bond in data:
                if bond != "average":
                    self.assertEqual(data[bond]["length"],
                                     lengths_sites_bise[bond][0])
                    self.assertEqual(data[bond]["sites"],
                                     lengths_sites_bise[bond][1])
                    self.assertEqual(len(data[bond]["COHP"][Spin.up]), 241)
                    self.assertEqual(len(data[bond]["ICOHP"][Spin.up]), 241)
        for data in [self.cohp_fe.cohp_data, self.coop_fe.cohp_data]:
            for bond in data:
                if bond != "average":
                    self.assertEqual(data[bond]["length"],
                                     lengths_sites_fe[bond][0])
                    self.assertEqual(data[bond]["sites"],
                                     lengths_sites_fe[bond][1])
                    self.assertEqual(len(data[bond]["COHP"][Spin.up]), 301)
                    self.assertEqual(len(data[bond]["ICOHP"][Spin.up]), 301)

        # Lobster 3.1
        for data in [self.cohp_KF.cohp_data, self.coop_KF.cohp_data]:
            for bond in data:
                if bond != "average":
                    self.assertEqual(data[bond]["length"],
                                     lengths_sites_KF[bond][0])
                    self.assertEqual(data[bond]["sites"],
                                     lengths_sites_KF[bond][1])
                    self.assertEqual(len(data[bond]["COHP"][Spin.up]), 6)
                    self.assertEqual(len(data[bond]["ICOHP"][Spin.up]), 6)

    def test_orbital_resolved_cohp(self):
        orbitals = [tuple((Orbital(i), Orbital(j))) for j in range(4)
                    for i in range(4)]
        self.assertIsNone(self.cohp_bise.orb_res_cohp)
        self.assertIsNone(self.coop_bise.orb_res_cohp)
        self.assertIsNone(self.cohp_fe.orb_res_cohp)
        self.assertIsNone(self.coop_fe.orb_res_cohp)
        self.assertIsNone(self.orb_notot.cohp_data["1"]["COHP"])
        self.assertIsNone(self.orb_notot.cohp_data["1"]["ICOHP"])
        for orbs in self.orb.orb_res_cohp["1"]:
            orb_set = self.orb.orb_res_cohp["1"][orbs]["orbitals"]
            # print(orb_set[0][0])
            self.assertEqual(orb_set[0][0], 4)
            self.assertEqual(orb_set[1][0], 4)
            self.assertIn(tuple((orb_set[0][1], orb_set[1][1])), orbitals)

        # test d and f orbitals
        comparelist = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6,
                       6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                       7, 7, 7, 7]
        comparelist2 = ["f0", "f0", "f0", "f0", "f1", "f1", "f1", "f1", "f2", "f2", "f2", "f2", "f3", "f3", "f3", "f3",
                        "f_1", "f_1", "f_1", "f_1", "f_2", "f_2", "f_2", "f_2", "f_3", "f_3", "f_3", "f_3", "dx2",
                        "dx2", "dx2", "dx2", "dxy", "dxy", "dxy", "dxy", "dxz", "dxz", "dxz", "dxz", "dyz", "dyz",
                        "dyz", "dyz", "dz2", "dz2", "dz2", "dz2", "px", "px", "px", "px", "py", "py", "py", "py", "pz",
                        "pz", "pz", "pz", "s", "s", "s", "s", "s", "s", "s", "s"]
        for iorb, orbs in enumerate(sorted(self.cohp_Na2UO4.orb_res_cohp["49"])):
            orb_set = self.cohp_Na2UO4.orb_res_cohp["49"][orbs]["orbitals"]
            self.assertEqual(orb_set[0][0], comparelist[iorb])
            self.assertEqual(str(orb_set[0][1]), comparelist2[iorb])

        # The sum of the orbital-resolved COHPs should be approximately
        # the total COHP. Due to small deviations in the LOBSTER calculation,
        # the precision is not very high though.
        cohp = self.orb.cohp_data["1"]["COHP"][Spin.up]
        icohp = self.orb.cohp_data["1"]["ICOHP"][Spin.up]
        tot = np.sum([self.orb.orb_res_cohp["1"][orbs]["COHP"][Spin.up]
                      for orbs in self.orb.orb_res_cohp["1"]], axis=0)
        self.assertArrayAlmostEqual(tot, cohp, decimal=3)
        tot = np.sum([self.orb.orb_res_cohp["1"][orbs]["ICOHP"][Spin.up]
                      for orbs in self.orb.orb_res_cohp["1"]], axis=0)
        self.assertArrayAlmostEqual(tot, icohp, decimal=3)

        # Lobster 3.1
        cohp_KF = self.cohp_KF.cohp_data["1"]["COHP"][Spin.up]
        icohp_KF = self.cohp_KF.cohp_data["1"]["ICOHP"][Spin.up]
        tot_KF = np.sum([self.cohp_KF.orb_res_cohp["1"][orbs]["COHP"][Spin.up]
                         for orbs in self.cohp_KF.orb_res_cohp["1"]], axis=0)
        self.assertArrayAlmostEqual(tot_KF, cohp_KF, decimal=3)
        tot_KF = np.sum([self.cohp_KF.orb_res_cohp["1"][orbs]["ICOHP"][Spin.up]
                         for orbs in self.cohp_KF.orb_res_cohp["1"]], axis=0)
        self.assertArrayAlmostEqual(tot_KF, icohp_KF, decimal=3)

        # d and f orbitals
        cohp_Na2UO4 = self.cohp_Na2UO4.cohp_data["49"]["COHP"][Spin.up]
        icohp_Na2UO4 = self.cohp_Na2UO4.cohp_data["49"]["ICOHP"][Spin.up]
        tot_Na2UO4 = np.sum([self.cohp_Na2UO4.orb_res_cohp["49"][orbs]["COHP"][Spin.up]
                             for orbs in self.cohp_Na2UO4.orb_res_cohp["49"]], axis=0)
        self.assertArrayAlmostEqual(tot_Na2UO4, cohp_Na2UO4, decimal=3)
        tot_Na2UO4 = np.sum([self.cohp_Na2UO4.orb_res_cohp["49"][orbs]["ICOHP"][Spin.up]
                             for orbs in self.cohp_Na2UO4.orb_res_cohp["49"]], axis=0)
        self.assertArrayAlmostEqual(tot_Na2UO4, icohp_Na2UO4, decimal=3)


class IcohplistTest(unittest.TestCase):
    def setUp(self):
        self.icohp_bise = Icohplist(filename=os.path.join(test_dir, "ICOHPLIST.lobster.BiSe"))
        self.icoop_bise = Icohplist(filename=os.path.join(test_dir, "ICOOPLIST.lobster.BiSe"),
                                    are_coops=True)
        self.icohp_fe = Icohplist(filename=os.path.join(test_dir, "ICOHPLIST.lobster"))
        self.icoop_fe = Icohplist(filename=os.path.join(test_dir, "ICOHPLIST.lobster"), are_coops=True)

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
        icohplist_bise = {"1": {"length": 2.88231, "number_of_bonds": 3,
                                "icohp": {Spin.up: -2.18042},
                                "translation": [0, 0, 0]},
                          "2": {"length": 3.10144, "number_of_bonds": 3,
                                "icohp": {Spin.up: -1.14347},
                                "translation": [0, 0, 0]},
                          "3": {"length": 2.88231, "number_of_bonds": 3,
                                "icohp": {Spin.up: -2.18042},
                                "translation": [0, 0, 0]},
                          "4": {"length": 3.10144, "number_of_bonds": 3,
                                "icohp": {Spin.up: -1.14348},
                                "translation": [0, 0, 0]},
                          "5": {"length": 3.05001, "number_of_bonds": 3,
                                "icohp": {Spin.up: -1.30006},
                                "translation": [0, 0, 0]},
                          "6": {"length": 2.91676, "number_of_bonds": 3,
                                "icohp": {Spin.up: -1.96843},
                                "translation": [0, 0, 0]},
                          "7": {"length": 3.05001, "number_of_bonds": 3,
                                "icohp": {Spin.up: -1.30006}, "translation": [0, 0, 0]},
                          "8": {"length": 2.91676, "number_of_bonds": 3,
                                "icohp": {Spin.up: -1.96843}, "translation": [0, 0, 0]},
                          "9": {"length": 3.37522, "number_of_bonds": 3,
                                "icohp": {Spin.up: -0.47531}, "translation": [0, 0, 0]},
                          "10": {"length": 3.07294, "number_of_bonds": 3,
                                 "icohp": {Spin.up: -2.38796}, "translation": [0, 0, 0]},
                          "11": {"length": 3.37522, "number_of_bonds": 3,
                                 "icohp": {Spin.up: -0.47531}, "translation": [0, 0, 0]}, }
        icooplist_bise = {"1": {"length": 2.88231, "number_of_bonds": 3,
                                "icohp": {Spin.up: 0.14245}, "translation": [0, 0, 0]},
                          "2": {"length": 3.10144, "number_of_bonds": 3,
                                "icohp": {Spin.up: -0.04118}, "translation": [0, 0, 0]},
                          "3": {"length": 2.88231, "number_of_bonds": 3,
                                "icohp": {Spin.up: 0.14245}, "translation": [0, 0, 0]},
                          "4": {"length": 3.10144, "number_of_bonds": 3,
                                "icohp": {Spin.up: -0.04118}, "translation": [0, 0, 0]},
                          "5": {"length": 3.05001, "number_of_bonds": 3,
                                "icohp": {Spin.up: -0.03516}, "translation": [0, 0, 0]},
                          "6": {"length": 2.91676, "number_of_bonds": 3,
                                "icohp": {Spin.up: 0.10745}, "translation": [0, 0, 0]},
                          "7": {"length": 3.05001, "number_of_bonds": 3,
                                "icohp": {Spin.up: -0.03516}, "translation": [0, 0, 0]},
                          "8": {"length": 2.91676, "number_of_bonds": 3,
                                "icohp": {Spin.up: 0.10745}, "translation": [0, 0, 0]},
                          "9": {"length": 3.37522, "number_of_bonds": 3,
                                "icohp": {Spin.up: -0.12395}, "translation": [0, 0, 0]},
                          "10": {"length": 3.07294, "number_of_bonds": 3,
                                 "icohp": {Spin.up: 0.24714}, "translation": [0, 0, 0]},
                          "11": {"length": 3.37522, "number_of_bonds": 3,
                                 "icohp": {Spin.up: -0.12395}, "translation": [0, 0, 0]}}
        icooplist_fe = {'1': {'length': 2.83189, 'number_of_bonds': 2,
                              'icohp': {Spin.up: -0.10218, Spin.down: -0.19701},
                              'translation': [0, 0, 0]},
                        '2': {'length': 2.45249, 'number_of_bonds': 1,
                              'icohp': {Spin.up: -0.28485, Spin.down: -0.58279},
                              'translation': [0, 0, 0]}}

        self.assertEqual(icohplist_bise, self.icohp_bise.icohplist)
        self.assertEqual(icooplist_fe, self.icoop_fe.icohplist)


class DoscarTest(unittest.TestCase):
    def setUp(self):
        # first for spin polarized version
        doscar = os.path.join(test_dir_doscar, "DOSCAR.lobster.spin")
        vasprun = os.path.join(test_dir_doscar, "vasprun.xml.lobster.spin")
        doscar2 = os.path.join(test_dir_doscar, "DOSCAR.lobster.nonspin")
        vasprun2 = os.path.join(test_dir_doscar, "vasprun.xml.lobster.nonspin")

        self.DOSCAR_spin_pol = Doscar(doscar=doscar, vasprun=vasprun)
        self.DOSCAR_nonspin_pol = Doscar(doscar=doscar2, vasprun=vasprun2)

        with open(os.path.join(test_dir_doscar, 'structure_KF.json'), 'r') as f:
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
        self.assertDictEqual(self.DOSCAR_spin_pol.completedos.structure.as_dict(), self.structure.as_dict())
        self.assertListEqual(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]['2s'][Spin.up].tolist(),
                             PDOS_F_2s_up)
        self.assertListEqual(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]['2s'][Spin.down].tolist(),
                             PDOS_F_2s_down)
        self.assertListEqual(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]['2p_y'][Spin.up].tolist(),
                             PDOS_F_2py_up)
        self.assertListEqual(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]['2p_y'][Spin.down].tolist(),
                             PDOS_F_2py_down)
        self.assertListEqual(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]['2p_z'][Spin.up].tolist(),
                             PDOS_F_2pz_up)
        self.assertListEqual(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]['2p_z'][Spin.down].tolist(),
                             PDOS_F_2pz_down)
        self.assertListEqual(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]['2p_x'][Spin.up].tolist(),
                             PDOS_F_2px_up)
        self.assertListEqual(self.DOSCAR_spin_pol.completedos.pdos[self.structure[0]]['2p_x'][Spin.down].tolist(),
                             PDOS_F_2px_down)

        energies_nonspin = [-11.25000, -7.50000, -3.75000, 0.00000, 3.75000, 7.50000]
        tdos_nonspin = [0.00000, 1.60000, 0.00000, 1.60000, 0.00000, 0.02418]
        PDOS_F_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        PDOS_F_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        PDOS_F_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        PDOS_F_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        self.assertListEqual(energies_nonspin, self.DOSCAR_nonspin_pol.completedos.energies.tolist())

        self.assertListEqual(tdos_nonspin, self.DOSCAR_nonspin_pol.completedos.densities[Spin.up].tolist())

        self.assertAlmostEqual(fermi, self.DOSCAR_nonspin_pol.completedos.efermi)
        self.assertDictEqual(self.DOSCAR_nonspin_pol.completedos.structure.as_dict(), self.structure.as_dict())

        self.assertListEqual(self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]['2s'][Spin.up].tolist(),
                             PDOS_F_2s)
        self.assertListEqual(self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]['2p_y'][Spin.up].tolist(),
                             PDOS_F_2py)
        self.assertListEqual(self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]['2p_z'][Spin.up].tolist(),
                             PDOS_F_2pz)
        self.assertListEqual(self.DOSCAR_nonspin_pol.completedos.pdos[self.structure[0]]['2p_x'][Spin.up].tolist(),
                             PDOS_F_2px)

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

        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]['2s'][Spin.up].tolist(),
                             PDOS_F_2s_up)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]['2s'][Spin.down].tolist(),
                             PDOS_F_2s_down)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]['2p_y'][Spin.up].tolist(),
                             PDOS_F_2py_up)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]['2p_y'][Spin.down].tolist(),
                             PDOS_F_2py_down)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]['2p_z'][Spin.up].tolist(),
                             PDOS_F_2pz_up)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]['2p_z'][Spin.down].tolist(),
                             PDOS_F_2pz_down)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]['2p_x'][Spin.up].tolist(),
                             PDOS_F_2px_up)
        self.assertListEqual(self.DOSCAR_spin_pol.pdos[0]['2p_x'][Spin.down].tolist(),
                             PDOS_F_2px_down)

        # non spin
        PDOS_F_2s = [0.00000, 0.00320, 0.00000, 0.00017, 0.00000, 0.00060]
        PDOS_F_2py = [0.00000, 0.00322, 0.00000, 0.51635, 0.00000, 0.00037]
        PDOS_F_2pz = [0.00000, 0.00322, 0.00000, 0.51636, 0.00000, 0.00037]
        PDOS_F_2px = [0.00000, 0.00322, 0.00000, 0.51634, 0.00000, 0.00037]

        self.assertListEqual(self.DOSCAR_nonspin_pol.pdos[0]['2s'][Spin.up].tolist(), PDOS_F_2s)
        self.assertListEqual(self.DOSCAR_nonspin_pol.pdos[0]['2p_y'][Spin.up].tolist(), PDOS_F_2py)
        self.assertListEqual(self.DOSCAR_nonspin_pol.pdos[0]['2p_z'][Spin.up].tolist(), PDOS_F_2pz)
        self.assertListEqual(self.DOSCAR_nonspin_pol.pdos[0]['2p_x'][Spin.up].tolist(), PDOS_F_2px)

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

    def test_is_spin_polarized(self):
        # first for spin polarized version
        self.assertTrue(self.DOSCAR_spin_pol.is_spin_polarized)

        self.assertFalse(self.DOSCAR_nonspin_pol.is_spin_polarized)


class ChargeTest(PymatgenTest):
    def setUp(self):
        self.charge2 = Charge(filename=os.path.join(test_dir, "CHARGE.lobster.MnO"))

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
        structure_dict2 = {'lattice': {'c': 3.198244, 'volume': 23.132361565928807, 'b': 3.1982447183003364,
                                       'gamma': 60.00000011873414, 'beta': 60.00000401737447,
                                       'alpha': 60.00000742944491,
                                       'matrix': [[2.769761, 0.0, 1.599122], [0.923254, 2.611356, 1.599122],
                                                  [0.0, 0.0, 3.198244]], 'a': 3.1982443884113985},
                           '@class': 'Structure', 'sites': [{'xyz': [1.846502883732, 1.305680611356, 3.198248797366],
                                                             'properties': {'Loewdin Charges': -1.25,
                                                                            'Mulliken Charges': -1.3},
                                                             'abc': [0.499998, 0.500001, 0.500002],
                                                             'species': [{'occu': 1, 'element': 'O'}], 'label': 'O'},
                                                            {'xyz': [0.0, 0.0, 0.0],
                                                             'properties': {'Loewdin Charges': 1.25,
                                                                            'Mulliken Charges': 1.3},
                                                             'abc': [0.0, 0.0, 0.0],
                                                             'species': [{'occu': 1, 'element': 'Mn'}], 'label': 'Mn'}],
                           'charge': None, '@module': 'pymatgen.core.structure'}
        s2 = Structure.from_dict(structure_dict2)
        self.assertEqual(s2, self.charge2.get_structure_with_charges(os.path.join(this_dir, "POSCAR.MnO")))


if __name__ == "__main__":
    unittest.main()
