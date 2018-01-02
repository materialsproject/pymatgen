# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest
import os
import numpy as np
from pymatgen.io.lobster import Cohpcar, Icohplist
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.util.testing import PymatgenTest

__author__ = "Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.2"
__email__ = "esters@uoregon.edu"
__date__ = "Dec 10, 2017"

test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "..", "test_files", "cohp")
this_dir = os.path.dirname(os.path.abspath(__file__))


class CohpcarTest(PymatgenTest):
    def setUp(self):
        os.chdir(test_dir)
        self.cohp_bise = Cohpcar(filename="COHPCAR.lobster.BiSe")
        self.coop_bise = Cohpcar(filename="COOPCAR.lobster.BiSe",
                                 are_coops=True)
        self.cohp_fe = Cohpcar()
        self.coop_fe = Cohpcar(are_coops=True)
        self.orb = Cohpcar(filename="COHPCAR.lobster.orbitalwise")
        self.orb_notot = Cohpcar(filename="COHPCAR.lobster.notot.orbitalwise")

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

    def test_energies(self):
        efermi_bise = 5.90043
        elim_bise = (-0.124679, 11.9255)
        efermi_fe = 9.75576
        elim_fe = (-0.277681, 14.7725)
        self.assertEqual(self.cohp_bise.efermi, efermi_bise)
        self.assertEqual(self.coop_bise.efermi, efermi_bise)
        self.assertEqual(self.cohp_fe.efermi, efermi_fe)
        self.assertEqual(self.coop_fe.efermi, efermi_fe)
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

    def test_cohp_data(self):
        lengths_sites_bise = {"Bi1-Se7": (2.882308829886294, (0, 6)),
                              "Bi1-Se10": (3.1014396233274444, (0, 9)),
                              "Bi2-Se8": (2.8823088298862083, (1, 7)),
                              "Bi2-Se9": (3.1014396233275434, (1, 8)),
                              "Bi3-Se10": (3.0500070394403904, (2, 9)),
                              "Bi3-Se11": (2.9167594580335807, (2, 10)),
                              "Bi4-Se9": (3.05000703944039, (3, 8)),
                              "Bi4-Se12": (2.9167594580335803, (3, 11)),
                              "Bi5-Se12": (3.3752173204052101, (4, 11)),
                              "Bi5-Bi6": (3.0729354518345948, (4, 5)),
                              "Bi6-Se11": (3.3752173204052101, (5, 10))}
        lengths_sites_fe = {"Fe8-Fe7": (2.8318907764979082, (7, 6)),
                            "Fe8-Fe9": (2.4524893531900283, (7, 8))}
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

    def test_orbital_resolved_cohp(self):
        orbitals = [tuple((Orbital(i), Orbital(j))) for j in range(4)
                    for i in range(4)]
        self.assertIsNone(self.cohp_bise.orb_res_cohp)
        self.assertIsNone(self.coop_bise.orb_res_cohp)
        self.assertIsNone(self.cohp_fe.orb_res_cohp)
        self.assertIsNone(self.coop_fe.orb_res_cohp)
        self.assertIsNone(self.orb_notot.cohp_data["Ga1-As2"]["COHP"])
        self.assertIsNone(self.orb_notot.cohp_data["Ga1-As2"]["ICOHP"])
        for orbs in self.orb.orb_res_cohp["Ga1-As2"]:
            orb_set = self.orb.orb_res_cohp["Ga1-As2"][orbs]["orbitals"]
            self.assertEqual(orb_set[0][0], 4)
            self.assertEqual(orb_set[1][0], 4)
            self.assertIn(tuple((orb_set[0][1], orb_set[1][1])), orbitals)

        # The sum of the orbital-resolved COHPs should be approximately
        # the total COHP. Due to small deviations in the LOBSTER calculation,
        # the precision is not very high though.
        cohp = self.orb.cohp_data["Ga1-As2"]["COHP"][Spin.up]
        icohp = self.orb.cohp_data["Ga1-As2"]["ICOHP"][Spin.up]
        tot = np.sum([self.orb.orb_res_cohp["Ga1-As2"][orbs]["COHP"][Spin.up]
                      for orbs in self.orb.orb_res_cohp["Ga1-As2"]], axis=0)
        self.assertArrayAlmostEqual(tot, cohp, decimal=3)
        tot = np.sum([self.orb.orb_res_cohp["Ga1-As2"][orbs]["ICOHP"][Spin.up]
                      for orbs in self.orb.orb_res_cohp["Ga1-As2"]], axis=0)
        self.assertArrayAlmostEqual(tot, icohp, decimal=3)

    def tearDown(self):
        os.chdir(this_dir)


class IcohplistTest(unittest.TestCase):
    def setUp(self):
        os.chdir(test_dir)
        self.icohp_bise = Icohplist(filename="ICOHPLIST.lobster.BiSe")
        self.icoop_bise = Icohplist(filename="ICOOPLIST.lobster.BiSe",
                                    are_coops=True)
        self.icohp_fe = Icohplist()
        self.icoop_fe = Icohplist(are_coops=True)

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
        icohplist_bise = {"Bi1-Se7": {"length": 2.88231, "number_of_bonds": 3,
                                      "icohp": {Spin.up: -2.18042}},
                          "Bi1-Se10": {"length": 3.10144, "number_of_bonds": 3,
                                       "icohp": {Spin.up: -1.14347}},
                          "Bi2-Se8": {"length": 2.88231, "number_of_bonds": 3,
                                      "icohp": {Spin.up: -2.18042}},
                          "Bi2-Se9": {"length": 3.10144, "number_of_bonds": 3,
                                      "icohp": {Spin.up: -1.14348}},
                          "Bi3-Se10": {"length": 3.05001, "number_of_bonds": 3,
                                       "icohp": {Spin.up: -1.30006}},
                          "Bi3-Se11": {"length": 2.91676, "number_of_bonds": 3,
                                       "icohp": {Spin.up: -1.96843}},
                          "Bi4-Se9": {"length": 3.05001, "number_of_bonds": 3,
                                      "icohp": {Spin.up: -1.30006}},
                          "Bi4-Se12": {"length": 2.91676, "number_of_bonds": 3,
                                       "icohp": {Spin.up: -1.96843}},
                          "Bi5-Se12": {"length": 3.37522, "number_of_bonds": 3,
                                       "icohp": {Spin.up: -0.47531}},
                          "Bi5-Bi6": {"length": 3.07294, "number_of_bonds": 3,
                                      "icohp": {Spin.up: -2.38796}},
                          "Bi6-Se11": {"length": 3.37522, "number_of_bonds": 3,
                                       "icohp": {Spin.up: -0.47531}}}
        icooplist_bise = {"Bi1-Se7": {"length": 2.88231, "number_of_bonds": 3,
                                      "icohp": {Spin.up: 0.14245}},
                          "Bi1-Se10": {"length": 3.10144, "number_of_bonds": 3,
                                       "icohp": {Spin.up: -0.04118}},
                          "Bi2-Se8": {"length": 2.88231, "number_of_bonds": 3,
                                      "icohp": {Spin.up: 0.14245}},
                          "Bi2-Se9": {"length": 3.10144, "number_of_bonds": 3,
                                      "icohp": {Spin.up: -0.04118}},
                          "Bi3-Se10": {"length": 3.05001, "number_of_bonds": 3,
                                       "icohp": {Spin.up: -0.03516}},
                          "Bi3-Se11": {"length": 2.91676, "number_of_bonds": 3,
                                       "icohp": {Spin.up: 0.10745}},
                          "Bi4-Se9": {"length": 3.05001, "number_of_bonds": 3,
                                      "icohp": {Spin.up: -0.03516}},
                          "Bi4-Se12": {"length": 2.91676, "number_of_bonds": 3,
                                       "icohp": {Spin.up: 0.10745}},
                          "Bi5-Se12": {"length": 3.37522, "number_of_bonds": 3,
                                       "icohp": {Spin.up: -0.12395}},
                          "Bi5-Bi6": {"length": 3.07294, "number_of_bonds": 3,
                                      "icohp": {Spin.up: 0.24714}},
                          "Bi6-Se11": {"length": 3.37522, "number_of_bonds": 3,
                                       "icohp": {Spin.up: -0.12395}}}
        icohplist_fe = {"Fe8-Fe7": {"length": 2.83189, "number_of_bonds": 2,
                                    "icohp": {Spin.up: -0.10218,
                                              Spin.down: -0.19701}},
                        "Fe8-Fe9": {"length": 2.45249, "number_of_bonds": 1,
                                    "icohp": {Spin.up: -0.28485,
                                              Spin.down: -0.58279}}}
        icooplist_fe = {"Fe8-Fe7": {"length": 2.83189, "number_of_bonds": 2,
                                    "icohp": {Spin.up: -0.11389,
                                              Spin.down: -0.20828}},
                        "Fe8-Fe9": {"length": 2.45249, "number_of_bonds": 1,
                                    "icohp": {Spin.up: -0.04087,
                                              Spin.down: -0.05756}}}

        self.assertEqual(icohplist_bise, self.icohp_bise.icohplist)
        self.assertEqual(icooplist_bise, self.icoop_bise.icohplist)
        self.assertEqual(icooplist_fe, self.icoop_fe.icohplist)
        self.assertEqual(icohplist_fe, self.icohp_fe.icohplist)

    def tearDown(self):
        os.chdir(this_dir)

if __name__ == "__main__":
    unittest.main()
