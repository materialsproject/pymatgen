# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import unittest
import json
import os
from pymatgen.electronic_structure.cohp import CompleteCohp, Cohp
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "cohp")


class CohpTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(test_dir, "cohp.json"), "r") as f:
            self.cohp = Cohp.from_dict(json.load(f))
        self.cohp_only = Cohp(self.cohp.efermi,
                              self.cohp.energies,
                              self.cohp.cohp)
        with open(os.path.join(test_dir, "coop.json"), "r") as f:
            self.coop = Cohp.from_dict(json.load(f))

    def test_as_from_dict(self):
        with open(os.path.join(test_dir, "cohp.json"), "r") as f:
            cohp_dict = json.load(f)
        self.assertEqual(self.cohp.as_dict(), cohp_dict)

    def test_attributes(self):
        self.assertEqual(len(self.cohp.energies), 301)
        self.assertEqual(self.cohp.efermi, 9.75576)
        self.assertEqual(self.coop.efermi, 5.90043)
        self.assertFalse(self.cohp.are_coops)
        self.assertTrue(self.coop.are_coops)

    def test_get_icohp(self):
        self.assertEqual(self.cohp.get_icohp(),
                         self.cohp.get_cohp(integrated=True))
        self.assertEqual(None, self.cohp_only.get_icohp())

    def test_get_interpolated_value(self):
        # icohp_ef are the ICHOP(Ef) values taken from
        # the ICOHPLIST.lobster file.
        icohp_ef_dict = {Spin.up: -0.10218, Spin.down: -0.19701}
        icoop_ef_dict = {Spin.up: 0.24714}
        icohp_ef = self.cohp.get_interpolated_value(self.cohp.efermi,
                                                    integrated=True)
        icoop_ef = self.coop.get_interpolated_value(self.coop.efermi,
                                                    integrated=True)
        self.assertAlmostEqual(icohp_ef_dict, icohp_ef)
        self.assertAlmostEqual(icoop_ef_dict, icoop_ef)
        with self.assertRaises(ValueError):
            self.cohp_only.get_interpolated_value(5.0, integrated=True)

    def test_str(self):
        with open(os.path.join(test_dir, "cohp.str"), "rt") as f:
            str_cohp = f.read()
        with open(os.path.join(test_dir, "coop.str"), "rt") as f:
            str_coop = f.read()
        self.assertEqual(self.cohp.__str__(), str_cohp)
        self.assertEqual(self.coop.__str__(), str_coop)


class CompleteCohpTest(PymatgenTest):
    def setUp(self):
        filepath = os.path.join(test_dir, "complete_cohp_lobster.json")
        with open(filepath, "r") as f:
            self.cohp_lobster_dict = CompleteCohp.from_dict(json.load(f))
        filepath = os.path.join(test_dir, "complete_coop_lobster.json")
        with open(filepath, "r") as f:
            self.coop_lobster_dict = CompleteCohp.from_dict(json.load(f))
        filepath = os.path.join(test_dir, "complete_cohp_lmto.json")
        with open(filepath, "r") as f:
            self.cohp_lmto_dict = CompleteCohp.from_dict(json.load(f))
        filepath = os.path.join(test_dir, "complete_cohp_orbitalwise.json")
        with open(filepath, "r") as f:
            self.cohp_orb_dict = CompleteCohp.from_dict(json.load(f))

        filepath = os.path.join(test_dir, "COPL.BiSe")
        structure = os.path.join(test_dir, "CTRL.BiSe")
        self.cohp_lmto = CompleteCohp.from_file("lmto", filename=filepath,
                                                structure_file=structure)
        filepath = os.path.join(test_dir, "COHPCAR.lobster")
        structure = os.path.join(test_dir, "POSCAR")
        self.cohp_lobster = CompleteCohp.from_file("lobster",
                                                   filename=filepath,
                                                   structure_file=structure)
        filepath = os.path.join(test_dir, "COOPCAR.lobster.BiSe")
        structure = os.path.join(test_dir, "POSCAR.BiSe")
        self.coop_lobster = CompleteCohp.from_file("lobster",
                                                   filename=filepath,
                                                   structure_file=structure,
                                                   are_coops=True)
        filepath = os.path.join(test_dir, "COHPCAR.lobster.orbitalwise")
        structure = os.path.join(test_dir, "POSCAR.orbitalwise")
        self.cohp_orb = CompleteCohp.from_file("lobster",
                                               filename=filepath,
                                               structure_file=structure)
        filepath = os.path.join(test_dir, "COHPCAR.lobster.notot.orbitalwise")
        self.cohp_notot = CompleteCohp.from_file("lobster",
                                                 filename=filepath,
                                                 structure_file=structure)

    def test_attiributes(self):
        self.assertFalse(self.cohp_lobster.are_coops)
        self.assertFalse(self.cohp_lobster_dict.are_coops)
        self.assertFalse(self.cohp_lmto.are_coops)
        self.assertFalse(self.cohp_lmto_dict.are_coops)
        self.assertTrue(self.coop_lobster.are_coops)
        self.assertTrue(self.coop_lobster_dict.are_coops)
        self.assertEqual(len(self.cohp_lobster.energies), 301)
        self.assertEqual(len(self.cohp_lmto.energies), 801)
        self.assertEqual(len(self.coop_lobster.energies), 241)
        self.assertEqual(self.cohp_lobster.efermi, 9.75576)
        self.assertEqual(self.cohp_lmto.efermi, -2.3433)
        self.assertEqual(self.coop_lobster.efermi, 5.90043)

    def test_dict(self):
        # The json files are dict representations of the COHPs from the LMTO
        # and LOBSTER calculations and should thus be the same.

        self.assertEqual(self.cohp_lobster.as_dict(),
                         self.cohp_lobster_dict.as_dict())
        self.assertEqual(self.cohp_orb.as_dict(),
                         self.cohp_orb_dict.as_dict())

        # Testing the LMTO dicts will be more involved. Since the average
        # is calculated and not read, there may be differences in rounding
        # with a very small number of matrix elements, which would cause the
        # test to fail
        for key in ["COHP", "ICOHP"]:
            self.assertArrayAlmostEqual(
                self.cohp_lmto.as_dict()[key]["average"]["1"],
                self.cohp_lmto_dict.as_dict()[key]["average"]["1"], 5)
        for key in self.cohp_lmto.as_dict():
            if key not in ["COHP", "ICOHP"]:
                self.assertEqual(self.cohp_lmto.as_dict()[key],
                                 self.cohp_lmto_dict.as_dict()[key])
            else:
                for bond in self.cohp_lmto.as_dict()[key]:
                    if bond != "average":
                        self.assertEqual(self.cohp_lmto.as_dict()[key][bond],
                                         self.cohp_lmto_dict.as_dict()[key][bond])

    def test_icohp_values(self):
        # icohp_ef are the ICHOP(Ef) values taken from
        # the ICOHPLIST.lobster file.
        icohp_ef_dict = {"Fe8-Fe7": {Spin.up: -0.10218, Spin.down: -0.19701},
                         "Fe8-Fe9": {Spin.up: -0.28485, Spin.down: -0.58279}}
        all_cohps_lobster = self.cohp_lobster.all_cohps
        for bond in icohp_ef_dict:
            icohp_ef = all_cohps_lobster[bond].get_interpolated_value(
                           self.cohp_lobster.efermi, integrated=True)
            self.assertEqual(icohp_ef_dict[bond], icohp_ef)

        icoop_ef_dict = {"Bi1-Se7": {Spin.up: 0.14245},
                         "Bi1-Se10": {Spin.up: -0.04118},
                         "Bi2-Se8": {Spin.up: 0.14245},
                         "Bi2-Se9": {Spin.up: -0.04118},
                         "Bi3-Se10": {Spin.up: -0.03516},
                         "Bi3-Se11": {Spin.up: 0.10745},
                         "Bi4-Se9": {Spin.up: -0.03516},
                         "Bi4-Se12": {Spin.up: 0.10745},
                         "Bi5-Se12": {Spin.up: -0.12395},
                         "Bi5-Bi6": {Spin.up: 0.24714},
                         "Bi6-Se11": {Spin.up: -0.12395}}
        all_coops_lobster = self.coop_lobster.all_cohps
        for bond in icoop_ef_dict:
            icoop_ef = all_coops_lobster[bond].get_interpolated_value(
                           self.coop_lobster.efermi, integrated=True)
            self.assertEqual(icoop_ef_dict[bond], icoop_ef)

    def test_orbital_resolved_cohp(self):
        # When read from a COHPCAR file, total COHPs are calculated from
        # the orbital-resolved COHPs if the total is missing. This may be
        # case for LOBSTER version 2.2.0 and earlier due to a bug with the
        # cohpgenerator keyword. The calculated total should be approximately
        # the total COHP calculated by LOBSTER. Due to numerical errors in
        # the LOBSTER calculation, the precision is not very high though.
        self.assertArrayAlmostEqual(
            self.cohp_orb.all_cohps["Ga1-As2"].cohp[Spin.up],
            self.cohp_notot.all_cohps["Ga1-As2"].cohp[Spin.up], decimal=3)
        self.assertArrayAlmostEqual(
            self.cohp_orb.all_cohps["Ga1-As2"].icohp[Spin.up],
            self.cohp_notot.all_cohps["Ga1-As2"].icohp[Spin.up], decimal=3)

        # Tests different methods for getting orbital-resolved COHPs
        ref = self.cohp_orb.orb_res_cohp["Ga1-As2"]["4s-4px"]
        cohp_label = self.cohp_orb.get_orbital_resolved_cohp("Ga1-As2",
                                                             "4s-4px")
        self.assertEqual(cohp_label.cohp, ref["COHP"])
        self.assertEqual(cohp_label.icohp, ref["ICOHP"])
        orbitals = [[Orbital.s, Orbital.px], ["s", "px"], [0, 3]]
        cohps = [self.cohp_orb.get_orbital_resolved_cohp("Ga1-As2",
                 [[4, orb[0]], [4, orb[1]]]) for orb in orbitals]
        for cohp in cohps:
            self.assertEqual(cohp.as_dict(), cohp_label.as_dict())


if __name__ == "__main__":
    unittest.main()
