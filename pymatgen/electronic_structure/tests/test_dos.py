# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest
import os
import json

from pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType
from pymatgen.electronic_structure.dos import CompleteDos

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

import scipy


class DosTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "complete_dos.json"), "r") as f:
            self.dos = CompleteDos.from_dict(json.load(f))

    def test_get_gap(self):
        dos = self.dos
        self.assertAlmostEqual(dos.get_gap(), 2.0589, 4)
        self.assertEqual(len(dos.energies), 301)
        self.assertAlmostEqual(dos.get_interpolated_gap(tol=0.001,
                                                        abs_tol=False,
                                                        spin=None)[0],
                               2.16815942458015, 7)
        self.assertAlmostEqual(dos.get_cbm_vbm(), (3.8729, 1.8140000000000001))

        self.assertAlmostEqual(dos.get_interpolated_value(9.9)[Spin.up],
                               1.744588888888891, 7)
        self.assertAlmostEqual(dos.get_interpolated_value(9.9)[Spin.down],
                               1.756888888888886, 7)
        self.assertRaises(ValueError, dos.get_interpolated_value, 1000)

    def test_get_smeared_densities(self):
        dos = self.dos
        smeared = dos.get_smeared_densities(0.2)
        dens = dos.densities
        for spin in Spin:
            self.assertAlmostEqual(sum(dens[spin]), sum(smeared[spin]))


class CompleteDosTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "complete_dos.json"), "r") as f:
            self.dos = CompleteDos.from_dict(json.load(f))

    def test_get_gap(self):
        dos = self.dos
        self.assertAlmostEqual(dos.get_gap(), 2.0589, 4, "Wrong gap from dos!")
        self.assertEqual(len(dos.energies), 301)
        self.assertAlmostEqual(dos.get_interpolated_gap(tol=0.001,
                                                        abs_tol=False,
                                                        spin=None)[0],
                               2.16815942458015, 7)
        spd_dos = dos.get_spd_dos()
        self.assertEqual(len(spd_dos), 3)
        el_dos = dos.get_element_dos()
        self.assertEqual(len(el_dos), 4)
        sum_spd = spd_dos[OrbitalType.s] + spd_dos[OrbitalType.p] + spd_dos[OrbitalType.d]
        sum_element = None
        for pdos in el_dos.values():
            if sum_element is None:
                sum_element = pdos
            else:
                sum_element += pdos

        #The sums of the SPD or the element doses should be the same.
        self.assertTrue((abs(sum_spd.energies
                             - sum_element.energies) < 0.0001).all())
        self.assertTrue((abs(sum_spd.densities[Spin.up]
                             - sum_element.densities[Spin.up])
                         < 0.0001).all())
        self.assertTrue((abs(sum_spd.densities[Spin.down]
                             - sum_element.densities[Spin.down])
                         < 0.0001).all())

        site = dos.structure[0]
        self.assertIsNotNone(dos.get_site_dos(site))
        self.assertAlmostEqual(sum(dos.get_site_dos(site).get_densities(
            Spin.up)), 2.0391)
        self.assertAlmostEqual(sum(dos.get_site_dos(site).get_densities(
            Spin.down)), 2.0331999999999995)
        self.assertIsNotNone(dos.get_site_orbital_dos(site, Orbital.s))
        egt2g = dos.get_site_t2g_eg_resolved_dos(site)
        self.assertAlmostEqual(sum(egt2g["e_g"].get_densities(Spin.up)),
                               0.0)
        self.assertAlmostEqual(sum(egt2g["t2g"].get_densities(Spin.up)),
                               0.0)
        egt2g = dos.get_site_t2g_eg_resolved_dos(dos.structure[4])
        self.assertAlmostEqual(sum(egt2g["e_g"].get_densities(Spin.up)),
                               15.004399999999997)
        self.assertAlmostEqual(sum(egt2g["t2g"].get_densities(Spin.up)),
                               22.910399999999999)
        self.assertAlmostEqual(dos.get_cbm_vbm(), (3.8729, 1.8140000000000001))

        self.assertAlmostEqual(dos.get_interpolated_value(9.9)[Spin.up],
                               1.744588888888891, 7)
        self.assertAlmostEqual(dos.get_interpolated_value(9.9)[Spin.down],
                               1.756888888888886, 7)
        self.assertRaises(ValueError, dos.get_interpolated_value, 1000)

    def test_to_from_dict(self):
        d = self.dos.as_dict()
        dos = CompleteDos.from_dict(d)
        el_dos = dos.get_element_dos()
        self.assertEqual(len(el_dos), 4)
        spd_dos = dos.get_spd_dos()
        sum_spd = spd_dos[OrbitalType.s] + spd_dos[OrbitalType.p] + spd_dos[OrbitalType.d]
        sum_element = None
        for pdos in el_dos.values():
            if sum_element is None:
                sum_element = pdos
            else:
                sum_element += pdos

        #The sums of the SPD or the element doses should be the same.
        self.assertTrue((abs(sum_spd.energies
                             - sum_element.energies) < 0.0001).all())

    def test_str(self):
        self.assertIsNotNone(str(self.dos))

if __name__ == '__main__':
    unittest.main()
