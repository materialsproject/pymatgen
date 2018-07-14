# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os
import json

from monty.serialization import loadfn

from pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType
from pymatgen.electronic_structure.dos import CompleteDos, DOS, FermiDos
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


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

class FermiDosTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(test_dir, "complete_dos.json"), "r") as f:
            self.dos = CompleteDos.from_dict(json.load(f))
        self.dos = FermiDos(self.dos)

    def test_doping_fermi(self):
        T = 300
        fermi0 = self.dos.efermi
        frange = [fermi0-0.5, fermi0, fermi0+2.0, fermi0+2.2]
        dopings = [self.dos.get_doping(fermi=f, T=T) for f in frange]
        ref_dopings = [3.48077e+21, 1.9235e+18, -2.6909e+16, -4.8723e+19]
        for i, c_ref in enumerate(ref_dopings):
            self.assertLessEqual(abs(dopings[i]/c_ref-1.0), 0.01)

        calc_fermis = [self.dos.get_fermi(c=c, T=T) for c in ref_dopings]
        for j, f_ref in enumerate(frange):
            self.assertAlmostEqual(calc_fermis[j], f_ref, 4)

        sci_dos = FermiDos(self.dos, bandgap=3.0)
        for i, c_ref in enumerate(ref_dopings):
            if c_ref < 0:
                self.assertAlmostEqual(
                    sci_dos.get_fermi(c_ref, T=T) - frange[i], 0.47, places=2)
            else:
                self.assertAlmostEqual(
                    sci_dos.get_fermi(c_ref, T=T) - frange[i], -0.47, places=2)

        self.assertAlmostEqual(sci_dos.get_fermi_interextrapolated(-1e26, 300),
                               7.5108, 4)
        self.assertAlmostEqual(sci_dos.get_fermi_interextrapolated(1e26, 300),
                               -1.6884, 4)
        self.assertAlmostEqual(sci_dos.get_fermi_interextrapolated(0.0, 300),
                               3.2382, 4)

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


class DOSTest(PymatgenTest):

    def setUp(self):
        with open(os.path.join(test_dir, "complete_dos.json"), "r") as f:
            d = json.load(f)
            y = list(zip(d["densities"]["1"], d["densities"]["-1"]))
            self.dos = DOS(d["energies"], y, d["efermi"])

    def test_get_gap(self):
        dos = self.dos
        self.assertAlmostEqual(dos.get_gap(), 2.0589, 4)
        self.assertEqual(len(dos.x), 301)
        self.assertAlmostEqual(dos.get_interpolated_gap(tol=0.001,
                                                        abs_tol=False,
                                                        spin=None)[0],
                               2.16815942458015, 7)
        self.assertArrayAlmostEqual(dos.get_cbm_vbm(),
                                    (3.8729, 1.8140000000000001))

        self.assertAlmostEqual(dos.get_interpolated_value(9.9)[0],
                               1.744588888888891, 7)
        self.assertAlmostEqual(dos.get_interpolated_value(9.9)[1],
                               1.756888888888886, 7)
        self.assertRaises(ValueError, dos.get_interpolated_value, 1000)

        self.assertArrayAlmostEqual(dos.get_cbm_vbm(spin=Spin.up),
                                    (3.8729, 1.2992999999999999))

        self.assertArrayAlmostEqual(dos.get_cbm_vbm(spin=Spin.down),
                                    (4.645, 1.8140000000000001))


class SpinPolarizationTest(unittest.TestCase):

    def test_spin_polarization(self):

        dos_path = os.path.join(test_dir, "dos_spin_polarization_mp-865805.json")
        dos = loadfn(dos_path)
        self.assertAlmostEqual(dos.spin_polarization, 0.6460514663341762)

if __name__ == '__main__':
    unittest.main()
