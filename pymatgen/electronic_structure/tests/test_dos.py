#!/usr/bin/python

import unittest
import os
import json

import numpy as np

from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.electronic_structure.dos import Dos, PDos, CompleteDos

import pymatgen

test_dir = os.path.join(os.path.dirname(os.path.abspath(pymatgen.__file__)), '..', 'test_files')

class DosTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "complete_dos.json"), "r") as f:
            self.dos = Dos.from_dict(json.load(f))

    def test_get_gap(self):
        self.assertAlmostEqual(self.dos.get_gap(), 2.0589, 4)
        self.assertEqual(len(self.dos.energies), 301)
        self.assertAlmostEqual(self.dos.get_interpolated_gap(tol=0.001, abs_tol=False, spin=None)[0], 2.16815942458015, 7)
        self.assertAlmostEqual(self.dos.get_cbm_vbm(), (3.8729, 1.8140000000000001))

        self.assertAlmostEqual(self.dos.get_interpolated_value(9.9)[Spin.up], 1.744588888888891, 7)
        self.assertAlmostEqual(self.dos.get_interpolated_value(9.9)[Spin.down], 1.756888888888886, 7)
        self.assertRaises(ValueError, self.dos.get_interpolated_value, 1000)

    def test_get_smeared_densities(self):
        smeared = self.dos.get_smeared_densities(0.2)
        dens = self.dos.densities
        for spin in Spin.all_spins:
            self.assertAlmostEqual(sum(dens[spin]), sum(smeared[spin]))

class PDosTest(unittest.TestCase):

    def setUp(self):
        self.pdos = PDos(5, [0, 5, 10], {Spin.up: [10, 20, 30]}, Orbital.s)

    def test_to_from_dict(self):
        d = self.pdos.to_dict
        pdos = PDos.from_dict(d)
        self.assertTrue(np.allclose(pdos.energies, [0, 5, 10]))

class CompleteDosTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, "complete_dos.json"), "r") as f:
            self.dos = CompleteDos.from_dict(json.load(f))

    def test_get_gap(self):
        self.assertAlmostEqual(self.dos.get_gap(), 2.0589, 4, "Wrong gap from dos!")
        self.assertEqual(len(self.dos.energies), 301)
        self.assertAlmostEqual(self.dos.get_interpolated_gap(tol=0.001, abs_tol=False, spin=None)[0], 2.16815942458015, 7)
        spd_dos = self.dos.get_spd_dos()
        self.assertEqual(len(spd_dos), 3)
        el_dos = self.dos.get_element_dos()
        self.assertEqual(len(el_dos), 4)
        sum_spd = spd_dos['S'] + spd_dos['P'] + spd_dos['D']
        sum_element = None
        for dos in el_dos.values():
            if sum_element == None:
                sum_element = dos
            else:
                sum_element += dos

        #The sums of the SPD or the element doses should be the same.
        self.assertTrue((abs(sum_spd.energies - sum_element.energies) < 0.0001).all())
        self.assertTrue((abs(sum_spd.densities[Spin.up] - sum_element.densities[Spin.up]) < 0.0001).all())
        self.assertTrue((abs(sum_spd.densities[Spin.down] - sum_element.densities[Spin.down]) < 0.0001).all())

        self.assertIsNotNone(self.dos.get_site_dos(self.dos.structure[0]))
        self.assertIsNotNone(self.dos.get_site_orbital_dos(self.dos.structure[0], Orbital.s))
        self.assertAlmostEqual(self.dos.get_cbm_vbm(), (3.8729, 1.8140000000000001))

        self.assertAlmostEqual(self.dos.get_interpolated_value(9.9)[Spin.up], 1.744588888888891, 7)
        self.assertAlmostEqual(self.dos.get_interpolated_value(9.9)[Spin.down], 1.756888888888886, 7)
        self.assertRaises(ValueError, self.dos.get_interpolated_value, 1000)


if __name__ == '__main__':
    unittest.main()

