#!/usr/bin/python

import unittest
import pickle
import os

from pymatgen.core.electronic_structure import Orbital, Spin

module_dir = os.path.dirname(os.path.abspath(__file__))

class SpinTest(unittest.TestCase):
    
    def test_init(self):
        self.assertEquals(int(Spin.up), 1)
        self.assertEquals(int(Spin.down), -1)
    
class OrbitalTest(unittest.TestCase):
    
    def test_init(self):
        self.assertEqual(Orbital.from_vasp_index(1), Orbital.py)

class DosTest(unittest.TestCase):
    
    def setUp(self):
        with open(os.path.join(module_dir,"dos_test_file.pkl"), "rb") as f:
            self.dos = pickle.load(f)        
                
    def test_get_gap(self):
        self.assertAlmostEqual(self.dos.get_gap(),2.0589,4,"Wrong gap from dos!")
        self.assertEqual(len(self.dos.energies),301)
        self.assertAlmostEqual(self.dos.get_interpolated_gap(tol=0.001,abs_tol=False,spin=None)[0], 2.16815942458015, 7)
        spd_dos = self.dos.get_spd_dos()
        self.assertEqual(len(spd_dos),3)
        el_dos = self.dos.get_element_dos()
        self.assertEqual(len(el_dos),4)
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
        
        self.assertAlmostEqual(self.dos.get_interpolated_value(9.9)[Spin.up],1.744588888888891, 7)
        self.assertAlmostEqual(self.dos.get_interpolated_value(9.9)[Spin.down],1.756888888888886, 7)
        self.assertRaises(ValueError, self.dos.get_interpolated_value, 1000)

if __name__ == '__main__':
    unittest.main()

