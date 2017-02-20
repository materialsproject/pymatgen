from __future__ import unicode_literals

import unittest
import os
import json

from pymatgen.core.periodic_table import Element
from pymatgen.phonon.dos import PhononDos, CompletePhononDos
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

import scipy


class DosTest(PymatgenTest):

    def setUp(self):
        with open(os.path.join(test_dir, "NaCl_ph_dos.json"), "r") as f:
            self.dos = PhononDos.from_dict(json.load(f))

    def test_properties(self):
        self.assertAlmostEqual(self.dos.densities[15], 0.0001665998)
        self.assertAlmostEqual(self.dos.frequencies[20], 0.0894965119)
        self.assertAlmostEqual(self.dos.get_interpolated_value(3.), 1.2915532670115628)
        self.assertEqual(len(self.dos.frequencies), 201)
        self.assertEqual(len(self.dos.densities), 201)

    def test_get_smeared_densities(self):
        smeared = self.dos.get_smeared_densities(0.01)
        self.assertAlmostEqual(smeared[20], 0.00084171007635058825)
        dens = self.dos.densities
        self.assertAlmostEqual(sum(dens), sum(smeared))

    def test_dict_methods(self):
        s = json.dumps(self.dos.as_dict())
        self.assertIsNotNone(s)
        self.assertMSONable(self.dos)


class CompleteDosTest(PymatgenTest):

    def setUp(self):
        with open(os.path.join(test_dir, "NaCl_complete_ph_dos.json"), "r") as f:
            self.cdos = CompletePhononDos.from_dict(json.load(f))

    def test_properties(self):
        site_Na = self.cdos.structure[0]
        site_Cl = self.cdos.structure[1]

        self.assertEqual(len(self.cdos.frequencies), 201)
        self.assertAlmostEqual(self.cdos.pdos[site_Na][30],  0.008058208)
        self.assertAlmostEqual(self.cdos.get_site_dos(site_Na).densities[30],  0.008058208)
        self.assertAlmostEqual(self.cdos.pdos[site_Cl][30],  0.0119040783)

        self.assertIn(Element.Na, self.cdos.get_element_dos())
        self.assertIn(Element.Cl, self.cdos.get_element_dos())

        sum_dos = self.cdos.get_element_dos()[Element.Na] + self.cdos.get_element_dos()[Element.Cl]
        self.assertArrayAlmostEqual(sum_dos.frequencies, self.cdos.frequencies)

    def test_dict_methods(self):
        s = json.dumps(self.cdos.as_dict())
        self.assertIsNotNone(s)
        self.assertMSONable(self.cdos)

    def test_str(self):
        self.assertIsNotNone(str(self.cdos))

if __name__ == '__main__':
    unittest.main()
