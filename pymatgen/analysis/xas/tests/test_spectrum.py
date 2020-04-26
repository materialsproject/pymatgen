# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.analysis.xas.spectrum import XAS, site_weighted_spectrum
from pymatgen.util.testing import PymatgenTest
from pymatgen.core import Element
from monty.json import MontyDecoder
import numpy as np
import unittest
import json
import os

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files/spectrum_test")

with open(os.path.join(test_dir, 'LiCoO2_k_xanes.json')) as fp:
    k_xanes_dict = json.load(fp, cls=MontyDecoder)
with open(os.path.join(test_dir, 'LiCoO2_k_exafs.json')) as fp:
    k_exafs_dict = json.load(fp, cls=MontyDecoder)
with open(os.path.join(test_dir, 'ZnO_l2_xanes.json')) as fp:
    l2_xanes_dict = json.load(fp, cls=MontyDecoder)
with open(os.path.join(test_dir, 'ZnO_l3_xanes.json')) as fp:
    l3_xanes_dict = json.load(fp, cls=MontyDecoder)
with open(os.path.join(test_dir, 'site1_k_xanes.json')) as fp:
    site1_xanes_dict = json.load(fp, cls=MontyDecoder)
with open(os.path.join(test_dir, 'site2_k_xanes.json')) as fp:
    site2_xanes_dict = json.load(fp, cls=MontyDecoder)


class XASTest(PymatgenTest):
    def setUp(self):
        self.k_xanes = XAS.from_dict(k_xanes_dict)
        self.k_exafs = XAS.from_dict(k_exafs_dict)
        self.l2_xanes = XAS.from_dict(l2_xanes_dict)
        self.l3_xanes = XAS.from_dict(l3_xanes_dict)
        self.site1_xanes = XAS.from_dict(site1_xanes_dict)
        self.site2_xanes = XAS.from_dict(site2_xanes_dict)

    def test_e0(self):
        self.assertAlmostEqual(7728.565, self.k_xanes.e0)

    def test_k(self):
        self.assertEqual(len(self.k_xanes.x), len(self.k_xanes.k))
        self.assertAlmostEqual(self.k_xanes.e0,
                               self.k_xanes.x[self.k_xanes.k.index(0)])

    def test_normalization(self):
        self.k_xanes.normalize(mode="sum")
        self.assertAlmostEqual(1.0, np.sum(self.k_xanes.y))

    def test_add_mul(self):
        scaled_spect = self.k_xanes + self.k_xanes
        scaled_spect2 = self.k_xanes * 3
        self.assertTrue(np.allclose(scaled_spect.y, 2 * self.k_xanes.y))
        self.assertTrue(np.allclose(scaled_spect2.y, 3 * self.k_xanes.y))
        self.assertAlmostEqual(0.274302,
                               self.k_xanes.get_interpolated_value(7720.422), 3)

    def test_to_from_dict(self):
        s = XAS.from_dict(self.k_xanes.as_dict())
        self.assertArrayAlmostEqual(s.y, self.k_xanes.y)

    def test_attributes(self):
        self.assertArrayEqual(self.k_xanes.energy, self.k_xanes.x)
        self.assertArrayEqual(self.k_xanes.intensity, self.k_xanes.y)

    def test_str(self):
        self.assertIsNotNone(str(self.k_xanes))

    def test_stitch_xafs(self):
        self.assertRaises(ValueError, XAS.stitch, self.k_xanes, self.k_exafs,
                          mode="invalid")
        xafs = XAS.stitch(self.k_xanes, self.k_exafs, mode="XAFS")
        self.assertIsInstance(xafs, XAS)
        self.assertEqual("XAFS", xafs.spectrum_type)
        self.assertEqual(len(xafs.x), 500)
        self.assertAlmostEqual(min(xafs.x), min(self.k_xanes.x), 2)
        self.assertAlmostEqual(max(xafs.y), max(self.k_xanes.y), 2)
        self.assertAlmostEqual(xafs.x[np.argmax(np.gradient(xafs.y) /
                                                np.gradient(xafs.x))],
                               self.k_xanes.e0, 2)
        self.assertRaises(ValueError, XAS.stitch,
                          self.k_xanes, self.l2_xanes, mode="XAFS")
        self.k_xanes.x = np.zeros(100)
        self.assertRaises(ValueError, XAS.stitch,
                          self.k_xanes, self.k_exafs)
        self.k_xanes.absorbing_element = Element("Pt")
        self.assertRaises(ValueError, XAS.stitch,
                          self.k_xanes, self.k_exafs, mode="XAFS")

    def test_stitch_l23(self):
        l23 = XAS.stitch(self.l2_xanes, self.l3_xanes, 100, mode="L23")
        self.assertIsInstance(l23, XAS)
        self.assertEqual("L23", l23.edge)
        self.assertAlmostEqual(min(l23.x), min(self.l3_xanes.x), 3)
        self.assertAlmostEqual(max(l23.x), max(self.l2_xanes.x), 3)
        self.assertTrue(np.greater_equal(l23.y, self.l2_xanes.y).all())
        self.assertEqual(len(l23.x), 100)
        self.l2_xanes.spectrum_type = "EXAFS"
        self.assertRaises(ValueError,  XAS.stitch,
                          self.l2_xanes, self.l3_xanes, mode="L23")
        self.l2_xanes.absorbing_element = Element("Pt")
        self.assertRaises(ValueError, XAS.stitch,
                          self.l2_xanes, self.l3_xanes, mode="L23")
        self.assertRaises(ValueError, XAS.stitch,
                          self.k_xanes, self.l3_xanes, mode="L23")

    def test_site_weighted_spectrum(self):
        weighted_spectrum = site_weighted_spectrum([self.site1_xanes,
                                                    self.site2_xanes])
        self.assertIsInstance(weighted_spectrum, XAS)
        self.assertTrue(len(weighted_spectrum.x), 500)
        # The site multiplicities for site1 and site2 are 4 and 2, respectively.
        self.assertAlmostEqual(weighted_spectrum.y[0], (4*self.site1_xanes.y[0] +
                               2*self.site2_xanes.y[0])/6, 2)
        self.assertEqual(min(weighted_spectrum.x),
                         max(min(self.site1_xanes.x), min(self.site2_xanes.x)))
        self.site2_xanes.absorbing_index = self.site1_xanes.absorbing_index
        self.assertRaises(ValueError, site_weighted_spectrum,
                          [self.site1_xanes, self.site2_xanes])


if __name__ == '__main__':
    unittest.main()
