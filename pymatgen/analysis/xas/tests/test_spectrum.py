# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.analysis.xas.spectrum import *
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


class XANESTest(PymatgenTest):
    def setUp(self):
        self.xanes = XANES.from_dict(k_xanes_dict)

    def test_e0(self):
        self.assertAlmostEqual(7728.565, self.xanes.e0)

    def test_k(self):
        self.assertEqual(len(self.xanes.x), len(self.xanes.k))
        self.assertAlmostEqual(self.xanes.e0,
                               self.xanes.x[self.xanes.k.index(0)])

    def test_normalization(self):
        self.xanes.normalize(mode="sum")
        self.assertAlmostEqual(1.0, np.sum(self.xanes.y))

    def test_add_mul(self):
        scaled_spect = self.xanes + self.xanes
        scaled_spect2 = self.xanes * 3
        self.assertTrue(np.allclose(scaled_spect.y, 2 * self.xanes.y))
        self.assertTrue(np.allclose(scaled_spect2.y, 3 * self.xanes.y))
        self.assertAlmostEqual(0.274302,
                               self.xanes.get_interpolated_value(7720.422), 4)

    def test_to_from_dict(self):
        s = XANES.from_dict(self.xanes.as_dict())
        self.assertArrayAlmostEqual(s.y, self.xanes.y)

    def test_attributes(self):
        self.assertArrayEqual(self.xanes.energy, self.xanes.x)
        self.assertArrayEqual(self.xanes.intensity, self.xanes.y)

    def test_str(self):
        self.assertIsNotNone(str(self.xanes))


class StitchTest(PymatgenTest):
    def setup(self):
        pass

    def test_stitch_xanes_exafs(self):
        k_xanes = XANES.from_dict(k_xanes_dict)
        k_exafs = EXAFS.from_dict(k_exafs_dict)
        l2_xanes = XANES.from_dict(l2_xanes_dict)
        xas_x, xas_y = stitch_xanes_exafs(k_xanes, k_exafs)
        self.assertEqual(len(xas_x), 500)
        self.assertAlmostEqual(min(xas_x), min(k_xanes.x), 2)
        self.assertAlmostEqual(max(xas_y), max(k_xanes.y), 2)
        self.assertAlmostEqual(xas_x[np.argmax(np.gradient(xas_y) / np.gradient(xas_x))],
                               k_xanes.e0, 2)
        self.assertRaises(ValueError, stitch_xanes_exafs, k_xanes, l2_xanes)
        k_xanes.x = np.zeros(100)
        self.assertRaises(ValueError, stitch_xanes_exafs, k_xanes, k_exafs)

    def test_stitch_l23(self):
        l2_xanes = XANES.from_dict(l2_xanes_dict)
        l3_xanes = XANES.from_dict(l3_xanes_dict)
        l23_x, l23_y = stitch_l23(l2_xanes, l3_xanes, 100)
        self.assertAlmostEqual(min(l23_x), min(l3_xanes.x), 3)
        self.assertAlmostEqual(max(l23_x), max(l2_xanes.x), 3)
        self.assertTrue(np.greater_equal(l23_y, l2_xanes.y).all())
        self.assertEqual(len(l23_x), 100)
        l2_xanes.absorption_element = Element("Pt")
        self.assertRaises(ValueError, stitch_l23, l2_xanes, l3_xanes)





if __name__ == '__main__':
    unittest.main()
