# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Test for the piezo tensor class
"""


__author__ = "Shyam Dwaraknath"
__version__ = "0.1"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "4/1/16"

import os
import unittest2 as unittest
import numpy as np
from pymatgen.analysis.piezo import PiezoTensor
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class PiezoTest(PymatgenTest):
    def setUp(self):
        self.piezo_struc = self.get_structure('BaNiO3')
        self.voigt_matrix = np.array([[0., 0., 0., 0., 0.03839, 0.],
                                      [0., 0., 0., 0.03839, 0., 0.],
                                      [6.89822, 6.89822, 27.46280, 0., 0., 0.]])
        self.full_tensor_array = [[[0., 0., 0.03839],
                                   [0., 0., 0.],
                                   [0.03839, 0., 0.]],
                                  [[0., 0., 0.],
                                   [0., 0., 0.03839],
                                   [0.,  0.03839,  0.]],
                                  [[6.89822, 0., 0.],
                                   [0.,  6.89822, 0.],
                                   [0.,  0.,  27.4628]]]
    def test_new(self):
        pt = PiezoTensor(self.full_tensor_array)
        self.assertArrayAlmostEqual(pt, self.full_tensor_array)
        bad_dim_array = np.zeros((3, 3))
        self.assertRaises(ValueError, PiezoTensor, bad_dim_array)

        bad_pt = PiezoTensor.from_voigt([[0., 0., 1., 0., 0.03839, 2.],
                                         [0., 0., 0., 0.03839, 0., 0.],
                                         [6.89822, 6.89822, 27.4628, 0., 0., 0.]])

        sym_pt = bad_pt.symmeterized(piezo_struc)
        alt_tensor = pt(full_tensor)

        self.assertArrayAlmostEqual(pt.voigt, full_tensor)
        self.assertArrayEqual(pt, alt_tensor)
        #TODO Recheck this test. Commented out for now to enable Py3k testing.
        self.assertTrue(pt.is_fit_to_structure(piezo_struc))

        self.assertTrue(sym_pt.is_fit_to_structure(piezo_struc))
        self.assertArrayAlmostEqual(sym_pt, pt)

    def test_from_voigt(self):

if __name__ == '__main__':
    unittest.main()
