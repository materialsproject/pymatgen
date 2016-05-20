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

    def test_from_voigt(self):
        bad_voigt = np.zeros((3, 7))
        pt = PiezoTensor.from_voigt(self.voigt_matrix)
        self.assertArrayEqual(pt, self.full_tensor_array)
        self.assertRaises(ValueError, PiezoTensor.from_voigt, bad_voigt)
        self.assertArrayEqual(self.voigt_matrix, pt.voigt)

if __name__ == '__main__':
    unittest.main()
