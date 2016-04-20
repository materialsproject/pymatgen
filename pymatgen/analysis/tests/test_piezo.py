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

    def runTest(self):

        piezo_struc = self.get_structure('BaNiO3')
        tensor = np.asarray([[0., 0., 0., 0., 0.03839, 0.],
                             [0., 0., 0., 0.03839, 0., 0.],
                             [6.89822, 6.89822, 27.46280, 0., 0., 0.]]).reshape(3, 6)
        pt = PiezoTensor(tensor)
        full_tensor = [[[0., 0., 0.03839],
                        [0., 0., 0.],
                        [0.03839, 0., 0.]],
                       [[0., 0., 0.],
                        [0., 0., 0.03839],
                        [0.,  0.03839,  0.]],
                       [[6.89822, 0., 0.],
                        [0.,  6.89822, 0.],
                        [0.,  0.,  27.4628]]]
        bad_pt = PiezoTensor([[0., 0., 1., 0., 0.03839, 2.],
                              [0., 0., 0., 0.03839, 0., 0.],
                              [6.89822, 6.89822, 27.4628, 0., 0., 0.]])

        sym_pt = bad_pt.symmeterize(piezo_struc)
        alt_tensor = pt.from_full_tensor(full_tensor)

        self.assertArrayAlmostEqual(pt.full_tensor, full_tensor)
        self.assertArrayEqual(pt, alt_tensor)
        self.assertTrue(pt.is_valid(piezo_struc))

        self.assertTrue(sym_pt.is_valid(piezo_struc))
        self.assertArrayAlmostEqual(sym_pt, pt)

if __name__ == '__main__':
    unittest.main()
