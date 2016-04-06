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
import unittest
import numpy as np
from pymatgen.analysis.piezo import PiezoTensor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import PymatgenTest
from pymatgen.matproj.rest import MPRester

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class PiezoTest(PymatgenTest):

    def runTest(self):

        mat_id = 'mp-19241'
        mpr = MPRester()
        piezo_data =mpr.get_entry_by_material_id(mat_id,property_data=['piezo'])
        piezo_struc = mpr.get_structure_by_material_id(mat_id)
        piezo_tensor = PiezoTensor(piezo_data.data['piezo']['piezoelectric_tensor'])
        sg = SpacegroupAnalyzer(piezo_struc)
        self.assertTrue(piezo_tensor.is_valid(sg.get_conventional_standard_structure()))



if __name__ == '__main__':
    unittest.main()
