# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


__author__ = "Handong Ling"
__version__ = "0.1"
__maintainer__ = "Handong Ling"
__email__ = "handongling@berkeley.edu"
__status__ = "Development"
__date__ = "4/23/19"

import os
import unittest
import numpy as np
import pymatgen
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry import site_symmetries as ss

test_dir = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "test_files", "site_symmetries"
)


class SiteSymmetriesTest(PymatgenTest):
    def setUp(self):
        self.pointops = np.load(
            os.path.join(test_dir, "pointops.npy"), allow_pickle=True
        )
        self.sharedops = np.load(
            os.path.join(test_dir, "sharedops.npy"), allow_pickle=True
        )
        self.piezo_struc = self.get_structure("Pb2TiZrO6")

    def test_get_site_symmetries(self):
        pointops = ss.get_site_symmetries(self.piezo_struc)
        self.assertTrue(np.all(pointops == self.pointops))

    def test_get_shared_symmetries_operations(self):
        sharedops = ss.get_shared_symmetry_operations(self.piezo_struc, self.pointops)
        self.assertTrue(np.all(sharedops == self.sharedops))
