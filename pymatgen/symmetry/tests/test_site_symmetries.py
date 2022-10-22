# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


__author__ = "Handong Ling"
__version__ = "0.1"
__maintainer__ = "Handong Ling"
__email__ = "handongling@berkeley.edu"
__status__ = "Development"
__date__ = "4/23/19"

import os

import numpy as np

from pymatgen.symmetry import site_symmetries as ss
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "site_symmetries")


class SiteSymmetriesTest(PymatgenTest):
    def setUp(self):
        self.point_ops = np.load(os.path.join(test_dir, "point_ops.npy"), allow_pickle=True)
        self.shared_ops = np.load(os.path.join(test_dir, "shared_ops.npy"), allow_pickle=True)
        self.piezo_struc = self.get_structure("Pb2TiZrO6")

    def test_get_site_symmetries(self):
        point_ops = ss.get_site_symmetries(self.piezo_struc)
        assert np.all(point_ops == self.point_ops)

    def test_get_shared_symmetries_operations(self):
        shared_ops = ss.get_shared_symmetry_operations(self.piezo_struc, self.point_ops)
        assert np.all(shared_ops == self.shared_ops)
