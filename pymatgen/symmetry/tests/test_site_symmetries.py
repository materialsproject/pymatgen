# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import os
import pickle

import numpy as np

from pymatgen.symmetry import site_symmetries as ss
from pymatgen.util.testing import PymatgenTest

__author__ = "Handong Ling"
__version__ = "0.1"
__maintainer__ = "Handong Ling"
__email__ = "handongling@berkeley.edu"
__status__ = "Development"
__date__ = "4/23/19"

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "site_symmetries")


class SiteSymmetriesTest(PymatgenTest):
    def setUp(self):
        with open(os.path.join(test_dir, "point_ops.pkl"), "rb") as f:
            self.point_ops = pickle.load(f)
        with open(os.path.join(test_dir, "shared_ops.pkl"), "rb") as f:
            self.shared_ops = pickle.load(f)
        self.piezo_struc = self.get_structure("Pb2TiZrO6")

    def test_get_site_symmetries(self):
        point_ops = ss.get_site_symmetries(self.piezo_struc)
        # with open(os.path.join(test_dir, "point_ops.pkl"), "wb") as f:
        #     pickle.dump(point_ops, f)  # update test file
        assert np.all(point_ops == self.point_ops)

    def test_get_shared_symmetries_operations(self):
        shared_ops = ss.get_shared_symmetry_operations(self.piezo_struc, ss.get_site_symmetries(self.piezo_struc))
        # with open(os.path.join(test_dir, "shared_ops.pkl"), "wb") as f:
        #     pickle.dump(shared_ops, f)  # update test file
        assert np.all(shared_ops == self.shared_ops)
