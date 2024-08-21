from __future__ import annotations

import gzip
import json

from monty.json import MontyDecoder

from pymatgen.symmetry import site_symmetries as ss
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

__author__ = "Handong Ling"
__version__ = "0.1"
__maintainer__ = "Handong Ling"
__email__ = "handongling@berkeley.edu"
__status__ = "Development"
__date__ = "4/23/19"

TEST_DIR = f"{TEST_FILES_DIR}/symmetry/site_symmetries"


class TestSiteSymmetries(PymatgenTest):
    def setUp(self):
        with gzip.open(f"{TEST_DIR}/point_ops.json.gz", mode="rt") as file:
            self.point_ops = MontyDecoder().process_decoded(json.load(file))

        with gzip.open(f"{TEST_DIR}/shared_ops.json.gz", mode="rt") as file:
            self.shared_ops = MontyDecoder().process_decoded(json.load(file))

        self.piezo_struct = self.get_structure("Pb2TiZrO6")

        # following code can be used to update point_ops/shared_ops reference file
        # def handler(obj):
        #     if hasattr(obj, "as_dict"):
        #         return obj.as_dict()
        #     if isinstance(obj, np.ndarray):
        #         return obj.tolist()
        #     raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")

        # with gzip.open(f"{TEST_FILES}/point_ops.json.gz", mode="wt") as file:
        #     json.dump(self.point_ops, file, default=handler)

    def test_get_site_symmetries(self):
        point_ops = ss.get_site_symmetries(self.piezo_struct)

        assert point_ops == self.point_ops

    def test_get_shared_symmetries_operations(self):
        shared_ops = list(
            map(list, ss.get_shared_symmetry_operations(self.piezo_struct, ss.get_site_symmetries(self.piezo_struct)))
        )
        assert shared_ops == self.shared_ops
