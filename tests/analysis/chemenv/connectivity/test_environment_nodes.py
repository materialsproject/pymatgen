from __future__ import annotations

import os

from pymatgen.analysis.chemenv.connectivity.environment_nodes import EnvironmentNode
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

try:
    import bson
except ModuleNotFoundError:
    bson = None  # type: ignore

__author__ = "waroquiers"

json_files_dir = os.path.join(
    TEST_FILES_DIR,
    "chemenv",
    "json_test_files",
)


class TestEnvironmentNodes(PymatgenTest):
    def test_equal(self):
        struct = self.get_structure("SiO2")
        en = EnvironmentNode(central_site=struct[0], i_central_site=0, ce_symbol="T:4")

        en1 = EnvironmentNode(central_site=struct[2], i_central_site=0, ce_symbol="T:4")
        assert en == en1
        assert not en.everything_equal(en1)

        en2 = EnvironmentNode(central_site=struct[0], i_central_site=3, ce_symbol="T:4")
        assert en != en2
        assert not en.everything_equal(en2)

        en3 = EnvironmentNode(central_site=struct[0], i_central_site=0, ce_symbol="O:6")
        assert en == en3
        assert not en.everything_equal(en3)

        en4 = EnvironmentNode(central_site=struct[0], i_central_site=0, ce_symbol="T:4")
        assert en == en4
        assert en.everything_equal(en4)

    def test_as_dict(self):
        struct = self.get_structure("SiO2")
        en = EnvironmentNode(central_site=struct[2], i_central_site=2, ce_symbol="T:4")

        en_from_dict = EnvironmentNode.from_dict(en.as_dict())
        assert en.everything_equal(en_from_dict)

        if bson is not None:
            bson_data = bson.BSON.encode(en.as_dict())
            en_from_bson = EnvironmentNode.from_dict(bson_data.decode())
            assert en.everything_equal(en_from_bson)

    def test_str(self):
        struct = self.get_structure("SiO2")
        en = EnvironmentNode(central_site=struct[2], i_central_site=2, ce_symbol="T:4")
        assert str(en) == "Node #2 Si (T:4)"


if __name__ == "__main__":
    import unittest

    unittest.main()
