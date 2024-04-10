from __future__ import annotations

import json

from pymatgen.analysis.chemenv.connectivity.environment_nodes import EnvironmentNode
from pymatgen.util.testing import PymatgenTest

try:
    import bson
except ModuleNotFoundError:
    bson = None  # type: ignore[assignment]

__author__ = "waroquiers"


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
        env_node = EnvironmentNode(central_site=struct[2], i_central_site=2, ce_symbol="T:4")

        env_node_from_dict = EnvironmentNode.from_dict(env_node.as_dict())
        assert env_node.everything_equal(env_node_from_dict)

        json_str = self.assert_msonable(env_node)
        env_node_from_json = EnvironmentNode.from_dict(json.loads(json_str))
        assert env_node.everything_equal(env_node_from_json)

    def test_str(self):
        struct = self.get_structure("SiO2")
        env_node = EnvironmentNode(central_site=struct[2], i_central_site=2, ce_symbol="T:4")
        assert str(env_node) == "Node #2 Si (T:4)"
