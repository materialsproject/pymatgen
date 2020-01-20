#!/usr/bin/env python


__author__ = 'waroquiers'

import networkx as nx
import os
import shutil
from pymatgen.analysis.chemenv.connectivity.environment_nodes import get_environment_node, EnvironmentNode
from pymatgen.util.testing import PymatgenTest
try:
    import bson  # type: ignore  # Ignore bson import for mypy
except ModuleNotFoundError:
    bson = None

json_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "json_test_files")


class EnvironmentNodesTest(PymatgenTest):

    def test_equal(self):
        s = self.get_structure('SiO2')
        en = EnvironmentNode(central_site=s[0], i_central_site=0, ce_symbol='T:4')

        en1 = EnvironmentNode(central_site=s[2], i_central_site=0, ce_symbol='T:4')
        assert en == en1
        assert not en.everything_equal(en1)

        en2 = EnvironmentNode(central_site=s[0], i_central_site=3, ce_symbol='T:4')
        assert en != en2
        assert not en.everything_equal(en2)

        en3 = EnvironmentNode(central_site=s[0], i_central_site=0, ce_symbol='O:6')
        assert en == en3
        assert not en.everything_equal(en3)

        en4 = EnvironmentNode(central_site=s[0], i_central_site=0, ce_symbol='T:4')
        assert en == en4
        assert en.everything_equal(en4)

    def test_as_dict(self):
        s = self.get_structure('SiO2')
        en = EnvironmentNode(central_site=s[2], i_central_site=2, ce_symbol='T:4')

        en_from_dict = EnvironmentNode.from_dict(en.as_dict())
        assert en.everything_equal(en_from_dict)

        if bson is not None:
            bson_data = bson.BSON.encode(en.as_dict())
            en_from_bson = EnvironmentNode.from_dict(bson_data.decode())
            assert en.everything_equal(en_from_bson)

    def test_str(self):
        s = self.get_structure('SiO2')
        en = EnvironmentNode(central_site=s[2], i_central_site=2, ce_symbol='T:4')
        assert str(en) == 'Node #2 Si (T:4)'


if __name__ == "__main__":
    import unittest
    unittest.main()
