# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
FIXME: Proper module docstring
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Nov 9, 2012"

import unittest2 as unittest
import os
import json

from monty.json import MontyDecoder
from pymatgen.entries.entry_tools import group_entries_by_structure

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class FuncTest(unittest.TestCase):

    def test_group_entries_by_structure(self):
        with open(os.path.join(test_dir, "TiO2_entries.json"), "r") as f:
            entries = json.load(f, cls=MontyDecoder)
        groups = group_entries_by_structure(entries)
        self.assertEqual(sorted([len(g) for g in groups]),
                         [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4])
        self.assertLess(len(groups), len(entries))
        #Make sure no entries are left behind
        self.assertEqual(sum([len(g) for g in groups]), len(entries))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
