# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on Mar 18, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 18, 2012"

import unittest
import os
import warnings

from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class BorgQueenTest(unittest.TestCase):

    def test_get_data(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            drone = VaspToComputedEntryDrone()
            self.queen = BorgQueen(drone, test_dir, 1)
            data = self.queen.get_data()
            self.assertEqual(len(data), 8)

    def test_load_data(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            drone = VaspToComputedEntryDrone()
            queen = BorgQueen(drone)
            queen.load_data(os.path.join(test_dir, "assimilated.json"))
            self.assertEqual(len(queen.get_data()), 1)


if __name__ == "__main__":
    unittest.main()
