"""Created on Mar 18, 2012."""


from __future__ import annotations

import unittest
import warnings

from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.util.testing import TEST_FILES_DIR

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 18, 2012"


class TestBorgQueen(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_get_data(self):
        drone = VaspToComputedEntryDrone()
        self.queen = BorgQueen(drone, TEST_FILES_DIR, 1)
        data = self.queen.get_data()
        assert len(data) == 15

    def test_load_data(self):
        drone = VaspToComputedEntryDrone()
        queen = BorgQueen(drone)
        queen.load_data(f"{TEST_FILES_DIR}/assimilated.json")
        assert len(queen.get_data()) == 1
