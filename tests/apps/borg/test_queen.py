from __future__ import annotations

from pytest import approx

from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.util.testing import TEST_FILES_DIR

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__date__ = "Mar 18, 2012"

TEST_DIR = f"{TEST_FILES_DIR}/apps/borg"


class TestBorgQueen:
    def test_get_data(self):
        """Test get data from vasprun.xml.xe.gz file."""
        drone = VaspToComputedEntryDrone()
        queen = BorgQueen(drone, TEST_DIR, 1)
        data = queen.get_data()
        assert len(data) == 1
        assert data[0].energy == approx(0.5559329, 1e-6)

    def test_load_data(self):
        drone = VaspToComputedEntryDrone()
        queen = BorgQueen(drone)
        queen.load_data(f"{TEST_DIR}/assimilated.json")
        assert len(queen.get_data()) == 1
