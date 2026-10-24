from __future__ import annotations

import gzip
import json
import shutil

import pytest
from monty.json import MontyDecoder
from pytest import approx

from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen, order_assimilation
from pymatgen.entries.computed_entries import ComputedEntry
from tests.testing import TEST_FILES_DIR

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__date__ = "Mar 18, 2012"

TEST_DIR = f"{TEST_FILES_DIR}/apps/borg"
VASP_FIXTURES = f"{TEST_FILES_DIR}/io/vasp/fixtures"


@pytest.fixture
def run_tree(tmp_path):
    """Directory tree with two valid runs (Xe, LiH) and one unparsable vasprun."""
    for sub, src in (
        ("xe", f"{TEST_DIR}/vasprun.xml.xe.gz"),
        ("nested/lih", f"{VASP_FIXTURES}/scan_relaxation/vasprun.xml.gz"),
    ):
        (tmp_path / sub).mkdir(parents=True)
        shutil.copy(src, tmp_path / sub / "vasprun.xml.gz")
    (tmp_path / "bad").mkdir()
    with gzip.open(tmp_path / "bad" / "vasprun.xml.gz", "wt", encoding="utf-8") as file:
        file.write("<modeling><not valid xml")
    return tmp_path


class TestBorgQueen:
    def test_get_data(self):
        """Test get data from vasprun.xml.xe.gz file."""
        drone = VaspToComputedEntryDrone()
        queen = BorgQueen(drone, TEST_DIR, 1)
        data = queen.get_data()
        assert len(data) == 1
        assert data[0].energy == approx(0.5559329, 1e-6)

    def test_no_rootpath(self):
        queen = BorgQueen(VaspToComputedEntryDrone())
        assert queen.get_data() == []

    @pytest.mark.parametrize("number_of_drones", [1, 2])
    def test_assimilate_tree(self, run_tree, number_of_drones):
        """Serial and parallel assimilation give the same entries and skip failed runs."""
        queen = BorgQueen(VaspToComputedEntryDrone(), run_tree, number_of_drones)
        data = queen.get_data()
        assert all(isinstance(entry, ComputedEntry) for entry in data)
        energies = {entry.reduced_formula: entry.energy for entry in data}
        assert energies == approx({"Xe": 0.5559329, "LiH": -6.82078391})

    def test_load_data(self):
        drone = VaspToComputedEntryDrone()
        queen = BorgQueen(drone)
        queen.load_data(f"{TEST_DIR}/assimilated.json")
        assert len(queen.get_data()) == 1

    @pytest.mark.parametrize("filename", ["assimilated.json", "assimilated.json.gz"])
    def test_save_load_data(self, tmp_path, filename):
        queen = BorgQueen(VaspToComputedEntryDrone(), TEST_DIR)
        path = tmp_path / filename
        queen.save_data(path)
        if filename.endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8") as file:
                assert "ComputedEntry" in file.read()

        new_queen = BorgQueen(VaspToComputedEntryDrone())
        new_queen.load_data(path)
        (entry,) = new_queen.get_data()
        assert isinstance(entry, ComputedEntry)
        assert entry.reduced_formula == "Xe"
        assert entry.energy == approx(0.5559329)


def test_order_assimilation(tmp_path):
    """order_assimilation (run in worker processes) serializes entries and skips failures."""
    drone = VaspToComputedEntryDrone()
    data: list[str] = []
    status = {"count": 0, "total": 2}
    order_assimilation((TEST_DIR, drone, data, status))
    order_assimilation((tmp_path, drone, data, status))  # no vasprun.xml
    assert status["count"] == 2
    (json_str,) = data
    entry = json.loads(json_str, cls=MontyDecoder)
    assert isinstance(entry, ComputedEntry)
    assert entry.energy == approx(0.5559329)
