from __future__ import annotations

import os
import unittest
import warnings

from pytest import approx

from pymatgen.apps.borg.hive import (
    GaussianToComputedEntryDrone,
    SimpleVaspToComputedEntryDrone,
    VaspToComputedEntryDrone,
)
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.util.testing import TEST_FILES_DIR


class TestVaspToComputedEntryDrone(unittest.TestCase):
    def setUp(self):
        self.drone = VaspToComputedEntryDrone(data=["efermi"])
        self.structure_drone = VaspToComputedEntryDrone(True)

    def test_get_valid_paths(self):
        for path in os.walk(TEST_FILES_DIR):
            if path[0] == TEST_FILES_DIR:
                assert len(self.drone.get_valid_paths(path)) > 0

    def test_assimilate(self):
        entry = self.drone.assimilate(TEST_FILES_DIR)
        for p in ["hubbards", "is_hubbard", "potcar_spec", "run_type"]:
            assert p in entry.parameters
        assert entry.data["efermi"] == approx(-6.62148548)
        assert entry.composition.reduced_formula == "Xe"
        assert entry.energy == approx(0.5559329)
        entry = self.structure_drone.assimilate(TEST_FILES_DIR)
        assert entry.composition.reduced_formula == "Xe"
        assert entry.energy == approx(0.5559329)
        assert isinstance(entry, ComputedStructureEntry)
        assert entry.structure is not None
        # assert len(entry.parameters["history"]) == 2

    def tearDown(self):
        warnings.simplefilter("default")

    def test_to_from_dict(self):
        d = self.structure_drone.as_dict()
        drone = VaspToComputedEntryDrone.from_dict(d)
        assert isinstance(drone, VaspToComputedEntryDrone)


class TestSimpleVaspToComputedEntryDrone(unittest.TestCase):
    def setUp(self):
        self.drone = SimpleVaspToComputedEntryDrone()
        self.structure_drone = SimpleVaspToComputedEntryDrone(True)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_get_valid_paths(self):
        for path in os.walk(TEST_FILES_DIR):
            if path[0] == TEST_FILES_DIR:
                assert len(self.drone.get_valid_paths(path)) > 0

    def test_to_from_dict(self):
        d = self.structure_drone.as_dict()
        drone = SimpleVaspToComputedEntryDrone.from_dict(d)
        assert isinstance(drone, SimpleVaspToComputedEntryDrone)


class TestGaussianToComputedEntryDrone(unittest.TestCase):
    def setUp(self):
        self.drone = GaussianToComputedEntryDrone(data=["corrections"])
        self.structure_drone = GaussianToComputedEntryDrone(True)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_get_valid_paths(self):
        for path in os.walk(f"{TEST_FILES_DIR}/molecules"):
            if path[0] == f"{TEST_FILES_DIR}/molecules":
                assert len(self.drone.get_valid_paths(path)) > 0

    def test_assimilate(self):
        test_file = f"{TEST_FILES_DIR}/molecules/methane.log"
        entry = self.drone.assimilate(test_file)
        for p in [
            "functional",
            "basis_set",
            "charge",
            "spin_multiplicity",
            "route_parameters",
        ]:
            assert p in entry.parameters
        for p in ["corrections"]:
            assert p in entry.data

        assert entry.composition.reduced_formula == "H4C"
        assert entry.energy == approx(-39.9768775602)
        entry = self.structure_drone.assimilate(test_file)
        assert entry.composition.reduced_formula == "H4C"
        assert entry.energy == approx(-39.9768775602)
        assert isinstance(entry, ComputedStructureEntry)
        assert entry.structure is not None
        for p in ["properly_terminated", "stationary_type"]:
            assert p in entry.data

    def test_to_from_dict(self):
        d = self.structure_drone.as_dict()
        drone = GaussianToComputedEntryDrone.from_dict(d)
        assert isinstance(drone, GaussianToComputedEntryDrone)
