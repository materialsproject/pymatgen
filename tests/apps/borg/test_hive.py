from __future__ import annotations

import os
from unittest import TestCase

from pytest import approx

from pymatgen.apps.borg.hive import (
    GaussianToComputedEntryDrone,
    SimpleVaspToComputedEntryDrone,
    VaspToComputedEntryDrone,
)
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.util.testing import TEST_FILES_DIR, VASP_OUT_DIR

TEST_DIR = f"{TEST_FILES_DIR}/apps/borg"

MOL_TEST_DIR = f"{TEST_FILES_DIR}/io/gaussian"


class TestVaspToComputedEntryDrone(TestCase):
    def setUp(self):
        self.drone = VaspToComputedEntryDrone(data=["efermi"])
        self.structure_drone = VaspToComputedEntryDrone(inc_structure=True)

    def test_get_valid_paths(self):
        for path in os.walk(VASP_OUT_DIR):
            if path[0] == VASP_OUT_DIR:
                assert len(self.drone.get_valid_paths(path)) > 0

    def test_assimilate(self):
        """Test assimilate data from "vasprun.xml.xe.gz" file."""

        entry = self.drone.assimilate(f"{TEST_DIR}")

        for param in ("hubbards", "is_hubbard", "potcar_spec", "run_type"):
            assert param in entry.parameters
        assert entry.data["efermi"] == approx(-6.62148548)
        assert entry.reduced_formula == "Xe"
        assert entry.energy == approx(0.5559329)

        entry = self.structure_drone.assimilate(f"{TEST_DIR}")
        assert entry.reduced_formula == "Xe"
        assert entry.energy == approx(0.5559329)
        assert isinstance(entry, ComputedStructureEntry)
        assert entry.structure is not None

    def test_as_from_dict(self):
        dct = self.structure_drone.as_dict()
        drone = VaspToComputedEntryDrone.from_dict(dct)
        assert isinstance(drone, VaspToComputedEntryDrone)


class TestSimpleVaspToComputedEntryDrone(TestCase):
    def setUp(self):
        self.drone = SimpleVaspToComputedEntryDrone()
        self.structure_drone = SimpleVaspToComputedEntryDrone(inc_structure=True)

    def test_get_valid_paths(self):
        for path in os.walk(VASP_OUT_DIR):
            if path[0] == VASP_OUT_DIR:
                assert len(self.drone.get_valid_paths(path)) > 0

    def test_as_from_dict(self):
        dct = self.structure_drone.as_dict()
        drone = SimpleVaspToComputedEntryDrone.from_dict(dct)
        assert isinstance(drone, SimpleVaspToComputedEntryDrone)


class TestGaussianToComputedEntryDrone(TestCase):
    def setUp(self):
        self.drone = GaussianToComputedEntryDrone(data=["corrections"])
        self.structure_drone = GaussianToComputedEntryDrone(inc_structure=True)

    def test_get_valid_paths(self):
        for path in os.walk(MOL_TEST_DIR):
            if path[0] == MOL_TEST_DIR:
                assert len(self.drone.get_valid_paths(path)) > 0

    def test_assimilate(self):
        test_file = f"{MOL_TEST_DIR}/methane.log"
        entry = self.drone.assimilate(test_file)
        for param in [
            "functional",
            "basis_set",
            "charge",
            "spin_multiplicity",
            "route_parameters",
        ]:
            assert param in entry.parameters
        for param in ["corrections"]:
            assert param in entry.data

        assert entry.reduced_formula == "H4C"
        assert entry.energy == approx(-39.9768775602)
        entry = self.structure_drone.assimilate(test_file)
        assert entry.reduced_formula == "H4C"
        assert entry.energy == approx(-39.9768775602)
        assert isinstance(entry, ComputedStructureEntry)
        assert entry.structure is not None
        for param in ["properly_terminated", "stationary_type"]:
            assert param in entry.data

    def test_as_from_dict(self):
        dct = self.structure_drone.as_dict()
        drone = GaussianToComputedEntryDrone.from_dict(dct)
        assert isinstance(drone, GaussianToComputedEntryDrone)
