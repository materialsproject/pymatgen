# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest
import warnings

from pymatgen.apps.borg.hive import (
    GaussianToComputedEntryDrone,
    SimpleVaspToComputedEntryDrone,
    VaspToComputedEntryDrone,
)
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.util.testing import PymatgenTest


class VaspToComputedEntryDroneTest(unittest.TestCase):
    def setUp(self):
        self.drone = VaspToComputedEntryDrone(data=["efermi"])
        self.structure_drone = VaspToComputedEntryDrone(True)

    def test_get_valid_paths(self):
        for path in os.walk(PymatgenTest.TEST_FILES_DIR):
            if path[0] == PymatgenTest.TEST_FILES_DIR:
                self.assertTrue(len(self.drone.get_valid_paths(path)) > 0)

    def test_assimilate(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            entry = self.drone.assimilate(PymatgenTest.TEST_FILES_DIR)
            for p in ["hubbards", "is_hubbard", "potcar_spec", "run_type"]:
                self.assertIn(p, entry.parameters)
            self.assertAlmostEqual(entry.data["efermi"], -6.62148548)
            self.assertEqual(entry.composition.reduced_formula, "Xe")
            self.assertAlmostEqual(entry.energy, 0.5559329)
            entry = self.structure_drone.assimilate(PymatgenTest.TEST_FILES_DIR)
            self.assertEqual(entry.composition.reduced_formula, "Xe")
            self.assertAlmostEqual(entry.energy, 0.5559329)
            self.assertIsInstance(entry, ComputedStructureEntry)
            self.assertIsNotNone(entry.structure)
            # self.assertEqual(len(entry.parameters["history"]), 2)

    def tearDown(self):
        warnings.simplefilter("default")

    def test_to_from_dict(self):
        d = self.structure_drone.as_dict()
        drone = VaspToComputedEntryDrone.from_dict(d)
        self.assertEqual(type(drone), VaspToComputedEntryDrone)


class SimpleVaspToComputedEntryDroneTest(unittest.TestCase):
    def setUp(self):
        self.drone = SimpleVaspToComputedEntryDrone()
        self.structure_drone = SimpleVaspToComputedEntryDrone(True)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_get_valid_paths(self):
        for path in os.walk(PymatgenTest.TEST_FILES_DIR):
            if path[0] == PymatgenTest.TEST_FILES_DIR:
                self.assertTrue(len(self.drone.get_valid_paths(path)) > 0)

    def test_to_from_dict(self):
        d = self.structure_drone.as_dict()
        drone = SimpleVaspToComputedEntryDrone.from_dict(d)
        self.assertEqual(type(drone), SimpleVaspToComputedEntryDrone)


class GaussianToComputedEntryDroneTest(unittest.TestCase):
    def setUp(self):
        self.drone = GaussianToComputedEntryDrone(data=["corrections"])
        self.structure_drone = GaussianToComputedEntryDrone(True)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_get_valid_paths(self):
        for path in os.walk(os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules")):
            if path[0] == os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules"):
                self.assertTrue(len(self.drone.get_valid_paths(path)) > 0)

    def test_assimilate(self):
        test_file = os.path.join(PymatgenTest.TEST_FILES_DIR, "molecules", "methane.log")
        entry = self.drone.assimilate(test_file)
        for p in [
            "functional",
            "basis_set",
            "charge",
            "spin_multiplicity",
            "route_parameters",
        ]:
            self.assertIn(p, entry.parameters)
        for p in ["corrections"]:
            self.assertIn(p, entry.data)

        self.assertEqual(entry.composition.reduced_formula, "H4C")
        self.assertAlmostEqual(entry.energy, -39.9768775602)
        entry = self.structure_drone.assimilate(test_file)
        self.assertEqual(entry.composition.reduced_formula, "H4C")
        self.assertAlmostEqual(entry.energy, -39.9768775602)
        self.assertIsInstance(entry, ComputedStructureEntry)
        self.assertIsNotNone(entry.structure)
        for p in ["properly_terminated", "stationary_type"]:
            self.assertIn(p, entry.data)

    def test_to_from_dict(self):
        d = self.structure_drone.as_dict()
        drone = GaussianToComputedEntryDrone.from_dict(d)
        self.assertEqual(type(drone), GaussianToComputedEntryDrone)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
