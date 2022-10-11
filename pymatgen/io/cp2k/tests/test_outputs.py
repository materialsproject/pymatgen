# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
from pathlib import Path

from pymatgen.io.cp2k.outputs import Cp2kOutput
from pymatgen.util.testing import PymatgenTest


# TODO More comprehensive testing
class SetTest(PymatgenTest):
    def setUp(self):
        self.TEST_FILES_DIR = Path.joinpath(self.TEST_FILES_DIR, "cp2k")
        self.out = Cp2kOutput(Path.joinpath(self.TEST_FILES_DIR, "cp2k.out"), auto_load=True)

    def test_files(self):
        self.out.parse_files()
        self.assertEqual(len(self.out.filenames["PDOS"]), 1)

    def test_basic(self):
        self.assertEqual(self.out.spin_polarized, True)
        self.assertEqual(self.out.completed, True)
        self.assertEqual(self.out.num_warnings, [[2]])
        self.assertEqual(self.out.run_type.upper(), "ENERGY_FORCE")
        self.assertEqual(self.out.final_energy, -197.40000341992783)

    def test_band(self):
        self.assertTrue(self.out.band_structure)
        self.assertEqual(self.out.band_structure.get_band_gap().get("energy"), 0.27940141999999923)


if __name__ == "__main__":
    unittest.main()
