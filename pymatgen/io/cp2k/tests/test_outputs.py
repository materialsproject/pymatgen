# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
from pathlib import Path

from pymatgen.io.cp2k.outputs import Cp2kOutput
from pymatgen.util.testing import PymatgenTest


class SetTest(PymatgenTest):
    def setUp(self):
        self.TEST_FILES_DIR = Path.joinpath(self.TEST_FILES_DIR, "cp2k")
        self.out = Cp2kOutput(Path.joinpath(self.TEST_FILES_DIR, "cp2k.out"), auto_load=True)

    def test_files(self):
        self.out.parse_files()
        self.assertEqual(len(self.out.filenames["PDOS"]), 2)

    def test(self):
        self.assertEqual(self.out.spin_polarized, False)
        self.assertEqual(self.out.completed, True)
        self.assertEqual(self.out.num_warnings, [[1]])
        self.assertEqual(self.out.run_type.upper(), "GEO_OPT")


if __name__ == "__main__":
    unittest.main()
