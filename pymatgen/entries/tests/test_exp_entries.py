# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


from __future__ import annotations

import json
import os
import unittest

from monty.json import MontyDecoder

from pymatgen.entries.exp_entries import ExpEntry
from pymatgen.util.testing import PymatgenTest


class ExpEntryTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(PymatgenTest.TEST_FILES_DIR, "Fe2O3_exp.json")) as f:
            thermodata = json.load(f, cls=MontyDecoder)
        self.entry = ExpEntry("Fe2O3", thermodata)

    def test_energy(self):
        assert round(abs(self.entry.energy - -825.5), 7) == 0

    def test_to_from_dict(self):
        d = self.entry.as_dict()
        e = ExpEntry.from_dict(d)
        assert round(abs(e.energy - -825.5), 7) == 0

    def test_str(self):
        assert str(self.entry) is not None


if __name__ == "__main__":
    unittest.main()
