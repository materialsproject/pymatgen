# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import json
import os
import unittest

from monty.json import MontyDecoder

from pymatgen.entries.exp_entries import ExpEntry

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "test_files")


class ExpEntryTest(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(test_dir, "Fe2O3_exp.json"), "r") as f:
            thermodata = json.load(f, cls=MontyDecoder)
        self.entry = ExpEntry("Fe2O3", thermodata)

    def test_energy(self):
        self.assertAlmostEqual(self.entry.energy, -825.5)

    def test_to_from_dict(self):
        d = self.entry.as_dict()
        e = ExpEntry.from_dict(d)
        self.assertAlmostEqual(e.energy, -825.5)

    def test_str(self):
        self.assertIsNotNone(str(self.entry))


if __name__ == "__main__":
    unittest.main()
