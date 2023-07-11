from __future__ import annotations

import json
import os
import unittest

from monty.json import MontyDecoder
from monty.io import zopen

from pytest import approx

from pymatgen.entries.exp_entries import ExpEntry
from pymatgen.util.testing import PymatgenTest


class ExpEntryTest(unittest.TestCase):
    def setUp(self):
        with zopen(os.path.join(PymatgenTest.TEST_FILES_DIR, "Fe2O3_exp.json.gz"), 'r') as f:
            thermo_data = json.loads(f.read(), cls=MontyDecoder)
        self.entry = ExpEntry("Fe2O3", thermo_data)

    def test_energy(self):
        assert self.entry.energy == approx(-825.5)

    def test_to_from_dict(self):
        d = self.entry.as_dict()
        e = ExpEntry.from_dict(d)
        assert e.energy == approx(-825.5)

    def test_str(self):
        assert str(self.entry) is not None
