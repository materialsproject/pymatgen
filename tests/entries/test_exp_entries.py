from __future__ import annotations

import json
from unittest import TestCase

from monty.json import MontyDecoder
from pytest import approx

from pymatgen.entries.exp_entries import ExpEntry
from pymatgen.util.testing import TEST_FILES_DIR


class TestExpEntry(TestCase):
    def setUp(self):
        with open(f"{TEST_FILES_DIR}/entries/Fe2O3_exp.json") as file:
            thermo_data = json.load(file, cls=MontyDecoder)
        self.entry = ExpEntry("Fe2O3", thermo_data)

    def test_energy(self):
        assert self.entry.energy == approx(-825.5)

    def test_as_from_dict(self):
        dct = self.entry.as_dict()
        entry = ExpEntry.from_dict(dct)
        assert entry.energy == approx(-825.5)

    def test_str(self):
        assert str(self.entry) == "ExpEntry Fe2 O3, Energy = -825.5000"
