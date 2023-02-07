# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import unittest

from pytest import approx

from pymatgen.analysis.hhi import HHIModel


class HHIModelTest(unittest.TestCase):
    def test_hhi(self):
        hhi = HHIModel()
        assert hhi.get_hhi("He") == (3200, 3900)
        assert hhi.get_hhi_production("He") == 3200
        assert hhi.get_hhi_reserve("He") == 3900

        assert hhi.get_hhi_production("Li2O") == approx(1614.96, abs=1e-1)
        assert hhi.get_hhi_reserve("Li2O") == approx(2218.90, abs=1e-1)

        assert hhi.get_hhi_designation(1400) == "low"
        assert hhi.get_hhi_designation(1800) == "medium"
        assert hhi.get_hhi_designation(3000) == "high"
        assert hhi.get_hhi_designation(None) is None


if __name__ == "__main__":
    unittest.main()
