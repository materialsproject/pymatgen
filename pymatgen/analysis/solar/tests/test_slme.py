from __future__ import annotations

import os
import unittest
import warnings

import pytest

from pymatgen.analysis.solar.slme import optics, slme
from pymatgen.util.testing import PymatgenTest


class SolarTest(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_slme_from_vasprun(self):
        path = os.path.join(os.path.dirname(__file__), "vasprun.xml")
        en, abz, dirgap, indirgap = optics(path)
        abz = abz * 100.0
        eff = slme(en, abz, indirgap, indirgap, plot_current_voltage=False)
        assert eff == pytest.approx(27.728998512)


if __name__ == "__main__":
    unittest.main()
