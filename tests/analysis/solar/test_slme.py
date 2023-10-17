from __future__ import annotations

import os

from pytest import approx

from pymatgen.analysis.solar.slme import optics, slme
from pymatgen.util.testing import PymatgenTest


class TestSolar(PymatgenTest):
    def test_slme_from_vasprun(self):
        en, abz, dir_gap, indir_gap = optics(f"{os.path.dirname(__file__)}/vasprun.xml")
        abz = abz * 100.0
        eff = slme(en, abz, indir_gap, indir_gap, plot_current_voltage=False)
        assert eff == approx(27.729, abs=1e-2)
