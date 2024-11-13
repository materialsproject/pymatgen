from __future__ import annotations

from pytest import approx

from pymatgen.analysis.solar.slme import optics, slme
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/analysis/solar"


class TestSolar(PymatgenTest):
    def test_slme_from_vasprun(self):
        en, abz, dir_gap, indir_gap = optics(f"{TEST_DIR}/vasprun.xml")
        abz *= 100.0
        eff = slme(en, abz, indir_gap, indir_gap, plot_current_voltage=False)
        assert eff == approx(27.729, abs=1e-2)
        assert dir_gap == approx(0.85389999, abs=1e-6)
