from __future__ import annotations

import os

from pytest import approx

from pymatgen.analysis.topological.spillage import SOCSpillage
from pymatgen.util.testing import PymatgenTest


class TestSolar(PymatgenTest):
    def test_spillage_from_vasprun(self):
        wf_noso = os.path.join(os.path.dirname(__file__), "WAVECAR-NonSOC")
        wf_so = os.path.join(os.path.dirname(__file__), "WAVECAR-SOC")
        # JVASP-1044
        gamma_max = SOCSpillage(wf_noso=wf_noso, wf_so=wf_so).overlap_so_spinpol()

        assert gamma_max == approx(1.3634111271008775, abs=1e-5)
