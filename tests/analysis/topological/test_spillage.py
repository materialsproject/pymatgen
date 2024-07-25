from __future__ import annotations

from pymatgen.analysis.topological.spillage import SOCSpillage
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest
from pytest import approx

TEST_DIR = f"{TEST_FILES_DIR}/analysis/topological"


class TestSolar(PymatgenTest):
    def test_spillage_from_vasprun(self):
        wf_noso = f"{TEST_DIR}/WAVECAR-NonSOC"
        wf_so = f"{TEST_DIR}/WAVECAR-SOC"
        # JVASP-1044
        gamma_max = SOCSpillage(wf_noso=wf_noso, wf_so=wf_so).overlap_so_spinpol()

        assert gamma_max == approx(1.3634111271008775, abs=1e-5)
