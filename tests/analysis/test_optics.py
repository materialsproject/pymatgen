from __future__ import annotations

import pytest

from pymatgen.analysis.optics import DielectricAnalysis
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/io/vasp/outputs"


def test_dielectric_analysis():
    analysis = DielectricAnalysis.from_vasprun(Vasprun(f"{TEST_DIR}/vasprun.dielectric_6.0.8.xml.gz"))
    assert analysis.n.shape == (2000, 6)
    assert analysis.n[0, 1] == pytest.approx(1.96636721)
