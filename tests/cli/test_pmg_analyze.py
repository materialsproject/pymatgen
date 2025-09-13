from __future__ import annotations

import os

from pymatgen.cli import pmg
from pymatgen.util.testing import TEST_FILES_DIR


def test_pmg_analyze():
    pmg.main(
        ["analyze", f"{TEST_FILES_DIR}/io/vasp/fixtures/scan_relaxation"],
    )
    assert os.path.isfile("vasp_data.gz")
