from __future__ import annotations

import os
import subprocess

from pymatgen.util.testing import TEST_FILES_DIR


def test_pmg_analyze():
    subprocess.run(
        ["pmg", "analyze", f"{TEST_FILES_DIR}/io/vasp/fixtures/scan_relaxation"],
        check=True,
    )
    assert os.path.isfile("vasp_data.gz")
