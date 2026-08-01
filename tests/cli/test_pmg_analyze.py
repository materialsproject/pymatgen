from __future__ import annotations

import os

from pymatgen.cli import pmg
from tests.testing import TEST_FILES_DIR


def test_pmg_analyze():
    pmg.main(
        ["analyze", f"{TEST_FILES_DIR}/io/vasp/fixtures/scan_relaxation"],
    )
    assert os.path.isfile("vasp_data.gz")


def test_pmg_analyze_all_magnetizations(capsys):
    assert (
        pmg.main(
            ["analyze", "--mag", "All", f"{TEST_FILES_DIR}/io/vasp/fixtures/scan_relaxation"],
        )
        == 0
    )
    assert "OUTCAR" in capsys.readouterr().out
