from __future__ import annotations

import os
import subprocess
from typing import TYPE_CHECKING

import pytest

from pymatgen.util.testing import TEST_FILES_DIR, VASP_IN_DIR

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def cd_tmp_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_pmg_analyze(cd_tmp_path: Path):
    subprocess.run(
        ["pmg", "analyze", f"{TEST_FILES_DIR}/io/vasp/fixtures/scan_relaxation"],
        check=True,
    )
    assert os.path.isfile("vasp_data.gz")


def test_pmg_structure(cd_tmp_path: Path):
    subprocess.run(
        ["pmg", "structure", "--convert", "--filenames", f"{TEST_FILES_DIR}/cif/Li2O.cif", "POSCAR_Li2O_test"],
        check=True,
    )
    assert os.path.isfile("POSCAR_Li2O_test"), "Output file 'POSCAR_Li2O_test' not found"

    subprocess.run(
        ["pmg", "structure", "--symmetry", "0.1", "--filenames", f"{TEST_FILES_DIR}/cif/Li2O.cif"], check=True
    )

    subprocess.run(
        [
            "pmg",
            "structure",
            "--group",
            "element",
            "--filenames",
            f"{TEST_FILES_DIR}/cif/Li2O.cif",
            f"{TEST_FILES_DIR}/cif/Li.cif",
        ],
        check=True,
    )

    subprocess.run(
        ["pmg", "structure", "--localenv", "Li-O=3", "--filenames", f"{TEST_FILES_DIR}/cif/Li2O.cif"], check=True
    )


def test_pmg_diff(cd_tmp_path: Path):
    subprocess.run(
        ["pmg", "diff", "--incar", f"{VASP_IN_DIR}/INCAR", f"{VASP_IN_DIR}/INCAR_2"],
        check=True,
    )
