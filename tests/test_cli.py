from __future__ import annotations

import os
from typing import TYPE_CHECKING

import pytest

from pymatgen.util.testing import TEST_FILES_DIR, VASP_IN_DIR

if TYPE_CHECKING:
    from pathlib import Path

    from pytest import MonkeyPatch


@pytest.fixture()
def cd_tmp_path(tmp_path: Path, monkeypatch: MonkeyPatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_pmg_analyze(cd_tmp_path: Path):
    exit_status = os.system(f"pmg analyze {TEST_FILES_DIR}/io/vasp/fixtures/scan_relaxation")
    assert exit_status == 0
    assert os.path.isfile("vasp_data.gz")


def test_pmg_structure(cd_tmp_path: Path):
    exit_status = os.system(f"pmg structure --convert --filenames {TEST_FILES_DIR}/cif/Li2O.cif POSCAR_Li2O_test")
    assert exit_status == 0
    assert os.path.isfile("POSCAR_Li2O_test")

    exit_status = os.system(f"pmg structure --symmetry 0.1 --filenames {TEST_FILES_DIR}/cif/Li2O.cif")
    assert exit_status == 0

    exit_status = os.system(
        f"pmg structure --group element --filenames {TEST_FILES_DIR}/cif/Li2O.cif {TEST_FILES_DIR}/cif/Li.cif"
    )
    assert exit_status == 0

    exit_status = os.system(f"pmg structure --localenv Li-O=3 --filenames {TEST_FILES_DIR}/cif/Li2O.cif")
    assert exit_status == 0


def test_pmg_diff(cd_tmp_path: Path):
    exit_status = os.system(f"pmg diff --incar {VASP_IN_DIR}/INCAR {VASP_IN_DIR}/INCAR_2")
    assert exit_status == 0
