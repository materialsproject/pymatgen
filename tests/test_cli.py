from __future__ import annotations

import os
from typing import TYPE_CHECKING

from pymatgen.util.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from pathlib import Path

    from pytest import MonkeyPatch


def test_entrypoint(tmp_path: Path, monkeypatch: MonkeyPatch):
    monkeypatch.chdir(tmp_path)

    exit_status = os.system(f"pmg analyze {TEST_FILES_DIR}/scan_relaxation")
    assert exit_status == 0
    assert os.path.exists("vasp_data.gz")

    exit_status = os.system(f"pmg structure --convert --filenames {TEST_FILES_DIR}/Li2O.cif POSCAR.Li2O.test")
    assert exit_status == 0
    assert os.path.exists("POSCAR.Li2O.test")

    exit_status = os.system(f"pmg structure --symmetry 0.1 --filenames {TEST_FILES_DIR}/Li2O.cif")
    assert exit_status == 0

    exit_status = os.system(
        f"pmg structure --group element --filenames {TEST_FILES_DIR}/Li2O.cif {TEST_FILES_DIR}/Li.cif"
    )
    assert exit_status == 0

    exit_status = os.system(f"pmg structure --localenv Li-O=3 --filenames {TEST_FILES_DIR}/Li2O.cif")
    assert exit_status == 0

    exit_status = os.system(f"pmg diff --incar {TEST_FILES_DIR}/INCAR {TEST_FILES_DIR}/INCAR.2")
    assert exit_status == 0
