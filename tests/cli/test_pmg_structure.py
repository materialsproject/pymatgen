from __future__ import annotations

import os
import subprocess

from pymatgen.util.testing import TEST_FILES_DIR


def test_pmg_structure():
    subprocess.run(
        [
            "pmg",
            "structure",
            "--convert",
            "--filenames",
            f"{TEST_FILES_DIR}/cif/Li2O.cif",
            "POSCAR_Li2O_test",
        ],
        check=True,
    )
    assert os.path.isfile("POSCAR_Li2O_test"), "Output file 'POSCAR_Li2O_test' not found"

    subprocess.run(
        [
            "pmg",
            "structure",
            "--symmetry",
            "0.1",
            "--filenames",
            f"{TEST_FILES_DIR}/cif/Li2O.cif",
        ],
        check=True,
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
        [
            "pmg",
            "structure",
            "--localenv",
            "Li-O=3",
            "--filenames",
            f"{TEST_FILES_DIR}/cif/Li2O.cif",
        ],
        check=True,
    )
