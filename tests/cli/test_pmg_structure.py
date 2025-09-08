from __future__ import annotations

import os

from pymatgen.cli import pmg
from pymatgen.util.testing import TEST_FILES_DIR


def test_pmg_structure():
    pmg.main(
        [
            "structure",
            "--convert",
            "--filenames",
            f"{TEST_FILES_DIR}/cif/Li2O.cif",
            "POSCAR_Li2O_test",
        ],
    )
    assert os.path.isfile("POSCAR_Li2O_test"), "Output file 'POSCAR_Li2O_test' not found"

    pmg.main(
        [
            "structure",
            "--symmetry",
            "0.1",
            "--filenames",
            f"{TEST_FILES_DIR}/cif/Li2O.cif",
        ],
    )

    pmg.main(
        [
            "structure",
            "--group",
            "element",
            "--filenames",
            f"{TEST_FILES_DIR}/cif/Li2O.cif",
            f"{TEST_FILES_DIR}/cif/Li.cif",
        ],
    )

    pmg.main(
        [
            "structure",
            "--localenv",
            "Li-O=3",
            "--filenames",
            f"{TEST_FILES_DIR}/cif/Li2O.cif",
        ],
    )
