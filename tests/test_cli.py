from __future__ import annotations

import os


def test_entrypoint(TEST_FILES_DIR):
    # exit_status = os.system("pmg analyze .")

    exit_status = os.system(f"pmg structure --convert --filenames {TEST_FILES_DIR}/Li2O.cif POSCAR.Li2O.test")
    assert exit_status == 0
    assert os.path.exists("POSCAR.Li2O.test")
    os.remove("POSCAR.Li2O.test")

    exit_status = os.system(f"pmg structure --symmetry 0.1 --filenames {TEST_FILES_DIR}/Li2O.cif")
    assert exit_status == 0

    exit_status = os.system(
        f"pmg structure --group element --filenames {TEST_FILES_DIR}/Li2O.cif {TEST_FILES_DIR}/Li.cif"
    )
    assert exit_status == 0
