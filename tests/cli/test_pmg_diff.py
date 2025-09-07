from __future__ import annotations

import subprocess

from pymatgen.util.testing import VASP_IN_DIR


def test_pmg_diff():
    subprocess.run(
        ["pmg", "diff", "--incar", f"{VASP_IN_DIR}/INCAR", f"{VASP_IN_DIR}/INCAR_2"],
        check=True,
    )
