from __future__ import annotations

import subprocess
from typing import TYPE_CHECKING

from pymatgen.util.testing import VASP_IN_DIR

if TYPE_CHECKING:
    from pathlib import Path


def test_pmg_diff(cd_tmp_path: Path):
    subprocess.run(
        ["pmg", "diff", "--incar", f"{VASP_IN_DIR}/INCAR", f"{VASP_IN_DIR}/INCAR_2"],
        check=True,
    )
