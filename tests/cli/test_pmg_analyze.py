from __future__ import annotations

import os
import subprocess
from typing import TYPE_CHECKING

from pymatgen.util.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from pathlib import Path


def test_pmg_analyze(cd_tmp_path: Path):
    subprocess.run(
        ["pmg", "analyze", f"{TEST_FILES_DIR}/io/vasp/fixtures/scan_relaxation"],
        check=True,
    )
    assert os.path.isfile("vasp_data.gz")
