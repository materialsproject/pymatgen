from __future__ import annotations

import os
import subprocess
from typing import TYPE_CHECKING

from pymatgen.util.testing import VASP_IN_DIR, VASP_OUT_DIR

if TYPE_CHECKING:
    from pathlib import Path


def test_plot_xrd(cd_tmp_path: Path):
    subprocess.run(
        ["pmg", "plot", "--xrd", f"{VASP_IN_DIR}/POSCAR_Fe3O4", "--out_file", "xrd.png"],
        check=True,
    )
    assert os.path.isfile("xrd.png")
    assert os.path.getsize("xrd.png") > 1024


def test_plot_dos(cd_tmp_path: Path):
    subprocess.run(
        ["pmg", "plot", "--dos", f"{VASP_OUT_DIR}/vasprun_Li_no_projected.xml.gz", "--out_file", "dos.png"],
        check=True,
    )
    assert os.path.isfile("dos.png")
    assert os.path.getsize("dos.png") > 1024


def test_plot_chgint(cd_tmp_path: Path):
    subprocess.run(
        ["pmg", "plot", "--chgint", f"{VASP_OUT_DIR}/CHGCAR.Fe3O4.gz", "--out_file", "chg.png"],
        check=True,
    )
    assert os.path.isfile("chg.png")
    assert os.path.getsize("chg.png") > 1024
