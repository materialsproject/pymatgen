from __future__ import annotations

import os

import pytest

from pymatgen.cli import pmg
from pymatgen.util.testing import VASP_IN_DIR, VASP_OUT_DIR


def test_plot_xrd():
    pmg.main(
        [
            "plot",
            "--xrd",
            f"{VASP_IN_DIR}/POSCAR_Fe3O4",
            "--out_file",
            "xrd.png",
        ],
    )
    assert os.path.isfile("xrd.png")
    assert os.path.getsize("xrd.png") > 1024


def test_plot_dos():
    pmg.main(
        [
            "plot",
            "--dos",
            f"{VASP_OUT_DIR}/vasprun_Li_no_projected.xml.gz",
            "--out_file",
            "dos.png",
        ],
    )
    assert os.path.isfile("dos.png")
    assert os.path.getsize("dos.png") > 1024


def test_plot_chgint():
    pmg.main(
        [
            "plot",
            "--chgint",
            f"{VASP_OUT_DIR}/CHGCAR.Fe3O4.gz",
            "--out_file",
            "chg.png",
        ],
    )
    assert os.path.isfile("chg.png")
    assert os.path.getsize("chg.png") > 1024


def test_plot_wrong_arg(capsys):
    with pytest.raises(SystemExit) as exc_info:
        pmg.main(
            ["plot", "--wrong", f"{VASP_OUT_DIR}/CHGCAR.Fe3O4.gz"],
        )

    assert exc_info.value.code == 2
    captured = capsys.readouterr()
    assert "one of the arguments -d/--dos -c/--chgint -x/--xrd is required" in captured.err
