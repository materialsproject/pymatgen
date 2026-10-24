from __future__ import annotations

import os
from argparse import Namespace

import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure

from pymatgen.cli import pmg
from pymatgen.cli.pmg_plot import get_chgint_plot, get_dos_plot, plot
from pymatgen.util.testing import VASP_IN_DIR, VASP_OUT_DIR
from tests.testing import TEST_FILES_DIR

LIF_VASPRUN = str(TEST_FILES_DIR / "io" / "vasp" / "outputs" / "vasprun.LiF.xml.gz")  # projected DOS


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _dos_args(**kwargs) -> Namespace:
    return Namespace(**{"dos_file": LIF_VASPRUN, "site": None, "element": None, "orbital": None} | kwargs)


def _legend_labels(ax) -> list[str]:
    return [text.get_text() for text in ax.get_legend().get_texts()]


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


@pytest.mark.parametrize(
    ("kwargs", "labels"),
    [
        ({}, ["Total"]),
        ({"site": True}, ["Total", "Site 0 Li", "Site 1 F"]),
        ({"element": ["F, Na"]}, ["F"]),
        ({"orbital": True}, ["s", "p", "d"]),
    ],
)
def test_get_dos_plot_projections(kwargs, labels):
    ax = get_dos_plot(_dos_args(**kwargs))
    assert sorted(_legend_labels(ax)) == sorted(labels)


def test_get_chgint_plot_explicit_indices():
    args = Namespace(chgcar_file=f"{VASP_OUT_DIR}/CHGCAR.Fe3O4.gz", inds=["0,3"], radius=2.0)
    ax = get_chgint_plot(args)
    assert _legend_labels(ax) == ["Atom 0 - Fe", "Atom 3 - Fe"]
    assert ax.get_xlabel() == "Radius (A)"


@pytest.fixture
def shown_figures(monkeypatch) -> list[Figure]:
    """Record Figure.show() calls instead of opening a window."""
    shown: list[Figure] = []

    def fake_show(self):
        shown.append(self)

    monkeypatch.setattr(Figure, "show", fake_show)
    return shown


def test_plot_shows_figure_without_out_file(shown_figures):
    shown = shown_figures
    plot(Namespace(chgcar_file=None, xrd_structure_file=None, **vars(_dos_args()), out_file=None))
    assert len(shown) == 1
    assert isinstance(shown[0], Figure)


def test_plot_nothing_selected(shown_figures):
    shown = shown_figures
    plot(Namespace(chgcar_file=None, xrd_structure_file=None, dos_file=None, out_file="unused.png"))
    assert shown == []
    assert not os.path.isfile("unused.png")
