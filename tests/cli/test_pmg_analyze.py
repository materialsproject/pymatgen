from __future__ import annotations

import logging
import os
import shutil
from argparse import Namespace
from typing import TYPE_CHECKING

import pytest
from pytest import approx

from pymatgen.cli import pmg
from pymatgen.cli.pmg_analyze import SAVE_FILE, analyze, get_energies, get_magnetizations
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.outputs import Vasprun
from tests.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from pathlib import Path

VASP_FIXTURES = TEST_FILES_DIR / "io" / "vasp" / "fixtures"
VASP_IN_DIR = TEST_FILES_DIR / "io" / "vasp" / "inputs"
VASP_OUT_DIR = TEST_FILES_DIR / "io" / "vasp" / "outputs"


@pytest.fixture
def serial(monkeypatch: pytest.MonkeyPatch) -> None:
    """Force BorgQueen.serial_assimilate (no process pool) for speed and determinism."""
    monkeypatch.setattr("pymatgen.cli.pmg_analyze.multiprocessing.cpu_count", lambda: 1)


@pytest.fixture
def runs_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Two VASP run dirs: 'b_li' (Li, GGA) and 'a_lih' (LiH, r2SCAN). Cwd is tmp_path."""
    monkeypatch.chdir(tmp_path)
    root = tmp_path / "runs"
    for name, src in (("b_li", "relaxation"), ("a_lih", "scan_relaxation")):
        (root / name).mkdir(parents=True)
        shutil.copy(VASP_FIXTURES / src / "vasprun.xml.gz", root / name / "vasprun.xml.gz")
    return root


def _parse_tsv(out: str) -> list[list[str]]:
    """Parse the tsv table printed by `pmg analyze -f tsv` (header excluded)."""
    lines = [line for line in out.splitlines() if "\t" in line]
    return [[cell.strip() for cell in line.split("\t")] for line in lines[1:]]


def test_pmg_analyze():
    pmg.main(
        ["analyze", f"{TEST_FILES_DIR}/io/vasp/fixtures/scan_relaxation"],
    )
    assert os.path.isfile("vasp_data.gz")


def test_pmg_analyze_all_magnetizations(capsys):
    assert (
        pmg.main(
            ["analyze", "--mag", "All", f"{TEST_FILES_DIR}/io/vasp/fixtures/scan_relaxation"],
        )
        == 0
    )
    assert "OUTCAR" in capsys.readouterr().out


@pytest.mark.parametrize(
    ("sort", "sort_col"),
    [("energy_per_atom", 3), ("filename", 0)],
)
def test_get_energies_sorted(runs_dir, serial, capsys, sort, sort_col):
    assert pmg.main(["analyze", "-e", "-r", "-f", "tsv", "-s", sort, str(runs_dir)]) == 0

    out = capsys.readouterr().out
    rows = _parse_tsv(out)
    assert len(rows) == 2
    assert sorted(row[1] for row in rows) == ["Li1", "Li1H1"]
    keys = [float(row[sort_col]) if sort_col == 3 else row[sort_col] for row in rows]
    assert keys == sorted(keys)
    assert f"Analysis results saved to {SAVE_FILE}" in out
    assert os.path.isfile(SAVE_FILE)


def test_get_energies_values(runs_dir, serial, capsys):
    assert get_energies(str(runs_dir), reanalyze=True, verbose=False, quick=False, sort="filename", fmt="tsv") == 0
    rows = _parse_tsv(capsys.readouterr().out)
    assert [os.path.basename(os.path.dirname(row[0])) for row in rows] == ["a_lih", "b_li"]

    for row, name in zip(rows, ("a_lih", "b_li"), strict=True):
        vasprun = Vasprun(runs_dir / name / "vasprun.xml.gz", parse_dos=False, parse_eigen=False)
        n_atoms = len(vasprun.final_structure)
        assert float(row[2]) == approx(vasprun.final_energy, abs=1e-4)
        assert float(row[3]) == approx(vasprun.final_energy / n_atoms, abs=1e-4)
        vol_chg = 100 * (vasprun.final_structure.volume / vasprun.initial_structure.volume - 1)
        assert float(row[4]) == approx(vol_chg, abs=1e-2)


def test_get_energies_reuses_saved_data(runs_dir, serial, capsys):
    assert pmg.main(["analyze", str(runs_dir)]) == 0
    assert "Analysis results saved" in capsys.readouterr().out

    # Without -r the saved data is loaded, even if the run dirs are gone.
    shutil.rmtree(runs_dir)
    assert pmg.main(["analyze", str(runs_dir)]) == 0
    out = capsys.readouterr().out
    assert f"Using previously assimilated data from {SAVE_FILE}" in out
    assert "a_lih" in out
    assert "b_li" in out

    # With -r the (now empty) dir is re-analyzed and the stale save file removed.
    runs_dir.mkdir()
    assert pmg.main(["analyze", "-r", str(runs_dir)]) == 0
    assert "No valid vasp run found." in capsys.readouterr().out
    assert not os.path.isfile(SAVE_FILE)


def test_get_energies_verbose(runs_dir, serial, capsys):
    root_logger = logging.getLogger()
    handlers, level = root_logger.handlers[:], root_logger.level
    try:
        assert pmg.main(["analyze", "-v", "-r", str(runs_dir)]) == 0
    finally:
        root_logger.handlers[:] = handlers
        root_logger.setLevel(level)
    assert "a_lih" in capsys.readouterr().out


def test_get_energies_quick(tmp_path, monkeypatch, serial, capsys):
    """Quick mode parses INCAR/POTCAR/POSCAR/CONTCAR/OSZICAR instead of vasprun.xml."""
    monkeypatch.chdir(tmp_path)
    run_dir = tmp_path / "runs" / "lfp"
    run_dir.mkdir(parents=True)
    structure = Poscar.from_file(VASP_IN_DIR / "POSCAR_LiFePO4", check_for_potcar=False).structure
    order = {"Li": 0, "Fe": 1, "P": 2, "O": 3}  # match POTCAR_LiFePO4 ordering
    structure = structure.get_sorted_structure(key=lambda site: order[site.specie.symbol])
    Poscar(structure).write_file(run_dir / "POSCAR")
    Poscar(structure).write_file(run_dir / "CONTCAR")
    Incar({"ENCUT": 520}).write_file(run_dir / "INCAR")
    shutil.copy(VASP_IN_DIR / "POTCAR_LiFePO4.gz", run_dir / "POTCAR.gz")
    shutil.copy(VASP_OUT_DIR / "OSZICAR", run_dir / "OSZICAR")

    assert pmg.main(["analyze", "-q", "-f", "tsv", str(tmp_path / "runs")]) == 0

    rows = _parse_tsv(capsys.readouterr().out)
    assert len(rows) == 1
    assert os.path.basename(rows[0][0]) == "lfp"
    assert rows[0][1] == "Li4Fe4P4O16"
    assert rows[0][4] == "NA"


@pytest.fixture
def outcar_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.chdir(tmp_path)
    root = tmp_path / "outcars"
    for name, src in (("a", "OUTCAR.gz"), ("b", "OUTCAR.Al")):
        (root / name).mkdir(parents=True)
        shutil.copy(VASP_OUT_DIR / src, root / name / src)
    # Unparsable OUTCARs are silently skipped.
    (root / "c").mkdir()
    (root / "c" / "OUTCAR").write_text("not an OUTCAR", encoding="utf-8")
    return root


def test_get_magnetizations_all(outcar_dir, capsys):
    assert pmg.main(["analyze", "-m", "All", str(outcar_dir)]) == 0

    lines = capsys.readouterr().out.splitlines()
    header = lines[0].split()
    assert header == ["Filename", *map(str, range(7))]
    rows = {line.split()[0]: line.split()[1:] for line in lines[2:]}
    assert len(rows) == 2
    gz_row = next(vals for name, vals in rows.items() if name.endswith("OUTCAR.gz"))
    al_row = next(vals for name, vals in rows.items() if name.endswith("OUTCAR.Al"))
    assert gz_row == ["0.005", "0.813", "0.813", "-0.112", "-0.162", "-0.112", "-0.162"]
    assert al_row == ["0"]  # tabulate renders 0.0 as "0"; remaining cells are blank padding


def test_get_magnetizations_range(outcar_dir, capsys):
    # "1-2" selects 0-based ion indices 1 and 2; OUTCAR.Al (1 ion) raises IndexError and is skipped.
    assert analyze(Namespace(directories=[str(outcar_dir)], get_energies=False, ion_list=["1-2"])) == 0

    lines = capsys.readouterr().out.splitlines()
    assert lines[0].split() == ["Filename", "0", "1"]
    assert len(lines) == 3
    assert lines[2].split()[1:] == ["0.813", "0.813"]


def test_get_magnetizations_empty_dir(tmp_path, capsys):
    assert get_magnetizations(str(tmp_path), None) == 0
    assert capsys.readouterr().out.split() == ["Filename", "----------"]


def test_analyze_no_directories():
    args = Namespace(directories=[], get_energies=True, ion_list=None)
    assert analyze(args) == -1
