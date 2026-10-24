from __future__ import annotations

import os

import pytest

from pymatgen.cli import pmg
from pymatgen.core import Structure
from tests.testing import TEST_FILES_DIR

LI2O_CIF = str(TEST_FILES_DIR / "cif" / "Li2O.cif")  # oxidation-state decorated
LI2O_POSCAR = str(TEST_FILES_DIR / "io" / "vasp" / "inputs" / "POSCAR_Li2O")  # 2x2x2 supercell, no oxi states
LI_CIF = str(TEST_FILES_DIR / "cif" / "Li.cif")


@pytest.fixture(autouse=True)
def _in_tmp_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)


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


@pytest.mark.parametrize(
    ("out_name", "n_sites"),
    [("Li2O.cif", 12), ("POSCAR_prim", 3), ("Li2O_prim.json", 3)],
)
def test_convert(tmp_path, out_name, n_sites):
    """'prim' in the output name triggers a primitive-cell reduction on read."""
    pmg.main(["structure", "-c", "-f", LI2O_CIF, out_name])

    struct = Structure.from_file(tmp_path / out_name)
    assert len(struct) == n_sites
    assert struct.reduced_formula == "Li2O"
    assert struct.matches(Structure.from_file(LI2O_CIF))


def test_convert_extra_filenames(tmp_path, capsys):
    pmg.main(["structure", "-c", "-f", LI2O_CIF, "POSCAR", "ignored.cif"])

    assert "File format conversion takes in only two filenames." in capsys.readouterr().out
    # Conversion still proceeds with the first two filenames.
    assert (tmp_path / "POSCAR").is_file()
    assert not (tmp_path / "ignored.cif").exists()


def test_symmetry(capsys):
    pmg.main(["structure", "-s", "0.1", "-f", LI2O_CIF, LI_CIF])

    lines = capsys.readouterr().out.splitlines()
    assert lines[0].split() == ["Filename", "Int", "Symbol", "Int", "number", "Hall"]
    assert lines[2].split()[1:] == ["Fm-3m", "225", "-F", "4", "2", "3"]
    assert lines[3].split()[1:] == ["Im-3m", "229", "-I", "4", "2", "3"]


@pytest.mark.parametrize(
    ("mode", "expected_groups"),
    [
        # Element mode ignores oxidation states: the decorated CIF and the plain POSCAR match.
        ("element", [[LI2O_CIF, LI2O_POSCAR], [LI_CIF]]),
        # Species mode distinguishes Li+/O2- from Li/O.
        ("species", [[LI2O_POSCAR], [LI_CIF], [LI2O_CIF]]),
    ],
)
def test_group(capsys, mode, expected_groups):
    pmg.main(["structure", "-g", mode, "-f", LI2O_CIF, LI2O_POSCAR, LI_CIF])

    groups: list[list[str]] = []
    for line in capsys.readouterr().out.splitlines():
        if line.startswith("Group "):
            groups.append([])
        elif line.startswith("- "):
            groups[-1].append(line[2:].rsplit(" (", 1)[0])
    assert sorted(groups) == sorted(expected_groups)


def test_group_needs_two_files():
    with pytest.raises(SystemExit, match="You need more than one structure to compare!"):
        pmg.main(["structure", "-g", "element", "-f", LI2O_CIF])


def test_group_unreadable_file(tmp_path, capsys):
    missing = str(tmp_path / "does_not_exist.cif")
    with pytest.raises(SystemExit):
        pmg.main(["structure", "-g", "species", "-f", LI2O_CIF, missing])
    assert "Error converting file. Are they in the right format?" in capsys.readouterr().out


def test_localenv(capsys):
    pmg.main(["structure", "-l", "Li-O=2.5", "O-Li=2.5", "-f", LI2O_CIF])

    out = capsys.readouterr().out
    assert f"Analyzing {LI2O_CIF}..." in out
    rows = [line.split(None, 3) for line in out.splitlines()[3:]]
    assert len(rows) == 12
    li_rows = [row for row in rows if row[1] == "Li"]
    o_rows = [row for row in rows if row[1] == "O"]
    assert len(li_rows) == 8
    assert len(o_rows) == 4
    # Antifluorite: Li is tetrahedrally and O cubically coordinated at a*sqrt(3)/4.
    assert all(row[2] == "O" and row[3] == ", ".join(["1.996"] * 4) for row in li_rows)
    assert all(row[2] == "Li" and row[3] == ", ".join(["1.996"] * 8) for row in o_rows)
