from __future__ import annotations

import gzip
import shutil
import sys
from typing import TYPE_CHECKING

import pytest
from monty.serialization import dumpfn, loadfn

from pymatgen.cli import pmg, pmg_config

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def settings_files(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    """Redirect SETTINGS_FILE / OLD_SETTINGS_FILE to temp paths so the real ~/.config/.pmgrc.yaml
    is never touched.
    """
    new_file = tmp_path / "config" / ".pmgrc.yaml"
    new_file.parent.mkdir()
    old_file = tmp_path / ".pmgrc.yaml"
    monkeypatch.setattr(pmg_config, "SETTINGS_FILE", str(new_file))
    monkeypatch.setattr(pmg_config, "OLD_SETTINGS_FILE", str(old_file))
    monkeypatch.chdir(tmp_path)
    return new_file, old_file


def test_add_creates_new_settings_file(settings_files, capsys):
    new_file, old_file = settings_files
    pmg.main(
        [
            "config",
            "--add",
            "PMG_VASP_PSP_DIR",
            "/path/to/psp",
            "PMG_DEFAULT_FUNCTIONAL",
            "PBE_54",
            "PMG_FLAG_TRUE",
            "True",
            "PMG_FLAG_FALSE",
            "false",
            "PMG_NONE",
            "None",
            "PMG_NULL",
            "null",
        ]
    )

    assert loadfn(new_file) == {
        "PMG_VASP_PSP_DIR": "/path/to/psp",
        "PMG_DEFAULT_FUNCTIONAL": "PBE_54",
        "PMG_FLAG_TRUE": True,
        "PMG_FLAG_FALSE": False,
        "PMG_NONE": None,
        "PMG_NULL": None,
    }
    assert not old_file.exists()
    out = capsys.readouterr().out
    assert f"New {new_file} written!" in out
    assert "backed up" not in out


@pytest.mark.parametrize(
    ("backup_args", "suffix"),
    [([], ".bak"), (["-b", ".orig"], ".orig"), (["-b", ""], None)],
)
def test_add_updates_existing_settings_file(settings_files, capsys, backup_args, suffix):
    new_file, _old_file = settings_files
    dumpfn({"PMG_MAPI_KEY": "abc", "PMG_DEFAULT_FUNCTIONAL": "PBE"}, new_file)

    pmg.main(["config", "-a", "PMG_DEFAULT_FUNCTIONAL", "PBE_64", *backup_args])

    assert loadfn(new_file) == {"PMG_MAPI_KEY": "abc", "PMG_DEFAULT_FUNCTIONAL": "PBE_64"}
    backups = list(new_file.parent.glob(".pmgrc.yaml?*"))
    if suffix is None:
        assert backups == []
        assert "backed up" not in capsys.readouterr().out
    else:
        backup = new_file.with_name(new_file.name + suffix)
        assert backups == [backup]
        assert loadfn(backup, fmt="yaml") == {"PMG_MAPI_KEY": "abc", "PMG_DEFAULT_FUNCTIONAL": "PBE"}
        assert f"backed up to {new_file}{suffix}" in capsys.readouterr().out


def test_add_falls_back_to_old_settings_file(settings_files):
    new_file, old_file = settings_files
    dumpfn({"PMG_MAPI_KEY": "abc"}, old_file)

    pmg.main(["config", "-a", "PMG_VASP_PSP_DIR", "/psp"])

    assert not new_file.exists()
    assert loadfn(old_file) == {"PMG_MAPI_KEY": "abc", "PMG_VASP_PSP_DIR": "/psp"}


def test_add_odd_number_of_tokens(settings_files):
    new_file, _old_file = settings_files
    with pytest.raises(ValueError, match="Uneven number 3 of tokens"):
        pmg.main(["config", "-a", "KEY1", "VAL1", "KEY2"])
    assert not new_file.exists()


def test_config_requires_an_action(capsys):
    with pytest.raises(SystemExit) as exc:
        pmg.main(["config"])
    assert exc.value.code == 2
    assert "one of the arguments" in capsys.readouterr().err


@pytest.mark.skipif(
    sys.platform == "win32" or shutil.which("gzip") is None or shutil.which("gunzip") is None,
    reason="setup_potcars shells out to the gzip/gunzip executables",
)
def test_setup_potcars(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    psp_dir = tmp_path / "psp_in"
    for functional_dir, element in (("potpaw_PBE", "Fe"), ("potpaw_PBE", "Osmium"), ("potpaw_LDA.54", "Li_sv")):
        (psp_dir / functional_dir / element).mkdir(parents=True)
        (psp_dir / functional_dir / element / "POTCAR").write_text(f"fake {element}\n", encoding="utf-8")
    # Compressed source POTCARs are decompressed before being renamed and re-gzipped.
    (psp_dir / "potpaw_PBE" / "P").mkdir()
    with gzip.open(psp_dir / "potpaw_PBE" / "P" / "POTCAR.gz", mode="wt", encoding="utf-8") as file:
        file.write("fake P\n")
    # Directories without a POTCAR are ignored.
    (psp_dir / "potpaw_PBE" / "empty").mkdir()
    out_dir = tmp_path / "psp_out"

    pmg.main(["config", "-p", str(psp_dir), str(out_dir)])

    expected = {
        out_dir / "POT_GGA_PAW_PBE" / "POTCAR.Fe.gz": "fake Fe\n",
        out_dir / "POT_GGA_PAW_PBE" / "POTCAR.Os.gz": "fake Osmium\n",
        out_dir / "POT_GGA_PAW_PBE" / "POTCAR.P.gz": "fake P\n",
        out_dir / "POT_LDA_PAW_54" / "POTCAR.Li_sv.gz": "fake Li_sv\n",
    }
    assert sorted(out_dir.rglob("*")) == sorted({*expected, *(path.parent for path in expected)})
    for path, content in expected.items():
        with gzip.open(path, mode="rt", encoding="utf-8") as file:
            assert file.read() == content
    out = capsys.readouterr().out
    assert "Generating pymatgen resources directory..." in out
    assert f"run 'pmg config --add PMG_VASP_PSP_DIR {out_dir}'" in out


def test_setup_potcars_existing_target_abort(tmp_path, monkeypatch):
    monkeypatch.setattr("builtins.input", lambda _prompt="": "n")
    with pytest.raises(SystemExit, match="Exiting ..."):
        pmg_config.setup_potcars([str(tmp_path), str(tmp_path)])


def test_setup_cp2k_data_existing_target_abort(tmp_path, monkeypatch):
    pytest.importorskip("yaml")
    monkeypatch.setattr("builtins.input", lambda _prompt="": "n")
    with pytest.raises(SystemExit, match="Exiting ..."):
        pmg.main(["config", "--cp2k", str(tmp_path), str(tmp_path)])


@pytest.mark.parametrize(
    ("package", "builder"),
    [("enumlib", "build_enum"), ("bader", "build_bader")],
)
def test_install_software(monkeypatch, capsys, package, builder):
    """Only the dispatch is tested; the builders clone/download and compile external code."""
    compilers = []
    monkeypatch.setattr(pmg_config.subprocess, "call", lambda cmd: compilers.append(cmd[0]))
    built = []
    monkeypatch.setattr(pmg_config, builder, lambda fortran_command: built.append(fortran_command) or True)

    pmg.main(["config", "--install", package])

    assert compilers == ["ifort"]
    assert built == ["ifort"]
    out = capsys.readouterr().out
    assert "Found ifort" in out
    assert "to your PATH" in out


def test_install_software_no_compiler(monkeypatch, capsys):
    def no_compiler(cmd):
        raise FileNotFoundError(cmd[0])

    monkeypatch.setattr(pmg_config.subprocess, "call", no_compiler)
    with pytest.raises(SystemExit, match="No fortran compiler found."):
        pmg_config.install_software("enumlib")
    assert "gfortran" in capsys.readouterr().out


BASIS_FILE_STR = """
# comment lines are skipped
 H  SZV-MOLOPT-GTH SZV-MOLOPT-GTH-q1
 1
 2 0 0 7 1
     11.478000339908  0.024916243200
      3.700758562763  0.079825490000
      1.446884268432  0.128862675300
      0.716814589696  0.379448894600
      0.247918564176  0.324552432600
      0.066918004004  0.037148121400
      0.021708243634 -0.001125195500
"""

POTENTIAL_FILE_STR = """
H GTH-PBE-q1 GTH-PBE
    1
    0.20000000    2    -4.17890044     0.72446331
    0
"""


def test_setup_cp2k_data(tmp_path, monkeypatch, capsys):
    yaml = pytest.importorskip("yaml")
    monkeypatch.chdir(tmp_path)
    data_dir = tmp_path / "cp2k_data"
    data_dir.mkdir()
    (data_dir / "BASIS_MOLOPT").write_text(BASIS_FILE_STR, encoding="utf-8")
    (data_dir / "GTH_POTENTIALS").write_text(POTENTIAL_FILE_STR, encoding="utf-8")
    # Files whose first non-comment line is not an element header cannot be chunked and are skipped.
    (data_dir / "BAD_BASIS").write_text("1 2 3\n", encoding="utf-8")
    (data_dir / "BAD_POTENTIAL").write_text("1 2 3\n", encoding="utf-8")
    out_dir = tmp_path / "cp2k_out"

    pmg.main(["config", "--cp2k", str(data_dir), str(out_dir)])

    with open(out_dir / "H", encoding="utf-8") as file:
        h_settings = yaml.safe_load(file)
    assert len(h_settings["basis_sets"]) == 1
    assert len(h_settings["potentials"]) == 1
    (basis,) = h_settings["basis_sets"].values()
    (potential,) = h_settings["potentials"].values()
    assert basis["filename"] == "BASIS_MOLOPT"
    assert basis["name"] == "SZV-MOLOPT-GTH"
    assert basis["n"] == [2]
    assert potential["filename"] == "GTH_POTENTIALS"
    assert potential["r_loc"] == pytest.approx(0.2)
    with open(out_dir / "Fe", encoding="utf-8") as file:
        assert yaml.safe_load(file) == {"basis_sets": {}, "potentials": {}}
    out = capsys.readouterr().out
    assert "Generating pymatgen resource directory for CP2K..." in out
    assert f"'pmg config --add PMG_CP2K_DATA_DIR {out_dir}'" in out
