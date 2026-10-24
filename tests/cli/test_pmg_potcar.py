from __future__ import annotations

from argparse import Namespace

import pytest

from pymatgen.cli import pmg
from pymatgen.cli.pmg_potcar import generate_potcar
from pymatgen.core import SETTINGS
from pymatgen.io.vasp.inputs import Potcar, PotcarSingle
from tests.testing import TEST_FILES_DIR

FAKE_POTCAR_DIR = TEST_FILES_DIR / "io" / "vasp" / "inputs" / "fake_potcars"


@pytest.fixture(autouse=True)
def fake_psp_dir(tmp_path, monkeypatch):
    monkeypatch.setitem(SETTINGS, "PMG_VASP_PSP_DIR", str(FAKE_POTCAR_DIR))
    monkeypatch.delitem(SETTINGS, "PMG_DEFAULT_FUNCTIONAL", raising=False)
    monkeypatch.chdir(tmp_path)


def _assert_potcar_from(path, functional_dir: str, symbols: list[str]) -> None:
    """Check that the POTCAR at path concatenates the fake POTCARs from functional_dir."""
    potcar = Potcar.from_file(path)
    assert potcar.symbols == symbols
    for single, symbol in zip(potcar, symbols, strict=True):
        source = PotcarSingle.from_file(FAKE_POTCAR_DIR / functional_dir / f"POTCAR.{symbol}.gz")
        assert single.data == source.data


@pytest.mark.parametrize(
    ("functional_args", "functional_dir"),
    [([], "POT_GGA_PAW_PBE"), (["-f", "LDA"], "POT_LDA_PAW")],
)
def test_symbols(tmp_path, functional_args, functional_dir):
    pmg.main(["potcar", "-s", "Fe", "P", *functional_args])
    _assert_potcar_from(tmp_path / "POTCAR", functional_dir, ["Fe", "P"])


def test_symbols_default_functional_from_settings(tmp_path, monkeypatch):
    monkeypatch.setitem(SETTINGS, "PMG_DEFAULT_FUNCTIONAL", "LDA_64")
    pmg.main(["potcar", "-s", "Si"])
    _assert_potcar_from(tmp_path / "POTCAR", "POT_LDA_PAW_64", ["Si"])


@pytest.mark.parametrize(
    ("args", "msg"),
    [
        (["-s", "Fe", "-f", "NOT_A_FUNCTIONAL"], "Invalid functional 'NOT_A_FUNCTIONAL'"),
        (["-s", "Xx"], "An error has occurred"),
    ],
)
def test_symbols_error_is_reported(tmp_path, capsys, args, msg):
    pmg.main(["potcar", *args])
    assert msg in capsys.readouterr().out
    assert not (tmp_path / "POTCAR").exists()


def test_recursive(tmp_path):
    for name, symbols in (("a", ["Fe", "P", "O"]), ("b/c", ["Si"])):
        (tmp_path / name).mkdir(parents=True)
        (tmp_path / name / "POTCAR.spec").write_text("\n".join(symbols) + "\n\n", encoding="utf-8")
    (tmp_path / "b" / "other_file").write_text("ignored", encoding="utf-8")

    pmg.main(["potcar", "-r", str(tmp_path)])

    _assert_potcar_from(tmp_path / "a" / "POTCAR", "POT_GGA_PAW_PBE", ["Fe", "P", "O"])
    _assert_potcar_from(tmp_path / "b" / "c" / "POTCAR", "POT_GGA_PAW_PBE", ["Si"])
    assert not (tmp_path / "b" / "POTCAR").exists()


def test_no_valid_option(capsys):

    generate_potcar(Namespace(functional=None, recursive=None, symbols=None))
    assert "No valid options selected." in capsys.readouterr().out


def test_requires_symbols_or_recursive(capsys):
    with pytest.raises(SystemExit) as exc:
        pmg.main(["potcar"])
    assert exc.value.code == 2
    assert "one of the arguments -s/--symbols -r/--recursive is required" in capsys.readouterr().err
