from __future__ import annotations

import json
import runpy
import sys
from typing import TYPE_CHECKING

import pytest

from pymatgen.cli import get_environment
from tests.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path


@pytest.fixture
def chemenv_home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point the home directory (where ChemEnv keeps .chemenv/config.json) at a temp dir."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    monkeypatch.chdir(tmp_path)
    return tmp_path


def _feed_input(monkeypatch: pytest.MonkeyPatch, answers: Iterable[str]) -> list[str]:
    """Replace builtins.input with a scripted sequence of answers and record the prompts."""
    answers = iter(answers)
    prompts: list[str] = []

    def fake_input(prompt: str = "") -> str:
        prompts.append(prompt)
        return next(answers)

    monkeypatch.setattr("builtins.input", fake_input)
    return prompts


def test_main_quit_immediately(chemenv_home, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["get_environment"])
    _feed_input(monkeypatch, ["q"])

    assert get_environment.main() == 0

    out = capsys.readouterr().out
    assert "Unable to load configuration" in out
    assert "Chemical Environment package (ChemEnv)" in out
    assert "Maximum distance factor : 1.5000" in out
    assert 'Default strategy is "SimplestChemenvStrategy"' in out
    assert "Thank you for using the ChemEnv package" in out
    assert not (chemenv_home / ".chemenv").exists()


def test_dunder_main(chemenv_home, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["get_environment"])
    _feed_input(monkeypatch, ["q"])
    with pytest.raises(SystemExit) as exc:
        runpy.run_path(get_environment.__file__, run_name="__main__")
    assert exc.value.code == 0
    assert "Thank you for using the ChemEnv package" in capsys.readouterr().out


def test_main_setup_saves_and_reloads_config(chemenv_home, monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["get_environment", "--setup"])
    # setup(): "1" -> choose strategy options, "1" -> SimplestChemenvStrategy, then accept all
    # option defaults with "", "S" -> save, then "q" to quit compute_environments.
    from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy

    n_options = len(SimplestChemenvStrategy.STRATEGY_OPTIONS)
    _feed_input(monkeypatch, ["1", "1", *[""] * n_options, "S", "q"])

    assert get_environment.main() == 0

    out = capsys.readouterr().out
    assert "Configuration of the ChemEnv package" in out
    assert "Setup completed" in out
    config_file = chemenv_home / ".chemenv" / "config.json"
    assert config_file.is_file()
    with open(config_file, encoding="utf-8") as file:
        config = json.load(file)
    assert config["package_options"]["default_strategy"]["strategy"] == "SimplestChemenvStrategy"

    # A second run without --setup now loads the saved configuration.
    monkeypatch.setattr(sys, "argv", ["get_environment"])
    _feed_input(monkeypatch, ["q"])
    assert get_environment.main() == 0
    out = capsys.readouterr().out
    assert "Unable to load configuration" not in out
    assert "Thank you for using the ChemEnv package" in out


def test_main_compute_environments_from_cif(chemenv_home, monkeypatch, capsys):
    # Force the text-only branch so that no VTK window is ever requested.
    monkeypatch.setattr("pymatgen.analysis.chemenv.utils.scripts_utils.StructureVis", None)
    monkeypatch.setattr(sys, "argv", ["get_environment", "-m", "ERROR"])
    cif_path = str(TEST_FILES_DIR / "cif" / "Li2O.cif")
    # "bad" -> wrong key, cif path -> source, "n" -> do not list environments,
    # "y" -> next structure, "q" -> quit.
    # Listing ("y") is not exercised: SimplestChemenvStrategy raises TypeError on the
    # anion sites that LocalGeometryFinder skips (only_cations=True).
    prompts = _feed_input(monkeypatch, ["bad", cif_path, "n", "y", "q"])

    assert get_environment.main() == 0

    out = capsys.readouterr().out
    assert "Wrong key, try again" in out
    assert "Computing environments for Li2O" in out
    assert "Computing environments finished" in out
    assert "Thank you for using the ChemEnv package" in out
    assert len(prompts) == 5
