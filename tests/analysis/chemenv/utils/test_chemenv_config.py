from __future__ import annotations

import copy
import os

import pytest

import pymatgen.analysis.chemenv.utils.chemenv_config as chemenv_config_module
from pymatgen.analysis.chemenv.utils.chemenv_config import ChemEnvConfig
from pymatgen.core import SETTINGS
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

__author__ = "waroquiers"

config_file_dir = f"{TEST_FILES_DIR}/analysis/chemenv/config"


class TestChemenvConfig(MatSciTest):
    def test_chemenv_config(self):
        config = ChemEnvConfig()

        assert config.has_materials_project_access == bool(SETTINGS.get("PMG_MAPI_KEY"))

        package_options = {
            **ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS,
            "default_max_distance_factor": 1.8,
        }

        config = ChemEnvConfig(package_options=package_options)

        assert (
            config.package_options_description() == "Package options :\n"
            " - Maximum distance factor : 1.8000\n"
            ' - Default strategy is "SimplestChemenvStrategy" :\n'
            "Simplest ChemenvStrategy using fixed angle and distance parameters \n"
            "for the definition of neighbors in the Voronoi approach. \n"
            "The coordination environment is then given as the one with the \n"
            "lowest continuous symmetry measure.\n"
            "   with options :\n"
            "     - distance_cutoff : 1.4\n"
            "     - angle_cutoff : 0.3\n"
            "     - additional_condition : 1\n"
            "     - continuous_symmetry_measure_cutoff : 10\n"
        )

        config.save(root_dir=self.tmp_path)

        config = config.auto_load(root_dir=self.tmp_path)

        assert config.package_options == package_options


@pytest.fixture
def fake_home(tmp_path, monkeypatch):
    """Redirect the home directory to a temporary directory."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    return tmp_path


def _feed_input(monkeypatch, answers):
    answers_iter = iter(answers)
    monkeypatch.setattr("builtins.input", lambda *_args: next(answers_iter))


@pytest.mark.parametrize(("mapi_key", "expected"), [("fake-key", True), (None, False)])
def test_materials_project_access(monkeypatch, mapi_key, expected):
    monkeypatch.setitem(chemenv_config_module.SETTINGS, "PMG_MAPI_KEY", mapi_key)
    config = ChemEnvConfig()
    assert config.has_materials_project_access is expected
    assert config.materials_project_configuration == mapi_key


def test_save_and_auto_load_default_root_dir(fake_home, monkeypatch):
    package_options = {**ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS, "default_max_distance_factor": 2.0}
    config_file = ChemEnvConfig(package_options=package_options).save()
    assert config_file == f"{fake_home}/.chemenv/config.json"
    assert os.path.isfile(config_file)

    assert ChemEnvConfig.auto_load().package_options == package_options

    # Existing file: declining the overwrite keeps the previous content
    _feed_input(monkeypatch, ["n"])
    ChemEnvConfig(package_options={**package_options, "default_max_distance_factor": 3.0}).save()
    assert ChemEnvConfig.auto_load().package_options["default_max_distance_factor"] == 2.0

    # Confirming the overwrite replaces it
    _feed_input(monkeypatch, ["Y"])
    ChemEnvConfig(package_options={**package_options, "default_max_distance_factor": 3.0}).save()
    assert ChemEnvConfig.auto_load().package_options["default_max_distance_factor"] == 3.0


def test_auto_load_missing_file_falls_back_to_default(tmp_path, capsys):
    config = ChemEnvConfig.auto_load(root_dir=tmp_path / "does_not_exist")
    assert config.package_options == ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS
    assert "loading default configuration" in capsys.readouterr().out


def test_setup_quit_and_wrong_key(monkeypatch, capsys):
    _feed_input(monkeypatch, ["x", "q"])
    ChemEnvConfig().setup()
    out = capsys.readouterr().out
    assert "wrong key, try again" in out
    assert out.count("Configuration of the ChemEnv package") == 2


def test_setup_package_options_and_save(fake_home, monkeypatch):
    default_options = copy.deepcopy(ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS)
    # Menu "1" -> strategy 1 (SimplestChemenvStrategy); invalid then valid distance cutoff;
    # defaults for the remaining options; then "S" to save.
    _feed_input(monkeypatch, ["1", "1", "abc", "1.6", "", "", "", "S"])
    config = ChemEnvConfig()
    config.setup()

    strategy = config.package_options["default_strategy"]
    assert strategy["strategy"] == "SimplestChemenvStrategy"
    assert strategy["strategy_options"] == {
        "distance_cutoff": 1.6,
        "angle_cutoff": 0.3,
        "additional_condition": 1,
        "continuous_symmetry_measure_cutoff": 10,
    }
    # The interactive setup must not mutate the class-level defaults
    assert default_options == ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS
    loaded = ChemEnvConfig.auto_load()
    assert loaded.package_options["default_strategy"]["strategy_options"]["distance_cutoff"] == 1.6
