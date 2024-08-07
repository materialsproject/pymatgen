from __future__ import annotations

from typing import TYPE_CHECKING

from pymatgen.core import _load_pmg_settings

if TYPE_CHECKING:
    from pathlib import Path

    from pytest import MonkeyPatch


__author__ = "Janosh Riebesell"
__date__ = "2022-10-21"
__email__ = "janosh@lbl.gov"


def test_load_settings(tmp_path: Path, monkeypatch: MonkeyPatch) -> None:
    """Test .pmgrc.yaml file is loaded correctly and env vars take precedence."""
    monkeypatch.setattr("os.environ", {})  # reset outer env vars

    settings_file = tmp_path / ".pmgrc.yaml"
    monkeypatch.setattr("pymatgen.core.SETTINGS_FILE", settings_file)

    # should return empty dict if file doesn't exist
    assert _load_pmg_settings() == {}

    # should return empty dict if file is empty
    settings_file.write_text("")
    assert _load_pmg_settings() == {}

    settings_file.write_text("PMG_VASP_PSP_DIR: /path/to/psp")
    assert _load_pmg_settings() == {"PMG_VASP_PSP_DIR": "/path/to/psp"}

    settings_file.write_text("PMG_MAPI_KEY: FOOBAR")
    assert _load_pmg_settings() == {"PMG_MAPI_KEY": "FOOBAR"}

    # env vars should override .pmgrc.yaml
    with monkeypatch.context() as ctx:
        ctx.setenv("PMG_MAPI_KEY", "BAZ")
        assert _load_pmg_settings() == {"PMG_MAPI_KEY": "BAZ"}

    with monkeypatch.context() as ctx:
        ctx.setenv("HOME", "/home/fakeuser")
        ctx.setenv("PMG_VASP_PSP_DIR", "$HOME/psp")
        assert _load_pmg_settings()["PMG_VASP_PSP_DIR"] == "/home/fakeuser/psp"

    # should return empty dict if file is invalid
    settings_file.write_text("---")
    assert _load_pmg_settings() == {}


def test_env_var_pmg_config_file(tmp_path: Path, monkeypatch: MonkeyPatch) -> None:
    custom_config_file = tmp_path / "custom_config.yaml"
    custom_config_file.write_text("PMG_CUSTOM_SETTING: custom_value")

    with monkeypatch.context() as ctx:
        ctx.setenv("PMG_CONFIG_FILE", str(custom_config_file))
        settings = _load_pmg_settings()
        assert "PMG_CUSTOM_SETTING" in settings
        assert settings["PMG_CUSTOM_SETTING"] == "custom_value"

    # Test that PMG_CONFIG_FILE takes precedence over the default location
    settings_file = tmp_path / ".pmgrc.yaml"
    monkeypatch.setattr("pymatgen.core.SETTINGS_FILE", settings_file)
    settings_file.write_text("PMG_DEFAULT_SETTING: default_value")
    custom_config_file.write_text("PMG_CUSTOM_SETTING: custom_value")

    with monkeypatch.context() as ctx:
        ctx.setenv("PMG_CONFIG_FILE", str(custom_config_file))
        settings = _load_pmg_settings()
        assert "PMG_CUSTOM_SETTING" in settings
        assert "PMG_DEFAULT_SETTING" not in settings
        assert settings["PMG_CUSTOM_SETTING"] == "custom_value"

    # Test that env vars still take precedence over the values specified in PMG_CONFIG_FILE
    with monkeypatch.context() as ctx:
        ctx.setenv("PMG_CONFIG_FILE", str(custom_config_file))
        ctx.setenv("PMG_CUSTOM_SETTING", "env_value")
        settings = _load_pmg_settings()
        assert settings["PMG_CUSTOM_SETTING"] == "env_value"
