from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from pytest import MonkeyPatch

from pymatgen.core import _load_pmg_settings

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

    # should return empty dict if file is invalid
    settings_file.write_text("---")
    assert _load_pmg_settings() == {}
