from __future__ import annotations

import os
import tempfile
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from collections.abc import Generator


@pytest.fixture(autouse=True)
def setup_teardown() -> Generator:
    """Use tempdir for all tests."""
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        yield tmpdir
        os.chdir(cwd)


# Value at session start, before any test can modify it.
_SPGLIB_ERROR_HANDLING_ENV = "SPGLIB_OLD_ERROR_HANDLING"
_SPGLIB_ERROR_HANDLING_INITIAL = os.environ.get(_SPGLIB_ERROR_HANDLING_ENV)


@pytest.fixture(autouse=True)
def restore_spglib_error_handling() -> Generator:
    """Undo leaks of SPGLIB_OLD_ERROR_HANDLING between tests.

    ase.utils.spglib_new_errorhandling sets it to "false" and does not unset it when it
    was previously unset, which switches spglib to raising mode for all later tests.

    spglib.error.OLD_ERROR_HANDLING is deliberately not reset: phonopy and symfc set it to
    False once at import and rely on raising mode afterwards. Tests that depend on spglib's
    legacy (return None) behaviour should pin it themselves with monkeypatch.
    """
    yield
    if _SPGLIB_ERROR_HANDLING_INITIAL is None:
        os.environ.pop(_SPGLIB_ERROR_HANDLING_ENV, None)
    else:
        os.environ[_SPGLIB_ERROR_HANDLING_ENV] = _SPGLIB_ERROR_HANDLING_INITIAL
