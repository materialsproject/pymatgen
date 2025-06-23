from __future__ import annotations

import os
import tempfile
import typing

import pytest


@pytest.fixture(autouse=True)
def setup_teardown() -> typing.Generator:
    """Use tempdir for all tests."""
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        yield tmpdir
        os.chdir(cwd)
