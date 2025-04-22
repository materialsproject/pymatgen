from __future__ import annotations

import typing

import pytest
from monty.tempfile import ScratchDir


@pytest.fixture(autouse=True)
def setup_teardown() -> typing.Generator:
    """Use ScratchDir for all tests in this session."""
    with ScratchDir(".", copy_from_current_on_enter=False, copy_to_current_on_exit=False):
        yield
