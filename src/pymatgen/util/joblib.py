"""This module provides utility functions for getting progress bar with joblib."""

from __future__ import annotations

import contextlib
import os
from typing import TYPE_CHECKING, Any

import joblib

if TYPE_CHECKING:
    from collections.abc import Iterator

    from tqdm import tqdm


@contextlib.contextmanager
def tqdm_joblib(tqdm_object: tqdm) -> Iterator[None]:
    """Context manager to patch joblib to report into tqdm progress bar given
    as argument.
    """

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args: tuple, **kwargs: dict[str, Any]) -> None:
            """This will be called after each batch, to update the progress bar."""
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


@contextlib.contextmanager
def set_python_warnings(warnings):
    """Context manager to set the PYTHONWARNINGS environment variable to the
    given value. This is useful for preventing spam when using parallel processing.
    """
    original_warnings = os.environ.get("PYTHONWARNINGS")
    os.environ["PYTHONWARNINGS"] = warnings
    try:
        yield
    finally:
        if original_warnings is None:
            del os.environ["PYTHONWARNINGS"]
        else:
            os.environ["PYTHONWARNINGS"] = original_warnings
