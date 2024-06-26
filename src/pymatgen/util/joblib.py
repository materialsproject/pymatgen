import contextlib
from collections.abc import Iterator
from typing import Any
import os

import joblib
from tqdm import tqdm


@contextlib.contextmanager
def tqdm_joblib(tqdm_object: tqdm) -> Iterator[None]:
    """Context manager to patch joblib to report into tqdm progress bar given
    as argument.
    """

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args: tuple, **kwargs: dict[str, Any]) -> None:
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
    original_warnings = os.environ.get('PYTHONWARNINGS')
    os.environ['PYTHONWARNINGS'] = warnings
    try:
        yield
    finally:
        if original_warnings is None:
            del os.environ['PYTHONWARNINGS']
        else:
            os.environ['PYTHONWARNINGS'] = original_warnings
