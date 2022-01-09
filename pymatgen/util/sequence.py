# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides utilities to chunk large sequences and display progress
bars during processing.
"""

import math


def get_chunks(sequence, size=1):
    """
    Args:
        sequence ():
        size ():

    Returns:

    """
    chunks = int(math.ceil(len(sequence) / float(size)))
    return [sequence[i * size : (i + 1) * size] for i in range(chunks)]


try:
    # noinspection PyUnresolvedReferences
    if get_ipython().__class__.__name__ == "ZMQInteractiveShell":  # type: ignore
        from tqdm import tqdm_notebook as PBar  # noqa: F401
    else:  # likely 'TerminalInteractiveShell'
        from tqdm import tqdm as PBar  # noqa: F401
except (NameError, ImportError):
    from tqdm import tqdm as PBar  # noqa: F401
