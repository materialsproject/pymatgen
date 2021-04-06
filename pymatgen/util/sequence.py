# coding: utf-8
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


class PBarSafe:
    """
    Progress bar.
    """

    def __init__(self, total):
        """
        Args:
            total (): Total value.
        """
        self.total = total
        self.done = 0
        self.report()

    def update(self, amount):
        """
        Update progress bar by amount.

        Args:
            amount (float):
        """
        self.done += amount
        self.report()

    def report(self):
        """
        Print progress.
        """
        print("{} of {} done {:.1%}".format(self.done, self.total, self.done / self.total))


try:
    # noinspection PyUnresolvedReferences
    if get_ipython().__class__.__name__ == "ZMQInteractiveShell":  # type: ignore
        from tqdm import tqdm_notebook as PBar
    else:  # likely 'TerminalInteractiveShell'
        from tqdm import tqdm as PBar
except NameError:
    try:
        from tqdm import tqdm as PBar
    except ImportError:
        PBar = PBarSafe
except ImportError:
    PBar = PBarSafe
