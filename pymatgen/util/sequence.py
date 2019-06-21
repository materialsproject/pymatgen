# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import math

"""
This module provides utilities to chunk large sequences and display progress
bars during processing.
"""


def get_chunks(sequence, size=1):
    chunks = int(math.ceil(len(sequence) / float(size)))
    return [sequence[i * size:(i + 1) * size]
            for i in range(chunks)]

class PBarSafe:
    def __init__(self, total):
        self.total = total
        self.done = 0
        self.report()

    def update(self, amount):
        self.done += amount
        self.report()

    def report(self):
        print("{} of {} done {:.1%}".format(
            self.done, self.total, self.done / self.total))

try:
    # noinspection PyUnresolvedReferences
    if get_ipython().__class__.__name__ == 'ZMQInteractiveShell':
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
