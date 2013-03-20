"""
This module provides utilities for basic math operations.
"""
from __future__ import division, print_function

import itertools

def iterator_from_slice(s):
    """
    Constructs an iterator given a slice object s.

    :note: The function returns an infinite iterator if s.stop is None
    """
    start = s.start if s.start is not None else 0
    step  = s.step  if s.step is not None else 1

    if s.stop is None: 
        # Infinite iterator.
        return itertools.count(start=start, step=step)
    else:
        # xrange-like iterator.
        return iter(range(start, s.stop, step))
