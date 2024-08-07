"""This module provides a wrapper for numba such that no functionality
is lost if numba is not available. Numba is a just-in-time compiler
that can significantly accelerate the evaluation of certain functions
if installed.
"""

from __future__ import annotations

try:
    from numba import jit, njit
except ImportError:

    def njit(func):
        """Replacement for numba.njit when numba is not installed that does nothing."""
        return func

    def jit(func):
        """Replacement for numba.jit when numba is not installed that does nothing."""
        return func
