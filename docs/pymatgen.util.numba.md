---
layout: default
title: pymatgen.util.numba.md
nav_exclude: true
---

# pymatgen.util.numba module

This module provides a wrapper for numba such that no functionality
is lost if numba is not available. Numba is a just-in-time compiler
that can significantly accelerate the evaluation of certain functions
if installed.


### pymatgen.util.numba.jit(func)
Replacement for numba.jit when numba is not installed that does nothing.


### pymatgen.util.numba.njit(func)
Replacement for numba.njit when numba is not installed that does nothing.