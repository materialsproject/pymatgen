"""This module implements functions to perform various useful operations on
entries, such as grouping entries by structure.
"""

from __future__ import annotations

import warnings

from pymatgen.analysis.compatibility.entry_tools import *  # noqa: F403

warnings.warn(
    "All Entry tools have been moved to the pymatgen.analysis.compatibility.entry_tools package. "
    "This stub will be removed v2027.1.",
    DeprecationWarning,
    stacklevel=2,
)
