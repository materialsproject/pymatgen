"""This module defines Entry classes for containing experimental data."""

from __future__ import annotations

import warnings

from pymatgen.analysis.compatibility.exp_entries import *  # noqa: F403

warnings.warn(
    "All exp_entries have been moved to the pymatgen.analysis.compatibility.exp_entries module. "
    "This stub will be removed v2027.1.",
    DeprecationWarning,
    stacklevel=2,
)
