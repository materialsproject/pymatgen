"""This module implements equivalents of the basic ComputedEntry objects, which
is the basic entity that can be used to perform many analyses. ComputedEntries
contain calculated information, typically from VASP or other electronic
structure codes. For example, ComputedEntries can be used as inputs for phase
diagram analysis.
"""

from __future__ import annotations

import warnings

from pymatgen.core.entries import *  # noqa: F403

warnings.warn(
    "Entry, ComputedEntry and ComputedStructureEntry have been moved to pymatgen.core.entries module. "
    "This stub will be removed v2027.1.",
    DeprecationWarning,
    stacklevel=2,
)
