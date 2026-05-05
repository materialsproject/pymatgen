"""Stub computed_entries from pymatgen-core for backwards compatibility.

This module implements equivalents of the basic ComputedEntry objects, which
is the basic entity that can be used to perform many analyses. ComputedEntries
contain calculated information, typically from VASP or other electronic
structure codes. For example, ComputedEntries can be used as inputs for phase
diagram analysis.
"""

from __future__ import annotations

from pymatgen.analysis.compatibility.computed_entries import *  # noqa: F403
from pymatgen.core.entries import *  # noqa: F403
