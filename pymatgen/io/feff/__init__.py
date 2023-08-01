"""
This package provides the modules to perform FEFF IO.

FEFF: http://feffproject.org/feffproject-feff.html
"""

from __future__ import annotations

from .inputs import (
    VALID_FEFF_TAGS,
    Atoms,
    FeffParseError,
    Header,
    Paths,
    Potential,
    Tags,
    get_absorbing_atom_symbol_index,
    get_atom_map,
)
from .outputs import Eels, LDos, Xmu
