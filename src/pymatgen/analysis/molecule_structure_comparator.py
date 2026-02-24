"""
This module provides classes to comparison the structures of the two
molecule. As long as the two molecule have the same bond connection tables,
the molecules are deemed to be same. The atom in the two molecule must be
paired accordingly.
This module is supposed to perform rough comparisons with the atom order
correspondence prerequisite, while molecule_matcher is supposed to do exact
comparisons without the atom order correspondence prerequisite.
"""

from __future__ import annotations

import warnings

warnings.warn(
    "This module has been moved to the pymatgen.core.molecule_structure_comparator module. This stub will be removed"
    " v2027.1. ",
    DeprecationWarning,
    stacklevel=2,
)

from pymatgen.core.molecule_structure_comparator import *  # noqa: E402, F403
