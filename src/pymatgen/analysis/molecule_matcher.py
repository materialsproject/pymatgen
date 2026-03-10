"""
This module provides classes to perform fitting of molecule with arbitrary
atom orders.
This module is supposed to perform exact comparisons without the atom order
correspondence prerequisite, while molecule_structure_comparator is supposed
to do rough comparisons with the atom order correspondence prerequisite.

The implementation is based on an excellent python package called `rmsd` that
you can find at https://github.com/charnley/rmsd.
"""

from __future__ import annotations

import warnings

warnings.warn(
    "This module has been moved to the pymatgen.core.molecule_matcher module. This stub will be removed v2027.1. ",
    DeprecationWarning,
    stacklevel=2,
)


from pymatgen.core.molecule_matcher import *  # noqa: E402, F403
