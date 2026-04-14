"""Stub for backwards compatibility.

This module provides classes to perform fitting of molecule with arbitrary
atom orders.
This module is supposed to perform exact comparisons without the atom order
correspondence prerequisite, while molecule_structure_comparator is supposed
to do rough comparisons with the atom order correspondence prerequisite.

The implementation is based on an excellent python package called `rmsd` that
you can find at https://github.com/charnley/rmsd.
"""

from __future__ import annotations

from pymatgen.core.molecule_matcher import *  # noqa: F403
