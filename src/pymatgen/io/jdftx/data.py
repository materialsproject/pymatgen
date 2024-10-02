"""Module containing reference data for JDFTx.

This module contains reference data for JDFTx.
"""

from __future__ import annotations

import numpy as np

from pymatgen.core.periodic_table import Element


def get_atom_valence_electrons(el: str) -> int:
    """Return number of electrons in valence shell(s).

    Return number of valence electrons for an element. This is mostly a copy-paste
    of the valence property for the Element class in pymatgen, but supersedes
    the error raised by ambiguous valence shells (multiple partially filled shells).
    """
    pmg_el = Element(el)
    if pmg_el.group == 18:
        return 0  # The number of valence of noble gas is 0

    L_symbols = "SPDFGHIKLMNOQRTUVWXYZ"
    valence: list[tuple[int, int]] = []
    full_electron_config = pmg_el.full_electronic_structure
    last_orbital = full_electron_config[-1]
    for n, l_symbol, ne in full_electron_config:
        idx = L_symbols.lower().index(l_symbol)
        if ne < (2 * idx + 1) * 2 or (
            (n, l_symbol, ne) == last_orbital and ne == (2 * idx + 1) * 2 and len(valence) == 0
        ):  # check for full last shell (e.g. column 2)
            valence.append((idx, ne))
    # if len(valence) > 1:
    #     raise ValueError(f"{pmg_el} has ambiguous valence")
    return np.sum(np.array([v[1] for v in valence]))
