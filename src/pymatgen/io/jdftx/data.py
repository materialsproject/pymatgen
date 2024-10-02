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
    return np.sum(v[1] for v in valence)


# atom_valence_electrons = {
#     "H": 1,
#     "He": 2,
#     "Li": 1,
#     "Be": 2,
#     "B": 3,
#     "C": 4,
#     "N": 5,
#     "O": 6,
#     "F": 7,
#     "Ne": 8,
#     "Na": 1,
#     "Mg": 2,
#     "Al": 3,
#     "Si": 4,
#     "P": 5,
#     "S": 6,
#     "Cl": 7,
#     "Ar": 8,
#     "K": 1,
#     "Ca": 2,
#     "Sc": 3,
#     "Ti": 4,
#     "V": 5,
#     "Cr": 6,
#     "Mn": 7,
#     "Fe": 8,
#     "Co": 9,
#     "Ni": 10,
#     "Cu": 11,
#     "Zn": 12,
#     "Ga": 3,
#     "Ge": 4,
#     "As": 5,
#     "Se": 6,
#     "Br": 7,
#     "Kr": 8,
#     "Rb": 1,
#     "Sr": 2,
#     "Y": 3,
#     "Zr": 4,
#     "Nb": 5,
#     "Mo": 6,
#     "Tc": 7,
#     "Ru": 8,
#     "Rh": 9,
#     "Pd": 10,
#     "Ag": 11,
#     "Cd": 12,
#     "In": 3,
#     "Sn": 4,
#     "Sb": 5,
#     "Te": 6,
#     "I": 7,
#     "Xe": 8,
#     "Cs": 1,
#     "Ba": 2,
#     "Hf": 4,
#     "Ta": 5,
#     "W": 6,
#     "Re": 7,
#     "Os": 8,
#     "Ir": 9,
#     "Pt": 10,
#     "Au": 11,
#     "Hg": 12,
#     "Tl": 3,
#     "Pb": 4,
#     "Bi": 5,
#     "La": 3,
#     "Ce": 4,
#     "Gd": 10,
# }
