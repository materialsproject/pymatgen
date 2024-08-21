"""This module provides various methods to analyze order/disorder in materials."""

from __future__ import annotations

import itertools
from collections import defaultdict
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pymatgen.core import Structure


def get_warren_cowley_parameters(structure: Structure, r: float, dr: float) -> dict[tuple, float]:
    """
    Warren-Crowley parameters.

    Args:
        structure: Pymatgen Structure.
        r: Radius
        dr: Shell width

    Returns:
        Warren-Crowley parameters in the form of a dict, e.g. {(Element Mo, Element W): -1.0, ...}
    """
    comp = structure.composition

    n_ij: defaultdict = defaultdict(int)
    n_neighbors: defaultdict = defaultdict(int)
    for site in structure:
        for nn in structure.get_neighbors_in_shell(site.coords, r, dr):
            n_ij[site.specie, nn.specie] += 1
            n_neighbors[site.specie] += 1

    alpha_ij = {}
    for sp1, sp2 in itertools.product(comp, comp):
        pij = n_ij.get((sp1, sp2), 0) / n_neighbors[sp1]
        conc2 = comp.get_atomic_fraction(sp2)
        alpha_ij[sp1, sp2] = (pij - conc2) / ((1 if sp1 == sp2 else 0) - conc2)

    return alpha_ij
