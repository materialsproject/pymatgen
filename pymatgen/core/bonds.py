# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This class implements definitions for various kinds of bonds. Typically used in
Molecule analysis.
"""

from __future__ import annotations

import collections
import json
import os
import warnings

from pymatgen.core.periodic_table import Element
from pymatgen.core.sites import Site
from pymatgen.util.typing import SpeciesLike


def _load_bond_length_data():
    """Loads bond length data from json file"""
    with open(os.path.join(os.path.dirname(__file__), "bond_lengths.json")) as f:
        data = collections.defaultdict(dict)
        for row in json.load(f):
            els = sorted(row["elements"])
            data[tuple(els)][row["bond_order"]] = row["length"]
        return data


bond_lengths = _load_bond_length_data()


class CovalentBond:
    """
    Defines a covalent bond between two sites.
    """

    def __init__(self, site1: Site, site2: Site):
        """
        Initializes a covalent bond between two sites.

        Args:
            site1 (Site): First site.
            site2 (Site): Second site.
        """
        self.site1 = site1
        self.site2 = site2

    @property
    def length(self) -> float:
        """
        Length of the bond.
        """
        return self.site1.distance(self.site2)

    def get_bond_order(self, tol: float = 0.2, default_bl: float | None = None) -> float:
        """
        The bond order according the distance between the two sites
        Args:
            tol (float): Relative tolerance to test.
                (1 + tol) * the longest bond distance is considered
                to be the threshold length for a bond to exist.
                (1 - tol) * the shortest bond distance is considered
                to be the shortest possible bond length
                Defaults to 0.2.
            default_bl: If a particular type of bond does not exist,
                use this bond length as a default value
                (bond order = 1). If None, a ValueError will be thrown.

        Returns:
            Float value of bond order. For example, for C-C bond in
            benzene, return 1.7.
        """
        sp1 = list(self.site1.species)[0]
        sp2 = list(self.site2.species)[0]
        dist = self.site1.distance(self.site2)
        return get_bond_order(sp1, sp2, dist, tol, default_bl)

    @staticmethod
    def is_bonded(site1, site2, tol: float = 0.2, bond_order: float | None = None, default_bl: float | None = None):
        """
        Test if two sites are bonded, up to a certain limit.

        Args:
            site1 (Site): First site
            site2 (Site): Second site
            tol (float): Relative tolerance to test. Basically, the code
                checks if the distance between the sites is less than (1 +
                tol) * typical bond distances. Defaults to 0.2, i.e.,
                20% longer.
            bond_order: Bond order to test. If None, the code simply checks
                against all possible bond data. Defaults to None.
            default_bl: If a particular type of bond does not exist, use this
                bond length. If None, a ValueError will be thrown.

        Returns:
            Boolean indicating whether two sites are bonded.
        """
        sp1 = list(site1.species)[0]
        sp2 = list(site2.species)[0]
        dist = site1.distance(site2)
        syms = tuple(sorted([sp1.symbol, sp2.symbol]))
        if syms in bond_lengths:
            all_lengths = bond_lengths[syms]
            if bond_order:
                return dist < (1 + tol) * all_lengths[bond_order]
            return any(dist < (1 + tol) * v for v in all_lengths.values())
        if default_bl:
            return dist < (1 + tol) * default_bl
        raise ValueError(f"No bond data for elements {syms[0]} - {syms[1]}")

    def __repr__(self):
        return f"Covalent bond between {self.site1} and {self.site2}"

    def __str__(self):
        return self.__repr__()


def obtain_all_bond_lengths(sp1, sp2, default_bl: float | None = None):
    """
    Obtain bond lengths for all bond orders from bond length database

    Args:
        sp1 (Species): First specie.
        sp2 (Species): Second specie.
        default_bl: If a particular type of bond does not exist, use this
            bond length as a default value (bond order = 1).
            If None, a ValueError will be thrown.

    Return:
        A dict mapping bond order to bond length in angstrom
    """
    if isinstance(sp1, Element):
        sp1 = sp1.symbol
    if isinstance(sp2, Element):
        sp2 = sp2.symbol
    syms = tuple(sorted([sp1, sp2]))
    if syms in bond_lengths:
        return bond_lengths[syms].copy()
    if default_bl is not None:
        return {1: default_bl}
    raise ValueError(f"No bond data for elements {syms[0]} - {syms[1]}")


def get_bond_order(sp1, sp2, dist: float, tol: float = 0.2, default_bl: float | None = None):
    """
    Calculate the bond order given the distance of 2 species

    Args:
        sp1 (Species): First specie.
        sp2 (Species): Second specie.
        dist: Their distance in angstrom
        tol (float): Relative tolerance to test. Basically, the code
            checks if the distance between the sites is larger than
            (1 + tol) * the longest bond distance or smaller than
            (1 - tol) * the shortest bond distance to determine if
            they are bonded or the distance is too short.
            Defaults to 0.2.
        default_bl: If a particular type of bond does not exist, use this
            bond length (bond order = 1). If None, a ValueError will be thrown.

    Returns:
        Float value of bond order. For example, for C-C bond in benzene,
        return 1.7.
    """
    all_lengths = obtain_all_bond_lengths(sp1, sp2, default_bl)
    # Transform bond lengths dict to list assuming bond data is successive
    # and add an imaginary bond 0 length
    lengths_list = [all_lengths[1] * (1 + tol)] + [all_lengths[idx + 1] for idx in range(len(all_lengths))]
    trial_bond_order = 0
    while trial_bond_order < len(lengths_list):
        if lengths_list[trial_bond_order] < dist:
            if trial_bond_order == 0:
                return trial_bond_order
            low_bl = lengths_list[trial_bond_order]
            high_bl = lengths_list[trial_bond_order - 1]
            return trial_bond_order - (dist - low_bl) / (high_bl - low_bl)
        trial_bond_order += 1
    # Distance shorter than the shortest bond length stored,
    # check if the distance is too short
    if dist < lengths_list[-1] * (1 - tol):  # too short
        warnings.warn(f"{dist:.2f} angstrom distance is too short for {sp1} and {sp2}")
    # return the highest bond order
    return trial_bond_order - 1


def get_bond_length(sp1: SpeciesLike, sp2: SpeciesLike, bond_order: float = 1) -> float:
    """
    Get the bond length between two species.

    Args:
        sp1 (Species): First specie.
        sp2 (Species): Second specie.
        bond_order: For species with different possible bond orders,
            this allows one to obtain the bond length for a particular bond
            order. For example, to get the C=C bond length instead of the
            C-C bond length, this should be set to 2. Defaults to 1.

    Returns:
        Bond length in Angstrom. If no data is available, the sum of the atomic
        radius is used.
    """
    sp1 = Element(sp1) if isinstance(sp1, str) else sp1
    sp2 = Element(sp2) if isinstance(sp2, str) else sp2
    try:
        all_lengths = obtain_all_bond_lengths(sp1, sp2)
        return all_lengths[bond_order]
    # The ValueError is raised in `obtain_all_bond_lengths` where no bond
    # data for both elements is found. The KeyError is raised in
    # `__getitem__` method of `dict` builtin class where although bond data
    # for both elements is found, the data for specified bond order does
    # not exist. In both cases, sum of atomic radius is returned.
    except (ValueError, KeyError):
        warnings.warn(
            f"No order {bond_order} bond lengths between {sp1} and {sp2} found in "
            "database. Returning sum of atomic radius."
        )
        return sp1.atomic_radius + sp2.atomic_radius  # type: ignore
