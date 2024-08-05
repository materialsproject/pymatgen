"""This module implements definitions for various kinds of bonds. Typically used in
Molecule analysis.
"""

from __future__ import annotations

import json
import os
import warnings
from collections import defaultdict
from typing import TYPE_CHECKING

from pymatgen.core import Element

if TYPE_CHECKING:
    from pymatgen.core.sites import Site
    from pymatgen.util.typing import SpeciesLike


def _load_bond_length_data() -> dict[tuple[str, ...], dict[float, float]]:
    """Load bond length data from bond_lengths.json file."""
    with open(
        os.path.join(os.path.dirname(__file__), "bond_lengths.json"),
        encoding="utf-8",
    ) as file:
        data: dict[tuple, dict] = defaultdict(dict)
        for row in json.load(file):
            els = sorted(row["elements"])
            data[tuple(els)][row["bond_order"]] = row["length"]
        return data


bond_lengths = _load_bond_length_data()


class CovalentBond:
    """A covalent bond between two sites."""

    def __init__(self, site1: Site, site2: Site) -> None:
        """Initialize a covalent bond between two sites.

        Args:
            site1 (Site): First site.
            site2 (Site): Second site.
        """
        self.site1 = site1
        self.site2 = site2

    def __repr__(self) -> str:
        return f"Covalent bond between {self.site1} and {self.site2}"

    @property
    def length(self) -> float:
        """Length of the bond."""
        return self.site1.distance(self.site2)

    def get_bond_order(
        self,
        tol: float = 0.2,
        default_bl: float | None = None,
    ) -> float:
        """The bond order according the distance between the two sites.

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
            float: value of bond order. E.g. 1.7 for C-C bond in benzene.
        """
        sp1 = next(iter(self.site1.species))
        sp2 = next(iter(self.site2.species))
        dist = self.site1.distance(self.site2)
        return get_bond_order(sp1, sp2, dist, tol, default_bl)

    @staticmethod
    def is_bonded(
        site1: Site,
        site2: Site,
        tol: float = 0.2,
        bond_order: float | None = None,
        default_bl: float | None = None,
    ) -> bool:
        """Check if two sites are bonded, up to a certain limit.

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
            bool: True if two sites are bonded.
        """
        sp1 = next(iter(site1.species))
        sp2 = next(iter(site2.species))
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


def obtain_all_bond_lengths(
    sp1: SpeciesLike,
    sp2: SpeciesLike,
    default_bl: float | None = None,
) -> dict[float, float]:
    """Obtain bond lengths for all bond orders from bond length database.

    Args:
        sp1 (Species): First specie.
        sp2 (Species): Second specie.
        default_bl: If a particular type of bond does not exist, use this
            bond length as a default value (bond order = 1).
            If None, a ValueError will be thrown.

    Returns:
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
        return {1.0: default_bl}
    raise ValueError(f"No bond data for elements {syms[0]} - {syms[1]}")


def get_bond_order(
    sp1: SpeciesLike,
    sp2: SpeciesLike,
    dist: float,
    tol: float = 0.2,
    default_bl: float | None = None,
) -> float:
    """Calculate the bond order given the distance of 2 species.

    Args:
        sp1 (Species): First specie.
        sp2 (Species): Second specie.
        dist (float): Distance in angstrom
        tol (float): Relative tolerance to test. Basically, the code
            checks if the distance between the sites is larger than
            (1 + tol) * the longest bond distance or smaller than
            (1 - tol) * the shortest bond distance to determine if
            they are bonded or the distance is too short.
            Defaults to 0.2.
        default_bl: If a particular type of bond does not exist, use this
            bond length (bond order = 1). If None, a ValueError will be thrown.

    Returns:
        float: Bond order. For example, 1.7 for C-C bond in benzene.
    """
    all_lens = obtain_all_bond_lengths(sp1, sp2, default_bl)
    # Transform bond lengths dict to list assuming bond data is successive
    # and add an imaginary bond 0 length
    lens = [all_lens[1] * (1 + tol)] + [all_lens[idx + 1] for idx in range(len(all_lens))]
    trial_bond_order = 0
    while trial_bond_order < len(lens):
        if lens[trial_bond_order] < dist:
            if trial_bond_order == 0:
                return trial_bond_order
            low_bl = lens[trial_bond_order]
            high_bl = lens[trial_bond_order - 1]
            return trial_bond_order - (dist - low_bl) / (high_bl - low_bl)
        trial_bond_order += 1
    # Distance shorter than the shortest bond length stored,
    # check if the distance is too short
    if dist < lens[-1] * (1 - tol):  # too short
        warnings.warn(f"{dist:.2f} angstrom distance is too short for {sp1} and {sp2}")
    # return the highest bond order
    return trial_bond_order - 1


def get_bond_length(
    sp1: SpeciesLike,
    sp2: SpeciesLike,
    bond_order: float = 1,
) -> float:
    """Get the bond length between two species.

    Args:
        sp1 (Species): First specie.
        sp2 (Species): Second specie.
        bond_order: For species with different possible bond orders,
            this allows one to obtain the bond length for a particular bond
            order. For example, to get the C=C bond length instead of the
            C-C bond length, this should be set to 2. Defaults to 1.

    Returns:
        float: Bond length in Angstrom. If no data is available,
            the sum of the atomic radius is used.
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
        return sp1.atomic_radius + sp2.atomic_radius  # type: ignore[operator]
