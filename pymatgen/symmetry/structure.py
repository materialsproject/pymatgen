"""This module implements symmetry-related structure forms."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from tabulate import tabulate

from pymatgen.core.structure import PeriodicSite, Structure

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pymatgen.symmetry.analyzer import SpacegroupOperations


class SymmetrizedStructure(Structure):
    """This class represents a symmetrized structure, i.e. a structure
    where the spacegroup and symmetry operations are defined. This class is
    typically not called but instead is typically obtained by calling
    pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetrized_structure.

    Attributes:
        equivalent_indices (list[List[int]]): A list of lists of indices of the sites in the structure that are
            considered equivalent based on the symmetry operations of the space group.
    """

    def __init__(
        self,
        structure: Structure,
        spacegroup: SpacegroupOperations,
        equivalent_positions: Sequence[int],
        wyckoff_letters: Sequence[str],
    ) -> None:
        """
        Args:
            structure (Structure): Original structure
            spacegroup (SpacegroupOperations): An input SpacegroupOperations from SpacegroupAnalyzer.
            equivalent_positions (list[int]): Equivalent positions from SpacegroupAnalyzer.
            wyckoff_letters (list[str]): Wyckoff letters.
        """
        self.spacegroup = spacegroup
        uniq, inverse = np.unique(equivalent_positions, return_inverse=True)
        self.site_labels = equivalent_positions

        super().__init__(
            structure.lattice,
            [site.species for site in structure],
            structure.frac_coords,
            site_properties=structure.site_properties,
            properties=structure.properties,
        )

        equivalent_indices: list[list[int]] = [[] for _ in range(len(uniq))]
        equivalent_sites: list[list[PeriodicSite]] = [[] for _ in range(len(uniq))]
        wyckoff_symbols: list[list[str]] = [[] for _ in range(len(uniq))]
        for idx, inv_ in enumerate(inverse):
            equivalent_indices[inv_].append(idx)
            equivalent_sites[inv_].append(self.sites[idx])
            wyckoff_symbols[inv_].append(wyckoff_letters[idx])
        self.equivalent_indices = equivalent_indices
        self.equivalent_sites = equivalent_sites
        self.wyckoff_letters = wyckoff_letters
        self.wyckoff_symbols = [f"{len(symb)}{symb[0]}" for symb in wyckoff_symbols]

    def copy(self):
        """Copy of structure."""
        return SymmetrizedStructure(
            self,
            spacegroup=self.spacegroup,
            equivalent_positions=self.site_labels,
            wyckoff_letters=self.wyckoff_letters,
        )

    def find_equivalent_sites(self, site: PeriodicSite) -> list[PeriodicSite]:
        """Finds all symmetrically equivalent sites for a particular site.

        Args:
            site (PeriodicSite): A site in the structure

        Raises:
            ValueError: if site is not in the structure.

        Returns:
            ([PeriodicSite]): List of all symmetrically equivalent sites.
        """
        for sites in self.equivalent_sites:
            if site in sites:
                return sites

        raise ValueError("Site not in structure")

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        outs = [
            "SymmetrizedStructure",
            f"Full Formula ({self.composition.formula})",
            f"Reduced Formula: {self.composition.reduced_formula}",
            f"Spacegroup: {self.spacegroup.int_symbol} ({self.spacegroup.int_number})",
        ]

        def to_str(x):
            return f"{x:>10.6f}"

        outs.append(f"abc   : {' '.join(to_str(val) for val in self.lattice.abc)}")
        outs.append(f"angles: {' '.join(to_str(val) for val in self.lattice.angles)}")
        if self._charge:
            outs.append(f"Overall Charge: {self._charge:+}")
        outs.append(f"Sites ({len(self)})")
        data = []
        props = self.site_properties
        keys = sorted(props)
        for idx, sites in enumerate(self.equivalent_sites):
            site = sites[0]
            row = [str(idx), site.species_string]
            row.extend([to_str(j) for j in site.frac_coords])
            row.append(self.wyckoff_symbols[idx])
            for k in keys:
                row.append(props[k][idx])
            data.append(row)
        outs.append(tabulate(data, headers=["#", "SP", "a", "b", "c", "Wyckoff", *keys]))
        return "\n".join(outs)

    def as_dict(self):
        """MSONable dict."""
        structure = Structure.from_sites(self.sites)
        return {
            "structure": structure.as_dict(),
            "spacegroup": self.spacegroup,
            "equivalent_positions": self.site_labels,
            "wyckoff_letters": self.wyckoff_letters,
        }

    @classmethod
    def from_dict(cls, dct):
        """:param d: Dict representation

        Returns:
            SymmetrizedStructure
        """
        return cls(
            Structure.from_dict(dct["structure"]),
            spacegroup=dct["spacegroup"],
            equivalent_positions=dct["equivalent_positions"],
            wyckoff_letters=dct["wyckoff_letters"],
        )
