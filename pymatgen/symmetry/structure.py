# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements symmetry-related structure forms.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
from tabulate import tabulate

from pymatgen.core.structure import PeriodicSite, Structure


class SymmetrizedStructure(Structure):
    """
    This class represents a symmetrized structure, i.e. a structure
    where the spacegroup and symmetry operations are defined. This class is
    typically not called but instead is typically obtained by calling
    pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetrized_structure.

    .. attribute: equivalent_indices

        indices of structure grouped by equivalency
    """

    def __init__(
        self,
        structure: Structure,
        spacegroup,
        equivalent_positions: Sequence[int],
        wyckoff_letters: Sequence[str],
    ):
        """
        Args:
            structure (Structure): Original structure
            spacegroup (SpacegroupOperations): An input SpacegroupOperations from
                SpacegroupAnalyzer.
            equivalent_positions: Equivalent positions from SpacegroupAnalyzer.
            wyckoff_letters: Wyckoff letters
        """
        self.spacegroup = spacegroup
        u, inv = np.unique(equivalent_positions, return_inverse=True)
        self.site_labels = equivalent_positions

        super().__init__(
            structure.lattice,
            [site.species for site in structure],
            structure.frac_coords,
            site_properties=structure.site_properties,
        )

        equivalent_indices = [[] for _ in range(len(u))]  # type: ignore
        equivalent_sites = [[] for _ in range(len(u))]  # type: ignore
        wyckoff_symbols = [[] for _ in range(len(u))]  # type: ignore
        for i, inv_ in enumerate(inv):
            equivalent_indices[inv_].append(i)
            equivalent_sites[inv_].append(self.sites[i])
            wyckoff_symbols[inv_].append(wyckoff_letters[i])
        self.equivalent_indices: list[int] = equivalent_indices  # type: ignore
        self.equivalent_sites: list[PeriodicSite] = equivalent_sites  # type: ignore
        self.wyckoff_letters = wyckoff_letters
        self.wyckoff_symbols = [f"{len(w)}{w[0]}" for w in wyckoff_symbols]

    def copy(self):
        """
        :return: Copy of structure.
        """
        return SymmetrizedStructure(
            self,
            spacegroup=self.spacegroup,
            equivalent_positions=self.site_labels,
            wyckoff_letters=self.wyckoff_letters,
        )

    def find_equivalent_sites(self, site) -> list[PeriodicSite]:
        """
        Finds all symmetrically equivalent sites for a particular site

        Args:
            site (PeriodicSite): A site in the structure

        Returns:
            ([PeriodicSite]): List of all symmetrically equivalent sites.
        """
        for sites in self.equivalent_sites:
            if site in sites:
                return sites

        raise ValueError("Site not in structure")

    def __repr__(self):
        return str(self)

    def __str__(self):
        outs = [
            "SymmetrizedStructure",
            f"Full Formula ({self.composition.formula})",
            f"Reduced Formula: {self.composition.reduced_formula}",
            f"Spacegroup: {self.spacegroup.int_symbol} ({self.spacegroup.int_number})",
        ]

        def to_s(x):
            return f"{x:0.6f}"

        outs.append("abc   : " + " ".join(to_s(i).rjust(10) for i in self.lattice.abc))
        outs.append("angles: " + " ".join(to_s(i).rjust(10) for i in self.lattice.angles))
        if self._charge:
            if self._charge >= 0:
                outs.append(f"Overall Charge: +{self._charge}")
            else:
                outs.append(f"Overall Charge: -{self._charge}")
        outs.append(f"Sites ({len(self)})")
        data = []
        props = self.site_properties
        keys = sorted(props)
        for i, sites in enumerate(self.equivalent_sites):
            site = sites[0]
            row = [str(i), site.species_string]
            row.extend([to_s(j) for j in site.frac_coords])
            row.append(self.wyckoff_symbols[i])
            for k in keys:
                row.append(props[k][i])
            data.append(row)
        outs.append(
            tabulate(
                data,
                headers=["#", "SP", "a", "b", "c", "Wyckoff"] + keys,
            )
        )
        return "\n".join(outs)

    def as_dict(self):
        """
        :return: MSONAble dict
        """
        structure = Structure.from_sites(self.sites)
        return {
            "structure": structure.as_dict(),
            "spacegroup": self.spacegroup,
            "equivalent_positions": self.site_labels,
            "wyckoff_letters": self.wyckoff_letters,
        }

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: SymmetrizedStructure
        """
        return cls(
            Structure.from_dict(d["structure"]),
            spacegroup=d["spacegroup"],
            equivalent_positions=d["equivalent_positions"],
            wyckoff_letters=d["wyckoff_letters"],
        )
