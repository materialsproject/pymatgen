# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements symmetry-related structure forms.
"""

from typing import Sequence, List

import numpy as np
from pymatgen.core.structure import Structure, PeriodicSite


class SymmetrizedStructure(Structure):
    """
    This class represents a symmetrized structure, i.e. a structure
    where the spacegroup and symmetry operations are defined. This class is
    typically not called but instead is typically obtained by calling
    pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetrized_structure.

    .. attribute: equivalent_indices

        indices of structure grouped by equivalency
    """

    def __init__(self, structure: Structure,
                 spacegroup,
                 equivalent_positions: Sequence[int],
                 wyckoff_letters: Sequence[str]):
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
            structure.lattice, [site.species for site in structure],
            structure.frac_coords, site_properties=structure.site_properties)

        equivalent_indices = [[] for _ in range(len(u))]  # type: ignore
        equivalent_sites = [[] for _ in range(len(u))]  # type: ignore
        wyckoff_symbols = [[] for _ in range(len(u))]  # type: ignore
        for i, inv in enumerate(inv):
            equivalent_indices[inv].append(i)
            equivalent_sites[inv].append(self.sites[i])
            wyckoff_symbols[inv].append(wyckoff_letters[i])
        self.equivalent_indices: List[int] = equivalent_indices  # type: ignore
        self.equivalent_sites: List[PeriodicSite] = equivalent_sites  # type: ignore
        self.wyckoff_letters = wyckoff_letters
        self.wyckoff_symbols = ["%d%s" % (len(w), w[0])
                                for w in wyckoff_symbols]

    def copy(self):
        """
        :return: Copy of structure.
        """
        return self.__class__(self, spacegroup=self.spacegroup,
                              equivalent_positions=self.site_labels,
                              wyckoff_letters=self.wyckoff_letters)

    def find_equivalent_sites(self, site):
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
