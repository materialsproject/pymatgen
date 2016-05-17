# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module implements symmetry-related structure forms.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 9, 2012"

import numpy as np
from pymatgen.core.structure import Structure


class SymmetrizedStructure(Structure):
    """
    This class represents a symmetrized structure, i.e. a structure
    where the spacegroup and symmetry operations are defined. This class is
    typically not called but instead is typically obtained by calling
    pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetrized_structure.

    Args:
        structure (Structure): Original structure
        spacegroup (Spacegroup): An input spacegroup from SpacegroupAnalyzer.
        equivalent_positions: Equivalent positions from SpacegroupAnalyzer.

    .. attribute: equivalent_indices

        indices of structure grouped by equivalency
    """

    def __init__(self, structure, spacegroup, equivalent_positions):
        super(SymmetrizedStructure, self).__init__(
            structure.lattice, [site.species_and_occu for site in structure],
            structure.frac_coords, site_properties=structure.site_properties)
        self._spacegroup = spacegroup
        u, inv = np.unique(equivalent_positions, return_inverse=True)
        self.site_labels = equivalent_positions
        self.equivalent_indices = [[] for i in range(len(u))]
        self._equivalent_sites = [[] for i in range(len(u))]
        for i, inv in enumerate(inv):
            self.equivalent_indices[inv].append(i)
            self._equivalent_sites[inv].append(self.sites[i])

    @property
    def spacegroup(self):
        return self._spacegroup

    @property
    def equivalent_sites(self):
        """
        All the sites grouped by symmetry equivalence in the form of [[sites
        in group1], [sites in group2], ...]
        """
        return self._equivalent_sites

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
