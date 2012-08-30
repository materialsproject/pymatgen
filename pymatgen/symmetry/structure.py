#!/usr/bin/env python

"""
This module implements symmetry-related structure forms.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 9, 2012"

import itertools

from pymatgen.core.structure import Structure


class SymmetrizedStructure(Structure):
    """
    This class represents a symmetrized structure, i.e. a structure
    where the spacegroup and symmetry operations are defined.
    """

    def __init__(self, refined_structure, spacegroup, equivalent_positions):
        Structure.__init__(self, refined_structure.lattice,
                           [site.species_and_occu
                            for site in refined_structure],
                           refined_structure.frac_coords)

        self._spacegroup = spacegroup
        site_map = zip(self._sites, equivalent_positions)
        site_map = sorted(site_map, key=lambda x: x[1])
        self._equivalent_sites = [[x[0] for x in g] \
                                  for k, g in itertools.groupby(site_map,
                                                        key=lambda x: x[1])]

    @property
    def equivalent_sites(self):
        return self._equivalent_sites

    def find_equivalent_sites(self, site):
        """
        Finds all symmetrically equivalent sites for a particular site

        Args:
            site:
                A site in the structure

        Returns:
            A list of all symmetrically equivalent sites.
        """
        for sites in self.equivalent_sites:
            if site in sites:
                return sites

        raise ValueError("Site not in structure")
