#!/usr/bin/env python

"""
This module implements symmetry-related structure forms.
"""

from __future__ import division

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
    pymatgen.symmetry.SymmetryFinder.get_symmetrized_structure.
    
    .. attribute: equivalent_indices
        
        indices of structure grouped by equivalency
    """

    def __init__(self, structure, spacegroup, equivalent_positions):
        """
        Args:
            structure:
                Original structure
            spacegroup:
                An input spacegroup from SymmetryFinder.
            equivalent_positions:
                Equivalent positions from SymmetryFinder.
        """
        Structure.__init__(self, structure.lattice,
                           [site.species_and_occu
                            for site in structure],
                           structure.frac_coords,
                           site_properties=structure.site_properties)
        
        self._spacegroup = spacegroup
        u, inv = np.unique(equivalent_positions, return_inverse = True)
        self.equivalent_indices = [[] for i in xrange(len(u))]
        self._equivalent_sites = [[] for i in xrange(len(u))]
        for i, inv in enumerate(inv):
            self.equivalent_indices[inv].append(i)
            self._equivalent_sites[inv].append(self.sites[i])

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
            site:
                A site in the structure

        Returns:
            A list of all symmetrically equivalent sites.
        """
        for sites in self.equivalent_sites:
            if site in sites:
                return sites

        raise ValueError("Site not in structure")
