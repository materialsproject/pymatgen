#!/usr/bin/env python

'''
This module implements a basic Spacegroup class
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 9, 2012"

from pymatgen.core.sites import PeriodicSite


class Spacegroup(object):
    """
    Represents a space group, which is a collection of symmetry operations
    """

    def __init__(self, int_symbol, int_number, symmops):
        """
        Args:
            int_symbol:
                The international symbol of the spacegroup.
            int_number:
                The international number of the spacegroup.
            symmops:
                The symmetry operations associated with the spacegroup.
        """
        self.int_symbol = int_symbol
        self.int_number = int_number
        self.symmops = symmops

    def are_symmetrically_equivalent(self, sites1, sites2, symm_prec=1e-3):
        """
        Given two sets of PeriodicSites, test if they are actually
        symmetrically equivalent under this space group.  Useful, for example,
        if you want to test if selecting atoms 1 and 2 out of a set of 4 atoms
        are symmetrically the same as selecting atoms 3 and 4, etc.

        One use is in PartialRemoveSpecie transformation to return only
        symmetrically distinct arrangements of atoms.

        Args:
            sites1:
                1st set of sites
            sites2:
                2nd set of sites
            symm_prec:
                The tolerance in atomic distance to test if atoms are
                symmetrically similar.

        Returns:
            Boolean indicating whether the two sets of sites are symmetrically
            equivalent.
        """
        def in_sites(site):
            for test_site in sites1:
                if test_site.is_periodic_image(site, symm_prec, False):
                    return True
            return False
        for op in self.symmops:
            newsites2 = [PeriodicSite(site.species_and_occu,
                                      op.operate(site.frac_coords),
                                      site.lattice) for site in sites2]
            ismapping = True
            for site in newsites2:
                if not in_sites(site):
                    ismapping = False
                    break
            if ismapping:
                return True
        return False

    def __str__(self):
        return "{} ({}) spacegroup".format(self._symbol, self._number)
