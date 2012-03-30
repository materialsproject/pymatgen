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

import os
import re
import glob
import numpy as np

from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import PeriodicSite

class Spacegroup(object):
    """
    Represents a space group, which is a collection of symmetry operations
    """

    def __init__(self, int_symbol, int_number, symmops):
        self._symbol = int_symbol
        self._number = int_number
        self._symmops = symmops

    @property
    def international_symbol(self):
        return self._symbol

    @property
    def international_number(self):
        return self._number

    def are_symmetrically_equivalent(self, sites1, sites2, symprec=1e-8):
        """
        Given two sets of PeriodicSites, test if they are actually symmetrically
        equivalent under this spacegroup.  Useful, for example, if you want to
        test if selecting atoms 1 and 2 out of a set of 4 atoms are symmetrically
        the same as selecting atoms 3 and 4, etc.
        
        One use is in PartialRemoveSpecie transformation to return only symmetrically
        distinct arrangements of atoms.
        
        Args:
            sites1:
                1st set of sites
            sites2:
                2nd set of sites
            symprec:
                The tolerance in atomic distance to test if atoms are 
                symmetrically similar.
        
        Returns:
            Boolean indicating whether the two sets of sites are symmetrically
            equivalent.
        """
        def in_sites(site):
            for test_site in sites1:
                if test_site.is_periodic_image(site, symprec):
                    return True
            return False
        for op in self._symmops:
            newsites2 = [PeriodicSite(site.species_and_occu, op.operate(site.frac_coords), site.lattice) for site in sites2]
            ismapping = True
            for site in newsites2:
                if not in_sites(site):
                    ismapping = False
                    break
            if ismapping:
                return True
        return False

    @staticmethod
    def from_spacegroup_number(sgnum):
        datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sg_data')
        filename = str(sgnum).zfill(3) + "*"
        files = sorted(glob.glob(os.path.join(datadir, filename)))
        with open(files[0], "r") as fid:
            symmops = []
            rots = []
            lines = fid.readlines()
            sgname = lines[0].strip()
            for i in xrange(1, len(lines)):
                toks = re.split(",", lines[i].strip())
                if len(toks) == 3:
                    rot = np.zeros((3, 3))
                    trans = [0, 0, 0]
                    for j in xrange(3):
                        tok = toks[j]
                        m = re.search("([\+\-]*)([xyz])", tok)
                        if m:
                            factor = -1 if m.group(1) == "-" else 1
                            loc = ord(m.group(2)) - 120
                            rot[j, loc] = factor
                            tok = re.sub("([\+\-]*)([xyz])", "", tok)
                            if tok.strip() != '':
                                trans[j] = eval(tok)
                        rots.append(rot)
                    symmops.append(SymmOp.from_rotation_matrix_and_translation_vector(rot, trans))
            return Spacegroup(sgname, sgnum, symmops)

    def __str__(self):
        return "{} ({}) spacegroup".format(self._symbol, self._number)

if __name__ == "__main__":
    print Spacegroup.from_spacegroup_number(230)
