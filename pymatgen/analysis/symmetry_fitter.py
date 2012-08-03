#!/usr/bin/env python

'''
This module implements an experimental symmetry fitter class.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 5, 2012"

import itertools
import logging

from pymatgen.symmetry.finder import SymmetryFinder

logger = logging.getLogger(__name__)


class SymmetryFitter(object):
    """
    The SymmetryFitter class is used to compare a sequence of structures and
    group them according to whether they are similar under a particular
    spacegroup.

    This class is still *EXPERIMENTAL*. Use with care.

    Attributes:

    .. attribute:: unique_groups

        All unique groups of structures under the spacegroup, given as
        [[struct1, struct2,...], ...]
    """

    def __init__(self, structures, spacegroup, symm_prec=0.1):
        """
        Args:
            structures:
                Sequence of structures to test.
            spacegroup:
                A spacegroup to test the structures.
            symm_prec:
                The symmetry precision to test with.
        """
        structure_symm = {}
        self.symm_prec = symm_prec
        self.spacegroup = spacegroup
        logger.debug("Computing spacegroups...")
        for i, s in enumerate(structures):
            finder = SymmetryFinder(s, symm_prec)
            structure_symm[s] = finder.get_spacegroup_number()
            logger.debug("Structure {} has spacegroup {}"
                         .format(i, structure_symm[s]))

        sorted_structures = sorted(structures, key=lambda s:-structure_symm[s])
        unique_groups = []
        for i, group in itertools.groupby(sorted_structures,
                                          key=lambda s: structure_symm[s]):
            logger.debug("Processing group of structures with sg number {}"
                         .format(i))
            subgroups = self._fit_group(group)
            unique_groups.extend(subgroups)
        self.unique_groups = unique_groups

    def num_groups(self):
        return len(self.unique_groups)

    def get_unique_structures(self):
        return (s[0] for s in self.unique_groups)

    def _fit_group(self, group):
        subgroups = []
        all_structures = list(group)
        logger.debug("{} structures in group".format(len(all_structures)))
        while len(all_structures) > 0:
            fixed = all_structures[0]
            subgroup = [fixed]
            for to_fit in all_structures[1:]:
                if self.spacegroup.are_symmetrically_equivalent(fixed, to_fit,
                                                    symm_prec=self.symm_prec):
                    subgroup.append(to_fit)
            all_structures = [s for s in all_structures if s not in subgroup]
            subgroups.append(subgroup)
        logger.debug("{} subgroups".format(len(subgroups)))
        return subgroups
