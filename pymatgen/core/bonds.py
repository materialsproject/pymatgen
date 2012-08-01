#!/usr/bin/env python

'''
This class implements definitions for various kinds of bonds. Typically used in
Molecule analysis.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 26, 2012"


import os
import json
import collections


def _load_bond_length_data():
    """Loads element data from json file"""
    module_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(module_dir, "bond_lengths.json")) as f:
        data = collections.defaultdict(dict)
        for row in json.load(f):
            els = sorted(row['elements'])
            data[tuple(els)][row['bond_order']] = row['length']
        return data

bond_lengths = _load_bond_length_data()


class CovalentBond(object):

    def __init__(self, site1, site2):
        self.site1 = site1
        self.site2 = site2

    @property
    def length(self):
        return self.site1.distance(self.site2)

    @staticmethod
    def is_bonded(site1, site2, tol=0.2, bond_order=None):
        """
        Test if two sites are bonded, up to a certain limit.

        Args:
            site1:
                First site
            site2:
                Second site
            tol:
                Relative tolerance to test. Basically, the code checks if the
                distance between the sites is less than (1 + tol) * typical
                bond distances. Defaults to 0.2, i.e., 20% longer.
            bond_order:
                Bond order to test. If None, the code simply checks against all
                possible bond data. Defaults to None.
        """
        sp1 = site1.species_and_occu.keys()[0]
        sp2 = site2.species_and_occu.keys()[0]
        dist = site1.distance(site2)
        syms = tuple(sorted([sp1.symbol, sp2.symbol]))
        if syms in bond_lengths:
            all_lengths = bond_lengths[syms]
            if bond_order:
                return dist < (1 + tol) * all_lengths[bond_order]
            for v in all_lengths.values():
                if dist < (1 + tol) * v:
                    return True
            return False
        raise ValueError("No bond data for elements {} - {}".format(*syms))

    def __repr__(self):
        output = ["Covalent bond"]
        output.append("between {}".format(self.site1))
        output.append("and {}".format(self.site2))
        return " ".join(output)
