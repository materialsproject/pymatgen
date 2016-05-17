# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This class implements definitions for various kinds of bonds. Typically used in
Molecule analysis.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 26, 2012"


import os
import json
import collections

from pymatgen.core.periodic_table import get_el_sp


def _load_bond_length_data():
    """Loads bond length data from json file"""
    with open(os.path.join(os.path.dirname(__file__),
                           "bond_lengths.json")) as f:
        data = collections.defaultdict(dict)
        for row in json.load(f):
            els = sorted(row['elements'])
            data[tuple(els)][row['bond_order']] = row['length']
        return data

bond_lengths = _load_bond_length_data()


class CovalentBond(object):
    """
    Defines a covalent bond between two sites.
    """

    def __init__(self, site1, site2):
        """
        Initializes a covalent bond between two sites.

        Args:
            site1 (Site): First site.
            site2 (Site): Second site.
        """
        self.site1 = site1
        self.site2 = site2

    @property
    def length(self):
        """
        Length of the bond.
        """
        return self.site1.distance(self.site2)

    @staticmethod
    def is_bonded(site1, site2, tol=0.2, bond_order=None):
        """
        Test if two sites are bonded, up to a certain limit.

        Args:
            site1 (Site): First site
            site2 (Site): Second site
            tol (float): Relative tolerance to test. Basically, the code
                checks if the distance between the sites is less than (1 +
                tol) * typical bond distances. Defaults to 0.2, i.e.,
                20% longer.
            bond_order: Bond order to test. If None, the code simply checks
                against all possible bond data. Defaults to None.

        Returns:
            Boolean indicating whether two sites are bonded.
        """
        sp1 = list(site1.species_and_occu.keys())[0]
        sp2 = list(site2.species_and_occu.keys())[0]
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
        return "Covalent bond between {} and {}".format(self.site1,
                                                        self.site2)

    def __str__(self):
        return self.__repr__()


def get_bond_length(sp1, sp2, bond_order=1):
    """
    Get the bond length between two species.

    Args:
        sp1 (Specie): First specie.
        sp2 (Specie): Second specie.
        bond_order: For species with different possible bond orders,
            this allows one to obtain the bond length for a particular bond
            order. For example, to get the C=C bond length instead of the
            C-C bond length, this should be set to 2. Defaults to 1.

    Returns:
        Bond length in Angstrom. None if no data is available.
    """

    syms = tuple(sorted([get_el_sp(sp1).symbol,
                         get_el_sp(sp2).symbol]))
    if syms in bond_lengths:
        all_lengths = bond_lengths[syms]
        if bond_order:
            return all_lengths.get(bond_order)
        else:
            return all_lengths.get(1)
    return None
