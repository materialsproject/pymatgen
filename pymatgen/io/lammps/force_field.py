# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines classes that set the force field parameters for the bonds,
angles and dihedrals.
"""

from collections import defaultdict
import yaml

from monty.json import MSONable

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'
__credits__ = 'Brandon Wood'


class ForceField(MSONable):
    """
    Stores force field information.

    Args:
        atoms (Dict): store atomic mass for each atom name.
            { "custom atom name": atom name, ... }
        bonds (Dict): store the bond distance and spring constant for each bond.
            { ("atom name1", "atom name2"): [spring const, distance], ... }
        angles (Dict): store the bond angle and spring constant.
            { ("atom name1", "atom name2", "atom name3"): [spring const, angle], ... }
        dihedrals (Dict): store dihedral paramters.
            { ("atom name1", "atom name2", "atom name3", "atom name4"):
            [val1, val2, ...], ... }
        imdihedrals (Dict): store improper dihedral information. Similar to
            dihedrals.
        pairs (Dict): store pair coefficient info.
            { ("atom name1", "atom name2"): [val1, val2, ..], ... }
    """

    def __init__(self, atoms, bonds, angles, dihedrals=None, imdihedrals=None,
                 pairs=None):
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.imdihedrals = imdihedrals
        self.pairs = pairs

    @staticmethod
    def from_file(filename):
        """
        Read in forcefield parameters from yaml file and return Forcefield
        object. Basically read in atom name mappings(key='Atoms'), bond
        coefficients(key='Bond Coeffs'), angle coefficients(key='Angle
        Coeffs'), pair coefficients(key='Pair Coeffs'),
        dihedral coefficients(key='Dihedral Coeffs') and the
        improper dihedral coefficients(key='Improper Coeffs').

        Args:
            filename (string)

        Returns:
            ForceField object
        """
        with open(filename, 'r') as f:
            d = yaml.load(f)
        ff_data = defaultdict(dict)
        for coeff_key, coeff in d.items():
            for k, v in coeff.items():
                tokens = k.split("-")
                key = tuple(tokens) if len(tokens) > 1 else k
                ff_data[coeff_key][key] = v
        pairs = ff_data.get("Pair Coeffs", None)
        if pairs:
            if len(ff_data["Atoms"]) != len(pairs):
                raise ValueError("Number of pairs coefficient parmaters > "
                                 "the number of atome types. Parameters i != j "
                                 "pairs cannot be set in the data file")
        return ForceField(ff_data["Atoms"],
                          ff_data["Bond Coeffs"],
                          ff_data["Angle Coeffs"],
                          dihedrals=ff_data.get("Dihedral Coeffs", None),
                          imdihedrals=ff_data.get("Improper Coeffs", None),
                          pairs=pairs)
