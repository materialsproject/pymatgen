# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module defines classes that set the force field parameters for the bonds,
angles and dihedrals.
"""

from collections import OrderedDict, defaultdict
import yaml

from monty.json import MSONable

__author__ = 'Kiran Mathew'


class ForceField(MSONable):
    """
    Stores force field information.

    Args:
        atoms (Dict): store atomic mass for each atom name.
            { "atom name": atom mass, ... }
        bonds (Dict): store the bond distance (A) and spring constant (
            Kcal/molA2) for each bond.
            { ("atom name1", "atom name2"): [spring const, distance], ... }
        angles (Dict): store the bond angle and spring constant
            (Kcal/mol*radian2).
            { ("atom name1", "atom name2", "atom name3"): [spring const, angle], ... }
        dihedrals (Dict): store the magnitude of torsion (Kcal/mol).
            { ("atom name1", "atom name2", "atom name3", "atom name4"): [
            function type, value, angle], ... }
        imdihedrals (Dict): store improper dihedral information.
            similar to dihedrals but the gaff atom name1 and gaff atom name2
            are marked 'X'
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
        coefficients(key='Bond Coeffs'), angle coefficients(key='Angle Coeffs')
        and the dihedral coefficients(key='Dihedral Coeffs').

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
        return ForceField(ff_data["Atoms"], ff_data["Bond Coeffs"],
                          ff_data["Angle Coeffs"], ff_data["Dihedral Coeffs"],
                          ff_data["Improper Coeffs"], ff_data["Pair Coeffs"])
