# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module defines classes that set the force field parameters for the bonds,
angles and dihedrals.

TODO: complete the class by defining interfaces to AMBER and CHARMM forcefields
"""

from collections import OrderedDict

from monty.json import MSONable

__author__ = 'Kiran Mathew, Navnidhi Rajput'


class ForceField(MSONable):
    """
    Stores force field information.

    Args:
        atoms (OrderedDict): store atomic mass for each atom name.
            { "atom name": atom mass, ... }
        bonds (OrderedDict): store the bond distance (A) and spring constant (
            Kcal/molA2) for each bond.
            { ("atom name1", "atom name2"): [spring const, distance], ... }
        angles (OrderedDict): store the bond angle and spring constant
            (Kcal/mol*radian2).
            { ("atom name1", "atom name2", "atom name3"): [spring const, angle], ... }
        dihedrals (OrderedDict): store the magnitude of torsion (Kcal/mol).
            { ("atom name1", "atom name2", "atom name3", "atom name4"): [
            function type, value, angle], ... }
        imdihedrals (OrderedDict): store improper dihedral information.
            similar to dihedrals but the gaff atom name1 and gaff atom name2 are marked 'X'
        vdws (OrderedDict): store the van der waal radius (A) and van der wall
            depth for a given atom (Kcal/mol). Lennard-Jones parameters.
            { "atom name": [sigma, epsilon], ... }
    """

    def __init__(self, atoms, bonds, angles, dihedrals=None, imdihedrals=None,
                 vdws=None, masses=None, charges=None):
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.imdihedrals = imdihedrals
        self.vdws = vdws
        self.masses = masses
        self.charges = OrderedDict() if charges is None else charges
        self.atomindex_to_atomname = OrderedDict()
        self.atomindex_to_gaffname = OrderedDict()
        self.atomname_to_gaffname = OrderedDict()
        self.gaffname_to_atomindex = OrderedDict()
        # self.set_atom_mappings()

    # TODO: add interfaces to AMBER and CHARMM force fields