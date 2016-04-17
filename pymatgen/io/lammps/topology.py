# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module defines classes that set the molecular topology i.e atoms, bonds,
angles and dihedrals

TODO: complete the class by defining interfaces to AMBER and CHARMM topologies
"""

__author__ = 'Kiran Mathew, Navnidhi Rajput'


class Topology(object):
    """
    Args:
        atoms (list): map atom names to force field(ff) atom name, [['c', 'c1'],...]
        charges (list): List of charges, [0.4, 0.7, ... ]
        bonds (list): List of bonds, [[i,j, bond_type], ... ] where i, j are integer(starts from 1)
            atom ids in the molecules and bond_type = (ff atomname_i, ff atomname_j)
        angles (list): List of angles, [[i,j,k, angle_type], ... ],
            angle_type = (ff atomname_i, ff atomname_j, ff atomname_k)
        dihedrals (list): List of dihedrals, [[i,j,k,l, dihedral_type], ... ]
            dihedral_type = (ff atomname_i, ff atomname_j, ff atomname_k, ff atomname_l)
        imdihedrals (list): List of improper dihedrals, [['i,j,k,l, dihedral_type], ... ]

        TODO: remove the rest of the params, dont see the need
    """

    def __init__(self, atoms, bonds, angles, charges=None, dihedrals=None,
                 imdihedrals=None):
        self.atoms = atoms
        self.charges = dict() if charges is None else charges
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.imdihedrals = imdihedrals

    # TODO: add from_file interfaces to AMBER and CHARMM topologies


