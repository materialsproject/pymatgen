# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

"""
This module defines classes that set the molecular topology i.e atoms, bonds,
angles and dihedrals
"""

import itertools
from pymatgen.core.bonds import CovalentBond

__author__ = 'Kiran Mathew'


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

    """

    def __init__(self, atoms, bonds, angles, charges=None, dihedrals=None,
                 imdihedrals=None):
        self.atoms = atoms
        self.charges = dict() if charges is None else charges
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.imdihedrals = imdihedrals

    @staticmethod
    def from_molecule(molecule, tol=0.1):
        """
        Return Topology object from molecule.

        Args:
            molecule (Molecule)
            tol (float): Relative tolerance to test. Basically, the code
                checks if the distance between the sites is less than (1 +
                tol) * typical bond distances. Defaults to 0.1, i.e.,
                10% longer.

        Returns:
            Topology object

        """
        bonds = molecule.get_covalent_bonds(tol=tol)
        angles = []
        dihedrals = []
        for bond1, bond2 in itertools.combinations(bonds, 2):
            bond1_sites = [bond1.site1, bond1.site2]
            bond2_sites = [bond2.site1, bond2.site2]
            s_bond1 = set(bond1_sites)
            s_bond2 = set(bond2_sites)
            common_site = s_bond1.intersection(s_bond2)
            if common_site:
                site1 = (s_bond1 - s_bond2).pop()
                site2 = common_site.pop()
                site3 = (s_bond2 - s_bond1).pop()
                angle = [molecule.index(site1), molecule.index(site2),
                         molecule.index(site3),
                         (str(site1.specie), str(site2.specie),
                          str(site3.specie))]
                angles.append(angle)
            else:
                for site1, site2 in itertools.product(bond1_sites,
                                                      bond2_sites):
                    if CovalentBond.is_bonded(site1, site2):
                        bond1_sites.remove(site1)
                        bond2_sites.remove(site2)
                        dihedral = bond1_sites + [site1,
                                                  site2] + bond2_sites
                        dihedral = [molecule.index(dihedral[0]),
                                    molecule.index(dihedral[1]),
                                    molecule.index(dihedral[2]),
                                    molecule.index(dihedral[3]),
                                    (str(dihedral[0].specie),
                                     str(dihedral[1].specie),
                                     str(dihedral[2].specie),
                                     str(dihedral[3].specie))]
                        dihedrals.append(dihedral)
                        break
        atoms = [[str(site.specie), str(site.specie)] for site in molecule]
        bonds = [[molecule.index(b.site1), molecule.index(b.site2),
                  (str(b.site1.specie), str(b.site2.specie))] for b in
                 bonds]
        return Topology(atoms, bonds, angles, dihedrals=dihedrals)
