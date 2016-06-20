# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines classes that set the molecular topology i.e atoms, bonds,
angles and dihedrals
"""

import itertools
from pymatgen.core.bonds import CovalentBond

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'
__credits__ = 'Brandon Wood'


class Topology(object):
    """
    Args:
        atoms (list): map atom names to force field(ff) atom name,
            [['c', 'c1'],...]
        bonds (list): List of bonds,
            [[i,j, bond_type], ... ] where i, j are integer(starts from 1)
            atom ids in the molecules and bond_type = (ff atomname_i, ff atomname_j)
        angles (list): List of angles,
            [[i,j,k, angle_type], ... ],
            angle_type = (ff atomname_i, ff atomname_j, ff atomname_k)
        charges (list): List of charges, [0.4, 0.7, ... ]
        dihedrals (list): List of dihedrals,
            [[i,j,k,l, dihedral_type], ... ]
            dihedral_type = (ff atomname_i, ff atomname_j, ff atomname_k, ff atomname_l)
        imdihedrals (list): List of improper dihedrals,
            [['i,j,k,l, dihedral_type], ... ]
    """

    def __init__(self, atoms, bonds, angles, charges=None, dihedrals=None,
                 imdihedrals=None):
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.charges = [] if charges is None else charges
        self.dihedrals = dihedrals
        self.imdihedrals = imdihedrals

    @staticmethod
    def from_molecule(molecule, tol=0.1, ff_map="ff_map"):
        """
        Return Topology object from molecule. Charges are also set if the
        molecule has 'charge' site property.

        Args:
            molecule (Molecule)
            tol (float): Relative tolerance to test in determining the bonds
                in the molecule. Basically, the code checks if the distance
                between the sites is less than (1 + tol) * typical bond
                distances. Defaults to 0.1, i.e., 10% longer.
            ff_map (string): Ensure this site property is set for each site if
                atoms need to be mapped to its forcefield name.
                eg: Carbon atom, 'C', on different sites can be mapped to
                either 'Ce' or 'Cm' depending on how the forcefield parameters
                are set.

        Returns:
            Topology object

        """
        type_attrib = ff_map if hasattr(molecule[0], ff_map) else "specie"
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
                         (str(getattr(site1, type_attrib)),
                          str(getattr(site2, type_attrib)),
                          str(getattr(site3, type_attrib)))]
                angles.append(angle)
            else:
                for site1, site2 in itertools.product(bond1_sites,
                                                      bond2_sites):
                    if CovalentBond.is_bonded(site1, site2, tol=tol):
                        bond1_sites.remove(site1)
                        bond2_sites.remove(site2)
                        dihedral = bond1_sites + [site1,
                                                  site2] + bond2_sites
                        dihedral = [molecule.index(dihedral[0]),
                                    molecule.index(dihedral[1]),
                                    molecule.index(dihedral[2]),
                                    molecule.index(dihedral[3]),
                                    (str(getattr(dihedral[0], type_attrib)),
                                     str(getattr(dihedral[1], type_attrib)),
                                     str(getattr(dihedral[2], type_attrib)),
                                     str(getattr(dihedral[3], type_attrib)))]
                        dihedrals.append(dihedral)
                        break
        atoms = [[str(site.specie), str(getattr(site, type_attrib))]
                 for site in molecule]
        bonds = [[molecule.index(b.site1), molecule.index(b.site2),
                  (str(getattr(b.site1, type_attrib)),
                   str(getattr(b.site2, type_attrib)))] for b in bonds]
        charges = None
        if hasattr(molecule[0], "charge"):
            charges = [site.charge for site in molecule]
        return Topology(atoms, bonds, angles,
                        charges=charges, dihedrals=dihedrals)
