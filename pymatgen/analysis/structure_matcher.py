#!/usr/bin/env python

"""
This module provides classes to perform fitting of structures.
"""

from __future__ import division

__author__ = "Stephen Dacek, William Davidson Richards, Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Stephen Dacek"
__email__ = "sdacek@mit.edu"
__status__ = "Beta"
__date__ = "Dec 3, 2012"

from pymatgen.core.structure_modifier import BasisChange
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.core.lattice import Lattice
from pymatgen.transformations.standard_transformations import\
    PrimitiveCellTransformation
from pymatgen.util.coord_utils import find_in_coord_list_pbc

import numpy as np
import itertools


class StructureMatcher(object):
    """
    Class to match structures by similarity.

    Algorithm:

    1. Given two structures: s1 and s2
    2. Optional: Reduce to primitive cells.
    3. If the number of sites do not match, return False
    4. Reduce to s1 and s2 to Niggli Cells
    5. Optional: Scale s1 and s2 to same volume.
    6. Find all possible lattice vectors for s2 within shell of ltol.
    7. For s1, translate an atom in the smallest set to the origin
    8. For s2: find all valid lattices from permutations of the list
       of lattice vectors (invalid if: det(Lattice Matrix) < half
       volume of original s2 lattice)
    9. For each valid lattice:

        a. If the lattice angles of are within tolerance of s1,
           basis change s2 into new lattice.
        b. For each atom in the smallest set of s2:

            i. Translate to origin and compare sites in structure within
               stol.
            ii. If true: break and return true
    """

    def __init__(self, ltol=0.2, stol=.4, angle_tol=5, primitive_cell=True,
                 scale=True, comparison_function=None):
        """
        Args:
            ltol:
                Fractional length tolerance. Default is 0.2.
            stol:
                Site tolerance in Angstrom. Default is 0.4 Angstrom.
            angle_tol:
                Angle tolerance in degrees. Default is 5 degrees.
            primitive_cell:
                If true: input structures will be reduced to primitive
                cells prior to matching. Default to True.
            scale:
                Input structures are scaled to equivalent volume if true;
                For exact matching, set to False.
            comparison_function:
                function declaring equivalency of sites.
                Default is None, implies rigid species mapping.
        """

        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol

        self._comparison_function = comparison_function if\
        comparison_function is not None else lambda sp1, sp2: sp1 == sp2

        self._primitive_cell = primitive_cell
        self._scale = scale

    def _basis_change(self, coords, new_lattice):
        #Convert cartesian coordinates to new basis
        frac_coords = []
        for i in range(len(coords)):
            frac_coords.append(new_lattice.get_fractional_coords(coords[i]))
        return frac_coords

    def _cmp_struct(self, s1, s2, frac_tol):
        #compares the fractional coordinates
        for s1_coords, s2_coords in zip(s1, s2):
            #Available vectors
            avail = [1] * len(s1_coords)
            for coord in s1_coords:
                ind = find_in_coord_list_pbc(s2_coords, coord, frac_tol)
                #if more than one match found, take closest
                if len(ind) > 1:
                    fcoords = np.tile(coord, (len(s2_coords[ind]), 1))
                    fdist = np.mod(s2_coords[ind], 1) - np.mod(fcoords, 1)
                    fdist -= np.round(fdist)
                    avail[np.where(
                        np.all(np.abs(fdist) == np.min(np.abs(fdist)), axis=1))[
                          0]] = 0
                elif len(ind) and avail[ind]:
                    avail[ind] = 0
                else:
                    return False
        return True

    def fit(self, struct1, struct2):
        """
        Fit two structures.

        Args:
            struct1:
                1st structure
            struct2:
                2nd structure

        Returns:
            True if the structures are the same, else False.
        """
        ltol = self.ltol
        stol = self.stol
        angle_tol = self.angle_tol

        #primitive cell transformation - Needs work
        if self._primitive_cell:
            prim = PrimitiveCellTransformation()
            struct1 = prim.apply_transformation(struct1)
            struct2 = prim.apply_transformation(struct2)
        # Same number of sites
        if struct1.num_sites != struct2.num_sites:
            return False

        #compute niggli lattices,
        nl1 = struct1.lattice.get_niggli_reduced_lattice()
        nl2 = struct2.lattice.get_niggli_reduced_lattice()
        struct1 = BasisChange(struct1, nl1).modified_structure
        struct2 = BasisChange(struct2, nl2).modified_structure

        #rescale lattice to same volume
        if self._scale:
            se1 = StructureEditor(struct1)
            nl1 = Lattice(nl1.matrix * ((nl2.volume / nl1.volume) ** (1.0 / 6)))
            se1.modify_lattice(nl1)
            struct1 = se1.modified_structure
            se2 = StructureEditor(struct2)
            nl2 = Lattice(nl2.matrix * ((nl1.volume / nl2.volume) ** (1.0 / 6)))
            se2.modify_lattice(nl2)
            struct2 = se2.modified_structure
            #Volume to determine invalid lattices
        halfs2vol = nl2.volume / 2

        #fractional tolerance of atomic positions
        frac_tol = np.array([stol / i for i in struct1.lattice.abc])

        #get possible new lattice vectors
        ds = Structure(struct2.lattice, ['X'], [[0, 0, 0]])
        nv = []
        for i in range(3):
            l = struct1.lattice.abc[i]
            nvi = []
            vs = ds.get_neighbors_in_shell([0, 0, 0], l, ltol * l)
            for site, dist in vs:
                nvi.append(site.coords)
            nv.append(nvi)
            #generate structure coordinate lists
        species_list = []
        s1 = []
        for site in struct1:
            ind = None
            for i, species in enumerate(species_list):
                if self._comparison_function(site.species_and_occu, species):
                    ind = i
                    break

            if ind is not None:
                s1[ind] = np.append(s1[ind], [site.frac_coords], axis=0)
            else:
                s1.append(np.array([site.frac_coords]))
                species_list.append(site.species_and_occu)

        zipped = sorted(zip(s1, species_list), key=lambda x: len(x[0]))
        s1 = [x[0] for x in zipped]
        species_list = [x[1] for x in zipped]
        s2_cart = [[] for i in s1]

        for site in struct2:
            ind = None
            for i, species in enumerate(species_list):
                if self._comparison_function(site.species_and_occu, species):
                    ind = i
                    break
            #get cartesian coords for s2, if no site match found return false
            if ind is None:
                return False
            if s2_cart[ind] != []:
                s2_cart[ind] = np.append(s2_cart[ind], [site.coords], axis=0)
            else:
                s2_cart[ind] = np.array([site.coords])

        #translate s1
        s1_translation = s1[0][0]

        for i in range(len(species_list)):
            s1[i] = (s1[i] - s1_translation) % 1

        #do permutations of vectors, check for equality
        for a, b, c in itertools.product(nv[0], nv[1], nv[2]):
            #invalid lattice
            if np.linalg.det([a, b, c]) < halfs2vol:
                continue

            nl = Lattice([a, b, c])
            if np.allclose(nl.angles, struct1.lattice.angles, rtol=0,
                           atol=angle_tol):
                #Basis Change into new lattice
                s2 = self._basis_change(s2_cart, nl)

                for coord in s2[0]:
                    t_s2 = []
                    for coords in s2:
                        t_s2.append((coords - coord) % 1)
                    if self._cmp_struct(s1, t_s2, frac_tol):
                        return True
        return False

    def find_indexes(self, s_list, group_list):
        """
        Given a list of structures, return list of indices where each
        structure appears in group_list.

        Args:
            s_list:
                list of structures to check
            group_list:
                list to find structures in
        """
        inds = [-1] * len(s_list)
        for j in range(len(s_list)):
            for i in range(len(group_list)):
                if len(np.where(s_list[j] in group_list[i])[0]):
                    inds[j] = i
                    break
        return inds

    def group_structures(self, s_list):
        """
        Given a list of structures, use fit to group
        them by structural equality.

        Args:
            s_list: List of structures to be grouped

        Returns:
            A list of lists of matched structures
            Assumption: if s1=s2 and s2=s3, then s1=s3
            This may not be true for small tolerances.
        """
        group_list = [[s_list[0]]]

        for i, j in itertools.combinations(range(len(s_list)), 2):
            s1_ind, s2_ind = self.find_indexes([s_list[i], s_list[j]],
                                               group_list)

            if s2_ind == -1 and self.fit(s_list[i], s_list[j]):
                group_list[s1_ind].append(s_list[j])
            elif (j - i) == 1 and s2_ind == -1:
                group_list.append([s_list[j]])

        return group_list
