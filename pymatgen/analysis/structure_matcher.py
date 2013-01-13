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


import numpy as np
import itertools
import abc

from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.core.lattice import Lattice
from pymatgen.util.coord_utils import find_in_coord_list_pbc
from pymatgen.util.coord_utils import pbc_diff
from pymatgen.core.composition import Composition


class AbstractComparator(object):
    """
    Abstract Comparator class. A Comparator defines how sites are compared in
    a structure.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def are_equal(self, sp1, sp2):
        """
        Defines how the species of two sites are considered equal. For
        example, one can consider sites to have the same species only when
        the species are exactly the same, i.e., Fe2+ matches Fe2+ but not
        Fe3+. Or one can define that only the element matters,
        and all oxidation state information are ignored.

        Args:
            sp1:
                First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2:
                Second species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.

        Returns:
            Boolean indicating whether species are considered equal.
        """
        return

    @abc.abstractmethod
    def get_structure_hash(self, structure):
        """
        Defines a hash for structures. This allows structures to be grouped
        efficiently for comparison. For example, in exact matching,
        you should only try to perform matching if structures have the same
        reduced formula (structures with different formulas can't possibly
        match). So the reduced_formula is a good hash. The hash function
        should be relatively fast to compute relative to the actual matching.

        Args:
            structure:
                A structure

        Returns:
            A hashable object. Examples can be string formulas, etc.
        """
        return


class SpeciesComparator(AbstractComparator):
    """
    A Comparator that matches species exactly. The default used in
    StructureMatcher.
    """

    def are_equal(self, sp1, sp2):
        """
        True if species are exactly the same, i.e., Fe2+ == Fe2+ but not Fe3+.

        Args:
            sp1:
                First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2:
                Second species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.

        Returns:
            Boolean indicating whether species are equal.
        """
        return sp1 == sp2

    def get_structure_hash(self, structure):
        """
        Hash for structure.

        Args:
            structure:
                A structure

        Returns:
            Reduced formula for the structure is used as a hash for the
            SpeciesComparator.
        """
        return structure.composition.reduced_formula


class ElementComparator(AbstractComparator):
    """
    A Comparator that matches elements. i.e. oxidation states are
    ignored.
    """

    def are_equal(self, sp1, sp2):
        """
        True if element:amounts are exactly the same, i.e.,
        oxidation state is not considered.

        Args:
            sp1:
                First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2:
                Second species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.

        Returns:
            Boolean indicating whether species are the same based on element
            and amounts.
        """
        comp1 = Composition(sp1)
        comp2 = Composition(sp2)
        return comp1.get_el_amt_dict() == comp2.get_el_amt_dict()

    def get_structure_hash(self, structure):
        """
        Hash for structure.

        Args:
            structure:
                A structure

        Returns:
            Reduced formula for the structure is used as a hash for the
            SpeciesComparator.
        """
        return structure.composition.reduced_formula


class FrameworkComparator(AbstractComparator):
    """
    A Comparator that matches sites, regardless of species.
    """

    def are_equal(self, sp1, sp2):
        """
        True if there are atoms on both sites.

        Args:
            sp1:
                First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2:
                Second species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.

        Returns:
            True always
        """
        return True

    def get_structure_hash(self, structure):
        """
        Hash for structure.

        Args:
            structure:
                A structure

        Returns:
            Number of atoms is a good hash for simple framework matching.
        """
        return len(structure)


class StructureMatcher(object):
    """
    Class to match structures by similarity.

    Algorithm:

    1. Given two structures: s1 and s2
    2. Optional: Reduce to primitive cells.
    3. If the number of sites do not match, return False
    4. Reduce to s1 and s2 to Niggli Cells
    5. Optional: Scale s1 and s2 to same volume.
    6. Optional: Remove oxidation states associated with sites
    7. Find all possible lattice vectors for s2 within shell of ltol.
    8. For s1, translate an atom in the smallest set to the origin
    9. For s2: find all valid lattices from permutations of the list
       of lattice vectors (invalid if: det(Lattice Matrix) < half
       volume of original s2 lattice)
    10. For each valid lattice:

        a. If the lattice angles of are within tolerance of s1,
           basis change s2 into new lattice.
        b. For each atom in the smallest set of s2:

            i. Translate to origin and compare sites in structure within
               stol.
            ii. If true: break and return true
    """

    def __init__(self, ltol=0.2, stol=0.5, angle_tol=5, primitive_cell=True,
                 scale=True, comparator=SpeciesComparator()):
        """
        Args:
            ltol:
                Fractional length tolerance. Default is 0.2.
            stol:
                Site tolerance in Angstrom. Default is 0.5 Angstrom.
            angle_tol:
                Angle tolerance in degrees. Default is 5 degrees.
            primitive_cell:
                If true: input structures will be reduced to primitive
                cells prior to matching. Default to True.
            scale:
                Input structures are scaled to equivalent volume if true;
                For exact matching, set to False.
            comparator:
                A comparator object implementing an equals method that declares
                declaring equivalency of sites. Default is
                SpeciesComparator, which implies rigid species
                mapping, i.e., Fe2+ only matches Fe2+ and not Fe3+.

                Other comparators are provided, e.g., ElementComparator which
                matches only the elements and not the species.

                The reason why a comparator object is used instead of
                supplying a comparison function is that it is not possible to
                pickle a function, which makes it otherwise difficult to use
                StructureMatcher with Python's multiprocessing.
        """
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol
        self._comparator = comparator
        self._primitive_cell = primitive_cell
        self._scale = scale

    def _cmp_struct(self, s1, s2, nl, frac_tol):
        #compares the fractional coordinates
        for s1_coords, s2_coords in zip(s1, s2):
            #Available vectors
            avail = [1] * len(s1_coords)
            for coord in s1_coords:
                ind = find_in_coord_list_pbc(s2_coords, coord, frac_tol)
                #if more than one match found, take closest
                if len(ind) > 1:
                    #only check against available vectors
                    ind = [i for i in ind if avail[i]]
                    if len(ind) > 1:
                        #get cartesian distances from periodic distances
                        pb_dists = np.array([pbc_diff(s2_coords[i], coord)
                                             for i in ind])
                        carts = nl.get_cartesian_coords(pb_dists)
                        dists = np.array([np.linalg.norm(carts[i])
                                          for i in range(len(ind))])
                        #use smallest distance
                        ind = np.where(dists == np.min(dists))[0][0]
                        avail[ind] = 0
                    elif len(ind):
                        avail[ind[0]] = 0
                    else:
                        return False
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
            True if the structures are the equivalent, else False.
        """
        ltol = self.ltol
        stol = self.stol
        angle_tol = self.angle_tol
        comparator = self._comparator

        if comparator.get_structure_hash(struct1) != \
            comparator.get_structure_hash(struct2):
            return False

        #primitive cell transformation
        if self._primitive_cell and struct1.num_sites != struct2.num_sites:
            struct1 = struct1.get_primitive_structure()
            struct2 = struct2.get_primitive_structure()

        # Same number of sites
        if struct1.num_sites != struct2.num_sites:
            return False

        # Get niggli reduced cells. Though technically not necessary, this
        # minimizes cell lengths and speeds up the matching of skewed
        # cells considerably.
        struct1 = struct1.get_reduced_structure(reduction_algo="niggli")
        struct2 = struct2.get_reduced_structure(reduction_algo="niggli")

        nl1 = struct1.lattice
        nl2 = struct2.lattice

        #rescale lattice to same volume
        if self._scale:
            scale_vol = (nl2.volume / nl1.volume) ** (1 / 6)
            se1 = StructureEditor(struct1)
            nl1 = Lattice(nl1.matrix * scale_vol)
            se1.modify_lattice(nl1)
            struct1 = se1.modified_structure
            se2 = StructureEditor(struct2)
            nl2 = Lattice(nl2.matrix / scale_vol)
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
            vs = ds.get_neighbors_in_shell([0, 0, 0], l, ltol * l)
            nvi = [site.coords for site, dist in vs]
            nv.append(nvi)

        #generate structure coordinate lists
        species_list = []
        s1 = []
        for site in struct1:
            found = False
            for i, species in enumerate(species_list):
                if comparator.are_equal(site.species_and_occu, species):
                    found = True
                    s1[i].append(site.frac_coords)
                    break
            if not found:
                s1.append([site.frac_coords])
                species_list.append(site.species_and_occu)

        zipped = sorted(zip(s1, species_list), key=lambda x: len(x[0]))

        s1 = [x[0] for x in zipped]
        species_list = [x[1] for x in zipped]
        s2_cart = [[] for i in s1]

        for site in struct2:
            found = False
            for i, species in enumerate(species_list):
                if comparator.are_equal(site.species_and_occu, species):
                    found = True
                    s2_cart[i].append(site.coords)
                    break
            #if no site match found return false
            if not found:
                return False

        #translate s1
        s1_translation = s1[0][0]

        for i in range(len(species_list)):
            s1[i] = np.mod(s1[i] - s1_translation, 1)
        #do permutations of vectors, check for equality
        s1_angles = struct1.lattice.angles
        for a, b, c in itertools.product(nv[0], nv[1], nv[2]):
            #invalid lattice
            if np.abs(np.linalg.det([a, b, c])) < halfs2vol:
                continue

            nl = Lattice([a, b, c])
            if np.allclose(nl.angles, s1_angles, rtol=0, atol=angle_tol):
                #Basis Change into new lattice
                s2 = [nl.get_fractional_coords(c) for c in s2_cart]
                for coord in s2[0]:
                    t_s2 = [np.mod(coords - coord, 1) for coords in s2]
                    if self._cmp_struct(s1, t_s2, nl, frac_tol):
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
            s_list:
                List of structures to be grouped

        Returns:
            A list of lists of matched structures
            Assumption: if s1=s2 and s2=s3, then s1=s3
            This may not be true for small tolerances.
        """
        #Use structure hash to pre-group structures.
        sorted_s_list = sorted(s_list, key=self._comparator.get_structure_hash)
        all_groups = []

        #For each pre-grouped list of structures, perform actual matching.
        for k, g in itertools.groupby(sorted_s_list,
                                      key=self._comparator.get_structure_hash):
            g = list(g)
            group_list = [[g[0]]]
            for i, j in itertools.combinations(range(len(g)), 2):
                s1_ind, s2_ind = self.find_indexes([g[i], g[j]], group_list)
                if s2_ind == -1 and self.fit(g[i], g[j]):
                    group_list[s1_ind].append(g[j])
                elif (j - i) == 1 and s2_ind == -1:
                    group_list.append([g[j]])
            all_groups.extend(group_list)

        return all_groups
