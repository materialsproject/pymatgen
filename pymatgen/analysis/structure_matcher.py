#!/usr/bin/env python

"""
This module provides classes to perform fitting of structures.
"""

from __future__ import division

__author__ = "Stephen Dacek, William Davidson Richards, Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Stephen Dacek"
__email__ = "sdacek@mit.edu"
__status__ = "Production"
__date__ = "Dec 3, 2012"

import numpy as np
import itertools
import abc

from pymatgen.serializers.json_coders import MSONable
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.composition import Composition
from pymatgen.optimization.linear_assignment import LinearAssignment
from pymatgen.util.coord_utils import get_points_in_sphere_pbc, \
    pbc_shortest_vectors


class AbstractComparator(MSONable):
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

    @classmethod
    def from_dict(cls, d):
        for trans_modules in ['structure_matcher']:
            mod = __import__('pymatgen.analysis.' + trans_modules,
                             globals(), locals(), [d['@class']], -1)
            if hasattr(mod, d['@class']):
                trans = getattr(mod, d['@class'])
                return trans()
        raise ValueError("Invalid Comparator dict")

    @property
    def to_dict(self):
        return {"version": __version__, "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}


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


class StructureMatcher(MSONable):
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

            i. Translate to origin and compare fractional sites in
            structure within a fractional tolerance.
            ii. If true:

                ia. Convert both lattices to cartesian and place
                both structures on an average lattice
                ib. Compute and return the average and max rms
                displacement between the two structures normalized
                by the average free length per atom

                if fit function called:
                    if normalized max rms displacement is less than
                    stol. Return True

                if get_rms_dist function called:
                    if normalized average rms displacement is less
                    than the stored rms displacement, store and
                    continue. (This function will search all possible
                    lattices for the smallest average rms displacement
                    between the two structures)

    """

    def __init__(self, ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True,
                 scale=True, attempt_supercell=False,
                 comparator=SpeciesComparator()):
        """
        Args:
            ltol:
                Fractional length tolerance. Default is 0.2.
            stol:
                Site tolerance. Defined as the fraction of the
                average free length per atom := ( V / Nsites ) ** (1/3)
                Default is 0.3.
            angle_tol:
                Angle tolerance in degrees. Default is 5 degrees.
            primitive_cell:
                If true: input structures will be reduced to primitive
                cells prior to matching. Default to True.
            scale:
                Input structures are scaled to equivalent volume if true;
                For exact matching, set to False.
            attempt_supercell:
                If set to True and number of sites in cells differ
                after a primitive cell reduction (divisible by an integer)
                attempts to generate a supercell transformation of the
                smaller cell which is equivalent to the larger structure.
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
        self._supercell = attempt_supercell

    def _get_lattices(self, s1, s2, vol_tol):
        s1_lengths, s1_angles = s1.lattice.lengths_and_angles
        all_nn = get_points_in_sphere_pbc(
            s2.lattice, [[0, 0, 0]], [0, 0, 0],
            (1 + self.ltol) * max(s1_lengths))[:, [0, 1]]

        nv = []
        for l in s1_lengths:
            nvi = all_nn[np.where((all_nn[:, 1] < (1 + self.ltol) * l)
                                  &
                                  (all_nn[:, 1] > (1 - self.ltol) * l))][:, 0]
            if not len(nvi):
                return
            nvi = [np.array(site) for site in nvi]
            nvi = np.dot(nvi, s2.lattice.matrix)
            nv.append(nvi)
            #The vectors are broadcast into a 5-D array containing
        #all permutations of the entries in nv[0], nv[1], nv[2]
        #produces the same result as three nested loops over the
        #same variables and calculating determinants individually
        bfl = (np.array(nv[0])[None, None, :, None, :] *
               np.array([1, 0, 0])[None, None, None, :, None] +
               np.array(nv[1])[None, :, None, None, :] *
               np.array([0, 1, 0])[None, None, None, :, None] +
               np.array(nv[2])[:, None, None, None, :] *
               np.array([0, 0, 1])[None, None, None, :, None])

        #Compute volume of each array
        vol = np.sum(bfl[:, :, :, 0, :] * np.cross(bfl[:, :, :, 1, :],
                                                   bfl[:, :, :, 2, :]), 3)
        #Find valid lattices
        valid = np.where(abs(vol) >= vol_tol)
        if not len(valid[0]):
            return
            #loop over valid lattices to compute the angles for each
        lengths = np.sum(bfl[valid] ** 2, axis=2) ** 0.5
        angles = np.zeros((len(bfl[valid]), 3), float)
        for i in xrange(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            angles[:, i] = \
                np.sum(bfl[valid][:, j, :] * bfl[valid][:, k, :], 1) \
                / (lengths[:, j] * lengths[:, k])
        angles = np.arccos(angles) * 180. / np.pi
        #Check angles are within tolerance
        valid_angles = np.where(np.all(np.abs(angles - s1_angles) <
                                       self.angle_tol, axis=1))
        if len(valid_angles[0]) == 0:
            return
            #yield valid lattices
        for lat in bfl[valid][valid_angles]:
            nl = Lattice(lat)
            yield nl

    def _cmp_fractional_struct(self, s1, s2, frac_tol):
        #compares the fractional coordinates
        for s1_coords, s2_coords in zip(s1, s2):
            dist = s1_coords[:, None] - s2_coords[None, :]
            dist = abs(dist - np.round(dist))
            dist[np.where(dist > frac_tol[None, None, :])] = 3 * len(dist)
            cost = np.sum(dist, axis=-1)
            if np.max(np.min(cost, axis=0)) >= 3 * len(dist):
                return False
            lin = LinearAssignment(cost)
            if lin.min_cost >= 3 * len(dist):
                return False
        return True

    def _cmp_cartesian_struct(self, s1, s2, l1, l2):
        """
        Once a fit is found, a rms minimizing fit is done to
        ensure the fit is correct. To do this,

        1) The structures are placed into an average lattice
        2) All sites are shifted by the mean
            displacement vector between matched sites.
        3) calculate distances
        4) return rms distance normalized by (V/Natom) ^ 1/3
            and the maximum distance found
        """
        nsites = sum(map(len, s1))

        avg_params = (np.array(l1.lengths_and_angles) +
                      np.array(l2.lengths_and_angles)) / 2

        avg_lattice = Lattice.from_lengths_and_angles(avg_params[0],
                                                      avg_params[1])
        dist = np.zeros([nsites, nsites]) + 100 * nsites
        vec_matrix = np.zeros([nsites, nsites, 3])
        i = 0
        for s1_coords, s2_coords in zip(s1, s2):
            j = len(s1_coords)
            vecs = pbc_shortest_vectors(avg_lattice, s1_coords, s2_coords)
            distances = (np.sum(vecs ** 2, axis=-1)) ** 0.5
            dist[i: i + j, i: i + j] = distances
            vec_matrix[i: i + j, i: i + j] = vecs
            i += j
        lin = LinearAssignment(dist)
        inds = np.arange(nsites)

        shortest_vecs = vec_matrix[inds, lin.solution, :]
        shortest_vec_square = np.sum(
            (shortest_vecs - np.average(shortest_vecs, axis=0)) ** 2, -1)

        norm_length = (avg_lattice.volume / nsites) ** (1 / 3)

        rms = np.average(shortest_vec_square) ** 0.5 / norm_length

        max_dist = np.max(shortest_vec_square) ** 0.5 / norm_length

        return rms, max_dist

    def _supercell_fit(self, struct1, struct2, break_on_match):
        """
        Calculate RMS displacement between two structures
        where one is a potential supercell of the other

        Args:
            struct1:
                1st structure
            struct2:
                2nd structure
            break_on_match:
                True or False. Will break if the maximum
                    distance found is less than the
                    provided stol

        Returns:
            rms displacement normalized by (Vol / nsites) ** (1/3) and
            maximum distance found between two paired sites
        """
        struct1 = Structure.from_sites(struct1)
        struct2 = Structure.from_sites(struct2)

        stored_rms = None

        #reset struct1 to supercell
        if struct2.num_sites > struct1.num_sites:
            [struct2, struct1] = [struct1, struct2]

        #number of formula units
        fu = struct1.num_sites / struct2.num_sites

        struct1 = struct1.get_reduced_structure(reduction_algo="niggli")
        struct2 = struct2.get_reduced_structure(reduction_algo="niggli")

        nl1 = struct1.lattice
        nl2 = struct2.lattice

        #Volume to determine invalid lattices
        vol_tol = fu * nl2.volume / 2

        #fractional tolerance of atomic positions (2x for initial fitting)
        frac_tol = np.array([self.stol / ((1 - self.ltol) * np.pi) * i
                             for i in struct1.lattice.reciprocal_lattice.abc])
        frac_tol *= ((nl1.volume + fu * nl2.volume)
                     / (2 * struct1.num_sites)) ** (1 / 3)

        #generate structure coordinate lists
        species_list = []
        s1 = []
        for site in struct1:
            found = False
            for i, species in enumerate(species_list):
                if self._comparator.are_equal(site.species_and_occu,
                                              species):
                    found = True
                    s1[i].append(site.frac_coords)
                    break
            if not found:
                s1.append([site.frac_coords])
                species_list.append(site.species_and_occu)

        zipped = sorted(zip(s1, species_list), key=lambda x: len(x[0]))

        s1 = [x[0] for x in zipped]
        species_list = [x[1] for x in zipped]

        #translate s1
        s1_translation = s1[0][0]
        for i in range(len(species_list)):
            s1[i] = np.mod(s1[i] - s1_translation, 1)

        #do permutations of vectors, check for equality
        for nl in self._get_lattices(struct1, struct2, vol_tol):

            s2_cart = [[] for i in s1]

            scale_matrix = np.round(np.dot(nl.matrix, nl2.inv_matrix))

            supercell = struct2.copy()
            supercell.make_supercell(scale_matrix.astype('int'))

            for site in supercell:
                found = False
                for i, species in enumerate(species_list):
                    if self._comparator.are_equal(site.species_and_occu,
                                                  species):
                        found = True
                        s2_cart[i].append(site.coords)
                        break
                        #if no site match found return None
                if not found:
                    return None

            #check that sizes of the site groups are identical
            for f1, c2 in zip(s1, s2_cart):
                if len(f1) != len(c2):
                    return None

            s2 = [nl.get_fractional_coords(c) for c in s2_cart]
            for coord in s2[0]:
                t_s2 = [np.mod(coords - coord, 1) for coords in s2]
                if self._cmp_fractional_struct(s1, t_s2, frac_tol):
                    rms, max_dist = self._cmp_cartesian_struct(s1, t_s2, nl,
                                                               nl1)
                    if break_on_match and max_dist < self.stol:
                        return max_dist
                    elif stored_rms is None or rms < stored_rms[0]:
                        stored_rms = rms, max_dist

        if break_on_match:
            return None
        else:
            return stored_rms

    def fit(self, struct1, struct2):
        """
        Fit two structures.

        Args:
            struct1:
                1st structure
            struct2:
                2nd structure

        Returns:
            True or False.
        """

        fit_dist = self._calc_rms(struct1, struct2, break_on_match=True)

        if fit_dist is None:
            return False
        else:
            return fit_dist <= self.stol

    def get_rms_dist(self, struct1, struct2):
        """
        Calculate RMS displacement between two structures

        Args:
            struct1:
                1st structure
            struct2:
                2nd structure

        Returns:
            rms displacement normalized by (Vol / nsites) ** (1/3)
            and maximum distance between paired sites. If no matching
            lattice is found None is returned.
        """

        return self._calc_rms(struct1, struct2, break_on_match=False)

    def _calc_rms(self, struct1, struct2, break_on_match):
        """
        Calculate RMS displacement between two structures

        Args:
            struct1:
                1st structure
            struct2:
                2nd structure
            break_on_match:
                True or False. Will break if the maximum
                    distance found is less than the
                    provided stol

        Returns:
            rms displacement normalized by (Vol / nsites) ** (1/3) and
            maximum distance found between two paired sites
        """
        struct1 = Structure.from_sites(struct1.sites)
        struct2 = Structure.from_sites(struct2.sites)

        stol = self.stol
        comparator = self._comparator
        #initial stored rms
        stored_rms = None

        if comparator.get_structure_hash(struct1) != \
                comparator.get_structure_hash(struct2):
            return None

        #primitive cell transformation
        if self._primitive_cell and struct1.num_sites != struct2.num_sites:
            struct1 = struct1.get_primitive_structure()
            struct2 = struct2.get_primitive_structure()

        # Same number of sites
        if struct1.num_sites != struct2.num_sites:
            #if mismatch try to fit a supercell or return None
            if self._supercell and not (
                    (struct1.num_sites % struct2.num_sites)
                    and (struct2.num_sites % struct1.num_sites)):
                return self._supercell_fit(struct1, struct2, break_on_match)
            else:
                return None

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
            nl1 = Lattice(nl1.matrix * scale_vol)
            struct1.modify_lattice(nl1)
            nl2 = Lattice(nl2.matrix / scale_vol)
            struct2.modify_lattice(nl2)

        #Volume to determine invalid lattices
        vol_tol = nl2.volume / 2

        #fractional tolerance of atomic positions (2x for initial fitting)
        frac_tol = \
            np.array([stol / ((1 - self.ltol) * np.pi) * i for
                      i in struct1.lattice.reciprocal_lattice.abc]) * \
            ((nl1.volume + nl2.volume) /
             (2 * struct1.num_sites)) ** (1.0 / 3)

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
                    #if no site match found return None
            if not found:
                return None

        #check that sizes of the site groups are identical
        for f1, c2 in zip(s1, s2_cart):
            if len(f1) != len(c2):
                return None

        #translate s1
        s1_translation = s1[0][0]
        for i in range(len(species_list)):
            s1[i] = np.mod(s1[i] - s1_translation, 1)
            #do permutations of vectors, check for equality
        for nl in self._get_lattices(struct1, struct2, vol_tol):
            s2 = [nl.get_fractional_coords(c) for c in s2_cart]
            for coord in s2[0]:
                t_s2 = [np.mod(coords - coord, 1) for coords in s2]
                if self._cmp_fractional_struct(s1, t_s2, frac_tol):
                    rms, max_dist = self._cmp_cartesian_struct(s1, t_s2, nl,
                                                               nl1)
                    if break_on_match and max_dist < stol:
                        return max_dist
                    elif stored_rms is None or rms < stored_rms[0]:
                        stored_rms = rms, max_dist

        if break_on_match:
            return None
        else:
            return stored_rms

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
                if s_list[j] in group_list[i]:
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

    @property
    def to_dict(self):
        return {"version": __version__, "@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "comparator": self._comparator.to_dict,
                "stol": self.stol,
                "ltol": self.ltol,
                "angle_tol": self.angle_tol,
                "primitive_cell": self._primitive_cell,
                "scale": self._scale}

    @classmethod
    def from_dict(cls, d):
        return StructureMatcher(
            ltol=d["ltol"], stol=d["stol"], angle_tol=d["angle_tol"],
            primitive_cell=d["primitive_cell"], scale=d["scale"],
            comparator=AbstractComparator.from_dict(d["comparator"]))

    def get_minimax_rms_anonymous(self, struct1, struct2):
        """
        Performs an anonymous fitting, which allows distinct species in one
        structure to map to another. E.g., to compare if the Li2O and Na2O
        structures are similar.

        Args:
            struct1:
                1st structure
            struct2:
                2nd structure

        Returns:
            (minimax_rms, min_mapping)
            min_rms is the minimum of the max rms calculated, and min_mapping
            is the corresponding minimal species mapping that would map
            struct1 to struct2. (None, None) is returned if the minimax_rms
            exceeds the threshold.
        """
        sp1 = list(set(struct1.species_and_occu))
        sp2 = list(set(struct2.species_and_occu))

        if len(sp1) != len(sp2):
            return None, None

        latt1 = struct1.lattice
        fcoords1 = struct1.frac_coords
        min_rms = float("inf")
        min_mapping = None
        for perm in itertools.permutations(sp2):
            sp_mapping = dict(zip(sp1, perm))
            mapped_sp = [sp_mapping[site.species_and_occu] for site in struct1]
            transformed_structure = Structure(latt1, mapped_sp, fcoords1)
            rms = self.get_rms_dist(transformed_structure, struct2)
            if rms is not None:
                if min_rms > rms[1]:
                    min_rms = rms[1]
                    min_mapping = {k: v for k, v in sp_mapping.items()
                                   if k != v}
        if min_mapping is None:
            return None, None
        else:
            return min_rms, min_mapping

    def fit_anonymous_all_mapping(self, struct1, struct2):
        """
        Performs an anonymous fitting, which allows distinct species in one
        structure to map to another. Preferentially pairs together species
        with closer oxidation states over minimizing rms_distance

        Args:
            struct1:
                1st oxidation state decorated structure
            struct2:
                2nd oxidation state decorated structure

        Returns:
            (mappings)
            mappings is a list of possible species mappings that
            would map struct1 to struct2.
        """
        sp1 = list(set(struct1.species_and_occu))
        sp2 = list(set(struct2.species_and_occu))

        if len(sp1) != len(sp2):
            return None

        latt1 = struct1.lattice
        fcoords1 = struct1.frac_coords
        mappings = []
        for perm in itertools.permutations(sp2):
            sp_mapping = dict(zip(sp1, perm))
            mapped_sp = [sp_mapping[site.species_and_occu] for site in struct1]
            transformed_structure = Structure(latt1, mapped_sp, fcoords1)
            rms = self.get_rms_dist(transformed_structure, struct2)
            if rms is not None:

                possible_mapping = {k: v for k, v in sp_mapping.items()
                                    if k != v}

                if rms[1] < self.stol:
                    #check if mapping already found
                    for k, v in possible_mapping.iteritems():
                        if {k: v} not in mappings:
                                mappings.append({k: v})
        if not mappings:
            return None
        else:
            return mappings

    def fit_anonymous(self, struct1, struct2):
        """
        Performs an anonymous fitting, which allows distinct species in one
        structure to map to another. E.g., to compare if the Li2O and Na2O
        structures are similar.

        Args:
            struct1:
                1st structure
            struct2:
                2nd structure

        Returns:
            A minimal species mapping that would map struct1 to struct2 in
            terms of structure similarity, or None if no fit is found. For
            example, to map the cubic Li2O to cubic Na2O,
            we need a Li->Na mapping. This method will return
            [({Element("Li"): 1}, {Element("Na"): 1})]. Since O is the same
            in both structures, there is no O to O mapping required.
            Note that the return form is a list of pairs of species and
            occupancy dicts. This complicated return for is necessary because
            species and occupancy dicts are non-hashable.
        """
        min_rms, min_mapping = self.get_minimax_rms_anonymous(struct1, struct2)

        if min_rms is None or min_rms > self.stol:
            return None
        else:
            return min_mapping