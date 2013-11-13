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


class OrderDisorderElementComparator(AbstractComparator):
    """
    A Comparator that matches sites, given some overlap in the element
    composition
    """

    def are_equal(self, sp1, sp2):
        """
        True if there is some overlap in composition between the species

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
        set1 = set(sp1.element_composition.keys())
        set2 = set(sp2.element_composition.keys())
        if set1.intersection(set2):
            return True
        return False

    def get_structure_hash(self, structure):
        """
        Hash for structure.

        Args:
            structure:
                A structure

        Returns:
            TODO No good hash yet
        """
        return 1


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
                 scale=True, attempt_supercell=False, allow_subset=False,
                 comparator=SpeciesComparator(), supercell_size='num_sites'):
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
            allow_subset:
                Allow one structure to match to the subset of another
                structure. Eg. Matching of an ordered structure onto a
                disordered one, or matching a delithiated to a lithiated
                structure. This option cannot be combined with
                attempt_supercell, or with structure grouping.
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
            supercell_size:
                Method to use for determining the size of a supercell (if
                applicable). Possible values are num_sites, num_atoms or
                volume.
        """
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol
        self._comparator = comparator
        self._primitive_cell = primitive_cell
        self._scale = scale
        self._supercell = attempt_supercell
        self._supercell_size = supercell_size
        self._subset = allow_subset

    def _get_supercell_size(self, target_s, s):
        """
        Returns the unrounded number of units of s in target_s
        """
        if self._supercell_size == 'num_sites':
            fu = target_s.num_sites / s.num_sites
        elif self._supercell_size == 'num_atoms':
            fu = target_s.composition.num_atoms / s.composition.num_atoms
        elif self._supercell_size == 'volume':
            fu = target_s.volume / s.volume
        else:
            raise ValueError('invalid argument for supercell_size')
        return fu

    def _get_lattices(self, target_s, s, supercell_size=1):
        """
        Yields lattices for s with lengths and angles close to the
        lattice of target_s. If supercell_size is specified, the
        returned lattice will have that number of primitive cells
        in it

        Args:
            s, target_s: Structure objects
        """
        t_l, t_a = target_s.lattice.lengths_and_angles
        r = (1 + self.ltol) * max(t_l)
        fpts, dists, i = get_points_in_sphere_pbc(
            lattice=s.lattice, frac_points=[[0, 0, 0]], center=[0, 0, 0],
            r=r).T
        #get possible vectors for a, b, and c
        new_v = []
        for l in t_l:
            max_r = (1 + self.ltol) * l
            min_r = (1 - self.ltol) * l
            vi = fpts[np.where((dists < max_r) & (dists > min_r))]
            if len(vi) == 0:
                return
            cart_vi = np.dot(np.array([i for i in vi]), s.lattice.matrix)
            new_v.append(cart_vi)

        #The vectors are broadcast into a 5-D array containing
        #all permutations of the entries in new_v[0], new_v[1], new_v[2]
        #Produces the same result as three nested loops over the
        #same variables and calculating determinants individually
        bfl = (np.array(new_v[0])[None, None, :, None, :] *
               np.array([1, 0, 0])[None, None, None, :, None] +
               np.array(new_v[1])[None, :, None, None, :] *
               np.array([0, 1, 0])[None, None, None, :, None] +
               np.array(new_v[2])[:, None, None, None, :] *
               np.array([0, 0, 1])[None, None, None, :, None])

        #Compute volume of each lattice
        vol = np.abs(np.sum(bfl[:, :, :, 0, :] *
                            np.cross(bfl[:, :, :, 1, :],
                                     bfl[:, :, :, 2, :]), 3))
        #valid lattices must not change volume
        min_vol = s.volume * 0.999 * supercell_size
        max_vol = s.volume * 1.001 * supercell_size
        bfl = bfl[np.where((vol > min_vol) & (vol < max_vol))]
        if len(bfl) == 0:
            return

        #compute angles
        lengths = np.sum(bfl ** 2, axis=2) ** 0.5
        angles = np.zeros((len(bfl), 3), float)
        for i in xrange(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            angles[:, i] = \
                np.sum(bfl[:, j, :] * bfl[:, k, :], 1) \
                / (lengths[:, j] * lengths[:, k])
        angles = np.arccos(angles) * 180. / np.pi
        #Check angles are within tolerance
        valid_angles = np.where(np.all(np.abs(angles - t_a) <
                                       self.angle_tol, axis=1))
        for lat in bfl[valid_angles]:
            nl = Lattice(lat)
            yield nl

    def _cmp_fractional_struct(self, s1, s2, frac_tol, mask):
        #ensure that we always calculate distances from the subset
        #to the superset
        if len(s1) > len(s2):
            s_superset, s_subset = s1, s2
        else:
            s_superset, s_subset = s2, s1
            mask = mask.T
        #compares the fractional coordinates
        mask_val = 3 * len(s_superset)
        #distance from subset to superset
        dist = s_superset[None, :] - s_subset[:, None]
        dist = abs(dist - np.round(dist))
        dist[np.where(dist > frac_tol[None, None, :])] = mask_val
        cost = np.sum(dist, axis=-1)
        cost[mask] = mask_val
        if np.max(np.min(cost, axis=1)) >= mask_val:
            return False
        if self._subset:
            n = len(s_superset)
            square_cost = np.zeros((n, n))
            square_cost[:cost.shape[0], :cost.shape[1]] = cost
            cost = square_cost
        lin = LinearAssignment(cost)
        if lin.min_cost >= mask_val:
            return False
        return True

    def _cart_dists(self, s1, s2, l1, l2, mask):
        """
        Finds the cartesian distances normalized by (V/Natom) ^ 1/3
        between two structures on the average lattice of l1 and l2
        s_superset and s_subset are lists of fractional coordinates.
        Minimizes the RMS distance of the matching with an additional
        translation (but doesn't change the mapping)
        returns distances, fractional_translation vector
        """
        #ensure that we always calculate distances from the subset
        #to the superset
        if len(s1) > len(s2):
            s_superset, s_subset, mult = s1, s2, 1
        else:
            s_superset, s_subset, mult = s2, s1, -1
            mask = mask.T
        #create the average lattice
        avg_params = (np.array(l1.lengths_and_angles) +
                      np.array(l2.lengths_and_angles)) / 2
        avg_lattice = Lattice.from_lengths_and_angles(*avg_params)
        norm_length = (avg_lattice.volume / len(s_superset)) ** (1 / 3)
        mask_val = 1e20 * norm_length * self.stol

        all_d_2 = np.zeros([len(s_superset), len(s_superset)])
        vec_matrix = np.zeros([len(s_superset), len(s_superset), 3])

        #vectors from subset to superset
        #1st index subset, 2nd index superset
        vecs = pbc_shortest_vectors(avg_lattice, s_subset, s_superset)
        vec_matrix[:len(s_subset), :len(s_superset)] = vecs
        vec_matrix[mask] = mask_val
        d_2 = (np.sum(vecs ** 2, axis=-1))
        all_d_2[:len(s_subset), :len(s_superset)] = d_2
        all_d_2[mask] = mask_val
        lin = LinearAssignment(all_d_2)
        inds = np.arange(len(s_subset))
        #shortest vectors from the subset to the superset
        shortest_vecs = vec_matrix[inds, lin.solution[:len(s_subset)], :]
        translation = np.average(shortest_vecs, axis=0)
        f_translation = avg_lattice.get_fractional_coords(translation)
        shortest_distances = np.sum((shortest_vecs - translation) ** 2,
                                    -1) ** 0.5
        return shortest_distances / norm_length, f_translation * mult

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

        match = self._find_match(struct1, struct2, break_on_match=True)

        if match is None:
            return False
        else:
            return match[0] <= self.stol

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
        match = self._find_match(struct1, struct2, break_on_match=False,
                                 use_rms=True)
        if match is None:
            return None
        else:
            return match[0], max(match[1])

    def _find_match(self, struct1, struct2, break_on_match=False,
                    use_rms=False, niggli=True):
        """
        Finds the best match between two structures.
        Typically, 'best' is determined by minimax cartesian distance
        on the average lattice

        Args:
            struct1:
                1st structure
            struct2:
                2nd structure
            break_on_match:
                If true, breaks once the max distance is below
                the stol (RMS distance if use_rms is true)
            use_rms:
                If True, finds the match that minimizes
                RMS instead of minimax
            niggli:
                whether to compute the niggli cells of the input
                structures

        Returns:
            the value, distances, s2 lattice, and s2 translation
            vector for the best match
        """
        struct1 = Structure.from_sites(struct1.sites)
        struct2 = Structure.from_sites(struct2.sites)

        if (self._comparator.get_structure_hash(struct1) !=
                self._comparator.get_structure_hash(struct2)
                and not self._subset):
            return None

        #primitive cell transformation
        if self._primitive_cell and struct1.num_sites != struct2.num_sites:
            struct1 = struct1.get_primitive_structure()
            struct2 = struct2.get_primitive_structure()

        if self._supercell:
            fu = self._get_supercell_size(struct1, struct2)
            #force struct1 to be the larger one
            if fu < 1:
                struct2, struct1 = struct1, struct2
                fu = 1 / fu
            fu = int(round(fu))
        else:
            fu = 1

        #can't do the check until we group with the comparator
        if (not self._subset) and struct1.num_sites != struct2.num_sites * fu:
            return None

        # Get niggli reduced cells. Though technically not necessary, this
        # minimizes cell lengths and speeds up the matching of skewed
        # cells considerably.
        if niggli:
            struct1 = struct1.get_reduced_structure(reduction_algo="niggli")
            struct2 = struct2.get_reduced_structure(reduction_algo="niggli")

        nl1 = struct1.lattice
        nl2 = struct2.lattice

        #rescale lattice to same volume
        if self._scale:
            ratio = (fu * nl2.volume / nl1.volume) ** (1 / 6)
            nl1 = Lattice(nl1.matrix * ratio)
            struct1.modify_lattice(nl1)
            nl2 = Lattice(nl2.matrix / ratio)
            struct2.modify_lattice(nl2)

        #fractional tolerance of atomic positions (2x for initial fitting)
        normalization = ((2 * struct2.num_sites * fu) /
                         (struct1.volume + struct2.volume * fu)) ** (1 / 3)
        frac_tol = np.array(struct1.lattice.reciprocal_lattice.abc) * \
            self.stol / ((1 - self.ltol) * np.pi) / normalization

        #make array mask
        mask = np.zeros((len(struct2) * fu, len(struct1)), dtype=np.bool)
        i = 0
        for site2 in struct2:
            for repeat in range(fu):
                for j, site1 in enumerate(struct1):
                    mask[i, j] = not self._comparator.are_equal(
                        site2.species_and_occu, site1.species_and_occu)
                i += 1

        #check that there is some valid mapping between sites
        nmax = max(mask.shape)
        sq_mask = np.zeros((nmax, nmax))
        sq_mask[mask] = 10000
        if LinearAssignment(sq_mask).min_cost > 0:
            return None

        #find the best sites for the translation vector
        num_s1_invalid_matches = np.sum(mask, axis=1)
        s2_translation_index = np.argmax(num_s1_invalid_matches)
        s1_translation_indices = np.argwhere(
            mask[s2_translation_index] == 0).flatten()

        s1fc = np.array(struct1.frac_coords)
        s2cc = np.array(struct2.cart_coords)

        best_match = None
        for nl in self._get_lattices(struct1, struct2, fu):
            #if supercell needs to be created, update s2_cart
            if self._supercell and fu > 1:
                scale_matrix = np.round(np.dot(nl.matrix, nl2.inv_matrix))
                supercell = struct2.copy()
                supercell.make_supercell(scale_matrix.astype('int'))
                s2fc = np.array(supercell.frac_coords)
            else:
                s2fc = nl.get_fractional_coords(s2cc)
                #loop over possible translations
            for s1i in s1_translation_indices:
                translation = s1fc[s1i] - s2fc[s2_translation_index]
                t_s2fc = s2fc + translation
                if self._cmp_fractional_struct(s1fc, t_s2fc, frac_tol, mask):
                    distances, t = self._cart_dists(s1fc, t_s2fc, nl, nl1,
                                                    mask)
                    if use_rms:
                        val = np.linalg.norm(distances) / len(distances) ** 0.5
                    else:
                        val = max(distances)
                    if best_match is None or val < best_match[0]:
                        total_translation = translation + t
                        total_translation -= np.round(total_translation)
                        best_match = val, distances, nl, total_translation
                    if break_on_match and val < self.stol:
                        return best_match
        if best_match and best_match[0] < self.stol:
            return best_match

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
        if self._subset:
            raise ValueError("allow_subset cannot be used with"
                             " group_structures")

        #Use structure hash to pre-group structures.
        shash = self._comparator.get_structure_hash
        sorted_s_list = sorted(s_list, key=shash)
        all_groups = []

        #For each pre-grouped list of structures, perform actual matching.
        for k, g in itertools.groupby(sorted_s_list, key=shash):
            unmatched = list(g)
            while len(unmatched) > 0:
                ref = unmatched.pop(0)
                matches = [ref]
                inds = filter(lambda i: self.fit(ref, unmatched[i]),
                              xrange(len(unmatched)))
                matches.extend([unmatched[i] for i in inds])
                unmatched = [unmatched[i] for i in xrange(len(unmatched))
                             if i not in inds]
                all_groups.append(matches)
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

    def fit_with_electronegativity(self, struct1, struct2):
        """
        Performs an anonymous fitting, which allows distinct species in one
        structure to map to another. E.g., to compare if the Li2O and Na2O
        structures are similar. If multiple substitutions are within tolerance
        this will return the one which minimizes the difference in
        electronegativity between the matches species.

        Args:
            struct1:
                1st structure
            struct2:
                2nd structure

        Returns:
            min_mapping
            min_rms is the minimum of the max rms calculated, and min_mapping
            is the corresponding minimal species mapping that would map
            struct1 to struct2. None is returned if no fit is found.
        """

        sp1 = list(set(struct1.species_and_occu))
        sp2 = list(set(struct2.species_and_occu))
        if len(sp1) != len(sp2):
            return None, None

        latt1 = struct1.lattice
        fcoords1 = struct1.frac_coords
        min_X_diff = np.inf
        min_mapping = None
        for perm in itertools.permutations(sp2):
            sp_mapping = dict(zip(sp1, perm))
            mapped_sp = [sp_mapping[site.species_and_occu] for site in struct1]
            transformed_structure = Structure(latt1, mapped_sp, fcoords1)
            if self.fit(transformed_structure, struct2):
                #Calculate electronegativity difference
                X_diff = np.average(
                    [(host_sp.elements[0].X - map_sp.elements[0].X) *
                     struct1.composition.get(host_sp.elements[0]) for
                     host_sp, map_sp in sp_mapping.iteritems()])

                if min_X_diff == 0:
                    return {}

                if min_X_diff > X_diff:
                    min_X_diff = X_diff
                    min_mapping = {k: v for k, v in sp_mapping.items()
                                   if k != v}
        if min_mapping is None:
            return None
        else:
            return min_mapping

    def fit_anonymous_all_mapping(self, struct1, struct2):
        """
        Performs an anonymous fitting, which allows distinct species in one
        structure to map to another. Returns a dictionary of species
        substitutions that are within tolerance

        Args:
            struct1:
                1st structure
            struct2:
                2nd structure

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

                possible_mapping = {k: v for k, v in sp_mapping.items()}

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

    def get_supercell_matrix(self, supercell, struct):
        """
        Returns the supercell matrix for transforming struct to supercell. This
        can be used for very distorted 'supercells' where the primitive cell
        is impossible to find
        """
        if self._primitive_cell:
            raise ValueError("get_supercell_matrix cannot be used with the "
                             "primitive cell option")
        if self._supercell \
                and self._get_supercell_size(supercell, struct) <= 1:
            raise ValueError("The non-supercell must be put onto the basis"
                             " of the supercell, not the other way around")
        match = self._find_match(supercell, struct, break_on_match=False,
                                 use_rms=True, niggli=False)
        if match is None:
            return None

        return np.round(np.dot(match[2].matrix,
                               struct.lattice.inv_matrix)).astype('int')

    def get_s2_like_s1(self, struct1, struct2):
        """
        Performs transformations on struct2 to put it in a basis similar to
        struct1 (without changing any of the inter-site distances)
        """
        if self._primitive_cell:
            raise ValueError("get_s2_like_s1 cannot be used with the primitive"
                             " cell option")
        if self._supercell and self._get_supercell_size(struct1, struct2) <= 1:
            raise ValueError("The non-supercell must be put onto the basis"
                             " of the supercell, not the other way around")
        if self._subset and struct2.num_sites > struct1.num_sites:
            raise ValueError("The smaller structure must be transformed onto"
                             " the larger one, not the other way around")
        match = self._find_match(struct1, struct2, break_on_match=False,
                                 use_rms=True, niggli=False)
        if match is None:
            return None
        scale_matrix = np.round(
            np.dot(match[2].matrix, struct2.lattice.inv_matrix)).astype('int')

        temp = struct2.copy()
        temp.make_supercell(scale_matrix)
        fc = temp.frac_coords + match[3]
        return Structure(temp.lattice, temp.species_and_occu, fc,
                         to_unit_cell=True)
