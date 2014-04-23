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
    pbc_shortest_vectors, lattice_points_in_supercell
from pymatgen.symmetry.finder import SymmetryFinder


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
            sp1: First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2: Second species. A dict of {specie/element: amt} as per the
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
            structure: A structure

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
            sp1: First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2: Second species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.

        Returns:
            Boolean indicating whether species are equal.
        """
        return sp1 == sp2

    def get_structure_hash(self, structure):
        """
        Hash for structure.

        Args:
            structure: A structure

        Returns:
            Reduced formula for the structure is used as a hash for the
            SpeciesComparator.
        """
        return structure.composition.reduced_formula


class SpinComparator(AbstractComparator):
    """
    A Comparator that matches magnetic structures to their inverse spins.
    This comparator is primarily used to filter magnetically ordered
    structures with opposite spins, which are equivalent.
    """

    def are_equal(self, sp1, sp2):
        """
        True if species are exactly the same, i.e., Fe2+ == Fe2+ but not
        Fe3+. and the spins are reversed. i.e., spin up maps to spin down,
        and vice versa.

        Args:
            sp1: First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2: Second species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.

        Returns:
            Boolean indicating whether species are equal.
        """
        for s1 in sp1.keys():
            spin1 = getattr(s1, "spin", 0)
            oxi1 = getattr(s1, "oxi_state", 0)
            found = False
            for s2 in sp2.keys():
                spin2 = getattr(s2, "spin", 0)
                oxi2 = getattr(s2, "oxi_state", 0)
                if (s1.symbol == s2.symbol and oxi1 == oxi2
                        and spin2 == -spin1):
                    found = True
                    break
            if not found:
                return False
        return True

    def get_structure_hash(self, structure):
        """
        Hash for structure.

        Args:
            structure: A structure

        Returns:
            Reduced formula for the structure is used as a hash for the
            SpeciesComparator.
        """
        f = SymmetryFinder(structure, 0.1)
        return "{} {}".format(f.get_spacegroup_number(),
                              structure.composition.reduced_formula)


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
            sp1: First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2: Second species. A dict of {specie/element: amt} as per the
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
            structure: A structure

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
            sp1: First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2: Second species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.

        Returns:
            True always
        """
        return True

    def get_structure_hash(self, structure):
        """
        Hash for structure.

        Args:
            structure: A structure

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
            sp1: First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2: Second species. A dict of {specie/element: amt} as per the
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
            structure: A structure

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

    Args:
        ltol (float): Fractional length tolerance. Default is 0.2.
        stol (float): Site tolerance. Defined as the fraction of the
            average free length per atom := ( V / Nsites ) ** (1/3)
            Default is 0.3.
        angle_tol (float): Angle tolerance in degrees. Default is 5 degrees.
        primitive_cell (bool): If true: input structures will be reduced to
            primitive cells prior to matching. Default to True.
        scale (bool): Input structures are scaled to equivalent volume if
           true; For exact matching, set to False.
        attempt_supercell (bool): If set to True and number of sites in
            cells differ after a primitive cell reduction (divisible by an
            integer) attempts to generate a supercell transformation of the
            smaller cell which is equivalent to the larger structure.
        allow_subset (bool): Allow one structure to match to the subset of
            another structure. Eg. Matching of an ordered structure onto a
            disordered one, or matching a delithiated to a lithiated
            structure. This option cannot be combined with
            attempt_supercell, or with structure grouping.
        comparator (Comparator): A comparator object implementing an equals
            method that declares declaring equivalency of sites. Default is
            SpeciesComparator, which implies rigid species
            mapping, i.e., Fe2+ only matches Fe2+ and not Fe3+.

            Other comparators are provided, e.g., ElementComparator which
            matches only the elements and not the species.

            The reason why a comparator object is used instead of
            supplying a comparison function is that it is not possible to
            pickle a function, which makes it otherwise difficult to use
            StructureMatcher with Python's multiprocessing.
        supercell_size: Method to use for determining the size of a
            supercell (if applicable). Possible values are num_sites,
            num_atoms or volume.
    """

    def __init__(self, ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True,
                 scale=True, attempt_supercell=False, allow_subset=False,
                 comparator=SpeciesComparator(), supercell_size='num_sites'):

        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol
        self._comparator = comparator
        self._primitive_cell = primitive_cell
        self._scale = scale
        self._supercell = attempt_supercell
        self._supercell_size = supercell_size
        self._subset = allow_subset

    def _get_supercell_size(self, s1, s2):
        """
        Returns the supercell size, and whether the supercell should
        be applied to s1. If fu == 1, s1_supercell is returned as
        true, to avoid ambiguity.
        """
        if self._supercell_size == 'num_sites':
            fu = s2.num_sites / s1.num_sites
        elif self._supercell_size == 'num_atoms':
            fu = s2.composition.num_atoms / s1.composition.num_atoms
        elif self._supercell_size == 'volume':
            fu = s2.volume / s1.volume
        else:
            raise ValueError('invalid argument for supercell_size')

        if fu < 2/3:
            return round(1/fu), False
        else:
            return round(fu), True

    def _get_lattices(self, target_lattice, s, supercell_size=1):
        """
        Yields lattices for s with lengths and angles close to the
        lattice of target_s. If supercell_size is specified, the
        returned lattice will have that number of primitive cells
        in it

        Args:
            s, target_s: Structure objects
        """
        t_l, t_a = target_lattice.lengths_and_angles
        r = (1 + self.ltol) * max(t_l)
        fpts, dists, i = s.lattice.get_points_in_sphere(
            frac_points=[[0, 0, 0]], center=[0, 0, 0],
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

    def _get_supercells(self, struct1, struct2, fu, s1_supercell):
        """
        Computes all supercells of one structure close to the lattice of the
        other
        if s1_supercell == True, it makes the supercells of struct1, otherwise
        it makes them of s2

        yields: s1, s2, supercell_matrix, average_lattice, supercell_matrix
        """
        def av_lat(l1, l2):
            params = (np.array(l1.lengths_and_angles) + \
                      np.array(l2.lengths_and_angles)) / 2
            return Lattice.from_lengths_and_angles(*params)

        def sc_generator(s1, s2):
            s2_fc = np.array(s2.frac_coords)
            if fu == 1:
                cc = np.array(s1.cart_coords)
                for l in self._get_lattices(s2.lattice, s1, fu):
                    supercell_matrix = np.round(np.dot(l.matrix,
                        s1.lattice.inv_matrix)).astype('int')
                    fc = l.get_fractional_coords(cc)
                    fc -= np.floor(fc)
                    yield fc, s2_fc, av_lat(l, s2.lattice), supercell_matrix
            else:
                fc_init = np.array(s1.frac_coords)
                for l in self._get_lattices(s2.lattice, s1, fu):
                    supercell_matrix = np.round(np.dot(l.matrix,
                        s1.lattice.inv_matrix)).astype('int')
                    fc = np.dot(fc_init, np.linalg.inv(supercell_matrix))
                    lp = lattice_points_in_supercell(supercell_matrix)
                    fc = (fc[:, None, :] + lp[None, :, :]).reshape((-1, 3))
                    fc -= np.floor(fc)
                    yield fc, s2_fc, av_lat(l, s2.lattice), supercell_matrix
        if s1_supercell:
            for x in sc_generator(struct1, struct2):
                yield x
        else:
            for x in sc_generator(struct2, struct1):
                #reorder generator output so s1 is still first
                yield x[1], x[0], x[2], x[3]

    def _cmp_fstruct(self, s1, s2, frac_tol, mask):
        """
        Returns true if a matching exists between s2 and s2
        under frac_tol. s2 should be a subset of s1
        """
        if len(s2) > len(s1):
            raise ValueError("s1 must be larger than s2")
        if mask.shape != (len(s2), len(s1)):
            raise ValueError("mask has incorrect shape")

        mask_val = 3 * len(s1)
        #distance from subset to superset
        dist = s1[None, :] - s2[:, None]
        dist = abs(dist - np.round(dist))

        dist[np.where(dist > frac_tol[None, None, :])] = mask_val
        cost = np.sum(dist, axis=-1)
        cost[mask] = mask_val

        #maximin is a lower bound on linear assignment
        #(and faster to compute)
        if np.max(np.min(cost, axis=1)) >= mask_val:
            return False

        return LinearAssignment(cost).min_cost < mask_val

    def _cart_dists(self, s1, s2, avg_lattice, mask):
        """
        Finds a matching in cartesian space. Finds an additional
        fractional translation vector to minimize RMS distance

        Args:
            s1, s2: numpy arrays of fractional coordinates.
                len(s1) >= len(s2)
            avg_lattice: Lattice on which to calculate distances
            mask: numpy array of booleans. mask[i, j] = True indicates
                that s2[i] cannot be matched to s1[j]

        Returns:
            Distances from s2 to s1, normalized by (V/Natom) ^ 1/3
            Fractional translation vector to apply to s2.
            Mapping from s1 to s2, i.e. with numpy slicing, s1[mapping] => s2
        """
        if len(s2) > len(s1):
            raise ValueError("s1 must be larger than s2")
        if mask.shape != (len(s2), len(s1)):
            raise ValueError("mask has incorrect shape")

        norm_length = (avg_lattice.volume / len(s1)) ** (1 / 3)
        mask_val = 1e10 * norm_length * self.stol
        #vectors are from s2 to s1
        vecs = pbc_shortest_vectors(avg_lattice, s2, s1)
        vecs[mask] = mask_val
        d_2 = np.sum(vecs ** 2, axis=-1)
        lin = LinearAssignment(d_2)
        s = lin.solution
        short_vecs = vecs[np.arange(len(s)), s]
        translation = np.average(short_vecs, axis=0)
        f_translation = avg_lattice.get_fractional_coords(translation)
        new_d2 = np.sum((short_vecs - translation) ** 2, axis=-1)

        return new_d2 ** 0.5 / norm_length, f_translation, s

    def _get_mask(self, struct1, struct2, fu, s1_supercell):
        """
        Returns mask for matching struct2 to struct1. If struct1 has sites
        a b c, and fu = 2, assumes supercells of struct2 will be ordered
        aabbcc (rather than abcabc)

        Returns:
        mask, struct1 translation indices, struct2 translation index
        """
        mask = np.zeros((len(struct2), len(struct1), fu), dtype=np.bool)
        for i, site2 in enumerate(struct2):
            for j, site1 in enumerate(struct1):
                mask[i, j, :] = not self._comparator.are_equal(
                    site2.species_and_occu, site1.species_and_occu)
        if s1_supercell:
            mask = mask.reshape((len(struct2), -1))
        else:
            #supercell is of struct2, roll fu axis back to preserve
            #correct ordering
            mask = np.rollaxis(mask, 2, 1)
            mask = mask.reshape((-1, len(struct1)))

        #find the best translation index
        i = np.argmax(np.sum(mask, axis=-1))
        return mask, np.where(np.invert(mask[i]))[0], i

    def fit(self, struct1, struct2):
        """
        Fit two structures.

        Args:
            struct1 (Structure): 1st structure
            struct2 (Structure): 2nd structure

        Returns:
            True or False.
        """
        if not self._subset and self._comparator.get_structure_hash(struct1) \
                != self._comparator.get_structure_hash(struct2):
            return None

        struct1, struct2, fu, s1_supercell = self._preprocess(struct1, struct2)
        ratio = fu if s1_supercell else 1/fu

        if len(struct1) * ratio >= len(struct2):
            match = self._match(struct1, struct2, fu, s1_supercell=s1_supercell,
                                break_on_match=True)
        else:
            match = self._match(struct2, struct1, fu,
                                s1_supercell=(not s1_supercell),
                                break_on_match=True)

        if match is None:
            return False
        else:
            return match[0] <= self.stol

    def get_rms_dist(self, struct1, struct2):
        """
        Calculate RMS displacement between two structures

        Args:
            struct1 (Structure): 1st structure
            struct2 (Structure): 2nd structure

        Returns:
            rms displacement normalized by (Vol / nsites) ** (1/3)
            and maximum distance between paired sites. If no matching
            lattice is found None is returned.
        """
        struct1, struct2, fu, s1_supercell = self._preprocess(struct1, struct2)
        ratio = fu if s1_supercell else 1/fu

        if len(struct1) * ratio >= len(struct2):
            match = self._match(struct1, struct2, fu, s1_supercell=s1_supercell,
                                break_on_match=False, use_rms=True)
        else:
            match = self._match(struct2, struct1, fu,
                                s1_supercell=(not s1_supercell),
                                break_on_match=False, use_rms=True)
        if match is None:
            return None
        else:
            return match[0], max(match[1])

    def _preprocess(self, struct1, struct2, niggli=True):
        """
        Rescales, finds the reduced structures (primitive and niggli),
        and finds fu, the supercell size to make struct1 comparable to
        s2
        """
        struct1 = struct1.copy()
        struct2 = struct2.copy()

        if niggli:
            struct1 = struct1.get_reduced_structure(reduction_algo="niggli")
            struct2 = struct2.get_reduced_structure(reduction_algo="niggli")

        #primitive cell transformation
        if self._primitive_cell:
            struct1 = struct1.get_primitive_structure()
            struct2 = struct2.get_primitive_structure()

        if self._supercell:
            fu, s1_supercell = self._get_supercell_size(struct1, struct2)
        else:
            fu, s1_supercell = 1, True
        mult = fu if s1_supercell else 1/fu

        #rescale lattice to same volume
        if self._scale:
            ratio = (struct2.volume / (struct1.volume * mult)) ** (1 / 6)
            nl1 = Lattice(struct1.lattice.matrix * ratio)
            struct1.modify_lattice(nl1)
            nl2 = Lattice(struct2.lattice.matrix / ratio)
            struct2.modify_lattice(nl2)

        return struct1, struct2, fu, s1_supercell

    def _match(self, struct1, struct2, fu, s1_supercell=True, use_rms=False,
                   break_on_match=False):
        """
        Matches struct2 onto struct1 (which should contain all sites in
        struct2).
        Args:
            struct1, struct2 (Structure): structures to be matched
            fu (int): size of supercell to create
            s1_supercell (bool): whether to create the supercell of
                struct1 (vs struct2)
            use_rms (bool): whether to minimize the rms of the matching
            break_on_match (bool): whether to stop search at first
                valid match
        """
        if fu < 1:
            raise ValueError("fu cannot be less than 1")

        mask, s1_t_inds, s2_t_ind = self._get_mask(struct1, struct2,
                                                   fu, s1_supercell)

        if mask.shape[0] > mask.shape[1]:
            raise ValueError('after supercell creation, struct1 must '
                             'have more sites than struct2')

        #check that a valid mapping exists
        if not self._subset and mask.shape[1] != mask.shape[0]:
            return None
        if LinearAssignment(mask).min_cost > 0:
            return None

        best_match = None
        #loop over all lattices
        for s1fc, s2fc, avg_l, sc_m in \
                self._get_supercells(struct1, struct2, fu, s1_supercell):
            #compute fractional tolerance
            normalization = (len(s1fc) / avg_l.volume) ** (1/3)
            inv_abc = np.array(avg_l.reciprocal_lattice.abc)
            frac_tol = inv_abc * self.stol / (np.pi * normalization)
            #loop over all translations
            for s1i in s1_t_inds:
                t = s1fc[s1i] - s2fc[s2_t_ind]
                t_s2fc = s2fc + t
                if self._cmp_fstruct(s1fc, t_s2fc, frac_tol, mask):
                    dist, t_adj, mapping = self._cart_dists(s1fc, t_s2fc,
                                                        avg_l, mask)
                    if use_rms:
                        val = np.linalg.norm(dist) / len(dist) ** 0.5
                    else:
                        val = max(dist)
                    if best_match is None or val < best_match[0]:
                        total_t = t + t_adj
                        total_t -= np.round(total_t)
                        best_match = val, dist, sc_m, total_t, mapping
                        if (break_on_match or val < 1e-5) and val < self.stol:
                            return best_match
        if best_match and best_match[0] < self.stol:
            return best_match

    def group_structures(self, s_list):
        """
        Given a list of structures, use fit to group
        them by structural equality.

        Args:
            s_list ([Structure]): List of structures to be grouped

        Returns:
            A list of lists of matched structures
            Assumption: if s1 == s2 but s1 != s3, than s2 and s3 will be put
            in different groups without comparison.
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
            struct1 (Structure): 1st structure
            struct2 (Structure): 2nd structure

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
            struct1 (Structure): 1st structure
            struct2 (Structure): 2nd structure

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
            struct1 (Structure): 1st structure
            struct2 (Structure): 2nd structure

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
            possible_mapping = {k: v for k, v in sp_mapping.items()}
            if self.fit(transformed_structure, struct2):
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
            struct1 (Structure): 1st structure
            struct2 (Structure): 2nd structure

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
        Returns the matrix for transforming struct to supercell. This
        can be used for very distorted 'supercells' where the primitive cell
        is impossible to find
        """
        if self._primitive_cell:
            raise ValueError("get_supercell_matrix cannot be used with the "
                             "primitive cell option")
        struct, supercell, fu, s1_supercell = self._preprocess(struct,
                                                               supercell, False)
        ratio = fu if s1_supercell else 1/fu

        if not s1_supercell:
            raise ValueError("The non-supercell must be put onto the basis"
                             " of the supercell, not the other way around")

        if len(supercell) >= len(struct) * ratio:
            match = self._match(supercell, struct, fu, s1_supercell=False,
                                use_rms=True, break_on_match=False)
        else:
            match = self._match(struct, supercell, fu, s1_supercell=True,
                                use_rms=True, break_on_match=False)

        if match is None:
            return None

        return match[2]

    def get_s2_like_s1(self, struct1, struct2):
        """
        Performs transformations on struct2 to put it in a basis similar to
        struct1 (without changing any of the inter-site distances)
        """
        if self._primitive_cell:
            raise ValueError("get_s2_like_s1 cannot be used with the primitive"
                             " cell option")

        s1, s2, fu, s1_supercell = self._preprocess(struct1, struct2, False)
        ratio = fu if s1_supercell else 1/fu
        if s1_supercell and fu > 1:
            raise ValueError("Struct1 must be the supercell, "
                             "not the other way around")

        if len(s1) * ratio >= len(s2):
            #s1 is superset
            match = self._match(s1, s2, fu=fu, s1_supercell=False,
                                use_rms=True, break_on_match=False)
        else:
            #s2 is superset
            match = self._match(s2, s1, fu=fu, s1_supercell=True,
                                use_rms=True, break_on_match=False)

        if match is None:
            return None

        temp = struct2.copy()
        temp.make_supercell(match[2])

        if len(struct1) >= len(temp):
            #invert the mapping, since it needs to be from s2 to s1
            mapping = np.argsort(match[4])
            tvec = match[3]
        else:
            #add sites not included in the mapping
            not_included = range(len(temp))
            for i in match[4]:
                not_included.remove(i)
            mapping = list(match[4]) + not_included
            tvec = -match[3]

        temp.translate_sites(range(len(temp)), tvec)
        return Structure.from_sites([temp.sites[i] for i in mapping])

    def get_mapping(self, superset, subset):
        """
        Calculate the mapping from superset to subset, i.e.
        superset[mapping] = subset
        """
        if self._supercell:
            raise ValueError("cannot compute mapping to supercell")
        if self._primitive_cell:
            raise ValueError("cannot compute mapping with primitive cell "
                             "option")
        if len(subset) > len(superset):
            raise ValueError("subset is larger than superset")

        superset, subset, _, _ = self._preprocess(superset, subset, True)
        match = self._match(superset, subset, 1, break_on_match=False)

        if match[0] > self.stol:
            return None
        return match[4]
