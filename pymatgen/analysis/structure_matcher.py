# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to perform fitting of structures.
"""
import abc
import itertools

import numpy as np
from monty.json import MSONable

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.structure import Structure
from pymatgen.optimization.linear_assignment import LinearAssignment
from pymatgen.util.coord import lattice_points_in_supercell
from pymatgen.util.coord_cython import is_coord_subset_pbc, pbc_shortest_vectors

__author__ = "William Davidson Richards, Stephen Dacek, Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "William Davidson Richards"
__email__ = "wrichard@mit.edu"
__status__ = "Production"
__date__ = "Dec 3, 2012"


class AbstractComparator(MSONable, metaclass=abc.ABCMeta):
    """
    Abstract Comparator class. A Comparator defines how sites are compared in
    a structure.
    """

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
    def get_hash(self, composition):
        """
        Defines a hash to group structures. This allows structures to be
        grouped efficiently for comparison. The hash must be invariant under
        supercell creation. (e.g. composition is not a good hash, but
        fractional_composition might be). Reduced formula is not a good formula,
        due to weird behavior with fractional occupancy.

        Composition is used here instead of structure because for anonymous
        matches it is much quicker to apply a substitution to a composition
        object than a structure object.

        Args:
            composition (Composition): composition of the structure

        Returns:
            A hashable object. Examples can be string formulas, integers etc.
        """
        return

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: Comparator.
        """
        for trans_modules in ["structure_matcher"]:
            mod = __import__(
                "pymatgen.analysis." + trans_modules,
                globals(),
                locals(),
                [d["@class"]],
                0,
            )
            if hasattr(mod, d["@class"]):
                trans = getattr(mod, d["@class"])
                return trans()
        raise ValueError("Invalid Comparator dict")

    def as_dict(self):
        """
        :return: MSONable dict
        """
        return {
            "version": __version__,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }


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

    def get_hash(self, composition):
        """
        Returns: Fractional composition
        """
        return composition.fractional_composition


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
        for s1 in sp1:
            spin1 = getattr(s1, "spin", 0)
            oxi1 = getattr(s1, "oxi_state", 0)
            for s2 in sp2:
                spin2 = getattr(s2, "spin", 0)
                oxi2 = getattr(s2, "oxi_state", 0)
                if s1.symbol == s2.symbol and oxi1 == oxi2 and spin2 == -spin1:
                    break
            else:
                return False
        return True

    def get_hash(self, composition):
        """
        Returns: Fractional composition
        """
        return composition.fractional_composition


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

    def get_hash(self, composition):
        """
        Returns: Fractional element composition
        """
        return composition.element_composition.fractional_composition


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

    def get_hash(self, composition):
        """
        No hash possible
        """
        return 1


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
        set1 = set(sp1.elements)
        set2 = set(sp2.elements)
        return set1.issubset(set2) or set2.issubset(set1)

    def get_hash(self, composition):
        """
        Returns: Fractional composition
        """
        return composition.fractional_composition


class OccupancyComparator(AbstractComparator):
    """
    A Comparator that matches occupancies on sites,
    irrespective of the species of those sites.
    """

    def are_equal(self, sp1, sp2):
        """
        Args:
            sp1: First species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.
            sp2: Second species. A dict of {specie/element: amt} as per the
                definition in Site and PeriodicSite.

        Returns:
            True if sets of occupancies (amt) are equal on both sites.
        """
        return set(sp1.element_composition.values()) == set(sp2.element_composition.values())

    def get_hash(self, composition):
        """
        :param composition: Composition.
        :return: 1. Difficult to define sensible hash
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

                ia. Convert both lattices to Cartesian and place
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

    def __init__(
        self,
        ltol=0.2,
        stol=0.3,
        angle_tol=5,
        primitive_cell=True,
        scale=True,
        attempt_supercell=False,
        allow_subset=False,
        comparator=None,
        supercell_size="num_sites",
        ignored_species=None,
    ):
        """
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
            supercell_size (str or list): Method to use for determining the
                size of a supercell (if applicable). Possible values are
                num_sites, num_atoms, volume, or an element or list of elements
                present in both structures.
            ignored_species (list): A list of ions to be ignored in matching.
                Useful for matching structures that have similar frameworks
                except for certain ions, e.g., Li-ion intercalation frameworks.
                This is more useful than allow_subset because it allows better
                control over what species are ignored in the matching.
        """

        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol
        self._comparator = comparator or SpeciesComparator()
        self._primitive_cell = primitive_cell
        self._scale = scale
        self._supercell = attempt_supercell
        self._supercell_size = supercell_size
        self._subset = allow_subset
        self._ignored_species = [] if ignored_species is None else ignored_species[:]

    def _get_supercell_size(self, s1, s2):
        """
        Returns the supercell size, and whether the supercell should
        be applied to s1. If fu == 1, s1_supercell is returned as
        true, to avoid ambiguity.
        """
        if self._supercell_size == "num_sites":
            fu = s2.num_sites / s1.num_sites
        elif self._supercell_size == "num_atoms":
            fu = s2.composition.num_atoms / s1.composition.num_atoms
        elif self._supercell_size == "volume":
            fu = s2.volume / s1.volume
        elif not isinstance(self._supercell_size, str):
            s1comp, s2comp = 0, 0
            for el in self._supercell_size:
                el = get_el_sp(el)
                s1comp += s1.composition[el]
                s2comp += s2.composition[el]
            fu = s2comp / s1comp
        else:
            el = get_el_sp(self._supercell_size)
            if (el in s2.composition) and (el in s1.composition):
                fu = s2.composition[el] / s1.composition[el]
            else:
                raise ValueError("Invalid argument for supercell_size.")

        if fu < 2 / 3:
            return int(round(1 / fu)), False

        return int(round(fu)), True

    def _get_lattices(self, target_lattice, s, supercell_size=1):
        """
        Yields lattices for s with lengths and angles close to the lattice of target_s. If
        supercell_size is specified, the returned lattice will have that number of primitive
        cells in it

        Args:
            s, target_s: Structure objects
        """
        lattices = s.lattice.find_all_mappings(
            target_lattice,
            ltol=self.ltol,
            atol=self.angle_tol,
            skip_rotation_matrix=True,
        )
        for l, _, scale_m in lattices:
            if abs(abs(np.linalg.det(scale_m)) - supercell_size) < 0.5:
                yield l, scale_m

    def _get_supercells(self, struct1, struct2, fu, s1_supercell):
        """
        Computes all supercells of one structure close to the lattice of the
        other
        if s1_supercell is True, it makes the supercells of struct1, otherwise
        it makes them of s2

        yields: s1, s2, supercell_matrix, average_lattice, supercell_matrix
        """

        def av_lat(l1, l2):
            params = (np.array(l1.parameters) + np.array(l2.parameters)) / 2
            return Lattice.from_parameters(*params)

        def sc_generator(s1, s2):
            s2_fc = np.array(s2.frac_coords)
            if fu == 1:
                cc = np.array(s1.cart_coords)
                for l, sc_m in self._get_lattices(s2.lattice, s1, fu):
                    fc = l.get_fractional_coords(cc)
                    fc -= np.floor(fc)
                    yield fc, s2_fc, av_lat(l, s2.lattice), sc_m
            else:
                fc_init = np.array(s1.frac_coords)
                for l, sc_m in self._get_lattices(s2.lattice, s1, fu):
                    fc = np.dot(fc_init, np.linalg.inv(sc_m))
                    lp = lattice_points_in_supercell(sc_m)
                    fc = (fc[:, None, :] + lp[None, :, :]).reshape((-1, 3))
                    fc -= np.floor(fc)
                    yield fc, s2_fc, av_lat(l, s2.lattice), sc_m

        if s1_supercell:
            for x in sc_generator(struct1, struct2):
                yield x
        else:
            for x in sc_generator(struct2, struct1):
                # reorder generator output so s1 is still first
                yield x[1], x[0], x[2], x[3]

    @classmethod
    def _cmp_fstruct(cls, s1, s2, frac_tol, mask):
        """
        Returns true if a matching exists between s2 and s2
        under frac_tol. s2 should be a subset of s1
        """
        if len(s2) > len(s1):
            raise ValueError("s1 must be larger than s2")
        if mask.shape != (len(s2), len(s1)):
            raise ValueError("mask has incorrect shape")

        return is_coord_subset_pbc(s2, s1, frac_tol, mask)

    @classmethod
    def _cart_dists(cls, s1, s2, avg_lattice, mask, normalization, lll_frac_tol=None):
        """
        Finds a matching in Cartesian space. Finds an additional
        fractional translation vector to minimize RMS distance

        Args:
            s1, s2: numpy arrays of fractional coordinates. len(s1) >= len(s2)
            avg_lattice: Lattice on which to calculate distances
            mask: numpy array of booleans. mask[i, j] = True indicates
                that s2[i] cannot be matched to s1[j]
            normalization (float): inverse normalization length

        Returns:
            Distances from s2 to s1, normalized by (V/Natom) ^ 1/3
            Fractional translation vector to apply to s2.
            Mapping from s1 to s2, i.e. with numpy slicing, s1[mapping] => s2
        """
        if len(s2) > len(s1):
            raise ValueError("s1 must be larger than s2")
        if mask.shape != (len(s2), len(s1)):
            raise ValueError("mask has incorrect shape")

        # vectors are from s2 to s1
        vecs, d_2 = pbc_shortest_vectors(avg_lattice, s2, s1, mask, return_d2=True, lll_frac_tol=lll_frac_tol)
        lin = LinearAssignment(d_2)
        s = lin.solution  # pylint: disable=E1101
        short_vecs = vecs[np.arange(len(s)), s]
        translation = np.average(short_vecs, axis=0)
        f_translation = avg_lattice.get_fractional_coords(translation)
        new_d2 = np.sum((short_vecs - translation) ** 2, axis=-1)

        return new_d2**0.5 * normalization, f_translation, s

    def _get_mask(self, struct1, struct2, fu, s1_supercell):
        """
        Returns mask for matching struct2 to struct1. If struct1 has sites
        a b c, and fu = 2, assumes supercells of struct2 will be ordered
        aabbcc (rather than abcabc)

        Returns:
        mask, struct1 translation indices, struct2 translation index
        """
        mask = np.zeros((len(struct2), len(struct1), fu), dtype=bool)

        inner = []
        for sp2, i in itertools.groupby(enumerate(struct2.species_and_occu), key=lambda x: x[1]):
            i = list(i)
            inner.append((sp2, slice(i[0][0], i[-1][0] + 1)))

        for sp1, j in itertools.groupby(enumerate(struct1.species_and_occu), key=lambda x: x[1]):
            j = list(j)
            j = slice(j[0][0], j[-1][0] + 1)
            for sp2, i in inner:
                mask[i, j, :] = not self._comparator.are_equal(sp1, sp2)

        if s1_supercell:
            mask = mask.reshape((len(struct2), -1))
        else:
            # supercell is of struct2, roll fu axis back to preserve
            # correct ordering
            mask = np.rollaxis(mask, 2, 1)
            mask = mask.reshape((-1, len(struct1)))

        # find the best translation indices
        i = np.argmax(np.sum(mask, axis=-1))
        inds = np.where(np.invert(mask[i]))[0]
        if s1_supercell:
            # remove the symmetrically equivalent s1 indices
            inds = inds[::fu]
        return np.array(mask, dtype=int), inds, i

    def fit(
        self, struct1: Structure, struct2: Structure, symmetric: bool = False, skip_structure_reduction: bool = False
    ) -> bool:
        """
        Fit two structures.

        Args:
            struct1 (Structure): 1st structure
            struct2 (Structure): 2nd structure
            symmetric (Bool): Defaults to False
                If True, check the equality both ways.
                This only impacts a small percentage of structures
            skip_structure_reduction (Bool): Defaults to False
                If True, skip to get a primitive structure and perform Niggli reduction for struct1 and struct2

        Returns:
            True or False.
        """
        struct1, struct2 = self._process_species([struct1, struct2])

        if not self._subset and self._comparator.get_hash(struct1.composition) != self._comparator.get_hash(
            struct2.composition
        ):
            return False

        if not symmetric:
            struct1, struct2, fu, s1_supercell = self._preprocess(
                struct1, struct2, skip_structure_reduction=skip_structure_reduction
            )
            match = self._match(struct1, struct2, fu, s1_supercell, break_on_match=True)
            if match is None:
                return False

            return match[0] <= self.stol

        struct1, struct2, fu, s1_supercell = self._preprocess(
            struct1, struct2, skip_structure_reduction=skip_structure_reduction
        )
        match1 = self._match(struct1, struct2, fu, s1_supercell, break_on_match=True)
        struct1, struct2 = struct2, struct1
        struct1, struct2, fu, s1_supercell = self._preprocess(
            struct1, struct2, skip_structure_reduction=skip_structure_reduction
        )
        match2 = self._match(struct1, struct2, fu, s1_supercell, break_on_match=True)

        if match1 is None or match2 is None:
            return False

        return max(match1[0], match2[0]) <= self.stol

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
        struct1, struct2 = self._process_species([struct1, struct2])
        struct1, struct2, fu, s1_supercell = self._preprocess(struct1, struct2)
        match = self._match(struct1, struct2, fu, s1_supercell, use_rms=True, break_on_match=False)

        if match is None:
            return None

        return match[0], max(match[1])

    def _process_species(self, structures):
        copied_structures = []
        for s in structures:
            # We need the copies to be actual Structure to work properly, not
            # subclasses. So do type(s) == Structure.
            ss = Structure.from_sites(s)
            if self._ignored_species:
                ss.remove_species(self._ignored_species)
            copied_structures.append(ss)
        return copied_structures

    def _preprocess(self, struct1, struct2, niggli=True, skip_structure_reduction: bool = False):
        """
        Rescales, finds the reduced structures (primitive and niggli),
        and finds fu, the supercell size to make struct1 comparable to
        s2.
        If skip_structure_reduction is True, skip to get reduced structures (by primitive transformation and
        niggli reduction). This option is useful for fitting a set of structures several times.
        """
        if skip_structure_reduction:
            # Need to copy original structures to rescale lattices later
            struct1 = struct1.copy()
            struct2 = struct2.copy()
        else:
            struct1 = self._get_reduced_structure(struct1, self._primitive_cell, niggli)
            struct2 = self._get_reduced_structure(struct2, self._primitive_cell, niggli)

        if self._supercell:
            fu, s1_supercell = self._get_supercell_size(struct1, struct2)
        else:
            fu, s1_supercell = 1, True
        mult = fu if s1_supercell else 1 / fu

        # rescale lattice to same volume
        if self._scale:
            ratio = (struct2.volume / (struct1.volume * mult)) ** (1 / 6)
            nl1 = Lattice(struct1.lattice.matrix * ratio)
            struct1.lattice = nl1
            nl2 = Lattice(struct2.lattice.matrix / ratio)
            struct2.lattice = nl2

        return struct1, struct2, fu, s1_supercell

    def _match(
        self,
        struct1,
        struct2,
        fu,
        s1_supercell=True,
        use_rms=False,
        break_on_match=False,
    ):
        """
        Matches one struct onto the other
        """
        ratio = fu if s1_supercell else 1 / fu
        if len(struct1) * ratio >= len(struct2):
            return self._strict_match(
                struct1,
                struct2,
                fu,
                s1_supercell=s1_supercell,
                break_on_match=break_on_match,
                use_rms=use_rms,
            )
        return self._strict_match(
            struct2,
            struct1,
            fu,
            s1_supercell=(not s1_supercell),
            break_on_match=break_on_match,
            use_rms=use_rms,
        )

    def _strict_match(
        self,
        struct1,
        struct2,
        fu,
        s1_supercell=True,
        use_rms=False,
        break_on_match=False,
    ):
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

        mask, s1_t_inds, s2_t_ind = self._get_mask(struct1, struct2, fu, s1_supercell)

        if mask.shape[0] > mask.shape[1]:
            raise ValueError("after supercell creation, struct1 must have more sites than struct2")

        # check that a valid mapping exists
        if (not self._subset) and mask.shape[1] != mask.shape[0]:
            return None

        if LinearAssignment(mask).min_cost > 0:  # pylint: disable=E1101
            return None

        best_match = None
        # loop over all lattices
        for s1fc, s2fc, avg_l, sc_m in self._get_supercells(struct1, struct2, fu, s1_supercell):
            # compute fractional tolerance
            normalization = (len(s1fc) / avg_l.volume) ** (1 / 3)
            inv_abc = np.array(avg_l.reciprocal_lattice.abc)
            frac_tol = inv_abc * self.stol / (np.pi * normalization)
            # loop over all translations
            for s1i in s1_t_inds:
                t = s1fc[s1i] - s2fc[s2_t_ind]
                t_s2fc = s2fc + t
                if self._cmp_fstruct(s1fc, t_s2fc, frac_tol, mask):
                    inv_lll_abc = np.array(avg_l.get_lll_reduced_lattice().reciprocal_lattice.abc)
                    lll_frac_tol = inv_lll_abc * self.stol / (np.pi * normalization)
                    dist, t_adj, mapping = self._cart_dists(s1fc, t_s2fc, avg_l, mask, normalization, lll_frac_tol)
                    if use_rms:
                        val = np.linalg.norm(dist) / len(dist) ** 0.5
                    else:
                        val = max(dist)
                    # pylint: disable=E1136
                    if best_match is None or val < best_match[0]:
                        total_t = t + t_adj
                        total_t -= np.round(total_t)
                        best_match = val, dist, sc_m, total_t, mapping
                        if (break_on_match or val < 1e-5) and val < self.stol:
                            return best_match

        if best_match and best_match[0] < self.stol:
            return best_match

        return None

    def group_structures(self, s_list, anonymous=False):
        """
        Given a list of structures, use fit to group
        them by structural equality.

        Args:
            s_list ([Structure]): List of structures to be grouped
            anonymous (bool): Whether to use anonymous mode.

        Returns:
            A list of lists of matched structures
            Assumption: if s1 == s2 but s1 != s3, than s2 and s3 will be put
            in different groups without comparison.
        """
        if self._subset:
            raise ValueError("allow_subset cannot be used with group_structures")

        original_s_list = list(s_list)
        s_list = self._process_species(s_list)
        # Prepare reduced structures beforehand
        s_list = [self._get_reduced_structure(s, self._primitive_cell, niggli=True) for s in s_list]

        # Use structure hash to pre-group structures
        if anonymous:

            def c_hash(c):
                return c.anonymized_formula

        else:
            c_hash = self._comparator.get_hash

        def s_hash(s):
            return c_hash(s[1].composition)

        sorted_s_list = sorted(enumerate(s_list), key=s_hash)
        all_groups = []

        # For each pre-grouped list of structures, perform actual matching.
        for _, g in itertools.groupby(sorted_s_list, key=s_hash):
            unmatched = list(g)
            while len(unmatched) > 0:
                i, refs = unmatched.pop(0)
                matches = [i]
                if anonymous:
                    inds = filter(
                        lambda i: self.fit_anonymous(refs, unmatched[i][1], skip_structure_reduction=True),
                        list(range(len(unmatched))),
                    )
                else:
                    inds = filter(
                        lambda i: self.fit(refs, unmatched[i][1], skip_structure_reduction=True),
                        list(range(len(unmatched))),
                    )
                inds = list(inds)
                matches.extend([unmatched[i][0] for i in inds])
                unmatched = [unmatched[i] for i in range(len(unmatched)) if i not in inds]
                all_groups.append([original_s_list[i] for i in matches])

        return all_groups

    def as_dict(self):
        """
        :return: MSONable dict
        """
        return {
            "version": __version__,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "comparator": self._comparator.as_dict(),
            "stol": self.stol,
            "ltol": self.ltol,
            "angle_tol": self.angle_tol,
            "primitive_cell": self._primitive_cell,
            "scale": self._scale,
            "attempt_supercell": self._supercell,
            "allow_subset": self._subset,
            "supercell_size": self._supercell_size,
            "ignored_species": self._ignored_species,
        }

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: StructureMatcher
        """
        return cls(
            ltol=d["ltol"],
            stol=d["stol"],
            angle_tol=d["angle_tol"],
            primitive_cell=d["primitive_cell"],
            scale=d["scale"],
            attempt_supercell=d["attempt_supercell"],
            allow_subset=d["allow_subset"],
            comparator=AbstractComparator.from_dict(d["comparator"]),
            supercell_size=d["supercell_size"],
            ignored_species=d["ignored_species"],
        )

    def _anonymous_match(
        self,
        struct1,
        struct2,
        fu,
        s1_supercell=True,
        use_rms=False,
        break_on_match=False,
        single_match=False,
    ):
        """
        Tries all permutations of matching struct1 to struct2.
        Args:
            struct1, struct2 (Structure): Preprocessed input structures
        Returns:
            List of (mapping, match)
        """
        if not isinstance(self._comparator, SpeciesComparator):
            raise ValueError("Anonymous fitting currently requires SpeciesComparator")

        # check that species lists are comparable
        sp1 = struct1.composition.elements
        sp2 = struct2.composition.elements
        if len(sp1) != len(sp2):
            return None

        ratio = fu if s1_supercell else 1 / fu
        swapped = len(struct1) * ratio < len(struct2)

        s1_comp = struct1.composition
        s2_comp = struct2.composition
        matches = []
        for perm in itertools.permutations(sp2):
            sp_mapping = dict(zip(sp1, perm))

            # do quick check that compositions are compatible
            mapped_comp = Composition({sp_mapping[k]: v for k, v in s1_comp.items()})
            if (not self._subset) and (self._comparator.get_hash(mapped_comp) != self._comparator.get_hash(s2_comp)):
                continue

            mapped_struct = struct1.copy()
            mapped_struct.replace_species(sp_mapping)
            if swapped:
                m = self._strict_match(
                    struct2,
                    mapped_struct,
                    fu,
                    (not s1_supercell),
                    use_rms,
                    break_on_match,
                )
            else:
                m = self._strict_match(mapped_struct, struct2, fu, s1_supercell, use_rms, break_on_match)
            if m:
                matches.append((sp_mapping, m))
                if single_match:
                    break
        return matches

    @classmethod
    def _get_reduced_structure(cls, struct: Structure, primitive_cell: bool = True, niggli: bool = True) -> Structure:
        """
        Helper method to find a reduced structure
        """
        reduced = struct.copy()
        if niggli:
            reduced = reduced.get_reduced_structure(reduction_algo="niggli")
        if primitive_cell:
            reduced = reduced.get_primitive_structure()
        return reduced

    def get_rms_anonymous(self, struct1, struct2):
        """
        Performs an anonymous fitting, which allows distinct species in one
        structure to map to another. E.g., to compare if the Li2O and Na2O
        structures are similar.

        Args:
            struct1 (Structure): 1st structure
            struct2 (Structure): 2nd structure

        Returns:
            (min_rms, min_mapping)
            min_rms is the minimum rms distance, and min_mapping is the
            corresponding minimal species mapping that would map
            struct1 to struct2. (None, None) is returned if the minimax_rms
            exceeds the threshold.
        """
        struct1, struct2 = self._process_species([struct1, struct2])
        struct1, struct2, fu, s1_supercell = self._preprocess(struct1, struct2)

        matches = self._anonymous_match(struct1, struct2, fu, s1_supercell, use_rms=True, break_on_match=False)
        if matches:
            best = sorted(matches, key=lambda x: x[1][0])[0]
            return best[1][0], best[0]

        return None, None

    def get_best_electronegativity_anonymous_mapping(self, struct1, struct2):
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
            min_mapping (dict): Mapping of struct1 species to struct2 species
        """
        struct1, struct2 = self._process_species([struct1, struct2])
        struct1, struct2, fu, s1_supercell = self._preprocess(struct1, struct2)

        matches = self._anonymous_match(struct1, struct2, fu, s1_supercell, use_rms=True, break_on_match=True)

        if matches:
            min_X_diff = np.inf
            for m in matches:
                X_diff = 0
                for k, v in m[0].items():
                    X_diff += struct1.composition[k] * (k.X - v.X) ** 2
                if X_diff < min_X_diff:
                    min_X_diff = X_diff
                    best = m[0]
            return best

        return None

    def get_all_anonymous_mappings(self, struct1, struct2, niggli=True, include_dist=False):
        """
        Performs an anonymous fitting, which allows distinct species in one
        structure to map to another. Returns a dictionary of species
        substitutions that are within tolerance

        Args:
            struct1 (Structure): 1st structure
            struct2 (Structure): 2nd structure
            niggli (bool): Find niggli cell in preprocessing
            include_dist (bool): Return the maximin distance with each mapping

        Returns:
            list of species mappings that map struct1 to struct2.
        """
        struct1, struct2 = self._process_species([struct1, struct2])
        struct1, struct2, fu, s1_supercell = self._preprocess(struct1, struct2, niggli)

        matches = self._anonymous_match(struct1, struct2, fu, s1_supercell, break_on_match=not include_dist)
        if matches:
            if include_dist:
                return [(m[0], m[1][0]) for m in matches]

            return [m[0] for m in matches]

        return None

    def fit_anonymous(
        self, struct1: Structure, struct2: Structure, niggli: bool = True, skip_structure_reduction: bool = False
    ):
        """
        Performs an anonymous fitting, which allows distinct species in one
        structure to map to another. E.g., to compare if the Li2O and Na2O
        structures are similar.

        Args:
            struct1 (Structure): 1st structure
            struct2 (Structure): 2nd structure
            niggli (Bool): If true, perform Niggli reduction for struct1 and struct2
            skip_structure_reduction (Bool): Defaults to False
                If True, skip to get a primitive structure and perform Niggli reduction for struct1 and struct2

        Returns:
            True/False: Whether a species mapping can map struct1 to stuct2
        """
        struct1, struct2 = self._process_species([struct1, struct2])
        struct1, struct2, fu, s1_supercell = self._preprocess(struct1, struct2, niggli, skip_structure_reduction)

        matches = self._anonymous_match(struct1, struct2, fu, s1_supercell, break_on_match=True, single_match=True)

        return bool(matches)

    def get_supercell_matrix(self, supercell, struct):
        """
        Returns the matrix for transforming struct to supercell. This
        can be used for very distorted 'supercells' where the primitive cell
        is impossible to find
        """
        if self._primitive_cell:
            raise ValueError("get_supercell_matrix cannot be used with the primitive cell option")
        struct, supercell, fu, s1_supercell = self._preprocess(struct, supercell, niggli=False)

        if not s1_supercell:
            raise ValueError("The non-supercell must be put onto the basis of the supercell, not the other way around")

        match = self._match(struct, supercell, fu, s1_supercell, use_rms=True, break_on_match=False)

        if match is None:
            return None

        return match[2]

    def get_transformation(self, struct1, struct2):
        """
        Returns the supercell transformation, fractional translation vector,
        and a mapping to transform struct2 to be similar to struct1.

        Args:
            struct1 (Structure): Reference structure
            struct2 (Structure): Structure to transform.

        Returns:
            supercell (numpy.ndarray(3, 3)): supercell matrix
            vector (numpy.ndarray(3)): fractional translation vector
            mapping (list(int or None)):
                The first len(struct1) items of the mapping vector are the
                indices of struct1's corresponding sites in struct2 (or None
                if there is no corresponding site), and the other items are
                the remaining site indices of struct2.
        """
        if self._primitive_cell:
            raise ValueError("get_transformation cannot be used with the primitive cell option")

        struct1, struct2 = self._process_species((struct1, struct2))

        s1, s2, fu, s1_supercell = self._preprocess(struct1, struct2, niggli=False)
        ratio = fu if s1_supercell else 1 / fu
        if s1_supercell and fu > 1:
            raise ValueError("Struct1 must be the supercell, not the other way around")

        if len(s1) * ratio >= len(s2):
            # s1 is superset
            match = self._strict_match(s1, s2, fu=fu, s1_supercell=False, use_rms=True, break_on_match=False)
            if match is None:
                return None
            # invert the mapping, since it needs to be from s1 to s2
            mapping = [list(match[4]).index(i) if i in match[4] else None for i in range(len(s1))]
            return match[2], match[3], mapping
        # s2 is superset
        match = self._strict_match(s2, s1, fu=fu, s1_supercell=True, use_rms=True, break_on_match=False)
        if match is None:
            return None
        # add sites not included in the mapping
        not_included = list(range(len(s2) * fu))
        for i in match[4]:
            not_included.remove(i)
        mapping = list(match[4]) + not_included
        return match[2], -match[3], mapping

    def get_s2_like_s1(self, struct1, struct2, include_ignored_species=True):
        """
        Performs transformations on struct2 to put it in a basis similar to
        struct1 (without changing any of the inter-site distances)

        Args:
            struct1 (Structure): Reference structure
            struct2 (Structure): Structure to transform.
            include_ignored_species (bool): Defaults to True,
                the ignored_species is also transformed to the struct1
                lattice orientation, though obviously there is no direct
                matching to existing sites.

        Returns:
            A structure object similar to struct1, obtained by making a
            supercell, sorting, and translating struct2.
        """
        s1, s2 = self._process_species([struct1, struct2])
        trans = self.get_transformation(s1, s2)
        if trans is None:
            return None
        sc, t, mapping = trans
        sites = list(s2)
        # Append the ignored sites at the end.
        sites.extend([site for site in struct2 if site not in s2])
        temp = Structure.from_sites(sites)

        temp.make_supercell(sc)
        temp.translate_sites(list(range(len(temp))), t)
        # translate sites to correct unit cell
        for i, j in enumerate(mapping[: len(s1)]):
            if j is not None:
                vec = np.round(struct1[i].frac_coords - temp[j].frac_coords)
                temp.translate_sites(j, vec, to_unit_cell=False)

        sites = [temp.sites[i] for i in mapping if i is not None]

        if include_ignored_species:
            start = int(round(len(temp) / len(struct2) * len(s2)))
            sites.extend(temp.sites[start:])

        return Structure.from_sites(sites)

    def get_mapping(self, superset, subset):
        """
        Calculate the mapping from superset to subset.

        Args:
            superset (Structure): Structure containing at least the sites in
                subset (within the structure matching tolerance)
            subset (Structure): Structure containing some of the sites in
                superset (within the structure matching tolerance)

        Returns:
            numpy array such that superset.sites[mapping] is within matching
            tolerance of subset.sites or None if no such mapping is possible
        """
        if self._supercell:
            raise ValueError("cannot compute mapping to supercell")
        if self._primitive_cell:
            raise ValueError("cannot compute mapping with primitive cell option")
        if len(subset) > len(superset):
            raise ValueError("subset is larger than superset")

        superset, subset, _, _ = self._preprocess(superset, subset, niggli=True)
        match = self._strict_match(superset, subset, 1, break_on_match=False)

        if match is None or match[0] > self.stol:
            return None

        return match[4]
