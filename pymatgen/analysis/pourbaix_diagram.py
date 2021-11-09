# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module is intended to be used to compute Pourbaix diagrams
of arbitrary compositions and formation energies.  If you use
this module in your work, please consider citing the following:

General formalism for solid-aqueous equilibria from DFT:
    Persson et al., DOI: 10.1103/PhysRevB.85.235438
Decomposition maps, or Pourbaix hull diagrams
    Singh et al., DOI: 10.1021/acs.chemmater.7b03980
Fast computation of many-element Pourbaix diagrams:
    Patel et al., https://arxiv.org/abs/1909.00035 (submitted)
"""

import itertools
import logging
import re
import warnings
from copy import deepcopy
from functools import cmp_to_key, lru_cache, partial
from multiprocessing import Pool
from typing import Optional, Union, List, Dict

import numpy as np
from monty.json import MontyDecoder, MSONable
from scipy.spatial import ConvexHull, HalfspaceIntersection

try:
    from scipy.special import comb
except ImportError:
    from scipy.misc import comb

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError
from pymatgen.core.composition import Composition
from pymatgen.core.ion import Ion
from pymatgen.core.periodic_table import Element
from pymatgen.entries.compatibility import MU_H2O
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.coord import Simplex
from pymatgen.util.plotting import pretty_plot
from pymatgen.util.string import Stringify
from pymatgen.util.sequence import PBar


__author__ = "Sai Jayaraman"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.4"
__maintainer__ = "Joseph Montoya"
__credits__ = "Arunima Singh, Joseph Montoya, Anjli Patel"
__email__ = "joseph.montoya@tri.global"
__status__ = "Production"
__date__ = "Nov 1, 2012"

logger = logging.getLogger(__name__)

PREFAC = 0.0591


# TODO: Revise to more closely reflect PDEntry, invoke from energy/composition
# TODO: PourbaixEntries depend implicitly on having entry energies be
#       formation energies, should be a better way to get from raw energies
# TODO: uncorrected_energy is a bit of a misnomer, but not sure what to rename
class PourbaixEntry(MSONable, Stringify):
    """
    An object encompassing all data relevant to a solid or ion
    in a pourbaix diagram.  Each bulk solid/ion has an energy
    g of the form: e = e0 + 0.0591 log10(conc) - nO mu_H2O
    + (nH - 2nO) pH + phi (-nH + 2nO + q)

    Note that the energies corresponding to the input entries
    should be formation energies with respect to hydrogen and
    oxygen gas in order for the pourbaix diagram formalism to
    work. This may be changed to be more flexible in the future.
    """

    def __init__(self, entry, entry_id=None, concentration=1e-6):
        """
        Args:
            entry (ComputedEntry/ComputedStructureEntry/PDEntry/IonEntry): An
                entry object
            entry_id ():
            concentration ():
        """
        self.entry = entry
        if isinstance(entry, IonEntry):
            self.concentration = concentration
            self.phase_type = "Ion"
            self.charge = entry.ion.charge
        else:
            self.concentration = 1.0
            self.phase_type = "Solid"
            self.charge = 0.0
        self.uncorrected_energy = entry.energy
        if entry_id is not None:
            self.entry_id = entry_id
        elif hasattr(entry, "entry_id") and entry.entry_id:
            self.entry_id = entry.entry_id
        else:
            self.entry_id = None

    @property
    def npH(self):
        """
        Returns:
        """
        return self.entry.composition.get("H", 0.0) - 2 * self.entry.composition.get("O", 0.0)

    @property
    def nH2O(self):
        """
        Returns: Number of H2O.
        """
        return self.entry.composition.get("O", 0.0)

    @property
    def nPhi(self):
        """
        Returns: Number of H2O.
        """
        return self.npH - self.charge

    @property
    def name(self):
        """
        Returns: Name for entry
        """
        if self.phase_type == "Solid":
            return self.entry.composition.reduced_formula + "(s)"

        return self.entry.name

    @property
    def energy(self):
        """
        returns energy

        Returns (float): total energy of the pourbaix
            entry (at pH, V = 0 vs. SHE)
        """
        # Note: this implicitly depends on formation energies as input
        return self.uncorrected_energy + self.conc_term - (MU_H2O * self.nH2O)

    @property
    def energy_per_atom(self):
        """
        energy per atom of the pourbaix entry

        Returns (float): energy per atom
        """
        return self.energy / self.composition.num_atoms

    def energy_at_conditions(self, pH, V):
        """
        Get free energy for a given pH and V

        Args:
            pH (float): pH at which to evaluate free energy
            V (float): voltage at which to evaluate free energy

        Returns:
            free energy at conditions
        """
        return self.energy + self.npH * PREFAC * pH + self.nPhi * V

    def get_element_fraction(self, element):
        """
        Gets the elemental fraction of a given non-OH element

        Args:
            element (Element or str): string or element corresponding
                to element to get from composition

        Returns:
            fraction of element / sum(all non-OH elements)

        """
        return self.composition.get(element) * self.normalization_factor

    @property
    def normalized_energy(self):
        """
        Returns:
             energy normalized by number of non H or O atoms, e. g.
             for Zn2O6, energy / 2 or for AgTe3(OH)3, energy / 4
        """
        return self.energy * self.normalization_factor

    def normalized_energy_at_conditions(self, pH, V):
        """
        Energy at an electrochemical condition, compatible with
        numpy arrays for pH/V input

        Args:
            pH (float): pH at condition
            V (float): applied potential at condition

        Returns:
            energy normalized by number of non-O/H atoms at condition
        """
        return self.energy_at_conditions(pH, V) * self.normalization_factor

    @property
    def conc_term(self):
        """
        Returns the concentration contribution to the free energy,
        and should only be present when there are ions in the entry
        """
        return PREFAC * np.log10(self.concentration)

    # TODO: not sure if these are strictly necessary with refactor
    def as_dict(self):
        """
        Returns dict which contains Pourbaix Entry data.
        Note that the pH, voltage, H2O factors are always calculated when
        constructing a PourbaixEntry object.
        """
        d = {"@module": self.__class__.__module__, "@class": self.__class__.__name__}
        if isinstance(self.entry, IonEntry):
            d["entry_type"] = "Ion"
        else:
            d["entry_type"] = "Solid"
        d["entry"] = self.entry.as_dict()
        d["concentration"] = self.concentration
        d["entry_id"] = self.entry_id
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Invokes a PourbaixEntry from a dictionary
        """
        entry_type = d["entry_type"]
        if entry_type == "Ion":
            entry = IonEntry.from_dict(d["entry"])
        else:
            entry = MontyDecoder().process_decoded(d["entry"])
        entry_id = d["entry_id"]
        concentration = d["concentration"]
        return PourbaixEntry(entry, entry_id, concentration)

    @property
    def normalization_factor(self):
        """
        Sum of number of atoms minus the number of H and O in composition
        """
        return 1.0 / (self.num_atoms - self.composition.get("H", 0) - self.composition.get("O", 0))

    @property
    def composition(self):
        """
        Returns composition
        """
        return self.entry.composition

    @property
    def num_atoms(self):
        """
        Return number of atoms in current formula. Useful for normalization
        """
        return self.composition.num_atoms

    def to_pretty_string(self) -> str:
        """
        :return: A pretty string representation.
        """
        if self.phase_type == "Solid":
            return self.entry.composition.reduced_formula + "(s)"

        return self.entry.name

    def __repr__(self):
        return "Pourbaix Entry : {} with energy = {:.4f}, npH = {}, nPhi = {}, nH2O = {}, entry_id = {} ".format(
            self.entry.composition,
            self.energy,
            self.npH,
            self.nPhi,
            self.nH2O,
            self.entry_id,
        )

    def __str__(self):
        return self.__repr__()


class MultiEntry(PourbaixEntry):
    """
    PourbaixEntry-like object for constructing multi-elemental Pourbaix
    diagrams.
    """

    def __init__(self, entry_list, weights=None):
        """
        Initializes a MultiEntry.

        Args:
            entry_list ([PourbaixEntry]): List of component PourbaixEntries
            weights ([float]): Weights associated with each entry. Default is None
        """
        if weights is None:
            self.weights = [1.0] * len(entry_list)
        else:
            self.weights = weights
        self.entry_list = entry_list

    @lru_cache()
    def __getattr__(self, item):
        """
        Because most of the attributes here are just weighted
        averages of the entry_list, we save some space by
        having a set of conditionals to define the attributes
        """
        # Attributes that are weighted averages of entry attributes
        if item in [
            "energy",
            "npH",
            "nH2O",
            "nPhi",
            "conc_term",
            "composition",
            "uncorrected_energy",
        ]:
            # TODO: Composition could be changed for compat with sum
            if item == "composition":
                start = Composition({})
            else:
                start = 0
            return sum(
                [getattr(e, item) * w for e, w in zip(self.entry_list, self.weights)],
                start,
            )

        # Attributes that are just lists of entry attributes
        if item in ["entry_id", "phase_type"]:
            return [getattr(e, item) for e in self.entry_list]

        # normalization_factor, num_atoms should work from superclass
        return self.__getattribute__(item)

    @property
    def name(self):
        """
        MultiEntry name, i. e. the name of each entry joined by ' + '
        """
        return " + ".join([e.name for e in self.entry_list])

    def __repr__(self):
        return (
            "Multiple Pourbaix Entry: energy = {:.4f}, npH = {}, nPhi = {}, "
            "nH2O = {}, entry_id = {}, species: {}".format(
                self.energy, self.npH, self.nPhi, self.nH2O, self.entry_id, self.name
            )
        )

    def __str__(self):
        return self.__repr__()

    def as_dict(self):
        """
        Returns: MSONable dict
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "entry_list": [e.as_dict() for e in self.entry_list],
            "weights": self.weights,
        }

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (): Dict representation

        Returns:
            MultiEntry
        """
        entry_list = [PourbaixEntry.from_dict(e) for e in d.get("entry_list")]
        return cls(entry_list, d.get("weights"))


# TODO: this class isn't particularly useful in its current form, could be
#       refactored to include information about the reference solid
class IonEntry(PDEntry):
    """
    Object similar to PDEntry, but contains an Ion object instead of a
    Composition object.

    .. attribute:: name

        A name for the entry. This is the string shown in the phase diagrams.
        By default, this is the reduced formula for the composition, but can be
        set to some other string for display purposes.
    """

    def __init__(self, ion, energy, name=None, attribute=None):
        """
        Args:
            ion: Ion object
            energy: Energy for composition.
            name: Optional parameter to name the entry. Defaults to the
                chemical formula.
        """
        self.ion = ion
        # Auto-assign name
        name = name if name else self.ion.reduced_formula
        super().__init__(composition=ion.composition, energy=energy, name=name, attribute=attribute)

    @classmethod
    def from_dict(cls, d):
        """
        Returns an IonEntry object from a dict.
        """
        return IonEntry(Ion.from_dict(d["ion"]), d["energy"], d.get("name"), d.get("attribute"))

    def as_dict(self):
        """
        Creates a dict of composition, energy, and ion name
        """
        d = {"ion": self.ion.as_dict(), "energy": self.energy, "name": self.name}
        return d

    def __repr__(self):
        return "IonEntry : {} with energy = {:.4f}".format(self.composition, self.energy)

    def __str__(self):
        return self.__repr__()


def ion_or_solid_comp_object(formula):
    """
    Returns either an ion object or composition object given
    a formula.

    Args:
        formula: String formula. Eg. of ion: NaOH(aq), Na[+];
            Eg. of solid: Fe2O3(s), Fe(s), Na2O

    Returns:
        Composition/Ion object
    """
    m = re.search(r"\[([^\[\]]+)\]|\(aq\)", formula)
    if m:
        comp_obj = Ion.from_formula(formula)
    elif re.search(r"\(s\)", formula):
        comp_obj = Composition(formula[:-3])
    else:
        comp_obj = Composition(formula)
    return comp_obj


ELEMENTS_HO = {Element("H"), Element("O")}


# TODO: the solids filter breaks some of the functionality of the
#       heatmap plotter, because the reference states for decomposition
#       don't include oxygen/hydrogen in the OER/HER regions

# TODO: create a from_phase_diagram class method for non-formation energy
#       invocation
# TODO: invocation from a MultiEntry entry list could be a bit more robust
# TODO: serialization is still a bit rough around the edges
class PourbaixDiagram(MSONable):
    """
    Class to create a Pourbaix diagram from entries
    """

    def __init__(
        self,
        entries: Union[List[PourbaixEntry], List[MultiEntry]],
        comp_dict: Optional[Dict[str, float]] = None,
        conc_dict: Optional[Dict[str, float]] = None,
        filter_solids: bool = True,
        nproc: Optional[int] = None,
    ):
        """
        Args:
            entries ([PourbaixEntry] or [MultiEntry]): Entries list
                containing Solids and Ions or a list of MultiEntries
            comp_dict ({str: float}): Dictionary of compositions,
                defaults to equal parts of each elements
            conc_dict ({str: float}): Dictionary of ion concentrations,
                defaults to 1e-6 for each element
            filter_solids (bool): applying this filter to a Pourbaix
                diagram ensures all included solid phases are filtered by
                stability on the compositional phase diagram. Defaults to True.
                The practical consequence of this is that highly oxidized or reduced
                phases that might show up in experiments due to kinetic limitations
                on oxygen/hydrogen evolution won't appear in the diagram, but they are
                not actually "stable" (and are frequently overstabilized from DFT errors).
                Hence, including only the stable solid phases generally leads to the
                most accurate Pourbaix diagrams.
            nproc (int): number of processes to generate multientries with
                in parallel.  Defaults to None (serial processing)
        """
        entries = deepcopy(entries)
        self.filter_solids = filter_solids

        # Get non-OH elements
        self.pbx_elts = list(
            set(itertools.chain.from_iterable([entry.composition.elements for entry in entries])) - ELEMENTS_HO
        )
        self.dim = len(self.pbx_elts) - 1

        # Process multientry inputs
        if isinstance(entries[0], MultiEntry):
            self._processed_entries = entries
            # Extract individual entries
            single_entries = list(set(itertools.chain.from_iterable([e.entry_list for e in entries])))
            self._unprocessed_entries = single_entries
            self._filtered_entries = single_entries
            self._conc_dict = None
            self._elt_comp = {k: v for k, v in entries[0].composition.items() if k not in ELEMENTS_HO}
            self._multielement = True

        # Process single entry inputs
        else:
            # Set default conc/comp dicts
            if not comp_dict:
                comp_dict = {elt.symbol: 1.0 / len(self.pbx_elts) for elt in self.pbx_elts}
            if not conc_dict:
                conc_dict = {elt.symbol: 1e-6 for elt in self.pbx_elts}
            self._conc_dict = conc_dict

            self._elt_comp = comp_dict
            self.pourbaix_elements = self.pbx_elts

            solid_entries = [entry for entry in entries if entry.phase_type == "Solid"]
            ion_entries = [entry for entry in entries if entry.phase_type == "Ion"]

            # If a conc_dict is specified, override individual entry concentrations
            for entry in ion_entries:
                ion_elts = list(set(entry.composition.elements) - ELEMENTS_HO)
                # TODO: the logic here for ion concentration setting is in two
                #       places, in PourbaixEntry and here, should be consolidated
                if len(ion_elts) == 1:
                    entry.concentration = conc_dict[ion_elts[0].symbol] * entry.normalization_factor
                elif len(ion_elts) > 1 and not entry.concentration:
                    raise ValueError("Elemental concentration not compatible with multi-element ions")

            self._unprocessed_entries = solid_entries + ion_entries

            if not len(solid_entries + ion_entries) == len(entries):
                raise ValueError("All supplied entries must have a phase type of " 'either "Solid" or "Ion"')

            if self.filter_solids:
                # O is 2.46 b/c pbx entry finds energies referenced to H2O
                entries_HO = [ComputedEntry("H", 0), ComputedEntry("O", 2.46)]
                solid_pd = PhaseDiagram(solid_entries + entries_HO)
                solid_entries = list(set(solid_pd.stable_entries) - set(entries_HO))

            self._filtered_entries = solid_entries + ion_entries
            if len(comp_dict) > 1:
                self._multielement = True
                self._processed_entries = self._preprocess_pourbaix_entries(self._filtered_entries, nproc=nproc)
            else:
                self._processed_entries = self._filtered_entries
                self._multielement = False

        self._stable_domains, self._stable_domain_vertices = self.get_pourbaix_domains(self._processed_entries)

    def _convert_entries_to_points(self, pourbaix_entries):
        """
        Args:
            pourbaix_entries ([PourbaixEntry]): list of pourbaix entries
                to process into vectors in nph-nphi-composition space

        Returns:
            list of vectors, [[nph, nphi, e0, x1, x2, ..., xn-1]]
            corresponding to each entry in nph-nphi-composition space

        """
        vecs = [
            [entry.npH, entry.nPhi, entry.energy] + [entry.composition.get(elt) for elt in self.pbx_elts[:-1]]
            for entry in pourbaix_entries
        ]
        vecs = np.array(vecs)
        norms = np.transpose([[entry.normalization_factor for entry in pourbaix_entries]])
        vecs *= norms
        return vecs

    def _get_hull_in_nph_nphi_space(self, entries):
        """
        Generates convex hull of pourbaix diagram entries in composition,
        npH, and nphi space.  This enables filtering of multi-entries
        such that only compositionally stable combinations of entries
        are included.

        Args:
            entries ([PourbaixEntry]): list of PourbaixEntries to construct
                the convex hull

        Returns: list of entries and stable facets corresponding to that
            list of entries

        """
        ion_entries = [entry for entry in entries if entry.phase_type == "Ion"]
        solid_entries = [entry for entry in entries if entry.phase_type == "Solid"]

        # Pre-filter solids based on min at each composition
        logger.debug("Pre-filtering solids by min energy at each composition")
        sorted_entries = sorted(
            solid_entries,
            key=lambda x: (x.composition.reduced_composition, x.entry.energy_per_atom),
        )
        grouped_by_composition = itertools.groupby(sorted_entries, key=lambda x: x.composition.reduced_composition)
        min_entries = [list(grouped_entries)[0] for comp, grouped_entries in grouped_by_composition]
        min_entries += ion_entries

        logger.debug("Constructing nph-nphi-composition points for qhull")

        vecs = self._convert_entries_to_points(min_entries)
        maxes = np.max(vecs[:, :3], axis=0)
        extra_point = np.concatenate([maxes, np.ones(self.dim) / self.dim], axis=0)

        # Add padding for extra point
        pad = 1000
        extra_point[2] += pad
        points = np.concatenate([vecs, np.array([extra_point])], axis=0)
        logger.debug("Constructing convex hull in nph-nphi-composition space")
        hull = ConvexHull(points, qhull_options="QJ i")

        # Create facets and remove top
        facets = [facet for facet in hull.simplices if not len(points) - 1 in facet]

        if self.dim > 1:
            logger.debug("Filtering facets by pourbaix composition")
            valid_facets = []
            for facet in facets:
                comps = vecs[facet][:, 3:]
                full_comps = np.concatenate([comps, 1 - np.sum(comps, axis=1).reshape(len(comps), 1)], axis=1)
                # Ensure an compositional interior point exists in the simplex
                if np.linalg.matrix_rank(full_comps) > self.dim:
                    valid_facets.append(facet)
        else:
            valid_facets = facets

        return min_entries, valid_facets

    def _preprocess_pourbaix_entries(self, entries, nproc=None):
        """
        Generates multi-entries for pourbaix diagram

        Args:
            entries ([PourbaixEntry]): list of PourbaixEntries to preprocess
                into MultiEntries
            nproc (int): number of processes to be used in parallel
                treatment of entry combos

        Returns:
            ([MultiEntry]) list of stable MultiEntry candidates

        """
        # Get composition
        tot_comp = Composition(self._elt_comp)

        min_entries, valid_facets = self._get_hull_in_nph_nphi_space(entries)

        combos = []
        for facet in valid_facets:
            for i in range(1, self.dim + 2):
                these_combos = []
                for combo in itertools.combinations(facet, i):
                    these_entries = [min_entries[i] for i in combo]
                    these_combos.append(frozenset(these_entries))
                combos.append(these_combos)

        all_combos = set(itertools.chain.from_iterable(combos))

        list_combos = []
        for i in all_combos:
            list_combos.append(list(i))
        all_combos = list_combos

        multi_entries = []

        # Parallel processing of multi-entry generation
        if nproc is not None:
            f = partial(self.process_multientry, prod_comp=tot_comp)
            with Pool(nproc) as p:
                multi_entries = list(PBar(p.imap(f, all_combos), total=len(all_combos)))
            multi_entries = list(filter(bool, multi_entries))
        else:
            # Serial processing of multi-entry generation
            for combo in PBar(all_combos):
                multi_entry = self.process_multientry(combo, prod_comp=tot_comp)
                if multi_entry:
                    multi_entries.append(multi_entry)

        return multi_entries

    def _generate_multielement_entries(self, entries, nproc=None):
        """
        Create entries for multi-element Pourbaix construction.

        This works by finding all possible linear combinations
        of entries that can result in the specified composition
        from the initialized comp_dict.

        Args:
            entries ([PourbaixEntries]): list of pourbaix entries
                to process into MultiEntries
            nproc (int): number of processes to be used in parallel
                treatment of entry combos
        """

        N = len(self._elt_comp)  # No. of elements
        total_comp = Composition(self._elt_comp)

        # generate all combinations of compounds that have all elements
        entry_combos = [itertools.combinations(entries, j + 1) for j in range(N)]
        entry_combos = itertools.chain.from_iterable(entry_combos)

        entry_combos = filter(lambda x: total_comp < MultiEntry(x).composition, entry_combos)

        # Generate and filter entries
        processed_entries = []
        total = sum([comb(len(entries), j + 1) for j in range(N)])
        if total > 1e6:
            warnings.warn(
                "Your pourbaix diagram includes {} entries and may take a long time to generate.".format(total)
            )

        # Parallel processing of multi-entry generation
        if nproc is not None:
            f = partial(self.process_multientry, prod_comp=total_comp)
            with Pool(nproc) as p:
                processed_entries = list(PBar(p.imap(f, entry_combos), total=total))
            processed_entries = list(filter(bool, processed_entries))
        # Serial processing of multi-entry generation
        else:
            for entry_combo in entry_combos:
                processed_entry = self.process_multientry(entry_combo, total_comp)
                if processed_entry is not None:
                    processed_entries.append(processed_entry)

        return processed_entries

    @staticmethod
    def process_multientry(entry_list, prod_comp, coeff_threshold=1e-4):
        """
        Static method for finding a multientry based on
        a list of entries and a product composition.
        Essentially checks to see if a valid aqueous
        reaction exists between the entries and the
        product composition and returns a MultiEntry
        with weights according to the coefficients if so.

        Args:
            entry_list ([Entry]): list of entries from which to
                create a MultiEntry
            prod_comp (Composition): composition constraint for setting
                weights of MultiEntry
            coeff_threshold (float): threshold of stoichiometric
                coefficients to filter, if weights are lower than
                this value, the entry is not returned
        """
        dummy_oh = [Composition("H"), Composition("O")]
        try:
            # Get balanced reaction coeffs, ensuring all < 0 or conc thresh
            # Note that we get reduced compositions for solids and non-reduced
            # compositions for ions because ions aren't normalized due to
            # their charge state.
            entry_comps = [e.composition for e in entry_list]
            rxn = Reaction(entry_comps + dummy_oh, [prod_comp])
            react_coeffs = [-rxn.get_coeff(comp) for comp in entry_comps]
            all_coeffs = react_coeffs + [rxn.get_coeff(prod_comp)]

            # Check if reaction coeff threshold met for pourbaix compounds
            # All reactant/product coefficients must be positive nonzero
            if all(coeff > coeff_threshold for coeff in all_coeffs):
                return MultiEntry(entry_list, weights=react_coeffs)

            return None
        except ReactionError:
            return None

    @staticmethod
    def get_pourbaix_domains(pourbaix_entries, limits=None):
        """
        Returns a set of pourbaix stable domains (i. e. polygons) in
        pH-V space from a list of pourbaix_entries

        This function works by using scipy's HalfspaceIntersection
        function to construct all of the 2-D polygons that form the
        boundaries of the planes corresponding to individual entry
        gibbs free energies as a function of pH and V. Hyperplanes
        of the form a*pH + b*V + 1 - g(0, 0) are constructed and
        supplied to HalfspaceIntersection, which then finds the
        boundaries of each pourbaix region using the intersection
        points.

        Args:
            pourbaix_entries ([PourbaixEntry]): Pourbaix entries
                with which to construct stable pourbaix domains
            limits ([[float]]): limits in which to do the pourbaix
                analysis

        Returns:
            Returns a dict of the form {entry: [boundary_points]}.
            The list of boundary points are the sides of the N-1
            dim polytope bounding the allowable ph-V range of each entry.
        """
        if limits is None:
            limits = [[-2, 16], [-4, 4]]

        # Get hyperplanes
        hyperplanes = [
            np.array([-PREFAC * entry.npH, -entry.nPhi, 0, -entry.energy]) * entry.normalization_factor
            for entry in pourbaix_entries
        ]
        hyperplanes = np.array(hyperplanes)
        hyperplanes[:, 2] = 1

        max_contribs = np.max(np.abs(hyperplanes), axis=0)
        g_max = np.dot(-max_contribs, [limits[0][1], limits[1][1], 0, 1])

        # Add border hyperplanes and generate HalfspaceIntersection
        border_hyperplanes = [
            [-1, 0, 0, limits[0][0]],
            [1, 0, 0, -limits[0][1]],
            [0, -1, 0, limits[1][0]],
            [0, 1, 0, -limits[1][1]],
            [0, 0, -1, 2 * g_max],
        ]
        hs_hyperplanes = np.vstack([hyperplanes, border_hyperplanes])
        interior_point = np.average(limits, axis=1).tolist() + [g_max]
        hs_int = HalfspaceIntersection(hs_hyperplanes, np.array(interior_point))

        # organize the boundary points by entry
        pourbaix_domains = {entry: [] for entry in pourbaix_entries}
        for intersection, facet in zip(hs_int.intersections, hs_int.dual_facets):
            for v in facet:
                if v < len(pourbaix_entries):
                    this_entry = pourbaix_entries[v]
                    pourbaix_domains[this_entry].append(intersection)

        # Remove entries with no pourbaix region
        pourbaix_domains = {k: v for k, v in pourbaix_domains.items() if v}
        pourbaix_domain_vertices = {}

        for entry, points in pourbaix_domains.items():
            points = np.array(points)[:, :2]
            # Initial sort to ensure consistency
            points = points[np.lexsort(np.transpose(points))]
            center = np.average(points, axis=0)
            points_centered = points - center

            # Sort points by cross product of centered points,
            # isn't strictly necessary but useful for plotting tools
            points_centered = sorted(points_centered, key=cmp_to_key(lambda x, y: x[0] * y[1] - x[1] * y[0]))
            points = points_centered + center

            # Create simplices corresponding to pourbaix boundary
            simplices = [Simplex(points[indices]) for indices in ConvexHull(points).simplices]
            pourbaix_domains[entry] = simplices
            pourbaix_domain_vertices[entry] = points

        return pourbaix_domains, pourbaix_domain_vertices

    def find_stable_entry(self, pH, V):
        """
        Finds stable entry at a pH,V condition
        Args:
            pH (float): pH to find stable entry
            V (float): V to find stable entry

        Returns:

        """
        energies_at_conditions = [e.normalized_energy_at_conditions(pH, V) for e in self.stable_entries]
        return self.stable_entries[np.argmin(energies_at_conditions)]

    def get_decomposition_energy(self, entry, pH, V):
        """
        Finds decomposition to most stable entries in eV/atom,
        supports vectorized inputs for pH and V

        Args:
            entry (PourbaixEntry): PourbaixEntry corresponding to
                compound to find the decomposition for
            pH (float, [float]): pH at which to find the decomposition
            V (float, [float]): voltage at which to find the decomposition

        Returns:
            Decomposition energy for the entry, i. e. the energy above
                the "pourbaix hull" in eV/atom at the given conditions
        """

        # Check composition consistency between entry and Pourbaix diagram:
        pbx_comp = Composition(self._elt_comp).fractional_composition
        entry_pbx_comp = Composition(
            {elt: coeff for elt, coeff in entry.composition.items() if elt not in ELEMENTS_HO}
        ).fractional_composition
        if entry_pbx_comp != pbx_comp:
            raise ValueError("Composition of stability entry does not match Pourbaix Diagram")
        entry_normalized_energy = entry.normalized_energy_at_conditions(pH, V)
        hull_energy = self.get_hull_energy(pH, V)
        decomposition_energy = entry_normalized_energy - hull_energy

        # Convert to eV/atom instead of eV/normalized formula unit
        decomposition_energy /= entry.normalization_factor
        decomposition_energy /= entry.composition.num_atoms
        return decomposition_energy

    def get_hull_energy(self, pH, V):
        """
        Gets the minimum energy of the pourbaix "basin" that is formed
        from the stable pourbaix planes.  Vectorized.

        Args:
            pH (float or [float]): pH at which to find the hull energy
            V (float or [float]): V at which to find the hull energy

        Returns:
            (float or [float]) minimum pourbaix energy at conditions

        """
        all_gs = np.array([e.normalized_energy_at_conditions(pH, V) for e in self.stable_entries])
        base = np.min(all_gs, axis=0)
        return base

    def get_stable_entry(self, pH, V):
        """
        Gets the stable entry at a given pH, V condition

        Args:
            pH (float): pH at a given condition
            V (float): V at a given condition

        Returns:
            (PourbaixEntry or MultiEntry): pourbaix or multi-entry
                corresponding ot the minimum energy entry at a given
                pH, V condition

        """
        all_gs = np.array([e.normalized_energy_at_conditions(pH, V) for e in self.stable_entries])
        return self.stable_entries[np.argmin(all_gs)]

    @property
    def stable_entries(self):
        """
        Returns the stable entries in the Pourbaix diagram.
        """
        return list(self._stable_domains.keys())

    @property
    def unstable_entries(self):
        """
        Returns all unstable entries in the Pourbaix diagram
        """
        return [e for e in self.all_entries if e not in self.stable_entries]

    @property
    def all_entries(self):
        """
        Return all entries used to generate the pourbaix diagram
        """
        return self._processed_entries

    @property
    def unprocessed_entries(self):
        """
        Return unprocessed entries
        """
        return self._unprocessed_entries

    def as_dict(self, include_unprocessed_entries=None):
        """
        Args:
            include_unprocessed_entries (): DEPRECATED. Whether to include unprocessed
                entries (equivalent to filter_solids=False). Serialization now includes
                all unprocessed entries by default. Set filter_solids=False before
                serializing to include unstable solids from the generated Pourbaix Diagram.

        Returns:
            MSONable dict.
        """
        if include_unprocessed_entries:
            warnings.warn(
                DeprecationWarning(
                    "The include_unprocessed_entries kwarg is deprecated! "
                    "Set filter_solids=True / False before serializing instead."
                )
            )
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "entries": [e.as_dict() for e in self._unprocessed_entries],
            "comp_dict": self._elt_comp,
            "conc_dict": self._conc_dict,
            "filter_solids": self.filter_solids,
        }
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (): Dict representation.

        Returns:
            PourbaixDiagram
        """
        decoded_entries = MontyDecoder().process_decoded(d["entries"])
        return cls(decoded_entries, d.get("comp_dict"), d.get("conc_dict"), d.get("filter_solids"))


class PourbaixPlotter:
    """
    A plotter class for phase diagrams.
    """

    def __init__(self, pourbaix_diagram):
        """
        Args:
            pourbaix_diagram (PourbaixDiagram): A PourbaixDiagram object.
        """
        self._pbx = pourbaix_diagram

    def show(self, *args, **kwargs):
        """
        Shows the pourbaix plot

        Args:
            *args: args to get_pourbaix_plot
            **kwargs: kwargs to get_pourbaix_plot

        Returns:
            None
        """
        plt = self.get_pourbaix_plot(*args, **kwargs)
        plt.show()

    def get_pourbaix_plot(self, limits=None, title="", label_domains=True, plt=None):
        """
        Plot Pourbaix diagram.

        Args:
            limits: 2D list containing limits of the Pourbaix diagram
                of the form [[xlo, xhi], [ylo, yhi]]
            title (str): Title to display on plot
            label_domains (bool): whether to label pourbaix domains
            plt (pyplot): Pyplot instance for plotting

        Returns:
            plt (pyplot) - matplotlib plot object with pourbaix diagram
        """
        if limits is None:
            limits = [[-2, 16], [-3, 3]]

        plt = plt or pretty_plot(16)

        xlim = limits[0]
        ylim = limits[1]

        h_line = np.transpose([[xlim[0], -xlim[0] * PREFAC], [xlim[1], -xlim[1] * PREFAC]])
        o_line = np.transpose([[xlim[0], -xlim[0] * PREFAC + 1.23], [xlim[1], -xlim[1] * PREFAC + 1.23]])
        neutral_line = np.transpose([[7, ylim[0]], [7, ylim[1]]])
        V0_line = np.transpose([[xlim[0], 0], [xlim[1], 0]])

        ax = plt.gca()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        lw = 3
        plt.plot(h_line[0], h_line[1], "r--", linewidth=lw)
        plt.plot(o_line[0], o_line[1], "r--", linewidth=lw)
        plt.plot(neutral_line[0], neutral_line[1], "k-.", linewidth=lw)
        plt.plot(V0_line[0], V0_line[1], "k-.", linewidth=lw)

        for entry, vertices in self._pbx._stable_domain_vertices.items():
            center = np.average(vertices, axis=0)
            x, y = np.transpose(np.vstack([vertices, vertices[0]]))
            plt.plot(x, y, "k-", linewidth=lw)

            if label_domains:
                plt.annotate(
                    generate_entry_label(entry),
                    center,
                    ha="center",
                    va="center",
                    fontsize=20,
                    color="b",
                ).draggable()

        plt.xlabel("pH")
        plt.ylabel("E (V)")
        plt.title(title, fontsize=20, fontweight="bold")
        return plt

    def plot_entry_stability(
        self,
        entry,
        pH_range=None,
        pH_resolution=100,
        V_range=None,
        V_resolution=100,
        e_hull_max=1,
        cmap="RdYlBu_r",
        **kwargs,
    ):
        """
        Args:
            entry ():
            pH_range ():
            pH_resolution ():
            V_range ():
            V_resolution ():
            e_hull_max ():
            cmap ():
            **kwargs ():

        Returns:

        """
        if pH_range is None:
            pH_range = [-2, 16]
        if V_range is None:
            V_range = [-3, 3]

        # plot the Pourbaix diagram
        plt = self.get_pourbaix_plot(**kwargs)
        pH, V = np.mgrid[
            pH_range[0] : pH_range[1] : pH_resolution * 1j,
            V_range[0] : V_range[1] : V_resolution * 1j,
        ]

        stability = self._pbx.get_decomposition_energy(entry, pH, V)

        # Plot stability map
        plt.pcolor(pH, V, stability, cmap=cmap, vmin=0, vmax=e_hull_max)
        cbar = plt.colorbar()
        cbar.set_label("Stability of {} (eV/atom)".format(generate_entry_label(entry)))

        # Set ticklabels
        # ticklabels = [t.get_text() for t in cbar.ax.get_yticklabels()]
        # ticklabels[-1] = '>={}'.format(ticklabels[-1])
        # cbar.ax.set_yticklabels(ticklabels)

        return plt

    def domain_vertices(self, entry):
        """
        Returns the vertices of the Pourbaix domain.

        Args:
            entry: Entry for which domain vertices are desired

        Returns:
            list of vertices
        """
        return self._pbx._stable_domain_vertices[entry]


def generate_entry_label(entry):
    """
    Generates a label for the pourbaix plotter

    Args:
        entry (PourbaixEntry or MultiEntry): entry to get a label for
    """
    if isinstance(entry, MultiEntry):
        return " + ".join([e.name for e in entry.entry_list])

    # TODO - a more elegant solution could be added later to Stringify
    # for example, the pattern re.sub(r"([-+][\d\.]*)", r"$^{\1}$", )
    # will convert B(OH)4- to B(OH)$_4^-$.
    # for this to work, the ion's charge always must be written AFTER
    # the sign (e.g., Fe+2 not Fe2+)
    string = entry.to_latex_string()
    return re.sub(r"()\[([^)]*)\]", r"\1$^{\2}$", string)
