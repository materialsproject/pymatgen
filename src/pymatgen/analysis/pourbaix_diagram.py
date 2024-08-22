"""
This module is intended to be used to compute Pourbaix diagrams of arbitrary compositions
and formation energies.
"""

from __future__ import annotations

import itertools
import logging
import re
import warnings
from copy import deepcopy
from functools import cmp_to_key, partial
from multiprocessing import Pool
from typing import TYPE_CHECKING, no_type_check

import numpy as np
from monty.json import MontyDecoder, MSONable
from scipy.spatial import ConvexHull, HalfspaceIntersection
from scipy.special import comb

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError
from pymatgen.core import Composition, Element
from pymatgen.core.ion import Ion
from pymatgen.entries.compatibility import MU_H2O
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.util.coord import Simplex
from pymatgen.util.due import Doi, due
from pymatgen.util.plotting import pretty_plot
from pymatgen.util.string import Stringify

if TYPE_CHECKING:
    from typing import Any

    import matplotlib.pyplot as plt
    from typing_extensions import Self

__author__ = "Sai Jayaraman"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.4"
__maintainer__ = "Joseph Montoya"
__credits__ = "Arunima Singh, Joseph Montoya, Anjli Patel"
__email__ = "joseph.montoya@tri.global"
__status__ = "Production"
__date__ = "Nov 1, 2012"

# If you use this module in your work, consider citing:
due.cite(
    Doi("10.1103/PhysRevB.85.235438"),
    description="Prediction of solid-aqueous equilibria: Scheme to combine first-principles "
    "calculations of solids with experimental aqueous states",
)
due.cite(
    Doi("10.1021/acs.chemmater.7b03980"),
    description="Electrochemical Stability of Metastable Materials",
)
due.cite(
    Doi("10.1021/acs.chemmater.7b03980"),
    description="Fast computation of many-element Pourbaix diagrams",
)

logger = logging.getLogger(__name__)

PREFAC = 0.0591


# TODO: Revise to more closely reflect PDEntry, invoke from energy/composition
# TODO: PourbaixEntries depend implicitly on having entry energies be
# formation energies, should be a better way to get from raw energies
# TODO: uncorrected_energy is a bit of a misnomer, but not sure what to rename
class PourbaixEntry(MSONable, Stringify):
    """
    An object encompassing all data relevant to a solid or ion
    in a Pourbaix diagram. Each bulk solid/ion has an energy
    g of the form: e = e0 + 0.0591 log10(conc) - nO mu_H2O
    + (nH - 2nO) pH + phi (-nH + 2nO + q).

    Note that the energies corresponding to the input entries
    should be formation energies with respect to hydrogen and
    oxygen gas in order for the Pourbaix diagram formalism to
    work. This may be changed to be more flexible in the future.
    """

    def __init__(self, entry, entry_id=None, concentration=1e-6):
        """
        Args:
            entry (ComputedEntry | ComputedStructureEntry | PDEntry | IonEntry): An entry object
            entry_id (str): A string id for the entry
            concentration (float): Concentration of the entry in M. Defaults to 1e-6.
        """
        self.entry = entry
        if isinstance(entry, IonEntry):
            self.concentration = concentration
            self.phase_type = "Ion"
            self.charge = entry.ion.charge
        else:
            self.concentration = 1.0
            self.phase_type = "Solid"
            self.charge = 0
        self.uncorrected_energy = entry.energy
        if entry_id is not None:
            self.entry_id = entry_id
        elif getattr(entry, "entry_id", None):
            self.entry_id = entry.entry_id
        else:
            self.entry_id = None

    @property
    def npH(self):
        """The number of H."""
        return self.entry.composition.get("H", 0) - 2 * self.entry.composition.get("O", 0)

    @property
    def nH2O(self):
        """The number of H2O."""
        return self.entry.composition.get("O", 0)

    @property
    def nPhi(self):
        """The number of electrons."""
        return self.npH - self.charge

    @property
    def name(self):
        """The entry's name."""
        if self.phase_type == "Solid":
            return f"{self.entry.reduced_formula}(s)"

        return self.entry.name

    @property
    def energy(self):
        """Total energy of the Pourbaix entry (at pH, V = 0 vs. SHE)."""
        # Note: this implicitly depends on formation energies as input
        return self.uncorrected_energy + self.conc_term - (MU_H2O * self.nH2O)

    @property
    def energy_per_atom(self):
        """Energy per atom of the Pourbaix entry."""
        return self.energy / self.composition.num_atoms

    @property
    def elements(self):
        """Elements in the entry."""
        return self.entry.elements

    def energy_at_conditions(self, pH, V):
        """Get free energy for a given pH and V.

        Args:
            pH (float): pH at which to evaluate free energy
            V (float): voltage at which to evaluate free energy

        Returns:
            free energy at conditions
        """
        return self.energy + self.npH * PREFAC * pH + self.nPhi * V

    def get_element_fraction(self, element):
        """Get the elemental fraction of a given non-OH element.

        Args:
            element (Element or str): string or element corresponding
                to element to get from composition

        Returns:
            fraction of element / sum(all non-OH elements)
        """
        return self.composition.get(element) * self.normalization_factor

    @property
    def normalized_energy(self):
        """Energy normalized by number of non H or O atoms, e.g.
        for Zn2O6, energy / 2 or for AgTe3(OH)3, energy / 4.
        """
        return self.energy * self.normalization_factor

    def normalized_energy_at_conditions(self, pH, V):
        """Energy at an electrochemical condition, compatible with
        numpy arrays for pH/V input.

        Args:
            pH (float): pH at condition
            V (float): applied potential at condition

        Returns:
            energy normalized by number of non-O/H atoms at condition
        """
        return self.energy_at_conditions(pH, V) * self.normalization_factor

    @property
    def conc_term(self):
        """The concentration contribution to the free energy. Should only be present
        when there are ions in the entry.
        """
        return PREFAC * np.log10(self.concentration)

    # TODO: not sure if these are strictly necessary with refactor
    def as_dict(self):
        """Get dict which contains Pourbaix Entry data.
        Note that the pH, voltage, H2O factors are always calculated when
        constructing a PourbaixEntry object.
        """
        dct = {"@module": type(self).__module__, "@class": type(self).__name__}
        if isinstance(self.entry, IonEntry):
            dct["entry_type"] = "Ion"
        else:
            dct["entry_type"] = "Solid"
        dct["entry"] = self.entry.as_dict()
        dct["concentration"] = self.concentration
        dct["entry_id"] = self.entry_id
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Invokes a PourbaixEntry from a dictionary."""
        entry_type = dct["entry_type"]
        entry = (
            IonEntry.from_dict(dct["entry"]) if entry_type == "Ion" else MontyDecoder().process_decoded(dct["entry"])
        )
        entry_id = dct["entry_id"]
        concentration = dct["concentration"]
        return cls(entry, entry_id, concentration)

    @property
    def normalization_factor(self):
        """Sum of number of atoms minus the number of H and O in composition."""
        return 1.0 / (self.num_atoms - self.composition.get("H", 0) - self.composition.get("O", 0))

    @property
    def composition(self):
        """Composition."""
        return self.entry.composition

    @property
    def num_atoms(self):
        """Number of atoms in current formula. Useful for normalization."""
        return self.composition.num_atoms

    def to_pretty_string(self) -> str:
        """A pretty string representation."""
        if self.phase_type == "Solid":
            return f"{self.entry.reduced_formula}(s)"

        return self.entry.name

    def __repr__(self):
        energy, npH, nPhi, nH2O, entry_id = self.energy, self.npH, self.nPhi, self.nH2O, self.entry_id
        return (
            f"{type(self).__name__}({self.entry.composition} with {energy=:.4f}, {npH=}, "
            f"{nPhi=}, {nH2O=}, {entry_id=})"
        )


class MultiEntry(PourbaixEntry):
    """PourbaixEntry-like object for constructing multi-elemental Pourbaix diagrams."""

    def __init__(self, entry_list, weights=None):
        """Initialize a MultiEntry.

        Args:
            entry_list ([PourbaixEntry]): List of component PourbaixEntries
            weights ([float]): Weights associated with each entry. Default is None
        """
        self.weights = weights or [1.0] * len(entry_list)
        self.entry_list = entry_list

    def __getattr__(self, attr):
        """
        Because most of the attributes here are just weighted averages of the entry_list,
        we save some space by having a set of conditionals to define the attributes.
        """
        # Attributes that are weighted averages of entry attributes
        if attr in ["energy", "npH", "nH2O", "nPhi", "conc_term", "composition", "uncorrected_energy", "elements"]:
            # TODO: Composition could be changed for compat with sum
            start = Composition() if attr == "composition" else 0
            weighted_values = (
                getattr(entry, attr) * weight for entry, weight in zip(self.entry_list, self.weights, strict=True)
            )
            return sum(weighted_values, start)

        # Attributes that are just lists of entry attributes
        if attr in ["entry_id", "phase_type"]:
            return [getattr(entry, attr) for entry in self.entry_list]

        # normalization_factor, num_atoms should work from superclass
        return self.__getattribute__(attr)

    @property
    def name(self):
        """MultiEntry name, i.e. the name of each entry joined by ' + '."""
        return " + ".join(entry.name for entry in self.entry_list)

    def __repr__(self):
        energy, npH, nPhi, nH2O, entry_id = self.energy, self.npH, self.nPhi, self.nH2O, self.entry_id
        cls_name, species = type(self).__name__, self.name
        return f"Pourbaix{cls_name}({energy=:.4f}, {npH=}, {nPhi=}, {nH2O=}, {entry_id=}, {species=})"

    def as_dict(self):
        """Get MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "entry_list": [entry.as_dict() for entry in self.entry_list],
            "weights": self.weights,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            MultiEntry
        """
        entry_list = [PourbaixEntry.from_dict(entry) for entry in dct.get("entry_list", ())]
        return cls(entry_list, dct.get("weights"))


# TODO: this class isn't particularly useful in its current form, could be
# refactored to include information about the reference solid
class IonEntry(PDEntry):
    """
    Object similar to PDEntry, but contains an Ion object instead of a
    Composition object.

    Attributes:
        name (str): A name for the entry. This is the string shown in the phase diagrams.
            By default, this is the reduced formula for the composition, but can be
            set to some other string for display purposes.
    """

    def __init__(self, ion: Ion, energy: float, name: str | None = None, attribute=None):
        """
        Args:
            ion: Ion object
            energy: Energy for composition.
            name: Optional parameter to name the entry. Defaults to the
                chemical formula.
            attribute: Optional attribute of the entry, e.g. band gap.
        """
        self.ion = ion
        # Auto-assign name
        name = name or self.ion.reduced_formula
        super().__init__(composition=ion.composition, energy=energy, name=name, attribute=attribute)

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Get an IonEntry object from a dict."""
        return cls(Ion.from_dict(dct["ion"]), dct["energy"], dct.get("name"), dct.get("attribute"))

    def as_dict(self):
        """Create a dict of composition, energy, and ion name."""
        return {"ion": self.ion.as_dict(), "energy": self.energy, "name": self.name}

    def __repr__(self):
        return f"IonEntry : {self.composition} with energy = {self.energy:.4f}"


def ion_or_solid_comp_object(formula):
    """Get an Ion or Composition object given a formula.

    Args:
        formula: String formula. Eg. of ion: NaOH(aq), Na[+];
            Eg. of solid: Fe2O3(s), Fe(s), Na2O

    Returns:
        Composition/Ion object
    """
    if re.match(r"\[([^\[\]]+)\]|\(aq\)", formula):
        comp_obj = Ion.from_formula(formula)
    elif re.search(r"\(s\)", formula):
        comp_obj = Composition(formula[:-3])
    else:
        comp_obj = Composition(formula)
    return comp_obj


ELEMENTS_HO = {Element("H"), Element("O")}


# TODO: the solids filter breaks some of the functionality of the
# heatmap plotter, because the reference states for decomposition
# don't include oxygen/hydrogen in the OER/HER regions


# TODO: create a from_phase_diagram class method for non-formation energy invocation
# TODO: invocation from a MultiEntry entry list could be a bit more robust
# TODO: serialization is still a bit rough around the edges
class PourbaixDiagram(MSONable):
    """Create a Pourbaix diagram from entries."""

    def __init__(
        self,
        entries: list[PourbaixEntry] | list[MultiEntry],
        comp_dict: dict[str, float] | None = None,
        conc_dict: dict[str, float] | None = None,
        filter_solids: bool = True,
        nproc: int | None = None,
    ):
        """
        Args:
            entries ([PourbaixEntry] or [MultiEntry]): Entries list
                containing Solids and Ions or a list of MultiEntries
            comp_dict (dict[str, float]): Dictionary of compositions,
                defaults to equal parts of each elements
            conc_dict (dict[str, float]): Dictionary of ion concentrations,
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
            nproc (int): number of processes to generate multi-entries with
                in parallel. Defaults to None (serial processing).
        """
        entries = deepcopy(entries)
        self.filter_solids = filter_solids

        # Get non-OH elements
        self.pbx_elts = list(
            set(itertools.chain.from_iterable([entry.composition.elements for entry in entries])) - ELEMENTS_HO
        )
        self.dim = len(self.pbx_elts) - 1

        # Process multi-entry inputs
        if isinstance(entries[0], MultiEntry):
            self._processed_entries = entries
            # Extract individual entries
            single_entries = list(set(itertools.chain.from_iterable([entry.entry_list for entry in entries])))
            self._unprocessed_entries = single_entries
            self._filtered_entries = single_entries
            self._conc_dict = None
            self._elt_comp = {k: v for k, v in entries[0].composition.items() if k not in ELEMENTS_HO}
            self._multi_element = True

        # Process single entry inputs
        else:
            # Set default conc/comp dicts
            if not comp_dict:
                comp_dict = {elt.symbol: 1 / len(self.pbx_elts) for elt in self.pbx_elts}
            if not conc_dict:
                conc_dict = {elt.symbol: 1e-6 for elt in self.pbx_elts}
            self._conc_dict = conc_dict

            self._elt_comp = comp_dict
            self.pourbaix_elements = self.pbx_elts

            solid_entries = [entry for entry in entries if entry.phase_type == "Solid"]
            ion_entries = [entry for entry in entries if entry.phase_type == "Ion"]

            # If a conc_dict is specified, override individual entry concentrations
            for entry in ion_entries:
                ion_elts = list(set(entry.elements) - ELEMENTS_HO)
                # TODO: the logic here for ion concentration setting is in two
                # places, in PourbaixEntry and here, should be consolidated
                if len(ion_elts) == 1:
                    entry.concentration = conc_dict[ion_elts[0].symbol] * entry.normalization_factor
                elif len(ion_elts) > 1 and not entry.concentration:
                    raise ValueError("Elemental concentration not compatible with multi-element ions")

            self._unprocessed_entries = solid_entries + ion_entries

            if len(solid_entries + ion_entries) != len(entries):
                raise ValueError('All supplied entries must have a phase type of either "Solid" or "Ion"')

            if self.filter_solids:
                # O is 2.46 b/c pbx entry finds energies referenced to H2O
                entries_HO = [ComputedEntry("H", 0), ComputedEntry("O", 2.46)]
                solid_pd = PhaseDiagram(solid_entries + entries_HO)
                solid_entries = list(set(solid_pd.stable_entries) - set(entries_HO))

            self._filtered_entries = solid_entries + ion_entries
            if len(comp_dict) > 1:
                self._multi_element = True
                self._processed_entries = self._preprocess_pourbaix_entries(self._filtered_entries, nproc=nproc)
            else:
                self._processed_entries = self._filtered_entries
                self._multi_element = False

        self._stable_domains, self._stable_domain_vertices = self.get_pourbaix_domains(self._processed_entries)

    def _convert_entries_to_points(self, pourbaix_entries):
        """
        Args:
            pourbaix_entries ([PourbaixEntry]): list of Pourbaix entries
                to process into vectors in nph-nphi-composition space.

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

    def _get_hull_in_nph_nphi_space(self, entries) -> tuple[list[PourbaixEntry], list[Simplex]]:
        """Generate convex hull of Pourbaix diagram entries in composition,
        npH, and nphi space. This enables filtering of multi-entries
        such that only compositionally stable combinations of entries
        are included.

        Args:
            entries ([PourbaixEntry]): list of PourbaixEntries to construct
                the convex hull

        Returns:
            tuple[list[PourbaixEntry], list[Simplex]]: PourbaixEntry list and stable
                facets corresponding to that list
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
        min_entries = [next(iter(grouped_entries)) for comp, grouped_entries in grouped_by_composition]
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
        facets = [facet for facet in hull.simplices if len(points) - 1 not in facet]

        if self.dim > 1:
            logger.debug("Filtering facets by Pourbaix composition")
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
        """Generate multi-entries for Pourbaix diagram.

        Args:
            entries ([PourbaixEntry]): list of PourbaixEntries to preprocess
                into MultiEntries
            nproc (int): number of processes to be used in parallel
                treatment of entry combos

        Returns:
            list[MultiEntry]: stable MultiEntry candidates
        """
        # Get composition
        tot_comp = Composition(self._elt_comp)

        min_entries, valid_facets = self._get_hull_in_nph_nphi_space(entries)

        combos = []
        for facet in valid_facets:
            for idx in range(1, self.dim + 2):
                these_combos = []
                for combo in itertools.combinations(facet, idx):
                    these_entries = [min_entries[i] for i in combo]
                    these_combos.append(frozenset(these_entries))
                combos.append(these_combos)

        all_combos = set(itertools.chain.from_iterable(combos))

        list_combos = []
        for idx in all_combos:
            list_combos.append(list(idx))
        all_combos = list_combos

        multi_entries = []

        # Parallel processing of multi-entry generation
        if nproc is not None:
            func = partial(self.process_multientry, prod_comp=tot_comp)
            with Pool(nproc) as proc_pool:
                multi_entries = list(proc_pool.imap(func, all_combos))
            multi_entries = list(filter(bool, multi_entries))
        else:
            # Serial processing of multi-entry generation
            for combo in all_combos:
                if multi_entry := self.process_multientry(combo, prod_comp=tot_comp):
                    multi_entries.append(multi_entry)

        return multi_entries

    def _generate_multielement_entries(self, entries, nproc=None):
        """
        Create entries for multi-element Pourbaix construction.

        This works by finding all possible linear combinations
        of entries that can result in the specified composition
        from the initialized comp_dict.

        Args:
            entries ([PourbaixEntries]): list of Pourbaix entries
                to process into MultiEntries
            nproc (int): number of processes to be used in parallel
                treatment of entry combos
        """
        n_elems = len(self._elt_comp)  # No. of elements
        total_comp = Composition(self._elt_comp)

        # generate all combinations of compounds that have all elements
        entry_combos = [itertools.combinations(entries, idx + 1) for idx in range(n_elems)]
        entry_combos = itertools.chain.from_iterable(entry_combos)

        entry_combos = filter(lambda x: total_comp < MultiEntry(x).composition, entry_combos)

        # Generate and filter entries
        processed_entries = []
        total = sum(comb(len(entries), idx + 1) for idx in range(n_elems))
        if total > 1e6:
            warnings.warn(f"Your Pourbaix diagram includes {total} entries and may take a long time to generate.")

        # Parallel processing of multi-entry generation
        if nproc is not None:
            func = partial(self.process_multientry, prod_comp=total_comp)
            with Pool(nproc) as proc_pool:
                processed_entries = list(proc_pool.imap(func, entry_combos))
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
        """Static method for finding a multientry based on
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
            entry_comps = [entry.composition for entry in entry_list]
            rxn = Reaction(entry_comps + dummy_oh, [prod_comp])
            react_coeffs = [-coeff for coeff in rxn.coeffs[: len(entry_list)]]
            all_coeffs = [*react_coeffs, rxn.get_coeff(prod_comp)]

            # Check if reaction coeff threshold met for Pourbaix compounds
            # All reactant/product coefficients must be positive nonzero
            if all(coeff > coeff_threshold for coeff in all_coeffs):
                return MultiEntry(entry_list, weights=react_coeffs)

            return None
        except ReactionError:
            return None

    @staticmethod
    def get_pourbaix_domains(pourbaix_entries, limits=None):
        """Get a set of Pourbaix stable domains (i.e. polygons) in
        pH-V space from a list of pourbaix_entries.

        This function works by using scipy's HalfspaceIntersection
        function to construct all of the 2-D polygons that form the
        boundaries of the planes corresponding to individual entry
        gibbs free energies as a function of pH and V. Hyperplanes
        of the form a*pH + b*V + 1 - g(0, 0) are constructed and
        supplied to HalfspaceIntersection, which then finds the
        boundaries of each Pourbaix region using the intersection
        points.

        Args:
            pourbaix_entries ([PourbaixEntry]): Pourbaix entries
                with which to construct stable Pourbaix domains
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
        interior_point = [*np.mean(limits, axis=1).tolist(), g_max]
        hs_int = HalfspaceIntersection(hs_hyperplanes, np.array(interior_point))

        # organize the boundary points by entry
        pourbaix_domains = {entry: [] for entry in pourbaix_entries}
        for intersection, facet in zip(hs_int.intersections, hs_int.dual_facets, strict=True):
            for v in facet:
                if v < len(pourbaix_entries):
                    this_entry = pourbaix_entries[v]
                    pourbaix_domains[this_entry].append(intersection)

        # Remove entries with no Pourbaix region
        pourbaix_domains = {k: v for k, v in pourbaix_domains.items() if v}
        pourbaix_domain_vertices = {}

        for entry, points in pourbaix_domains.items():
            points = np.array(points)[:, :2]
            # Initial sort to ensure consistency
            points = points[np.lexsort(np.transpose(points))]
            center = np.mean(points, axis=0)
            points_centered = points - center

            # Sort points by cross product of centered points,
            # isn't strictly necessary but useful for plotting tools
            points_centered = sorted(points_centered, key=cmp_to_key(lambda x, y: x[0] * y[1] - x[1] * y[0]))
            points = points_centered + center

            # Create simplices corresponding to Pourbaix boundary
            simplices = [Simplex(points[indices]) for indices in ConvexHull(points).simplices]
            pourbaix_domains[entry] = simplices
            pourbaix_domain_vertices[entry] = points

        return pourbaix_domains, pourbaix_domain_vertices

    def find_stable_entry(self, pH, V):
        """Find stable entry at a pH,V condition.

        Args:
            pH (float): pH to find stable entry
            V (float): V to find stable entry.

        Returns:
            PourbaixEntry: stable entry at pH, V
        """
        energies_at_conditions = [entry.normalized_energy_at_conditions(pH, V) for entry in self.stable_entries]
        return self.stable_entries[np.argmin(energies_at_conditions)]

    def get_decomposition_energy(self, entry, pH, V):
        """Find decomposition to most stable entries in eV/atom,
        supports vectorized inputs for pH and V.

        Args:
            entry (PourbaixEntry): PourbaixEntry corresponding to
                compound to find the decomposition for
            pH (float, list[float]): pH at which to find the decomposition
            V (float, list[float]): voltage at which to find the decomposition

        Returns:
            Decomposition energy for the entry, i.e. the energy above
                the "Pourbaix hull" in eV/atom at the given conditions
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

    def get_hull_energy(self, pH: float | list[float], V: float | list[float]) -> np.ndarray:
        """Get the minimum energy of the Pourbaix "basin" that is formed
        from the stable Pourbaix planes. Vectorized.

        Args:
            pH (float | list[float]): pH at which to find the hull energy
            V (float | list[float]): V at which to find the hull energy

        Returns:
            np.array: minimum Pourbaix energy at conditions
        """
        all_gs = np.array([entry.normalized_energy_at_conditions(pH, V) for entry in self.stable_entries])
        return np.min(all_gs, axis=0)

    def get_stable_entry(self, pH, V):
        """Get the stable entry at a given pH, V condition.

        Args:
            pH (float): pH at a given condition
            V (float): V at a given condition

        Returns:
            PourbaixEntry | MultiEntry: Pourbaix or multi-entry
                corresponding to the minimum energy entry at a given pH, V condition
        """
        all_gs = np.array([entry.normalized_energy_at_conditions(pH, V) for entry in self.stable_entries])
        return self.stable_entries[np.argmin(all_gs)]

    @property
    def stable_entries(self):
        """The stable entries in the Pourbaix diagram."""
        return list(self._stable_domains)

    @property
    def unstable_entries(self):
        """All unstable entries in the Pourbaix diagram."""
        return [entry for entry in self.all_entries if entry not in self.stable_entries]

    @property
    def all_entries(self):
        """All entries used to generate the Pourbaix diagram."""
        return self._processed_entries

    @property
    def unprocessed_entries(self):
        """Unprocessed entries."""
        return self._unprocessed_entries

    def as_dict(self):
        """Get MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "entries": [entry.as_dict() for entry in self._unprocessed_entries],
            "comp_dict": self._elt_comp,
            "conc_dict": self._conc_dict,
            "filter_solids": self.filter_solids,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            PourbaixDiagram
        """
        decoded_entries = MontyDecoder().process_decoded(dct["entries"])
        return cls(
            decoded_entries,
            comp_dict=dct.get("comp_dict"),
            conc_dict=dct.get("conc_dict"),
            filter_solids=bool(dct.get("filter_solids")),
        )


class PourbaixPlotter:
    """A plotter class for phase diagrams."""

    def __init__(self, pourbaix_diagram):
        """
        Args:
            pourbaix_diagram (PourbaixDiagram): A PourbaixDiagram object.
        """
        self._pbx = pourbaix_diagram

    def show(self, *args, **kwargs):
        """Show the Pourbaix plot.

        Args:
            *args: args to get_pourbaix_plot
            **kwargs: kwargs to get_pourbaix_plot
        """
        plt = self.get_pourbaix_plot(*args, **kwargs)
        plt.show()

    @no_type_check
    def get_pourbaix_plot(
        self,
        limits: tuple[float, float] | None = None,
        title: str = "",
        label_domains: bool = True,
        label_fontsize: int = 20,
        show_water_lines: bool = True,
        show_neutral_axes: bool = True,
        ax: plt.Axes = None,
    ) -> plt.Axes:
        """
        Plot Pourbaix diagram.

        Args:
            limits: 2D list containing limits of the Pourbaix diagram
                of the form [[xlo, xhi], [ylo, yhi]]
            title (str): Title to display on plot
            label_domains (bool): whether to label Pourbaix domains
            label_fontsize: font size for domain labels
            show_water_lines: whether to show dashed lines indicating the region
                of water stability.
            show_neutral_axes; whether to show dashed horizontal and vertical lines
                at 0 V and pH 7, respectively.
            ax (Axes): Matplotlib Axes instance for plotting

        Returns:
            Axes: matplotlib Axes object with Pourbaix diagram
        """
        if limits is None:
            limits = [[-2, 16], [-3, 3]]

        ax = ax or pretty_plot(16)

        xlim, ylim = limits
        lw = 3

        if show_water_lines:
            h_line = np.transpose([[xlim[0], -xlim[0] * PREFAC], [xlim[1], -xlim[1] * PREFAC]])
            o_line = np.transpose([[xlim[0], -xlim[0] * PREFAC + 1.23], [xlim[1], -xlim[1] * PREFAC + 1.23]])
            ax.plot(h_line[0], h_line[1], "r--", linewidth=lw)
            ax.plot(o_line[0], o_line[1], "r--", linewidth=lw)

        if show_neutral_axes:
            neutral_line = np.transpose([[7, ylim[0]], [7, ylim[1]]])
            V0_line = np.transpose([[xlim[0], 0], [xlim[1], 0]])
            ax.plot(neutral_line[0], neutral_line[1], "k-.", linewidth=lw)
            ax.plot(V0_line[0], V0_line[1], "k-.", linewidth=lw)

        for entry, vertices in self._pbx._stable_domain_vertices.items():
            center = np.mean(vertices, axis=0)
            x, y = np.transpose(np.vstack([vertices, vertices[0]]))
            ax.plot(x, y, "k-", linewidth=lw)

            if label_domains:
                ax.annotate(
                    generate_entry_label(entry),
                    center,
                    ha="center",
                    va="center",
                    fontsize=label_fontsize,
                    color="b",
                ).draggable()

        ax.set_title(title, fontsize=20, fontweight="bold")
        ax.set(xlabel="pH", ylabel="E (V)", xlim=xlim, ylim=ylim)
        return ax

    @no_type_check
    def plot_entry_stability(
        self,
        entry: Any,
        pH_range: tuple[float, float] = (-2, 16),
        pH_resolution: int = 100,
        V_range: tuple[float, float] = (-3, 3),
        V_resolution: int = 100,
        e_hull_max: float = 1,
        cmap: str = "RdYlBu_r",
        ax: plt.Axes | None = None,
        **kwargs: Any,
    ) -> plt.Axes:
        """
        Plots the stability of an entry in the Pourbaix diagram.

        Args:
            entry (Any): The entry to plot stability for.
            pH_range (tuple[float, float], optional): pH range for the plot. Defaults to (-2, 16).
            pH_resolution (int, optional): pH resolution. Defaults to 100.
            V_range (tuple[float, float], optional): Voltage range for the plot. Defaults to (-3, 3).
            V_resolution (int, optional): Voltage resolution. Defaults to 100.
            e_hull_max (float, optional): Maximum energy above the hull. Defaults to 1.
            cmap (str, optional): Colormap for the plot. Defaults to "RdYlBu_r".
            ax (Axes, optional): Existing matplotlib Axes object for plotting. Defaults to None.
            **kwargs (Any): Additional keyword arguments passed to `get_pourbaix_plot`.

        Returns:
            plt.Axes: Matplotlib Axes object with the plotted stability.
        """
        # Plot the Pourbaix diagram
        ax = self.get_pourbaix_plot(ax=ax, **kwargs)
        pH, V = np.mgrid[
            pH_range[0] : pH_range[1] : pH_resolution * 1j,
            V_range[0] : V_range[1] : V_resolution * 1j,
        ]

        stability = self._pbx.get_decomposition_energy(entry, pH, V)

        # Plot stability map
        cax = ax.pcolor(pH, V, stability, cmap=cmap, vmin=0, vmax=e_hull_max)
        cbar = ax.figure.colorbar(cax)
        cbar.set_label(f"Stability of {generate_entry_label(entry)} (eV/atom)")

        # Set ticklabels
        # ticklabels = [t.get_text() for t in cbar.ax.get_yticklabels()]
        # ticklabels[-1] = f">={ticklabels[-1]}"
        # cbar.ax.set_yticklabels(ticklabels)

        return ax

    def domain_vertices(self, entry):
        """Get the vertices of the Pourbaix domain.

        Args:
            entry: Entry for which domain vertices are desired

        Returns:
            list of vertices
        """
        return self._pbx._stable_domain_vertices[entry]


def generate_entry_label(entry):
    """
    Generates a label for the Pourbaix plotter.

    Args:
        entry (PourbaixEntry or MultiEntry): entry to get a label for
    """
    if isinstance(entry, MultiEntry):
        return " + ".join(entry.name for entry in entry.entry_list)

    # TODO - a more elegant solution could be added later to Stringify
    # for example, the pattern re.sub(r"([-+][\d\.]*)", r"$^{\1}$", )
    # will convert B(OH)4- to B(OH)$_4^-$.
    # for this to work, the ion's charge always must be written AFTER
    # the sign (e.g., Fe+2 not Fe2+)
    string = entry.to_latex_string()
    return re.sub(r"()\[([^)]*)\]", r"\1$^{\2}$", string)
