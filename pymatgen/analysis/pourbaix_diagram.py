# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import logging
import numpy as np
import itertools
import re
from copy import deepcopy
from functools import cmp_to_key, partial, lru_cache
from monty.json import MSONable, MontyDecoder
from six.moves import zip
from multiprocessing import Pool
import warnings

from scipy.spatial import ConvexHull, HalfspaceIntersection

try:
    from scipy.special import comb
except ImportError:
    from scipy.misc import comb
from pymatgen.util.coord import Simplex
from pymatgen.util.string import latexify
from pymatgen.util.plotting import pretty_plot
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.core.ion import Ion
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry

try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda x: x

__author__ = "Sai Jayaraman"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.3"
__maintainer__ = "Joseph Montoya"
__credits__ = "Arunima Singh, Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Production"
__date__ = "Nov 1, 2012"

logger = logging.getLogger(__name__)

MU_H2O = -2.4583
PREFAC = 0.0591


# TODO: Revise to more closely reflect PDEntry, invoke from energy/composition
# TODO: PourbaixEntries depend implicitly on having entry energies be
#       formation energies, should be a better way to get from raw energies
# TODO: uncorrected_energy is a bit of a misnomer, but not sure what to rename
class PourbaixEntry(MSONable):
    """
    An object encompassing all data relevant to a solid or ion
    in a pourbaix diagram.  Each bulk solid/ion has an energy
    g of the form: e = e0 + 0.0591 log10(conc) - nO mu_H2O
    + (nH - 2nO) pH + phi (-nH + 2nO + q)

    Note that the energies corresponding to the input entries
    should be formation energies with respect to hydrogen and
    oxygen gas in order for the pourbaix diagram formalism to
    work. This may be changed to be more flexible in the future.

    Args:
        entry (ComputedEntry/ComputedStructureEntry/PDEntry/IonEntry): An
            entry object
    """

    def __init__(self, entry, entry_id=None, concentration=1e-6):
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
        return self.entry.composition.get("H", 0.) \
               - 2 * self.entry.composition.get("O", 0.)

    @property
    def nH2O(self):
        return self.entry.composition.get("O", 0.)

    @property
    def nPhi(self):
        return self.npH - self.charge

    @property
    def name(self):
        if self.phase_type == "Solid":
            return self.entry.composition.reduced_formula + "(s)"
        elif self.phase_type == "Ion":
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
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
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
        Invokes
        """
        entry_type = d["entry_type"]
        if entry_type == "Ion":
            entry = IonEntry.from_dict(d["entry"])
        else:
            entry = PDEntry.from_dict(d["entry"])
        entry_id = d["entry_id"]
        concentration = d["concentration"]
        return PourbaixEntry(entry, entry_id, concentration)

    @property
    def normalization_factor(self):
        """
        Sum of number of atoms minus the number of H and O in composition
        """
        return 1.0 / (self.num_atoms - self.composition.get('H', 0)
                      - self.composition.get('O', 0))

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

    def __repr__(self):
        return "Pourbaix Entry : {} with energy = {:.4f}, npH = {}, " \
               "nPhi = {}, nH2O = {}, entry_id = {} ".format(
            self.entry.composition, self.energy, self.npH,
            self.nPhi, self.nH2O, self.entry_id)

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
        if item in ["energy", "npH", "nH2O", "nPhi", "conc_term",
                    "composition", "uncorrected_energy"]:
            # TODO: Composition could be changed for compat with sum
            if item == "composition":
                start = Composition({})
            else:
                start = 0
            return sum([getattr(e, item) * w
                        for e, w in zip(self.entry_list, self.weights)], start)
        # Attributes that are just lists of entry attributes
        elif item in ["entry_id", "phase_type"]:
            return [getattr(e, item) for e in self.entry_list]
        # normalization_factor, num_atoms should work from superclass
        return self.__getattribute__(item)

    @property
    def name(self):
        """
        Multientry name, i. e. the name of each entry joined by ' + '
        """
        return " + ".join([e.name for e in self.entry_list])

    def __repr__(self):
        return "Multiple Pourbaix Entry : with energy = {:.4f}, npH = {}, " \
               "nPhi = {}, nH2O = {}, entry_id = {}, species: {}".format(
            self.energy, self.npH, self.nPhi, self.nH2O,
            self.entry_id, self.name)

    def __str__(self):
        return self.__repr__()

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "entry_list": [e.as_dict() for e in self.entry_list],
                "weights": self.weights}

    @classmethod
    def from_dict(cls, d):
        entry_list = [PourbaixEntry.from_dict(e) for e in d.get("entry_list")]
        return cls(entry_list, d.get("weights"))


# TODO: this class isn't particularly useful in its current form, could be
#       refactored to include information about the reference solid
class IonEntry(PDEntry):
    """
    Object similar to PDEntry, but contains an Ion object instead of a
    Composition object.

    Args:
        comp: Ion object
        energy: Energy for composition.
        name: Optional parameter to name the entry. Defaults to the
            chemical formula.

    .. attribute:: name

        A name for the entry. This is the string shown in the phase diagrams.
        By default, this is the reduced formula for the composition, but can be
        set to some other string for display purposes.
    """

    def __init__(self, ion, energy, name=None):
        self.energy = energy
        self.ion = ion
        self.composition = ion.composition
        self.name = name if name else self.ion.reduced_formula

    @classmethod
    def from_dict(cls, d):
        """
        Returns an IonEntry object from a dict.
        """
        return IonEntry(Ion.from_dict(d["ion"]), d["energy"], d.get("name", None))

    def as_dict(self):
        """
        Creates a dict of composition, energy, and ion name
        """
        d = {"ion": self.ion.as_dict(), "energy": self.energy,
             "name": self.name}
        return d

    def __repr__(self):
        return "IonEntry : {} with energy = {:.4f}".format(self.composition,
                                                           self.energy)

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


elements_HO = {Element('H'), Element('O')}


# TODO: There's a lot of functionality here that diverges
#   based on whether or not the pbx diagram is multielement
#   or not.  Could be a more elegant way to treat the two distinct modes,
#   for example by slicing Phase diagram at comp ratio

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

    Args:
        entries ([PourbaixEntry] or [MultiEntry]): Entries list
            containing Solids and Ions or a list of MultiEntries
        comp_dict ({str: float}): Dictionary of compositions,
            defaults to equal parts of each elements
        conc_dict ({str: float}): Dictionary of ion concentrations,
            defaults to 1e-6 for each element
        filter_solids (bool): applying this filter to a pourbaix
            diagram ensures all included phases are filtered by
            stability on the compositional phase diagram.  This
            breaks some of the functionality of the analysis,
            though, so use with caution.
        nproc (int): number of processes to generate multientries with
            in parallel.  Defaults to None (serial processing)
    """

    def __init__(self, entries, comp_dict=None, conc_dict=None,
                 filter_solids=False, nproc=None):
        entries = deepcopy(entries)

        # Get non-OH elements
        pbx_elts = set(itertools.chain.from_iterable(
            [entry.composition.elements for entry in entries]))
        pbx_elts = list(pbx_elts - elements_HO)

        # Process multientry inputs
        if isinstance(entries[0], MultiEntry):
            self._processed_entries = entries
            # Extract individual entries
            single_entries = list(set(itertools.chain.from_iterable(
                [e.entry_list for e in entries])))
            self._unprocessed_entries = single_entries
            self._filtered_entries = single_entries
            self._conc_dict = None
            self._elt_comp = {k: v for k, v in entries[0].composition.items()
                              if not k in elements_HO}
            self._multielement = True

        # Process single entry inputs
        else:
            # Set default conc/comp dicts
            if not comp_dict:
                comp_dict = {elt.symbol: 1. / len(pbx_elts) for elt in pbx_elts}
            if not conc_dict:
                conc_dict = {elt.symbol: 1e-6 for elt in pbx_elts}
            self._conc_dict = conc_dict

            self._elt_comp = comp_dict
            self.pourbaix_elements = pbx_elts

            solid_entries = [entry for entry in entries
                             if entry.phase_type == "Solid"]
            ion_entries = [entry for entry in entries
                           if entry.phase_type == "Ion"]

            # If a conc_dict is specified, override individual entry concentrations
            for entry in ion_entries:
                ion_elts = list(set(entry.composition.elements) - elements_HO)
                # TODO: the logic here for ion concentration setting is in two
                #       places, in PourbaixEntry and here, should be consolidated
                if len(ion_elts) == 1:
                    entry.concentration = conc_dict[ion_elts[0].symbol] \
                                          * entry.normalization_factor
                elif len(ion_elts) > 1 and not entry.concentration:
                    raise ValueError("Elemental concentration not compatible "
                                     "with multi-element ions")

            self._unprocessed_entries = solid_entries + ion_entries

            if not len(solid_entries + ion_entries) == len(entries):
                raise ValueError("All supplied entries must have a phase type of "
                                 "either \"Solid\" or \"Ion\"")

            if filter_solids:
                # O is 2.46 b/c pbx entry finds energies referenced to H2O
                entries_HO = [ComputedEntry('H', 0), ComputedEntry('O', 2.46)]
                solid_pd = PhaseDiagram(solid_entries + entries_HO)
                solid_entries = list(set(solid_pd.stable_entries) - set(entries_HO))

            self._filtered_entries = solid_entries + ion_entries

            if len(comp_dict) > 1:
                self._multielement = True
                self._processed_entries = self._generate_multielement_entries(
                    self._filtered_entries, nproc=nproc)
            else:
                self._processed_entries = self._filtered_entries
                self._multielement = False

        self._stable_domains, self._stable_domain_vertices = \
            self.get_pourbaix_domains(self._processed_entries)

    def _generate_multielement_entries(self, entries, forced_include=None,
                                       nproc=None):
        """
        Create entries for multi-element Pourbaix construction.

        This works by finding all possible linear combinations
        of entries that can result in the specified composition
        from the initialized comp_dict.

        Args:
            entries ([PourbaixEntries]): list of pourbaix entries
                to process into MultiEntries
            forced_include ([PourbaixEntries]) list of pourbaix entries
                that must be included in multielement entries
            nproc (int): number of processes to be used in parallel
                treatment of entry combos
        """
        N = len(self._elt_comp)  # No. of elements
        total_comp = Composition(self._elt_comp)
        forced_include = forced_include or []

        # generate all combinations of compounds that have all elements
        entry_combos = [itertools.combinations(
            entries, j + 1 - len(forced_include)) for j in range(N)]
        entry_combos = itertools.chain.from_iterable(entry_combos)
        if forced_include:
            entry_combos = [forced_include + list(ec) for ec in entry_combos]
        entry_combos = filter(lambda x: total_comp < MultiEntry(x).composition,
                              entry_combos)

        # Generate and filter entries
        processed_entries = []
        total = sum([comb(len(entries), j + 1 - len(forced_include))
                     for j in range(N)])
        if total > 1e6:
            warnings.warn("Your pourbaix diagram includes {} entries and may "
                          "take a long time to generate.".format(total))

        # Parallel processing of multi-entry generation
        if nproc is not None:
            f = partial(self.process_multientry, prod_comp=total_comp)
            with Pool(nproc) as p:
                processed_entries = list(tqdm(p.imap(f, entry_combos),
                                              total=total))
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
            coeffs = -np.array([rxn.get_coeff(comp) for comp in entry_comps])
            # Return None if reaction coeff threshold is not met
            # TODO: this filtration step might be put somewhere else
            if (coeffs > coeff_threshold).all():
                return MultiEntry(entry_list, weights=coeffs.tolist())
            else:
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
        hyperplanes = [np.array([-PREFAC * entry.npH, -entry.nPhi,
                                 0, -entry.energy]) * entry.normalization_factor
                       for entry in pourbaix_entries]
        hyperplanes = np.array(hyperplanes)
        hyperplanes[:, 2] = 1

        max_contribs = np.max(np.abs(hyperplanes), axis=0)
        g_max = np.dot(-max_contribs, [limits[0][1], limits[1][1], 0, 1])

        # Add border hyperplanes and generate HalfspaceIntersection
        border_hyperplanes = [[-1, 0, 0, limits[0][0]],
                              [1, 0, 0, -limits[0][1]],
                              [0, -1, 0, limits[1][0]],
                              [0, 1, 0, -limits[1][1]],
                              [0, 0, -1, 2 * g_max]]
        hs_hyperplanes = np.vstack([hyperplanes, border_hyperplanes])
        interior_point = np.average(limits, axis=1).tolist() + [g_max]
        hs_int = HalfspaceIntersection(hs_hyperplanes, np.array(interior_point))

        # organize the boundary points by entry
        pourbaix_domains = {entry: [] for entry in pourbaix_entries}
        for intersection, facet in zip(hs_int.intersections,
                                       hs_int.dual_facets):
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
            point_comparator = lambda x, y: x[0] * y[1] - x[1] * y[0]
            points_centered = sorted(points_centered,
                                     key=cmp_to_key(point_comparator))
            points = points_centered + center

            # Create simplices corresponding to pourbaix boundary
            simplices = [Simplex(points[indices])
                         for indices in ConvexHull(points).simplices]
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
        energies_at_conditions = [e.normalized_energy_at_conditions(pH, V)
                                  for e in self.stable_entries]
        return self.stable_entries[np.argmin(energies_at_conditions)]

    def get_decomposition_energy(self, entry, pH, V):
        """
        Finds decomposition to most stable entry

        Args:
            entry (PourbaixEntry): PourbaixEntry corresponding to
                compound to find the decomposition for
            pH (float): pH at which to find the decomposition
            V (float): voltage at which to find the decomposition

        Returns:
            reaction corresponding to the decomposition
        """
        # Find representative multientry
        if self._multielement and not isinstance(entry, MultiEntry):
            possible_entries = self._generate_multielement_entries(
                self._filtered_entries, forced_include=[entry])

            # Filter to only include materials where the entry is only solid
            if entry.phase_type == "solid":
                possible_entries = [e for e in possible_entries
                                    if e.phase_type.count("Solid") == 1]
            possible_energies = [e.normalized_energy_at_conditions(pH, V)
                                 for e in possible_entries]
        else:
            possible_energies = [entry.normalized_energy_at_conditions(pH, V)]

        min_energy = np.min(possible_energies, axis=0)

        # Find entry and take the difference
        hull = self.get_hull_energy(pH, V)
        return min_energy - hull

    def get_hull_energy(self, pH, V):
        all_gs = np.array([e.normalized_energy_at_conditions(
            pH, V) for e in self.stable_entries])
        base = np.min(all_gs, axis=0)
        return base

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

    def as_dict(self, include_unprocessed_entries=False):
        if include_unprocessed_entries:
            entries = [e.as_dict() for e in self._unprocessed_entries]
        else:
            entries = [e.as_dict() for e in self._processed_entries]
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "entries": entries,
             "comp_dict": self._elt_comp,
             "conc_dict": self._conc_dict}
        return d

    @classmethod
    def from_dict(cls, d):
        decoded_entries = MontyDecoder().process_decoded(d['entries'])
        return cls(decoded_entries, d.get('comp_dict'),
                   d.get('conc_dict'))


class PourbaixPlotter(object):
    """
    A plotter class for phase diagrams.

    Args:
        phasediagram: A PhaseDiagram object.
        show_unstable: Whether unstable phases will be plotted as well as
            red crosses. Defaults to False.
    """

    def __init__(self, pourbaix_diagram):
        self._pd = pourbaix_diagram

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

    def get_pourbaix_plot(self, limits=None, title="",
                          label_domains=True, plt=None):
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

        h_line = np.transpose([[xlim[0], -xlim[0] * PREFAC],
                               [xlim[1], -xlim[1] * PREFAC]])
        o_line = np.transpose([[xlim[0], -xlim[0] * PREFAC + 1.23],
                               [xlim[1], -xlim[1] * PREFAC + 1.23]])
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

        for entry, vertices in self._pd._stable_domain_vertices.items():
            center = np.average(vertices, axis=0)
            x, y = np.transpose(np.vstack([vertices, vertices[0]]))
            plt.plot(x, y, 'k-', linewidth=lw)
            if label_domains:
                plt.annotate(generate_entry_label(entry), center, ha='center',
                             va='center', fontsize=20, color="b")

        plt.xlabel("pH")
        plt.ylabel("E (V)")
        plt.title(title, fontsize=20, fontweight='bold')
        return plt

    def plot_entry_stability(self, entry, pH_range=None, pH_resolution=100,
                             V_range=None, V_resolution=100, e_hull_max=1,
                             cmap='RdYlBu_r', **kwargs):
        if pH_range is None:
            pH_range = [-2, 16]
        if V_range is None:
            V_range = [-3, 3]
        # plot the Pourbaix diagram
        plt = self.get_pourbaix_plot(**kwargs)
        pH, V = np.mgrid[pH_range[0]:pH_range[1]:pH_resolution * 1j,
                V_range[0]:V_range[1]:V_resolution * 1j]

        stability = self._pd.get_decomposition_energy(entry, pH, V)
        # Plot stability map
        plt.pcolor(pH, V, stability, cmap=cmap, vmin=0, vmax=e_hull_max)
        cbar = plt.colorbar()
        cbar.set_label("Stability of {} (eV/atom)".format(
            generate_entry_label(entry)))

        # Set ticklabels
        ticklabels = [t.get_text() for t in cbar.ax.get_yticklabels()]
        ticklabels[-1] = '>={}'.format(ticklabels[-1])
        cbar.ax.set_yticklabels(ticklabels)

        return plt

    def domain_vertices(self, entry):
        """
        Returns the vertices of the Pourbaix domain.

        Args:
            entry: Entry for which domain vertices are desired

        Returns:
            list of vertices
        """
        return self._pd._pourbaix_domain_vertices[entry]


def generate_entry_label(entry):
    """
    Generates a label for the pourbaix plotter

    Args:
        entry (PourbaixEntry or MultiEntry): entry to get a label for
    """
    if isinstance(entry, MultiEntry):
        return " + ".join([latexify_ion(e.name) for e in entry.entry_list])
    else:
        return latexify_ion(latexify(entry.name))


def latexify_ion(formula):
    return re.sub(r"()\[([^)]*)\]", r"\1$^{\2}$", formula)
