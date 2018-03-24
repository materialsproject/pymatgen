# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import logging
import numpy as np
import itertools
from copy import deepcopy
from functools import cmp_to_key

from scipy.spatial import ConvexHull, HalfspaceIntersection
from pymatgen.util.coord import Simplex
from pymatgen.analysis.pourbaix.entry import MultiEntry
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError
from pymatgen.analysis.phase_diagram import PhaseDiagram

"""
Module containing analysis classes which compute a pourbaix diagram given a
target compound/element.
"""

from six.moves import zip

__author__ = "Sai Jayaraman"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.0"
__maintainer__ = "Sai Jayaraman"
__credits__ = "Arunima Singh, Joseph Montoya"
__email__ = "sjayaram@mit.edu"
__status__ = "Development"
__date__ = "Nov 1, 2012"


logger = logging.getLogger(__name__)

PREFAC = 0.0591
MU_H2O = -2.4583
elements_HO = {Element('H'), Element('O')}
# TODO: There's a lot of functionality here that diverges
#   based on whether or not the pbx diagram is multielement
#   or not.  Could be a more elegant way to
#   treat the two distinct modes.

# TODO: the solids filter breaks some of the functionality of the
#       heatmap plotter, because the reference states for decomposition
#       don't include oxygen/hydrogen in the OER/HER regions

class PourbaixDiagram(object):
    """
    Class to create a Pourbaix diagram from entries

    Args:
        entries [Entry]: Entries list containing both Solids and Ions
        comp_dict {str: float}: Dictionary of compositions, defaults to
            equal parts of each elements
        conc_dict {str: float}: Dictionary of ion concentrations, defaults
            to 1e-6 for each element
        filter_multielement (bool): applying this filter to a multi-
            element pourbaix diagram makes generates it a bit more
            efficiently by filtering the entries used to generate
            the hull.  This breaks some of the functionality of
            the analyzer, though, so use with caution.
    """
    def __init__(self, entries, comp_dict=None, conc_dict=None,
                 filter_solids=True):

        entries = deepcopy(entries)
        # Get non-OH elements
        pbx_elts = set(itertools.chain.from_iterable(
            [entry.composition.elements for entry in entries]))
        pbx_elts = list(pbx_elts - elements_HO)

        # Set default conc/comp dicts
        if not comp_dict:
            comp_dict = {elt.symbol : 1. / len(pbx_elts) for elt in pbx_elts}
        if not conc_dict:
            conc_dict = {elt.symbol : 1e-6 for elt in pbx_elts}

        self._elt_comp = comp_dict
        self.pourbaix_elements = pbx_elts

        solid_entries = [entry for entry in entries
                         if entry.phase_type == "Solid"]
        ion_entries = [entry for entry in entries
                       if entry.phase_type == "Ion"]

        for entry in ion_entries:
            ion_elts = list(set(entry.composition.elements) - elements_HO)
            if len(ion_elts) != 1:
                raise ValueError("Elemental concentration not compatible "
                                 "with multi-element ions")
            entry.conc = conc_dict[ion_elts[0].symbol]

        if not len(solid_entries + ion_entries) == len(entries):
            raise ValueError("All supplied entries must have a phase type of "
                             "either \"Solid\" or \"Ion\"")

        self._unprocessed_entries = entries

        if len(comp_dict) > 1:
            self._multielement = True
            if filter_multielement:
                # Add two high-energy H/O entries that ensure the hull
                # includes all stable solids.
                entries_HO = [ComputedEntry('H', 10000), ComputedEntry('O', 10000)]
                solid_pd = PhaseDiagram(solid_entries + entries_HO)
                solid_entries = list(set(solid_pd.stable_entries) - set(entries_HO))
            self._processed_entries = self._generate_multielement_entries(
                    solid_entries + ion_entries)
        else:
            self._multielement = False

            self._processed_entries = solid_entries + ion_entries
        self._stable_domains, self._stable_domain_vertices = \
            self.get_pourbaix_domains(self._processed_entries)


    def _generate_multielement_entries(self, entries):
        """
        Create entries for multi-element Pourbaix construction.

        This works by finding all possible linear combinations
        of entries that can result in the specified composition
        from the initialized comp_dict.

        Args:
            entries ([PourbaixEntries]): list of pourbaix entries
                to process into MultiEntries
        """
        N = len(self._elt_comp)  # No. of elements
        total_comp = Composition(self._elt_comp)

        # generate all possible combinations of compounds that have all elts
        entry_combos = [itertools.combinations(entries, j+1) for j in range(N)]
        entry_combos = itertools.chain.from_iterable(entry_combos)
        entry_combos = filter(lambda x: total_comp < MultiEntry(x).total_composition,
                              entry_combos)

        # Generate and filter entries
        processed_entries = []
        for entry_combo in entry_combos:
            processed_entry = self.process_multientry(entry_combo, total_comp)
            if processed_entry is not None:
                processed_entries.append(processed_entry)

        return processed_entries

    @staticmethod
    def process_multientry(entry_list, prod_comp):
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
            comp (Composition): composition constraint for setting
                weights of MultiEntry
        """
        dummy_oh = [Composition("H"), Composition("O")]
        try:
            # Get balanced reaction coeffs, ensuring all < 0 or conc thresh
            # Note that we get reduced compositions for solids and non-reduced
            # compositions for ions because ions aren't normalized due to
            # their charge state.
            entry_comps = [e.composition if e.phase_type=='Ion'
                           else e.composition.reduced_composition
                           for e in entry_list]
            rxn = Reaction(entry_comps + dummy_oh, [prod_comp])
            thresh = np.array([pe.conc if pe.phase_type == "Ion"
                               else 1e-3 for pe in entry_list])
            coeffs = -np.array([rxn.get_coeff(comp) for comp in entry_comps])
            if (coeffs > thresh).all():
                weights = coeffs / coeffs[0]
                return MultiEntry(entry_list, weights=weights.tolist())
            else:
                return None
        except ReactionError:
            return None

    @staticmethod
    def get_pourbaix_domains(pourbaix_entries, limits=[[-2, 16], [-4, 4]]):
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
        # Get hyperplanes
        hyperplanes = [[entry.npH, entry.nPhi, 1, -entry.g0]
                       for entry in pourbaix_entries]

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
            point_comparator = lambda x, y: x[0]*y[1] - x[1]*y[0]
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
        for entry, simplex in self._stable_domains.items():
            if simplex.in_simplex(pH, V):
                return entry

    def get_decomposition(self, entry, pH, V):
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
               self._preprocessed_entries, forced_include=[single_entry])
            # Filter to only include materials where the entry is only solid
            possible_entries = [e for e in possible_entries
                                if e.phases.count("Solid") == 1]
            entry = sorted(possible_entries, key=lambda x: x.g(pH, V))[0]

        # Find stable entry and take the difference
        stable_entry = self.find_stable_entry(pH, V)
        return entry.g(pH, V) - stable_entry.g(pH, V)

    def get_stability_map(self, entry, pH_range=[-2, 16], pH_resolution=100,
                          V_range=[-3, 3], V_resolution=100):
        """
        Finds hull energies

        Args:
            entry (PourbaixEntry): entry for which to find pourbaix stability
            pH_range (list): range of pHs for which to find stability
            pH_resolution (int): resolution for pH
            V_range (list): range of Vs for which to find stability
            V_resolution (int): resolution for V

        Returns:
            pHs, Vs, and pourbaix hull energies for the specified entry

        """
        # Get pourbaix "base" from all stable entries

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
        return [e for e in self._stable_domains.keys()
                if e not in self.stable_entries]

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
