# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import logging
import numpy as np
import itertools

from scipy.spatial import ConvexHull
from pymatgen.analysis.pourbaix.entry import MultiEntry, ion_or_solid_comp_object
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
                 filter_multielement=False):
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
        self._make_pourbaix_diagram()

    def _create_conv_hull_data(self):
        """
        Make data conducive to convex hull generator.
        """
        entries_to_process = list()
        for entry in self._processed_entries:
            entry.scale(entry.normalization_factor)
            entry.correction += (- MU_H2O * entry.nH2O + entry.conc_term)
            entries_to_process.append(entry)
        self._qhull_entries = entries_to_process
        return self._process_conv_hull_data(entries_to_process)

    def _process_conv_hull_data(self, entries_to_process):
        """
        From a sequence of ion+solid entries, generate the necessary data
        for generation of the convex hull.
        """
        data = []
        for entry in entries_to_process:
            row = [entry.npH, entry.nPhi, entry.g0]
            data.append(row)
        temp = sorted(zip(data, self._qhull_entries),
                      key=lambda x: x[0][2])
        [data, self._qhull_entries] = list(zip(*temp))
        return data

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

    def _make_pourbaix_diagram(self):
        """
        Calculates entries on the convex hull in the dual space.
        """
        stable_entries = set()
        self._qhull_data = self._create_conv_hull_data()
        dim = len(self._qhull_data[0])
        if len(self._qhull_data) < dim:
            # TODO: might want to lift this restriction and
            # supply a warning instead, should work even if it's slow.
            raise NotImplementedError("Can only do elements with at-least "
                                      "3 entries for now")
        if len(self._qhull_data) == dim:
            self._facets = [list(range(dim))]
        else:
            facets_hull = np.array(ConvexHull(self._qhull_data).simplices)
            self._facets = np.sort(np.array(facets_hull))
            logger.debug("Final facets are\n{}".format(self._facets))

            logger.debug("Removing vertical facets...")
            vert_facets_removed = list()
            for facet in self._facets:
                facetmatrix = np.zeros((len(facet), len(facet)))
                count = 0
                for vertex in facet:
                    facetmatrix[count] = np.array(self._qhull_data[vertex])
                    facetmatrix[count, dim - 1] = 1
                    count += 1
                if abs(np.linalg.det(facetmatrix)) > 1e-8:
                    vert_facets_removed.append(facet)
                else:
                    logger.debug("Removing vertical facet : {}".format(facet))

            logger.debug("Removing UCH facets by eliminating normal.z >0 ...")

            # Find center of hull
            vertices = set()
            for facet in vert_facets_removed:
                for vertex in facet:
                    vertices.add(vertex)
            c = [0.0, 0.0, 0.0]
            c[0] = np.average([self._qhull_data[vertex][0]
                               for vertex in vertices])
            c[1] = np.average([self._qhull_data[vertex][1]
                               for vertex in vertices])
            c[2] = np.average([self._qhull_data[vertex][2]
                               for vertex in vertices])

            # Shift origin to c
            new_qhull_data = np.array(self._qhull_data)
            for vertex in vertices:
                new_qhull_data[vertex] -= c

            # For each facet, find normal n, find dot product with P, and
            # check if this is -ve
            final_facets = list()
            for facet in vert_facets_removed:
                a = new_qhull_data[facet[1]] - new_qhull_data[facet[0]]
                b = new_qhull_data[facet[2]] - new_qhull_data[facet[0]]
                n = np.cross(a, b)
                val = np.dot(n, new_qhull_data[facet[0]])
                if val < 0:
                    n = -n
                if n[2] <= 0:
                    final_facets.append(facet)
                else:
                    logger.debug("Removing UCH facet : {}".format(facet))
            final_facets = np.array(final_facets)
            self._facets = final_facets

        stable_vertices = set()
        for facet in self._facets:
            for vertex in facet:
                stable_vertices.add(vertex)
                stable_entries.add(self._qhull_entries[vertex])
        self._stable_entries = stable_entries
        self._vertices = stable_vertices

    @property
    def facets(self):
        """
        Facets of the convex hull in the form of  [[1,2,3],[4,5,6]...]
        """
        return self._facets

    @property
    def qhull_data(self):
        """
        Data used in the convex hull operation. This is essentially a matrix of
        composition data and energy per atom values created from qhull_entries.
        """
        return self._qhull_data

    @property
    def qhull_entries(self):
        """
        Return qhull entries
        """
        return self._qhull_entries

    @property
    def stable_entries(self):
        """
        Returns the stable entries in the Pourbaix diagram.
        """
        return list(self._stable_entries)
    
    @property
    def unstable_entries(self):
        """
        Returns all unstable entries in the Pourbaix diagram
        """
        return [e for e in self.qhull_entries if e not in self.stable_entries]

    @property
    def all_entries(self):
        """
        Return all entries used to generate the pourbaix diagram
        """
        return self._processed_entries

    @property
    def vertices(self):
        """
        Return vertices of the convex hull
        """
        return self._vertices

    @property
    def unprocessed_entries(self):
        """
        Return unprocessed entries
        """
        return self._unprocessed_entries
