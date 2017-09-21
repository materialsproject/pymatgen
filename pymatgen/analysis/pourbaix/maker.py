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
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError

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

# TODO: There's a lot of functionality here that diverges
#   based on whether or not the pbx diagram is multielement
#   or not.  Could be a more elegant way to
#   treat the two distinct modes.

class PourbaixDiagram(object):
    """
    Class to create a Pourbaix diagram from entries
    Args:
        entries: Entries list containing both Solids and Ions
        comp_dict: Dictionary of compositions
    """
    def __init__(self, entries, comp_dict=None):
        self._solid_entries = list()
        self._ion_entries = list()
        for entry in entries:
            if entry.phase_type == "Solid":
                self._solid_entries.append(entry)
            elif entry.phase_type == "Ion":
                self._ion_entries.append(entry)
            else:
                raise Exception("Incorrect Phase type - needs to be \
                Pourbaix entry of phase type Ion/Solid")
        if len(self._ion_entries) == 0:
            raise Exception("No ion phase. Equilibrium between ion/solid "
                            "is required to make a Pourbaix Diagram")
        self._unprocessed_entries = self._solid_entries + self._ion_entries
        self._elt_comp = comp_dict

        if comp_dict and len(comp_dict) > 1:
            self._multielement = True
            pbx_elements = set()
            for comp in comp_dict.keys():
                for el in [el for el in
                           ion_or_solid_comp_object(comp).elements
                           if el not in ["H", "O"]]:
                    pbx_elements.add(el.symbol)
            self.pourbaix_elements = pbx_elements
            w = [comp_dict[key] for key in comp_dict]
            A = []
            for comp in comp_dict:
                comp_obj = ion_or_solid_comp_object(comp)
                Ai = []
                for elt in self.pourbaix_elements:
                    Ai.append(comp_obj[elt])
                A.append(Ai)
            A = np.array(A).T.astype(float)
            w = np.array(w)
            A /= np.dot([a.sum() for a in A], w)
            x = np.linalg.solve(A, w)
            self._elt_comp = dict(zip(self.pourbaix_elements, x))

        else:
            self._multielement = False
            self.pourbaix_elements = [el.symbol
                                      for el in entries[0].composition.elements
                                      if el.symbol not in ["H", "O"]]
            # TODO: document the physical meaning
            # of ternary oxide of pbx diagrams in both modes
            self._elt_comp = {self.pourbaix_elements[0]: 1.0}
        self._make_pourbaix_diagram()

    def _create_conv_hull_data(self):
        """
        Make data conducive to convex hull generator.
        """
        if self._multielement:
            self._all_entries = self._process_multielement_entries()
        else:
            self._all_entries = self._unprocessed_entries
        entries_to_process = list()
        for entry in self._all_entries:
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

    def _process_multielement_entries(self):
        """
        Create entries for multi-element Pourbaix construction.

        This works by finding all possible linear combinations
        of entries that can result in the specified composition
        from the initialized comp_dict.
        """
        N = len(self._elt_comp)  # No. of elements
        entries = self._unprocessed_entries
        dummy_prod = Composition(self._elt_comp) 
        total_comp = Composition(self._elt_comp)

        # generate all possible combinations of compounds that have all elts
        entry_combos = [itertools.combinations(entries, j+1) for j in range(N)]
        entry_combos = itertools.chain.from_iterable(entry_combos)
        entry_combos = filter(lambda x: dummy_prod < MultiEntry(x).total_composition, 
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
        Return all entries
        """
        return self._all_entries

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
