#!/usr/bin/env python

"""
This module provides classes for analyzing phase diagrams.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "May 16, 2012"

import numpy as np
import itertools

from pymatgen.core.structure import Composition
from pymatgen.phasediagram.pdmaker import PhaseDiagram, \
    GrandPotentialPhaseDiagram
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.comp_geometry.simplex import Simplex


class PDAnalyzer(object):
    """
    A class for performing analyses on Phase Diagrams.

    The algorithm is based on the work in the following papers:

    1. S. P. Ong, L. Wang, B. Kang, and G. Ceder, Li-Fe-P-O2 Phase Diagram from
       First Principles Calculations. Chem. Mater., 2008, 20(5), 1798-1807.
       doi:10.1021/cm702327g

    2. S. P. Ong, A. Jain, G. Hautier, B. Kang, G. Ceder, Thermal stabilities
       of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first
       principles calculations. Electrochem. Comm., 2010, 12(3), 427-430.
       doi:10.1016/j.elecom.2010.01.010
    """

    numerical_tol = 1e-8

    def __init__(self, pd):
        """
        Args:
            pd - Phase Diagram to analyze.
        """
        self._pd = pd

    def _make_comp_matrix(self, complist):
        """
        Helper function to generates a normalized composition matrix from a
        list of compositions.
        """
        return np.array([[comp.get_atomic_fraction(el)
                          for el in self._pd.elements] for comp in complist])

    def _in_facet(self, facet, comp):
        """
        Checks if a composition is in a facet.

        Args:
            facet:
                facet to test.
            comp:
                Composition to test.
        """
        dim = len(self._pd.elements)
        if dim > 1:
            coords = [np.array(self._pd.qhull_data[facet[i]][0:dim - 1])
                      for i in xrange(len(facet))]
            simplex = Simplex(coords)
            comp_point = [comp.get_atomic_fraction(self._pd.elements[i])
                          for i in xrange(1, len(self._pd.elements))]
            return simplex.in_simplex(comp_point, PDAnalyzer.numerical_tol)
        else:
            return True

    def _get_facets(self, comp):
        """
        Get the facets that a composition falls into.
        """
        memberfacets = list()
        for facet in self._pd.facets:
            if self._in_facet(facet, comp):
                memberfacets.append(facet)
        return memberfacets

    def _get_facet(self, comp):
        """
        Get the facets that a composition falls into.
        """
        for facet in self._pd.facets:
            if self._in_facet(facet, comp):
                return facet
        raise RuntimeError("No facet found for comp = {}".format(comp))

    def get_decomposition(self, comp):
        """
        Provides the decomposition at a particular composition

        Args:
            comp:
                A composition

        Returns:
            Decomposition as a dict of {PDEntry: amount}
        """
        facet = self._get_facet(comp)
        complist = [self._pd.qhull_entries[i].composition for i in facet]
        m = self._make_comp_matrix(complist)
        compm = self._make_comp_matrix([comp])
        decompamts = np.dot(np.linalg.inv(m.transpose()), compm.transpose())
        decomp = dict()
        #Scrub away zero amounts
        for i in xrange(len(decompamts)):
            if abs(decompamts[i][0]) > PDAnalyzer.numerical_tol:
                decomp[self._pd.qhull_entries[facet[i]]] = decompamts[i][0]
        return decomp

    def get_decomp_and_e_above_hull(self, entry):
        """
        Provides the decomposition and energy above convex hull for an entry

        Args:
            entry:
                A PDEntry like object

        Returns:
            (decomp, energy above convex hull)  Stable entries should have
            energy above hull of 0.
        """
        comp = entry.composition
        eperatom = entry.energy_per_atom
        decomp = self.get_decomposition(comp)
        hullenergy = sum([entry.energy_per_atom * amt
                          for entry, amt in decomp.items()])
        if abs(eperatom) < PDAnalyzer.numerical_tol:
            return (decomp, 0)
        return (decomp, eperatom - hullenergy)

    def get_e_above_hull(self, entry):
        """
        Provides the energy above convex hull for an entry

        Args:
            entry - A PDEntry like object

        Returns:
            Energy above convex hull of entry. Stable entries should have
            energy above hull of 0.
        """
        return self.get_decomp_and_e_above_hull(entry)[1]

    def get_equilibrium_reaction_energy(self, entry):
        """
        Provides the reaction energy of a stable entry from the neighboring
        equilibrium stable entries (also known as the inverse distance to
        hull).

        Args:
            entry:
                A PDEntry like object

        Returns:
            Equilibrium reaction energy of entry. Stable entries should have
            equilibrium reaction energy <= 0.
        """
        if entry not in self._pd.stable_entries:
            raise ValueError("Equilibrium reaction energy is available only "
                             "for stable entries.")
        if entry.is_element:
            return 0
        entries = [e for e in self._pd.stable_entries if e != entry]
        modpd = PhaseDiagram(entries, self._pd.elements)
        analyzer = PDAnalyzer(modpd)
        return analyzer.get_decomp_and_e_above_hull(entry)[1]

    def get_facet_chempots(self, facet):
        """
        Calculates the chemical potentials for each element within a facet.

        Args:
            facet:
                Facet of the phase diagram.

        Returns:
            { element: chempot } for all elements in the phase diagram.
        """
        complist = [self._pd.qhull_entries[i].composition for i in facet]
        energylist = [self._pd.qhull_entries[i].energy_per_atom for i in facet]
        m = self._make_comp_matrix(complist)
        chempots = np.dot(np.linalg.inv(m), energylist)
        return dict(zip(self._pd.elements, chempots))

    def get_transition_chempots(self, element):
        """
        Get the critical chemical potentials for an element in the Phase
        Diagram.

        Args:
            element:
                An element. Has to be in the PD in the first place.

        Returns:
            A sorted sequence of critical chemical potentials, from less
            negative to more negative.
        """
        if element not in self._pd.elements:
            raise ValueError("get_transition_chempots can only be called with "
                             "elements in the phase diagram.")

        critical_chempots = []
        for facet in self._pd.facets:
            chempots = self.get_facet_chempots(facet)
            critical_chempots.append(chempots[element])

        clean_pots = []
        for c in sorted(critical_chempots):
            if len(clean_pots) == 0:
                clean_pots.append(c)
            else:
                if abs(c - clean_pots[-1]) > PDAnalyzer.numerical_tol:
                    clean_pots.append(c)
        clean_pots.reverse()
        return tuple(clean_pots)

    def get_element_profile(self, element, comp, comp_tol=1e-5):
        """
        Provides the element evolution data for a composition.
        For example, can be used to analyze Li conversion voltages by varying
        uLi and looking at the phases formed. Also can be used to analyze O2
        evolution by varying uO2.

        Args:
            element:
                An element. Must be in the phase diagram.
            comp:
                A Composition
            comp_tol:
                The tolerance to use when calculating decompositions. Phases
                with amounts less than this tolerance are excluded. Defaults to
                1e-5.

        Returns:
            Evolution data as a list of dictionaries of the following format:
            [ {'chempot': -10.487582010000001, 'evolution': -2.0,
            'reaction': Reaction Object], ...]
        """
        if element not in self._pd.elements:
            raise ValueError("get_transition_chempots can only be called with"
                             " elements in the phase diagram.")
        chempots = self.get_transition_chempots(element)
        stable_entries = self._pd.stable_entries
        gccomp = Composition({el: amt for el, amt in comp.items()
                              if el != element})
        elref = self._pd.el_refs[element]
        elcomp = Composition.from_formula(element.symbol)
        prev_decomp = []
        evolution = []

        def are_same_decomp(decomp1, decomp2):
            for comp in decomp2:
                if comp not in decomp1:
                    return False
            return True

        for c in chempots:
            gcpd = GrandPotentialPhaseDiagram(stable_entries,
                                              {element: c - 0.01},
                                              self._pd.elements)
            analyzer = PDAnalyzer(gcpd)
            decomp = [gcentry.original_entry.composition for gcentry,
                      amt in analyzer.get_decomposition(gccomp).items()
                      if amt > comp_tol]
            decomp_entries = [gcentry.original_entry for gcentry,
                              amt in analyzer.get_decomposition(gccomp).items()
                              if amt > comp_tol]

            if not are_same_decomp(prev_decomp, decomp):
                if elcomp not in decomp:
                    decomp.insert(0, elcomp)
                rxn = Reaction([comp], decomp)
                rxn.normalize_to(comp)
                prev_decomp = decomp
                amt = -rxn.coeffs[rxn.all_comp.index(elcomp)]
                evolution.append({'chempot': c,
                                  'evolution': amt,
                                  'element_reference': elref,
                                  'reaction': rxn, 'entries': decomp_entries})

        return evolution

    def get_chempot_range_map(self, elements):
        """
        Returns a chemical potential range map for each stable entry.

        Args:
            elements:
                Sequence of elements to be considered as independent variables.
                E.g., if you want to show the stability ranges of all Li-Co-O
                phases wrt to uLi and uO, you will supply
                [Element("Li"), Element("O")]

        Returns:
            Returns a dict of the form {entry: [simplices]}. The list of
            simplices are the sides of the N-1 dim polytope bounding the
            allowable chemical potential range of each entry.
        """
        elrefs = self._pd.el_refs
        chempot_ranges = {}
        for entry in self._pd.stable_entries:
            all_facets = self._get_facets(entry.composition)
            simplices = []
            # For each entry, go through all possible combinations of 2 facets.
            for facets in itertools.combinations(all_facets, 2):
                # Get the intersection of the 2 facets.
                inter = set(facets[0]).intersection(set(facets[1]))

                #Check if the intersection has N-1 vertices. if so, add the
                #line to the list of simplices.
                if len(inter) == self._pd.dim - 1:
                    coords = []
                    for facet in facets:
                        chempots = self.get_facet_chempots(facet)
                        coords.append([chempots[el]
                                       - elrefs[el].energy_per_atom
                                       for el in elements])
                    sim = Simplex(coords)
                    simplices.append(sim)

            if len(simplices) > 0:
                chempot_ranges[entry] = simplices

        return chempot_ranges
