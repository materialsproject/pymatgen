#!/usr/bin/env python

"""
This module provides classes for analyzing phase diagrams.
"""

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import numpy as np
from pymatgen.core.structure import Composition
from pymatgen.phasediagram.pdmaker import PhaseDiagram, GrandPotentialPhaseDiagram
from pymatgen.analysis.reaction_calculator import Reaction

class PDAnalyzer(object):
    """
    A class for performing analyses on Phase Diagrams.
    """
    
    numerical_tol = 1e-8
    
    def __init__(self, pd):
        """
        Arguments:
            pd - Phase Diagram to analyze.
        """
        self._pd = pd

    def _make_comp_matrix(self,complist):
        """
        Helper function to generates a normalized composition matrix from a list of composition.
        """
        return np.array([[comp.get_atomic_fraction(el) for el in self._pd.elements] for comp in complist])

    def _in_facet(self, facet, comp):
        """
        Checks if a composition is in a facet using the standard barycentric coordinate system algorithm.
        The general concept is that the facets from a convex hull form a simplex in n-dimensional space.
        Taking an arbitrary vertex as an origin, we compute the basis for the simplex from this origin by
        subtracting all other vertices from the origin. We then project the composition into this 
        coordinate system and determine the coefficients in this coordinate system.  If the coeffs satisfy
        that all coeffs >= 0 and sum(coeffs) <= 1, the composition is in the facet.
        For example, take a 4-comp PD where the facets are tetrahedrons. For a tetrahedron, let's label
        the vertices as O, A, B anc C.  Let's call our comp coordinate as X.
        We form the composition matrix M with vectors OA, OB and OB, transponse it, and solve for 
            M'.a = OX
        where a are the coefficients.
        if (a >= 0).all() and sum(a) <= 1, X is in the tetrahedron.
        Note that in reality, the test needs to provide a tolerance (set to 1e-8 by default) for numerical errors.
        """
        dim = len(self._pd.elements)
        if dim > 1:
            origin = np.array(self._pd.qhull_data[facet[0]][0:dim-1])
            cm = np.array([np.array(self._pd.qhull_data[facet[i]][0:dim-1]) - origin for i in xrange(1, len(facet))])
            row = [comp.get_atomic_fraction(self._pd.elements[i]) for i in xrange(1,len(self._pd.elements))]
            compm = np.array(row) - origin
            coeffs = np.linalg.solve(cm.transpose(), compm)
            return (coeffs >= -PDAnalyzer.numerical_tol).all() and sum(coeffs) <= (1 + PDAnalyzer.numerical_tol)
        else:
            return True

    def _get_facets(self,comp):
        """
        Get the facets that a composition falls into.
        """
        memberfacets = list()
        for facet in self._pd.facets:
            if self._in_facet(facet, comp):
                memberfacets.append(facet)
        return memberfacets
    
    def _get_facet(self,comp):
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
        Arguments:
            comp - A composition
        Returns:
            Decomposition as a dict of {PDEntry: amount}
        """
        facet = self._get_facet(comp)
        complist = [self._pd.qhull_entries[i].composition for i in facet]
        m = self._make_comp_matrix(complist)
        compm = self._make_comp_matrix([comp])
        decompamts = np.dot(np.linalg.inv(m.transpose()),compm.transpose())
        decomp = dict()
        #Scrub away zero amounts
        for i in xrange(len(decompamts)):
            if abs(decompamts[i][0]) > PDAnalyzer.numerical_tol:
                decomp[self._pd.qhull_entries[facet[i]]] = decompamts[i][0]
        return decomp
    
    def get_decomp_and_e_above_hull(self, entry):
        """
        Provides the decomposition and energy above convex hull for an entry
        Arguments:
            entry - A PDEntry like object
        Returns:
            (decomp, energy above convex hull)  Stable entries should have energy above hull of 0.
        """
        comp = entry.composition
        eperatom = entry.energy_per_atom
        decomp = self.get_decomposition(comp)
        hullenergy = sum([entry.energy_per_atom*amt for entry, amt in decomp.items()])
        if abs(eperatom) < PDAnalyzer.numerical_tol:
            return (decomp, 0)
        return (decomp, eperatom - hullenergy)

    def get_e_above_hull(self, entry):
        """
        Provides the energy above convex hull for an entry
        Arguments:
            entry - A PDEntry like object
        Returns:
            Energy above convex hull of entry.  Stable entries should have energy above hull of 0.
        """
        return self.get_decomp_and_e_above_hull(entry)[1]


    def get_equilibrium_reaction_energy(self, entry):
        """
        Provides the reaction energy of a stable entry from the neighboring equilibrium stable entries.
        (also known as the inverse distance to hull).
        Arguments:
            entry - A PDEntry like object
        Returns:
            Equilibrium reaction energy of entry.  Stable entries should have equilibrium reaction energy <= 0.
        """
        if entry not in self._pd.stable_entries:
            raise ValueError("Equilibrium reaction energy is available only for stable entries.")
        entries = [e for e in self._pd.all_entries if e != entry]
        modpd = PhaseDiagram(entries, self._pd.elements)
        analyzer = PDAnalyzer(modpd)
        return analyzer.get_decomp_and_e_above_hull(entry)[1]
    
    def get_transition_chempots(self, element):
        """
        Get the critical chemical potentials for an element in the Phase Diagram.
        Arguments:
            element - An element.  Has to be in the PD in the first place.
        Returns:
            A sorted sequence of critical chemical potentials, from less negative to more negative.
        """
        if element not in self._pd.elements:
            raise ValueError("get_transition_chempots can only be called with elements in the phase diagram.")
        ind = self._pd.elements.index(element)
        critical_chempots = []
        for facet in self._pd.facets:
            complist = [self._pd.qhull_entries[i].composition for i in facet]
            energylist = [self._pd.qhull_entries[i].energy_per_atom for i in facet]
            m = self._make_comp_matrix(complist)
            chempots = np.dot(np.linalg.inv(m),energylist)
            critical_chempots.append(chempots[ind])
                
        clean_pots = []
        for c in sorted(critical_chempots):
            if len(clean_pots) == 0:
                clean_pots.append(c)
            else:
                if abs(c-clean_pots[-1]) > PDAnalyzer.numerical_tol:
                    clean_pots.append(c)
        clean_pots.reverse()
        return tuple(clean_pots)
    
    def get_element_profile(self, element, comp):
        """
        Provides the element evolution data for a composition.
        For example, can be used to analyze Li conversion voltages by varying uLi and looking at the phases formed.
        Also can be used to analyze O2 evolution by varying uO2.
        Arguments:
            element - An element. Must be in the phase diagram.
            comp - a Composition
        Returns:
            Evolution data as a list of dictionaries of the following format: [ {'chempot': -10.487582010000001, 'evolution': -2.0, 'reaction': Reaction Object], ...]
        """
        if element not in self._pd.elements:
            raise ValueError("get_transition_chempots can only be called with elements in the phase diagram.")
        chempots = self.get_transition_chempots(element)
        stable_entries = self._pd.stable_entries
        gccomp = Composition({el:amt for el, amt in comp.items() if el != element})
        elref = self._pd.el_refs[element]
        elcomp = Composition.from_formula(element.symbol)
        prev_decomp = [];
        evolution = []
        def are_same_decomp(decomp1, decomp2):
            for comp in decomp2:
                if comp not in decomp1:
                    return False
            return True
        
        for c in chempots:
            gcpd = GrandPotentialPhaseDiagram(stable_entries, {element:c-0.01}, self._pd.elements)
            analyzer = PDAnalyzer(gcpd)
            decomp = [gcentry.original_entry.composition for gcentry, amt in analyzer.get_decomposition(gccomp).items() if amt > 1e-5]
            decomp_entries = [gcentry.original_entry for gcentry, amt in analyzer.get_decomposition(gccomp).items() if amt > 1e-5]

            if not are_same_decomp(prev_decomp, decomp):
                if elcomp not in decomp:
                    decomp.insert(0,elcomp)
                rxn = Reaction([comp], decomp)
                rxn.normalize_to(comp)
                prev_decomp = decomp
                evolution.append({'chempot':c, 'evolution' : - rxn.coeffs[rxn.all_comp.index(elcomp)], 'element_reference': elref, 'reaction':rxn, 'entries':decomp_entries})
            
        return evolution
