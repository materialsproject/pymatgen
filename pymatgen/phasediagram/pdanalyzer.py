#!/usr/bin/env python

from __future__ import division

"""
This module provides classes for analyzing phase diagrams.
"""

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import numpy as np
from pymatgen.core.structure import Composition
from pymatgen.phasediagram.pdmaker import GrandPotentialPhaseDiagram
from pymatgen.analysis.reaction_calculator import Reaction

class PDAnalyzer(object):
    """
    A class for performing analyses on Phase Diagrams.
    """

    def __init__(self, pd):
        """
        Arguments:
            pd - Phase Diagram to analyze.
        """
        self._pd = pd

    def _make_comp_matrix(self,complist):
        m = list()
        for comp in complist:
            sumform = comp.num_atoms
            row = list()
            for el in self._pd.elements:
                row.append(comp[el]/sumform)
            m.append(row)
        return np.array(m)

    def _in_facet(self, facet, comp):
        dim = len(self._pd.elements)
        m = np.array([self._pd.qhull_data[i][0:dim-1] for i in facet])
        cm = np.array([m[i] - m[0] for i in xrange(1, len(m))])
        sumform = comp.num_atoms
        row = list()
        for i in xrange(1,len(self._pd.elements)):
            row.append(comp[self._pd.elements[i]]/sumform)
        compm = np.array(row) - m[0]
        
        coeffs = np.linalg.solve(cm.transpose(), compm)
        return (coeffs >= -1e-8).all() and sum(coeffs) <= (1 + 1e-8)

    def _get_facets(self,comp):
        memberfacets = list()
        for facet in self._pd.facets:
            if self._in_facet(facet,comp):
                memberfacets.append(facet)
        return memberfacets

    def get_decomposition(self,comp):
        """
        Provides the decomposition at a particular composition
        Arguments:
            comp - A composition
        Returns:
            Decomposition as a dict of {PDEntry: amount}
        """
        memberfacets = self._get_facets(comp)
        facet = memberfacets[0]
        complist = [self._pd.qhull_entries[i].composition for i in facet]
        m = self._make_comp_matrix(complist)
        compm = self._make_comp_matrix([comp])
        decompamts = np.dot(np.linalg.inv(m.transpose()),compm.transpose())
        decomp = dict()
        #Scrub away zero amounts
        for i in xrange(len(decompamts)):
            if abs(decompamts[i][0]) > 1e-8:
                decomp[self._pd.qhull_entries[facet[i]]] = decompamts[i][0]
        return decomp

    def get_e_above_hull(self,entry):
        """
        Provides the energy above convex hull for an entry
        Arguments:
            entry - A PDEntry like object
        Returns:
            Energy above convex hull of entry.  Stable entries should have energy above hull of 0.
        """
        comp = entry.composition
        eperatom = entry.energy_per_atom
        decomp = self.get_decomposition(comp)
        hullenergy = sum([entry.energy_per_atom*amt for entry, amt in decomp.items()])
        if abs(eperatom) < 1e-8:
            return 0
        return eperatom - hullenergy
    
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
        if ind >= 0:
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
                if abs(c-clean_pots[-1]) > 1e-8:
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
        
        elcomp = Composition.from_formula(element.symbol)
        evolution = []
        for c in chempots:
            gcpd = GrandPotentialPhaseDiagram(stable_entries, {element:c-0.01}, self._pd.elements)
            analyzer = PDAnalyzer(gcpd)
            decomp = [gcentry.original_entry.composition for gcentry, amt in analyzer.get_decomposition(gccomp).items() if amt > 1e-5]
            if elcomp not in decomp:
                decomp.insert(0,elcomp)
            rxn = Reaction([comp], decomp)
            rxn.normalize_to(comp)
            evolution.append({'chempot':c, 'evolution' : - rxn.coeffs[rxn.all_comp.index(elcomp)], 'reaction':rxn})
            
        return evolution
