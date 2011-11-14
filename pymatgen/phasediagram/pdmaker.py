#!/usr/bin/env python

"""
This module provides classes to create phase diagrams.
"""

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import numpy as np
import warnings

from entries import GrandPotPDEntry
from scipy.spatial import Delaunay

class PhaseDiagram (object):
    '''
    Simple phase diagram class taking in elements and entries as inputs.
    '''
    FORMATION_ENERGY_TOLERANCE = 1e-11
    
    def __init__(self, entries,elements = None):
        """
        Standard constructor for phase diagram.
        Arguments:
            entries - a list of PDEntry-like objects having an energy, energy_per_atom and composition.
            elements - Optional list of elements in the phase diagram. If set to None, the elements are determined from the
                       the entries themselves.
        """
        if elements == None:
            elements = set()
            map(elements.update, [entry.composition.elements for entry in entries])
        self._all_entries = entries
        self._elements = sorted(tuple(elements))
        self._qhull_data = None
        self._facets = None
        self._qhull_entries = None
        self._stable_entries = None
        self._all_entries_hulldata = None
        self._el_refs = None
        self.make_phasediagram()

    @property
    def all_entries(self):
        """
        All entries submitted to construct the PhaseDiagram.
        Note that this does not mean that all these entries are actually used in the phase diagram.
        """
        return self._all_entries

    @property
    def elements(self):
        """
        Elements in the phase diagram.
        """
        return self._elements

    @property
    def facets(self):
        """
        Facets of the phase diagram in the form of  [[1,2,3],[4,5,6]...]
        """
        return self._facets

    @property
    def qhull_data(self):
        """
        Data used in the convex hull operation. This is essentially a matrix of
        composition data and energy per atom values created from the qhull_entries.
        """
        return self._qhull_data

    @property
    def qhull_entries(self):
        """
        Actual entries used in convex hull. Excludes all positive formation energy entries.
        """
        return self._qhull_entries
    
    @property
    def unstable_entries(self):
        """
        Entries that are unstable in the phase diagram. Includes positive formation energy entries.
        """
        return [e for e in self.all_entries if e not in self.stable_entries]
    
    @property
    def stable_entries(self):
        '''
        Returns the stable entries.
        '''
        return self._stable_entries

    @property
    def all_entries_hulldata(self):
        """
        Same as qhull_data, but for all entries rather than just positive formation energy ones.
        """
        return self._all_entries_hulldata

    @property
    def el_refs(self):
        """
        List of elemental references for the phase diagrams. 
        These are typically entries corresponding to the lowest energy element entries.
        """
        return self._el_refs


    def get_form_energy(self,entry):
        '''
        Returns the formation energy for an entry (NOT normalized) from the
        elemental references.
        Arguments:
            entry - A PDEntry
        Returns:
            formation energy from the elementals references.
        '''
        comp = entry.composition
        energy = entry.energy - sum([comp[el]*self.el_ref[el].energy_per_atom for el in comp.elements])
        return energy
        
    def get_form_energy_per_atom(self,entry):
        '''
        Returns the formation energy per atom for an entry from the
        elemental references.
        '''
        comp = entry.composition
        return self.get_form_energy(entry)/comp.num_atoms

    def _process_entries_qhulldata(self,entries_to_process):
        data = list()
        for entry in entries_to_process:
            comp = entry.composition

            energy_per_atom = entry.energy_per_atom
            row = list()
            for i in xrange(1,len(self._elements)):
                row.append(comp.get_atomic_fraction(self._elements[i]))
            row.append(energy_per_atom)
            data.append(row)
        return data

    def _create_convhull_data(self):
        '''
        Make data suitable for convex hull procedure from the list of entries.
        '''
        #Determine the elemental references based on lowest energy for each.
        self.el_ref = dict()
        for entry in self._all_entries:
            if entry.composition.is_element:
                el = entry.composition.elements[0]
                e_per_atom = entry.energy_per_atom
                if el not in self.el_ref:
                    self.el_ref[el] = entry
                elif self.el_ref[el].energy_per_atom > e_per_atom:
                    self.el_ref[el] = entry

        # Remove positive formation energy entries and duplicate entries
        def in_list(entries_list, entry):
            for test_entry in entries_list:
                if abs(entry.energy - test_entry.energy) <= self.FORMATION_ENERGY_TOLERANCE and entry.composition == test_entry.composition:
                    warnings.warn("Duplicate entry found!  Discarding...")
                    return True
            return False
        
        entries_to_process = list()
        for entry in self._all_entries:
            if self.get_form_energy(entry) <= self.FORMATION_ENERGY_TOLERANCE:# and (not in_list(entries_to_process, entry)):
                entries_to_process.append(entry)
                
        self._qhull_entries = entries_to_process
        return self._process_entries_qhulldata(entries_to_process)

    def make_phasediagram(self):
        stable_entries = set()
        dim = len(self._elements)
        self._qhull_data = self._create_convhull_data()
        if len(self._qhull_data) == dim:
            self._facets = [range(len(self._elements))]
        else:
            self._facets = Delaunay(self._qhull_data).convex_hull
            finalfacets = list()

            # Remove vertical facets and elemental facets
            for facet in self._facets:
                facetmatrix = np.zeros((len(facet),len(facet)))
                count = 0
                is_element_facet = True

                for vertex in facet:
                    facetmatrix[count] = np.array(self._qhull_data[vertex])
                    facetmatrix[count, dim-1] = 1
                    count += 1
                    if len(self._qhull_entries[vertex].composition) > 1:
                        is_element_facet = False
                if abs(np.linalg.det(facetmatrix)) > 1e-8 and (not is_element_facet):
                    finalfacets.append(facet)
            self._facets = finalfacets

        for facet in self._facets:
            for vertex in facet:
                stable_entries.add(self._qhull_entries[vertex])
        self._stable_entries = stable_entries
        self._all_entries_hulldata = self._process_entries_qhulldata(self._all_entries)


class GrandPotentialPhaseDiagram (PhaseDiagram):
    '''
    Grand potential phase diagram class taking in elements and entries as inputs.
    '''

    def __init__(self, entries, chempots, elements):
        """
        Standard constructor for phase diagram.
        Arguments:
            entries - a list of PDEntry-like objects having an energy, energy_per_atom and composition.
            chempots - a dict of {element: float} to specify the chemical potentials of the open elements.
            elements - Optional list of elements in the phase diagram. If set to None, the elements are determined from the
                       the entries themselves.
        """
        allentries = list()
        for entry in entries:
            if not (entry.is_element and (entry.composition.elements[0] in chempots)):
                allentries.append(GrandPotPDEntry(entry,chempots))
        self.chempots = chempots
        filteredels = list()
        for el in elements:
            if el not in chempots:
                filteredels.append(el)
        elements = sorted(filteredels)
        super(GrandPotentialPhaseDiagram,self).__init__(allentries, elements)


