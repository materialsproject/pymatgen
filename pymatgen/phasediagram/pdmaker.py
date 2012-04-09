#!/usr/bin/env python

"""
This module provides classes to create phase diagrams.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import numpy as np
import logging
import itertools

from scipy.spatial import Delaunay

from pymatgen.core.structure import Composition
from pymatgen.command_line.qhull_caller import qconvex
from pymatgen.phasediagram.entries import GrandPotPDEntry, TransformedPDEntry
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.reaction_calculator import Reaction

logger = logging.getLogger(__name__)

class PhaseDiagram (object):
    '''
    Simple phase diagram class taking in elements and entries as inputs.
    '''
    FORMATION_ENERGY_TOLERANCE = 1e-11

    def __init__(self, entries, elements=None, use_external_qhull=False):
        """
        Standard constructor for phase diagram.
        
        Args:
            entries:
                A list of PDEntry-like objects having an energy, energy_per_atom and composition.
            elements:
                Optional list of elements in the phase diagram. If set to None, the elements are determined from the
                the entries themselves.
        """
        if elements == None:
            elements = set()
            map(elements.update, [entry.composition.elements for entry in entries])
        self._all_entries = entries
        self._elements = tuple(elements)
        self._qhull_data = None
        self._facets = None
        self._qhull_entries = None
        self._stable_entries = None
        self._all_entries_hulldata = None
        self._use_external_qhull = use_external_qhull
        self._make_phasediagram()

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

    def get_form_energy(self, entry):
        '''
        Returns the formation energy for an entry (NOT normalized) from the
        elemental references.
        
        Args:
            entry:
                A PDEntry
        
        Returns:
            Formation energy from the elementals references.
        '''
        comp = entry.composition
        energy = entry.energy - sum([comp[el] * self._el_refs[el].energy_per_atom for el in comp.elements])
        return energy

    def get_form_energy_per_atom(self, entry):
        '''
        Returns the formation energy per atom for an entry from the
        elemental references.
        '''
        comp = entry.composition
        return self.get_form_energy(entry) / comp.num_atoms

    def _process_entries_qhulldata(self, entries_to_process):
        data = list()
        for entry in entries_to_process:
            comp = entry.composition
            energy_per_atom = entry.energy_per_atom
            row = list()
            for i in xrange(1, len(self._elements)):
                row.append(comp.get_atomic_fraction(self._elements[i]))
            row.append(energy_per_atom)
            data.append(row)
        return data

    def _create_convhull_data(self):
        '''
        Make data suitable for convex hull procedure from the list of entries.
        '''
        logger.debug("Creating convex hull data...")
        #Determine the elemental references based on lowest energy for each.
        self._el_refs = dict()
        for entry in self._all_entries:
            if entry.composition.is_element:
                for el in entry.composition.elements:
                    if entry.composition[el] > Composition.amount_tolerance:
                        break
                e_per_atom = entry.energy_per_atom
                if el not in self._el_refs:
                    self._el_refs[el] = entry
                elif self._el_refs[el].energy_per_atom > e_per_atom:
                    self._el_refs[el] = entry
        # Remove positive formation energy entries
        entries_to_process = list()
        for entry in self._all_entries:
            if self.get_form_energy(entry) <= self.FORMATION_ENERGY_TOLERANCE:
                entries_to_process.append(entry)
            else:
                logger.debug("Removing positive formation energy entry {}".format(entry))

        self._qhull_entries = entries_to_process
        return self._process_entries_qhulldata(entries_to_process)

    def _make_phasediagram(self):
        """
        Make the phase diagram.
        """
        stable_entries = set()
        dim = len(self._elements)
        self._qhull_data = self._create_convhull_data()
        if len(self._qhull_data) == dim:
            self._facets = [range(len(self._elements))]
        else:
            if self._use_external_qhull:
                logger.debug("> 4D hull encountered. Computing hull using external qconvex call.")
                self._facets = qconvex(self._qhull_data)
            else:
                logger.debug("Computing hull using scipy.spatial.delaunay")
                delau = Delaunay(self._qhull_data)
                self._facets = delau.convex_hull
            logger.debug("Final facets are\n{}".format(self._facets))

            logger.debug("Removing vertical facets...")
            finalfacets = list()
            for facet in self._facets:
                facetmatrix = np.zeros((len(facet), len(facet)))
                count = 0
                is_element_facet = True
                for vertex in facet:
                    facetmatrix[count] = np.array(self._qhull_data[vertex])
                    facetmatrix[count, dim - 1] = 1
                    count += 1
                    if len(self._qhull_entries[vertex].composition) > 1:
                        is_element_facet = False
                if abs(np.linalg.det(facetmatrix)) > 1e-8 and (not is_element_facet):
                    finalfacets.append(facet)
                else:
                    logger.debug("Removing vertical facet : {}".format(facet))
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

    def __init__(self, entries, chempots, elements, use_external_qhull=False):
        """
        Standard constructor for phase diagram.
        
        Args:
            entries:
                A list of PDEntry-like objects having an energy, energy_per_atom and composition.
            chempots:
                A dict of {element: float} to specify the chemical potentials of the open elements.
            elements:
                Optional list of elements in the phase diagram. If set to None, the elements are determined from the
                the entries themselves.
        """
        allentries = list()
        for entry in entries:
            if not (entry.is_element and (entry.composition.elements[0] in chempots)):
                allentries.append(GrandPotPDEntry(entry, chempots))
        self.chempots = chempots
        filteredels = list()
        for el in elements:
            if el not in chempots:
                filteredels.append(el)
        elements = sorted(filteredels)
        super(GrandPotentialPhaseDiagram, self).__init__(allentries, elements, use_external_qhull)


class CompoundPhaseDiagram(PhaseDiagram):
    """
    Experimental feature. Generates phase diagrams from compounds as termninations
    instead of elements.
    """

    def __init__(self, entries, terminal_compositions, use_external_qhull=False):
        entries = get_entries_within_compositional_space(entries, terminal_compositions)
        elset = get_non_coplanar_element_set(entries)
        els = list(elset)
        pentries = get_transformed_entries(entries, els)
        super(CompoundPhaseDiagram, self).__init__(pentries, use_external_qhull=use_external_qhull)


def get_comp_matrix_from_comp(compositions, elements, normalize_row=True):
    """
    Helper function to generates a normalized composition matrix from a list of 
    composition.
    """
    comp_matrix = np.array([[comp.get_atomic_fraction(el) for el in elements] for comp in compositions])
    if not normalize_row:
        return comp_matrix
    factor = np.tile(np.sum(comp_matrix, 1), (len(elements), 1)).transpose()
    return comp_matrix / factor

def get_comp_matrix(entries, elements, normalize_row=True):
    """
    Helper function to generates a normalized composition matrix from a list of 
    composition.
    """
    return get_comp_matrix_from_comp([entry.composition for entry in entries], elements, normalize_row)

def is_coplanar(entries, elements):
    comp_matrix = get_comp_matrix(entries, elements)
    for submatrix in itertools.combinations(comp_matrix, min(len(elements), len(entries))):
        if abs(np.linalg.det(submatrix)) > 1e-5:
            return False
    return True

def get_non_coplanar_element_set(entries):
    elements = set()
    map(elements.update, [entry.composition.elements for entry in entries])
    for i in xrange(len(elements), 1, -1):
        for elset in itertools.combinations(elements, i):
            if not is_coplanar(entries, elset):
                return elset
    return None

def get_transformed_entries(entries, elements):
    comp_matrix = get_comp_matrix(entries, elements)
    newmat = []
    energies = []
    for i in xrange(len(elements)):
        col = comp_matrix[:, i]
        maxval = max(col)
        maxind = list(col).index(maxval)
        newmat.append(comp_matrix[maxind])
        energies.append(entries[i].energy_per_atom)
    invm = np.linalg.inv(np.array(newmat).transpose())
    newentries = []
    for i in xrange(len(entries)):
        entry = entries[i]
        lincomp = np.dot(invm, comp_matrix[i])
        lincomp = np.around(lincomp, 5)
        comp = Composition({Element.from_Z(j + 1):lincomp[j] for j in xrange(len(elements))})
        scaled_energy = entry.energy_per_atom - sum(lincomp * energies)
        newentries.append(TransformedPDEntry(comp, scaled_energy, entry))
    return newentries

def get_entries_within_compositional_space(entries, terminal_compositions):
    newentries = []
    for entry in entries:
        try:
            rxn = Reaction(terminal_compositions, [entry.composition])
            newentries.append(entry)
        except:
            pass
    return newentries


