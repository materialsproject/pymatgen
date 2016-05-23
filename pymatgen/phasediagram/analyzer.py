# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module provides classes for analyzing phase diagrams.
"""

from six.moves import zip

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "May 16, 2012"

import numpy as np
import itertools
import collections

from monty.functools import lru_cache

from pymatgen.core.composition import Composition
from pymatgen.phasediagram.maker import PhaseDiagram, \
    GrandPotentialPhaseDiagram, get_facets
from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.util.coord_utils import Simplex


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
        Initializes analyzer with a PhaseDiagram.

        Args:
            pd: Phase Diagram to analyze.
        """
        self._pd = pd

    def _make_comp_matrix(self, complist):
        """
        Helper function to generates a normalized composition matrix from a
        list of compositions.
        """
        return np.array([[comp.get_atomic_fraction(el)
                          for el in self._pd.elements] for comp in complist])

    @lru_cache(1)
    def _get_facet(self, comp):
        """
        Get any facet that a composition falls into. Cached so successive
        calls at same composition are fast.
        """
        if set(comp.elements).difference(self._pd.elements):
            raise ValueError('{} has elements not in the phase diagram {}'
                             ''.format(comp, self._pd.elements))
        c = [comp.get_atomic_fraction(e) for e in self._pd.elements[1:]]
        for f, s in zip(self._pd.facets, self._pd.simplices):
            if Simplex(s).in_simplex(c, PDAnalyzer.numerical_tol / 10):
                return f
        raise RuntimeError("No facet found for comp = {}".format(comp))

    def get_decomposition(self, comp):
        """
        Provides the decomposition at a particular composition.

        Args:
            comp: A composition

        Returns:
            Decomposition as a dict of {Entry: amount}
        """
        facet = self._get_facet(comp)
        comp_list = [self._pd.qhull_entries[i].composition for i in facet]
        m = self._make_comp_matrix(comp_list)
        compm = self._make_comp_matrix([comp])
        decomp_amts = np.linalg.solve(m.T, compm.T)
        return {self._pd.qhull_entries[f]: amt[0]
                for f, amt in zip(facet, decomp_amts)
                if abs(amt[0]) > PDAnalyzer.numerical_tol}

    def get_hull_energy(self, comp):
        """
        Args:
            comp (Composition): Input composition

        Returns:
            Energy of lowest energy equilibrium at desired composition. Not
            normalized by atoms, i.e. E(Li4O2) = 2 * E(Li2O)
        """
        e = 0
        for k, v in self.get_decomposition(comp).items():
            e += k.energy_per_atom * v
        return e * comp.num_atoms

    def get_decomp_and_e_above_hull(self, entry, allow_negative=False):
        """
        Provides the decomposition and energy above convex hull for an entry.
        Due to caching, can be much faster if entries with the same composition
        are processed together.

        Args:
            entry: A PDEntry like object
            allow_negative: Whether to allow negative e_above_hulls. Used to
                calculate equilibrium reaction energies. Defaults to False.

        Returns:
            (decomp, energy above convex hull)  Stable entries should have
            energy above hull of 0. The decomposition is provided as a dict of
            {Entry: amount}.
        """
        if entry in self._pd.stable_entries:
            return {entry: 1}, 0

        facet = self._get_facet(entry.composition)
        comp_list = [self._pd.qhull_entries[i].composition for i in facet]
        m = self._make_comp_matrix(comp_list)
        compm = self._make_comp_matrix([entry.composition])
        decomp_amts = np.linalg.solve(m.T, compm.T)[:, 0]
        decomp = {self._pd.qhull_entries[facet[i]]: decomp_amts[i]
                  for i in range(len(decomp_amts))
                  if abs(decomp_amts[i]) > PDAnalyzer.numerical_tol}
        energies = [self._pd.qhull_entries[i].energy_per_atom for i in facet]
        ehull = entry.energy_per_atom - np.dot(decomp_amts, energies)
        if allow_negative or ehull >= -PDAnalyzer.numerical_tol:
            return decomp, ehull
        raise ValueError("No valid decomp found!")

    def get_e_above_hull(self, entry):
        """
        Provides the energy above convex hull for an entry

        Args:
            entry: A PDEntry like object

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
            entry: A PDEntry like object

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
        return analyzer.get_decomp_and_e_above_hull(entry,
                                                    allow_negative=True)[1]

    def get_facet_chempots(self, facet):
        """
        Calculates the chemical potentials for each element within a facet.

        Args:
            facet: Facet of the phase diagram.

        Returns:
            { element: chempot } for all elements in the phase diagram.
        """
        complist = [self._pd.qhull_entries[i].composition for i in facet]
        energylist = [self._pd.qhull_entries[i].energy_per_atom for i in facet]
        m = self._make_comp_matrix(complist)
        chempots = np.linalg.solve(m, energylist)
        return dict(zip(self._pd.elements, chempots))

    def get_composition_chempots(self, comp):
        facet = self._get_facet(comp)
        return self.get_facet_chempots(facet)

    def get_transition_chempots(self, element):
        """
        Get the critical chemical potentials for an element in the Phase
        Diagram.

        Args:
            element: An element. Has to be in the PD in the first place.

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
            element: An element. Must be in the phase diagram.
            comp: A Composition
            comp_tol: The tolerance to use when calculating decompositions.
                Phases with amounts less than this tolerance are excluded.
                Defaults to 1e-5.

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
        elcomp = Composition(element.symbol)
        prev_decomp = []
        evolution = []

        def are_same_decomp(decomp1, decomp2):
            for comp in decomp2:
                if comp not in decomp1:
                    return False
            return True

        for c in chempots:
            gcpd = GrandPotentialPhaseDiagram(
                stable_entries, {element: c - 1e-5}, self._pd.elements
            )
            analyzer = PDAnalyzer(gcpd)
            gcdecomp = analyzer.get_decomposition(gccomp)
            decomp = [gcentry.original_entry.composition
                      for gcentry, amt in gcdecomp.items()
                      if amt > comp_tol]
            decomp_entries = [gcentry.original_entry
                              for gcentry, amt in gcdecomp.items()
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

    def get_chempot_range_map(self, elements, referenced=True, joggle=True,
                              force_use_pyhull=False):
        """
        Returns a chemical potential range map for each stable entry.

        Args:
            elements: Sequence of elements to be considered as independent
                variables. E.g., if you want to show the stability ranges
                of all Li-Co-O phases wrt to uLi and uO, you will supply
                [Element("Li"), Element("O")]
            referenced: If True, gives the results with a reference being the
                energy of the elemental phase. If False, gives absolute values.
            joggle (boolean): Whether to joggle the input to avoid precision
                errors.
            force_use_pyhull (boolean): Whether the pyhull algorithm is always
                used, even when scipy is present.

        Returns:
            Returns a dict of the form {entry: [simplices]}. The list of
            simplices are the sides of the N-1 dim polytope bounding the
            allowable chemical potential range of each entry.
        """
        all_chempots = []
        pd = self._pd
        facets = pd.facets
        for facet in facets:
            chempots = self.get_facet_chempots(facet)
            all_chempots.append([chempots[el] for el in pd.elements])
        inds = [pd.elements.index(el) for el in elements]
        el_energies = {el: 0.0 for el in elements}
        if referenced:
            el_energies = {el: pd.el_refs[el].energy_per_atom
                           for el in elements}
        chempot_ranges = collections.defaultdict(list)
        vertices = [list(range(len(self._pd.elements)))]
        if len(all_chempots) > len(self._pd.elements):
            vertices = get_facets(all_chempots, joggle=joggle,
                                  force_use_pyhull=force_use_pyhull)
        for ufacet in vertices:
            for combi in itertools.combinations(ufacet, 2):
                data1 = facets[combi[0]]
                data2 = facets[combi[1]]
                common_ent_ind = set(data1).intersection(set(data2))
                if len(common_ent_ind) == len(elements):
                    common_entries = [pd.qhull_entries[i]
                                      for i in common_ent_ind]
                    data = np.array([[all_chempots[i][j]
                                      - el_energies[pd.elements[j]]
                                      for j in inds] for i in combi])
                    sim = Simplex(data)
                    for entry in common_entries:
                        chempot_ranges[entry].append(sim)

        return chempot_ranges

    def getmu_vertices_stability_phase(self, target_comp, dep_elt, tol_en=1e-2):
        """
        returns a set of chemical potentials corresponding to the vertices of the simplex
        in the chemical potential phase diagram.
        The simplex is built using all elements in the target_composition except dep_elt.
        The chemical potential of dep_elt is computed from the target composition energy.
        This method is useful to get the limiting conditions for
        defects computations for instance.

        Args:
            target_comp: A Composition object
            dep_elt: the element for which the chemical potential is computed from the energy of
            the stable phase at the target composition
            tol_en: a tolerance on the energy to set

        Returns:
             [{Element:mu}]: An array of conditions on simplex vertices for
             which each element has a chemical potential set to a given
             value. "absolute" values (i.e., not referenced to element energies)
        """
        muref = np.array([self._pd.el_refs[e].energy_per_atom
                          for e in self._pd.elements if e != dep_elt])
        chempot_ranges = self.get_chempot_range_map(
            [e for e in self._pd.elements if e != dep_elt])

        for e in self._pd.elements:
            if not e in target_comp.elements:
                target_comp = target_comp + Composition({e: 0.0})
        coeff = [-target_comp[e] for e in self._pd.elements if e != dep_elt]
        for e in chempot_ranges.keys():
            if e.composition.reduced_composition == \
                    target_comp.reduced_composition:
                multiplicator = e.composition[dep_elt] / target_comp[dep_elt]
                ef = e.energy / multiplicator
                all_coords = []
                for s in chempot_ranges[e]:
                    for v in s._coords:
                        elts = [e for e in self._pd.elements if e != dep_elt]
                        res = {}
                        for i in range(len(elts)):
                            res[elts[i]] = v[i] + muref[i]
                        res[dep_elt]=(np.dot(v+muref, coeff)+ef)/target_comp[dep_elt]
                        already_in = False
                        for di in all_coords:
                            dict_equals = True
                            for k in di:
                                if abs(di[k]-res[k]) > tol_en:
                                    dict_equals = False
                                    break
                            if dict_equals:
                                already_in = True
                                break
                        if not already_in:
                            all_coords.append(res)
        return all_coords

    def get_chempot_range_stability_phase(self, target_comp, open_elt):
        """
        returns a set of chemical potentials correspoding to the max and min
        chemical potential of the open element for a given composition. It is
        quite common to have for instance a ternary oxide (e.g., ABO3) for
        which you want to know what are the A and B chemical potential leading
        to the highest and lowest oxygen chemical potential (reducing and
        oxidizing conditions). This is useful for defect computations.

        Args:
            target_comp: A Composition object
            open_elt: Element that you want to constrain to be max or min

        Returns:
             {Element:(mu_min,mu_max)}: Chemical potentials are given in
             "absolute" values (i.e., not referenced to 0)
        """
        muref = np.array([self._pd.el_refs[e].energy_per_atom
                          for e in self._pd.elements if e != open_elt])
        chempot_ranges = self.get_chempot_range_map(
            [e for e in self._pd.elements if e != open_elt])
        for e in self._pd.elements:
            if not e in target_comp.elements:
                target_comp = target_comp + Composition({e: 0.0})
        coeff = [-target_comp[e] for e in self._pd.elements if e != open_elt]
        max_open = -float('inf')
        min_open = float('inf')
        max_mus = None
        min_mus = None
        for e in chempot_ranges.keys():
            if e.composition.reduced_composition == \
                    target_comp.reduced_composition:
                multiplicator = e.composition[open_elt] / target_comp[open_elt]
                ef = e.energy / multiplicator
                all_coords = []
                for s in chempot_ranges[e]:
                    for v in s._coords:
                        all_coords.append(v)
                        if (np.dot(v + muref, coeff) + ef) / target_comp[
                            open_elt] > max_open:
                            max_open = (np.dot(v + muref, coeff) + ef) / \
                                       target_comp[open_elt]
                            max_mus = v
                        if (np.dot(v + muref, coeff) + ef) / target_comp[
                            open_elt] < min_open:
                            min_open = (np.dot(v + muref, coeff) + ef) / \
                                       target_comp[open_elt]
                            min_mus = v
        elts = [e for e in self._pd.elements if e != open_elt]
        res = {}
        for i in range(len(elts)):
            res[elts[i]] = (min_mus[i] + muref[i], max_mus[i] + muref[i])
        res[open_elt] = (min_open, max_open)
        return res
