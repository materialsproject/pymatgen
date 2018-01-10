# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import numpy as np

from pymatgen.util.coord import Simplex
from functools import cmp_to_key
from scipy.spatial import HalfspaceIntersection, ConvexHull
from pymatgen.analysis.pourbaix.entry import MultiEntry
from six.moves import zip
import warnings

"""
Class for analyzing Pourbaix Diagrams. Similar to PDAnalyzer
"""

__author__ = "Sai Jayaraman"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.0"
__maintainer__ = "Sai Jayaraman"
__credits__ = "Arunima Singh, Joseph Montoya"
__email__ = "sjayaram@mit.edu"
__status__ = "Development"
__date__ = "Nov 7, 2012"


class PourbaixAnalyzer(object):
    """
    Class for performing analysis on Pourbaix Diagrams

    Args:
        pd: Pourbaix Diagram to analyze.
    """
    numerical_tol = 1e-8

    def __init__(self, pd):
        self._pd = pd
        self._keys = ['H+', 'V', '1']
        self.chempot_limits = None

    def get_facet_chempots(self, facet):
        """
        Calculates the chemical potentials for each element within a facet.

        Args:
            facet: Facet of the phase diagram.

        Returns:
            { element: chempot } for all elements in the phase diagram.
        """
        entrylist = [self._pd.qhull_entries[i] for i in facet]
        energylist = [self._pd.qhull_entries[i].g0 for i in facet]
        m = self._make_comp_matrix(entrylist)
        chempots = np.dot(np.linalg.inv(m), energylist)

        return dict(zip(self._keys, chempots))

    def _make_comp_matrix(self, entrylist):
        """
        Helper function to generates a normalized composition matrix from a
        list of Pourbaix Entries
        """
        return np.array([[entry.npH, entry.nPhi, 1] for entry in entrylist])

    def get_chempot_range_map(self, limits=[[-2,16], [-4,4]]):
        """
        Returns a chemical potential range map for each stable entry.

        This function works by using scipy's HalfspaceIntersection
        function to construct all of the 2-D polygons that form the
        boundaries of the planes corresponding to individual entry
        gibbs free energies as a function of pH and V. Hyperplanes
        of the form a*pH + b*V + 1 - g(0, 0) are constructed and
        supplied to HalfspaceIntersection, which then finds the
        boundaries of each pourbaix region using the intersection
        points.

        Args:
            limits ([[float]]): limits in which to do the pourbaix
                analysis

        Returns:
            Returns a dict of the form {entry: [boundary_points]}. 
            The list of boundary points are the sides of the N-1 
            dim polytope bounding the allowable ph-V range of each entry.
        """
        tol = PourbaixAnalyzer.numerical_tol
        all_chempots = []
        facets = self._pd.facets
        for facet in facets:
            chempots = self.get_facet_chempots(facet)
            chempots["H+"] /= -0.0591
            chempots["V"] = -chempots["V"]
            chempots["1"] = chempots["1"]
            all_chempots.append([chempots[el] for el in self._keys])

        # Get hyperplanes corresponding to G as function of pH and V
        halfspaces = []
        qhull_data = np.array(self._pd._qhull_data)
        stable_entries = self._pd.stable_entries
        stable_indices = [self._pd.qhull_entries.index(e)
                          for e in stable_entries]
        qhull_data = np.array(self._pd._qhull_data)
        hyperplanes = np.vstack([-0.0591 * qhull_data[:, 0], -qhull_data[:, 1],
                                 np.ones(len(qhull_data)), -qhull_data[:, 2]])
        hyperplanes = np.transpose(hyperplanes)
        max_contribs = np.max(np.abs(hyperplanes), axis=0)
        g_max = np.dot(-max_contribs, [limits[0][1], limits[1][1], 0, 1])

        # Add border hyperplanes and generate HalfspaceIntersection
        border_hyperplanes = [[-1, 0, 0, limits[0][0]],
                              [1, 0, 0, -limits[0][1]],
                              [0, -1, 0, limits[1][0]],
                              [0, 1, 0, -limits[1][1]],
                              [0, 0, -1, 2 * g_max]]
        hs_hyperplanes = np.vstack([hyperplanes[stable_indices],
                                    border_hyperplanes])
        interior_point = np.average(limits, axis=1).tolist() + [g_max]
        hs_int = HalfspaceIntersection(hs_hyperplanes, np.array(interior_point))

        # organize the boundary points by entry
        pourbaix_domains = {entry: [] for entry in stable_entries}
        for intersection, facet in zip(hs_int.intersections, 
                                       hs_int.dual_facets):
            for v in facet:
                if v < len(stable_entries):
                    pourbaix_domains[stable_entries[v]].append(intersection)

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

        self.pourbaix_domains = pourbaix_domains
        self.pourbaix_domain_vertices = pourbaix_domain_vertices
        return pourbaix_domains

    def _in_facet(self, facet, entry):
        """
        Checks if a Pourbaix Entry is in a facet.

        Args:
            facet: facet to test.
            entry: Pourbaix Entry to test.
        """
        dim = len(self._keys)
        if dim > 1:
            coords = [np.array(self._pd.qhull_data[facet[i]][0:dim - 1])
                      for i in range(len(facet))]
            simplex = Simplex(coords)
            comp_point = [entry.npH, entry.nPhi]
            return simplex.in_simplex(comp_point,
                                      PourbaixAnalyzer.numerical_tol)
        else:
            return True

    def _get_facets(self, entry):
        """
        Get the facets that an entry falls into.
        """
        memberfacets = list()
        for facet in self._pd.facets:
            if self._in_facet(facet, entry):
                memberfacets.append(facet)
        return memberfacets

    def _get_facet(self, entry):
        """
        Get any facet that a Pourbaix Entry falls into.
        """
        for facet in self._pd.facets:
            if self._in_facet(facet, entry):
                return facet
        raise RuntimeError("No facet found for comp = {}".format(entry.name))

    def _get_all_facets(self, entry):
        """
        Get all the facets that a Pourbaix Entry falls into
        """
        all_facets = []
        for facet in self._pd.facets:
            if self._in_facet(facet,entry):
               all_facets.append(facet)
        return all_facets
        raise RuntimeError("No facet found for comp = {}".format(entry.name))


    def _get_facet_entries(self, facet):
        """
        Get the entries corresponding to a facet
        """
        entries = []
        for vertex in facet:
            entries.append(self._pd.qhull_entries[vertex])
        return entries

    def g(self, entry, pH, V):
        """
        Get free energy for a given pH, and V.
        """
        g0 = entry.g0
        npH = -entry.npH * 0.0591
        nPhi = -entry.nPhi
        return g0 - npH * pH - nPhi * V

    def get_all_decomp_and_e_above_hull(self, single_entry):
        """
        Computes the decomposition entries, species and hull energies 
        for all the multi-entries which have the "material" as the only solid.  
        
        Args:
            single_entry: single entry for which to find all of the
                decompositions

        Returns:
            (decomp_entries, hull_energies, decomp_species, entries)
            for all multi_entries containing the single_entry as the
            only solid
        """
        decomp_entries, hull_energies, decomp_species, entries = [], [], [], []

        # for all entries where the material is the only solid
        if not self._pd._multielement:
           possible_entries = [e for e in self._pd.all_entries
                               if single_entry == e]
        else:
           possible_entries = [e for e in self._pd.all_entries
                               if e.phases.count("Solid") == 1
                               and single_entry in e.entrylist]
        
        for possible_entry in possible_entries:
            # Find the decomposition details if the material
            # is in the Pourbaix Multi Entry or Pourbaix Entry
            facets = self._get_all_facets(possible_entry)
            for facet in facets:
                entrylist = [self._pd.qhull_entries[i] for i in facet]
                m = self._make_comp_matrix(entrylist)
                compm = self._make_comp_matrix([possible_entry])
                decomp_amts = np.dot(np.linalg.inv(m.transpose()), compm.transpose())
                decomp, decomp_names = {}, {}
                for i, decomp_amt in enumerate(decomp_amts):
                    if abs(decomp_amt[0]) > PourbaixAnalyzer.numerical_tol:
                        decomp[self._pd.qhull_entries[facet[i]]] = decomp_amt[0]
                decomp_entries.append(decomp)
                hull_energy = sum([entry.g0 * amt for entry, amt in decomp.items()])
                hull_energies.append(possible_entry.g0 - hull_energy)
                entries.append(possible_entry)
        
        return decomp_entries, hull_energies, entries

    def get_decomposition(self, entry):
        """
        Provides the decomposition at a particular composition

        Args:
            comp: A composition

        Returns:
            Decomposition as a dict of {PourbaixEntry: amount}
        """
        facet = self._get_facet(entry)
        entrylist = [self._pd.qhull_entries[i] for i in facet]
        m = self._make_comp_matrix(entrylist)
        compm = self._make_comp_matrix([entry])
        decomp_amts = np.dot(np.linalg.inv(m.transpose()), compm.transpose())
        decomp = dict()
        self.decomp_names = dict()
        #Scrub away zero amounts
        for i in range(len(decomp_amts)):
            if abs(decomp_amts[i][0]) > PourbaixAnalyzer.numerical_tol:
                decomp[self._pd.qhull_entries[facet[i]]] = decomp_amts[i][0]
                self.decomp_names[self._pd.qhull_entries[facet[i]].name] = decomp_amts[i][0]
        return decomp

    def get_decomp_and_e_above_hull(self, entry):
        """
        Provides the decomposition and energy above convex hull for an entry

        Args:
            entry: A PourbaixEntry

        Returns:
            (decomp, energy above convex hull)  Stable entries should have
            energy above hull of 0.
        """
        g0 = entry.g0
        decomp = self.get_decomposition(entry)
        hull_energy = sum([entry.g0 * amt
                          for entry, amt in decomp.items()])
        return decomp, g0 - hull_energy, self.decomp_names

    def get_e_above_hull(self, entry):
        """
        Provides the energy above convex hull for an entry

        Args:
            entry: A PourbaixEntry object

        Returns:
            Energy above convex hull of entry. Stable entries should have
            energy above hull of 0.
        """
        return self.get_decomp_and_e_above_hull(entry)[1]

    # TODO: we might want to rename this, still a bit ambiguous
    def get_gibbs_free_energy(self, pH, V):
        """
        Provides the gibbs free energy of the Pourbaix stable entry
        at a given pH and V

        Args:
            pH: pH
             V: potential vs SHE

        Returns:
             gibbs free energy (eV/atom) 
        """
        data = {}

        for entry in self._pd.stable_entries:
            data.update({entry.name: self.g(entry, pH, V)})
        gibbs_energy = min(data.values())
        stable_entry = [k for k, v in data.items() if v == gibbs_energy]
        return (gibbs_energy, stable_entry)

    def _min_multientry_from_single_entry(self, single_entry):
        """
        Gives lowest energy multi-entry from single entry

        Args:
            single_entry (PourbaixEntry): pourbaix entry to find valid
                multientries from
        """
        de, ehulls, entries = self.get_all_decomp_and_e_above_hull(single_entry)
        if not ehulls:
            raise ValueError("No entries where {} is the only solid".format(
                             single_entry.name))
        return entries[np.argmin(ehulls)]

    def get_entry_stability(self, entry, pH, V):
        """
        Get the energy difference between an entry and the
        most stable decomposition product (i.e. the pourbaix-stable
        entry) at a given pH and voltage.

        Args:
            entry (PourbaixEntry): Pourbaix entry or MultiEntry
                corresponding to the stability to be calculated
            pH (float): pH at which to calculate stability of entry
            V (float): voltage at which to calculate stability of entry
        """
        if self._pd._multielement and not isinstance(entry, MultiEntry):
            _, _, entries = self.get_all_decomp_and_e_above_hull(entry)
            warnings.warn("{} is not a multi-entry, calculating stability of "
                          "representative {} multientry")
            gs = [self.g(e, pH, V) for e in entries]
            return min(gs) - self.get_gibbs_free_energy(pH, V)[0]
        return self.g(entry, pH, V) - self.get_gibbs_free_energy(pH, V)[0]
