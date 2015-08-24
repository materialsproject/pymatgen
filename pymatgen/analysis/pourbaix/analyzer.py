# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

"""
Class for analyzing Pourbaix Diagrams. Similar to PDAnalyzer
"""

__author__ = "Sai Jayaraman"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.0"
__maintainer__ = "Sai Jayaraman"
__email__ = "sjayaram@mit.edu"
__status__ = "Development"
__date__ = "Nov 7, 2012"

import numpy as np

from pyhull.simplex import Simplex
from functools import cmp_to_key
from pyhull.halfspace import Halfspace, HalfspaceIntersection
from pyhull.convex_hull import ConvexHull
from six.moves import zip


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

        Args:
            elements: Sequence of elements to be considered as independent
                variables. E.g., if you want to show the stability ranges of
                all Li-Co-O phases wrt to uLi and uO, you will supply
                [Element("Li"), Element("O")]

        Returns:
            Returns a dict of the form {entry: [simplices]}. The list of
            simplices are the sides of the N-1 dim polytope bounding the
            allowable chemical potential range of each entry.
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

        basis_vecs = []
        on_plane_points = []
        # Create basis vectors
        for entry in self._pd.stable_entries:
            ie = self._pd.qhull_entries.index(entry)
            row = self._pd._qhull_data[ie]
            on_plane_points.append([0, 0, row[2]])
            this_basis_vecs = []
            norm_vec = [-0.0591 * row[0], -1 * row[1], 1]
            if abs(norm_vec[0]) > tol:
                this_basis_vecs.append([-norm_vec[2]/norm_vec[0], 0, 1])
            if abs(norm_vec[1]) > tol:
                this_basis_vecs.append([0, -norm_vec[2]/norm_vec[1], 1])
            if len(this_basis_vecs) == 0:
                basis_vecs.append([[1, 0, 0], [0, 1, 0]])
            elif len(this_basis_vecs) == 1:
                if abs(this_basis_vecs[0][0]) < tol:
                    this_basis_vecs.append([1, 0, 0])
                else:
                    this_basis_vecs.append([0, 1, 0])
                basis_vecs.append(this_basis_vecs)
            else:
                basis_vecs.append(this_basis_vecs)

        # Find point in half-space in which optimization is desired
        ph_max_contrib = -1 * max([abs(0.0591 * row[0])
                                    for row in self._pd._qhull_data]) * limits[0][1]
        V_max_contrib = -1 * max([abs(row[1]) for row in self._pd._qhull_data]) * limits[1][1]
        g_max = (-1 * max([abs(pt[2]) for pt in on_plane_points])
                  + ph_max_contrib + V_max_contrib) - 10
        point_in_region = [7, 0, g_max]

        # Append border hyperplanes along limits
        for i in range(len(limits)):
            for j in range(len(limits[i])):
                basis_vec_1 = [0.0] * 3
                basis_vec_2 = [0.0] * 3
                point = [0.0] * 3
                basis_vec_1[2] = 1.0
                basis_vec_2[2] = 0.0
                for axis in range(len(limits)):
                    if axis is not i:
                        basis_vec_1[axis] = 0.0
                        basis_vec_2[axis] = 1.0
                basis_vecs.append([basis_vec_1, basis_vec_2])
                point[i] = limits[i][j]
                on_plane_points.append(point)

        # Hyperplane enclosing the very bottom
        basis_vecs.append([[1, 0, 0], [0, 1, 0]])
        on_plane_points.append([0, 0, 2 * g_max])
        hyperplane_list = [Halfspace.from_hyperplane(basis_vecs[i], on_plane_points[i], point_in_region)
                            for i in range(len(basis_vecs))]
        hs_int = HalfspaceIntersection(hyperplane_list, point_in_region)
        int_points = hs_int.vertices
        pourbaix_domains = {}
        self.pourbaix_domain_vertices = {}

        for i in range(len(self._pd.stable_entries)):
            vertices = [[int_points[vert][0], int_points[vert][1]] for vert in
                         hs_int.facets_by_halfspace[i]]
            if len(vertices) < 1:
                continue
            pourbaix_domains[self._pd.stable_entries[i]] = ConvexHull(vertices).simplices

            # Need to order vertices for highcharts area plot
            cx = sum([vert[0] for vert in vertices]) / len(vertices)
            cy = sum([vert[1] for vert in vertices]) / len(vertices)
            point_comp = lambda x, y: x[0]*y[1] - x[1]*y[0]
            vert_center = [[v[0] - cx, v[1] - cy] for v in vertices]
            vert_center.sort(key=cmp_to_key(point_comp))
            self.pourbaix_domain_vertices[self._pd.stable_entries[i]] =\
             [[v[0] + cx, v[1] + cy] for v in vert_center]

        self.pourbaix_domains = pourbaix_domains
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
        Get any facet that a composition falls into.
        """
        for facet in self._pd.facets:
            if self._in_facet(facet, entry):
                return facet
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
        decompamts = np.dot(np.linalg.inv(m.transpose()), compm.transpose())
        decomp = dict()
        #Scrub away zero amounts
        for i in range(len(decompamts)):
            if abs(decompamts[i][0]) > PourbaixAnalyzer.numerical_tol:
                decomp[self._pd.qhull_entries[facet[i]]] = decompamts[i][0]
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
        hullenergy = sum([entry.g0 * amt
                          for entry, amt in decomp.items()])
        return decomp, g0 - hullenergy

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
