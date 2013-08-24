#!/usr/bin/env python

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
import itertools
import collections

from pyhull.simplex import Simplex


class PourbaixAnalyzer(object):
    """
    Class for performing analysis on Pourbaix Diagrams
    """
    numerical_tol = 1e-8

    def __init__(self, pd):
        """
        Args:
            pd:
                Pourbaix Diagram to analyze.
        """
        self._pd = pd
        self._keys = ['H+', 'V', '1']
        self.chempot_limits = None

    def get_facet_chempots(self, facet):
        """
        Calculates the chemical potentials for each element within a facet.

        Args:
            facet:
                Facet of the phase diagram.

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
        tol = PourbaixAnalyzer.numerical_tol
        all_chempots = []
        facets = self._pd.facets
        entries = self._pd.qhull_entries
        for facet in facets:
            chempots = self.get_facet_chempots(facet)
            chempots["H+"] /= -0.0591
            chempots["V"] = -chempots["V"]
            chempots["1"] = chempots["1"]
            all_chempots.append([chempots[el] for el in self._keys])

        chempots_np = np.array(all_chempots)
        chempot_ranges = collections.defaultdict(list)
        edges_are = []
        for vertex in sorted(self._pd.vertices):
            entry0 = entries[vertex]
            where_vertex = np.where(facets == vertex)
            facets_thisvertex = [facets[where_vertex[0][i]]
                                 for i in xrange(len(where_vertex[0]))]
            vertices_other = set()
            for facet in facets_thisvertex:
                for vert in facet:
                    if (vert > vertex) & (vert != vertex):
                        vertices_other.add(vert)
            for vert in sorted(vertices_other):
                entry1 = entries[vert]
                facets_of_edge = [facet for facet in facets_thisvertex
                                  if ((vert in facet) & (vertex in facet))]
                data = []
                if len(facets_of_edge) < 2:
                    edges_are.append(facets_of_edge[0])
                    continue
                for facet in facets_of_edge:
                    which_facet = facets.tolist().index(facet.tolist())
                    chempot = chempots_np[which_facet]
                    data.append([chempot[0], chempot[1]])
                sim1 = Simplex(data)
                chempot_ranges[entry0].append(sim1)
                chempot_ranges[entry1].append(sim1)
        edges = []
        combos = []
        points_on_border = collections.defaultdict(list)
        edge_vertices_set = set()

        for facet in facets:
            for combi in itertools.combinations(facet, 2):
                combos.append(combi)
        ccombo = collections.Counter(combos)
        for combi in ccombo:
            if ccombo[combi] == 1:
                edges.append(combi)
                for vertex in combi:
                    edge_vertices_set.add(vertex)
        edge_vertices = list(edge_vertices_set)
        self.edge_vertices = edge_vertices

        min_pH = min([chempot[0] for chempot in all_chempots])
        max_pH = max([chempot[0] for chempot in all_chempots])
        minV = min([chempot[1] for chempot in all_chempots])
        maxV = max([chempot[0] for chempot in all_chempots])
        self.chempot_limits = [[min_pH, max_pH], [minV, maxV]]
        xlo = min(min_pH - 5.0, -10.0)
        xhi = max(max_pH + 5.0, 20.0)
        ylo = min(minV - 3.0, -10.0)
        yhi = max(maxV + 3.0, 10.0)

        max_dist = np.sqrt((yhi - ylo) ** 2 + (xhi - xlo) ** 2)

        dist_edge = {}
        bound_coords_edge = {}
        dist_dof = {}  # Number of possible boundary choices for each edge
        coord_facet_edge = {}
        for edge in edges:
            ifacet = divmod(combos.index(edge), 3)[0]
            entry1 = entries[edge[0]]
            entry2 = entries[edge[1]]

            del_npH = (entry1.npH - entry2.npH) * (- 0.0591)
            del_nPhi = (entry1.nPhi - entry2.nPhi) * -1
            del_g0 = self._pd.qhull_data[edge[0]][2] -\
                self._pd.qhull_data[edge[1]][2]

            mu_pH_facet = all_chempots[ifacet][0]
            mu_V_facet = all_chempots[ifacet][1]
            coord_facet = np.array([mu_pH_facet, mu_V_facet])
            coord_facet_edge[edge] = coord_facet

            # Find intersection of coexistence line with bounding box
            bound_coords = []
            if del_nPhi != 0:
                Vlow = (del_g0 - del_npH * xlo) / del_nPhi
                bound_coords.append([xlo, Vlow])
                Vhigh = (del_g0 - del_npH * xhi) / del_nPhi
                bound_coords.append([xhi, Vhigh])
            else:
                bound_coords.append([mu_pH_facet, -100000000])
                bound_coords.append([mu_pH_facet, 100000000])
            if del_npH != 0:
                pHlow = (del_g0 - del_nPhi * ylo) / del_npH
                pHhigh = (del_g0 - del_nPhi * yhi) / del_npH
                bound_coords.append([pHlow, ylo])
                bound_coords.append([pHhigh, yhi])
            else:
                bound_coords.append([-100000000, mu_V_facet])
                bound_coords.append([100000000, mu_V_facet])
            dist = [np.linalg.norm(np.array(coord_facet - coord))
                    for coord in bound_coords]
            dist = (lambda x: [100000000 if x[i] > max_dist else x[i]
                               for i in xrange(len(x))])(dist)
            for border_line in itertools.combinations(xrange(len(dist)), 2):
                if (dist[border_line[0]] == 100000000) | \
                        (dist[border_line[1]] == 100000000):
                    continue
                dist_btw_pts = np.linalg.norm(
                    np.array(bound_coords[border_line[0]]) -
                    np.array(bound_coords[border_line[1]]))
                if abs(abs(dist[border_line[0]] - dist[border_line[1]])
                        - dist_btw_pts) < tol:
                    if dist[border_line[0]] < dist[border_line[1]]:
                        dist[border_line[0]] = 100000000
                    else:
                        dist[border_line[1]] = 100000000
            for entry in self._pd.stable_entries:
                for line in chempot_ranges[entry]:
                    r = [line.coords[1][0] - line.coords[0][0],
                         line.coords[1][1] - line.coords[0][1]]
                    p = [line.coords[0][0], line.coords[0][1]]
                    for i in xrange(len(bound_coords)):
                        coord = bound_coords[i]
                        s = coord - coord_facet
                        q = coord_facet
                        if abs(np.cross(r, s)) < tol:
                            break
                        else:
                            t = np.cross(q - p, s) / np.cross(r, s)
                            u = np.cross(q - p, r) / np.cross(r, s)
                            if ((t > tol) & (t - 1 < tol)) & \
                                    ((u > tol) & (u - 1 < tol)):
                                dist[i] = 100000000
            bound_coords = np.array(bound_coords)
            dist_edge[edge] = dist
            dist_dof[edge] = len(filter(lambda x: x < 10000, dist))
            bound_coords_edge[edge] = bound_coords

        dist_dof_sorted = collections.OrderedDict(
            sorted(dist_dof.items(), key=lambda t: t[1]))

        for edge in [key for key in dist_dof_sorted.keys()
                     if dist_dof_sorted[key] == 1]:
            border_coord = bound_coords_edge[edge][
                np.argmin(np.array(dist_edge[edge]))]
            coord_facet = coord_facet_edge[edge]
            line = np.array([border_coord, coord_facet])
            sim = Simplex(line)

            for vertex in edge:
                entry = entries[vertex]
                chempot_ranges[entry].append(sim)
                points_on_border[vertex].append(border_coord)

        # And now for all dof == 2
        for edge in [key for key in dist_dof_sorted.keys()
                     if dist_dof_sorted[key] == 2]:
            # 1. Evaluate g at all boundary points
            for ic in xrange(len(bound_coords_edge[edge])):
                coord = bound_coords_edge[edge][ic]
                if dist_edge[edge][ic] > 100000:
                    continue
                else:
                    g_edge = []
                    for vertex in edge_vertices:
                        entry = entries[vertex]
                        npH = entry.npH * (-0.0591)
                        nPhi = entry.nPhi * -1
                        g0 = entry.g0
                        pH = coord[0]
                        V = coord[1]
                        g_edge.append(g0 - npH * pH - nPhi * V)
                if len(np.where((np.array(g_edge) - min(g_edge) *
                       np.ones(len(g_edge))) < tol)[0]) == 2:
                    # For correct boundary coord, there will be 2 min. values
                    border_coord = coord
                    break
            coord_facet = coord_facet_edge[edge]
            line = np.array([border_coord, coord_facet])
            sim = Simplex(line)
            for vertex in edge:
                entry = entries[vertex]
                chempot_ranges[entry].append(sim)
                points_on_border[vertex].append(border_coord)

        # Now add vertical and horizontal bounding boxes
        corner_points = [[xlo, ylo], [xhi, ylo], [xhi, yhi], [xlo, yhi]]
        for vertex in edge_vertices:
            corner_pts_add = []
            for point in corner_points:
                g_others = [self.g(entry, point[0], point[1])
                            for entry in entries
                            if entry is not entries[vertex]]
                g_this = self.g(entries[vertex], point[0], point[1])
                if g_this < min(g_others):
                    corner_pts_add.append(point)
            if len(corner_pts_add) > 1:
                corner_pts_indices = [i for i in xrange(len(corner_pts_add))]
                for combi in itertools.combinations(corner_pts_indices, 2):
                    if ((abs(corner_pts_add[combi[0]][0] -
                             corner_pts_add[combi[1]][0]) < tol) |
                        (abs(corner_pts_add[combi[0]][1] -
                             corner_pts_add[combi[1]][1]) < tol)):
                        line = [corner_pts_add[combi[0]],
                                corner_pts_add[combi[1]]]
                        sim = Simplex(line)
                        chempot_ranges[entries[vertex]].append(sim)
            elif len(corner_pts_add) < 1:
                if len(points_on_border[vertex]) > 2:
                    raise StandardError("Undefined region. More than 2 points"
                                        " lie on same bounding line!")
                line = [point for point in points_on_border[vertex]]
                sim = Simplex(line)
                chempot_ranges[entries[vertex]].append(sim)
                continue

            for point in points_on_border[vertex]:
                for pt_corner in corner_pts_add:
                    if ((abs(point[0] - pt_corner[0]) < tol) |
                            (abs(point[1] - pt_corner[1]) < tol)):
                        line = [point, pt_corner]
                        sim = Simplex(line)
                        chempot_ranges[entries[vertex]].append(sim)
        chempot_ranges_cleaned = {}
        for entry in self._pd.stable_entries:
            chempot_ranges_cleaned[entry] = \
                self.check_regions(entry, chempot_ranges[entry])
        self.chempot_ranges = chempot_ranges_cleaned
        return chempot_ranges_cleaned

    def _in_facet(self, facet, entry):
        """
        Checks if a Pourbaix Entry is in a facet.

        Args:
            facet:
                facet to test.
            entry:
                Pourbaix Entry to test.
        """
        dim = len(self._keys)
        if dim > 1:
            coords = [np.array(self._pd.qhull_data[facet[i]][0:dim - 1])
                      for i in xrange(len(facet))]
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

    def check_regions(self, this_entry, lines):
        """
        Checks if a simplex region is enclosed. Computes intersections of each
        line with the other and clips any extraneous line segments. If not,
        returns null

        Args:
            this_entry:
                Pourbaix entry whose domain needs to be cleaned up
            lines:
                list of lines

        Returns:
            region:
                simplex of points
        """
        tol = PourbaixAnalyzer.numerical_tol
        coord_indx = -1
        coords_list = []
        line_coord = collections.defaultdict(list)
        new_lines = []
        for line in lines:
            if line not in new_lines:
                new_lines.append(line)
        lines = new_lines
        nlines = len(lines)
        line_indices = [i for i in xrange(nlines)]
        for line in lines:
            for coord in line.coords:
                if list(coord) not in coords_list:
                    coords_list.append(list(coord))
                    coord_indx += 1
                else:
                    coord_indx = coords_list.index(list(coord))
                line_coord[coord_indx].append(line)
        coords_count = list()
        for ic in xrange(len(coords_list)):
            coords_count.append(len(line_coord[ic]))
        line_indices_trunc = []
        for coord in [coords_list[ic] for ic in xrange(len(coords_list))
                      if coords_count[ic] == 1]:
            coord_indx = coords_list.index(coord)
            line_indices_trunc.append(lines.index(line_coord[coord_indx][0]))
        replacement_lines = list()
        replaced_lines = list()
        new_region = list()
        for combi in itertools.combinations(line_indices, 2):
            line1 = lines[combi[0]]
            line2 = lines[combi[1]]
            r = np.array([line1.coords[1][0] - line1.coords[0][0],
                          line1.coords[1][1] - line1.coords[0][1]])
            p = np.array([line1.coords[0][0], line1.coords[0][1]])
            s = np.array([line2.coords[1][0] - line2.coords[0][0],
                          line2.coords[1][1] - line2.coords[0][1]])
            q = np.array([line2.coords[0][0], line2.coords[0][1]])
            if abs(np.cross(r, s)) < tol:
                continue
            else:
                t = np.cross(q - p, s) / np.cross(r, s)
                u = np.cross(q - p, r) / np.cross(r, s)
                if ((t >= tol) & (abs(t - 1) <= tol)) & \
                        ((u >= tol) & (abs(u - 1) <= tol)):
                    if ((t < tol) | (abs(t - 1.0) < tol)) & \
                            ((u > tol) & (abs(u - 1.0) > tol)):
                        intersect_coords = q + u * s
                        intersect_coords_plus = q + u * 1.1 * s
                        intersect_coords_minus = q + u * 0.9 * s
                        line_to_be_clipped = line2
                    elif ((u < tol) | (abs(u - 1.0) < tol)) & \
                            ((t > tol) & (abs(t - 1.0) > tol)):
                        intersect_coords = p + t * r
                        intersect_coords_plus = p + t * 1.1 * r
                        intersect_coords_minus = p + t * 0.9 * r
                        line_to_be_clipped = line1
                    elif ((u < tol) | (abs(u - 1.0) < tol)) & \
                            ((t < tol) | (abs(t - 1.0) < tol)):
                        continue
                    else:
                        raise StandardError("Lines intersect, "
                                            "domain ill-defined")
                    g_other_plus = [self.g(entry, intersect_coords_plus[0],
                                           intersect_coords_plus[1])
                                    for entry in self._get_edge_entries()
                                    if entry is not this_entry]
                    g_this_plus = self.g(this_entry, intersect_coords_plus[0],
                                         intersect_coords_plus[1])
                    g_other_minus = [self.g(entry, intersect_coords_minus[0],
                                            intersect_coords_minus[1])
                                     for entry in self._get_edge_entries()
                                     if entry is not this_entry]
                    g_this_minus = self.g(this_entry, intersect_coords_minus[0],
                                          intersect_coords_minus[1])
                    line_new = []
                    if g_this_plus < min(g_other_plus):
                        line_new.append([line_to_be_clipped.coords[1][0],
                                         line_to_be_clipped.coords[1][1]])
                        line_new.append(intersect_coords)
                    elif g_this_minus < min(g_other_minus):
                        line_new.append([line_to_be_clipped.coords[0][0],
                                         line_to_be_clipped.coords[0][1]])
                        line_new.append(intersect_coords)
                    sim = Simplex(line_new)
                    replacement_lines.append(sim)
                    replaced_lines.append(line_to_be_clipped)
        for line in [line for line in lines if line not in replaced_lines]:
            new_region.append(line)
        for line in replacement_lines:
            new_region.append(line)
        return new_region

    def g(self, entry, pH, V):
        """
        Get free energy for a given pH, and V.
        """
        g0 = entry.g0
        npH = -entry.npH * 0.0591
        nPhi = -entry.nPhi
        return g0 - npH * pH - nPhi * V

    def _get_edge_entries(self):
        """
        Get all entries on edges of convex hull
        """
        entries = []
        for vertex in self.edge_vertices:
            entries.append(self._pd.qhull_entries[vertex])
        return entries

    def get_decomposition(self, entry):
        """
        Provides the decomposition at a particular composition

        Args:
            comp:
                A composition

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
        for i in xrange(len(decompamts)):
            if abs(decompamts[i][0]) > PourbaixAnalyzer.numerical_tol:
                decomp[self._pd.qhull_entries[facet[i]]] = decompamts[i][0]
        return decomp

    def get_decomp_and_e_above_hull(self, entry):
        """
        Provides the decomposition and energy above convex hull for an entry

        Args:
            entry:
                A PourbaixEntry

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
            entry - A PourbaixEntry object

        Returns:
            Energy above convex hull of entry. Stable entries should have
            energy above hull of 0.
        """
        return self.get_decomp_and_e_above_hull(entry)[1]