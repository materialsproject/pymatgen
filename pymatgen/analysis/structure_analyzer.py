#!/usr/bin/env python

"""
This module provides classes to perform topological analyses of structures.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Geoffroy Hautier"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import math
import numpy as np
import itertools
import collections

from pyhull.voronoi import VoronoiTess


class VoronoiCoordFinder(object):
    """
    Uses a Voronoi algorithm to determine the coordination for each site in a
    structure.
    """

    """Radius in Angstrom cutoff to look for coordinating atoms"""
    default_cutoff = 10.0

    def __init__(self, structure, target=None):
        """
        Args:
            structure:
                Input structure
            target:
                A list of target species to determine coordination for.
        """
        self._structure = structure
        if target is None:
            self._target = structure.composition.elements
        else:
            self._target = target

    def get_voronoi_polyhedra(self, n):
        """
        Gives a weighted polyhedra around a site. This uses the voronoi
        construction with solid angle weights.
        See ref: A Proposed Rigorous Definition of Coordination Number,
        M. O'Keeffe, Acta Cryst. (1979). A35, 772-775

        Args:
            n:
                site index

        Returns:
            A dictionary of sites sharing a common Voronoi facet with the site
            n and their solid angle weights
        """

        localtarget = self._target
        center = self._structure[n]
        neighbors = self._structure.get_sites_in_sphere(center.coords,
                                            VoronoiCoordFinder.default_cutoff)
        neighbors = [i[0] for i in sorted(neighbors, key=lambda s: s[1])]
        qvoronoi_input = [s.coords for s in neighbors]
        voro = VoronoiTess(qvoronoi_input)
        all_vertices = voro.vertices

        results = {}
        for nn, vind in voro.ridges.items():
            if 0 in nn:
                if 0 in vind:
                    raise RuntimeError("This structure is pathological,"
                                       " infinite vertex in the voronoi "
                                       "construction")

                facets = [all_vertices[i] for i in vind]
                results[neighbors[nn[1]]] = solid_angle(center.coords, facets)

        maxangle = max(results.values())

        resultweighted = {}
        for nn, angle in results.items():
            if nn.specie in localtarget:
                resultweighted[nn] = angle / maxangle

        return resultweighted

    def get_coordination_number(self, n):
        """
        Returns the coordination number of site with index n.

        Args:
            n:
                site index
        """
        return sum(self.get_voronoi_polyhedra(n).values())

    def get_coordinated_sites(self, n, tol=0, target=None):
        """
        Returns the sites that are in the coordination radius of site with
        index n.

        Args:
            n:
                Site number.
            tol:
                Weight tolerance to determine if a particular pair is
                considered a neighbor.
            Target:
                Target element

        Returns:
            Sites coordinating input site.
        """
        coordinated_sites = []
        for site, weight in self.get_voronoi_polyhedra(n).items():
            if weight > tol and (target == None or site.specie == target):
                coordinated_sites.append(site)
        return coordinated_sites


class RelaxationAnalyzer(object):
    """
    This class analyzes the relaxation in a calculation.
    """

    def __init__(self, initial_structure, final_structure):
        """
        Please note that the input and final structures should have the same
        ordering of sites. This is typically the case for most computational
        codes.

        Args:
            initial_structure:
                Initial input structure to calculation.
            final_structure:
                Final output structure from calculation.
        """
        if final_structure.formula != initial_structure.formula:
            raise ValueError("Initial and final structures have different " +
                             "formulas!")
        self.initial = initial_structure
        self.final = final_structure

    def get_percentage_volume_change(self):
        """
        Returns the percentage volume change.

        Returns:
            Volume change in percentage, e.g., 0.055 implies a 5.5% increase.
        """
        initial_vol = self.initial.lattice.volume
        final_vol = self.final.lattice.volume
        return final_vol / initial_vol - 1

    def get_percentage_lattice_parameter_changes(self):
        """
        Returns the percentage lattice parameter changes.

        Returns:
            A dict of the percentage change in lattice parameter, e.g.,
            {'a': 0.012, 'b': 0.021, 'c': -0.031} implies a change of 1.2%,
            2.1% and -3.1% in the a, b and c lattice parameters respectively.
        """
        initial_latt = self.initial.lattice
        final_latt = self.final.lattice
        d = {l: getattr(final_latt, l) / getattr(initial_latt, l) - 1
             for l in ["a", "b", "c"]}
        return d

    def get_percentage_bond_dist_changes(self, max_radius=3.0):
        """
        Returns the percentage bond distance changes for each site up to a
        maximum radius for nearest neighbors.

        Args:
            max_radius:
                Maximum radius to search for nearest neighbors. This radius is
                applied to the initial structure, not the final structure.

        Returns:
            Bond distance changes as a dict of dicts. E.g.,
            {index1: {index2: 0.011, ...}}. For economy of representation, the
            index1 is always less than index2, i.e., since bonding between
            site1 and siten is the same as bonding between siten and site1,
            there is no reason to duplicate the information or computation.
        """
        data = collections.defaultdict(dict)
        for inds in itertools.combinations(xrange(len(self.initial)), 2):
            (i, j) = sorted(inds)
            initial_dist = self.initial[i].distance(self.initial[j])
            if initial_dist < max_radius:
                final_dist = self.final[i].distance(self.final[j])
                data[i][j] = final_dist / initial_dist - 1
        return data


def solid_angle(center, coords):
    """
    Helper method to calculate the solid angle of a set of coords from the
    center.

    Args:
        center:
            Center to measure solid angle from.
        coords:
            List of coords to determine solid angle.

    Returns:
        The solid angle.
    """
    o = np.array(center)
    r = [np.array(c) - o for c in coords]
    r.append(r[0])
    n = [np.cross(r[i + 1], r[i]) for i in range(len(r) - 1)]
    n.append(np.cross(r[1], r[0]))
    phi = sum([math.acos(-np.dot(n[i], n[i + 1])
                         / (np.linalg.norm(n[i]) * np.linalg.norm(n[i + 1])))
               for i in range(len(n) - 1)])
    return phi + (3 - len(r)) * math.pi


def contains_peroxide(structure, relative_cutoff=1.2):
    """
    Determines if a structure contains peroxide anions.

    Args:
        structure:
            Input structure.
        relative_cutoff:
            The peroxide bond distance is 1.49 Angstrom. Relative_cutoff * 1.49
            stipulates the maximum distance two O atoms must be to each other
            to be considered a peroxide.

    Returns:
        Boolean indicating if structure contains a peroxide anion.
    """
    max_dist = relative_cutoff * 1.49
    o_sites = []
    for site in structure:
        syms = [sp.symbol for sp in site.species_and_occu.keys()]
        if "O" in syms:
            o_sites.append(site)

    for i, j in itertools.combinations(o_sites, 2):
        if i.distance(j) < max_dist:
            return True

    return False
