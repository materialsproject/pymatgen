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

from pymatgen.command_line.qhull_caller import qvoronoi, qvertex


class VoronoiCoordFinder:
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
        if target == None:
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
        closest = qvoronoi(qvoronoi_input)
        all_vertices = qvertex(qvoronoi_input)

        result = []
        angle = []
        for i in range(len(all_vertices)):
            if closest[i][1] == 0:
                facets = []
                result.append(neighbors[closest[i][2]])
                for j in range(3, len(closest[i])):
                    if closest[i][j] == 0:
                        raise RuntimeError("This structure is pathological,"
                                           " infinite vertex in the voronoi "
                                           "construction")
                    facets.append(all_vertices[closest[i][j] - 1])

                angle.append(solid_angle(center.coords, facets))

        maxangle = max(angle)

        resultweighted = {}
        for j in range(len(angle)):
            if result[j].specie in localtarget:
                resultweighted[result[j]] = angle[j] / maxangle

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
