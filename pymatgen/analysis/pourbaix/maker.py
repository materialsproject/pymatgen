#!/usr/bin/env python

"""
Module containing analysis classes which compute a pourbaix diagram given a target
compound/element.
TODO :1 Extension to compute alloy/compound/multiple-element pourbaix diagrams
"""

from __future__ import division

__author__ = "Sai Jayaraman"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.0"
__maintainer__ = "Sai Jayaraman"
__email__ = "sjayaram@mit.edu"
__status__ = "Development"
__date__ = "Nov 1, 2012"

PREFAC = 0.0591

import math
import logging
import numpy as np

MU_H2O = -2.4583
from pyhull.convex_hull import ConvexHull

logger = logging.getLogger(__name__)


class PourbaixDiagram(object):
    """
    Class to create a Pourbaix diagram from entries
    """
    def __init__(self, entries):
        """
        Args:
            entries:
                Entries list containing both Solids and Ions
        """
        self._all_entries = entries
        self._solid_entries = list()
        self._ion_entries = list()
        for entry in entries:
            if entry.phase_type == "Solid":
                self._solid_entries.append(entry)
            elif entry.phase_type == "Ion":
                self._ion_entries.append(entry)
            else:
                raise StandardError("Incorrect Phase type - needs to be Ion/Solid")
        self._make_pourbaixdiagram()

    def _create_conv_hull_data(self):
        """
        Make data conducive to convex hull generator.
        """
        entries_to_process = list()
        for entry in self._all_entries:
            entry.normalize(entry.normalization_factor)
            entry.g0_add(- MU_H2O * entry.nH2O + PREFAC * math.log10(entry.conc))
            entries_to_process.append(entry)
        self._qhull_entries = entries_to_process
        return self._process_conv_hull_data(entries_to_process)

    def _process_conv_hull_data(self, entries_to_process):
        """
        From a sequence of ion+solid entries, generate the necessary data
        for generation of the convex hull.
        """
        data = []
        for entry in entries_to_process:
            row = [entry.npH, entry.nPhi, entry.g0]
            data.append(row)
        return data

    def _make_pourbaixdiagram(self):
        """
        Calculates entries on the convex hull in the dual space.
        """
        stable_entries = set()
        self._qhull_data = self._create_conv_hull_data()
        dim = len(self._qhull_data[0])
        if len(self._qhull_data) < dim:
            raise StandardError("Can only do elements with atleast 3 entries for now")
        if len(self._qhull_data) == dim:
            self._facets = [range(dim)]
            vertices = set()
            for facet in self._facets:
                for vertex in facet:
                    vertices.add(vertex)

        else:
            facets_qhull = np.array(ConvexHull(self._qhull_data).vertices)
            self._facets = np.sort(np.array(facets_qhull))
            logger.debug("Final facets are\n{}".format(self._facets))

            logger.debug("Removing vertical facets...")
            vert_facets_removed = list()
            for facet in self._facets:
                facetmatrix = np.zeros((len(facet), len(facet)))
                count = 0
                for vertex in facet:
                    facetmatrix[count] = np.array(self._qhull_data[vertex])
                    facetmatrix[count, dim - 1] = 1
                    count += 1
                if abs(np.linalg.det(facetmatrix)) > 1e-8:
                    vert_facets_removed.append(facet)
                else:
                    print "removed facet", facet
                    logger.debug("Removing vertical facet : {}".format(facet))

            logger.debug("Removing UCH facets by eliminating normal.z >0 ...")

            # Find center of hull
            vertices = set()
            for facet in vert_facets_removed:
                for vertex in facet:
                    vertices.add(vertex)
            c = [0.0, 0.0, 0.0]
            c[0] = np.average([self._qhull_data[vertex][0] for vertex in vertices])
            c[1] = np.average([self._qhull_data[vertex][1] for vertex in vertices])
            c[2] = np.average([self._qhull_data[vertex][2] for vertex in vertices])
            # Shift origin to c
            new_qhull_data = np.array(self._qhull_data)
            for vertex in vertices:
                new_qhull_data[vertex] -= c

            # For each facet, find normal n, find dot product with P, and check if this is -ve
            final_facets = list()
            for facet in vert_facets_removed:
                a = new_qhull_data[facet[1]] - new_qhull_data[facet[0]]
                b = new_qhull_data[facet[2]] - new_qhull_data[facet[0]]
                n = np.cross(a, b)
                val = np.dot(n, new_qhull_data[facet[0]])
                if val < 0:
                    n = -n
                if n[2] <= 0:
                    final_facets.append(facet)
                else:
                    print "removed UCH facet", facet
                    logger.debug("Removing UCH facet : {}".format(facet))

            self._facets = final_facets

        for facet in self._facets:
            for vertex in facet:
                stable_entries.add(self._qhull_entries[vertex])
        self._stable_entries = stable_entries
        self._vertices = vertices

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
        composition data and energy per atom values created from qhull_entries.
        """
        return self._qhull_data

    @property
    def qhull_entries(self):
        """
        Return qhull entries
        """
        return self._qhull_entries

    @property
    def stable_entries(self):
        """
        Returns the stable entries in the phase diagram.
        """
        return self._stable_entries

    @property
    def all_entries(self):
        """
        Return all entries
        """
        return self._all_entries

    @property
    def vertices(self):
        """
        Return vertices of the convex hull
        """
        return self._vertices
