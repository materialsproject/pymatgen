#!/usr/bin/env python

'''
Created on May 15, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "May 15, 2012"

import itertools

import numpy as np


class Simplex(object):
    """
    A generalized simplex object. See http://en.wikipedia.org/wiki/Simplex.
    """

    def __init__(self, coords):
        """
        Args:
            coords:
                Coords of the vertices of the simplex. E.g.,
                [[1,2,3],[2,4,5],[6,7,8]]. Note that the length of each coord
                must equal to the number of vertices -1 
        """
        self.simplex_dim = len(coords) - 1
        self.space_dim = len(coords[0])
        for c in coords:
            if len(c) != self.space_dim:
                raise ValueError("All coords must have the same space dimension.")
        self.coords = coords

    def __eq__(self, other):
        for p in itertools.permutations(self.coords):
            if np.allclose(p, other.coords):
                return True
        return False

    def __hash__(self):
        return 7

    def get_plot_coords(self):
        coords = np.array(self.coords)
        return coords.transpose()
