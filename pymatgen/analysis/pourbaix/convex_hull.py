#!/usr/bin/env python

"""
This module implements a ConvexHull class. Note that this passes QJ as a
command line option to qconvex. This can handle co-planar points, which are
generated in pourbaix.analyzer.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Nov 19 2012"

from pyhull import qconvex
from pyhull.simplex import Simplex


class ConvexHull(object):
    """
    Convex hull for a set of points.

    .. attribute: dim

        Dimension of the points.

    .. attribute: points

        Original points supplied.

    .. attribute: vertices

        The vertices as a list of list of integer indices. E.g., [[0, 2], [1,
        0], [2, 3], [3, 1]]
    """

    def __init__(self, points):
        """
        Args:
            points:
                All the points as a sequence of sequences. e.g., [[-0.5, -0.5],
                [-0.5, 0.5], [0.5, -0.5], [0.5, 0.5]]
        """
        self.points = points
        dim = map(len, self.points)
        if max(dim) != min(dim):
            raise ValueError("Input points must all have the same dimension!")
        self.dim = dim[0]
        output = qconvex("i QJ", points)
        output.pop(0)
        self.vertices = [map(int, row.strip().split()) for row in output]

    @property
    def simplices(self):
        """
        Returns the simplices of the convex hull.
        """
        return [Simplex([self.points[i] for i in v]) for v in self.vertices]
