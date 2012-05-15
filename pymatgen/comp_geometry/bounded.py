#!/usr/bin/env python

'''
This module defines various bounded computational geometry classes. E.g., a
bounded line is a line that has a start and end point.

Note that this module is work in progress.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "May 15, 2012"

import numpy as np


class BoundedLine(object):
    """
    Defines a bounded line.
    
    TODO:
        Move this code to a generalized com
    """

    def __init__(self, coord1, coord2):
        """
        Args:
            coord1, coord2:
                Start and end points of the bounded line. It does not matter
                which one is defined as the start or end. The bounded line has
                no directionality.
        """
        self.coord1 = coord1
        self.coord2 = coord2

    def __eq__(self, other):
        if np.allclose(self.coord1, other.coord1) and np.allclose(self.coord2, other.coord2):
            return True
        if np.allclose(self.coord1, other.coord2) and np.allclose(self.coord2, other.coord1):
            return True
        return False

    def __hash__(self):
        return 7

    def get_plot_coords(self):
        return [[self.coord1[i], self.coord2[i]] for i in xrange(len(self.coord1))]
