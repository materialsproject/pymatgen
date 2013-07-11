#!/usr/bin/env python
"""
Useful coordinate geometry routines
"""

from __future__ import division

__author__ = "Sai Jayaraman" 
__copyright__ = "Copyright 2013, The Materials Project" 
__version__ = "1.0"
__maintainer__ = "Sai Jayaraman"
__email__ = "sjayaram@mit.edu"
__status__ = "Development"
__date__ = "Jul 10, 2013"

import numpy as np
#from pyhull.simplex import Simplex

numerical_tol = 1e-8


def intersect(l1, l2, tol=numerical_tol):
    """
    Intersection of two line segments.
    args:
        l1, l2: Simplex of 2 line segments whose intersection is desired.
    returns:
        parameters t, and u, assuming the vectors denoting the line segments
        are p + tr and q + us
    """
    
    try:
        r = [l1.coords[1][0] - l1.coords[0][0], \
         l1.coords[1][1] - l1.coords[0][1]]
        p = [l1.coords[0][0], l1.coords[0][1]]
    except AttributeError:
        r = [l1[1][0] - l1[0][0], \
         l1[1][1] - l1[0][1]]
        p = [l1[0][0], l1[0][1]]
    try:
        s = [l2.coords[1][0] - l2.coords[0][0], \
         l2.coords[1][1] - l2.coords[0][1]]
        q = [l2.coords[0][0], l2.coords[0][1]]
    except AttributeError:
        s = [l2[1][0] - l2[0][0], \
         l2[1][1] - l2[0][1]]
        q = [l2[0][0], l2[0][1]]
    p = np.array(p)
    r = np.array(r)
    q = np.array(q)
    s = np.array(s)
    if abs(np.cross(r, s)) < tol:
        if abs(np.cross(q - p, r)) < tol:
            print "The cross product", q - p, r, np.cross(q-p,r), l1.coords
            return None, None
        else:
            return None
    else:
        t = np.cross(q - p, s) / np.cross(r, s)
        u = np.cross(q - p, r) / np.cross(r, s)
    return (t, u), (np.array(p) + t * np.array(r))
