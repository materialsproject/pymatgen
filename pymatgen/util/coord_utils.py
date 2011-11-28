#!/usr/bin/env python

'''
Utilities for manipulating coordinates or list of coordinates
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 27, 2011"

import numpy as np

def in_coord_list(coord_list, coord):
    for test in coord_list:
        if np.allclose(test, coord):
            return True
    return False