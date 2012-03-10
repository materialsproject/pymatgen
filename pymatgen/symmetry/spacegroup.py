#!/usr/bin/env python

'''
This module implements a basic Spacegroup class
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 9, 2012"


class Spacegroup(object):
    """
    Represents a space group, which is a collection of symmetry operations
    """
    
    def __init__(self, int_symbol, int_number, symmops):
        self._symbol = int_symbol
        self._number = int_number
        self._symmops = symmops
    
    @property
    def international_symbol(self):
        return self._symbol
    
    @property
    def international_number(self):
        return self._number
    
    
    
    