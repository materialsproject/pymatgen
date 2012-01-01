#!/usr/bin/env python

"""
This module defines some useful design patterns.
"""

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="Sep 23, 2011"

class Enum(set):
    """
    Creates an enum out of a set.
    """
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError