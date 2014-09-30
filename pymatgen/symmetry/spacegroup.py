# coding: utf-8

from __future__ import division, unicode_literals

"""
This module implements a basic Spacegroup class
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 9, 2012"

from pymatgen.symmetry.analyzer import SpacegroupOperations
from monty.dev import deprecated


@deprecated(replacement=SpacegroupOperations)
class Spacegroup(SpacegroupOperations):
    pass