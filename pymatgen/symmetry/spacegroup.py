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

import warnings

from pymatgen.symmetry.analyzer import SpacegroupOperations
from monty.dev import deprecated


warnings.warn("The spacegroup module has been deprecated. All classes and "
              "functions have been moved to pymatgen.symmetry.analyzer. "
              "Please change your code. If you are looking for point group "
              "representations, check out pymatgen.symmetry.groups. This "
              "module will be removed in pmg v3.1.")

@deprecated(replacement=SpacegroupOperations,
            message="This class will be removed in pmg v3.1.")
class Spacegroup(SpacegroupOperations):
    pass