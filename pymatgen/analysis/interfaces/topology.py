# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals
from __future__ import absolute_import
from fractions import gcd
import numpy as np
from monty.serialization import loadfn
from pymatgen.core.surface import get_symmetrically_distinct_miller_indices
from pymatgen.core.surface import SlabGenerator
from monty.json import MSONable

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__credits__ = "Hong Ding"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Feb 3, 2016"

"""
This module provides a class used to identify matched super lattice solid-solid
interfaces using the algorithm of Zur and McGill 1984
"""
