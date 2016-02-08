# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals
from __future__ import absolute_import

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__credits__ = "Hong Ding"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Feb 3, 2016"

"""
This module provides a class used to determine the optimal substrates for
the epitaxial growth a given film.
"""

from pymatgen.analysis.substrate.cal import structure_interface_matching

from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen.core.structure import Structure

from pymatgen.matproj.rest import MPRester

from fractions import gcd


import numpy as np

class ZurSuperLatticeGenerator(object):
    """
    This class generate interface super lattices based on the methodology
    of lattice vector matching for heterostructural interfaces proposed by
    Zur and McGill:
    Journal of Applied Physics 55 (1984), 378 ; doi: 10.1063/1.333084

    The process of generating all possible interaces is as such:

    1.) Generate all slabs for the film and substrate for different orientations
    2.) For each film/substrate orientation pair:
        1.) Reduce lattice vectors and calculate area
        2.) Generate all super lattice transformations
        3.) For each superlattice set:
            1.) Reduce super lattice vectors
            2.) Check length and angle between film and substrate super lattice vectors to determine if the super lattices are the same and therefore coincident


    """

    def __init__(self,substrate,film):
        self.substrate = substrate
        self.film = film

    def generate(self,film_miller = None, substrate_miller = None):
        """
            Generates the film/substrate combinations for either set miller indicies or all possible miller indices up to a max miller index
        """

        structures = structure_interface_matching(self.film,self.substrate)

        return structures
