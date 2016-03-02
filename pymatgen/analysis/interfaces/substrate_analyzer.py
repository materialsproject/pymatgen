# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module provides classes to identify optimal substrates for film growth
"""

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Production"
__date__ = "Feb, 2016"
import numpy as np
import scipy.constants as const

from monty.json import MSONable
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.analysis.interfaces.topology import ZSLGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.elasticity.elastic import ElasticTensor
from pymatgen.analysis.elasticity.strain import Deformation



class SubstrateAnalyzer(MSONable):



    def __init__(self, film, elasticity_tensor = None):

        self.film = film
        self.elasticity_tensor = elasticity_tensor
        self.susbtrates = {}
        self.matches = {}


    def calculate(self,substrate):


        z = ZSLGenerator(self.film,substrate)

        for match in z.generate():

            energy = self.calculate_elastic_energy(match)




    def calculate_elastic_energy(self,match):

        if self.elasticity_tensor is None:
            return 9999

        film_matrix = match.film_sl_vectors
        film_matrix.append(np.cross(film_matrix[0],film_matrix[1]))

        substrate_matrix = match.substrate_sl_vectors
        substrate_matrix.append(np.cross(substrate_matrix[0],
            substrate_matrix[1]))

        transform_matrix = np.linalg.solve(film_matrix,substrate_matrix)

        energy_density = self.elasticity_tensor.energy_density(
            Deformation(transform_matrix).green_lagrange_strain)

        return self.film.volume*energy_density
