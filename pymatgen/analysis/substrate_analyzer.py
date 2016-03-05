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

from fractions import gcd
import numpy as np

from monty.json import MSONable
from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.core.surface import get_symmetrically_distinct_miller_indices
from pymatgen.core.surface import SlabGenerator


class ZSLMatch(MSONable):
    def __init__(self, film_miller, substrate_miller, film_sl_vectors,
                 substrate_sl_vectors, film_vectors, susbtrate_vectors,
                 match_area):
        self.film_miller = np.array(film_miller)
        self.substrate_miller = np.array(substrate_miller)
        self.film_sl_vectors = list(film_sl_vectors)
        self.substrate_sl_vectors = list(substrate_sl_vectors)
        self.film_vectors = list(film_vectors)
        self.substrate_vectors = list(susbtrate_vectors)
        self.match_area = float(match_area)

    def as_dict(self):
        """
        Returns dict which contains ZSL match
        """
        d = {}
        d["film miller"] = self.film_miller
        d["susbtrate miller"] = self.substrate_miller
        d["film super lattice vectors"] = self.film_sl_vectors
        d["substrate super lattice vectors"] = self.substrate_sl_vectors
        d["matching area"] = self.match_area
        d["film vectors"] = self.film_vectors
        d["susbtrate vectors"] = self.substrate_vectors

        return d


class ZSLGenerator(object):
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
            2.) Check length and angle between film and substrate super lattice
                vectors to determine if the super lattices are the same and
                therefore coincident
    """

    def __init__(self, film, substrate, max_area_ratio_tol=0.09,
                 max_area=400, film_max_miller=1, substrate_max_miller=1,
                 max_length_tol=0.03, max_angle_tol=0.01):
        """
        Intialize a Zur Super Lattice Generator for a specific film and
            substrate
        """
        self.substrate = substrate
        self.film = film
        self.max_area_ratio_tol = max_area_ratio_tol
        self.max_area = max_area
        self.film_max_miller = film_max_miller
        self.substrate_max_miller = substrate_max_miller
        self.max_length_tol = max_length_tol
        self.max_angle_tol = max_angle_tol

    def rel_strain(self, vec1, vec2):
        """
        Calculate relative strain between two vectors
        """
        return fast_norm(vec2) / fast_norm(vec1) - 1

    def rel_angle(self, vec_set1, vec_set2):
        """
        Calculate the relative angle between two vectors
        """
        return vec_angle(vec_set2[0], vec_set2[1]) / vec_angle(
            vec_set1[0], vec_set1[1]) - 1

    def is_same_vectors(self, vec_set1, vec_set2):
        """
        Determine if two sets of vectors are the same within length and
        angle tolerances
        """
        if (np.absolute(self.rel_strain(vec_set1[0], vec_set2[0])) >
                self.max_length_tol):
            return False
        elif (np.absolute(self.rel_strain(vec_set1[1], vec_set2[1])) >
                  self.max_length_tol):
            return False
        elif (np.absolute(self.rel_angle(vec_set1, vec_set2)) >
                  self.max_angle_tol):
            return False
        else:
            return True

    def generate_sl_transformation(self, area_multiple):
        """
        Generates the transformation matricies that convert a set of 2D
        vectors into a super lattice of integer area multiple as proven
        in Cassels:

        Cassels, John William Scott. An introduction to the geometry of
        numbers. Springer Science & Business Media, 2012.

        Args:
            area_multiple: integer multiple of unit cell area for super
            lattice area

        Returns:
            matrix_list: transformation matricies to covert unit vectors to
            super lattice vectors
        """

        for i in get_factors(area_multiple):
            for j in range(area_multiple // i):
                yield np.matrix(((i, j), (0, area_multiple / i)))

    def generate_sl_transformations(self, film_area, substrate_area):
        """
        Generates transformation sets for film/substrate pair given the
        area of the unit cell area for the film and substrate. The
        transformation sets map the film and substrate unit cells to super
        lattices with a maximum area

        Args:
            film_area: the unit cell area for the film
            substrate_area: the unit cell area for the substrate

        Returns:
            transformation_sets: a set of transformation_sets defined as:
                1.) the (i,j) pair corresponding to the integer multiple of
                the film area (i) and substrate area (j) that makes the two
                equal within tolerance
                2.) the transformation matricies for the film to create a
                super lattice of area i*film area
                3.) the tranformation matricies for the substrate to create
                a super lattice of area j*film area
        """

        for i in range(1, int(self.max_area / film_area)):
            for j in range(1, int(self.max_area / substrate_area)):
                if (gcd(i, j) == 1 and
                            np.absolute(film_area / substrate_area - float(
                                j) / i) <
                            self.max_area_ratio_tol):
                    yield [(i, j), self.generate_sl_transformation(i),
                           self.generate_sl_transformation(j)]

    def check_transformations(self, transformation_sets, film_vectors,
                              substrate_vectors):
        """
        Applies the transformation_sets to the film and substrate vectors
        to generate super-lattices and checks if they matches.
        Returns all matching vectors sets.
        """

        for [ij_pair, film_transformations, substrate_transformations] in \
                transformation_sets:

            films = []
            substrates = []
            # Apply transformations and reduce using Zur reduce methodology
            for f in film_transformations:
                films.append(reduce_vectors(*np.squeeze(np.asarray(
                    f * film_vectors))))
            for s in substrate_transformations:
                substrates.append(reduce_vectors(*np.squeeze(np.asarray(
                    s * substrate_vectors))))
            # Check if equivelant super lattices
            for f in films:
                for s in substrates:
                    if self.is_same_vectors(f, s):
                        yield [f, s]

    def generate_slabs(self, film_millers, substrate_millers):
        """
        Generates the film/substrate slab combinations for a set of given
        miller indicies
        """

        for f in film_millers:
            film_slab = SlabGenerator(self.film, f, 20, 15,
                                      primitive=False).get_slab()
            film_vectors = reduce_vectors(film_slab.lattice_vectors()[0],
                                               film_slab.lattice_vectors()[1])
            film_area = vec_area(*film_vectors)

            for s in substrate_millers:
                substrate_slab = SlabGenerator(self.substrate, s, 20, 15,
                                               primitive=False).get_slab()
                substrate_vectors = reduce_vectors(
                    substrate_slab.lattice_vectors()[0],
                    substrate_slab.lattice_vectors()[1])
                substrate_area = vec_area(*substrate_vectors)

                yield [film_area, substrate_area, film_vectors,
                       substrate_vectors, f, s]

    def generate(self, film_millers=None, substrate_millers=None):
        """
        Generates the film/substrate combinations for either set miller
        indicies or all possible miller indices up to a max miller index
        """

        # Generate miller indicies if none specified for film
        if film_millers is None:
            film_millers = sorted(get_symmetrically_distinct_miller_indices(
                self.film, self.film_max_miller))

        # Generate miller indicies if none specified for substrate
        if substrate_millers is None:
            substrate_millers = sorted(
                get_symmetrically_distinct_miller_indices(self.substrate,
                                                          self.substrate_max_miller))

        # Check each miller index combination
        for [film_area, substrate_area, film_vectors, substrate_vectors,
             film_miller, substrate_miller] in self.generate_slabs(film_millers,
                                                                   substrate_millers):
            # Generate all super lattice comnbinations for a given set of miller
            # indicies
            transformations = self.generate_sl_transformations(
                film_area, substrate_area)
            # Check each super-lattice pair to see if they match
            for match in self.check_transformations(transformations,
                                                    film_vectors,
                                                    substrate_vectors):
                # Yield the match area, the miller indicies,
                yield ZSLMatch(film_miller, substrate_miller, match[0],
                               match[1], film_vectors, substrate_vectors,
                               vec_area(*match[0]))


class SubstrateAnalyzer(MSONable):
    def __init__(self, film, elasticity_tensor=None):

        self.film = film
        self.elasticity_tensor = elasticity_tensor

    def calculate(self, substrate):

        z = ZSLGenerator(self.film, substrate)

        for match in z.generate():
            d = match.as_dict()
            energy = self.calculate_3D_elastic_energy(match)
            d["elastic energy"] = energy
            yield d

    def calculate_3D_elastic_energy(self, match):
        """
            Calculates the multi-plane elastic energy
        """
        if self.elasticity_tensor is None:
            return 9999

        # Generate 3D lattice vectors for film super lattice
        film_matrix = list(match.film_sl_vectors)
        film_matrix.append(np.cross(film_matrix[0], film_matrix[1]))

        # Generate 3D lattice vectors for substrate super lattice
        # Out of place substrate super lattice has to be same length as
        # Film out of plane vector to ensure no extra deformation in that
        # direction
        substrate_matrix = list(match.substrate_sl_vectors)
        temp_sub = np.cross(substrate_matrix[0], substrate_matrix[1])
        temp_sub = temp_sub * np.linalg.norm(film_matrix[2]) / np.linalg.norm(
            temp_sub)
        substrate_matrix.append(temp_sub)

        transform_matrix = np.transpose(np.linalg.solve(film_matrix,
                                                        substrate_matrix))

        dfm = Deformation(transform_matrix)

        energy_density = self.elasticity_tensor.energy_density(
            dfm.green_lagrange_strain)

        return self.film.volume * energy_density / len(self.film.sites)


def fast_norm(a):
    """
    Much faster variant of numpy linalg norm
    """
    return np.sqrt(np.dot(a, a))


def vec_angle(a, b):
    """
    Calculate angle between two vectors
    """
    cosang = np.dot(a, b)
    sinang = fast_norm(np.cross(a, b))
    return np.arctan2(sinang, cosang)


def vec_area(a, b):
    """
    Area of lattice plane defined by two vectors
    """
    return fast_norm(np.cross(a, b))


def reduce_vectors(a, b):
    """
    Generate independent and unique basis vectors based on the
    methodology of Zur and McGill
    """
    if np.dot(a, b) < 0:
        return reduce_vectors(a, -b)

    if fast_norm(a) > fast_norm(b):
        return reduce_vectors(b, a)

    if fast_norm(b) > fast_norm(np.add(b, a)):
        return reduce_vectors(a, np.add(b, a))

    if fast_norm(b) > fast_norm(np.subtract(b, a)):
        return reduce_vectors(a, np.subtract(b, a))

    return [a, b]


def get_factors(n):
    """
    Generate all factors of n
    """
    for x in range(1, n + 1):
        if n % x == 0:
            yield x