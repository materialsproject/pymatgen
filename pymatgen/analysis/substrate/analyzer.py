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
            2.) Check length and angle between film and substrate super lattice
                vectors to determine if the super lattices are the same and
                therefore coincident
    """

    def __init__(self, film, substrate, criteria=None):
        """
            Intialize a Zur Super Lattice Generator for a specific film and
                substrate
        """
        self.substrate = substrate
        self.film = film

        if criteria is None:
            self.criteria = loadfn('matching_criteria.json')
        else:
            self.criteria = criteria


    def reduce_vectors(self, a, b):
        """
            Generate independent and unique basis vectors based on the
            methodology of Zur and McGill
        """
        if np.dot(a, b) < 0:
            b = -b
            return self.reduce_vectors(a, b)

        if np.linalg.norm(a) > np.linalg.norm(b):
            a, b = b, a
            return self.reduce_vectors(a, b)

        if np.linalg.norm(b) > np.linalg.norm(np.add(b, a)):
            b = np.add(b, a)
            return self.reduce_vectors(a, b)

        if np.linalg.norm(b) > np.linalg.norm(np.subtract(b, a)):
            b = np.subtract(b, a)
            return self.reduce_vectors(a, b)

        return [a, b]

    def area(self, a, b):
        """
            Area of lattice plane defined by two vectors
        """
        return np.linalg.norm(np.cross(a, b))

    def factor(self, n):
        """
            Generate all factors of n
        """
        for x in range(1, n+1):
            if n % x == 0:
                yield x

    def vec_angle(self, a, b):
        """
            Calculate angle between two vectors
        """
        cosang = np.dot(a, b)
        sinang = np.linalg.norm(np.cross(a, b))
        return np.arctan2(sinang, cosang)

    def rel_strain(self, vec1, vec2):
        """
            Calculate relative strain between two vectors
        """
        return np.linalg.norm(vec2)/np.linalg.norm(vec1)-1

    def rel_angle(self, vec_set1, vec_set2):
        """
            Calculate the relative angle between two vectors
        """
        return self.vec_angle(vec_set2[0], vec_set2[1])/self.vec_angle(
            vec_set1[0], vec_set1[1])-1

    def is_same_vectors(self, vec_set1, vec_set2):
        """
            Determine if two sets of vectors are the same within length and
            angle tolerances specified in the criteria dictionary
        """
        if (np.absolute(self.rel_strain(vec_set1[0], vec_set2[0])) >
                self.criteria["max_length_tol"]):
            return False
        elif (np.absolute(self.rel_strain(vec_set1[1], vec_set2[1])) >
            self.criteria["max_length_tol"]):
            return False
        elif (np.absolute(self.rel_angle(vec_set1, vec_set2)) >
            self.criteria["max_angle_tol"]):
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

        for i in self.factor(area_multiple):
            for j in range(area_multiple//i):
                yield np.matrix(((i, j), (0, area_multiple/i)))

    def generate_sl_transformations(self, film_area, substrate_area):
        """
            Generates transformation sets for film/substrate pair given the
            area of the unit cell area for the film and substrate. The
            transformation sets map the film and substrate unit cells to super
            lattices with a maximum area given in criteria["max_area"]

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

        for i in range(1, int(self.criteria["max_area"]/film_area)):
            for j in range (1, int(self.criteria["max_area"]/substrate_area) ):
                if (gcd(i, j) == 1 and
                        np.absolute(film_area/substrate_area - float(j)/i) <
                        self.criteria["max_area_ratio_tol"]):

                    yield [(i, j), self.generate_sl_transformation(i),
                        self.generate_sl_transformation(j)]

    def check_transformations(self, transformation_sets, film_vectors,
        substrate_vectors):
        """
            Applies the transformation_sets to the film and substrate vectors
            to generate super-lattices and checks if they matches.
            Returns all matching vectors sets.
        """

        for t in transformation_sets:
            films = []
            substrates = []
            # Apply transformations and reduce using Zur reduce methodology
            for f in t[1]:
                films.append(self.reduce_vectors(*np.squeeze(np.asarray(
                    f*film_vectors))))
            for s in t[2]:
                substrates.append(self.reduce_vectors(*np.squeeze(np.asarray(
                    s*substrate_vectors))))
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
            film_vectors = self.reduce_vectors(film_slab.lattice_vectors()[0],
                film_slab.lattice_vectors()[1])
            film_area = self.area(*film_vectors)

            for s in substrate_millers:
                substrate_slab = SlabGenerator(self.substrate, s, 20, 15,
                    primitive=False).get_slab()
                substrate_vectors = self.reduce_vectors(
                    substrate_slab.lattice_vectors()[0],
                    substrate_slab.lattice_vectors()[1])
                substrate_area = self.area(*substrate_vectors)

                yield [film_area, substrate_area, film_vectors,
                    substrate_vectors, f]

    def generate(self, film_millers = None, substrate_millers = None):
        """
            Generates the film/substrate combinations for either set miller
            indicies or all possible miller indices up to a max miller index
        """

        # Generate miller indicies if none specified for film
        if film_millers is None:
            film_millers = sorted(get_symmetrically_distinct_miller_indices(
                self.film, self.criteria["film_max_miller_index"]))

        # Generate miller indicies if none specified for substrate
        if substrate_millers is None:
            substrate_millers = sorted(
                get_symmetrically_distinct_miller_indices(self.substrate,
                    self.criteria["sub_max_miller_index"]))


        # Check each miller index combination
        for slab_combination in self.generate_slabs(film_millers,
            substrate_millers):
            # Generate all super lattice comnbinations for a given set of miller
            # indicies
            transformations = self.generate_sl_transformations(
                slab_combination[0], slab_combination[1])
            # Check each super-lattice pair to see if they match
            for match in self.check_transformations(transformations,
                slab_combination[2], slab_combination[3]):
                yield [slab_combination[4], self.area(*match[0]), match]


    def find_min(self, matches):
        """
            Find the minimum area match for each miller index in list
        """
        area_result = dict()

        for m in matches:
            if m[0] in area_result.keys():
                if m[1] < area_result[m[0]][1]:
                    area_result[m[0]] = m
            else:
                area_result[m[0]] = m

        return area_result
