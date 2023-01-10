# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements the Zur and McGill lattice matching algorithm
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Iterator

import numpy as np
from monty.json import MSONable


@dataclass
class ZSLMatch(MSONable):
    """
    A match from the Zur and McGill Algorithm. The super_lattice vectors are listed
    as _sl_vectors. These are reduced according to the algorithm in the paper which
    effectively a rotation in 3D space. Use the match_transformation property to get
    the appropriate transformation matrix
    """

    film_sl_vectors: list
    substrate_sl_vectors: list
    film_vectors: list
    substrate_vectors: list
    film_transformation: list
    substrate_transformation: list

    @property
    def match_area(self):
        """The area of the match between the substrate and film super lattice vectors"""
        return vec_area(*self.film_sl_vectors)

    @property
    def match_transformation(self):
        """The transformation matrix to convert the film super lattice vectors to the substrate"""
        # Generate 3D lattice vectors for film super lattice
        film_matrix = list(self.film_sl_vectors)
        film_matrix.append(np.cross(film_matrix[0], film_matrix[1]))

        # Generate 3D lattice vectors for substrate super lattice
        # Out of plane substrate super lattice has to be same length as
        # Film out of plane vector to ensure no extra deformation in that
        # direction
        substrate_matrix = list(self.substrate_sl_vectors)
        temp_sub = np.cross(substrate_matrix[0], substrate_matrix[1])
        temp_sub = temp_sub * fast_norm(film_matrix[2]) / fast_norm(temp_sub)
        substrate_matrix.append(temp_sub)

        transform_matrix = np.transpose(np.linalg.solve(film_matrix, substrate_matrix))

        return transform_matrix


class ZSLGenerator(MSONable):
    """
    This class generate matching interface super lattices based on the methodology
    of lattice vector matching for heterostructural interfaces proposed by
    Zur and McGill:
    Journal of Applied Physics 55 (1984), 378 ; doi: 10.1063/1.333084
    The process of generating all possible matching super lattices is:
    1.) Reduce the surface lattice vectors and calculate area for the surfaces
    2.) Generate all super lattice transformations within a maximum allowed area
        limit that give nearly equal area super-lattices for the two
        surfaces - generate_sl_transformation_sets
    3.) For each superlattice set:
        1.) Reduce super lattice vectors
        2.) Check length and angle between film and substrate super lattice
            vectors to determine if the super lattices are the nearly same
            and therefore coincident - get_equiv_transformations
    """

    def __init__(
        self,
        max_area_ratio_tol=0.09,
        max_area=400,
        max_length_tol=0.03,
        max_angle_tol=0.01,
        bidirectional=False,
    ):
        """
        Initialize a Zur Super Lattice Generator for a specific film and
            substrate
        Args:
            max_area_ratio_tol(float): Max tolerance on ratio of
                super-lattices to consider equal
            max_area(float): max super lattice area to generate in search
            max_length_tol: maximum length tolerance in checking if two
                vectors are of nearly the same length
            max_angle_tol: maximum angle tolerance in checking of two sets
                of vectors have nearly the same angle between them
        """
        self.max_area_ratio_tol = max_area_ratio_tol
        self.max_area = max_area
        self.max_length_tol = max_length_tol
        self.max_angle_tol = max_angle_tol
        self.bidirectional = bidirectional

    def is_same_vectors(self, vec_set1, vec_set2) -> bool:
        """
        Determine if two sets of vectors are the same within length and angle
        tolerances
        Args:
            vec_set1(array[array]): an array of two vectors
            vec_set2(array[array]): second array of two vectors
        """
        if self.bidirectional:
            return self._bidirectional_same_vectors(vec_set1, vec_set2)

        return self._unidirectional_is_same_vectors(vec_set1, vec_set2)

    def _unidirectional_is_same_vectors(self, vec_set1, vec_set2):
        """
        Determine if two sets of vectors are the same within length and angle
        tolerances
        Args:
            vec_set1(array[array]): an array of two vectors
            vec_set2(array[array]): second array of two vectors
        """
        if np.absolute(rel_strain(vec_set1[0], vec_set2[0])) > self.max_length_tol:
            return False
        if np.absolute(rel_strain(vec_set1[1], vec_set2[1])) > self.max_length_tol:
            return False
        if np.absolute(rel_angle(vec_set1, vec_set2)) > self.max_angle_tol:
            return False
        return True

    def _bidirectional_same_vectors(self, vec_set1, vec_set2):
        """Bidirectional version of above matching constraint check"""
        return self._unidirectional_is_same_vectors(vec_set1, vec_set2) or self._unidirectional_is_same_vectors(
            vec_set2, vec_set1
        )

    def generate_sl_transformation_sets(self, film_area, substrate_area):
        """
        Generates transformation sets for film/substrate pair given the
        area of the unit cell area for the film and substrate. The
        transformation sets map the film and substrate unit cells to super
        lattices with a maximum area
        Args:
            film_area(int): the unit cell area for the film
            substrate_area(int): the unit cell area for the substrate
        Returns:
            transformation_sets: a set of transformation_sets defined as:
                1.) the transformation matrices for the film to create a
                super lattice of area i*film area
                2.) the transformation matrices for the substrate to create
                a super lattice of area j*film area
        """
        transformation_indices = [
            (i, j)
            for i in range(1, int(self.max_area / film_area))
            for j in range(1, int(self.max_area / substrate_area))
            if np.absolute(film_area / substrate_area - float(j) / i) < self.max_area_ratio_tol
        ] + [
            (i, j)
            for i in range(1, int(self.max_area / film_area))
            for j in range(1, int(self.max_area / substrate_area))
            if np.absolute(substrate_area / film_area - float(i) / j) < self.max_area_ratio_tol
        ]
        transformation_indices = list(set(transformation_indices))

        # Sort sets by the square of the matching area and yield in order
        # from smallest to largest
        for i, j in sorted(transformation_indices, key=lambda x: x[0] * x[1]):
            yield (gen_sl_transform_matricies(i), gen_sl_transform_matricies(j))

    def get_equiv_transformations(self, transformation_sets, film_vectors, substrate_vectors):
        """
        Applies the transformation_sets to the film and substrate vectors
        to generate super-lattices and checks if they matches.
        Returns all matching vectors sets.
        Args:
            transformation_sets(array): an array of transformation sets:
                each transformation set is an array with the (i,j)
                indicating the area multiples of the film and substrate it
                corresponds to, an array with all possible transformations
                for the film area multiple i and another array for the
                substrate area multiple j.
            film_vectors(array): film vectors to generate super lattices
            substrate_vectors(array): substrate vectors to generate super
                lattices
        """
        for (film_transformations, substrate_transformations) in transformation_sets:
            # Apply transformations and reduce using Zur reduce methodology
            films = [reduce_vectors(*np.dot(f, film_vectors)) for f in film_transformations]

            substrates = [reduce_vectors(*np.dot(s, substrate_vectors)) for s in substrate_transformations]

            # Check if equivalent super lattices
            for (f_trans, s_trans), (f, s) in zip(
                product(film_transformations, substrate_transformations),
                product(films, substrates),
            ):
                if self.is_same_vectors(f, s):
                    yield [f, s, f_trans, s_trans]

    def __call__(self, film_vectors, substrate_vectors, lowest=False) -> Iterator[ZSLMatch]:
        """
        Runs the ZSL algorithm to generate all possible matching
        :return:
        """
        film_area = vec_area(*film_vectors)
        substrate_area = vec_area(*substrate_vectors)

        # Generate all super lattice combinations for a given set of miller
        # indices
        transformation_sets = self.generate_sl_transformation_sets(film_area, substrate_area)

        # Check each super-lattice pair to see if they match
        for match in self.get_equiv_transformations(transformation_sets, film_vectors, substrate_vectors):
            # Yield the match area, the miller indices,
            yield ZSLMatch(
                film_sl_vectors=match[0],
                substrate_sl_vectors=match[1],
                film_vectors=film_vectors,
                substrate_vectors=substrate_vectors,
                film_transformation=match[2],
                substrate_transformation=match[3],
            )

            # Just want lowest match per direction
            if lowest:
                break


def gen_sl_transform_matricies(area_multiple):
    """
    Generates the transformation matricies that convert a set of 2D
    vectors into a super lattice of integer area multiple as proven
    in Cassels:

    Cassels, John William Scott. An introduction to the geometry of
    numbers. Springer Science & Business Media, 2012.

    Args:
        area_multiple(int): integer multiple of unit cell area for super
        lattice area

    Returns:
        matrix_list: transformation matricies to convert unit vectors to
        super lattice vectors
    """
    return [
        np.array(((i, j), (0, area_multiple / i)))
        for i in get_factors(area_multiple)
        for j in range(area_multiple // i)
    ]


def rel_strain(vec1, vec2):
    """
    Calculate relative strain between two vectors
    """
    return fast_norm(vec2) / fast_norm(vec1) - 1


def rel_angle(vec_set1, vec_set2):
    """
    Calculate the relative angle between two vector sets

    Args:
        vec_set1(array[array]): an array of two vectors
        vec_set2(array[array]): second array of two vectors
    """
    return vec_angle(vec_set2[0], vec_set2[1]) / vec_angle(vec_set1[0], vec_set1[1]) - 1


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
