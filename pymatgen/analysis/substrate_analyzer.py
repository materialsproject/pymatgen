# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to identify optimal substrates for film growth
"""

from itertools import product

import numpy as np

from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.core.surface import (
    SlabGenerator,
    get_symmetrically_distinct_miller_indices,
)

__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Production"
__date__ = "Feb, 2016"


class ZSLGenerator:
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
    ):
        """
        Intialize a Zur Super Lattice Generator for a specific film and
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

    def is_same_vectors(self, vec_set1, vec_set2):
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
                1.) the transformation matricies for the film to create a
                super lattice of area i*film area
                2.) the tranformation matricies for the substrate to create
                a super lattice of area j*film area
        """
        transformation_indicies = [
            (i, j)
            for i in range(1, int(self.max_area / film_area))
            for j in range(1, int(self.max_area / substrate_area))
            if np.absolute(film_area / substrate_area - float(j) / i) < self.max_area_ratio_tol
        ]

        # Sort sets by the square of the matching area and yield in order
        # from smallest to largest
        for x in sorted(transformation_indicies, key=lambda x: x[0] * x[1]):
            yield (gen_sl_transform_matricies(x[0]), gen_sl_transform_matricies(x[1]))

    def get_equiv_transformations(self, transformation_sets, film_vectors, substrate_vectors):
        """
        Applies the transformation_sets to the film and substrate vectors
        to generate super-lattices and checks if they matches.
        Returns all matching vectors sets.
        Args:
            transformation_sets(array): an array of transformation sets:
                each transformation set is an array with the (i,j)
                indicating the area multipes of the film and subtrate it
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

            # Check if equivalant super lattices
            for (f_trans, s_trans), (f, s) in zip(
                product(film_transformations, substrate_transformations),
                product(films, substrates),
            ):
                if self.is_same_vectors(f, s):
                    yield [f, s, f_trans, s_trans]

    def __call__(self, film_vectors, substrate_vectors, lowest=False):
        """
        Runs the ZSL algorithm to generate all possible matching
        :return:
        """

        film_area = vec_area(*film_vectors)
        substrate_area = vec_area(*substrate_vectors)

        # Generate all super lattice comnbinations for a given set of miller
        # indicies
        transformation_sets = self.generate_sl_transformation_sets(film_area, substrate_area)

        # Check each super-lattice pair to see if they match
        for match in self.get_equiv_transformations(transformation_sets, film_vectors, substrate_vectors):
            # Yield the match area, the miller indicies,
            yield self.match_as_dict(
                match[0],
                match[1],
                film_vectors,
                substrate_vectors,
                vec_area(*match[0]),
                match[2],
                match[3],
            )

            # Just want lowest match per direction
            if lowest:
                break

    @staticmethod
    def match_as_dict(
        film_sl_vectors,
        substrate_sl_vectors,
        film_vectors,
        substrate_vectors,
        match_area,
        film_transformation,
        substrate_transformation,
    ):
        """
        Returns dict which contains ZSL match
        Args:
            film_miller(array)
            substrate_miller(array)
        """
        d = {}
        d["film_sl_vecs"] = np.asarray(film_sl_vectors)
        d["sub_sl_vecs"] = np.asarray(substrate_sl_vectors)
        d["match_area"] = match_area
        d["film_vecs"] = np.asarray(film_vectors)
        d["sub_vecs"] = np.asarray(substrate_vectors)
        d["film_transformation"] = np.asarray(film_transformation)
        d["substrate_transformation"] = np.asarray(substrate_transformation)

        return d


class SubstrateAnalyzer:
    """
    This class applies a set of search criteria to identify suitable
    substrates for film growth. It first uses a topoplogical search by Zur
    and McGill to identify matching super-lattices on various faces of the
    two materials. Additional criteria can then be used to identify the most
    suitable substrate. Currently, the only additional criteria is the
    elastic strain energy of the super-lattices
    """

    def __init__(self, zslgen=ZSLGenerator(), film_max_miller=1, substrate_max_miller=1):
        """
        Initializes the substrate analyzer
        Args:
            zslgen(ZSLGenerator): Defaults to a ZSLGenerator with standard
                tolerances, but can be fed one with custom tolerances
            film_max_miller(int): maximum miller index to generate for film
                surfaces
            substrate_max_miller(int): maximum miller index to generate for
                substrate surfaces
        """
        self.zsl = zslgen
        self.film_max_miller = film_max_miller
        self.substrate_max_miller = substrate_max_miller

    def generate_surface_vectors(self, film_millers, substrate_millers):
        """
        Generates the film/substrate slab combinations for a set of given
        miller indicies

        Args:
            film_millers(array): all miller indices to generate slabs for
                film
            substrate_millers(array): all miller indicies to generate slabs
                for substrate
        """
        vector_sets = []

        for f in film_millers:
            film_slab = SlabGenerator(self.film, f, 20, 15, primitive=False).get_slab()
            film_vectors = reduce_vectors(film_slab.lattice.matrix[0], film_slab.lattice.matrix[1])

            for s in substrate_millers:
                substrate_slab = SlabGenerator(self.substrate, s, 20, 15, primitive=False).get_slab()
                substrate_vectors = reduce_vectors(substrate_slab.lattice.matrix[0], substrate_slab.lattice.matrix[1])

                vector_sets.append((film_vectors, substrate_vectors, f, s))

        return vector_sets

    def calculate(
        self,
        film,
        substrate,
        elasticity_tensor=None,
        film_millers=None,
        substrate_millers=None,
        ground_state_energy=0,
        lowest=False,
    ):
        """
        Finds all topological matches for the substrate and calculates elastic
        strain energy and total energy for the film if elasticity tensor and
        ground state energy are provided:

        Args:
            film(Structure): conventional standard structure for the film
            substrate(Structure): conventional standard structure for the
                substrate
            elasticity_tensor(ElasticTensor): elasticity tensor for the film
                in the IEEE orientation
            film_millers(array): film facets to consider in search as defined by
                miller indicies
            substrate_millers(array): substrate facets to consider in search as
                defined by miller indicies
            ground_state_energy(float): ground state energy for the film
            lowest(bool): only consider lowest matching area for each surface
        """
        self.film = film
        self.substrate = substrate

        # Generate miller indicies if none specified for film
        if film_millers is None:
            film_millers = sorted(get_symmetrically_distinct_miller_indices(self.film, self.film_max_miller))

        # Generate miller indicies if none specified for substrate
        if substrate_millers is None:
            substrate_millers = sorted(
                get_symmetrically_distinct_miller_indices(self.substrate, self.substrate_max_miller)
            )

        # Check each miller index combination
        surface_vector_sets = self.generate_surface_vectors(film_millers, substrate_millers)
        for [
            film_vectors,
            substrate_vectors,
            film_miller,
            substrate_miller,
        ] in surface_vector_sets:
            for match in self.zsl(film_vectors, substrate_vectors, lowest):
                match["film_miller"] = film_miller
                match["sub_miller"] = substrate_miller
                if elasticity_tensor is not None:
                    energy, strain = self.calculate_3D_elastic_energy(
                        film, match, elasticity_tensor, include_strain=True
                    )
                    match["elastic_energy"] = energy
                    match["strain"] = strain
                if ground_state_energy != 0:
                    match["total_energy"] = match.get("elastic_energy", 0) + ground_state_energy

                yield match

    def calculate_3D_elastic_energy(self, film, match, elasticity_tensor=None, include_strain=False):
        """
        Calculates the multi-plane elastic energy. Returns 999 if no elastic
        tensor was given on init

        Args:
            film(Structure): conventional standard structure for the film
            match(dictionary) : match dictionary from substrate analyzer
            elasticity_tensor(ElasticTensor): elasticity tensor for the film
            include_strain(bool): include strain in the output or not; changes
             return from just the energy to a tuple with the energy and strain
             in voigt notation
        """
        if elasticity_tensor is None:
            return 9999

        # Get the appropriate surface structure
        struc = SlabGenerator(self.film, match["film_miller"], 20, 15, primitive=False).get_slab().oriented_unit_cell

        # Generate 3D lattice vectors for film super lattice
        film_matrix = list(match["film_sl_vecs"])
        film_matrix.append(np.cross(film_matrix[0], film_matrix[1]))

        # Generate 3D lattice vectors for substrate super lattice
        # Out of plane substrate super lattice has to be same length as
        # Film out of plane vector to ensure no extra deformation in that
        # direction
        substrate_matrix = list(match["sub_sl_vecs"])
        temp_sub = np.cross(substrate_matrix[0], substrate_matrix[1])
        temp_sub = temp_sub * fast_norm(film_matrix[2]) / fast_norm(temp_sub)
        substrate_matrix.append(temp_sub)

        transform_matrix = np.transpose(np.linalg.solve(film_matrix, substrate_matrix))

        dfm = Deformation(transform_matrix)

        strain = dfm.green_lagrange_strain.convert_to_ieee(struc, initial_fit=False)

        energy_density = elasticity_tensor.energy_density(strain)

        if include_strain:
            return (
                film.volume * energy_density / len(film.sites),
                strain.von_mises_strain,
            )
        return film.volume * energy_density / len(film.sites)


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
        matrix_list: transformation matricies to covert unit vectors to
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
