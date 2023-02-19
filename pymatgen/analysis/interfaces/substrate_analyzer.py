# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to identify optimal substrates for film growth
"""

from __future__ import annotations

from dataclasses import dataclass

from pymatgen.analysis.elasticity.strain import Deformation, Strain
from pymatgen.analysis.interfaces.zsl import ZSLGenerator, ZSLMatch, reduce_vectors
from pymatgen.core import Structure
from pymatgen.core.surface import (
    SlabGenerator,
    get_symmetrically_distinct_miller_indices,
)


@dataclass
class SubstrateMatch(ZSLMatch):
    """
    A substrate match building on the Zur and McGill algorithm. This match class includes the miller
    planes of the film and substrate the full strain tensor, the Von Mises strain, the ground state
    energy if provided, and the elastic energy
    """

    film_miller: tuple[int, int, int]
    substrate_miller: tuple[int, int, int]
    strain: Strain
    von_mises_strain: float
    ground_state_energy: float
    elastic_energy: float

    @classmethod
    def from_zsl(
        cls,
        match: ZSLMatch,
        film: Structure,
        film_miller,
        substrate_miller,
        elasticity_tensor=None,
        ground_state_energy=0,
    ):
        """Generate a substrate match from a ZSL match plus metadata"""
        # Get the appropriate surface structure
        struct = SlabGenerator(film, film_miller, 20, 15, primitive=False).get_slab().oriented_unit_cell

        dfm = Deformation(match.match_transformation)

        strain = dfm.green_lagrange_strain.convert_to_ieee(struct, initial_fit=False)
        von_mises_strain = strain.von_mises_strain

        if elasticity_tensor is not None:
            energy_density = elasticity_tensor.energy_density(strain)

            elastic_energy = film.volume * energy_density / len(film.sites)
        else:
            elastic_energy = 0

        return cls(
            film_miller=film_miller,
            substrate_miller=substrate_miller,
            strain=strain,
            von_mises_strain=von_mises_strain,
            elastic_energy=elastic_energy,
            ground_state_energy=ground_state_energy,
            **{
                k: getattr(match, k)
                for k in [
                    "film_sl_vectors",
                    "substrate_sl_vectors",
                    "film_vectors",
                    "substrate_vectors",
                    "film_transformation",
                    "substrate_transformation",
                ]
            },
        )

    @property
    def total_energy(self):
        """Total energy of this match"""
        return self.ground_state_energy + self.elastic_energy


class SubstrateAnalyzer(ZSLGenerator):
    """
    This class applies a set of search criteria to identify suitable
    substrates for film growth. It first uses a topoplogical search by Zur
    and McGill to identify matching super-lattices on various faces of the
    two materials. Additional criteria can then be used to identify the most
    suitable substrate. Currently, the only additional criteria is the
    elastic strain energy of the super-lattices
    """

    def __init__(self, film_max_miller=1, substrate_max_miller=1, **kwargs):
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
        self.film_max_miller = film_max_miller
        self.substrate_max_miller = substrate_max_miller
        self.kwargs = kwargs
        super().__init__(**kwargs)

    def generate_surface_vectors(self, film_millers, substrate_millers):
        """
        Generates the film/substrate slab combinations for a set of given
        miller indices

        Args:
            film_millers(array): all miller indices to generate slabs for
                film
            substrate_millers(array): all miller indices to generate slabs
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
                miller indices
            substrate_millers(array): substrate facets to consider in search as
                defined by miller indices
            ground_state_energy(float): ground state energy for the film
            lowest(bool): only consider lowest matching area for each surface
        """
        self.film = film
        self.substrate = substrate

        # Generate miller indices if none specified for film
        if film_millers is None:
            film_millers = sorted(get_symmetrically_distinct_miller_indices(self.film, self.film_max_miller))

        # Generate miller indices if none specified for substrate
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
            for match in self(film_vectors, substrate_vectors, lowest):
                sub_match = SubstrateMatch.from_zsl(
                    match=match,
                    film=film,
                    film_miller=film_miller,
                    substrate_miller=substrate_miller,
                    elasticity_tensor=elasticity_tensor,
                    ground_state_energy=ground_state_energy,
                )

                yield sub_match
