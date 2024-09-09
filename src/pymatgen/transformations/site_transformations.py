"""This module defines site transformations which transforms a structure into
another structure. Site transformations differ from standard transformations
in that they operate in a site-specific manner.
All transformations should inherit the AbstractTransformation ABC.
"""

from __future__ import annotations

import itertools
import logging
import math
import time
from typing import TYPE_CHECKING

import numpy as np
from monty.json import MSONable

from pymatgen.analysis.ewald import EwaldMinimizer, EwaldSummation
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.transformation_abc import AbstractTransformation

if TYPE_CHECKING:
    from pymatgen.core.sites import PeriodicSite
    from pymatgen.core.structure import Structure


class InsertSitesTransformation(AbstractTransformation):
    """This transformation substitutes certain sites with certain species."""

    def __init__(self, species, coords, coords_are_cartesian=False, validate_proximity=True):
        """
        Args:
            species: A list of species. e.g. ["Li", "Fe"]
            coords: A list of coords corresponding to those species. e.g.
                [[0,0,0],[0.5,0.5,0.5]].
            coords_are_cartesian (bool): Set to True if coords are given in
                Cartesian coords. Defaults to False.
            validate_proximity (bool): Set to False if you do not wish to ensure
                that added sites are not too close to other sites. Defaults to True.
        """
        if len(species) != len(coords):
            raise ValueError("Species and coords must be the same length!")
        self.species = species
        self.coords = coords
        self.coords_are_cartesian = coords_are_cartesian
        self.validate_proximity = validate_proximity

    def apply_transformation(self, structure: Structure):
        """Apply the transformation.

        Args:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.

        Returns:
            A copy of structure with sites inserted.
        """
        struct = structure.copy()
        for idx, sp in enumerate(self.species):
            struct.insert(
                idx,
                sp,
                self.coords[idx],
                coords_are_cartesian=self.coords_are_cartesian,
                validate_proximity=self.validate_proximity,
            )
        return struct.get_sorted_structure()

    def __repr__(self):
        return f"InsertSiteTransformation : species {self.species}, coords {self.coords}"


class ReplaceSiteSpeciesTransformation(AbstractTransformation):
    """This transformation substitutes certain sites with certain species."""

    def __init__(self, indices_species_map):
        """
        Args:
            indices_species_map: A dict containing the species mapping in
                int-string pairs. e.g. { 1:"Na"} or {2:"Mn2+"}. Multiple
                substitutions can be done. Overloaded to accept sp_and_occu
                dictionary. E.g. {1: {"Ge":0.75, "C":0.25} }, which
                substitutes a single species with multiple species to generate a
                disordered structure.
        """
        self.indices_species_map = indices_species_map

    def apply_transformation(self, structure: Structure):
        """Apply the transformation.

        Args:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.

        Returns:
            A copy of structure with sites replaced.
        """
        struct = structure.copy()
        for idx, sp in self.indices_species_map.items():
            struct[int(idx)] = sp
        return struct

    def __repr__(self):
        return "ReplaceSiteSpeciesTransformation :" + ", ".join(
            [f"{key}->{val}" + val for key, val in self.indices_species_map.items()]
        )


class RemoveSitesTransformation(AbstractTransformation):
    """Remove certain sites in a structure."""

    def __init__(self, indices_to_remove):
        """
        Args:
            indices_to_remove: List of indices to remove. e.g. [0, 1, 2].
        """
        self.indices_to_remove = indices_to_remove

    def apply_transformation(self, structure: Structure):
        """Apply the transformation.

        Args:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.

        Returns:
            A copy of structure with sites removed.
        """
        struct = structure.copy()
        struct.remove_sites(self.indices_to_remove)
        return struct

    def __repr__(self):
        return "RemoveSitesTransformation :" + ", ".join(map(str, self.indices_to_remove))


class TranslateSitesTransformation(AbstractTransformation):
    """This class translates a set of sites by a certain vector."""

    def __init__(self, indices_to_move, translation_vector, vector_in_frac_coords=True):
        """
        Args:
            indices_to_move: The indices of the sites to move
            translation_vector: Vector to move the sites. If a list of list or numpy
                array of shape, (len(indices_to_move), 3), is provided then each
                translation vector is applied to the corresponding site in the
                indices_to_move.
            vector_in_frac_coords: Set to True if the translation vector is in
                fractional coordinates, and False if it is in cartesian
                coordinations. Defaults to True.
        """
        self.indices_to_move = indices_to_move
        self.translation_vector = np.array(translation_vector)
        self.vector_in_frac_coords = vector_in_frac_coords

    def apply_transformation(self, structure: Structure):
        """Apply the transformation.

        Args:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.

        Returns:
            A copy of structure with sites translated.
        """
        struct = structure.copy()
        if self.translation_vector.shape == (len(self.indices_to_move), 3):
            for idx, idx in enumerate(self.indices_to_move):
                struct.translate_sites(idx, self.translation_vector[idx], self.vector_in_frac_coords)
        else:
            struct.translate_sites(self.indices_to_move, self.translation_vector, self.vector_in_frac_coords)
        return struct

    def __repr__(self):
        return (
            f"TranslateSitesTransformation for indices {self.indices_to_move}, "
            f"vect {self.translation_vector} and "
            f"vect_in_frac_coords = {self.vector_in_frac_coords}"
        )

    @property
    def inverse(self) -> TranslateSitesTransformation:
        """TranslateSitesTransformation with the reverse translation."""
        return TranslateSitesTransformation(self.indices_to_move, -self.translation_vector, self.vector_in_frac_coords)

    def as_dict(self):
        """JSON-serializable dict representation."""
        dct = MSONable.as_dict(self)
        dct["translation_vector"] = self.translation_vector.tolist()
        return dct


class PartialRemoveSitesTransformation(AbstractTransformation):
    """Remove fraction of specie from a structure.
    Requires an oxidation state decorated structure for Ewald sum to be
    computed.

    Given that the solution to selecting the right removals is NP-hard, there
    are several algorithms provided with varying degrees of accuracy and speed.
    The options are as follows:

    ALGO_FAST:
        This is a highly optimized algorithm to quickly go through the search
        tree. It is guaranteed to find the optimal solution, but will return
        only a single lowest energy structure. Typically, you will want to use
        this.

    ALGO_COMPLETE:
        The complete algo ensures that you get all symmetrically distinct
        orderings, ranked by the estimated Ewald energy. But this can be an
        extremely time-consuming process if the number of possible orderings is
        very large. Use this if you really want all possible orderings. If you
        want just the lowest energy ordering, ALGO_FAST is accurate and faster.

    ALGO_BEST_FIRST:
        This algorithm is for ordering the really large cells that defeats even
        ALGO_FAST. For example, if you have 48 sites of which you want to
        remove 16 of them, the number of possible orderings is around
        2 x 10^12. ALGO_BEST_FIRST shortcircuits the entire search tree by
        removing the highest energy site first, then followed by the next
        highest energy site, and so on. It is guaranteed to find a solution
        in a reasonable time, but it is also likely to be highly inaccurate.

    ALGO_ENUMERATE:
        This algorithm uses the EnumerateStructureTransformation to perform
        ordering. This algo returns *complete* orderings up to a single unit
        cell size. It is more robust than the ALGO_COMPLETE, but requires
        Gus Hart's enumlib to be installed.
    """

    ALGO_FAST = 0
    ALGO_COMPLETE = 1
    ALGO_BEST_FIRST = 2
    ALGO_ENUMERATE = 3

    def __init__(self, indices, fractions, algo=ALGO_COMPLETE):
        """
        Args:
            indices:
                A list of list of indices, e.g. [[0, 1], [2, 3, 4, 5]].
            fractions:
                The corresponding fractions to remove. Must be same length as
                indices. e.g. [0.5, 0.25]
            algo:
                This parameter allows you to choose the algorithm to perform
                ordering. Use one of PartialRemoveSpecieTransformation.ALGO_*
                variables to set the algo.
        """
        self.indices = indices
        self.fractions = fractions
        self.algo = algo
        self.logger = logging.getLogger(type(self).__name__)

    def _best_first_ordering(self, structure: Structure, num_remove_dict):
        self.logger.debug("Performing best first ordering")
        start_time = time.perf_counter()
        self.logger.debug("Performing initial Ewald sum...")
        ewald_sum = EwaldSummation(structure)
        self.logger.debug(f"Ewald sum took {time.perf_counter() - start_time} seconds.")
        start_time = time.perf_counter()

        e_matrix = ewald_sum.total_energy_matrix
        to_delete = []

        total_removals = sum(num_remove_dict.values())
        removed = dict.fromkeys(num_remove_dict, 0)
        for _ in range(total_removals):
            max_idx = None
            max_ene = float("-inf")
            max_indices = None
            for indices in num_remove_dict:
                if removed[indices] < num_remove_dict[indices]:
                    for ind in indices:
                        if ind not in to_delete:
                            energy = sum(e_matrix[:, ind]) + sum(e_matrix[:, ind]) - e_matrix[ind, ind]
                            if energy > max_ene:
                                max_idx = ind
                                max_ene = energy
                                max_indices = indices
            removed[max_indices] += 1
            to_delete.append(max_idx)
            e_matrix[:, max_idx] = 0
            e_matrix[max_idx, :] = 0
        struct = structure.copy()
        struct.remove_sites(to_delete)
        self.logger.debug(f"Minimizing Ewald took {time.perf_counter() - start_time} seconds.")
        return [{"energy": sum(e_matrix), "structure": struct.get_sorted_structure()}]

    def _complete_ordering(self, structure: Structure, num_remove_dict):
        self.logger.debug("Performing complete ordering...")
        all_structures: list[dict[str, float | Structure]] = []
        symprec = 0.2
        spg_analyzer = SpacegroupAnalyzer(structure, symprec=symprec)
        self.logger.debug(f"Symmetry of structure is determined to be {spg_analyzer.get_space_group_symbol()}.")
        sg = spg_analyzer.get_space_group_operations()
        tested_sites: list[list[PeriodicSite]] = []
        start_time = time.perf_counter()
        self.logger.debug("Performing initial Ewald sum...")
        ewald_sum = EwaldSummation(structure)
        self.logger.debug(f"Ewald sum took {time.perf_counter() - start_time} seconds.")
        start_time = time.perf_counter()

        all_combis = [list(itertools.combinations(ind, num)) for ind, num in num_remove_dict.items()]

        for idx, all_indices in enumerate(itertools.product(*all_combis), start=1):
            sites_to_remove = []
            indices_list = []
            for indices in all_indices:
                sites_to_remove.extend([structure[i] for i in indices])
                indices_list.extend(indices)
            s_new = structure.copy()
            s_new.remove_sites(indices_list)
            energy = ewald_sum.compute_partial_energy(indices_list)
            already_tested = False
            for ii, t_sites in enumerate(tested_sites):
                t_energy = all_structures[ii]["energy"]
                if abs((energy - t_energy) / len(s_new)) < 1e-5 and sg.are_symmetrically_equivalent(
                    set(sites_to_remove), set(t_sites), symm_prec=symprec
                ):
                    already_tested = True

            if not already_tested:
                tested_sites.append(sites_to_remove)
                all_structures.append({"structure": s_new, "energy": energy})

            if idx % 10 == 0:
                now = time.perf_counter()
                self.logger.debug(f"{idx} structures, {now - start_time:.2f} seconds.")
                self.logger.debug(f"Average time per combi = {(now - start_time) / idx} seconds")
                self.logger.debug(f"{len(all_structures)} symmetrically distinct structures found.")

        self.logger.debug(f"Total symmetrically distinct structures found = {len(all_structures)}")
        return sorted(all_structures, key=lambda s: s["energy"])

    def _fast_ordering(self, structure: Structure, num_remove_dict, num_to_return=1):
        """Use the matrix form of Ewald sum to calculate the Ewald
        sums of the potential structures. This is on the order of 4 orders of
        magnitude faster when there are large numbers of permutations to
        consider. There are further optimizations possible (doing a smarter
        search of permutations for example), but this won't make a difference
        until the number of permutations is on the order of 30,000.
        """
        self.logger.debug("Performing fast ordering")
        start_time = time.perf_counter()
        self.logger.debug("Performing initial Ewald sum...")

        ewald_matrix = EwaldSummation(structure).total_energy_matrix
        self.logger.debug(f"Ewald sum took {time.perf_counter() - start_time} seconds.")
        start_time = time.perf_counter()
        m_list = [[0, num, list(indices), None] for indices, num in num_remove_dict.items()]

        self.logger.debug("Calling EwaldMinimizer...")
        minimizer = EwaldMinimizer(ewald_matrix, m_list, num_to_return, PartialRemoveSitesTransformation.ALGO_FAST)
        self.logger.debug(f"Minimizing Ewald took {time.perf_counter() - start_time} seconds.")

        all_structures = []

        lowest_energy = minimizer.output_lists[0][0]
        num_atoms = sum(structure.composition.values())

        for output in minimizer.output_lists:
            struct = structure.copy()
            del_indices = []

            for manipulation in output[1]:
                if manipulation[1] is None:
                    del_indices.append(manipulation[0])
                else:
                    struct.replace(manipulation[0], manipulation[1])
            struct.remove_sites(del_indices)
            struct = struct.get_sorted_structure()
            e_above_min = (output[0] - lowest_energy) / num_atoms
            all_structures.append({"energy": output[0], "energy_above_minimum": e_above_min, "structure": struct})

        return all_structures

    def _enumerate_ordering(self, structure: Structure):
        # Generate the disordered structure first.
        struct = structure.copy()
        for indices, fraction in zip(self.indices, self.fractions, strict=True):
            for ind in indices:
                new_sp = {sp: occu * fraction for sp, occu in structure[ind].species.items()}
                struct[ind] = new_sp
        # Perform enumeration
        from pymatgen.transformations.advanced_transformations import EnumerateStructureTransformation

        trans = EnumerateStructureTransformation()
        return trans.apply_transformation(struct, 10000)

    def apply_transformation(self, structure: Structure, return_ranked_list: bool | int = False):
        """Apply the transformation.

        Args:
            structure: input structure
            return_ranked_list (bool | int): Whether or not multiple structures are returned.
                If return_ranked_list is int, that number of structures is returned.

        Returns:
            Depending on returned_ranked list, either a transformed structure
            or a list of dictionaries, where each dictionary is of the form
            {"structure" = .... , "other_arguments"}
            the key "transformation" is reserved for the transformation that
            was actually applied to the structure.
            This transformation is parsed by the alchemy classes for generating
            a more specific transformation history. Any other information will
            be stored in the transformation_parameters dictionary in the
            transmuted structure class.
        """
        num_remove_dict = {}
        total_combos = 0
        for idx, frac in zip(self.indices, self.fractions, strict=True):
            n_to_remove = len(idx) * frac
            if abs(n_to_remove - int(round(n_to_remove))) > 1e-3:
                raise ValueError("Fraction to remove must be consistent with integer amounts in structure.")
            n_to_remove = int(round(n_to_remove))
            num_remove_dict[tuple(idx)] = n_to_remove
            n = len(idx)
            total_combos += int(
                round(math.factorial(n) / math.factorial(n_to_remove) / math.factorial(n - n_to_remove))
            )

        self.logger.debug(f"Total combinations = {total_combos}")

        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1

        num_to_return = max(1, num_to_return)
        self.logger.debug(f"Will return {num_to_return} best structures.")

        if self.algo == PartialRemoveSitesTransformation.ALGO_FAST:
            all_structures = self._fast_ordering(structure, num_remove_dict, num_to_return)
        elif self.algo == PartialRemoveSitesTransformation.ALGO_COMPLETE:
            all_structures = self._complete_ordering(structure, num_remove_dict)
        elif self.algo == PartialRemoveSitesTransformation.ALGO_BEST_FIRST:
            all_structures = self._best_first_ordering(structure, num_remove_dict)
        elif self.algo == PartialRemoveSitesTransformation.ALGO_ENUMERATE:
            all_structures = self._enumerate_ordering(structure)
        else:
            raise ValueError("Invalid algo.")

        opt_s = all_structures[0]["structure"]
        return opt_s if not return_ranked_list else all_structures[:num_to_return]

    def __repr__(self):
        return f"PartialRemoveSitesTransformation : Indices and fraction to remove = {self.indices}, ALGO = {self.algo}"

    @property
    def is_one_to_many(self) -> bool:
        """Transform one structure to many."""
        return True


class AddSitePropertyTransformation(AbstractTransformation):
    """Simple transformation to add site properties to a given structure."""

    def __init__(self, site_properties):
        """
        Args:
            site_properties (dict): site properties to be added to a structure.
        """
        self.site_properties = site_properties

    def apply_transformation(self, structure: Structure):
        """Apply the transformation.

        Args:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.

        Returns:
            A copy of structure with sites properties added.
        """
        new_struct = structure.copy()
        for prop in self.site_properties:
            new_struct.add_site_property(prop, self.site_properties[prop])
        return new_struct


class RadialSiteDistortionTransformation(AbstractTransformation):
    """Radially perturbs atoms around a site. Can be used to create spherical distortion due to a
    point defect.
    """

    def __init__(self, site_index: int, displacement: float = 0.1, nn_only: bool = False) -> None:
        """
        Args:
            site_index (int): index of the site in structure to place at the center of the distortion (will
                not be distorted). This index must be provided before the structure is provided in
                apply_transformation in order to keep in line with the base class.
            displacement (float): distance to perturb the atoms around the objective site
            nn_only (bool): Whether to perturb beyond the nearest neighbors. If True, then only the
                nearest neighbors will be perturbed, leaving the other sites undisturbed. If False, then
                the nearest neighbors will receive the full displacement, and then subsequent sites will receive
                a displacement=0.1 / r, where r is the distance each site to the origin site. For small displacements,
                atoms beyond the NN environment will receive very small displacements, and these are almost equal.
                For large displacements, this difference is noticeable.
        """
        self.site_index = site_index
        self.displacement = displacement
        self.nn_only = nn_only

    def apply_transformation(self, structure: Structure):
        """Apply the transformation.

        Args:
            structure: Structure or Molecule to apply the transformation to

        Returns:
            the transformed structure
        """
        structure = structure.copy()
        site = structure[self.site_index]

        def displace_dist(x, r, r0):
            return x * r0 / r

        r0 = max(site.distance(_["site"]) for _ in MinimumDistanceNN().get_nn_info(structure, self.site_index))
        if hasattr(structure, "lattice"):
            latt_mat = structure.lattice.matrix
            latt_mat = (abs(latt_mat) > 1e-5) * latt_mat  # round small values to 0
            a, b, c = latt_mat[0], latt_mat[1], latt_mat[2]
            x = abs(np.dot(a, np.cross(b, c)) / np.linalg.norm(np.cross(b, c)))
            y = abs(np.dot(b, np.cross(a, c)) / np.linalg.norm(np.cross(a, c)))
            z = abs(np.dot(c, np.cross(a, b)) / np.linalg.norm(np.cross(a, b)))
            r_max = np.floor(min([x, y, z]) / 2)
        else:
            r_max = np.max(structure.distance_matrix)

        for vals in structure.get_neighbors(site, r=r0 if self.nn_only else r_max):
            site2, distance, index = vals[:3]
            vec = site2.coords - site.coords
            kwargs = {
                "indices": [index],
                "vector": vec * displace_dist(self.displacement, distance, r0) / np.linalg.norm(vec),
            }
            if hasattr(structure, "lattice"):
                kwargs["frac_coords"] = False
            structure.translate_sites(**kwargs)
        return structure

    @property
    def is_one_to_many(self) -> bool:
        """Determine if a Transformation is a one-to-many transformation. If a
        Transformation is a one-to-many transformation, the
        apply_transformation method should have a keyword arg
        "return_ranked_list" which allows for the transformed structures to be
        returned as a ranked list.
        """
        return False

    @property
    def use_multiprocessing(self):
        """Indicates whether the transformation can be applied by a
        subprocessing pool. This should be overridden to return True for
        transformations that the transmuter can parallelize.
        """
        return False
