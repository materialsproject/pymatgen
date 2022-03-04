# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines site transformations which transforms a structure into
another structure. Site transformations differ from standard transformations
in that they operate in a site-specific manner.
All transformations should inherit the AbstractTransformation ABC.
"""

import itertools
import logging
import math
import time

import numpy as np
from monty.json import MSONable

from pymatgen.analysis.ewald import EwaldMinimizer, EwaldSummation
from pymatgen.analysis.local_env import MinimumDistanceNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.transformation_abc import AbstractTransformation


class InsertSitesTransformation(AbstractTransformation):
    """
    This transformation substitutes certain sites with certain species.
    """

    def __init__(self, species, coords, coords_are_cartesian=False, validate_proximity=True):
        """
        Args:
            species: A list of species. e.g., ["Li", "Fe"]
            coords: A list of coords corresponding to those species. e.g.,
                [[0,0,0],[0.5,0.5,0.5]].
            coords_are_cartesian (bool): Set to True if coords are given in
                cartesian coords. Defaults to False.
            validate_proximity (bool): Set to False if you do not wish to ensure
                that added sites are not too close to other sites. Defaults to True.
        """
        if len(species) != len(coords):
            raise ValueError("Species and coords must be the same length!")
        self.species = species
        self.coords = coords
        self.coords_are_cartesian = coords_are_cartesian
        self.validate_proximity = validate_proximity

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Arg:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.

        Return:
            Returns a copy of structure with sites inserted.
        """
        s = structure.copy()
        for i, sp in enumerate(self.species):
            s.insert(
                i,
                sp,
                self.coords[i],
                coords_are_cartesian=self.coords_are_cartesian,
                validate_proximity=self.validate_proximity,
            )
        return s.get_sorted_structure()

    def __str__(self):
        return "InsertSiteTransformation : " + f"species {self.species}, coords {self.coords}"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """Return: None"""
        return None

    @property
    def is_one_to_many(self):
        """Return: False"""
        return False


class ReplaceSiteSpeciesTransformation(AbstractTransformation):
    """
    This transformation substitutes certain sites with certain species.
    """

    def __init__(self, indices_species_map):
        """
        Args:
            indices_species_map: A dict containing the species mapping in
                int-string pairs. E.g., { 1:"Na"} or {2:"Mn2+"}. Multiple
                substitutions can be done. Overloaded to accept sp_and_occu
                dictionary. E.g. {1: {"Ge":0.75, "C":0.25} }, which
                substitutes a single species with multiple species to generate a
                disordered structure.
        """
        self.indices_species_map = indices_species_map

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Arg:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.

        Return:
            Returns a copy of structure with sites replaced.
        """
        s = structure.copy()
        for i, sp in self.indices_species_map.items():
            s[int(i)] = sp
        return s

    def __str__(self):
        return "ReplaceSiteSpeciesTransformation :" + ", ".join(
            [f"{k}->{v}" + v for k, v in self.indices_species_map.items()]
        )

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """Return: None"""
        return None

    @property
    def is_one_to_many(self):
        """Return: False"""
        return False


class RemoveSitesTransformation(AbstractTransformation):
    """
    Remove certain sites in a structure.
    """

    def __init__(self, indices_to_remove):
        """
        Args:
            indices_to_remove: List of indices to remove. E.g., [0, 1, 2]
        """

        self.indices_to_remove = indices_to_remove

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Arg:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.

        Return:
            Returns a copy of structure with sites removed.
        """
        s = structure.copy()
        s.remove_sites(self.indices_to_remove)
        return s

    def __str__(self):
        return "RemoveSitesTransformation :" + ", ".join(map(str, self.indices_to_remove))

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """Return: None"""
        return None

    @property
    def is_one_to_many(self):
        """Return: False"""
        return False


class TranslateSitesTransformation(AbstractTransformation):
    """
    This class translates a set of sites by a certain vector.
    """

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

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Arg:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.

        Return:
            Returns a copy of structure with sites translated.
        """
        s = structure.copy()
        if self.translation_vector.shape == (len(self.indices_to_move), 3):
            for i, idx in enumerate(self.indices_to_move):
                s.translate_sites(idx, self.translation_vector[i], self.vector_in_frac_coords)
        else:
            s.translate_sites(self.indices_to_move, self.translation_vector, self.vector_in_frac_coords)
        return s

    def __str__(self):
        return (
            f"TranslateSitesTransformation for indices {self.indices_to_move}, "
            f"vect {self.translation_vector} and "
            f"vect_in_frac_coords = {self.vector_in_frac_coords}"
        )

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """
        Returns:
            TranslateSitesTransformation with the reverse translation.
        """
        return TranslateSitesTransformation(self.indices_to_move, -self.translation_vector, self.vector_in_frac_coords)

    @property
    def is_one_to_many(self):
        """Return: False"""
        return False

    def as_dict(self):
        """
        Json-serializable dict representation.
        """
        d = MSONable.as_dict(self)
        d["translation_vector"] = self.translation_vector.tolist()
        return d


class PartialRemoveSitesTransformation(AbstractTransformation):
    """
    Remove fraction of specie from a structure.
    Requires an oxidation state decorated structure for ewald sum to be
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
        ALGO_FAST.  For example, if you have 48 sites of which you want to
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
                A list of list of indices.
                e.g. [[0, 1], [2, 3, 4, 5]]
            fractions:
                The corresponding fractions to remove. Must be same length as
                indices. e.g., [0.5, 0.25]
            algo:
                This parameter allows you to choose the algorithm to perform
                ordering. Use one of PartialRemoveSpecieTransformation.ALGO_*
                variables to set the algo.
        """
        self.indices = indices
        self.fractions = fractions
        self.algo = algo
        self.logger = logging.getLogger(self.__class__.__name__)

    def _best_first_ordering(self, structure, num_remove_dict):
        self.logger.debug("Performing best first ordering")
        starttime = time.time()
        self.logger.debug("Performing initial ewald sum...")
        ewaldsum = EwaldSummation(structure)
        self.logger.debug(f"Ewald sum took {time.time() - starttime} seconds.")
        starttime = time.time()

        ematrix = ewaldsum.total_energy_matrix
        to_delete = []

        totalremovals = sum(num_remove_dict.values())
        removed = {k: 0 for k in num_remove_dict.keys()}
        for i in range(totalremovals):
            maxindex = None
            maxe = float("-inf")
            maxindices = None
            for indices in num_remove_dict.keys():
                if removed[indices] < num_remove_dict[indices]:
                    for ind in indices:
                        if ind not in to_delete:
                            energy = sum(ematrix[:, ind]) + sum(ematrix[:, ind]) - ematrix[ind, ind]
                            if energy > maxe:
                                maxindex = ind
                                maxe = energy
                                maxindices = indices
            removed[maxindices] += 1
            to_delete.append(maxindex)
            ematrix[:, maxindex] = 0
            ematrix[maxindex, :] = 0
        s = structure.copy()
        s.remove_sites(to_delete)
        self.logger.debug(f"Minimizing Ewald took {time.time() - starttime} seconds.")
        return [{"energy": sum(sum(ematrix)), "structure": s.get_sorted_structure()}]

    def _complete_ordering(self, structure, num_remove_dict):
        self.logger.debug("Performing complete ordering...")
        all_structures = []
        symprec = 0.2
        s = SpacegroupAnalyzer(structure, symprec=symprec)
        self.logger.debug(f"Symmetry of structure is determined to be {s.get_space_group_symbol()}.")
        sg = s.get_space_group_operations()
        tested_sites = []
        starttime = time.time()
        self.logger.debug("Performing initial ewald sum...")
        ewaldsum = EwaldSummation(structure)
        self.logger.debug(f"Ewald sum took {time.time() - starttime} seconds.")
        starttime = time.time()

        allcombis = []
        for ind, num in num_remove_dict.items():
            allcombis.append(itertools.combinations(ind, num))

        count = 0
        for allindices in itertools.product(*allcombis):
            sites_to_remove = []
            indices_list = []
            for indices in allindices:
                sites_to_remove.extend([structure[i] for i in indices])
                indices_list.extend(indices)
            s_new = structure.copy()
            s_new.remove_sites(indices_list)
            energy = ewaldsum.compute_partial_energy(indices_list)
            already_tested = False
            for i, tsites in enumerate(tested_sites):
                tenergy = all_structures[i]["energy"]
                if abs((energy - tenergy) / len(s_new)) < 1e-5 and sg.are_symmetrically_equivalent(
                    sites_to_remove, tsites, symm_prec=symprec
                ):
                    already_tested = True

            if not already_tested:
                tested_sites.append(sites_to_remove)
                all_structures.append({"structure": s_new, "energy": energy})

            count += 1
            if count % 10 == 0:
                timenow = time.time()
                self.logger.debug(f"{count} structures, {timenow - starttime:.2f} seconds.")
                self.logger.debug(f"Average time per combi = {(timenow - starttime) / count} seconds")
                self.logger.debug(f"{len(all_structures)} symmetrically distinct structures found.")

        self.logger.debug(f"Total symmetrically distinct structures found = {len(all_structures)}")
        all_structures = sorted(all_structures, key=lambda s: s["energy"])
        return all_structures

    def _fast_ordering(self, structure, num_remove_dict, num_to_return=1):
        """
        This method uses the matrix form of ewaldsum to calculate the ewald
        sums of the potential structures. This is on the order of 4 orders of
        magnitude faster when there are large numbers of permutations to
        consider. There are further optimizations possible (doing a smarter
        search of permutations for example), but this won't make a difference
        until the number of permutations is on the order of 30,000.
        """
        self.logger.debug("Performing fast ordering")
        starttime = time.time()
        self.logger.debug("Performing initial ewald sum...")

        ewaldmatrix = EwaldSummation(structure).total_energy_matrix
        self.logger.debug(f"Ewald sum took {time.time() - starttime} seconds.")
        starttime = time.time()
        m_list = []
        for indices, num in num_remove_dict.items():
            m_list.append([0, num, list(indices), None])

        self.logger.debug("Calling EwaldMinimizer...")
        minimizer = EwaldMinimizer(ewaldmatrix, m_list, num_to_return, PartialRemoveSitesTransformation.ALGO_FAST)
        self.logger.debug(f"Minimizing Ewald took {time.time() - starttime} seconds.")

        all_structures = []

        lowest_energy = minimizer.output_lists[0][0]
        num_atoms = sum(structure.composition.values())

        for output in minimizer.output_lists:
            s = structure.copy()
            del_indices = []

            for manipulation in output[1]:
                if manipulation[1] is None:
                    del_indices.append(manipulation[0])
                else:
                    s.replace(manipulation[0], manipulation[1])
            s.remove_sites(del_indices)
            struct = s.get_sorted_structure()
            all_structures.append(
                {
                    "energy": output[0],
                    "energy_above_minimum": (output[0] - lowest_energy) / num_atoms,
                    "structure": struct,
                }
            )

        return all_structures

    def _enumerate_ordering(self, structure):
        # Generate the disordered structure first.
        s = structure.copy()
        for indices, fraction in zip(self.indices, self.fractions):
            for ind in indices:
                new_sp = {sp: occu * fraction for sp, occu in structure[ind].species.items()}
                s[ind] = new_sp
        # Perform enumeration
        from pymatgen.transformations.advanced_transformations import (
            EnumerateStructureTransformation,
        )

        trans = EnumerateStructureTransformation()
        return trans.apply_transformation(s, 10000)

    def apply_transformation(self, structure, return_ranked_list=False):
        """
        Apply the transformation.

        Args:
            structure: input structure
            return_ranked_list (bool): Whether or not multiple structures are
                returned. If return_ranked_list is a number, that number of
                structures is returned.

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
        total_combis = 0
        for indices, frac in zip(self.indices, self.fractions):
            num_to_remove = len(indices) * frac
            if abs(num_to_remove - int(round(num_to_remove))) > 1e-3:
                raise ValueError("Fraction to remove must be consistent with integer amounts in structure.")
            num_to_remove = int(round(num_to_remove))
            num_remove_dict[tuple(indices)] = num_to_remove
            n = len(indices)
            total_combis += int(
                round(math.factorial(n) / math.factorial(num_to_remove) / math.factorial(n - num_to_remove))
            )

        self.logger.debug(f"Total combinations = {total_combis}")

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
        return opt_s if not return_ranked_list else all_structures[0:num_to_return]

    def __str__(self):
        return f"PartialRemoveSitesTransformation : Indices and fraction to remove = {self.indices}, ALGO = {self.algo}"

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """Return: None"""
        return None

    @property
    def is_one_to_many(self):
        """Return: True"""
        return True


class AddSitePropertyTransformation(AbstractTransformation):
    """
    Simple transformation to add site properties to a given structure
    """

    def __init__(self, site_properties):
        """
        Args:
            site_properties (dict): site properties to be added to a structure
        """
        self.site_properties = site_properties

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Arg:
            structure (Structure): A structurally similar structure in
                regards to crystal and site positions.

        Return:
            Returns a copy of structure with sites properties added.
        """
        new_structure = structure.copy()
        for prop in self.site_properties.keys():
            new_structure.add_site_property(prop, self.site_properties[prop])
        return new_structure

    @property
    def inverse(self):
        """Return: None"""
        return None

    @property
    def is_one_to_many(self):
        """Return: False"""
        return False


class RadialSiteDistortionTransformation(AbstractTransformation):
    """
    Radially perturbs atoms around a site. Can be used to create spherical distortion due to a
    point defect.
    """

    def __init__(self, site_index, displacement=0.1, nn_only=False):
        """
        Args:
            site_index (int): index of the site in structure to place at the center of the distortion (will
                not be distorted). This index must be provided before the structure is provided in
                apply_transformation in order to keep in line with the base class.
            displacement (float): distance to perturb the atoms around the objective site
            nn_only (bool): Whether or not to perturb beyond the nearest neighbors. If True, then only the
                nearest neighbors will be perturbed, leaving the other sites undisturbed. If False, then
                the nearest neighbors will receive the full displacement, and then subsequent sites will receive
                a displacement=0.1 / r, where r is the distance each site to the origin site. For small displacements,
                atoms beyond the NN environment will receive very small displacements, and these are almost equal.
                For large displacements, this difference is noticeable.
        """
        self.site_index = site_index
        self.displacement = displacement
        self.nn_only = nn_only

    def apply_transformation(self, structure):
        """
        Apply the transformation.

        Args:
            structure: Structure or Molecule to apply the transformation to

        Returns:
            the transformed structure
        """
        structure = structure.copy()
        site = structure[self.site_index]

        def f(x, r, r0):
            return x * r0 / r

        r0 = max(site.distance(_["site"]) for _ in MinimumDistanceNN().get_nn_info(structure, self.site_index))
        if hasattr(structure, "lattice"):
            m = structure.lattice.matrix
            m = (abs(m) > 1e-5) * m
            a, b, c = m[0], m[1], m[2]
            x = abs(np.dot(a, np.cross(b, c)) / np.linalg.norm(np.cross(b, c)))
            y = abs(np.dot(b, np.cross(a, c)) / np.linalg.norm(np.cross(a, c)))
            z = abs(np.dot(c, np.cross(a, b)) / np.linalg.norm(np.cross(a, b)))
            rmax = np.floor(min([x, y, z]) / 2)
        else:
            rmax = np.max(structure.distance_matrix)

        for vals in structure.get_neighbors(site, r=r0 if self.nn_only else rmax):
            site2, distance, index = vals[:3]
            v = site2.coords - site.coords
            kwargs = {"indices": [index], "vector": v * f(self.displacement, distance, r0) / np.linalg.norm(v)}
            if hasattr(structure, "lattice"):
                kwargs["frac_coords"] = False
            structure.translate_sites(**kwargs)
        return structure

    @property
    def inverse(self):
        """
        Returns the inverse transformation if available.
        Otherwise, should return None.
        """
        return False

    @property
    def is_one_to_many(self):
        """
        Determines if a Transformation is a one-to-many transformation. If a
        Transformation is a one-to-many transformation, the
        apply_transformation method should have a keyword arg
        "return_ranked_list" which allows for the transformed structures to be
        returned as a ranked list.
        """
        return False

    @property
    def use_multiprocessing(self):
        """
        Indicates whether the transformation can be applied by a
        subprocessing pool. This should be overridden to return True for
        transformations that the transmuter can parallelize.
        """
        return False
