# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

from six.moves import map
from six.moves import zip

import math
import itertools
import logging
import time

from monty.json import MSONable

import numpy as np

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.analysis.ewald import EwaldSummation, EwaldMinimizer

"""
This module defines site transformations which transforms a structure into
another structure. Site transformations differ from standard transformations
in that they operate in a site-specific manner.
All transformations should inherit the AbstractTransformation ABC.
"""

__author__ = "Shyue Ping Ong, Will Richards"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.2"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Sep 23, 2011"


class InsertSitesTransformation(AbstractTransformation):
    """
    This transformation substitutes certain sites with certain species.

    Args:
        species: A list of species. e.g., ["Li", "Fe"]
        coords: A list of coords corresponding to those species. e.g.,
            [[0,0,0],[0.5,0.5,0.5]].
        coords_are_cartesian (bool): Set to True if coords are given in
            cartesian coords. Defaults to False.
        validate_proximity (bool): Set to False if you do not wish to ensure that added sites are
            not too close to other sites. Defaults to True.
    """
    def __init__(self, species, coords, coords_are_cartesian=False,
                 validate_proximity=True):
        if len(species) != len(coords):
            raise ValueError("Species and coords must be the same length!")
        self.species = species
        self.coords = coords
        self.coords_are_cartesian = coords_are_cartesian
        self.validate_proximity = validate_proximity

    def apply_transformation(self, structure):
        s = structure.copy()
        for i, sp in enumerate(self.species):
            s.insert(i, sp, self.coords[i],
                     coords_are_cartesian=self.coords_are_cartesian,
                     validate_proximity=self.validate_proximity)
        return s.get_sorted_structure()

    def __str__(self):
        return "InsertSiteTransformation : " + \
            "species {}, coords {}".format(self.species, self.coords)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False


class ReplaceSiteSpeciesTransformation(AbstractTransformation):
    """
    This transformation substitutes certain sites with certain species.

    Args:
        indices_species_map: A dict containing the species mapping in
            int-string pairs. E.g., { 1:"Na"} or {2:"Mn2+"}. Multiple
            substitutions can be done. Overloaded to accept sp_and_occu
            dictionary. E.g. {"Si: {"Ge":0.75, "C":0.25} }, which
            substitutes a single species with multiple species to generate a disordered
            structure.
    """
    def __init__(self, indices_species_map):
        self.indices_species_map = indices_species_map

    def apply_transformation(self, structure):
        s = structure.copy()
        for i, sp in self.indices_species_map.items():
            s[int(i)] = sp
        return s

    def __str__(self):
        return "ReplaceSiteSpeciesTransformation :" + \
            ", ".join(["{}->{}".format(k, v) + v for k, v in
                       self.indices_species_map.items()])

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False


class RemoveSitesTransformation(AbstractTransformation):
    """
    Remove certain sites in a structure.

    Args:
        indices_to_remove: List of indices to remove. E.g., [0, 1, 2]
    """
    def __init__(self, indices_to_remove):
        self.indices_to_remove = indices_to_remove

    def apply_transformation(self, structure):
        s = structure.copy()
        s.remove_sites(self.indices_to_remove)
        return s

    def __str__(self):
        return "RemoveSitesTransformation :" + ", ".join(
            map(str, self.indices_to_remove))

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return False


class TranslateSitesTransformation(AbstractTransformation):
    """
    This class translates a set of sites by a certain vector.

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
    def __init__(self, indices_to_move, translation_vector,
                 vector_in_frac_coords=True):
        self.indices_to_move = indices_to_move
        self.translation_vector = np.array(translation_vector)
        self.vector_in_frac_coords = vector_in_frac_coords

    def apply_transformation(self, structure):
        s = structure.copy()
        if self.translation_vector.shape == (len(self.indices_to_move), 3):
            for i, idx in enumerate(self.indices_to_move):
                s.translate_sites(idx, self.translation_vector[i],
                                  self.vector_in_frac_coords)
        else:
            s.translate_sites(self.indices_to_move, self.translation_vector,
                              self.vector_in_frac_coords)
        return s

    def __str__(self):
        return "TranslateSitesTransformation for indices " + \
            "{}, vect {} and vect_in_frac_coords = {}".format(
                self.indices_to_move, self.translation_vector,
                self.vector_in_frac_coords)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return TranslateSitesTransformation(
            self.indices_to_move, -self.translation_vector,
            self.vector_in_frac_coords)

    @property
    def is_one_to_many(self):
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
        self.indices = indices
        self.fractions = fractions
        self.algo = algo
        self.logger = logging.getLogger(self.__class__.__name__)

    def best_first_ordering(self, structure, num_remove_dict):
        self.logger.debug("Performing best first ordering")
        starttime = time.time()
        self.logger.debug("Performing initial ewald sum...")
        ewaldsum = EwaldSummation(structure)
        self.logger.debug("Ewald sum took {} seconds."
                          .format(time.time() - starttime))
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
                            energy = sum(ematrix[:, ind]) + \
                                sum(ematrix[:, ind]) - ematrix[ind, ind]
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
        self.logger.debug("Minimizing Ewald took {} seconds."
                          .format(time.time() - starttime))
        return [{"energy": sum(sum(ematrix)),
                 "structure": s.get_sorted_structure()}]

    def complete_ordering(self, structure, num_remove_dict):
        self.logger.debug("Performing complete ordering...")
        all_structures = []
        symprec = 0.2
        s = SpacegroupAnalyzer(structure, symprec=symprec)
        self.logger.debug("Symmetry of structure is determined to be {}."
                          .format(s.get_space_group_symbol()))
        sg = s.get_space_group_operations()
        tested_sites = []
        starttime = time.time()
        self.logger.debug("Performing initial ewald sum...")
        ewaldsum = EwaldSummation(structure)
        self.logger.debug("Ewald sum took {} seconds."
                          .format(time.time() - starttime))
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
                if abs((energy - tenergy) / len(s_new)) < 1e-5 and \
                        sg.are_symmetrically_equivalent(sites_to_remove,
                                                        tsites,
                                                        symm_prec=symprec):
                    already_tested = True

            if not already_tested:
                tested_sites.append(sites_to_remove)
                all_structures.append({"structure": s_new, "energy": energy})

            count += 1
            if count % 10 == 0:
                timenow = time.time()
                self.logger.debug("{} structures, {:.2f} seconds."
                                  .format(count, timenow - starttime))
                self.logger.debug("Average time per combi = {} seconds"
                                  .format((timenow - starttime) / count))
                self.logger.debug("{} symmetrically distinct structures found."
                                  .format(len(all_structures)))

        self.logger.debug("Total symmetrically distinct structures found = {}"
                          .format(len(all_structures)))
        all_structures = sorted(all_structures, key=lambda s: s["energy"])
        return all_structures

    def fast_ordering(self, structure, num_remove_dict, num_to_return=1):
        """
        This method uses the matrix form of ewaldsum to calculate the ewald
        sums of the potential structures. This is on the order of 4 orders of
        magnitude faster when there are large numbers of permutations to
        consider. There are further optimizations possible (doing a smarter
        search of permutations for example), but this wont make a difference
        until the number of permutations is on the order of 30,000.
        """
        self.logger.debug("Performing fast ordering")
        starttime = time.time()
        self.logger.debug("Performing initial ewald sum...")

        ewaldmatrix = EwaldSummation(structure).total_energy_matrix
        self.logger.debug("Ewald sum took {} seconds."
                          .format(time.time() - starttime))
        starttime = time.time()
        m_list = []
        for indices, num in num_remove_dict.items():
            m_list.append([0, num, list(indices), None])

        self.logger.debug("Calling EwaldMinimizer...")
        minimizer = EwaldMinimizer(ewaldmatrix, m_list, num_to_return,
                                   PartialRemoveSitesTransformation.ALGO_FAST)
        self.logger.debug("Minimizing Ewald took {} seconds."
                          .format(time.time() - starttime))

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
                {"energy": output[0],
                 "energy_above_minimum": (output[0] - lowest_energy)
                    / num_atoms,
                 "structure": struct})

        return all_structures

    def enumerate_ordering(self, structure):
        # Generate the disordered structure first.
        s = structure.copy()
        for indices, fraction in zip(self.indices, self.fractions):
            for ind in indices:
                new_sp = {sp: occu * fraction
                          for sp, occu
                          in structure[ind].species_and_occu.items()}
                s[ind] = new_sp
        # Perform enumeration
        from pymatgen.transformations.advanced_transformations import \
            EnumerateStructureTransformation
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
                raise ValueError("Fraction to remove must be consistent with "
                                 "integer amounts in structure.")
            else:
                num_to_remove = int(round(num_to_remove))
            num_remove_dict[tuple(indices)] = num_to_remove
            n = len(indices)
            total_combis += int(round(math.factorial(n) /
                                      math.factorial(num_to_remove) /
                                      math.factorial(n - num_to_remove)))

        self.logger.debug("Total combinations = {}".format(total_combis))

        try:
            num_to_return = int(return_ranked_list)
        except ValueError:
            num_to_return = 1

        num_to_return = max(1, num_to_return)
        self.logger.debug("Will return {} best structures."
                          .format(num_to_return))

        if self.algo == PartialRemoveSitesTransformation.ALGO_FAST:
            all_structures = self.fast_ordering(structure, num_remove_dict,
                                                num_to_return)
        elif self.algo == PartialRemoveSitesTransformation.ALGO_COMPLETE:
            all_structures = self.complete_ordering(structure, num_remove_dict)
        elif self.algo == PartialRemoveSitesTransformation.ALGO_BEST_FIRST:
            all_structures = self.best_first_ordering(structure,
                                                      num_remove_dict)
        elif self.algo == PartialRemoveSitesTransformation.ALGO_ENUMERATE:
            all_structures = self.enumerate_ordering(structure)
        else:
            raise ValueError("Invalid algo.")

        opt_s = all_structures[0]["structure"]
        return opt_s if not return_ranked_list \
            else all_structures[0:num_to_return]

    def __str__(self):
        return "PartialRemoveSitesTransformation : Indices and fraction" + \
               " to remove = {}, ALGO = {}".format(self.indices, self.algo)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        return None

    @property
    def is_one_to_many(self):
        return True
