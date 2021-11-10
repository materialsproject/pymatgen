# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes for predicting new structures from existing ones.
"""

import functools
import itertools
import logging
from operator import mul

from monty.json import MSONable

from pymatgen.alchemy.filters import RemoveDuplicatesFilter, RemoveExistingFilter
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.alchemy.transmuters import StandardTransmuter
from pymatgen.analysis.structure_prediction.substitution_probability import (
    SubstitutionProbability,
)
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.transformations.standard_transformations import SubstitutionTransformation

__author__ = "Will Richards, Geoffroy Hautier"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.2"
__maintainer__ = "Will Richards"
__email__ = "wrichard@mit.edu"
__date__ = "Aug 31, 2012"


class Substitutor(MSONable):
    """
    This object uses a data mined ionic substitution approach to propose
    compounds likely to be stable. It relies on an algorithm presented in
    Hautier, G., Fischer, C., Ehrlacher, V., Jain, A., and Ceder, G. (2011).
    Data Mined Ionic Substitutions for the Discovery of New Compounds.
    Inorganic Chemistry, 50(2), 656-663. doi:10.1021/ic102031h
    """

    def __init__(self, threshold=1e-3, symprec=0.1, **kwargs):
        """
        This substitutor uses the substitution probability class to
        find good substitutions for a given chemistry or structure.

        Args:
            threshold:
                probability threshold for predictions
            symprec:
                symmetry precision to determine if two structures
                are duplicates
            kwargs:
                kwargs for the SubstitutionProbability object
                lambda_table, alpha
        """
        self._kwargs = kwargs
        self._sp = SubstitutionProbability(**kwargs)
        self._threshold = threshold
        self._symprec = symprec

    def get_allowed_species(self):
        """
        returns the species in the domain of the probability function
        any other specie will not work
        """
        return self._sp.species

    def pred_from_structures(
        self,
        target_species,
        structures_list,
        remove_duplicates=True,
        remove_existing=False,
    ):
        """
        performs a structure prediction targeting compounds containing all of
        the target_species, based on a list of structure (those structures
        can for instance come from a database like the ICSD). It will return
        all the structures formed by ionic substitutions with a probability
        higher than the threshold

        Notes:
        If the default probability model is used, input structures must
        be oxidation state decorated. See AutoOxiStateDecorationTransformation

        This method does not change the number of species in a structure. i.e
        if the number of target species is 3, only input structures containing
        3 species will be considered.

        Args:
            target_species:
                a list of species with oxidation states
                e.g., [Species('Li',1),Species('Ni',2), Species('O',-2)]

            structures_list:
                a list of dictionnary of the form {'structure':Structure object
                ,'id':some id where it comes from}
                the id can for instance refer to an ICSD id.

            remove_duplicates:
                if True, the duplicates in the predicted structures will
                be removed

            remove_existing:
                if True, the predicted structures that already exist in the
                structures_list will be removed

        Returns:
            a list of TransformedStructure objects.
        """
        target_species = [get_el_sp(sp) for sp in target_species]
        result = []
        transmuter = StandardTransmuter([])
        if len(list(set(target_species) & set(self.get_allowed_species()))) != len(target_species):
            raise ValueError(
                "the species in target_species are not allowed " + "for the probability model you are using"
            )

        for permut in itertools.permutations(target_species):
            for s in structures_list:
                # check if: species are in the domain,
                # and the probability of subst. is above the threshold
                els = s["structure"].composition.elements
                if (
                    len(els) == len(permut)
                    and len(list(set(els) & set(self.get_allowed_species()))) == len(els)
                    and self._sp.cond_prob_list(permut, els) > self._threshold
                ):

                    clean_subst = {els[i]: permut[i] for i in range(0, len(els)) if els[i] != permut[i]}

                    if len(clean_subst) == 0:
                        continue

                    transf = SubstitutionTransformation(clean_subst)

                    if Substitutor._is_charge_balanced(transf.apply_transformation(s["structure"])):
                        ts = TransformedStructure(
                            s["structure"],
                            [transf],
                            history=[{"source": s["id"]}],
                            other_parameters={
                                "type": "structure_prediction",
                                "proba": self._sp.cond_prob_list(permut, els),
                            },
                        )
                        result.append(ts)
                        transmuter.append_transformed_structures([ts])

        if remove_duplicates:
            transmuter.apply_filter(RemoveDuplicatesFilter(symprec=self._symprec))
        if remove_existing:
            # Make the list of structures from structures_list that corresponds to the
            # target species
            chemsys = {sp.symbol for sp in target_species}
            structures_list_target = [
                st["structure"]
                for st in structures_list
                if Substitutor._is_from_chemical_system(chemsys, st["structure"])
            ]
            transmuter.apply_filter(RemoveExistingFilter(structures_list_target, symprec=self._symprec))
        return transmuter.transformed_structures

    @staticmethod
    def _is_charge_balanced(struct):
        """
        checks if the structure object is charge balanced
        """
        return sum([s.specie.oxi_state for s in struct.sites]) == 0.0

    @staticmethod
    def _is_from_chemical_system(chemical_system, struct):
        """
        checks if the structure object is from the given chemical system
        """
        return {sp.symbol for sp in struct.composition} == set(chemical_system)

    def pred_from_list(self, species_list):
        """
        There are an exceptionally large number of substitutions to
        look at (260^n), where n is the number of species in the
        list. We need a more efficient than brute force way of going
        through these possibilities. The brute force method would be::

            output = []
            for p in itertools.product(self._sp.species_list
                                       , repeat = len(species_list)):
                if self._sp.conditional_probability_list(p, species_list)
                                       > self._threshold:
                    output.append(dict(zip(species_list,p)))
            return output

        Instead of that we do a branch and bound.

        Args:
            species_list:
                list of species in the starting structure

        Returns:
            list of dictionaries, each including a substitutions
            dictionary, and a probability value
        """
        species_list = [get_el_sp(sp) for sp in species_list]
        # calculate the highest probabilities to help us stop the recursion
        max_probabilities = []
        for s2 in species_list:
            max_p = 0
            for s1 in self._sp.species:
                max_p = max([self._sp.cond_prob(s1, s2), max_p])
            max_probabilities.append(max_p)
        output = []

        def _recurse(output_prob, output_species):
            best_case_prob = list(max_probabilities)
            best_case_prob[: len(output_prob)] = output_prob
            if functools.reduce(mul, best_case_prob) > self._threshold:
                if len(output_species) == len(species_list):
                    odict = {
                        "substitutions": dict(zip(species_list, output_species)),
                        "probability": functools.reduce(mul, best_case_prob),
                    }
                    output.append(odict)
                    return
                for sp in self._sp.species:
                    i = len(output_prob)
                    prob = self._sp.cond_prob(sp, species_list[i])
                    _recurse(output_prob + [prob], output_species + [sp])

        _recurse([], [])
        logging.info("{} substitutions found".format(len(output)))
        return output

    def pred_from_comp(self, composition):
        """
        Similar to pred_from_list except this method returns a list after
        checking that compositions are charge balanced.
        """
        output = []
        predictions = self.pred_from_list(composition.elements)
        for p in predictions:
            subs = p["substitutions"]
            charge = 0
            for i_el in composition.elements:
                f_el = subs[i_el]
                charge += f_el.oxi_state * composition[i_el]
            if charge == 0:
                output.append(p)
        logging.info("{} charge balanced compositions found".format(len(output)))
        return output

    def as_dict(self):
        """
        Returns: MSONable dict
        """
        return {
            "name": self.__class__.__name__,
            "version": __version__,
            "kwargs": self._kwargs,
            "threshold": self._threshold,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
        }

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): Dict representation

        Returns:
            Class
        """
        t = d["threshold"]
        kwargs = d["kwargs"]
        return cls(threshold=t, **kwargs)
