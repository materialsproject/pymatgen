"""This module provides classes for representing species substitution probabilities."""

from __future__ import annotations

import functools
import itertools
import json
import logging
import math
import os
from collections import defaultdict
from operator import mul
from typing import TYPE_CHECKING

from monty.design_patterns import cached_class

from pymatgen.core import Species, get_el_sp
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing_extensions import Self

    from pymatgen.util.typing import SpeciesLike

__author__ = "Will Richards, Geoffroy Hautier"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.2"
__maintainer__ = "Will Richards"
__email__ = "wrichard@mit.edu"
__date__ = "Aug 31, 2012"


@due.dcite(
    Doi("10.1021/ic102031h"),
    description="Data Mined Ionic Substitutions for the Discovery of New Compounds",
)
@cached_class
class SubstitutionProbability:
    """
    This class finds substitution probabilities given lists of atoms
    to substitute. The inputs make more sense if you look through the
    from_defaults static method.

    The substitution prediction algorithm is presented in:
    Hautier, G., Fischer, C., Ehrlacher, V., Jain, A., and Ceder, G. (2011)
    Data Mined Ionic Substitutions for the Discovery of New Compounds.
    Inorganic Chemistry, 50(2), 656-663. doi:10.1021/ic102031h
    """

    def __init__(self, lambda_table=None, alpha=-5):
        """
        Args:
            lambda_table: JSON table of the weight functions lambda if None,
                will use the default lambda.json table
            alpha (float): weight function for never observed substitutions.
        """
        if lambda_table is not None:
            self._lambda_table = lambda_table
        else:
            module_dir = os.path.dirname(__file__)
            json_file = f"{module_dir}/data/lambda.json"
            with open(json_file) as file:
                self._lambda_table = json.load(file)

        # build map of specie pairs to lambdas
        self.alpha = alpha
        self._l = {}
        self.species = set()
        for row in self._lambda_table:
            if "D1+" not in row:
                s1 = Species.from_str(row[0])
                s2 = Species.from_str(row[1])
                self.species.add(s1)
                self.species.add(s2)
                self._l[frozenset([s1, s2])] = float(row[2])

        # create Z and px
        self.Z = 0
        self._px: dict[SpeciesLike, float] = defaultdict(float)
        for s1, s2 in itertools.product(self.species, repeat=2):
            value = math.exp(self.get_lambda(s1, s2))
            self._px[s1] += value / 2
            self._px[s2] += value / 2
            self.Z += value

    def get_lambda(self, s1, s2):
        """
        Args:
            s1 (SpeciesLike): Ion in 1st structure.
            s2 (SpeciesLike): Ion in 2nd structure.

        Returns:
            Lambda values
        """
        key = frozenset([get_el_sp(s1), get_el_sp(s2)])
        return self._l.get(key, self.alpha)

    def get_px(self, sp: SpeciesLike) -> float:
        """
        Args:
            sp (SpeciesLike): Species.

        Returns:
            float: Probability
        """
        return self._px[get_el_sp(sp)]

    def prob(self, s1, s2):
        """Get the probability of 2 species substitution. Not used by the
        structure predictor.

        Returns:
            Probability of s1 and s2 substitution.
        """
        return math.exp(self.get_lambda(s1, s2)) / self.Z

    def cond_prob(self, s1, s2):
        """
        Conditional probability of substituting s1 for s2.

        Args:
            s1:
                The *variable* specie
            s2:
                The *fixed* specie

        Returns:
            Conditional probability used by structure predictor.
        """
        return math.exp(self.get_lambda(s1, s2)) / self.get_px(s2)

    def pair_corr(self, s1, s2):
        """
        Pair correlation of two species.

        Returns:
            The pair correlation of 2 species
        """
        return math.exp(self.get_lambda(s1, s2)) * self.Z / (self.get_px(s1) * self.get_px(s2))

    def cond_prob_list(self, l1, l2):
        """Find the probabilities of 2 lists. These should include ALL species.
        This is the probability conditional on l2.

        Args:
            l1, l2:
                lists of species

        Returns:
            The conditional probability (assuming these species are in
            l2)
        """
        if len(l1) != len(l2):
            raise ValueError("lengths of l1 and l2 mismatch.")
        p = 1
        for s1, s2 in zip(l1, l2, strict=True):
            p *= self.cond_prob(s1, s2)
        return p

    def as_dict(self):
        """Get MSONable dict."""
        return {
            "name": type(self).__name__,
            "version": __version__,
            "init_args": {"lambda_table": self._l, "alpha": self.alpha},
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            Class
        """
        return cls(**dct["init_args"])


class SubstitutionPredictor:
    """
    Predicts likely substitutions either to or from a given composition
    or species list using the SubstitutionProbability.
    """

    def __init__(self, lambda_table=None, alpha=-5, threshold=1e-3):
        """
        Args:
            lambda_table (dict): Input lambda table.
            alpha (float): weight function for never observed substitutions
            threshold (float): Threshold to use to identify high probability structures.
        """
        self.p = SubstitutionProbability(lambda_table, alpha)
        self.threshold = threshold

    def list_prediction(self, species, to_this_composition=True):
        """
        Args:
            species:
                list of species
            to_this_composition:
                If true, substitutions with this as a final composition
                will be found. If false, substitutions with this as a
                starting composition will be found (these are slightly
                different).

        Returns:
            List of predictions in the form of dictionaries.
            If to_this_composition is true, the values of the dictionary
            will be from the list species. If false, the keys will be
            from that list.
        """
        for sp in species:
            if get_el_sp(sp) not in self.p.species:
                raise ValueError(f"the species {sp} is not allowed for the probability model you are using")
        max_probabilities = []
        for s1 in species:
            if to_this_composition:
                max_p = max(self.p.cond_prob(s2, s1) for s2 in self.p.species)
            else:
                max_p = max(self.p.cond_prob(s1, s2) for s2 in self.p.species)
            max_probabilities.append(max_p)

        output = []

        def _recurse(output_prob, output_species):
            best_case_prob = list(max_probabilities)
            best_case_prob[: len(output_prob)] = output_prob
            if functools.reduce(mul, best_case_prob) > self.threshold:
                if len(output_species) == len(species):
                    odict = {"probability": functools.reduce(mul, best_case_prob)}
                    if to_this_composition:
                        odict["substitutions"] = dict(zip(output_species, species, strict=True))
                    else:
                        odict["substitutions"] = dict(zip(species, output_species, strict=True))
                    if len(output_species) == len(set(output_species)):
                        output.append(odict)
                    return
                for sp in self.p.species:
                    i = len(output_prob)
                    prob = self.p.cond_prob(sp, species[i]) if to_this_composition else self.p.cond_prob(species[i], sp)
                    _recurse([*output_prob, prob], [*output_species, sp])

        _recurse([], [])
        logging.info(f"{len(output)} substitutions found")
        return output

    def composition_prediction(self, composition, to_this_composition=True):
        """Get charged balanced substitutions from a starting or ending
        composition.

        Args:
            composition:
                starting or ending composition
            to_this_composition:
                If true, substitutions with this as a final composition
                will be found. If false, substitutions with this as a
                starting composition will be found (these are slightly
                different)

        Returns:
            List of predictions in the form of dictionaries.
            If to_this_composition is true, the values of the dictionary
            will be from the list species. If false, the keys will be
            from that list.
        """
        preds = self.list_prediction(list(composition), to_this_composition)
        output = []
        for p in preds:
            subs = {v: k for k, v in p["substitutions"].items()} if to_this_composition else p["substitutions"]
            charge = 0
            for k, v in composition.items():
                charge += subs[k].oxi_state * v
            if abs(charge) < 1e-8:
                output.append(p)
        logging.info(f"{len(output)} charge balanced substitutions found")
        return output
