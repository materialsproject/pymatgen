#!/usr/bin/env python

"""
This module implements classes to perform bond valence analyses.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Oct 26, 2012"

import collections
import json
import itertools
import os
import operator
from math import exp, sqrt

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.core.structure_modifier import StructureEditor

#Let's initialize some module level properties.

#List of electronegative elements specified in M. O'Keefe, & N. Brese,
#JACS, 1991, 113(9), 3226-3229. doi:10.1021/ja00009a002.
ELECTRONEG = map(Element, ["H",
                           "B", "C", "Si",
                           "N", "P", "As", "Sb",
                           "O", "S", "Se", "Te",
                           "F", "Cl", "Br", "I"])

module_dir = os.path.dirname(os.path.abspath(__file__))

#Read in BV parameters.
BV_PARAMS = {}
with open(os.path.join(module_dir, "bvparam_1991.json"), "r") as f:
    for k, v in json.load(f).items():
        BV_PARAMS[Element(k)] = v

#Read in json containing data-mined ICSD BV data.
with open(os.path.join(module_dir, "icsd_bv.json"), "r") as f:
    all_data = json.load(f)
    ICSD_BV_DATA = {Specie.from_string(sp): data
                    for sp, data in all_data["bvsum"].items()}
    PRIOR_PROB = {Specie.from_string(sp): data
                  for sp, data in all_data["occurrence"].items()}


def calculate_bv_sum(site, nn_list, scale_factor=1):
    """
    Calculates the BV sum of a site.

    Args:
        site:
            The site
        nn_list:
            List of nearest neighbors in the format [(nn_site, dist), ...].
        anion_el:
            The most electronegative element in the structure.
        scale_factor:
            A scale factor to be applied. This is useful for scaling distance,
            esp in the case of calculation-relaxed structures which may tend
            to under (GGA) or over bind (LDA).
    """
    el1 = Element(site.specie.symbol)
    bvsum = 0
    for (nn, dist) in nn_list:
        el2 = Element(nn.specie.symbol)
        if (el1 in ELECTRONEG or el2 in ELECTRONEG) and el1 != el2:
            r1 = BV_PARAMS[el1]["r"]
            r2 = BV_PARAMS[el2]["r"]
            c1 = BV_PARAMS[el1]["c"]
            c2 = BV_PARAMS[el2]["c"]
            R = r1 + r2 - r1 * r2 * (sqrt(c1) - sqrt(c2)) ** 2 / \
                (c1 * r1 + c2 * r2)
            vij = exp((R - dist * scale_factor) / 0.31)
            bvsum += vij * (1 if el1.X < el2.X else -1)
    return bvsum


class BVAnalyzer(object):
    """
    This class implements a maximum a posteriori (MAP) estimation method to
    determine oxidation states in a structure. The algorithm is as follows:
    1) The bond valence sum of all symmetrically distinct sites in a structure
    is calculated using the element-based parameters in M. O'Keefe, & N. Brese,
    JACS, 1991, 113(9), 3226-3229. doi:10.1021/ja00009a002.
    2) The posterior probabilities of all oxidation states is then calculated
    using: P(oxi_state/BV) = K * P(BV/oxi_state) * P(oxi_state), where K is
    a constant factor for each element. P(BV/oxi_state) is calculated as a
    Gaussian with mean and std deviation determined from an analysis of
    the ICSD. The posterior P(oxi_state) is determined from a frequency
    analysis of the ICSD.
    3) The oxidation states are then ranked in order of decreasing probability
    and the oxidation state combination that result in a charge neutral cell
    is selected.
    """

    def __init__(self, symm_tol=0.1, max_radius=4, max_permutations=100000,
                 distance_scale_factor=1.015):
        """
        Args:
            symm_tol:
                Symmetry tolerance used to determine which sites are
                symmetrically equivalent. Set to 0 to turn off symmetry.
            max_radius:
                Maximum radius in Angstrom used to find nearest neighbors.
            max_permutations:
                The maximum number of permutations of oxidation states to test.
            distance_scale_factor:
                A scale factor to be applied. This is useful for scaling
                distances, esp in the case of calculation-relaxed structures
                which may tend to under (GGA) or over bind (LDA). The default
                of 1.015 works for GGA. For experimental structure, set this to
                1.
        """
        self.symm_tol = symm_tol
        self.max_radius = max_radius
        self.max_permutations = max_permutations
        self.dist_scale_factor = distance_scale_factor

    def _calc_site_probabilities(self, site, nn):
        el = site.specie.symbol
        bv_sum = calculate_bv_sum(site, nn,
                                  scale_factor=self.dist_scale_factor)
        prob = {}
        for sp, data in ICSD_BV_DATA.items():
            if sp.symbol == el and sp.oxi_state != 0 and data["std"] > 0:
                u = data["mean"]
                sigma = data["std"]
                #Calculate posterior probability. Note that constant
                #factors are ignored. They have no effect on the results.
                prob[sp.oxi_state] = exp(-(bv_sum - u) ** 2 /
                                         2 / (sigma ** 2)) \
                                         / sigma * PRIOR_PROB[sp]
        #Normalize the probabilities
        prob = {k: v / sum(prob.values()) for k, v in prob.items()}
        return prob

    def get_valences(self, structure):
        """
        Returns a list of valences for the structure. This currently works only
        for ordered structures only.

        Args:
            structure:
                Structure to analyze

        Returns:
            A list of valences for each site in the structure,
            e.g., [1, 1, -2].

        Raises:
            A ValueError is the valences cannot be determined.
        """
        els = [Element(el.symbol) for el in structure.composition.elements]

        if not set(els).issubset(set(BV_PARAMS.keys())):
            raise ValueError(
                "Structure contains elements not in set of BV parameters!"
            )

        #Perform symmetry determination and get sites grouped by symmetry.
        if self.symm_tol:
            finder = SymmetryFinder(structure, self.symm_tol)
            symm_structure = finder.get_symmetrized_structure()
            equi_sites = symm_structure.equivalent_sites
        else:
            equi_sites = [[site] for site in structure]

        #Sort the equivalent sites by decreasing electronegativity.
        equi_sites = sorted(equi_sites, key=lambda sites:-sites[0].specie.X)

        #Get a list of valences and probabilities for each symmetrically
        #distinct site.
        valences = []
        all_prob = []
        for sites in equi_sites:
            test_site = sites[0]
            nn = structure.get_neighbors(test_site, self.max_radius)
            prob = self._calc_site_probabilities(test_site, nn)
            all_prob.append(prob)
            val = list(prob.keys())
            #Sort valences in order of decreasing probability.
            val = sorted(val, key=lambda v:-prob[v])
            #Retain probabilities that are at least 1/100 of highest prob.
            valences.append(filter(lambda v: prob[v] > 0.01 * prob[val[0]],
                                   val))

        #Based on the max allowed permutations, determine the number of
        #candidates per site.
        num_perm = 0
        selected_valences = [[v.pop(0)] for v in valences]
        while num_perm < self.max_permutations:
            max_prob = 0
            ind = -1
            for i, v in enumerate(valences):
                if len(v) > 0 and all_prob[i][v[0]] > max_prob:
                    max_prob = v[0]
                    ind = i
            if ind == -1:
                break
            else:
                selected_valences[ind].append(valences[ind].pop(0))
            num_perm = reduce(operator.mul, map(len, selected_valences))

        scores = {}
        #Find valid valence combinations and score them by total probability.
        for v_set in itertools.product(*selected_valences):
            total = 0
            el_oxi = collections.defaultdict(list)
            for i, sites in enumerate(equi_sites):
                total += v_set[i] * len(sites)
                el_oxi[sites[0].specie.symbol].append(v_set[i])

            #Calculate the maximum range in oxidation states for each element.
            max_diff = max([max(v) - min(v) for v in el_oxi.values()])

            if total == 0 and max_diff <= 1:
                #Cell has to be charge neutral. And the maximum difference in
                #oxidation state for each element cannot exceed 1.
                score = reduce(operator.mul, [all_prob[i][v]
                                              for i, v in enumerate(v_set)])
                scores[tuple(v_set)] = score

        if scores:
            best = max(scores.keys(), key=lambda k: scores[k])
            assigned = {}
            for val, sites in zip(best, equi_sites):
                for site in sites:
                    assigned[site] = val

            return [int(assigned[site]) for site in structure]
        else:
            raise ValueError("Valences cannot be assigned!")

    def get_oxi_state_decorated_structure(self, structure):
        """
        Get an oxidation state decorated structure. This currently works only
        for ordered structures only.

        Args:
            structure:
                Structure to analyze

        Returns:
            A modified structure that is oxidation state decorated.

        Raises:
            A ValueError is the valences cannot be determined.
        """
        valences = self.get_valences(structure)
        editor = StructureEditor(structure)
        editor.add_oxidation_state_by_site(valences)
        return editor.modified_structure
