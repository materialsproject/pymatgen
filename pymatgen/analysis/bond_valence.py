#!/usr/bin/env python

"""
This module implements classes to perform bond valence analyses.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Oct 26, 2012"

import collections
import json
import numpy as np
import operator
import os

from math import exp, sqrt
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Structure
from pymatgen.symmetry.finder import SymmetryFinder

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
                prob[sp.oxi_state] = exp(-(bv_sum - u) ** 2 / 2 /
                                         (sigma ** 2)) \
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
        equi_sites = sorted(equi_sites,
                            key=lambda sites: -sites[0].species_and_occu
                            .average_electroneg)

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
            val = sorted(val, key=lambda v: -prob[v])
            #Retain probabilities that are at least 1/100 of highest prob.
            valences.append(filter(lambda v: prob[v] > 0.01 * prob[val[0]],
                                   val))

        #make variables needed for recursion
        nsites = np.array(map(len, equi_sites))
        vmin = np.array(map(min, valences))
        vmax = np.array(map(max, valences))

        self._n = 0
        self._best_score = 0
        self._best_vset = None

        def evaluate_assignment(v_set):
            el_oxi = collections.defaultdict(list)
            for i, sites in enumerate(equi_sites):
                el_oxi[sites[0].specie.symbol].append(v_set[i])
            max_diff = max([max(v) - min(v) for v in el_oxi.values()])
            if max_diff > 1:
                return
            score = reduce(operator.mul, [all_prob[i][v]
                                          for i, v in enumerate(v_set)])
            if score > self._best_score:
                self._best_vset = v_set
                self._best_score = score

        def _recurse(assigned=[]):
            #recurses to find permutations of valences based on whether a
            #charge balanced assignment can still be found
            if self._n > self.max_permutations:
                return

            i = len(assigned)
            highest = vmax.copy()
            highest[:i] = assigned
            highest *= nsites
            highest = np.sum(highest)

            lowest = vmin.copy()
            lowest[:i] = assigned
            lowest *= nsites
            lowest = np.sum(lowest)

            if highest < 0 or lowest > 0:
                self._n += 1
                return

            if i == len(valences):
                evaluate_assignment(assigned)
                self._n += 1
                return
            else:
                for v in valences[i]:
                    new_assigned = list(assigned)
                    _recurse(new_assigned + [v])

        _recurse()

        if self._best_vset:
            assigned = {}
            for val, sites in zip(self._best_vset, equi_sites):
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
        s = Structure.from_sites(structure.sites)
        s.add_oxidation_state_by_site(valences)
        return s
