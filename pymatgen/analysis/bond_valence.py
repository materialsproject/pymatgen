# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements classes to perform bond valence analyses.
"""

from __future__ import annotations

import collections
import functools
import operator
import os
from math import exp, sqrt

import numpy as np
from monty.serialization import loadfn

from pymatgen.core.periodic_table import Element, Species, get_el_sp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Let's initialize some module level properties.

# List of electronegative elements specified in M. O'Keefe, & N. Brese,
# JACS, 1991, 113(9), 3226-3229. doi:10.1021/ja00009a002.
ELECTRONEG = [
    Element(sym)
    for sym in [
        "H",
        "B",
        "C",
        "Si",
        "N",
        "P",
        "As",
        "Sb",
        "O",
        "S",
        "Se",
        "Te",
        "F",
        "Cl",
        "Br",
        "I",
    ]
]

module_dir = os.path.dirname(os.path.abspath(__file__))

# Read in BV parameters.
BV_PARAMS = {}
for k, v in loadfn(os.path.join(module_dir, "bvparam_1991.yaml")).items():
    BV_PARAMS[Element(k)] = v

# Read in yaml containing data-mined ICSD BV data.
all_data = loadfn(os.path.join(module_dir, "icsd_bv.yaml"))
ICSD_BV_DATA = {Species.from_string(sp): data for sp, data in all_data["bvsum"].items()}
PRIOR_PROB = {Species.from_string(sp): data for sp, data in all_data["occurrence"].items()}


def calculate_bv_sum(site, nn_list, scale_factor=1.0):
    """
    Calculates the BV sum of a site.

    Args:
        site (PeriodicSite): The central site to calculate the bond valence
        nn_list ([Neighbor]): A list of namedtuple Neighbors having "distance"
            and "site" attributes
        scale_factor (float): A scale factor to be applied. This is useful for
            scaling distance, esp in the case of calculation-relaxed structures
            which may tend to under (GGA) or over bind (LDA).
    """
    el1 = Element(site.specie.symbol)
    bvsum = 0
    for nn in nn_list:
        el2 = Element(nn.specie.symbol)
        if (el1 in ELECTRONEG or el2 in ELECTRONEG) and el1 != el2:
            r1 = BV_PARAMS[el1]["r"]
            r2 = BV_PARAMS[el2]["r"]
            c1 = BV_PARAMS[el1]["c"]
            c2 = BV_PARAMS[el2]["c"]
            R = r1 + r2 - r1 * r2 * (sqrt(c1) - sqrt(c2)) ** 2 / (c1 * r1 + c2 * r2)
            vij = exp((R - nn.nn_distance * scale_factor) / 0.31)
            bvsum += vij * (1 if el1.X < el2.X else -1)
    return bvsum


def calculate_bv_sum_unordered(site, nn_list, scale_factor=1):
    """
    Calculates the BV sum of a site for unordered structures.

    Args:
        site (PeriodicSite): The central site to calculate the bond valence
        nn_list ([Neighbor]): A list of namedtuple Neighbors having "distance"
            and "site" attributes
        scale_factor (float): A scale factor to be applied. This is useful for
            scaling distance, esp in the case of calculation-relaxed structures
            which may tend to under (GGA) or over bind (LDA).
    """
    # If the site "site" has N partial occupations as : f_{site}_0,
    # f_{site}_1, ... f_{site}_N of elements
    # X_{site}_0, X_{site}_1, ... X_{site}_N, and each neighbors nn_i in nn
    # has N_{nn_i} partial occupations as :
    # f_{nn_i}_0, f_{nn_i}_1, ..., f_{nn_i}_{N_{nn_i}}, then the bv sum of
    # site "site" is obtained as :
    # \sum_{nn} \sum_j^N \sum_k^{N_{nn}} f_{site}_j f_{nn_i}_k vij_full
    # where vij_full is the valence bond of the fully occupied bond
    bvsum = 0
    for specie1, occu1 in site.species.items():
        el1 = Element(specie1.symbol)
        for nn in nn_list:
            for specie2, occu2 in nn.species.items():
                el2 = Element(specie2.symbol)
                if (el1 in ELECTRONEG or el2 in ELECTRONEG) and el1 != el2:
                    r1 = BV_PARAMS[el1]["r"]
                    r2 = BV_PARAMS[el2]["r"]
                    c1 = BV_PARAMS[el1]["c"]
                    c2 = BV_PARAMS[el2]["c"]
                    R = r1 + r2 - r1 * r2 * (sqrt(c1) - sqrt(c2)) ** 2 / (c1 * r1 + c2 * r2)
                    vij = exp((R - nn.nn_distance * scale_factor) / 0.31)
                    bvsum += occu1 * occu2 * vij * (1 if el1.X < el2.X else -1)
    return bvsum


class BVAnalyzer:
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

    CHARGE_NEUTRALITY_TOLERANCE = 0.00001

    def __init__(
        self,
        symm_tol=0.1,
        max_radius=4,
        max_permutations=100000,
        distance_scale_factor=1.015,
        charge_neutrality_tolerance=CHARGE_NEUTRALITY_TOLERANCE,
        forbidden_species=None,
    ):
        """
        Initializes the BV analyzer, with useful defaults.

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
            charge_neutrality_tolerance:
                Tolerance on the charge neutrality when unordered structures
                are at stake.
            forbidden_species:
                List of species that are forbidden (example : ["O-"] cannot be
                used) It is used when e.g. someone knows that some oxidation
                state cannot occur for some atom in a structure or list of
                structures.
        """
        self.symm_tol = symm_tol
        self.max_radius = max_radius
        self.max_permutations = max_permutations
        self.dist_scale_factor = distance_scale_factor
        self.charge_neutrality_tolerance = charge_neutrality_tolerance
        forbidden_species = [get_el_sp(sp) for sp in forbidden_species] if forbidden_species else []
        self.icsd_bv_data = (
            {get_el_sp(specie): data for specie, data in ICSD_BV_DATA.items() if specie not in forbidden_species}
            if len(forbidden_species) > 0
            else ICSD_BV_DATA
        )

    def _calc_site_probabilities(self, site, nn):
        el = site.specie.symbol
        bv_sum = calculate_bv_sum(site, nn, scale_factor=self.dist_scale_factor)
        prob = {}
        for sp, data in self.icsd_bv_data.items():
            if sp.symbol == el and sp.oxi_state != 0 and data["std"] > 0:
                u = data["mean"]
                sigma = data["std"]
                # Calculate posterior probability. Note that constant
                # factors are ignored. They have no effect on the results.
                prob[sp.oxi_state] = exp(-((bv_sum - u) ** 2) / 2 / (sigma**2)) / sigma * PRIOR_PROB[sp]
        # Normalize the probabilities
        try:
            prob = {k: v / sum(prob.values()) for k, v in prob.items()}
        except ZeroDivisionError:
            prob = {k: 0.0 for k in prob}
        return prob

    def _calc_site_probabilities_unordered(self, site, nn):
        bv_sum = calculate_bv_sum_unordered(site, nn, scale_factor=self.dist_scale_factor)
        prob = {}
        for specie in site.species:
            el = specie.symbol

            prob[el] = {}
            for sp, data in self.icsd_bv_data.items():
                if sp.symbol == el and sp.oxi_state != 0 and data["std"] > 0:
                    u = data["mean"]
                    sigma = data["std"]
                    # Calculate posterior probability. Note that constant
                    # factors are ignored. They have no effect on the results.
                    prob[el][sp.oxi_state] = exp(-((bv_sum - u) ** 2) / 2 / (sigma**2)) / sigma * PRIOR_PROB[sp]
            # Normalize the probabilities
            try:
                prob[el] = {k: v / sum(prob[el].values()) for k, v in prob[el].items()}
            except ZeroDivisionError:
                prob[el] = {k: 0.0 for k in prob[el]}
        return prob

    def get_valences(self, structure):
        """
        Returns a list of valences for the structure. This currently works only
        for ordered structures only.

        Args:
            structure: Structure to analyze

        Returns:
            A list of valences for each site in the structure (for an ordered
            structure), e.g., [1, 1, -2] or a list of lists with the
            valences for each fractional element of each site in the
            structure (for an unordered structure),
            e.g., [[2, 4], [3], [-2], [-2], [-2]]

        Raises:
            A ValueError if the valences cannot be determined.
        """
        els = [Element(el.symbol) for el in structure.composition.elements]

        if not set(els).issubset(set(BV_PARAMS)):
            raise ValueError("Structure contains elements not in set of BV parameters!")

        # Perform symmetry determination and get sites grouped by symmetry.
        if self.symm_tol:
            finder = SpacegroupAnalyzer(structure, self.symm_tol)
            symm_structure = finder.get_symmetrized_structure()
            equi_sites = symm_structure.equivalent_sites
        else:
            equi_sites = [[site] for site in structure]

        # Sort the equivalent sites by decreasing electronegativity.
        equi_sites = sorted(equi_sites, key=lambda sites: -sites[0].species.average_electroneg)

        # Get a list of valences and probabilities for each symmetrically
        # distinct site.
        valences = []
        all_prob = []
        if structure.is_ordered:
            for sites in equi_sites:
                test_site = sites[0]
                nn = structure.get_neighbors(test_site, self.max_radius)
                prob = self._calc_site_probabilities(test_site, nn)
                all_prob.append(prob)
                val = list(prob)
                # Sort valences in order of decreasing probability.
                val = sorted(val, key=lambda v: -prob[v])
                # Retain probabilities that are at least 1/100 of highest prob.
                valences.append(list(filter(lambda v: prob[v] > 0.01 * prob[val[0]], val)))
        else:
            full_all_prob = []
            for sites in equi_sites:
                test_site = sites[0]
                nn = structure.get_neighbors(test_site, self.max_radius)
                prob = self._calc_site_probabilities_unordered(test_site, nn)
                all_prob.append(prob)
                full_all_prob.extend(prob.values())
                vals = []
                for (elsp, _) in get_z_ordered_elmap(test_site.species):
                    val = list(prob[elsp.symbol])
                    # Sort valences in order of decreasing probability.
                    val = sorted(val, key=lambda v: -prob[elsp.symbol][v])
                    # Retain probabilities that are at least 1/100 of highest
                    # prob.
                    vals.append(
                        list(
                            filter(
                                lambda v: prob[elsp.symbol][v] > 0.001 * prob[elsp.symbol][val[0]],
                                val,
                            )
                        )
                    )
                valences.append(vals)

        # make variables needed for recursion
        if structure.is_ordered:
            nsites = np.array([len(i) for i in equi_sites])
            vmin = np.array([min(i) for i in valences])
            vmax = np.array([max(i) for i in valences])

            self._n = 0
            self._best_score = 0
            self._best_vset = None

            def evaluate_assignment(v_set):
                el_oxi = collections.defaultdict(list)
                for i, sites in enumerate(equi_sites):
                    el_oxi[sites[0].specie.symbol].append(v_set[i])
                max_diff = max(max(v) - min(v) for v in el_oxi.values())
                if max_diff > 1:
                    return
                score = functools.reduce(operator.mul, [all_prob[i][v] for i, v in enumerate(v_set)])
                if score > self._best_score:
                    self._best_vset = v_set
                    self._best_score = score

            def _recurse(assigned=None):
                # recurses to find permutations of valences based on whether a
                # charge balanced assignment can still be found
                if self._n > self.max_permutations:
                    return None
                if assigned is None:
                    assigned = []

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
                    return None

                if i == len(valences):
                    evaluate_assignment(assigned)
                    self._n += 1
                    return None
                for v in valences[i]:
                    new_assigned = list(assigned)
                    _recurse(new_assigned + [v])
                return None

        else:
            nsites = np.array([len(i) for i in equi_sites])
            tmp = []
            attrib = []
            for insite, nsite in enumerate(nsites):
                for _ in valences[insite]:
                    tmp.append(nsite)
                    attrib.append(insite)
            new_nsites = np.array(tmp)
            fractions = []
            elements = []
            for sites in equi_sites:
                for sp, occu in get_z_ordered_elmap(sites[0].species):
                    elements.append(sp.symbol)
                    fractions.append(occu)
            fractions = np.array(fractions, np.float_)
            new_valences = []
            for vals in valences:
                for val in vals:
                    new_valences.append(val)
            vmin = np.array([min(i) for i in new_valences], np.float_)
            vmax = np.array([max(i) for i in new_valences], np.float_)

            self._n = 0
            self._best_score = 0
            self._best_vset = None

            def evaluate_assignment(v_set):
                el_oxi = collections.defaultdict(list)
                jj = 0
                for sites in equi_sites:
                    for specie, _ in get_z_ordered_elmap(sites[0].species):
                        el_oxi[specie.symbol].append(v_set[jj])
                        jj += 1
                max_diff = max(max(v) - min(v) for v in el_oxi.values())
                if max_diff > 2:
                    return

                score = functools.reduce(
                    operator.mul,
                    [all_prob[attrib[iv]][elements[iv]][vv] for iv, vv in enumerate(v_set)],
                )
                if score > self._best_score:
                    self._best_vset = v_set
                    self._best_score = score

            def _recurse(assigned=None):
                # recurses to find permutations of valences based on whether a
                # charge balanced assignment can still be found
                if self._n > self.max_permutations:
                    return None
                if assigned is None:
                    assigned = []

                i = len(assigned)
                highest = vmax.copy()
                highest[:i] = assigned
                highest *= new_nsites
                highest *= fractions
                highest = np.sum(highest)

                lowest = vmin.copy()
                lowest[:i] = assigned
                lowest *= new_nsites
                lowest *= fractions
                lowest = np.sum(lowest)

                if highest < -self.charge_neutrality_tolerance or lowest > self.charge_neutrality_tolerance:
                    self._n += 1
                    return None

                if i == len(new_valences):
                    evaluate_assignment(assigned)
                    self._n += 1
                    return None

                for v in new_valences[i]:
                    new_assigned = list(assigned)
                    _recurse(new_assigned + [v])

                return None

        _recurse()

        if self._best_vset:
            if structure.is_ordered:
                assigned = {}
                for val, sites in zip(self._best_vset, equi_sites):
                    for site in sites:
                        assigned[site] = val

                return [int(assigned[site]) for site in structure]
            assigned = {}
            new_best_vset = []
            for _ in equi_sites:
                new_best_vset.append([])
            for ival, val in enumerate(self._best_vset):
                new_best_vset[attrib[ival]].append(val)
            for val, sites in zip(new_best_vset, equi_sites):
                for site in sites:
                    assigned[site] = val

            return [[int(frac_site) for frac_site in assigned[site]] for site in structure]
        raise ValueError("Valences cannot be assigned!")

    def get_oxi_state_decorated_structure(self, structure):
        """
        Get an oxidation state decorated structure. This currently works only
        for ordered structures only.

        Args:
            structure: Structure to analyze

        Returns:
            A modified structure that is oxidation state decorated.

        Raises:
            ValueError if the valences cannot be determined.
        """
        s = structure.copy()
        if s.is_ordered:
            valences = self.get_valences(s)
            s.add_oxidation_state_by_site(valences)
        else:
            valences = self.get_valences(s)
            s = add_oxidation_state_by_site_fraction(s, valences)
        return s


def get_z_ordered_elmap(comp):
    """
    Arbitrary ordered elmap on the elements/species of a composition of a
    given site in an unordered structure. Returns a list of tuples (
    element_or_specie: occupation) in the arbitrary order.

    The arbitrary order is based on the Z of the element and the smallest
    fractional occupations first.
    Example : {"Ni3+": 0.2, "Ni4+": 0.2, "Cr3+": 0.15, "Zn2+": 0.34,
    "Cr4+": 0.11} will yield the species in the following order :
    Cr4+, Cr3+, Ni3+, Ni4+, Zn2+ ... or
    Cr4+, Cr3+, Ni4+, Ni3+, Zn2+
    """
    return sorted((elsp, comp[elsp]) for elsp in comp)


def add_oxidation_state_by_site_fraction(structure, oxidation_states):
    """
    Add oxidation states to a structure by fractional site.

    Args:
        oxidation_states (list): List of list of oxidation states for each
            site fraction for each site.
            E.g., [[2, 4], [3], [-2], [-2], [-2]]
    """
    try:
        for i, site in enumerate(structure):
            new_sp = collections.defaultdict(float)
            for j, (el, occu) in enumerate(get_z_ordered_elmap(site.species)):
                specie = Species(el.symbol, oxidation_states[i][j])
                new_sp[specie] += occu
            structure[i] = new_sp
        return structure
    except IndexError:
        raise ValueError("Oxidation state of all sites must be specified in the list.")
