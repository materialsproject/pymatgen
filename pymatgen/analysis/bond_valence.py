#!/usr/bin/env python

"""
Bond valence analysis module. Bond valence is calculated as:
vi = exp((R0 - Ri) / b)
and
V = sum(vi)
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Oct 24, 2012"

import itertools
import math
import os
import json
import collections

from pymatgen.serializers.json_coders import MSONable, PMGJSONDecoder
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure_modifier import StructureEditor


class BondValenceParam(MSONable):
    """
    A helper object to contain a particular set of bond valence parameters.
    """
    def __init__(self, species, R, b):
        """
        Args:
            species:
                A sequence of two species, e.g., Li+, O2- for the bond.
            R:
                The R parameter.
            b:
                The b parameter.
        """
        self.species = tuple(species)
        self.syms = set([sp.symbol for sp in species])
        self.R = R
        self.b = b

    def contains_element(self, el):
        """
        Whether parameter contains a particular element.

        Args:
            el:
                Element

        Returns:
            True if element is in parameter.
        """
        return el.symbol in self.syms

    def get_sign(self, el):
        """
        Returns sign of valence for element.

        Args:
            el:
                Element

        Returns:
            Sign for element.
        """
        for sp in self.species:
            if sp.symbol == el.symbol:
                return sp.oxi_state / abs(sp.oxi_state)
        raise ValueError("{} not in BV {}".format(el, self.__repr__()))

    def get_valence(self, el):
        """
        Valence for a particular element.

        Args:
            el:
                Element

        Returns:
            Valence of the element in the parameter.

        Raises:
            Value error if el is not found in parameter.
        """
        for sp in self.species:
            if sp.symbol == el.symbol:
                return sp.oxi_state
        raise ValueError("{} not in BV {}".format(el, self.__repr__()))

    def contains_element_pair(self, el1, el2):
        """
        Whether a particular element pair is in parameter.

        Args:
            el1:
                Element 1
            el2:
                Element 2

        Returns:
            True if elements are in the parameter.
        """
        for sp1, sp2 in itertools.permutations(self.species, 2):
            if sp1.symbol == el1.symbol and sp2.symbol == el2.symbol:
                return True
        return False

    def contains_species_pair(self, sp1, sp2):
        """
        Whether a particular species pair is in parameter.

        Args:
            sp1:
                Specie 1
            sp2:
                Specie 2

        Returns:
            True if species are in the parameter.
        """
        for s1, s2 in itertools.permutations(self.species, 2):
            if s1 == sp1 and s2 == sp2:
                return True
        return False

    def __repr__(self):
        return "{}: R = {}, b = {}".format(self.species, self.R, self.b)

    @property
    def to_dict(self):
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["species"] = [sp.to_dict for sp in self.species]
        d["R"] = self.R
        d["b"] = self.b
        return d

    @staticmethod
    def from_dict(d):
        species = [Specie.from_dict(sp_dict) for sp_dict in d["species"]]
        return BondValenceParam(species, d["R"], d["b"])


module_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(module_dir, "bvparm2011.json"), "r") as f:
    _BV_PARAM = json.load(f, cls=PMGJSONDecoder)


class BVAnalyzer(object):
    """
    Bond valence analyzer for assigning oxidation states to structures.
    """

    def __init__(self, symm_tol=0.1, max_radius=4):
        """
        Args:
            symm_tol:
                Symmetry tolerance used to determine which sites are
                symmetrically equivalent. Set to 0 to turn off symmetry.
            max_radius:
                Maximum radius in Angstrom used to find nearest neighbors.
        """
        self.symm_tol = symm_tol
        self.max_radius = max_radius

    def get_valences(self, structure):
        """
        Returns a list of valences for the structure.

        Args:
            structure:
                Structure to analyze

        Returns:
            A list of valences for each site in the structure,
            e.g., [1, 1, -2].

        Raises:
            A ValueError is the valences cannot be determined.
        """
        syms = [el.symbol for el in structure.composition.elements]
        relevant_bvs = [bv for bv in _BV_PARAM
                        if all([sp.symbol in syms for sp in bv.species])]
        if self.symm_tol:
            finder = SymmetryFinder(structure, self.symm_tol)
            symm_structure = finder.get_symmetrized_structure()
            equi_sites = symm_structure.equivalent_sites
        else:
            symm_structure = structure
            equi_sites = [[site] for site in structure]

        equi_sites = sorted(equi_sites, key=lambda sites:-sites[0].specie.X)

        all_nn = structure.get_all_neighbors(self.max_radius,
                                             include_index=True)

        el1 = Element(equi_sites[0][0].specie.symbol)
        filtered_bvs = []
        for bv in relevant_bvs:
            if (not bv.contains_element(el1)) or bv.get_sign(el1) < 0:
                filtered_bvs.append(bv)
        relevant_bvs = filtered_bvs

        def get_equi_index(site):
            for j, sites in enumerate(equi_sites):
                if site in sites:
                    return j
            raise ValueError("Can't find equi-index!")

        el_map = []
        bonds = collections.defaultdict(int)
        bond_dist = {}
        for i, sites in enumerate(equi_sites):
            test_site = sites[0]
            ind = structure.sites.index(test_site)
            el1 = Element(test_site.specie.symbol)
            el_map.append(el1)
            for (nn, dist, nn_index) in all_nn[ind]:
                equi_index = get_equi_index(structure[nn_index])
                bond_index = tuple(sorted([i, equi_index]))
                bonds[bond_index] += 1
                bond_dist[bond_index] = dist

        valence_list = []
        for el in el_map:
            valences = []
            for bv in relevant_bvs:
                if bv.contains_element(el):
                    valences.append(bv.get_valence(el))
            valence_list.append(set(valences))

        all_results = []
        for valences in itertools.product(*valence_list):
            sum_val = sum([k * len(v) for k, v in zip(valences, equi_sites)])
            if sum_val == 0:
                total_bvs = collections.defaultdict(float)
                for bond, num in bonds.items():
                    ind1 = bond[0]
                    ind2 = bond[1]
                    sp1 = Specie(el_map[ind1].symbol, valences[ind1])
                    sp2 = Specie(el_map[ind2].symbol, valences[ind2])
                    for bv in relevant_bvs:
                        if bv.contains_species_pair(sp1, sp2):
                            val = math.exp((bv.R - bond_dist[bond]) / bv.b)
                            total_bvs[ind1] += val * num * \
                                bv.get_sign(el_map[ind1])
                            total_bvs[ind2] += val * num * \
                                bv.get_sign(el_map[ind2])
                            break
                error = 0
                all_avg = []
                for i, sum_bv in total_bvs.items():
                    avg_bv = sum_bv / len(equi_sites[i])
                    error += math.pow(avg_bv - valences[i], 2)
                    all_avg.append(avg_bv)
                all_results.append((valences, error, all_avg))
        if all_results:
            all_results = sorted(all_results, key=lambda x: x[1])
            found_val = all_results[0][0]
            oxi_states = []
            for site in structure:
                for i, sites in enumerate(equi_sites):
                    if site in sites:
                        oxi_states.append(found_val[i])
                        break
            return oxi_states
        else:
            raise ValueError("Valences cannot be assigned!")

    def get_oxi_state_decorated_structure(self, structure):
        """
        Get an oxidation state decorated structure.

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
