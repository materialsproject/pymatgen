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

from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure_modifier import StructureEditor


class BondValenceParam(object):
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
        self.species = frozenset(species)
        self.syms = set([sp.symbol for sp in species])
        self.R = R
        self.b = b

    def contains_element(self, el):
        return el.symbol in self.syms

    def get_sign(self, el):
        for sp in self.species:
            if sp.symbol == el.symbol:
                return sp.oxi_state / abs(sp.oxi_state)
        raise ValueError("{} not in BV {}".format(el, self.__repr__()))

    def contains_element_pair(self, el1, el2):
        for sp1, sp2 in itertools.permutations(self.species, 2):
            if sp1.symbol == el1.symbol and sp2.symbol == el2.symbol:
                return True
        return False

    def __repr__(self):
        return "{}: R = {}, b = {}".format(self.species, self.R, self.b)

_BV_PARAM = []
_species_set = []
module_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(module_dir, 'bvparm2011'), 'r') as f:
    for l in f:
        toks = l.split()
        sp1 = Specie(toks[0], int(toks[1]))
        sp2 = Specie(toks[2], int(toks[3]))
        species = frozenset([sp1, sp2])
        if species not in _species_set:
            _BV_PARAM.append(BondValenceParam([sp1, sp2], float(toks[4]),
                                             float(toks[5])))
            _species_set.append(species)


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
        assigned_valences = []
        if self.symm_tol:
            finder = SymmetryFinder(structure, self.symm_tol)
            symm_structure = finder.get_symmetrized_structure()
            equi_sites = symm_structure.equivalent_sites
        else:
            symm_structure = structure
            equi_sites = [[site] for site in structure]
        to_test = []
        bond_types = {}
        for sites in equi_sites:
            data = {}
            data["nsites"] = len(sites)
            test_site = sites[0]
            data["test_site"] = test_site
            el1 = Element(test_site.specie.symbol)
            test_nn = []
            for nn, dist in symm_structure.get_neighbors(test_site,
                                                         self.max_radius):
                el2 = Element(nn.specie.symbol)
                bvs = [bv for bv in relevant_bvs
                       if bv.contains_element_pair(el1, el2)]
                if bvs:
                    test_nn.append((nn, dist))
                    bond_types[frozenset([el1, el2])] = bvs
            data["nn"] = test_nn
            to_test.append(data)

        found_val = None
        for bv_set in itertools.product(*bond_types.values()):
            assigned_valences = []
            total_valence = 0
            for data in to_test:
                el1 = Element(data["test_site"].specie.symbol)
                bvs = 0
                for nn, dist in data["nn"]:
                    nn_el = Element(nn.specie.symbol)
                    for bv in bv_set:
                        if bv.contains_element_pair(el1, nn_el):
                            bvs += math.exp((bv.R - dist) / bv.b) * \
                                bv.get_sign(el1)
                assigned_valences.append(bvs)
                total_valence += int(round(bvs)) * data['nsites']
            if total_valence == 0:
                found_val = [int(round(i)) for i in assigned_valences]
                break

        if found_val:
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
