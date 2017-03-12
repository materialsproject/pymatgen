# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import logging
import numpy as np
import re

from monty.json import MSONable
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from monty.json import MontyDecoder
from monty.fractions import gcd_float

import itertools

"""
This module provides classes that define a chemical reaction.
"""


__author__ = "Shyue Ping Ong, Anubhav Jain"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Jul 11 2012"


logger = logging.getLogger(__name__)


class BalancedReaction(MSONable):
    """
    An object representing a complete chemical reaction.
    """

    # Tolerance for determining if a particular component fraction is > 0.
    TOLERANCE = 1e-6

    def __init__(self, reactants_coeffs, products_coeffs):
        """
        Reactants and products to be specified as dict of {Composition: coeff}.

        Args:
            reactants_coeffs ({Composition: float}): Reactants as dict of
                {Composition: amt}.
            products_coeffs ({Composition: float}): Products as dict of
                {Composition: amt}.
        """
        # sum reactants and products
        all_reactants = sum([k * v for k, v in reactants_coeffs.items()],
                            Composition({}))
        all_products = sum([k * v for k, v in products_coeffs.items()],
                           Composition({}))

        if not all_reactants.almost_equals(all_products, rtol=0,
                                           atol=self.TOLERANCE):
            raise ReactionError("Reaction is unbalanced!")

        self._els = all_reactants.elements

        self.reactants_coeffs = reactants_coeffs
        self.products_coeffs = products_coeffs

        # calculate net reaction coefficients
        self._coeffs = []
        self._els = []
        self._all_comp = []
        for c in set(list(reactants_coeffs.keys()) +
                     list(products_coeffs.keys())):
            coeff = products_coeffs.get(c, 0) - reactants_coeffs.get(c, 0)

            if abs(coeff) > self.TOLERANCE:
                self._all_comp.append(c)
                self._coeffs.append(coeff)

    def calculate_energy(self, energies):
        """
        Calculates the energy of the reaction.

        Args:
            energies ({Composition: float}): Energy for each composition.
                E.g ., {comp1: energy1, comp2: energy2}.

        Returns:
            reaction energy as a float.
        """
        return sum([amt * energies[c] for amt, c in zip(self._coeffs,
                                                        self._all_comp)])

    def normalize_to(self, comp, factor=1):
        """
        Normalizes the reaction to one of the compositions.
        By default, normalizes such that the composition given has a
        coefficient of 1. Another factor can be specified.

        Args:
            comp (Composition): Composition to normalize to
            factor (float): Factor to normalize to. Defaults to 1.
        """
        scale_factor = abs(1 / self._coeffs[self._all_comp.index(comp)]
                           * factor)
        self._coeffs = [c * scale_factor for c in self._coeffs]

    def normalize_to_element(self, element, factor=1):
        """
        Normalizes the reaction to one of the elements.
        By default, normalizes such that the amount of the element is 1.
        Another factor can be specified.

        Args:
            element (Element/Specie): Element to normalize to.
            factor (float): Factor to normalize to. Defaults to 1.
        """
        all_comp = self._all_comp
        coeffs = self._coeffs
        current_el_amount = sum([all_comp[i][element] * abs(coeffs[i])
                                 for i in range(len(all_comp))]) / 2
        scale_factor = factor / current_el_amount
        self._coeffs = [c * scale_factor for c in coeffs]

    def get_el_amount(self, element):
        """
        Returns the amount of the element in the reaction.

        Args:
            element (Element/Specie): Element in the reaction

        Returns:
            Amount of that element in the reaction.
        """
        return sum([self._all_comp[i][element] * abs(self._coeffs[i])
                    for i in range(len(self._all_comp))]) / 2

    @property
    def elements(self):
        """
        List of elements in the reaction
        """
        return self._els[:]

    @property
    def coeffs(self):
        """
        Final coefficients of the calculated reaction
        """
        return self._coeffs[:]

    @property
    def all_comp(self):
        """
        List of all compositions in the reaction.
        """
        return self._all_comp

    @property
    def reactants(self):
        """
        List of reactants
        """
        return [self._all_comp[i] for i in range(len(self._all_comp))
                if self._coeffs[i] < 0]

    @property
    def products(self):
        """
        List of products
        """
        return [self._all_comp[i] for i in range(len(self._all_comp))
                if self._coeffs[i] > 0]

    def get_coeff(self, comp):
        """
        Returns coefficient for a particular composition
        """
        return self._coeffs[self._all_comp.index(comp)]

    def normalized_repr_and_factor(self):
        """
        Normalized representation for a reaction
        For example, ``4 Li + 2 O -> 2Li2O`` becomes ``2 Li + O -> Li2O``
        """
        return self._str_from_comp(self._coeffs, self._all_comp, True)

    @property
    def normalized_repr(self):
        """
        A normalized representation of the reaction. All factors are converted
        to lowest common factors.
        """
        return self.normalized_repr_and_factor()[0]

    def __eq__(self, other):
        if other is None:
            return False
        for comp in self._all_comp:
            coeff2 = other.get_coeff(comp) if comp in other._all_comp else 0
            if abs(self.get_coeff(comp) - coeff2) > self.TOLERANCE:
                return False
        return True

    def __hash__(self):
        return 7

    @classmethod
    def _str_from_formulas(cls, coeffs, formulas):
        reactant_str = []
        product_str = []
        for amt, formula in zip(coeffs, formulas):
            if abs(amt + 1) < cls.TOLERANCE:
                reactant_str.append(formula)
            elif abs(amt - 1) < cls.TOLERANCE:
                product_str.append(formula)
            elif amt < -cls.TOLERANCE:
                reactant_str.append("{:.4g} {}".format(-amt, formula))
            elif amt > cls.TOLERANCE:
                product_str.append("{:.4g} {}".format(amt, formula))

        return " + ".join(reactant_str) + " -> " + " + ".join(product_str)

    @classmethod
    def _str_from_comp(cls, coeffs, compositions, reduce=False):
        r_coeffs = np.zeros(len(coeffs))
        r_formulas = []
        for i, (amt, comp) in enumerate(zip(coeffs, compositions)):
            formula, factor = comp.get_reduced_formula_and_factor()
            r_coeffs[i] = amt * factor
            r_formulas.append(formula)
        if reduce:
            factor = 1 / gcd_float(np.abs(r_coeffs))
            r_coeffs *= factor
        else:
            factor = 1
        return cls._str_from_formulas(r_coeffs, r_formulas), factor

    def __str__(self):
        return self._str_from_comp(self._coeffs, self._all_comp)[0]

    __repr__ = __str__

    def as_entry(self, energies):
        """
        Returns a ComputedEntry representation of the reaction.
        :return:
        """
        relevant_comp = [comp * abs(coeff) for coeff, comp
                         in zip(self._coeffs, self._all_comp)]
        comp = sum(relevant_comp, Composition())
        entry = ComputedEntry(0.5 * comp, self.calculate_energy(energies))
        entry.name = self.__str__()
        return entry

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "reactants": {str(comp): coeff
                              for comp, coeff in self.reactants_coeffs.items()},
                "products": {str(comp): coeff
                             for comp, coeff in self.products_coeffs.items()}}

    @classmethod
    def from_dict(cls, d):
        reactants = {Composition(comp): coeff
                     for comp, coeff in d["reactants"].items()}
        products = {Composition(comp): coeff
                    for comp, coeff in d["products"].items()}
        return cls(reactants, products)

    @staticmethod
    def from_string(rxn_string):
        """
        Generates a balanced reaction from a string. The reaction must
        already be balanced.

        Args:
            rxn_string:
                The reaction string. For example, "4 Li + O2-> 2Li2O"

        Returns:
            BalancedReaction
        """
        rct_str, prod_str = rxn_string.split("->")

        def get_comp_amt(comp_str):
            return {Composition(m.group(2)): float(m.group(1) or 1)
                    for m in re.finditer(r"([\d\.]*(?:[eE]-?[\d\.]+)?)\s*([A-Z][\w\.\(\)]*)",
                                         comp_str)}
        return BalancedReaction(get_comp_amt(rct_str), get_comp_amt(prod_str))


class Reaction(BalancedReaction):
    """
    A more flexible class representing a Reaction. The reaction amounts will
    be automatically balanced. Reactants and products can swap sides so that
    all coefficients are positive. Normalizes so that the *FIRST* product
    (or products, if underdetermined) has a coefficient of one.
    """

    def __init__(self, reactants, products):
        """
        Reactants and products to be specified as list of
        pymatgen.core.structure.Composition.  e.g., [comp1, comp2]

        Args:
            reactants ([Composition]): List of reactants.
            products ([Composition]): List of products.
        """
        self._input_reactants = reactants
        self._input_products = products
        self._all_comp = reactants + products

        els = set()
        for c in self.all_comp:
            els.update(c.elements)
        els = sorted(els)

        # Solving:
        #          | 0  R |
        # [ x y ]  |      |  =  [ 1 .. 1 0 .. 0]
        #          | C  P |
        # x, y are the coefficients of the reactants and products
        # R, P the matrices of the element compositions of the reactants
        # and products
        # C is a constraint matrix that chooses which compositions to normalize to

        # try just normalizing to just the first product
        rp_mat = np.array([[c[el] for el in els] for c in self._all_comp])
        f_mat = np.concatenate([np.zeros((len(rp_mat), 1)), rp_mat], axis=1)
        f_mat[len(reactants), 0] = 1  # set normalization by the first product
        b = np.zeros(len(els) + 1)
        b[0] = 1
        coeffs, res, _, s = np.linalg.lstsq(f_mat.T, b)

        # for whatever reason the rank returned by lstsq isn't always correct
        # seems to be a problem with low-rank M but inconsistent system  M x = b.
        # the singular values seem ok, so checking based on those
        if sum(np.abs(s) > 1e-12) == len(f_mat):
            if res and res[0] > self.TOLERANCE ** 2:
                raise ReactionError("Reaction cannot be balanced.")
            else:
                ok = True
        else:
            # underdetermined, add product constraints to make non-singular
            ok = False
            n_constr = len(rp_mat) - np.linalg.matrix_rank(rp_mat)
            f_mat = np.concatenate([np.zeros((len(rp_mat), n_constr)),
                                    rp_mat], axis=1)
            b = np.zeros(f_mat.shape[1])
            b[:n_constr] = 1

            # try setting C to all n_constr combinations of products
            for inds in itertools.combinations(range(len(reactants),
                                                     len(f_mat)),
                                               n_constr):
                f_mat[:, :n_constr] = 0
                for j, i in enumerate(inds):
                    f_mat[i, j] = 1
                # try a solution
                coeffs, res, _, s = np.linalg.lstsq(f_mat.T, b)
                if sum(np.abs(s) > 1e-12) == len(self._all_comp) and \
                        (not res or res[0] < self.TOLERANCE ** 2):
                    ok = True
                    break

        if not ok:
            r_mat = np.array([[c[el] for el in els] for c in reactants])
            if np.linalg.lstsq(r_mat.T, np.zeros(len(els)))[2] != len(reactants):
                raise ReactionError("Reaction cannot be balanced. "
                                    "Reactants are underdetermined.")
            raise ReactionError("Reaction cannot be balanced. "
                                "Unknown error, please report.")

        self._els = els
        self._coeffs = coeffs

    def copy(self):
        """
        Returns a copy of the Reaction object.
        """
        return Reaction(self.reactants, self.products)

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "reactants": [comp.as_dict() for comp in self._input_reactants],
                "products": [comp.as_dict() for comp in self._input_products]}

    @classmethod
    def from_dict(cls, d):
        reactants = [Composition(sym_amt) for sym_amt in d["reactants"]]
        products = [Composition(sym_amt) for sym_amt in d["products"]]
        return cls(reactants, products)

class ReactionError(Exception):
    """
    Exception class for Reactions. Allows more information in exception
    messages to cover situations not covered by standard exception classes.
    """

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class ComputedReaction(Reaction):
    """
    Convenience class to generate a reaction from ComputedEntry objects, with
    some additional attributes, such as a reaction energy based on computed
    energies.
    """

    def __init__(self, reactant_entries, product_entries):
        """
        Args:
            reactant_entries ([ComputedEntry]): List of reactant_entries.
            product_entries ([ComputedEntry]): List of product_entries.
        """
        self._reactant_entries = reactant_entries
        self._product_entries = product_entries
        self._all_entries = reactant_entries + product_entries
        reactant_comp = set([e.composition
                             .get_reduced_composition_and_factor()[0]
                             for e in reactant_entries])
        product_comp = set([e.composition
                            .get_reduced_composition_and_factor()[0]
                            for e in product_entries])
        super(ComputedReaction, self).__init__(list(reactant_comp),
                                               list(product_comp))

    @property
    def all_entries(self):
        """
        Equivalent of all_comp but returns entries, in the same order as the
        coefficients.
        """
        entries = []
        for c in self._all_comp:
            for e in self._all_entries:
                if e.composition.reduced_formula == c.reduced_formula:
                    entries.append(e)
                    break
        return entries

    @property
    def calculated_reaction_energy(self):
        calc_energies = {}

        for entry in self._reactant_entries + self._product_entries:
            (comp, factor) = \
                entry.composition.get_reduced_composition_and_factor()
            calc_energies[comp] = min(calc_energies.get(comp, float('inf')),
                                      entry.energy / factor)

        return self.calculate_energy(calc_energies)

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "reactants": [e.as_dict() for e in self._reactant_entries],
                "products": [e.as_dict() for e in self._product_entries]}

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        reactants = [dec.process_decoded(e) for e in d["reactants"]]
        products = [dec.process_decoded(e) for e in d["products"]]
        return cls(reactants, products)
