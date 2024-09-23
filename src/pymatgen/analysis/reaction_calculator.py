"""This module provides classes that define a chemical reaction."""

from __future__ import annotations

import re
from itertools import chain, combinations
from typing import TYPE_CHECKING, no_type_check, overload

import numpy as np
from monty.fractions import gcd_float
from monty.json import MontyDecoder, MSONable
from uncertainties import ufloat

from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry

if TYPE_CHECKING:
    from collections.abc import Mapping

    from typing_extensions import Self

    from pymatgen.core import Element, Species
    from pymatgen.util.typing import CompositionLike

__author__ = "Shyue Ping Ong, Anubhav Jain"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Jul 11 2012"


class BalancedReaction(MSONable):
    """Represent a complete chemical reaction."""

    # Tolerance for determining if a particular component fraction is > 0.
    TOLERANCE = 1e-6

    @no_type_check
    def __init__(
        self,
        reactants_coeffs: Mapping[CompositionLike, int | float],
        products_coeffs: Mapping[CompositionLike, int | float],
    ) -> None:
        """
        Reactants and products to be specified as dict of {Composition: coeff}.

        Args:
            reactants_coeffs (dict[Composition, float]): Reactants as dict of {Composition: amt}.
            products_coeffs (dict[Composition, float]): Products as dict of {Composition: amt}.
        """
        # convert to Composition if necessary
        reactants_coeffs = {Composition(comp): coeff for comp, coeff in reactants_coeffs.items()}
        products_coeffs = {Composition(comp): coeff for comp, coeff in products_coeffs.items()}

        # sum reactants and products
        all_reactants = sum((comp * coeff for comp, coeff in reactants_coeffs.items()), Composition())

        all_products = sum((comp * coeff for comp, coeff in products_coeffs.items()), Composition())

        if not all_reactants.almost_equals(all_products, rtol=0, atol=self.TOLERANCE):
            raise ReactionError("Reaction is unbalanced!")

        self.reactants_coeffs: dict = reactants_coeffs
        self.products_coeffs: dict = products_coeffs

        # calculate net reaction coefficients
        self._coeffs: list[float] = []
        self._els: list[Element | Species] = []
        self._all_comp: list[Composition] = []
        for key in {*reactants_coeffs, *products_coeffs}:
            coeff = products_coeffs.get(key, 0) - reactants_coeffs.get(key, 0)

            if abs(coeff) > self.TOLERANCE:
                self._all_comp += [key]
                self._coeffs += [coeff]

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        for comp in self._all_comp:
            coeff2 = other.get_coeff(comp) if comp in other._all_comp else 0
            if abs(self.get_coeff(comp) - coeff2) > self.TOLERANCE:
                return False
        return True

    def __hash__(self) -> int:
        # Necessity for hash method is unclear (see gh-3673)
        return hash((frozenset(self.reactants_coeffs.items()), frozenset(self.products_coeffs.items())))

    def __str__(self):
        return self._str_from_comp(self._coeffs, self._all_comp)[0]

    __repr__ = __str__

    @overload
    def calculate_energy(self, energies: dict[Composition, ufloat]) -> ufloat:
        pass

    @overload
    def calculate_energy(self, energies: dict[Composition, float]) -> float:
        pass

    def calculate_energy(self, energies):
        """
        Calculates the energy of the reaction.

        Args:
            energies ({Composition: float}): Energy for each composition.
                E.g ., {comp1: energy1, comp2: energy2}.

        Returns:
            reaction energy as a float.
        """
        return sum(amt * energies[c] for amt, c in zip(self._coeffs, self._all_comp, strict=True))

    def normalize_to(self, comp: Composition, factor: float = 1) -> None:
        """
        Normalizes the reaction to one of the compositions.
        By default, normalizes such that the composition given has a
        coefficient of 1. Another factor can be specified.

        Args:
            comp (Composition): Composition to normalize to
            factor (float): Factor to normalize to. Defaults to 1.
        """
        scale_factor = abs(1 / self._coeffs[self._all_comp.index(comp)] * factor)
        self._coeffs = [c * scale_factor for c in self._coeffs]

    def normalize_to_element(self, element: Species | Element, factor: float = 1) -> None:
        """
        Normalizes the reaction to one of the elements.
        By default, normalizes such that the amount of the element is 1.
        Another factor can be specified.

        Args:
            element (SpeciesLike): Element to normalize to.
            factor (float): Factor to normalize to. Defaults to 1.
        """
        all_comp = self._all_comp
        coeffs = self._coeffs
        current_el_amount = sum(all_comp[i][element] * abs(coeffs[i]) for i in range(len(all_comp))) / 2
        scale_factor = factor / current_el_amount
        self._coeffs = [c * scale_factor for c in coeffs]

    def get_el_amount(self, element: Element | Species) -> float:
        """Get the amount of the element in the reaction.

        Args:
            element (SpeciesLike): Element in the reaction

        Returns:
            Amount of that element in the reaction.
        """
        return sum(self._all_comp[i][element] * abs(self._coeffs[i]) for i in range(len(self._all_comp))) / 2

    @property
    def elements(self) -> list[Element | Species]:
        """List of elements in the reaction."""
        return self._els

    @property
    def coeffs(self) -> list[float]:
        """Final coefficients of the calculated reaction."""
        return self._coeffs[:]

    @property
    def all_comp(self) -> list[Composition]:
        """List of all compositions in the reaction."""
        return self._all_comp

    @property
    def reactants(self) -> list[Composition]:
        """List of reactants."""
        return [self._all_comp[i] for i in range(len(self._all_comp)) if self._coeffs[i] < 0]

    @property
    def products(self) -> list[Composition]:
        """List of products."""
        return [self._all_comp[i] for i in range(len(self._all_comp)) if self._coeffs[i] > 0]

    def get_coeff(self, comp: Composition) -> float:
        """Get coefficient for a particular composition."""
        return self._coeffs[self._all_comp.index(comp)]

    def normalized_repr_and_factor(self) -> tuple[str, float]:
        """
        Normalized representation for a reaction
        For example, ``4 Li + 2 O -> 2Li2O`` becomes ``2 Li + O -> Li2O``.
        """
        return self._str_from_comp(self._coeffs, self._all_comp, reduce=True)

    @property
    def normalized_repr(self) -> str:
        """
        A normalized representation of the reaction. All factors are converted
        to lowest common factors.
        """
        return self.normalized_repr_and_factor()[0]

    @classmethod
    def _str_from_formulas(cls, coeffs, formulas) -> str:
        reactant_str = []
        product_str = []
        for amt, formula in zip(coeffs, formulas, strict=True):
            if abs(amt + 1) < cls.TOLERANCE:
                reactant_str.append(formula)
            elif abs(amt - 1) < cls.TOLERANCE:
                product_str.append(formula)
            elif amt < -cls.TOLERANCE:
                reactant_str.append(f"{-amt:.4g} {formula}")
            elif amt > cls.TOLERANCE:
                product_str.append(f"{amt:.4g} {formula}")

        return f"{' + '.join(reactant_str)} -> {' + '.join(product_str)}"

    @classmethod
    def _str_from_comp(cls, coeffs, compositions, reduce=False) -> tuple[str, float]:
        r_coeffs = np.zeros(len(coeffs))
        r_formulas = []
        for idx, (amt, comp) in enumerate(zip(coeffs, compositions, strict=True)):
            formula, factor = comp.get_reduced_formula_and_factor()
            r_coeffs[idx] = amt * factor
            r_formulas.append(formula)
        if reduce:
            factor = 1 / gcd_float(np.abs(r_coeffs))
            r_coeffs *= factor
        else:
            factor = 1
        return cls._str_from_formulas(r_coeffs, r_formulas), factor

    def as_entry(self, energies) -> ComputedEntry:
        """Get a ComputedEntry representation of the reaction."""
        relevant_comp = [comp * abs(coeff) for coeff, comp in zip(self._coeffs, self._all_comp, strict=True)]
        comp: Composition = sum(relevant_comp, Composition())  # type: ignore[assignment]

        entry = ComputedEntry(0.5 * comp, self.calculate_energy(energies))
        entry.name = str(self)
        return entry

    def as_dict(self) -> dict:
        """
        Returns:
            A dictionary representation of BalancedReaction.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "reactants": {str(comp): coeff for comp, coeff in self.reactants_coeffs.items()},
            "products": {str(comp): coeff for comp, coeff in self.products_coeffs.items()},
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): from as_dict().

        Returns:
            BalancedReaction
        """
        reactants = {Composition(comp): coeff for comp, coeff in dct["reactants"].items()}
        products = {Composition(comp): coeff for comp, coeff in dct["products"].items()}
        return cls(reactants, products)

    @classmethod
    def from_str(cls, rxn_str: str) -> Self:
        """Generate a balanced reaction from a string. The reaction must
        already be balanced.

        Args:
            rxn_string (str): The reaction string. For example, "4 Li + O2 -> 2Li2O"

        Returns:
            BalancedReaction
        """
        rct_str, prod_str = rxn_str.split("->")

        def get_comp_amt(comp_str):
            return {
                Composition(match[2]): float(match[1] or 1)
                for match in re.finditer(r"([\d\.]*(?:[eE]-?[\d\.]+)?)\s*([A-Z][\w\.\(\)]*)", comp_str)
            }

        return BalancedReaction(get_comp_amt(rct_str), get_comp_amt(prod_str))


class Reaction(BalancedReaction):
    """
    A more flexible class representing a Reaction. The reaction amounts will
    be automatically balanced. Reactants and products can swap sides so that
    all coefficients are positive, however this class will find the solution
    with the minimum number of swaps and coefficients of 0. Normalizes so that
    the *FIRST* product (or products, if underdetermined) has a coefficient of one.
    """

    def __init__(self, reactants: list[Composition], products: list[Composition]) -> None:
        """
        Reactants and products to be specified as list of
        pymatgen.core.structure.Composition. e.g. [comp1, comp2].

        Args:
            reactants ([Composition]): List of reactants.
            products ([Composition]): List of products.
        """
        self._input_reactants = reactants
        self._input_products = products
        self._all_comp = reactants + products
        self._num_comp = len(self.all_comp)

        all_elems = sorted({elem for c in self._all_comp for elem in c.elements})
        self._num_elems = len(all_elems)

        comp_matrix = np.array([[c[el] for el in all_elems] for c in self._all_comp]).T

        rank = np.linalg.matrix_rank(comp_matrix)
        diff = self._num_comp - rank
        num_constraints = diff if diff >= 2 else 1

        self._lowest_num_errors = np.inf  # an error = a component changing sides or disappearing

        self._coeffs = self._balance_coeffs(comp_matrix, num_constraints)
        self._els = all_elems

    def _balance_coeffs(self, comp_matrix, max_num_constraints):
        first_product_idx = len(self._input_reactants)

        # start with simplest product constraints, work towards most complex reactant constraints
        product_constraints = chain.from_iterable(
            [
                combinations(range(first_product_idx, self._num_comp), n_constr)
                for n_constr in range(max_num_constraints, 0, -1)
            ]
        )
        reactant_constraints = chain.from_iterable(
            [combinations(range(first_product_idx), n_constr) for n_constr in range(max_num_constraints, 0, -1)]
        )
        best_soln = None
        balanced = False

        for constraints in chain(product_constraints, reactant_constraints):
            n_constr = len(constraints)

            comp_and_constraints = np.append(comp_matrix, np.zeros((n_constr, self._num_comp)), axis=0)
            b = np.zeros((self._num_elems + n_constr, 1))
            b[-n_constr:] = 1 if min(constraints) >= first_product_idx else -1

            for num, idx in enumerate(constraints):
                comp_and_constraints[self._num_elems + num, idx] = 1
                # arbitrarily fix coeff to 1

            coeffs = np.matmul(np.linalg.pinv(comp_and_constraints), b)

            if np.allclose(np.matmul(comp_matrix, coeffs), np.zeros((self._num_elems, 1))):
                balanced = True
                expected_signs = np.array([-1] * len(self._input_reactants) + [+1] * len(self._input_products))
                num_errors = np.sum(np.multiply(expected_signs, coeffs.T) < self.TOLERANCE)

                if num_errors == 0:
                    self._lowest_num_errors = 0
                    return np.squeeze(coeffs)
                if num_errors < self._lowest_num_errors:
                    self._lowest_num_errors = num_errors
                    best_soln = coeffs

        if not balanced:
            raise ReactionError("Reaction cannot be balanced.")

        return np.squeeze(best_soln)

    def copy(self) -> Self:
        """Get a copy of the Reaction object."""
        return Reaction(self.reactants, self.products)

    def as_dict(self) -> dict:
        """
        Returns:
            A dictionary representation of Reaction.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "reactants": [comp.as_dict() for comp in self._input_reactants],
            "products": [comp.as_dict() for comp in self._input_products],
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): from as_dict().

        Returns:
            Reaction
        """
        reactants = [*map(Composition, dct["reactants"])]
        products = [*map(Composition, dct["products"])]
        return cls(reactants, products)


class ReactionError(Exception):
    """Exception class for Reactions. Allows more information in exception
    messages to cover situations not covered by standard exception classes.
    """

    def __init__(self, msg: str) -> None:
        """
        Create a ReactionError.

        Args:
            msg (str): More information about the ReactionError.
        """
        self.msg = msg

    def __str__(self) -> str:
        return self.msg


class ComputedReaction(Reaction):
    """
    Convenience class to generate a reaction from ComputedEntry objects, with
    some additional attributes, such as a reaction energy based on computed
    energies.
    """

    def __init__(self, reactant_entries: list[ComputedEntry], product_entries: list[ComputedEntry]) -> None:
        """
        Args:
            reactant_entries ([ComputedEntry]): List of reactant_entries.
            product_entries ([ComputedEntry]): List of product_entries.
        """
        self._reactant_entries = reactant_entries
        self._product_entries = product_entries
        self._all_entries = reactant_entries + product_entries
        reactant_comp = [entry.composition.reduced_composition for entry in reactant_entries]

        product_comp = [entry.composition.reduced_composition for entry in product_entries]

        super().__init__(list(reactant_comp), list(product_comp))

    @property
    def all_entries(self):
        """Equivalent of all_comp but returns entries, in the same order as the
        coefficients.
        """
        entries = []
        for comp in self._all_comp:
            for entry in self._all_entries:
                if entry.reduced_formula == comp.reduced_formula:
                    entries.append(entry)
                    break
        return entries

    @property
    def calculated_reaction_energy(self) -> float:
        """
        Returns:
            float: The calculated reaction energy.
        """
        calc_energies: dict[Composition, float] = {}

        for entry in self._reactant_entries + self._product_entries:
            comp, factor = entry.composition.get_reduced_composition_and_factor()
            calc_energies[comp] = min(calc_energies.get(comp, float("inf")), entry.energy / factor)

        return self.calculate_energy(calc_energies)

    @property
    def calculated_reaction_energy_uncertainty(self) -> float:
        """
        Calculates the uncertainty in the reaction energy based on the uncertainty in the
        energies of the products and reactants.
        """
        calc_energies: dict[Composition, float] = {}

        for entry in self._reactant_entries + self._product_entries:
            comp, factor = entry.composition.get_reduced_composition_and_factor()
            energy_ufloat = ufloat(entry.energy, entry.correction_uncertainty)
            calc_energies[comp] = min(calc_energies.get(comp, float("inf")), energy_ufloat / factor)

        return self.calculate_energy(calc_energies).std_dev

    def as_dict(self) -> dict:
        """
        Returns:
            A dictionary representation of ComputedReaction.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "reactants": [entry.as_dict() for entry in self._reactant_entries],
            "products": [entry.as_dict() for entry in self._product_entries],
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): from as_dict().

        Returns:
            A ComputedReaction object.
        """
        reactants = [MontyDecoder().process_decoded(entry) for entry in dct["reactants"]]
        products = [MontyDecoder().process_decoded(entry) for entry in dct["products"]]
        return cls(reactants, products)
