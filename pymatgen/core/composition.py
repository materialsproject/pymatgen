# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements a Composition class to represent compositions,
and a ChemicalPotential class to represent potentials.
"""

import collections
import numbers
import string
from itertools import combinations_with_replacement, product
import os
import re
from typing import Tuple, List
from functools import total_ordering

from monty.serialization import loadfn
from monty.fractions import gcd, gcd_float
from monty.json import MSONable

from pymatgen.core.periodic_table import get_el_sp, Element, Specie, DummySpecie
from pymatgen.util.string import formula_double_format
from pymatgen.core.units import Mass

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Nov 10, 2012"


@total_ordering
class Composition(collections.abc.Hashable, collections.abc.Mapping, MSONable):
    """
    Represents a Composition, which is essentially a {element:amount} mapping
    type. Composition is written to be immutable and hashable,
    unlike a standard Python dict.

    Note that the key can be either an Element or a Specie. Elements and Specie
    are treated differently. i.e., a Fe2+ is not the same as a Fe3+ Specie and
    would be put in separate keys. This differentiation is deliberate to
    support using Composition to determine the fraction of a particular Specie.

    Works almost completely like a standard python dictionary, except that
    __getitem__ is overridden to return 0 when an element is not found.
    (somewhat like a defaultdict, except it is immutable).

    Also adds more convenience methods relevant to compositions, e.g.,
    get_fraction.

    It should also be noted that many Composition related functionality takes
    in a standard string as a convenient input. For example,
    even though the internal representation of a Fe2O3 composition is
    {Element("Fe"): 2, Element("O"): 3}, you can obtain the amount of Fe
    simply by comp["Fe"] instead of the more verbose comp[Element("Fe")].

    >>> comp = Composition("LiFePO4")
    >>> comp.get_atomic_fraction(Element("Li"))
    0.14285714285714285
    >>> comp.num_atoms
    7.0
    >>> comp.reduced_formula
    'LiFePO4'
    >>> comp.formula
    'Li1 Fe1 P1 O4'
    >>> comp.get_wt_fraction(Element("Li"))
    0.04399794666951898
    >>> comp.num_atoms
    7.0
    """

    # Tolerance in distinguishing different composition amounts.
    # 1e-8 is fairly tight, but should cut out most floating point arithmetic
    # errors.
    amount_tolerance = 1e-8

    # Special formula handling for peroxides and certain elements. This is so
    # that formula output does not write LiO instead of Li2O2 for example.
    special_formulas = {"LiO": "Li2O2", "NaO": "Na2O2", "KO": "K2O2",
                        "HO": "H2O2", "CsO": "Cs2O2", "RbO": "Rb2O2",
                        "O": "O2", "N": "N2", "F": "F2", "Cl": "Cl2",
                        "H": "H2"}

    oxi_prob = None  # prior probability of oxidation used by oxi_state_guesses

    def __init__(self, *args, strict=False, **kwargs):  # allow_negative=False
        r"""
        Very flexible Composition construction, similar to the built-in Python
        dict(). Also extended to allow simple string init.

        Args:
            Any form supported by the Python built-in dict() function.

            1. A dict of either {Element/Specie: amount},

               {string symbol:amount}, or {atomic number:amount} or any mixture
               of these. E.g., {Element("Li"):2 ,Element("O"):1},
               {"Li":2, "O":1}, {3:2, 8:1} all result in a Li2O composition.
            2. Keyword arg initialization, similar to a dict, e.g.,

               Composition(Li = 2, O = 1)

            In addition, the Composition constructor also allows a single
            string as an input formula. E.g., Composition("Li2O").

            strict: Only allow valid Elements and Species in the Composition.

            allow_negative: Whether to allow negative compositions. This
                argument must be popped from the **kwargs due to *args
                ambiguity.
        """
        self.allow_negative = kwargs.pop('allow_negative', False)
        # it's much faster to recognize a composition and use the elmap than
        # to pass the composition to dict()
        if len(args) == 1 and isinstance(args[0], Composition):
            elmap = args[0]
        elif len(args) == 1 and isinstance(args[0], str):
            elmap = self._parse_formula(args[0])
        else:
            elmap = dict(*args, **kwargs)
        elamt = {}
        self._natoms = 0
        for k, v in elmap.items():
            if v < -Composition.amount_tolerance and not self.allow_negative:
                raise CompositionError("Amounts in Composition cannot be "
                                       "negative!")
            if abs(v) >= Composition.amount_tolerance:
                elamt[get_el_sp(k)] = v
                self._natoms += abs(v)
        self._data = elamt
        if strict and not self.valid:
            raise ValueError("Composition is not valid, contains: {}"
                             .format(", ".join(map(str, self.elements))))

    def __getitem__(self, item):
        try:
            sp = get_el_sp(item)
            return self._data.get(sp, 0)
        except ValueError as ex:
            raise TypeError("Invalid key {}, {} for Composition\n"
                            "ValueError exception:\n{}".format(item,
                                                               type(item), ex))

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return self._data.keys().__iter__()

    def __contains__(self, item):
        try:
            sp = get_el_sp(item)
            return sp in self._data
        except ValueError as ex:
            raise TypeError("Invalid key {}, {} for Composition\n"
                            "ValueError exception:\n{}".format(item,
                                                               type(item), ex))

    def __eq__(self, other):
        #  elements with amounts < Composition.amount_tolerance don't show up
        #  in the elmap, so checking len enables us to only check one
        #  compositions elements
        if len(self) != len(other):
            return False
        for el, v in self.items():
            if abs(v - other[el]) > Composition.amount_tolerance:
                return False
        return True

    def __ge__(self, other):
        """
        Defines >= for Compositions. Should ONLY be used for defining a sort
        order (the behavior is probably not what you'd expect)
        """
        for el in sorted(set(self.elements + other.elements)):
            if other[el] - self[el] >= Composition.amount_tolerance:
                return False
            if self[el] - other[el] >= Composition.amount_tolerance:
                return True
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        """
        Adds two compositions. For example, an Fe2O3 composition + an FeO
        composition gives a Fe3O4 composition.
        """
        new_el_map = collections.defaultdict(float)
        new_el_map.update(self)
        for k, v in other.items():
            new_el_map[get_el_sp(k)] += v
        return Composition(new_el_map, allow_negative=self.allow_negative)

    def __sub__(self, other):
        """
        Subtracts two compositions. For example, an Fe2O3 composition - an FeO
        composition gives an FeO2 composition.

        Raises:
            CompositionError if the subtracted composition is greater than the
            original composition in any of its elements, unless allow_negative
            is True
        """
        new_el_map = collections.defaultdict(float)
        new_el_map.update(self)
        for k, v in other.items():
            new_el_map[get_el_sp(k)] -= v
        return Composition(new_el_map, allow_negative=self.allow_negative)

    def __mul__(self, other):
        """
        Multiply a Composition by an integer or a float.
        Fe2O3 * 4 -> Fe8O12
        """
        if not isinstance(other, numbers.Number):
            return NotImplemented
        return Composition({el: self[el] * other for el in self},
                           allow_negative=self.allow_negative)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, numbers.Number):
            return NotImplemented
        return Composition({el: self[el] / other for el in self},
                           allow_negative=self.allow_negative)

    __div__ = __truediv__

    def __hash__(self):
        """
        Minimally effective hash function that just distinguishes between
        Compositions with different elements.
        """
        hashcode = 0
        for el, amt in self.items():
            if abs(amt) > Composition.amount_tolerance:
                hashcode += el.Z
        return hashcode

    @property
    def average_electroneg(self) -> float:
        """
        :return: Average electronegativity of the composition.
        """
        return sum((el.X * abs(amt) for el, amt in self.items())) / self.num_atoms

    @property
    def total_electrons(self) -> float:
        """
        :return: Total number of electrons in composition.
        """
        return sum((el.Z * abs(amt) for el, amt in self.items()))

    def almost_equals(self, other, rtol=0.1, atol=1e-8):
        """
        Returns true if compositions are equal within a tolerance.

        Args:
            other (Composition): Other composition to check
            rtol (float): Relative tolerance
            atol (float): Absolute tolerance
        """
        sps = set(self.elements + other.elements)
        for sp in sps:
            a = self[sp]
            b = other[sp]
            tol = atol + rtol * (abs(a) + abs(b)) / 2
            if abs(b - a) > tol:
                return False
        return True

    @property
    def is_element(self) -> bool:
        """
        True if composition is for an element.
        """
        return len(self) == 1

    def copy(self):
        """
        :return: A copy of the composition.
        """
        return Composition(self, allow_negative=self.allow_negative)

    @property
    def formula(self) -> str:
        """
        Returns a formula string, with elements sorted by electronegativity,
        e.g., Li4 Fe4 P4 O16.
        """
        sym_amt = self.get_el_amt_dict()
        syms = sorted(sym_amt.keys(), key=lambda sym: get_el_sp(sym).X)
        formula = [s + formula_double_format(sym_amt[s], False) for s in syms]
        return " ".join(formula)

    @property
    def alphabetical_formula(self) -> str:
        """
        Returns a formula string, with elements sorted by alphabetically
        e.g., Fe4 Li4 O16 P4.
        """
        sym_amt = self.get_el_amt_dict()
        syms = sorted(sym_amt.keys())
        formula = [s + formula_double_format(sym_amt[s], False) for s in syms]
        return " ".join(formula)

    @property
    def iupac_formula(self) -> str:
        """
        Returns a formula string, with elements sorted by the iupac
        electronegativity ordering defined in Table VI of "Nomenclature of
        Inorganic Chemistry (IUPAC Recommendations 2005)". This ordering
        effectively follows the groups and rows of the periodic table, except
        the Lanthanides, Actanides and hydrogen. Polyanions are still determined
        based on the true electronegativity of the elements.
        e.g. CH2(SO4)2
        """
        sym_amt = self.get_el_amt_dict()
        syms = sorted(sym_amt.keys(),
                      key=lambda s: get_el_sp(s).iupac_ordering)
        formula = [s + formula_double_format(sym_amt[s], False) for s in syms]
        return " ".join(formula)

    @property
    def element_composition(self) -> 'Composition':
        """
        Returns the composition replacing any species by the corresponding
        element.
        """
        return Composition(self.get_el_amt_dict(),
                           allow_negative=self.allow_negative)

    @property
    def fractional_composition(self) -> 'Composition':
        """
        Returns the normalized composition which the number of species sum to
        1.

        Returns:
            Normalized composition which the number of species sum to 1.
        """
        return self / self._natoms

    @property
    def reduced_composition(self) -> 'Composition':
        """
        Returns the reduced composition,i.e. amounts normalized by greatest
        common denominator. e.g., Composition("FePO4") for
        Composition("Fe4P4O16").
        """
        return self.get_reduced_composition_and_factor()[0]

    def get_reduced_composition_and_factor(self) -> Tuple['Composition', float]:
        """
        Calculates a reduced composition and factor.

        Returns:
            A normalized composition and a multiplicative factor, i.e.,
            Li4Fe4P4O16 returns (Composition("LiFePO4"), 4).
        """
        factor = self.get_reduced_formula_and_factor()[1]
        return self / factor, factor

    def get_reduced_formula_and_factor(self, iupac_ordering=False) -> Tuple[str, float]:
        """
        Calculates a reduced formula and factor.

        Args:
            iupac_ordering (bool, optional): Whether to order the
                formula by the iupac "electronegativity" series, defined in
                Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
                Recommendations 2005)". This ordering effectively follows
                the groups and rows of the periodic table, except the
                Lanthanides, Actanides and hydrogen. Note that polyanions
                will still be determined based on the true electronegativity of
                the elements.

        Returns:
            A pretty normalized formula and a multiplicative factor, i.e.,
            Li4Fe4P4O16 returns (LiFePO4, 4).
        """
        all_int = all(abs(x - round(x)) < Composition.amount_tolerance
                      for x in self.values())
        if not all_int:
            return self.formula.replace(" ", ""), 1
        d = {k: int(round(v)) for k, v in self.get_el_amt_dict().items()}
        (formula, factor) = reduce_formula(
            d, iupac_ordering=iupac_ordering)

        if formula in Composition.special_formulas:
            formula = Composition.special_formulas[formula]
            factor /= 2

        return formula, factor

    def get_integer_formula_and_factor(self, max_denominator=10000,
                                       iupac_ordering=False):
        """
        Calculates an integer formula and factor.

        Args:
            max_denominator (int): all amounts in the el:amt dict are
                first converted to a Fraction with this maximum denominator
            iupac_ordering (bool, optional): Whether to order the
                formula by the iupac "electronegativity" series, defined in
                Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
                Recommendations 2005)". This ordering effectively follows
                the groups and rows of the periodic table, except the
                Lanthanides, Actanides and hydrogen. Note that polyanions
                will still be determined based on the true electronegativity of
                the elements.

        Returns:
            A pretty normalized formula and a multiplicative factor, i.e.,
            Li0.5O0.25 returns (Li2O, 0.25). O0.25 returns (O2, 0.125)
        """
        el_amt = self.get_el_amt_dict()
        g = gcd_float(list(el_amt.values()), 1 / max_denominator)

        d = {k: round(v / g) for k, v in el_amt.items()}
        (formula, factor) = reduce_formula(
            d, iupac_ordering=iupac_ordering)
        if formula in Composition.special_formulas:
            formula = Composition.special_formulas[formula]
            factor /= 2
        return formula, factor * g

    @property
    def reduced_formula(self) -> str:
        """
        Returns a pretty normalized formula, i.e., LiFePO4 instead of
        Li4Fe4P4O16.
        """
        return self.get_reduced_formula_and_factor()[0]

    @property
    def hill_formula(self) -> str:
        """
        :return: Hill formula. The Hill system (or Hill notation) is a system
        of writing empirical chemical formulas, molecular chemical formulas and
        components of a condensed formula such that the number of carbon atoms
        in a molecule is indicated first, the number of hydrogen atoms next,
        and then the number of all other chemical elements subsequently, in
        alphabetical order of the chemical symbols. When the formula contains
        no carbon, all the elements, including hydrogen, are listed
        alphabetically.
        """
        c = self.element_composition
        elements = sorted([el.symbol for el in c.keys()])
        if "C" in elements:
            elements = ["C"] + [el for el in elements if el != "C"]

        formula = ["%s%s" % (el, formula_double_format(c[el]) if c[el] != 1 else "")
                   for el in elements]
        return " ".join(formula)

    @property
    def elements(self) -> List[Element]:
        """
        Returns view of elements in Composition.
        """
        return list(self.keys())

    def __str__(self):
        return " ".join([
            "{}{}".format(k, formula_double_format(v, ignore_ones=False))
            for k, v in self.as_dict().items()])

    @property
    def num_atoms(self):
        """
        Total number of atoms in Composition. For negative amounts, sum
        of absolute values
        """
        return self._natoms

    @property
    def weight(self):
        """
        Total molecular weight of Composition
        """
        return Mass(sum([amount * el.atomic_mass for el, amount in self.items()]), "amu")

    def get_atomic_fraction(self, el):
        """
        Calculate atomic fraction of an Element or Specie.

        Args:
            el (Element/Specie): Element or Specie to get fraction for.

        Returns:
            Atomic fraction for element el in Composition
        """
        return abs(self[el]) / self._natoms

    def get_wt_fraction(self, el):
        """
        Calculate weight fraction of an Element or Specie.

        Args:
            el (Element/Specie): Element or Specie to get fraction for.

        Returns:
            Weight fraction for element el in Composition
        """
        return get_el_sp(el).atomic_mass * abs(self[el]) / self.weight

    def contains_element_type(self, category):
        """
        Check if Composition contains any elements matching a given category.

        Args:
            category (str): one of "noble_gas", "transition_metal",
            "post_transition_metal", "rare_earth_metal", "metal", "metalloid",
            "alkali", "alkaline", "halogen", "chalcogen", "lanthanoid",
            "actinoid", "quadrupolar", "s-block", "p-block", "d-block", "f-block"


        Returns:
            True if any elements in Composition match category, otherwise False
        """

        allowed_categories = ("noble_gas", "transition_metal", "post_transition_metal",
                              "rare_earth_metal", "metal", "metalloid", "alkali",
                              "alkaline", "halogen", "chalcogen", "lanthanoid",
                              "actinoid", "quadrupolar", "s-block", "p-block",
                              "d-block", "f-block")

        if category not in allowed_categories:
            raise ValueError("Please pick a category from: {}".format(
                ", ".join(allowed_categories)))

        if "block" in category:
            return any([category[0] in el.block for el in self.elements])
        return any([getattr(el, "is_{}".format(category)) for el in self.elements])

    def _parse_formula(self, formula):
        """
        Args:
            formula (str): A string formula, e.g. Fe2O3, Li3Fe2(PO4)3

        Returns:
            Composition with that formula.

        Notes:
            In the case of Metallofullerene formula (e.g. Y3N@C80),
            the @ mark will be dropped and passed to parser.
        """
        # for Metallofullerene like "Y3N@C80"
        formula = formula.replace("@", "")

        def get_sym_dict(f, factor):
            sym_dict = collections.defaultdict(float)
            for m in re.finditer(r"([A-Z][a-z]*)\s*([-*\.e\d]*)", f):
                el = m.group(1)
                amt = 1
                if m.group(2).strip() != "":
                    amt = float(m.group(2))
                sym_dict[el] += amt * factor
                f = f.replace(m.group(), "", 1)
            if f.strip():
                raise CompositionError("{} is an invalid formula!".format(f))
            return sym_dict

        m = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
        if m:
            factor = 1
            if m.group(2) != "":
                factor = float(m.group(2))
            unit_sym_dict = get_sym_dict(m.group(1), factor)
            expanded_sym = "".join(["{}{}".format(el, amt)
                                    for el, amt in unit_sym_dict.items()])
            expanded_formula = formula.replace(m.group(), expanded_sym)
            return self._parse_formula(expanded_formula)
        return get_sym_dict(formula, 1)

    @property
    def anonymized_formula(self):
        """
        An anonymized formula. Unique species are arranged in ordering of
        increasing amounts and assigned ascending alphabets. Useful for
        prototyping formulas. For example, all stoichiometric perovskites have
        anonymized_formula ABC3.
        """
        reduced = self.element_composition
        if all(x == int(x) for x in self.values()):
            reduced /= gcd(*(int(i) for i in self.values()))

        anon = ""
        for e, amt in zip(string.ascii_uppercase, sorted(reduced.values())):
            if amt == 1:
                amt_str = ""
            elif abs(amt % 1) < 1e-8:
                amt_str = str(int(amt))
            else:
                amt_str = str(amt)
            anon += ("{}{}".format(e, amt_str))
        return anon

    @property
    def chemical_system(self) -> str:
        """
        Get the chemical system of a Composition, for example "O-Si" for
        SiO2. Chemical system is a string of a list of elements
        sorted alphabetically and joined by dashes, by convention for use
        in database keys.
        """
        return "-".join(sorted([str(el) for el in self.elements]))

    @property
    def valid(self):
        """
        Returns True if Composition contains valid elements or species and
        False if the Composition contains any dummy species.
        """
        return not any([isinstance(el, DummySpecie) for el in self.elements])

    def __repr__(self):
        return "Comp: " + self.formula

    @classmethod
    def from_dict(cls, d):
        """
        Creates a composition from a dict generated by as_dict(). Strictly not
        necessary given that the standard constructor already takes in such an
        input, but this method preserves the standard pymatgen API of having
        from_dict methods to reconstitute objects generated by as_dict(). Allows
        for easier introspection.

        Args:
            d (dict): {symbol: amount} dict.
        """
        return cls(d)

    def get_el_amt_dict(self):
        """
        Returns:
            Dict with element symbol and (unreduced) amount e.g.,
            {"Fe": 4.0, "O":6.0} or {"Fe3+": 4.0, "O2-":6.0}
        """
        d = collections.defaultdict(float)
        for e, a in self.items():
            d[e.symbol] += a
        return d

    def as_dict(self):
        """
        Returns:
            dict with species symbol and (unreduced) amount e.g.,
            {"Fe": 4.0, "O":6.0} or {"Fe3+": 4.0, "O2-":6.0}
        """
        d = collections.defaultdict(float)
        for e, a in self.items():
            d[str(e)] += a
        return d

    @property
    def to_reduced_dict(self):
        """
        Returns:
            Dict with element symbol and reduced amount e.g.,
            {"Fe": 2.0, "O":3.0}
        """
        return self.get_reduced_composition_and_factor()[0]

    @property
    def to_data_dict(self):
        """
        Returns:
            A dict with many keys and values relating to Composition/Formula,
            including reduced_cell_composition, unit_cell_composition,
            reduced_cell_formula, elements and nelements.
        """
        return {
            "reduced_cell_composition": self.get_reduced_composition_and_factor()[0],
            "unit_cell_composition": self.as_dict(),
            "reduced_cell_formula": self.reduced_formula,
            "elements": list(self.as_dict().keys()),
            "nelements": len(self.as_dict().keys())
        }

    def oxi_state_guesses(self, oxi_states_override=None, target_charge=0,
                          all_oxi_states=False, max_sites=None):
        """
        Checks if the composition is charge-balanced and returns back all
        charge-balanced oxidation state combinations. Composition must have
        integer values. Note that more num_atoms in the composition gives
        more degrees of freedom. e.g., if possible oxidation states of
        element X are [2,4] and Y are [-3], then XY is not charge balanced
        but X2Y2 is. Results are returned from most to least probable based
        on ICSD statistics. Use max_sites to improve performance if needed.

        Args:
            oxi_states_override (dict): dict of str->list to override an
                element's common oxidation states, e.g. {"V": [2,3,4,5]}
            target_charge (int): the desired total charge on the structure.
                Default is 0 signifying charge balance.
            all_oxi_states (bool): if True, an element defaults to
                all oxidation states in pymatgen Element.icsd_oxidation_states.
                Otherwise, default is Element.common_oxidation_states. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results.
            max_sites (int): if possible, will reduce Compositions to at most
                this many sites to speed up oxidation state guesses. If the
                composition cannot be reduced to this many sites a ValueError
                will be raised. Set to -1 to just reduce fully. If set to a
                number less than -1, the formula will be fully reduced but a
                ValueError will be thrown if the number of atoms in the reduced
                formula is greater than abs(max_sites).

        Returns:
            A list of dicts - each dict reports an element symbol and average
                oxidation state across all sites in that composition. If the
                composition is not charge balanced, an empty list is returned.
        """

        return self._get_oxid_state_guesses(all_oxi_states, max_sites, oxi_states_override, target_charge)[0]

    def add_charges_from_oxi_state_guesses(self,
                                           oxi_states_override=None,
                                           target_charge=0,
                                           all_oxi_states=False,
                                           max_sites=None):
        """
        Assign oxidation states basedon guessed oxidation states.

        See `oxi_state_guesses` for an explanation of how oxidation states are
        guessed. This operation uses the set of oxidation states for each site
        that were determined to be most likley from the oxidation state guessing
        routine.

        Args:
            oxi_states_override (dict): dict of str->list to override an
                element's common oxidation states, e.g. {"V": [2,3,4,5]}
            target_charge (int): the desired total charge on the structure.
                Default is 0 signifying charge balance.
            all_oxi_states (bool): if True, an element defaults to
                all oxidation states in pymatgen Element.icsd_oxidation_states.
                Otherwise, default is Element.common_oxidation_states. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results.
            max_sites (int): if possible, will reduce Compositions to at most
                this many sites to speed up oxidation state guesses. If the
                composition cannot be reduced to this many sites a ValueError
                will be raised. Set to -1 to just reduce fully. If set to a
                number less than -1, the formula will be fully reduced but a
                ValueError will be thrown if the number of atoms in the reduced
                formula is greater than abs(max_sites).

        Returns:
            Composition, where the elements are assigned oxidation states based
            on the results form guessing oxidation states. If no oxidation state
            is possible, returns a Composition where all oxidation states are 0.
        """

        _, oxidation_states = self._get_oxid_state_guesses(
            all_oxi_states, max_sites, oxi_states_override, target_charge)

        # Special case: No charged compound is possible
        if not oxidation_states:
            return Composition(dict((Specie(e, 0), f) for e, f in self.items()))

        # Generate the species
        species = []
        for el, charges in oxidation_states[0].items():
            species.extend([Specie(el, c) for c in charges])

        # Return the new object
        return Composition(collections.Counter(species))

    def remove_charges(self):
        """
        Removes the charges from any species in a Composition object.

        Returns:
            Composition object without charge decoration, for example
            {"Fe3+": 2.0, "O2-":3.0} becomes {"Fe": 2.0, "O":3.0}
        """
        d = collections.Counter()

        for e, f in self.items():
            e = re.findall(r"[A-z]+", str(e))[0]
            d[str(e)] += f

        return Composition(d)

    def _get_oxid_state_guesses(self, all_oxi_states, max_sites,
                                oxi_states_override, target_charge):
        """
        Utility operation for guessing oxidation states.

        See `oxi_state_guesses` for full details. This operation does the
        calculation of the most likely oxidation states

        Args:
            oxi_states_override (dict): dict of str->list to override an
                element's common oxidation states, e.g. {"V": [2,3,4,5]}
            target_charge (int): the desired total charge on the structure.
                Default is 0 signifying charge balance.
            all_oxi_states (bool): if True, an element defaults to
                all oxidation states in pymatgen Element.icsd_oxidation_states.
                Otherwise, default is Element.common_oxidation_states. Note
                that the full oxidation state list is *very* inclusive and
                can produce nonsensical results.
            max_sites (int): if possible, will reduce Compositions to at most
                this many sites to speed up oxidation state guesses. If the
                composition cannot be reduced to this many sites a ValueError
                will be raised. Set to -1 to just reduce fully. If set to a
                number less than -1, the formula will be fully reduced but a
                ValueError will be thrown if the number of atoms in the reduced
                formula is greater than abs(max_sites).

        Returns:
            A list of dicts - each dict reports an element symbol and average
                oxidation state across all sites in that composition. If the
                composition is not charge balanced, an empty list is returned.
            A list of dicts - each dict maps the element symbol to a list of
                oxidation states for each site of that element. For example, Fe3O4 could
                return a list of [2,2,2,3,3,3] for the oxidation states of If the composition
                is

            """
        comp = self.copy()
        # reduce Composition if necessary
        if max_sites and max_sites < 0:
            comp = self.reduced_composition

            if max_sites < -1 and comp.num_atoms > abs(max_sites):
                raise ValueError(
                    "Composition {} cannot accommodate max_sites "
                    "setting!".format(comp))

        elif max_sites and comp.num_atoms > max_sites:
            reduced_comp, reduced_factor = self. \
                get_reduced_composition_and_factor()
            if reduced_factor > 1:
                reduced_comp *= max(1, int(max_sites / reduced_comp.num_atoms))
                comp = reduced_comp  # as close to max_sites as possible
            if comp.num_atoms > max_sites:
                raise ValueError("Composition {} cannot accommodate max_sites "
                                 "setting!".format(comp))

        # Load prior probabilities of oxidation states, used to rank solutions
        if not Composition.oxi_prob:
            module_dir = os.path.join(os.path.
                                      dirname(os.path.abspath(__file__)))
            all_data = loadfn(os.path.join(module_dir, "..",
                                           "analysis", "icsd_bv.yaml"))
            Composition.oxi_prob = {Specie.from_string(sp): data
                                    for sp, data in
                                    all_data["occurrence"].items()}
        oxi_states_override = oxi_states_override or {}
        # assert: Composition only has integer amounts
        if not all(amt == int(amt) for amt in comp.values()):
            raise ValueError("Charge balance analysis requires integer "
                             "values in Composition!")

        # for each element, determine all possible sum of oxidations
        # (taking into account nsites for that particular element)
        el_amt = comp.get_el_amt_dict()
        els = el_amt.keys()
        el_sums = []  # matrix: dim1= el_idx, dim2=possible sums
        el_sum_scores = collections.defaultdict(set)  # dict of el_idx, sum -> score
        el_best_oxid_combo = {}  # dict of el_idx, sum -> oxid combo with best score
        for idx, el in enumerate(els):
            el_sum_scores[idx] = {}
            el_best_oxid_combo[idx] = {}
            el_sums.append([])
            if oxi_states_override.get(el):
                oxids = oxi_states_override[el]
            elif all_oxi_states:
                oxids = Element(el).oxidation_states
            else:
                oxids = Element(el).icsd_oxidation_states or \
                        Element(el).oxidation_states

            # get all possible combinations of oxidation states
            # and sum each combination
            for oxid_combo in combinations_with_replacement(oxids,
                                                            int(el_amt[el])):

                # List this sum as a possible option
                oxid_sum = sum(oxid_combo)
                if oxid_sum not in el_sums[idx]:
                    el_sums[idx].append(oxid_sum)

                # Determine how probable is this combo?
                score = sum([Composition.oxi_prob.get(Specie(el, o), 0) for
                             o in oxid_combo])

                # If it is the most probable combo for a certain sum,
                #   store the combination
                if oxid_sum not in el_sum_scores[idx] or score > el_sum_scores[idx].get(oxid_sum, 0):
                    el_sum_scores[idx][oxid_sum] = score
                    el_best_oxid_combo[idx][oxid_sum] = oxid_combo

        # Determine which combination of oxidation states for each element
        #    is the most probable
        all_sols = []  # will contain all solutions
        all_oxid_combo = []  # will contain the best combination of oxidation states for each site
        all_scores = []  # will contain a score for each solution
        for x in product(*el_sums):
            # each x is a trial of one possible oxidation sum for each element
            if sum(x) == target_charge:  # charge balance condition
                el_sum_sol = dict(zip(els, x))  # element->oxid_sum
                # normalize oxid_sum by amount to get avg oxid state
                sol = {el: v / el_amt[el] for el, v in el_sum_sol.items()}
                # add the solution to the list of solutions
                all_sols.append(sol)

                # determine the score for this solution
                score = 0
                for idx, v in enumerate(x):
                    score += el_sum_scores[idx][v]
                all_scores.append(score)

                # collect the combination of oxidation states for each site
                all_oxid_combo.append(
                    dict((e, el_best_oxid_combo[idx][v]) for idx, (e, v) in enumerate(zip(els, x))))

        # sort the solutions by highest to lowest score
        if all_scores:
            all_sols, all_oxid_combo = zip(*[(y, x) for (z, y, x) in sorted(zip(all_scores, all_sols, all_oxid_combo),
                                                                            key=lambda pair: pair[0],
                                                                            reverse=True)])
        return all_sols, all_oxid_combo

    @staticmethod
    def ranked_compositions_from_indeterminate_formula(fuzzy_formula,
                                                       lock_if_strict=True):
        """
        Takes in a formula where capitilization might not be correctly entered,
        and suggests a ranked list of potential Composition matches.
        Author: Anubhav Jain

        Args:
            fuzzy_formula (str): A formula string, such as "co2o3" or "MN",
                that may or may not have multiple interpretations
            lock_if_strict (bool): If true, a properly entered formula will
                only return the one correct interpretation. For example,
                "Co1" will only return "Co1" if true, but will return both
                "Co1" and "C1 O1" if false.

        Returns:
            A ranked list of potential Composition matches
        """

        # if we have an exact match and the user specifies lock_if_strict, just
        # return the exact match!
        if lock_if_strict:
            # the strict composition parsing might throw an error, we can ignore
            # it and just get on with fuzzy matching
            try:
                comp = Composition(fuzzy_formula)
                return [comp]
            except (CompositionError, ValueError):
                pass

        all_matches = Composition._comps_from_fuzzy_formula(fuzzy_formula)
        # remove duplicates
        all_matches = list(set(all_matches))
        # sort matches by rank descending
        all_matches = sorted(all_matches,
                             key=lambda match: (match[1], match[0]), reverse=True)
        all_matches = [m[0] for m in all_matches]
        return all_matches

    @staticmethod
    def _comps_from_fuzzy_formula(fuzzy_formula, m_dict=None, m_points=0,
                                  factor=1):
        """
        A recursive helper method for formula parsing that helps in
        interpreting and ranking indeterminate formulas.
        Author: Anubhav Jain

        Args:
            fuzzy_formula (str): A formula string, such as "co2o3" or "MN",
                that may or may not have multiple interpretations.
            m_dict (dict): A symbol:amt dictionary from the previously parsed
                formula.
            m_points: Number of points gained from the previously parsed
                formula.
            factor: Coefficient for this parse, e.g. (PO4)2 will feed in PO4
                as the fuzzy_formula with a coefficient of 2.

        Returns:
            A list of tuples, with the first element being a Composition and
            the second element being the number of points awarded that
            Composition intepretation.
        """
        m_dict = m_dict or {}

        def _parse_chomp_and_rank(m, f, m_dict, m_points):
            """
            A helper method for formula parsing that helps in interpreting and
            ranking indeterminate formulas
            Author: Anubhav Jain

            Args:
                m: A regex match, with the first group being the element and
                    the second group being the amount
                f: The formula part containing the match
                m_dict: A symbol:amt dictionary from the previously parsed
                    formula
                m_points: Number of points gained from the previously parsed
                    formula

            Returns:
                A tuple of (f, m_dict, points) where m_dict now contains data
                from the match and the match has been removed (chomped) from
                the formula f. The "goodness" of the match determines the
                number of points returned for chomping. Returns
                (None, None, None) if no element could be found...
            """

            points = 0
            # Points awarded if the first element of the element is correctly
            # specified as a capital
            points_first_capital = 100
            # Points awarded if the second letter of the element is correctly
            # specified as lowercase
            points_second_lowercase = 100

            # get element and amount from regex match
            el = m.group(1)
            if len(el) > 2 or len(el) < 1:
                raise CompositionError("Invalid element symbol entered!")
            amt = float(m.group(2)) if m.group(2).strip() != "" else 1

            # convert the element string to proper [uppercase,lowercase] format
            # and award points if it is already in that format
            char1 = el[0]
            char2 = el[1] if len(el) > 1 else ""

            if char1 == char1.upper():
                points += points_first_capital
            if char2 and char2 == char2.lower():
                points += points_second_lowercase

            el = char1.upper() + char2.lower()

            # if it's a valid element, chomp and add to the points
            if Element.is_valid_symbol(el):
                if el in m_dict:
                    m_dict[el] += amt * factor
                else:
                    m_dict[el] = amt * factor
                return f.replace(m.group(), "", 1), m_dict, m_points + points

            # else return None
            return None, None, None

        fuzzy_formula = fuzzy_formula.strip()

        if len(fuzzy_formula) == 0:
            # The entire formula has been parsed into m_dict. Return the
            # corresponding Composition and number of points
            if m_dict:
                yield (Composition.from_dict(m_dict), m_points)
        else:
            # if there is a parenthesis, remove it and match the remaining stuff
            # with the appropriate factor
            for mp in re.finditer(r"\(([^\(\)]+)\)([\.\d]*)", fuzzy_formula):
                mp_points = m_points
                mp_form = fuzzy_formula.replace(mp.group(), " ", 1)
                mp_dict = dict(m_dict)
                mp_factor = 1 if mp.group(2) == "" else float(mp.group(2))
                # Match the stuff inside the parenthesis with the appropriate
                # factor
                for match in \
                        Composition._comps_from_fuzzy_formula(mp.group(1),
                                                              mp_dict,
                                                              mp_points,
                                                              factor=mp_factor):
                    only_me = True
                    # Match the stuff outside the parentheses and return the
                    # sum.

                    for match2 in \
                            Composition._comps_from_fuzzy_formula(mp_form,
                                                                  mp_dict,
                                                                  mp_points,
                                                                  factor=1):
                        only_me = False
                        yield (match[0] + match2[0], match[1] + match2[1])
                    # if the stuff inside the parenthesis is nothing, then just
                    # return the stuff inside the parentheses
                    if only_me:
                        yield match
                return

            # try to match the single-letter elements
            m1 = re.match(r"([A-z])([\.\d]*)", fuzzy_formula)
            if m1:
                m_points1 = m_points
                m_form1 = fuzzy_formula
                m_dict1 = dict(m_dict)
                (m_form1, m_dict1, m_points1) = \
                    _parse_chomp_and_rank(m1, m_form1, m_dict1, m_points1)
                if m_dict1:
                    # there was a real match
                    for match in \
                            Composition._comps_from_fuzzy_formula(m_form1,
                                                                  m_dict1,
                                                                  m_points1,
                                                                  factor):
                        yield match

            # try to match two-letter elements
            m2 = re.match(r"([A-z]{2})([\.\d]*)", fuzzy_formula)
            if m2:
                m_points2 = m_points
                m_form2 = fuzzy_formula
                m_dict2 = dict(m_dict)
                (m_form2, m_dict2, m_points2) = \
                    _parse_chomp_and_rank(m2, m_form2, m_dict2, m_points2)
                if m_dict2:
                    # there was a real match
                    for match in \
                            Composition._comps_from_fuzzy_formula(m_form2, m_dict2,
                                                                  m_points2,
                                                                  factor):
                        yield match


def reduce_formula(sym_amt, iupac_ordering=False):
    """
    Helper method to reduce a sym_amt dict to a reduced formula and factor.

    Args:
        sym_amt (dict): {symbol: amount}.
        iupac_ordering (bool, optional): Whether to order the
            formula by the iupac "electronegativity" series, defined in
            Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
            Recommendations 2005)". This ordering effectively follows
            the groups and rows of the periodic table, except the
            Lanthanides, Actanides and hydrogen. Note that polyanions
            will still be determined based on the true electronegativity of
            the elements.

    Returns:
        (reduced_formula, factor).
    """
    syms = sorted(sym_amt.keys(), key=lambda x: [get_el_sp(x).X, x])

    syms = list(filter(
        lambda x: abs(sym_amt[x]) > Composition.amount_tolerance, syms))

    factor = 1
    # Enforce integers for doing gcd.
    if all((int(i) == i for i in sym_amt.values())):
        factor = abs(gcd(*(int(i) for i in sym_amt.values())))

    polyanion = []
    # if the composition contains a poly anion
    if len(syms) >= 3 and get_el_sp(syms[-1]).X - get_el_sp(syms[-2]).X < 1.65:
        poly_sym_amt = {syms[i]: sym_amt[syms[i]] / factor
                        for i in [-2, -1]}
        (poly_form, poly_factor) = reduce_formula(
            poly_sym_amt, iupac_ordering=iupac_ordering)

        if poly_factor != 1:
            polyanion.append("({}){}".format(poly_form, int(poly_factor)))

    syms = syms[:len(syms) - 2 if polyanion else len(syms)]

    if iupac_ordering:
        syms = sorted(syms,
                      key=lambda x: [get_el_sp(x).iupac_ordering, x])

    reduced_form = []
    for s in syms:
        normamt = sym_amt[s] * 1.0 / factor
        reduced_form.append(s)
        reduced_form.append(formula_double_format(normamt))

    reduced_form = "".join(reduced_form + polyanion)
    return reduced_form, factor


class CompositionError(Exception):
    """Exception class for composition errors"""


class ChemicalPotential(dict, MSONable):
    """
    Class to represent set of chemical potentials. Can be:
    multiplied/divided by a Number
    multiplied by a Composition (returns an energy)
    added/subtracted with other ChemicalPotentials.
    """

    def __init__(self, *args, **kwargs):
        """
        Args:
            *args, **kwargs: any valid dict init arguments
        """
        d = dict(*args, **kwargs)
        super().__init__((get_el_sp(k), v)
                         for k, v in d.items())
        if len(d) != len(self):
            raise ValueError("Duplicate potential specified")

    def __mul__(self, other):
        if isinstance(other, numbers.Number):
            return ChemicalPotential({k: v * other for k, v in self.items()})
        raise NotImplementedError()

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, numbers.Number):
            return ChemicalPotential({k: v / other for k, v in self.items()})
        raise NotImplementedError()

    __div__ = __truediv__

    def __sub__(self, other):
        if isinstance(other, ChemicalPotential):
            els = set(self.keys()).union(other.keys())
            return ChemicalPotential({e: self.get(e, 0) - other.get(e, 0)
                                      for e in els})
        raise NotImplementedError()

    def __add__(self, other):
        if isinstance(other, ChemicalPotential):
            els = set(self.keys()).union(other.keys())
            return ChemicalPotential({e: self.get(e, 0) + other.get(e, 0)
                                      for e in els})
        raise NotImplementedError()

    def get_energy(self, composition, strict=True):
        """
        Calculates the energy of a composition.

        Args:
            composition (Composition): input composition
            strict (bool): Whether all potentials must be specified
        """
        if strict and set(composition.keys()) > set(self.keys()):
            s = set(composition.keys()) - set(self.keys())
            raise ValueError("Potentials not specified for {}".format(s))
        return sum(self.get(k, 0) * v for k, v in composition.items())

    def __repr__(self):
        return "ChemPots: " + super().__repr__()


if __name__ == "__main__":
    import doctest

    doctest.testmod()
