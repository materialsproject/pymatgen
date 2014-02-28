#!/usr/bin/env python

"""
This module implements a Composition class to represent compositions.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Nov 10, 2012"

import re
import collections
import string
from fractions import gcd
from itertools import chain
from pymatgen.core.periodic_table import get_el_sp, Element
from pymatgen.util.string_utils import formula_double_format
from pymatgen.serializers.json_coders import MSONable
from pymatgen.core.units import unitized


class Composition(collections.Mapping, collections.Hashable, MSONable):
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

    """
    Tolerance in distinguishing different composition amounts.
    1e-8 is fairly tight, but should cut out most floating point arithmetic
    errors.
    """
    amount_tolerance = 1e-8

    """
    Special formula handling for peroxides and certain elements. This is so
    that formula output does not write LiO instead of Li2O2 for example.
    """
    special_formulas = {"LiO": "Li2O2", "NaO": "Na2O2", "KO": "K2O2",
                        "HO": "H2O2", "CsO": "Cs2O2", "RbO": "Rb2O2",
                        "O": "O2",  "N": "N2", "Cl": "Cl2", "H": "H2"}

    def __init__(self, *args, **kwargs):
        """
        Very flexible Composition construction, similar to the built-in Python
        dict(). Also extended to allow simple string init.

        Args:
            Any form supported by the Python built-in dict() function.

            1. A dict of either {Element/Specie: amount},

               {string symbol:amount}, or {atomic number:amount} or any mixture
               of these. E.g., {Element("Li"):2 ,Element("O"):1},
               {"Li":2, "O":1}, {3:2, 8:1} all result in a Li2O composition.
            2. Keyword arg initialization, similar to a dict, e.g.,

               Compostion(Li = 2, O = 1)

            In addition, the Composition constructor also allows a single
            string as an input formula. E.g., Composition("Li2O").
        """
        if len(args) == 1 and isinstance(args[0], basestring):
            elmap = self._parse_formula(args[0])
        else:
            elmap = dict(*args, **kwargs)
        for k, v in elmap.items():
            if v < -Composition.amount_tolerance:
                raise CompositionError("Amounts in Composition cannot be "
                                       "negative!")
            elif v < 0:
                del elmap[k]
        self._elmap = {get_el_sp(k): v for k, v in elmap.items()}
        self._natoms = sum(self._elmap.values())

    def __getitem__(self, el):
        """
        Get the amount for element.
        """
        return self._elmap.get(get_el_sp(el), 0)

    def __eq__(self, other):
        for el in chain(self.elements, other.elements):
            if abs(self[el] - other[el]) > Composition.amount_tolerance:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        """
        Adds two compositions. For example, an Fe2O3 composition + an FeO
        composition gives a Fe3O4 composition.
        """
        new_el_map = {el: self[el] for el in self}
        for k in other.keys():
            el = get_el_sp(k)
            if el in self:
                new_el_map[el] += other[k]
            else:
                new_el_map[el] = other[k]
        return Composition(new_el_map)

    def __sub__(self, other):
        """
        Subtracts two compositions. For example, an Fe2O3 composition - an FeO
        composition gives an FeO2 composition.

        Raises:
            CompositionError if the subtracted composition is greater than the
            original composition in any of its elements.
        """
        new_el_map = {el: self[el] for el in self}
        for k in other.keys():
            el = get_el_sp(k)
            if el in self and other[k] <= self[el]:
                new_el_map[el] -= other[k]
            else:
                raise CompositionError(
                    "All elements in subtracted composition must exist in "
                    "original composition in equal or lesser amount!")

            new_el_map = {sp: amt for sp, amt in new_el_map.iteritems()
                          if amt != 0}
        return Composition(new_el_map)

    def __mul__(self, other):
        """
        Multiply a Composition by an integer or a float.
        Fe2O3 * 4 -> Fe8O12
        """
        if not isinstance(other, (int, float)):
            raise ValueError("Multiplication can only be done for int/floats!")
        return Composition({el: self[el] * other for el in self})

    def __rmul__(self, other):
        """
        Multiply a Composition by an integer or a float. This provides for
        the reflected multiplication, e.g.,
        4 * Fe2O3 -> Fe8O12
        """
        if not isinstance(other, (int, float)):
            raise ValueError("Multiplication can only be done for int/floats!")
        return Composition({el: self[el] * other for el in self})

    def __hash__(self):
        """
        Minimally effective hash function that just distinguishes between
        Compositions with different elements.
        """
        hashcode = 0
        for el in self._elmap.keys():
            #Ignore elements with zero amounts.
            if self[el] > self.amount_tolerance:
                hashcode += el.Z
        return hashcode

    def __contains__(self, el):
        return el in self._elmap

    def __len__(self):
        return len(self._elmap)

    def __iter__(self):
        return self._elmap.__iter__()

    def arb_ordered_elmap(self):
        """
        Arbitrary ordered elmap on the elements/species of a composition of a given site
        in an unordered structure. Returns a list of tuples (element_or_specie: occupation)
        in the arbitrary order.
        The arbitrary order is based on the Z of the element and the smallest
        fractional occupations first.
        Example : {"Ni3+": 0.2, "Ni4+": 0.2, "Cr3+": 0.15, "Zn2+": 0.34, "Cr4+": 0.11}
                will yield the species in the following order :
                    Cr4+, Cr3+, Ni3+, Ni4+, Zn2+ ... or
                    Cr4+, Cr3+, Ni4+, Ni3+, Zn2+
        This method is ONLY used in the BVAnalyser for unordered structures
        """
        sorted_elmap_keys = sorted(self._elmap.keys(),
                                   key=lambda elsp: 2.0*float(elsp.Z) + self._elmap[elsp])
        return [(elsp, self._elmap[elsp]) for elsp in sorted_elmap_keys]

    @property
    def average_electroneg(self):
        return sum((el.X * amt for el, amt in self._elmap.items())) / \
            self.num_atoms

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
    def is_element(self):
        """
        True if composition is for an element.
        """
        positive_amts = [amt for amt in self._elmap.values()
                         if amt > self.amount_tolerance]
        return len(positive_amts) == 1

    def copy(self):
        return Composition(self._elmap)

    @property
    def formula(self):
        """
        Returns a formula string, with elements sorted by electronegativity,
        e.g., Li4 Fe4 P4 O16.
        """
        sym_amt = self.get_el_amt_dict()
        syms = sorted(sym_amt.keys(),
                      key=lambda s: get_el_sp(s).X)
        formula = []
        for s in syms:
            if sym_amt[s] != 0:
                formula.append(s + formula_double_format(sym_amt[s], False))
        return " ".join(formula)

    @property
    def alphabetical_formula(self):
        """
        Returns a formula string, with elements sorted by alphabetically
        e.g., Fe4 Li4 O16 P4.
        """
        sym_amt = self.get_el_amt_dict()
        syms = sorted(sym_amt.keys())
        formula = []
        for s in syms:
            if sym_amt[s] != 0:
                formula.append(s + formula_double_format(sym_amt[s], False))
        return " ".join(formula)

    @property
    def element_composition(self):
        """
        Returns the composition replacing any species by the corresponding
        element.
        """
        o = collections.defaultdict(float)
        for sp in self:
            o[Element(sp.symbol)] += self[sp]
        return Composition(o)

    @property
    def reduced_composition(self):
        """
        Returns the reduced composition,i.e. amounts normalized by greatest
        common denominator. e.g., Composition("FePO4") for
        Composition("Fe4P4O16").
        """
        return self.get_reduced_composition_and_factor()[0]

    def get_reduced_composition_and_factor(self):
        """
        Calculates a reduced composition and factor.

        Returns:
            A normalized composition and a multiplicative factor, i.e.,
            Li4Fe4P4O16 returns (Composition("LiFePO4"), 4).
        """
        factor = self.get_reduced_formula_and_factor()[1]
        reduced_comp = Composition({el: self[el] / factor for el in self})
        return reduced_comp, factor

    def get_reduced_formula_and_factor(self):
        """
        Calculates a reduced formula and factor.

        Returns:
            A pretty normalized formula and a multiplicative factor, i.e.,
            Li4Fe4P4O16 returns (LiFePO4, 4).
        """
        all_int = all([x == int(x) for x in self._elmap.values()])
        if not all_int:
            return self.formula.replace(" ", ""), 1
        d = self.get_el_amt_dict()
        (formula, factor) = reduce_formula(d)

        if formula in Composition.special_formulas:
            formula = Composition.special_formulas[formula]
            factor /= 2

        return formula, factor

    def get_fractional_composition(self):
        """
        Returns the normalized composition which the number of species sum to
        1.

        Returns:
            Normalized composition which the number of species sum to 1.
        """
        natoms = self._natoms
        frac_map = {k: v / natoms for k, v in self._elmap.items()}
        return Composition(frac_map)

    @property
    def reduced_formula(self):
        """
        Returns a pretty normalized formula, i.e., LiFePO4 instead of
        Li4Fe4P4O16.
        """
        return self.get_reduced_formula_and_factor()[0]

    @property
    def elements(self):
        """
        Returns view of elements in Composition.
        """
        return self._elmap.keys()

    def __str__(self):
        return " ".join([
            "{}{}".format(k, formula_double_format(v, ignore_ones=False))
            for k, v in self.to_dict.items()])

    @property
    def num_atoms(self):
        """
        Total number of atoms in Composition
        """
        return self._natoms

    @property
    @unitized("amu")
    def weight(self):
        """
        Total molecular weight of Composition
        """
        return sum([amount * el.atomic_mass
                    for el, amount in self._elmap.items()])

    def get_atomic_fraction(self, el):
        """
        Calculate atomic fraction of an Element or Specie.

        Args:
            el (Element/Specie): Element or Specie to get fraction for.

        Returns:
            Atomic fraction for element el in Composition
        """
        return self[el] / self._natoms

    def get_wt_fraction(self, el):
        """
        Calculate weight fraction of an Element or Specie.

        Args:
            el (Element/Specie): Element or Specie to get fraction for.

        Returns:
            Weight fraction for element el in Composition
        """
        return get_el_sp(el).atomic_mass * self[el] / self.weight

    def _parse_formula(self, formula):
        """
        Args:
            formula (str): A string formula, e.g. Fe2O3, Li3Fe2(PO4)3

        Returns:
            Composition with that formula.
        """
        def get_sym_dict(f, factor):
            sym_dict = {}
            for m in re.finditer(r"([A-Z][a-z]*)([\.\d]*)", f):
                el = m.group(1)
                amt = 1
                if m.group(2).strip() != "":
                    amt = float(m.group(2))
                if el in sym_dict:
                    sym_dict[el] += amt * factor
                else:
                    sym_dict[el] = amt * factor
                f = f.replace(m.group(), "", 1)
            if f.strip():
                raise CompositionError("{} is an invalid formula!".format(f))
            return sym_dict
        m = re.search(r"\(([^\(\)]+)\)([\.\d]*)", formula)
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

    @staticmethod
    def from_formula(formula):
        """
        .. deprecated:: 1.6.1

        Use Composition(formula) instead.
        """
        return Composition(formula)

    @property
    def anonymized_formula(self):
        """
        An anonymized formula. Unique species are arranged in ordering of
        increasing amounts and assigned ascending alphabets. Useful for
        prototyping formulas. For example, all stoichiometric perovskites have
        anonymized_formula ABC3.
        """
        reduced_comp = self.get_reduced_composition_and_factor()[0]
        els = sorted(reduced_comp.elements, key=lambda e: reduced_comp[e])
        anon_formula = []
        for anon, e in zip(string.ascii_uppercase, els):
            amt = reduced_comp[e]
            if amt > 0:
                if amt == 1:
                    amt_str = ""
                elif abs(amt % 1) < 1e-8:
                    amt_str = str(int(amt))
                else:
                    amt_str = str(amt)
                anon_formula.append("{}{}".format(anon, amt_str))
        return "".join(anon_formula)

    def __repr__(self):
        return "Comp: " + self.formula

    @classmethod
    def from_dict(cls, d):
        """
        Creates a composition from a dict generated by to_dict. Strictly not
        necessary given that the standard constructor already takes in such an
        input, but this method preserves the standard pymatgen API of having
        from_dict methods to reconstitute objects generated by to_dict. Allows
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

    @property
    def to_dict(self):
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
        c = Composition(self.reduced_formula)
        return c.to_dict

    @property
    def to_data_dict(self):
        """
        Returns:
            A dict with many keys and values relating to Composition/Formula,
            including reduced_cell_composition, unit_cell_composition,
            reduced_cell_formula, elements and nelements.
        """
        return {"reduced_cell_composition": self.to_reduced_dict,
                "unit_cell_composition": self.to_dict,
                "reduced_cell_formula": self.reduced_formula,
                "elements": self.to_dict.keys(),
                "nelements": len(self.to_dict.keys())}

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

        #if we have an exact match and the user specifies lock_if_strict, just
        #return the exact match!
        if lock_if_strict:
            #the strict composition parsing might throw an error, we can ignore
            #it and just get on with fuzzy matching
            try:
                comp = Composition.from_formula(fuzzy_formula)
                return [comp]
            except (CompositionError, ValueError):
                pass

        all_matches = Composition._comps_from_fuzzy_formula(fuzzy_formula)
        #remove duplicates
        all_matches = list(set(all_matches))
        #sort matches by rank descending
        all_matches = sorted(all_matches,
                             key=lambda match: match[1], reverse=True)
        all_matches = [m[0] for m in all_matches]
        return all_matches

    @staticmethod
    def _comps_from_fuzzy_formula(fuzzy_formula, m_dict={}, m_points=0,
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

            #get element and amount from regex match
            el = m.group(1)
            if len(el) > 2 or len(el) < 1:
                raise CompositionError("Invalid element symbol entered!")
            amt = float(m.group(2)) if m.group(2).strip() != "" else 1

            #convert the element string to proper [uppercase,lowercase] format
            #and award points if it is already in that format
            char1 = el[0]
            char2 = el[1] if len(el) > 1 else ""

            if char1 == char1.upper():
                points += points_first_capital
            if char2 and char2 == char2.lower():
                points += points_second_lowercase

            el = char1.upper() + char2.lower()

            #if it's a valid element, chomp and add to the points
            if Element.is_valid_symbol(el):
                if el in m_dict:
                    m_dict[el] += amt * factor
                else:
                    m_dict[el] = amt * factor
                return f.replace(m.group(), "", 1), m_dict, m_points + points

            #else return None
            return None, None, None

        fuzzy_formula = fuzzy_formula.strip()

        if len(fuzzy_formula) == 0:
            #The entire formula has been parsed into m_dict. Return the
            #corresponding Composition and number of points
            if m_dict:
                yield (Composition.from_dict(m_dict), m_points)
        else:
            #if there is a parenthesis, remove it and match the remaining stuff
            #with the appropriate factor
            for mp in re.finditer(r"\(([^\(\)]+)\)([\.\d]*)", fuzzy_formula):
                mp_points = m_points
                mp_form = fuzzy_formula.replace(mp.group(), " ", 1)
                mp_dict = dict(m_dict)
                mp_factor = 1 if mp.group(2) == "" else float(mp.group(2))
                #Match the stuff inside the parenthesis with the appropriate
                #factor
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
                    #if the stuff inside the parenthesis is nothing, then just
                    #return the stuff inside the parentheses
                    if only_me:
                        yield match
                return

            #try to match the single-letter elements
            m1 = re.match(r"([A-z])([\.\d]*)", fuzzy_formula)
            if m1:
                m_points1 = m_points
                m_form1 = fuzzy_formula
                m_dict1 = dict(m_dict)
                (m_form1, m_dict1, m_points1) = \
                    _parse_chomp_and_rank(m1, m_form1, m_dict1, m_points1)
                if m_dict1:
                    #there was a real match
                    for match in \
                        Composition._comps_from_fuzzy_formula(m_form1,
                                                              m_dict1,
                                                              m_points1,
                                                              factor):
                        yield match

            #try to match two-letter elements
            m2 = re.match(r"([A-z]{2})([\.\d]*)", fuzzy_formula)
            if m2:
                m_points2 = m_points
                m_form2 = fuzzy_formula
                m_dict2 = dict(m_dict)
                (m_form2, m_dict2, m_points2) = \
                    _parse_chomp_and_rank(m2, m_form2, m_dict2, m_points2)
                if m_dict2:
                    #there was a real match
                    for match in \
                        Composition._comps_from_fuzzy_formula(m_form2, m_dict2,
                                                              m_points2,
                                                              factor):
                        yield match


def reduce_formula(sym_amt):
    """
    Helper method to reduce a sym_amt dict to a reduced formula and factor.

    Args:
        sym_amt (dict): {symbol: amount}.

    Returns:
        (reduced_formula, factor).
    """
    syms = sorted(sym_amt.keys(),
                  key=lambda s: get_el_sp(s).X)

    syms = filter(lambda s: sym_amt[s] > Composition.amount_tolerance, syms)
    num_el = len(syms)
    contains_polyanion = (num_el >= 3 and
                          get_el_sp(syms[num_el - 1]).X
                          - get_el_sp(syms[num_el - 2]).X < 1.65)

    factor = reduce(gcd, sym_amt.values())
    reduced_form = []
    n = num_el - 2 if contains_polyanion else num_el
    for i in range(0, n):
        s = syms[i]
        normamt = sym_amt[s] * 1.0 / factor
        reduced_form.append(s)
        reduced_form.append(formula_double_format(normamt))

    if contains_polyanion:
        poly_sym_amt = {syms[i]: sym_amt[syms[i]] / factor
                        for i in range(n, num_el)}
        (poly_form, poly_factor) = reduce_formula(poly_sym_amt)

        if poly_factor != 1:
            reduced_form.append("({}){}".format(poly_form, int(poly_factor)))
        else:
            reduced_form.append(poly_form)

    reduced_form = "".join(reduced_form)

    return reduced_form, factor


class CompositionError(Exception):
    """Exception class for composition errors"""
    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()
