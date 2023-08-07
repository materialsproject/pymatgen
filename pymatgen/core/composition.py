"""
This module implements a Composition class to represent compositions,
and a ChemicalPotential class to represent potentials.
"""

from __future__ import annotations

import collections
import os
import re
import string
import warnings
from functools import total_ordering
from itertools import combinations_with_replacement, product
from typing import Generator, Iterator, Union, cast

from monty.fractions import gcd, gcd_float
from monty.json import MSONable
from monty.serialization import loadfn

from pymatgen.core.periodic_table import DummySpecies, Element, Species, get_el_sp
from pymatgen.core.units import Mass
from pymatgen.util.string import Stringify, formula_double_format

SpeciesLike = Union[str, Element, Species, DummySpecies]


@total_ordering
class Composition(collections.abc.Hashable, collections.abc.Mapping, MSONable, Stringify):
    """
    Represents a Composition, which is essentially a {element:amount} mapping
    type. Composition is written to be immutable and hashable,
    unlike a standard Python dict.

    Note that the key can be either an Element or a Species. Elements and Species
    are treated differently. i.e., a Fe2+ is not the same as a Fe3+ Species and
    would be put in separate keys. This differentiation is deliberate to
    support using Composition to determine the fraction of a particular Species.

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
    special_formulas = dict(
        LiO="Li2O2",
        NaO="Na2O2",
        KO="K2O2",
        HO="H2O2",
        CsO="Cs2O2",
        RbO="Rb2O2",
        O="O2",
        N="N2",
        F="F2",
        Cl="Cl2",
        H="H2",
    )

    oxi_prob = None  # prior probability of oxidation used by oxi_state_guesses

    def __init__(self, *args, strict: bool = False, **kwargs) -> None:
        """
        Very flexible Composition construction, similar to the built-in Python
        dict(). Also extended to allow simple string init.

        Takes any inputs supported by the Python built-in dict function.

        1. A dict of either {Element/Species: amount},

            {string symbol:amount}, or {atomic number:amount} or any mixture
            of these. E.g., {Element("Li"): 2, Element("O"): 1},
            {"Li":2, "O":1}, {3: 2, 8: 1} all result in a Li2O composition.
        2. Keyword arg initialization, similar to a dict, e.g.,

            Composition(Li = 2, O = 1)

        In addition, the Composition constructor also allows a single
        string as an input formula. E.g., Composition("Li2O").

        Args:
            *args: Any number of 2-tuples as key-value pairs.
            strict (bool): Only allow valid Elements and Species in the Composition. Defaults to False.
            allow_negative (bool): Whether to allow negative compositions. Defaults to False.
            **kwargs: Additional kwargs supported by the dict() constructor.
        """
        # allow_negative must be popped from **kwargs due to *args ambiguity
        self.allow_negative = kwargs.pop("allow_negative", False)
        # it's much faster to recognize a composition and use the el_map than
        # to pass the composition to {}
        if len(args) == 1 and isinstance(args[0], Composition):
            elem_map = args[0]
        elif len(args) == 1 and isinstance(args[0], str):
            elem_map = self._parse_formula(args[0])  # type: ignore
        else:
            elem_map = dict(*args, **kwargs)  # type: ignore
        elem_amt = {}
        self._natoms = 0
        for key, val in elem_map.items():
            if val < -Composition.amount_tolerance and not self.allow_negative:
                raise ValueError("Amounts in Composition cannot be negative!")
            if abs(val) >= Composition.amount_tolerance:
                elem_amt[get_el_sp(key)] = val
                self._natoms += abs(val)
        self._data = elem_amt
        if strict and not self.valid:
            raise ValueError(f"Composition is not valid, contains: {', '.join(map(str, self.elements))}")

    def __getitem__(self, key: SpeciesLike) -> float:
        try:
            sp = get_el_sp(key)
            return self._data.get(sp, 0)
        except ValueError as exc:
            raise KeyError(f"Invalid {key=}") from exc

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self) -> Iterator[Species | Element | DummySpecies]:
        return self._data.__iter__()

    def __contains__(self, key) -> bool:
        try:
            sp = get_el_sp(key)
            if isinstance(sp, Species):
                return sp in self._data
            # Element or str
            return any(sp.symbol == s.symbol for s in self._data)
        except ValueError as exc:
            raise TypeError(f"Invalid {key=} for Composition") from exc

    def __eq__(self, other: object) -> bool:
        """Defines == for Compositions."""
        if not isinstance(other, (Composition, dict)):
            return NotImplemented

        # elements with amounts < Composition.amount_tolerance don't show up
        # in the el_map, so checking len enables us to only check one
        # composition's elements
        if len(self) != len(other):
            return False

        return all(abs(amt - other[el]) <= Composition.amount_tolerance for el, amt in self.items())

    def __ge__(self, other: object) -> bool:
        """
        Defines >= for Compositions. Should ONLY be used for defining a sort
        order (the behavior is probably not what you'd expect).
        """
        if not isinstance(other, Composition):
            return NotImplemented

        for el in sorted(set(self.elements + other.elements)):
            if other[el] - self[el] >= Composition.amount_tolerance:
                return False
            if self[el] - other[el] >= Composition.amount_tolerance:
                return True
        return True

    def __add__(self, other: object) -> Composition:
        """
        Adds two compositions. For example, an Fe2O3 composition + an FeO
        composition gives a Fe3O4 composition.
        """
        if not isinstance(other, (Composition, dict)):
            return NotImplemented

        new_el_map: dict[SpeciesLike, float] = collections.defaultdict(float)
        new_el_map.update(self)
        for k, v in other.items():
            new_el_map[get_el_sp(k)] += v
        return Composition(new_el_map, allow_negative=self.allow_negative)

    def __sub__(self, other: object) -> Composition:
        """
        Subtracts two compositions. For example, an Fe2O3 composition - an FeO
        composition gives an FeO2 composition.

        Raises:
            ValueError if the subtracted composition is greater than the
            original composition in any of its elements, unless allow_negative
            is True
        """
        if not isinstance(other, (Composition, dict)):
            return NotImplemented

        new_el_map: dict[SpeciesLike, float] = collections.defaultdict(float)
        new_el_map.update(self)
        for k, v in other.items():
            new_el_map[get_el_sp(k)] -= v
        return Composition(new_el_map, allow_negative=self.allow_negative)

    def __mul__(self, other: object) -> Composition:
        """
        Multiply a Composition by an integer or a float.
        Fe2O3 * 4 -> Fe8O12.
        """
        if not isinstance(other, (int, float)):
            return NotImplemented
        return Composition({el: self[el] * other for el in self}, allow_negative=self.allow_negative)

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> Composition:
        if not isinstance(other, (int, float)):
            return NotImplemented
        return Composition({el: self[el] / other for el in self}, allow_negative=self.allow_negative)

    __div__ = __truediv__

    def __hash__(self) -> int:
        """Hash based on the chemical system."""
        return hash(frozenset(self._data))

    @property
    def average_electroneg(self) -> float:
        """Average electronegativity of the composition."""
        return sum((el.X * abs(amt) for el, amt in self.items())) / self.num_atoms

    @property
    def total_electrons(self) -> float:
        """Total number of electrons in composition."""
        return sum((el.Z * abs(amt) for el, amt in self.items()))

    def almost_equals(self, other: Composition, rtol: float = 0.1, atol: float = 1e-8) -> bool:
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
        """True if composition is an element."""
        return len(self) == 1

    def copy(self) -> Composition:
        """A copy of the composition."""
        return Composition(self, allow_negative=self.allow_negative)

    @property
    def formula(self) -> str:
        """
        Returns a formula string, with elements sorted by electronegativity,
        e.g., Li4 Fe4 P4 O16.
        """
        sym_amt = self.get_el_amt_dict()
        syms = sorted(sym_amt, key=lambda sym: get_el_sp(sym).X)
        formula = [f"{s}{formula_double_format(sym_amt[s], False)}" for s in syms]
        return " ".join(formula)

    @property
    def alphabetical_formula(self) -> str:
        """
        Returns a formula string, with elements sorted by alphabetically
        e.g., Fe4 Li4 O16 P4.
        """
        return " ".join(sorted(self.formula.split()))

    @property
    def iupac_formula(self) -> str:
        """
        Returns a formula string, with elements sorted by the iupac
        electronegativity ordering defined in Table VI of "Nomenclature of
        Inorganic Chemistry (IUPAC Recommendations 2005)". This ordering
        effectively follows the groups and rows of the periodic table, except
        the Lanthanides, Actinides and hydrogen. Polyanions are still determined
        based on the true electronegativity of the elements.
        e.g. CH2(SO4)2.
        """
        sym_amt = self.get_el_amt_dict()
        syms = sorted(sym_amt, key=lambda s: get_el_sp(s).iupac_ordering)
        formula = [f"{s}{formula_double_format(sym_amt[s], False)}" for s in syms]
        return " ".join(formula)

    @property
    def element_composition(self) -> Composition:
        """
        Returns the composition replacing any species by the corresponding
        element.
        """
        return Composition(self.get_el_amt_dict(), allow_negative=self.allow_negative)

    @property
    def fractional_composition(self) -> Composition:
        """
        Returns the normalized composition in which the amounts of each species sum to
        1.
        E.g. "Fe2 O3".fractional_composition = "Fe0.4 O0.6".
        """
        return self / self._natoms

    @property
    def reduced_composition(self) -> Composition:
        """
        Returns the reduced composition, i.e. amounts normalized by greatest common denominator.
        E.g. "Fe4 P4 O16".reduced_composition = "Fe P O4".
        """
        return self.get_reduced_composition_and_factor()[0]

    def get_reduced_composition_and_factor(self) -> tuple[Composition, float]:
        """
        Calculates a reduced composition and factor.

        Returns:
            A normalized composition and a multiplicative factor, i.e.,
            Li4Fe4P4O16 returns (Composition("LiFePO4"), 4).
        """
        factor = self.get_reduced_formula_and_factor()[1]
        return self / factor, factor

    def get_reduced_formula_and_factor(self, iupac_ordering: bool = False) -> tuple[str, float]:
        """
        Calculates a reduced formula and factor.

        Args:
            iupac_ordering (bool, optional): Whether to order the
                formula by the iupac "electronegativity" series, defined in
                Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
                Recommendations 2005)". This ordering effectively follows
                the groups and rows of the periodic table, except the
                Lanthanides, Actinides and hydrogen. Note that polyanions
                will still be determined based on the true electronegativity of
                the elements.

        Returns:
            A pretty normalized formula and a multiplicative factor, i.e.,
            Li4Fe4P4O16 returns (LiFePO4, 4).
        """
        all_int = all(abs(x - round(x)) < Composition.amount_tolerance for x in self.values())
        if not all_int:
            return self.formula.replace(" ", ""), 1
        d = {k: int(round(v)) for k, v in self.get_el_amt_dict().items()}
        (formula, factor) = reduce_formula(d, iupac_ordering=iupac_ordering)

        if formula in Composition.special_formulas:
            formula = Composition.special_formulas[formula]
            factor /= 2

        return formula, factor

    def get_integer_formula_and_factor(
        self, max_denominator: int = 10000, iupac_ordering: bool = False
    ) -> tuple[str, float]:
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
                Lanthanides, Actinides and hydrogen. Note that polyanions
                will still be determined based on the true electronegativity of
                the elements.

        Returns:
            A pretty normalized formula and a multiplicative factor, i.e.,
            Li0.5O0.25 returns (Li2O, 0.25). O0.25 returns (O2, 0.125)
        """
        el_amt = self.get_el_amt_dict()
        gcd = gcd_float(list(el_amt.values()), 1 / max_denominator)

        dct = {k: round(v / gcd) for k, v in el_amt.items()}
        formula, factor = reduce_formula(dct, iupac_ordering=iupac_ordering)
        if formula in Composition.special_formulas:
            formula = Composition.special_formulas[formula]
            factor /= 2
        return formula, factor * gcd

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
        elem_comp = self.element_composition
        elements = sorted(el.symbol for el in elem_comp)
        hill_elements = []
        if "C" in elements:
            hill_elements.append("C")
            elements.remove("C")
            if "H" in elements:
                hill_elements.append("H")
                elements.remove("H")
        hill_elements += elements

        formula = [f"{el}{formula_double_format(elem_comp[el]) if elem_comp[el] != 1 else ''}" for el in hill_elements]
        return " ".join(formula)

    @property
    def elements(self) -> list[Element | Species | DummySpecies]:
        """Returns list of elements in Composition."""
        return list(self)

    def __str__(self):
        return " ".join(f"{k}{formula_double_format(v, ignore_ones=False)}" for k, v in self.as_dict().items())

    def to_pretty_string(self) -> str:
        """
        Returns:
            str: Same as output __str__() but without spaces.
        """
        return re.sub(r"\s+", "", str(self))

    @property
    def num_atoms(self) -> float:
        """
        Total number of atoms in Composition. For negative amounts, sum
        of absolute values.
        """
        return self._natoms

    @property
    def weight(self) -> float:
        """Total molecular weight of Composition."""
        return Mass(sum(amount * el.atomic_mass for el, amount in self.items()), "amu")

    def get_atomic_fraction(self, el: SpeciesLike) -> float:
        """
        Calculate atomic fraction of an Element or Species.

        Args:
            el (Element/Species): Element or Species to get fraction for.

        Returns:
            Atomic fraction for element el in Composition
        """
        return abs(self[el]) / self._natoms

    def get_wt_fraction(self, el: SpeciesLike) -> float:
        """
        Calculate weight fraction of an Element or Species.

        Args:
            el (Element | Species): Element or Species to get fraction for.

        Returns:
            float: Weight fraction for element el in Composition.
        """
        el_mass = cast(float, get_el_sp(el).atomic_mass)
        return el_mass * abs(self[el]) / self.weight

    def contains_element_type(self, category: str) -> bool:
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
        allowed_categories = (
            "noble_gas",
            "transition_metal",
            "post_transition_metal",
            "rare_earth_metal",
            "metal",
            "metalloid",
            "alkali",
            "alkaline",
            "halogen",
            "chalcogen",
            "lanthanoid",
            "actinoid",
            "quadrupolar",
            "s-block",
            "p-block",
            "d-block",
            "f-block",
        )

        if category not in allowed_categories:
            raise ValueError(f"Please pick a category from: {allowed_categories}")

        if "block" in category:
            return any(category[0] in el.block for el in self.elements)
        return any(getattr(el, f"is_{category}") for el in self.elements)

    def _parse_formula(self, formula: str) -> dict[str, float]:
        """
        Args:
            formula (str): A string formula, e.g. Fe2O3, Li3Fe2(PO4)3.

        Returns:
            Composition with that formula.

        Notes:
            In the case of Metallofullerene formula (e.g. Y3N@C80),
            the @ mark will be dropped and passed to parser.
        """
        # for Metallofullerene like "Y3N@C80"
        formula = formula.replace("@", "")

        def get_sym_dict(form: str, factor: int | float) -> dict[str, float]:
            sym_dict: dict[str, float] = collections.defaultdict(float)
            for m in re.finditer(r"([A-Z][a-z]*)\s*([-*\.e\d]*)", form):
                el = m.group(1)
                amt = 1.0
                if m.group(2).strip() != "":
                    amt = float(m.group(2))
                sym_dict[el] += amt * factor
                form = form.replace(m.group(), "", 1)
            if form.strip():
                raise ValueError(f"{form} is an invalid formula!")
            return sym_dict

        m = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
        if m:
            factor = 1.0
            if m.group(2) != "":
                factor = float(m.group(2))
            unit_sym_dict = get_sym_dict(m.group(1), factor)
            expanded_sym = "".join(f"{el}{amt}" for el, amt in unit_sym_dict.items())
            expanded_formula = formula.replace(m.group(), expanded_sym)
            return self._parse_formula(expanded_formula)
        return get_sym_dict(formula, 1)

    @property
    def anonymized_formula(self) -> str:
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
            anon += f"{e}{amt_str}"
        return anon

    @property
    def chemical_system(self) -> str:
        """
        Get the chemical system of a Composition, for example "O-Si" for
        SiO2. Chemical system is a string of a list of elements
        sorted alphabetically and joined by dashes, by convention for use
        in database keys.
        """
        return "-".join(sorted(el.symbol for el in self.elements))

    @property
    def valid(self) -> bool:
        """
        Returns True if Composition contains valid elements or species and
        False if the Composition contains any dummy species.
        """
        return not any(isinstance(el, DummySpecies) for el in self.elements)

    def __repr__(self):
        formula = " ".join(f"{k}{':' if hasattr(k, 'oxi_state') else ''}{v:g}" for k, v in self.items())
        cls_name = type(self).__name__
        return f"{cls_name}({formula!r})"

    @classmethod
    def from_dict(cls, d) -> Composition:
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

    @classmethod
    def from_weight_dict(cls, weight_dict) -> Composition:
        """
        Creates a Composition based on a dict of atomic fractions calculated
        from a dict of weight fractions. Allows for quick creation of the class
        from weight-based notations commonly used in the industry, such as
        Ti6V4Al and Ni60Ti40.

        Args:
            weight_dict (dict): {symbol: weight_fraction} dict.

        Returns:
            Composition
        """
        weight_sum = sum(val / Element(el).atomic_mass for el, val in weight_dict.items())
        comp_dict = {el: val / Element(el).atomic_mass / weight_sum for el, val in weight_dict.items()}

        return cls(comp_dict)

    def get_el_amt_dict(self) -> dict[str, float]:
        """
        Returns:
            dict[str, float]: element symbol and (unreduced) amount. E.g.
                {"Fe": 4.0, "O":6.0} or {"Fe3+": 4.0, "O2-":6.0}.
        """
        dic: dict[str, float] = collections.defaultdict(float)
        for el, amt in self.items():
            dic[el.symbol] += amt
        return dic

    def as_dict(self) -> dict[str, float]:
        """
        Note: Subtly different from get_el_amt_dict in that they keys here are str(Element) instead of Element.symbol.

        Returns:
            dict[str, float]: element symbol and (unreduced) amount. E.g.
                {"Fe": 4.0, "O":6.0} or {"Fe3+": 4.0, "O2-":6.0}
        """
        dic: dict[str, float] = collections.defaultdict(float)
        for el, amt in self.items():
            dic[str(el)] += amt
        return dic

    @property
    def to_reduced_dict(self) -> dict[str, float]:
        """
        Returns:
            dict[str, float]: element symbols mapped to reduced amount e.g. {"Fe": 2.0, "O":3.0}.
        """
        return self.reduced_composition.as_dict()

    @property
    def to_weight_dict(self) -> dict[str, float]:
        """
        Returns:
            dict[str, float] with weight fraction of each component {"Ti": 0.90, "V": 0.06, "Al": 0.04}.
        """
        return {str(el): self.get_wt_fraction(el) for el in self.elements}

    @property
    def to_data_dict(self) -> dict:
        """
        Returns:
            A dict with many keys and values relating to Composition/Formula,
            including reduced_cell_composition, unit_cell_composition,
            reduced_cell_formula, elements and nelements.
        """
        return {
            "reduced_cell_composition": self.reduced_composition,
            "unit_cell_composition": self.as_dict(),
            "reduced_cell_formula": self.reduced_formula,
            "elements": list(map(str, self)),
            "nelements": len(self),
        }

    def oxi_state_guesses(
        self,
        oxi_states_override: dict | None = None,
        target_charge: float = 0,
        all_oxi_states: bool = False,
        max_sites: int | None = None,
    ) -> list[dict[str, float]]:
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
        return self._get_oxi_state_guesses(all_oxi_states, max_sites, oxi_states_override, target_charge)[0]

    def replace(self, elem_map: dict[str, str | dict[str, int | float]]) -> Composition:
        """
        Replace elements in a composition. Returns a new Composition, leaving the old one unchanged.

        Args:
            elem_map (dict[str, str | dict[str, int | float]]): dict of elements or species to swap. E.g.
                {"Li": "Na"} performs a Li for Na substitution. The target can be a {species: factor} dict. For
                example, in Fe2O3 you could map {"Fe": {"Mg": 0.5, "Cu":0.5}} to obtain MgCuO3.

        Returns:
            Composition: New object with elements remapped according to elem_map.
        """
        # drop inapplicable substitutions
        invalid_elems = [key for key in elem_map if key not in self]
        if invalid_elems:
            warnings.warn(
                "Some elements to be substituted are not present in composition. Please check your input. "
                f"Problematic element = {invalid_elems}; {self}"
            )
        for elem in invalid_elems:
            elem_map.pop(elem)

        # start with elements that remain unchanged (not in elem_map)
        new_comp = {elem: amount for elem, amount in self.as_dict().items() if elem not in elem_map}

        for old_elem, new_elem in elem_map.items():
            amount = self[old_elem]

            # build a dictionary of substitutions to be made
            subs = {}
            if isinstance(new_elem, dict):
                for el, factor in new_elem.items():
                    subs[el] = factor * amount
            else:
                subs = {new_elem: amount}

            # and apply the substitutions to the new composition
            for el, amt in subs.items():
                if el in new_comp:
                    new_comp[el] += amt
                else:
                    new_comp[el] = amt

                # check for ambiguous input (see issue #2553)
                if el in self:
                    warnings.warn(
                        f"Same element ({el}) in both the keys and values of the substitution!"
                        "This can be ambiguous, so be sure to check your result."
                    )

        return Composition(new_comp)

    def add_charges_from_oxi_state_guesses(
        self,
        oxi_states_override: dict | None = None,
        target_charge: float = 0,
        all_oxi_states: bool = False,
        max_sites: int | None = None,
    ) -> Composition:
        """
        Assign oxidation states based on guessed oxidation states.

        See `oxi_state_guesses` for an explanation of how oxidation states are
        guessed. This operation uses the set of oxidation states for each site
        that were determined to be most likely from the oxidation state guessing
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
        _, oxidation_states = self._get_oxi_state_guesses(all_oxi_states, max_sites, oxi_states_override, target_charge)

        # Special case: No charged compound is possible
        if not oxidation_states:
            return Composition({Species(e, 0): f for e, f in self.items()})

        # Generate the species
        species = []
        for el, charges in oxidation_states[0].items():
            species.extend([Species(el, c) for c in charges])

        # Return the new object
        return Composition(collections.Counter(species))

    def remove_charges(self) -> Composition:
        """
        Returns a new Composition with charges from each Species removed.

        Returns:
            Composition object without charge decoration, for example
            {"Fe3+": 2.0, "O2-":3.0} becomes {"Fe": 2.0, "O":3.0}
        """
        dct: dict[Element, float] = collections.defaultdict(float)
        for specie, amt in self.items():
            dct[Element(specie.symbol)] += amt
        return Composition(dct)

    def _get_oxi_state_guesses(self, all_oxi_states, max_sites, oxi_states_override, target_charge):
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
            list[dict]: Each dict maps the element symbol to a list of
                oxidation states for each site of that element. For example, Fe3O4 could
                return a list of [2,2,2,3,3,3] for the oxidation states of the 6 Fe sites.
                If the composition is not charge balanced, an empty list is returned.
        """
        comp = self.copy()
        # reduce Composition if necessary
        if max_sites and max_sites < 0:
            comp = self.reduced_composition

            if max_sites < -1 and comp.num_atoms > abs(max_sites):
                raise ValueError(f"Composition {comp} cannot accommodate max_sites setting!")

        elif max_sites and comp.num_atoms > max_sites:
            reduced_comp, reduced_factor = self.get_reduced_composition_and_factor()
            if reduced_factor > 1:
                reduced_comp *= max(1, int(max_sites / reduced_comp.num_atoms))
                comp = reduced_comp  # as close to max_sites as possible
            if comp.num_atoms > max_sites:
                raise ValueError(f"Composition {comp} cannot accommodate max_sites setting!")

        # Load prior probabilities of oxidation states, used to rank solutions
        if not Composition.oxi_prob:
            module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
            all_data = loadfn(f"{module_dir}/../analysis/icsd_bv.yaml")
            Composition.oxi_prob = {Species.from_str(sp): data for sp, data in all_data["occurrence"].items()}
        oxi_states_override = oxi_states_override or {}
        # assert: Composition only has integer amounts
        if not all(amt == int(amt) for amt in comp.values()):
            raise ValueError("Charge balance analysis requires integer values in Composition!")

        # for each element, determine all possible sum of oxidations
        # (taking into account nsites for that particular element)
        el_amt = comp.get_el_amt_dict()
        elements = list(el_amt)
        el_sums = []  # matrix: dim1= el_idx, dim2=possible sums
        el_sum_scores = collections.defaultdict(set)  # dict of el_idx, sum -> score
        el_best_oxid_combo = {}  # dict of el_idx, sum -> oxid combo with best score
        for idx, el in enumerate(elements):
            el_sum_scores[idx] = {}
            el_best_oxid_combo[idx] = {}
            el_sums.append([])
            if oxi_states_override.get(el):
                oxids = oxi_states_override[el]
            elif all_oxi_states:
                oxids = Element(el).oxidation_states
            else:
                oxids = Element(el).icsd_oxidation_states or Element(el).oxidation_states

            # get all possible combinations of oxidation states
            # and sum each combination
            for oxid_combo in combinations_with_replacement(oxids, int(el_amt[el])):
                # List this sum as a possible option
                oxid_sum = sum(oxid_combo)
                if oxid_sum not in el_sums[idx]:
                    el_sums[idx].append(oxid_sum)

                # Determine how probable is this combo?
                score = sum(Composition.oxi_prob.get(Species(el, o), 0) for o in oxid_combo)

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
                el_sum_sol = dict(zip(elements, x))  # element->oxid_sum
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
                all_oxid_combo.append({e: el_best_oxid_combo[idx][v] for idx, (e, v) in enumerate(zip(elements, x))})

        # sort the solutions by highest to lowest score
        if all_scores:
            all_sols, all_oxid_combo = zip(
                *(
                    (y, x)
                    for (z, y, x) in sorted(
                        zip(all_scores, all_sols, all_oxid_combo),
                        key=lambda pair: pair[0],
                        reverse=True,
                    )
                )
            )
        return all_sols, all_oxid_combo

    @staticmethod
    def ranked_compositions_from_indeterminate_formula(
        fuzzy_formula: str, lock_if_strict: bool = True
    ) -> list[Composition]:
        """
        Takes in a formula where capitalization might not be correctly entered,
        and suggests a ranked list of potential Composition matches.
        Author: Anubhav Jain.

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
            except ValueError:
                pass

        all_matches = Composition._comps_from_fuzzy_formula(fuzzy_formula)
        # remove duplicates
        uniq_matches = list(set(all_matches))
        # sort matches by rank descending
        ranked_matches = sorted(uniq_matches, key=lambda match: (match[1], match[0]), reverse=True)

        return [m[0] for m in ranked_matches]

    @staticmethod
    def _comps_from_fuzzy_formula(
        fuzzy_formula: str,
        m_dict: dict[str, float] | None = None,
        m_points: int = 0,
        factor: int | float = 1,
    ) -> Generator[tuple[Composition, int], None, None]:
        """
        A recursive helper method for formula parsing that helps in
        interpreting and ranking indeterminate formulas.
        Author: Anubhav Jain.

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
            list[tuple[Composition, int]]: A list of tuples, with the first element being a Composition
                and the second element being the number of points awarded that Composition interpretation.
        """
        m_dict = m_dict or {}

        def _parse_chomp_and_rank(m, f, m_dict, m_points):
            """
            A helper method for formula parsing that helps in interpreting and
            ranking indeterminate formulas
            Author: Anubhav Jain.

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
                raise ValueError("Invalid element symbol entered!")
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
                for match in Composition._comps_from_fuzzy_formula(mp.group(1), mp_dict, mp_points, factor=mp_factor):
                    only_me = True
                    # Match the stuff outside the parentheses and return the
                    # sum.

                    for match2 in Composition._comps_from_fuzzy_formula(mp_form, mp_dict, mp_points, factor=1):
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
                (m_form1, m_dict1, m_points1) = _parse_chomp_and_rank(m1, m_form1, m_dict1, m_points1)
                if m_dict1:
                    # there was a real match
                    for match in Composition._comps_from_fuzzy_formula(m_form1, m_dict1, m_points1, factor):
                        yield match

            # try to match two-letter elements
            m2 = re.match(r"([A-z]{2})([\.\d]*)", fuzzy_formula)
            if m2:
                m_points2 = m_points
                m_form2 = fuzzy_formula
                m_dict2 = dict(m_dict)
                (m_form2, m_dict2, m_points2) = _parse_chomp_and_rank(m2, m_form2, m_dict2, m_points2)
                if m_dict2:
                    # there was a real match
                    for match in Composition._comps_from_fuzzy_formula(m_form2, m_dict2, m_points2, factor):
                        yield match


def reduce_formula(sym_amt, iupac_ordering: bool = False) -> tuple[str, float]:
    """
    Helper method to reduce a sym_amt dict to a reduced formula and factor.

    Args:
        sym_amt (dict): {symbol: amount}.
        iupac_ordering (bool, optional): Whether to order the
            formula by the iupac "electronegativity" series, defined in
            Table VI of "Nomenclature of Inorganic Chemistry (IUPAC
            Recommendations 2005)". This ordering effectively follows
            the groups and rows of the periodic table, except the
            Lanthanides, Actinides and hydrogen. Note that polyanions
            will still be determined based on the true electronegativity of
            the elements.

    Returns:
        (reduced_formula, factor).
    """
    syms = sorted(sym_amt, key=lambda x: [get_el_sp(x).X, x])

    syms = list(filter(lambda x: abs(sym_amt[x]) > Composition.amount_tolerance, syms))

    factor = 1
    # Enforce integers for doing gcd.
    if all(int(i) == i for i in sym_amt.values()):
        factor = abs(gcd(*(int(i) for i in sym_amt.values())))

    polyanion = []
    # if the composition contains a poly anion
    if len(syms) >= 3 and get_el_sp(syms[-1]).X - get_el_sp(syms[-2]).X < 1.65:
        poly_sym_amt = {syms[i]: sym_amt[syms[i]] / factor for i in [-2, -1]}
        (poly_form, poly_factor) = reduce_formula(poly_sym_amt, iupac_ordering=iupac_ordering)

        if poly_factor != 1:
            polyanion.append(f"({poly_form}){poly_factor}")

    syms = syms[: len(syms) - 2 if polyanion else len(syms)]

    if iupac_ordering:
        syms = sorted(syms, key=lambda x: [get_el_sp(x).iupac_ordering, x])

    reduced_form = []
    for sym in syms:
        norm_amt = sym_amt[sym] * 1.0 / factor
        reduced_form.append(sym)
        reduced_form.append(str(formula_double_format(norm_amt)))

    return "".join([*reduced_form, *polyanion]), factor


class ChemicalPotential(dict, MSONable):
    """
    Class to represent set of chemical potentials. Can be: multiplied/divided by a Number
    multiplied by a Composition (returns an energy) added/subtracted with other ChemicalPotentials.
    """

    def __init__(self, *args, **kwargs):
        """
        Args:
            *args: any valid dict init arguments
            **kwargs: any valid dict init arguments.
        """
        dct = dict(*args, **kwargs)
        super().__init__((get_el_sp(k), v) for k, v in dct.items())
        if len(dct) != len(self):
            raise ValueError("Duplicate potential specified")

    def __mul__(self, other: object) -> ChemicalPotential:
        if isinstance(other, (int, float)):
            return ChemicalPotential({k: v * other for k, v in self.items()})
        return NotImplemented

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> ChemicalPotential:
        if isinstance(other, (int, float)):
            return ChemicalPotential({k: v / other for k, v in self.items()})
        return NotImplemented

    __div__ = __truediv__

    def __sub__(self, other: object) -> ChemicalPotential:
        if isinstance(other, ChemicalPotential):
            els = {*self} | {*other}
            return ChemicalPotential({e: self.get(e, 0) - other.get(e, 0) for e in els})
        return NotImplemented

    def __add__(self, other: object) -> ChemicalPotential:
        if isinstance(other, ChemicalPotential):
            els = {*self} | {*other}
            return ChemicalPotential({e: self.get(e, 0) + other.get(e, 0) for e in els})
        return NotImplemented

    def get_energy(self, composition: Composition, strict: bool = True) -> float:
        """
        Calculates the energy of a composition.

        Args:
            composition (Composition): input composition
            strict (bool): Whether all potentials must be specified
        """
        if strict and set(composition) > set(self):
            s = set(composition) - set(self)
            raise ValueError(f"Potentials not specified for {s}")
        return sum(self.get(k, 0) * v for k, v in composition.items())

    def __repr__(self):
        return "ChemPots: " + super().__repr__()


class CompositionError(Exception):
    """Exception class for composition errors."""
