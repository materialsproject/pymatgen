"""This module implements a Composition class to represent compositions,
and a ChemicalPotential class to represent potentials.
"""

from __future__ import annotations

import collections
import os
import re
import string
import warnings
from collections import defaultdict
from functools import total_ordering
from itertools import combinations_with_replacement, product
from math import isnan
from typing import TYPE_CHECKING, cast

from monty.fractions import gcd, gcd_float
from monty.json import MSONable
from monty.serialization import loadfn
from pymatgen.core.periodic_table import DummySpecies, Element, ElementType, Species, get_el_sp
from pymatgen.core.units import Mass
from pymatgen.util.string import Stringify, formula_double_format

if TYPE_CHECKING:
    from collections.abc import Generator, Iterator
    from typing import Any, ClassVar

    from pymatgen.util.typing import SpeciesLike
    from typing_extensions import Self

module_dir = os.path.dirname(os.path.abspath(__file__))


@total_ordering
class Composition(collections.abc.Hashable, collections.abc.Mapping, MSONable, Stringify):
    """Represents a Composition, which is essentially a {element:amount} mapping
    type. Composition is written to be immutable and hashable,
    unlike a standard Python dict.

    Note that the key can be either an Element or a Species. Elements and Species
    are treated differently. i.e., a Fe2+ is not the same as a Fe3+ Species and
    would be put in separate keys. This differentiation is deliberate to
    support using Composition to determine the fraction of a particular Species.

    Works almost completely like a standard python dictionary, except that
    __getitem__ is overridden to return 0 when an element is not found.
    (somewhat like a defaultdict, except it is immutable).

    Also adds more convenience methods relevant to compositions, e.g.
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
    charge_balanced_tolerance = 1e-8

    # Special formula handling for peroxides and certain elements. This is so
    # that formula output does not write LiO instead of Li2O2 for example.
    special_formulas: ClassVar[dict[str, str]] = {
        "LiO": "Li2O2",
        "NaO": "Na2O2",
        "KO": "K2O2",
        "HO": "H2O2",
        "CsO": "Cs2O2",
        "RbO": "Rb2O2",
        "O": "O2",
        "N": "N2",
        "F": "F2",
        "Cl": "Cl2",
        "H": "H2",
    }

    oxi_prob = None  # prior probability of oxidation used by oxi_state_guesses

    def __init__(self, *args, strict: bool = False, **kwargs) -> None:
        """Very flexible Composition construction, similar to the built-in Python
        dict(). Also extended to allow simple string init.

        Takes any inputs supported by the Python built-in dict function.

        1. A dict of either {Element/Species: amount},

            {string symbol:amount}, or {atomic number:amount} or any mixture
            of these. e.g. {Element("Li"): 2, Element("O"): 1},
            {"Li":2, "O":1}, {3: 2, 8: 1} all result in a Li2O composition.
        2. Keyword arg initialization, similar to a dict, e.g.

            Composition(Li = 2, O = 1)

        In addition, the Composition constructor also allows a single
        string as an input formula. e.g. Composition("Li2O").

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
        if len(args) == 1 and isinstance(args[0], type(self)):
            elem_map = args[0]
        elif len(args) == 1 and isinstance(args[0], str):
            elem_map = self._parse_formula(args[0])  # type: ignore[assignment]
        elif len(args) == 1 and isinstance(args[0], float) and isnan(args[0]):
            raise ValueError("float('NaN') is not a valid Composition, did you mean 'NaN'?")
        else:
            elem_map = dict(*args, **kwargs)  # type: ignore[assignment]
        elem_amt = {}
        self._n_atoms = 0
        for key, val in elem_map.items():
            if val < -type(self).amount_tolerance and not self.allow_negative:
                raise ValueError("Amounts in Composition cannot be negative!")
            if abs(val) >= type(self).amount_tolerance:
                elem_amt[get_el_sp(key)] = val
                self._n_atoms += abs(val)
        self._data = elem_amt
        if strict and not self.valid:
            raise ValueError(f"Composition is not valid, contains: {', '.join(map(str, self.elements))}")

    def __getitem__(self, key: SpeciesLike) -> float:
        try:
            sp = get_el_sp(key)
            if isinstance(sp, Species):
                return self._data.get(sp, 0)
            # sp is Element or str
            return sum(
                val for key, val in self._data.items() if getattr(key, "symbol", key) == getattr(sp, "symbol", sp)
            )
        except ValueError as exc:
            raise KeyError(f"Invalid {key=}") from exc

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self) -> Iterator[Species | Element | DummySpecies]:
        return iter(self._data)

    def __contains__(self, key) -> bool:
        try:
            sp = get_el_sp(key)
            if isinstance(sp, Species):
                return sp in self._data
            # key is Element or str
            return any(sp.symbol == s.symbol for s in self._data)
        except ValueError as exc:
            raise TypeError(f"Invalid {key=} for Composition") from exc

    def __eq__(self, other: object) -> bool:
        """Composition equality. We consider compositions equal if they have the
        same elements and the amounts are within Composition.amount_tolerance
        of each other.

        Args:
            other: Composition to compare to.
        """
        if not isinstance(other, (type(self), dict)):
            return NotImplemented

        # elements with amounts < Composition.amount_tolerance don't show up
        # in the el_map, so checking len enables us to only check one
        # composition's elements
        if len(self) != len(other):
            return False

        return all(abs(amt - other[el]) <= type(self).amount_tolerance for el, amt in self.items())

    def __ge__(self, other: object) -> bool:
        """Composition greater than or equal to. We sort the elements in the compositions in order of
        electronegativity. The amount of the most electropositive element that is not equal within a certain
        tolerance factor is used to make a comparison. Note that an element not present in a Composition has an implied
        amount of 0.

        Should ONLY be used for defining a sort order (the behavior is probably not what you'd expect).
        """
        if not isinstance(other, type(self)):
            return NotImplemented

        for el in sorted(set(self.elements + other.elements)):
            if other[el] - self[el] >= type(self).amount_tolerance:
                return False
            # TODO @janosh 2024-04-29: is this a bug? why would we return True early?
            if self[el] - other[el] >= type(self).amount_tolerance:
                return True
        return True

    def __add__(self, other: object) -> Self:
        """Add two compositions. For example, an Fe2O3 composition + an FeO
        composition gives a Fe3O4 composition.
        """
        if not isinstance(other, (type(self), dict)):
            return NotImplemented

        new_el_map: dict[SpeciesLike, float] = defaultdict(float)
        new_el_map.update(self)
        for key, val in other.items():
            new_el_map[get_el_sp(key)] += val
        return type(self)(new_el_map, allow_negative=self.allow_negative)

    def __sub__(self, other: object) -> Self:
        """Subtracts two compositions. For example, an Fe2O3 composition - an FeO
        composition gives an FeO2 composition.

        Raises:
            ValueError if the subtracted composition is greater than the
            original composition in any of its elements, unless allow_negative
            is True
        """
        if not isinstance(other, (type(self), dict)):
            return NotImplemented

        new_el_map: dict[SpeciesLike, float] = defaultdict(float)
        new_el_map.update(self)
        for key, val in other.items():
            new_el_map[get_el_sp(key)] -= val
        return type(self)(new_el_map, allow_negative=self.allow_negative)

    def __mul__(self, other: object) -> Self:
        """Multiply a Composition by an integer or a float.
        Fe2O3 * 4 -> Fe8O12.
        """
        if not isinstance(other, (int, float)):
            return NotImplemented
        return type(self)({el: self[el] * other for el in self}, allow_negative=self.allow_negative)

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> Self:
        if not isinstance(other, (int, float)):
            return NotImplemented
        return type(self)({el: self[el] / other for el in self}, allow_negative=self.allow_negative)

    __div__ = __truediv__

    def __hash__(self) -> int:
        """Hash based on the chemical system."""
        return hash(frozenset(self._data))

    def __repr__(self) -> str:
        formula = " ".join(f"{key}{':' if hasattr(key, 'oxi_state') else ''}{val:g}" for key, val in self.items())
        cls_name = type(self).__name__
        return f"{cls_name}({formula!r})"

    def __str__(self) -> str:
        return " ".join(f"{key}{formula_double_format(val, ignore_ones=False)}" for key, val in self.as_dict().items())

    def to_pretty_string(self) -> str:
        """
        Returns:
            str: Same output as __str__() but without spaces.
        """
        return re.sub(r"\s+", "", str(self))

    @property
    def average_electroneg(self) -> float:
        """Average electronegativity of the composition."""
        return sum((el.X * abs(amt) for el, amt in self.items())) / self.num_atoms

    @property
    def total_electrons(self) -> float:
        """Total number of electrons in composition."""
        return sum((el.Z * abs(amt) for el, amt in self.items()))

    def almost_equals(
        self,
        other: Composition,
        rtol: float = 0.1,
        atol: float = 1e-8,
    ) -> bool:
        """Get true if compositions are equal within a tolerance.

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

    def copy(self) -> Self:
        """A copy of the composition."""
        return type(self)(self, allow_negative=self.allow_negative)

    @property
    def formula(self) -> str:
        """A formula string, with elements sorted by electronegativity,
        e.g. Li4 Fe4 P4 O16.
        """
        sym_amt = self.get_el_amt_dict()
        syms = sorted(sym_amt, key=lambda sym: get_el_sp(sym).X)
        formula = [f"{s}{formula_double_format(sym_amt[s], ignore_ones= False)}" for s in syms]
        return " ".join(formula)

    @property
    def alphabetical_formula(self) -> str:
        """A formula string, with elements sorted by alphabetically
        e.g. Fe4 Li4 O16 P4.
        """
        return " ".join(sorted(self.formula.split()))

    @property
    def iupac_formula(self) -> str:
        """A formula string, with elements sorted by the IUPAC
        electronegativity ordering defined in Table VI of "Nomenclature of
        Inorganic Chemistry (IUPAC Recommendations 2005)". This ordering
        effectively follows the groups and rows of the periodic table, except
        the Lanthanides, Actinides and hydrogen. Polyanions are still determined
        based on the true electronegativity of the elements.
        e.g. CH2(SO4)2.
        """
        sym_amt = self.get_el_amt_dict()
        syms = sorted(sym_amt, key=lambda s: get_el_sp(s).iupac_ordering)
        formula = [f"{s}{formula_double_format(sym_amt[s], ignore_ones= False)}" for s in syms]
        return " ".join(formula)

    @property
    def element_composition(self) -> Self:
        """The composition replacing any species by the corresponding element."""
        return type(self)(self.get_el_amt_dict(), allow_negative=self.allow_negative)

    @property
    def fractional_composition(self) -> Self:
        """The normalized composition in which the amounts of each species sum to
        1.
        E.g. "Fe2 O3".fractional_composition = "Fe0.4 O0.6".
        """
        return self / self._n_atoms

    @property
    def reduced_composition(self) -> Self:
        """The reduced composition, i.e. amounts normalized by greatest common denominator.
        E.g. "Fe4 P4 O16".reduced_composition = "Fe P O4".
        """
        return self.get_reduced_composition_and_factor()[0]

    def get_reduced_composition_and_factor(self) -> tuple[Self, float]:
        """Calculate a reduced composition and factor.

        Returns:
            A normalized composition and a multiplicative factor, i.e.,
            Li4Fe4P4O16 returns (Composition("LiFePO4"), 4).
        """
        factor = self.get_reduced_formula_and_factor()[1]
        return self / factor, factor

    def get_reduced_formula_and_factor(self, iupac_ordering: bool = False) -> tuple[str, float]:
        """Calculate a reduced formula and factor.

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
        all_int = all(abs(val - round(val)) < type(self).amount_tolerance for val in self.values())
        if not all_int:
            return self.formula.replace(" ", ""), 1
        el_amt_dict = {key: int(round(val)) for key, val in self.get_el_amt_dict().items()}
        formula, factor = reduce_formula(el_amt_dict, iupac_ordering=iupac_ordering)

        if formula in type(self).special_formulas:
            formula = type(self).special_formulas[formula]
            factor /= 2

        return formula, factor

    def get_integer_formula_and_factor(
        self, max_denominator: int = 10000, iupac_ordering: bool = False
    ) -> tuple[str, float]:
        """Calculate an integer formula and factor.

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
        _gcd = gcd_float(list(el_amt.values()), 1 / max_denominator)

        dct = {key: round(val / _gcd) for key, val in el_amt.items()}
        formula, factor = reduce_formula(dct, iupac_ordering=iupac_ordering)
        if formula in type(self).special_formulas:
            formula = type(self).special_formulas[formula]
            factor /= 2
        return formula, factor * _gcd

    @property
    def reduced_formula(self) -> str:
        """A pretty normalized formula, i.e., LiFePO4 instead of
        Li4Fe4P4O16.
        """
        return self.get_reduced_formula_and_factor()[0]

    @property
    def hill_formula(self) -> str:
        """The Hill system (or Hill notation) is a system of writing empirical chemical
        formulas, molecular chemical formulas and components of a condensed formula such
        that the number of carbon atoms in a molecule is indicated first, the number of
        hydrogen atoms next, and then the number of all other chemical elements
        subsequently, in alphabetical order of the chemical symbols. When the formula
        contains no carbon, all the elements, including hydrogen, are listed
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
        """List of elements in Composition."""
        return list(self)

    @property
    def chemical_system_set(self) -> set[str]:
        """The set of elements in the Composition. E.g. {"O", "Si"} for SiO2."""
        return {el.symbol for el in self.elements}

    @property
    def chemical_system(self) -> str:
        """The chemical system of a Composition, for example "O-Si" for
        SiO2. Chemical system is a string of a list of elements
        sorted alphabetically and joined by dashes, by convention for use
        in database keys.
        """
        return "-".join(sorted(self.chemical_system_set))

    @property
    def num_atoms(self) -> float:
        """Total number of atoms in Composition. For negative amounts, sum
        of absolute values.
        """
        return self._n_atoms

    @property
    def weight(self) -> float:
        """Total molecular weight of Composition."""
        return Mass(sum(amount * el.atomic_mass for el, amount in self.items()), "amu")

    def get_atomic_fraction(self, el: SpeciesLike) -> float:
        """Calculate atomic fraction of an Element or Species.

        Args:
            el (Element/Species): Element or Species to get fraction for.

        Returns:
            Atomic fraction for element el in Composition
        """
        return abs(self[el]) / self._n_atoms

    def get_wt_fraction(self, el: SpeciesLike) -> float:
        """Calculate weight fraction of an Element or Species.

        Args:
            el (Element | Species): Element or Species to get fraction for.

        Returns:
            float: Weight fraction for element el in Composition.
        """
        el_mass = cast(float, get_el_sp(el).atomic_mass)
        return el_mass * abs(self[el]) / self.weight

    def contains_element_type(self, category: str) -> bool:
        """Check if Composition contains any elements matching a given category.

        Args:
            category (str): one of "noble_gas", "transition_metal",
                "post_transition_metal", "rare_earth_metal", "metal", "metalloid",
                "alkali", "alkaline", "halogen", "chalcogen", "lanthanoid",
                "actinoid", "radioactive", "quadrupolar", "s-block", "p-block", "d-block", "f-block".

        Returns:
            bool: Whether any elements in Composition match category.
        """
        allowed_categories = [element.value for element in ElementType]

        if category not in allowed_categories:
            raise ValueError(f"Invalid {category=}, pick from {allowed_categories}")

        if "block" in category:
            return category[0] in [el.block for el in self.elements]

        return any(getattr(el, f"is_{category}") for el in self.elements)

    def _parse_formula(self, formula: str, strict: bool = True) -> dict[str, float]:
        """
        Args:
            formula (str): A string formula, e.g. Fe2O3, Li3Fe2(PO4)3.
            strict (bool): Whether to throw an error if formula string is invalid (e.g. empty).
                Defaults to True.

        Returns:
            Composition with that formula.

        Notes:
            In the case of Metallofullerene formula (e.g. Y3N@C80),
            the @ mark will be dropped and passed to parser.
        """
        # Raise error if formula contains special characters or only spaces and/or numbers
        if strict and re.match(r"[\s\d.*/]*$", formula):
            raise ValueError(f"Invalid {formula=}")

        # For Metallofullerene like "Y3N@C80"
        formula = formula.replace("@", "")
        # Square brackets are used in formulas to denote coordination complexes (gh-3583)
        formula = formula.replace("[", "(")
        formula = formula.replace("]", ")")

        def get_sym_dict(form: str, factor: float) -> dict[str, float]:
            sym_dict: dict[str, float] = defaultdict(float)
            for match in re.finditer(r"([A-Z][a-z]*)\s*([-*\.e\d]*)", form):
                el = match[1]
                amt = 1.0
                if match[2].strip() != "":
                    amt = float(match[2])
                sym_dict[el] += amt * factor
                form = form.replace(match.group(), "", 1)
            if form.strip():
                raise ValueError(f"{form} is an invalid formula!")
            return sym_dict

        match = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
        while match:
            factor = 1.0
            if match[2] != "":
                factor = float(match[2])
            unit_sym_dict = get_sym_dict(match[1], factor)
            expanded_sym = "".join(f"{el}{amt}" for el, amt in unit_sym_dict.items())
            expanded_formula = formula.replace(match.group(), expanded_sym, 1)
            formula = expanded_formula
            match = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
        return get_sym_dict(formula, 1)

    @property
    def anonymized_formula(self) -> str:
        """An anonymized formula. Unique species are arranged in ordering of
        increasing amounts and assigned ascending alphabets. Useful for
        prototyping formulas. For example, all stoichiometric perovskites have
        anonymized_formula ABC3.
        """
        reduced = self.element_composition
        if all(val == int(val) for val in self.values()):
            reduced /= gcd(*(int(i) for i in self.values()))

        anon = ""
        for elem, amt in zip(string.ascii_uppercase, sorted(reduced.values())):
            if amt == 1:
                amt_str = ""
            elif abs(amt % 1) < 1e-8:
                amt_str = str(int(amt))
            else:
                amt_str = str(amt)
            anon += f"{elem}{amt_str}"
        return anon

    @property
    def valid(self) -> bool:
        """True if Composition contains valid elements or species and
        False if the Composition contains any dummy species.
        """
        return not any(isinstance(el, DummySpecies) for el in self.elements)

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Create a composition from a dict generated by as_dict(). Strictly not
        necessary given that the standard constructor already takes in such an
        input, but this method preserves the standard pymatgen API of having
        from_dict methods to reconstitute objects generated by as_dict(). Allows
        for easier introspection.

        Args:
            dct (dict): {symbol: amount} dict.
        """
        return cls(dct)

    @classmethod
    def from_weight_dict(cls, weight_dict: dict[SpeciesLike, float]) -> Self:
        """Create a Composition based on a dict of atomic fractions calculated
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
                {"Fe": 4.0, "O": 6.0}.
        """
        dct: dict[str, float] = defaultdict(float)
        for el, amt in self.items():
            dct[el.symbol] += amt
        return dict(dct)

    def as_dict(self) -> dict[str, float]:
        """Subtly different from get_el_amt_dict in that they keys here are str(Element)
        instead of Element.symbol.

        Returns:
            dict[str, float]: element symbol and (unreduced) amount. E.g.
                {"Fe": 4.0, "O": 6.0} or {"Fe3+": 4.0, "O2-": 6.0}
        """
        dct: dict[str, float] = defaultdict(float)
        for el, amt in self.items():
            dct[str(el)] += amt
        return dict(dct)

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
    def to_data_dict(self) -> dict[str, Any]:
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

    @property
    def charge(self) -> float | None:
        """Total charge based on oxidation states. If any oxidation states
        are None or they're all 0, returns None. Use add_charges_from_oxi_state_guesses to
        assign oxidation states to elements based on charge balancing.
        """
        warnings.warn(
            "Composition.charge is experimental and may produce incorrect results. Use with "
            "caution and open a GitHub issue pinging @janosh to report bad behavior."
        )
        oxi_states = [getattr(specie, "oxi_state", None) for specie in self]
        if {*oxi_states} <= {0, None}:
            # all oxidation states are None or 0
            return None
        return sum(oxi * amt for oxi, amt in zip(oxi_states, self.values()))

    @property
    def charge_balanced(self) -> bool | None:
        """True if composition is charge balanced, False otherwise. If any oxidation states
        are None, returns None. Use add_charges_from_oxi_state_guesses to assign oxidation
        states to elements.
        """
        warnings.warn(
            "Composition.charge_balanced is experimental and may produce incorrect results. "
            "Use with caution and open a GitHub issue pinging @janosh to report bad behavior."
        )
        if self.charge is None:
            if {getattr(el, "oxi_state", None) for el in self} == {0}:
                # all oxidation states are 0. this usually means no way of combining oxidation states
                # to get a zero charge was found, so the composition is not charge balanced
                return False
            return None
        return abs(self.charge) < type(self).charge_balanced_tolerance

    def oxi_state_guesses(
        self,
        oxi_states_override: dict | None = None,
        target_charge: float = 0,
        all_oxi_states: bool = False,
        max_sites: int | None = None,
    ) -> tuple[dict[str, float]]:
        """Check if the composition is charge-balanced and returns back all
        charge-balanced oxidation state combinations. Composition must have
        integer values. Note that more num_atoms in the composition gives
        more degrees of freedom. e.g. if possible oxidation states of
        element X are [2,4] and Y are [-3], then XY is not charge balanced
        but X2Y2 is. Results are returned from most to least probable based
        on ICSD statistics. Use max_sites to improve performance if needed.

        Args:
            oxi_states_override (dict): dict of str->list to override an
                element's common oxidation states, e.g. {"V": [2,3,4,5]}
            target_charge (int): the desired total charge on the structure.
                Default is 0 signifying charge balance.
            all_oxi_states (bool): if True, all oxidation states of an element, even rare ones, are used in the search
                for guesses. However, the full oxidation state list is *very* inclusive and can produce nonsensical
                results. If False, the icsd_oxidation_states list is used when present, or the common_oxidation_states
                is used when icsd_oxidation_states is not present. These oxidation states lists comprise more
                commonly occurring oxidation states and results in more reliable guesses, albeit at the cost of
                missing some uncommon situations. The default is False.
            max_sites (int): if possible, will reduce Compositions to at most
                this many sites to speed up oxidation state guesses. If the
                composition cannot be reduced to this many sites a ValueError
                will be raised. Set to -1 to just reduce fully. If set to a
                number less than -1, the formula will be fully reduced but a
                ValueError will be thrown if the number of atoms in the reduced
                formula is greater than abs(max_sites).

        Returns:
            list[dict]: each dict reports an element symbol and average
                oxidation state across all sites in that composition. If the
                composition is not charge balanced, an empty list is returned.
        """
        if len(self.elements) == 1:
            return ({self.elements[0].symbol: 0.0},)
        return self._get_oxi_state_guesses(all_oxi_states, max_sites, oxi_states_override, target_charge)[0]

    def replace(self, elem_map: dict[str, str | dict[str, float]]) -> Self:
        """Replace elements in a composition. Returns a new Composition, leaving the old one unchanged.

        Args:
            elem_map (dict[str, str | dict[str, float]]): dict of elements or species to swap. E.g.
                {"Li": "Na"} performs a Li for Na substitution. The target can be a {species: factor} dict. For
                example, in Fe2O3 you could map {"Fe": {"Mg": 0.5, "Cu":0.5}} to obtain MgCuO3.

        Returns:
            Composition: New object with elements remapped according to elem_map.
        """
        # Drop inapplicable substitutions
        invalid_elems = [key for key in elem_map if key not in self]
        if invalid_elems:
            warnings.warn(
                "Some elements to be substituted are not present in composition. Please check your input. "
                f"Problematic element = {invalid_elems}; {self}"
            )
        for elem in invalid_elems:
            elem_map.pop(elem)

        # Start with elements that remain unchanged (not in elem_map)
        new_comp = {elem: amount for elem, amount in self.as_dict().items() if elem not in elem_map}

        for old_elem, new_elem in elem_map.items():
            amount = self[old_elem]

            # Build a dictionary of substitutions to be made
            subs = {}
            if isinstance(new_elem, dict):
                for el, factor in new_elem.items():
                    subs[el] = factor * amount
            else:
                subs = {new_elem: amount}

            # Apply the substitutions to the new composition
            for el, amt in subs.items():
                if el in new_comp:
                    new_comp[el] += amt
                else:
                    new_comp[el] = amt

                # Check for ambiguous input (see issue #2553)
                if el in self:
                    warnings.warn(
                        f"Same element ({el}) in both the keys and values of the substitution!"
                        "This can be ambiguous, so be sure to check your result."
                    )

        return type(self)(new_comp)

    def add_charges_from_oxi_state_guesses(
        self,
        oxi_states_override: dict | None = None,
        target_charge: float = 0,
        all_oxi_states: bool = False,
        max_sites: int | None = None,
    ) -> Self:
        """Assign oxidation states based on guessed oxidation states.

        See `oxi_state_guesses` for an explanation of how oxidation states are
        guessed. This operation uses the set of oxidation states for each site
        that were determined to be most likely from the oxidation state guessing
        routine.

        Args:
            oxi_states_override (dict[str, list[float]]): Override an
                element's common oxidation states, e.g. {"V": [2, 3, 4, 5]}
            target_charge (float): the desired total charge on the structure.
                Default is 0 signifying charge balance.
            all_oxi_states (bool): if True, all oxidation states of an element, even rare ones, are used in the search
                for guesses. However, the full oxidation state list is *very* inclusive and can produce nonsensical
                results. If False, the icsd_oxidation_states list is used when present, or the common_oxidation_states
                is used when icsd_oxidation_states is not present. These oxidation states lists comprise more
                commonly occurring oxidation states and results in more reliable guesses, albeit at the cost of
                missing some uncommon situations. The default is False.
            max_sites (int): If possible, will reduce Compositions to at most
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
            return type(self)({Species(e, 0): f for e, f in self.items()})

        # Generate the species
        species = []
        for el, charges in oxidation_states[0].items():
            species.extend([Species(el, c) for c in charges])

        # Return the new object
        return type(self)(collections.Counter(species))

    def remove_charges(self) -> Self:
        """Get a new Composition with charges from each Species removed.

        Returns:
            Composition object without charge decoration, for example
            {"Fe3+": 2.0, "O2-":3.0} becomes {"Fe": 2.0, "O":3.0}
        """
        dct: dict[Element, float] = defaultdict(float)
        for specie, amt in self.items():
            dct[Element(specie.symbol)] += amt
        return type(self)(dct)

    def _get_oxi_state_guesses(
        self,
        all_oxi_states: bool,
        max_sites: int | None,
        oxi_states_override: dict[str, list] | None,
        target_charge: float,
    ) -> tuple[tuple, tuple]:
        """Utility operation for guessing oxidation states.

        See `oxi_state_guesses` for full details. This operation does the
        calculation of the most likely oxidation states

        Args:
            oxi_states_override (dict): dict of str->list to override an element's common oxidation states, e.g.
                {"V": [2,3,4,5]}.
            target_charge (float): the desired total charge on the structure. Default is 0 signifying charge balance.
            all_oxi_states (bool): if True, all oxidation states of an element, even rare ones, are used in the search
                for guesses. However, the full oxidation state list is *very* inclusive and can produce nonsensical
                results. If False, the icsd_oxidation_states list is used when present, or the common_oxidation_states
                is used when icsd_oxidation_states is not present. These oxidation states lists comprise more
                commonly occurring oxidation states and results in more reliable guesses, albeit at the cost of
                missing some uncommon situations. The default is False.
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
        # Reduce Composition if necessary
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
        if not type(self).oxi_prob:
            all_data = loadfn(f"{module_dir}/../analysis/icsd_bv.yaml")
            type(self).oxi_prob = {Species.from_str(sp): data for sp, data in all_data["occurrence"].items()}
        oxi_states_override = oxi_states_override or {}
        # Assert Composition only has integer amounts
        if not all(amt == int(amt) for amt in comp.values()):
            raise ValueError("Charge balance analysis requires integer values in Composition!")

        # For each element, determine all possible sum of oxidations
        # (taking into account nsites for that particular element)
        el_amt = comp.get_el_amt_dict()
        elements = list(el_amt)
        el_sums: list = []  # matrix: dim1= el_idx, dim2=possible sums
        el_sum_scores: defaultdict = defaultdict(set)  # dict of el_idx, sum -> score
        el_best_oxid_combo: dict = {}  # dict of el_idx, sum -> oxid combo with best score
        for idx, el in enumerate(elements):
            el_sum_scores[idx] = {}
            el_best_oxid_combo[idx] = {}
            el_sums.append([])
            if oxi_states_override.get(el):
                oxids: list | tuple = oxi_states_override[el]
            elif all_oxi_states:
                oxids = Element(el).oxidation_states
            else:
                oxids = Element(el).icsd_oxidation_states or Element(el).common_oxidation_states

            # Get all possible combinations of oxidation states
            # and sum each combination
            for oxid_combo in combinations_with_replacement(oxids, int(el_amt[el])):
                # List this sum as a possible option
                oxid_sum = sum(oxid_combo)
                if oxid_sum not in el_sums[idx]:
                    el_sums[idx].append(oxid_sum)

                # Determine how probable is this combo?
                score = sum(type(self).oxi_prob.get(Species(el, o), 0) for o in oxid_combo)  # type: ignore[union-attr]

                # If it is the most probable combo for a certain sum,
                # store the combination
                if oxid_sum not in el_sum_scores[idx] or score > el_sum_scores[idx].get(oxid_sum, 0):
                    el_sum_scores[idx][oxid_sum] = score
                    el_best_oxid_combo[idx][oxid_sum] = oxid_combo

        # Determine which combination of oxidation states for each element
        # is the most probable
        all_sols = []  # will contain all solutions
        all_oxid_combo = []  # will contain the best combination of oxidation states for each site
        all_scores = []  # will contain a score for each solution
        for x in product(*el_sums):
            # Each x is a trial of one possible oxidation sum for each element
            if sum(x) == target_charge:  # charge balance condition
                el_sum_sol = dict(zip(elements, x))  # element->oxid_sum
                # Normalize oxid_sum by amount to get avg oxid state
                sol = {el: v / el_amt[el] for el, v in el_sum_sol.items()}
                # Add the solution to the list of solutions
                all_sols.append(sol)

                # Determine the score for this solution
                score = 0
                for idx, v in enumerate(x):
                    score += el_sum_scores[idx][v]
                all_scores.append(score)

                # Collect the combination of oxidation states for each site
                all_oxid_combo.append({e: el_best_oxid_combo[idx][v] for idx, (e, v) in enumerate(zip(elements, x))})

        # Sort the solutions from highest to lowest score
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
        return tuple(all_sols), tuple(all_oxid_combo)

    @staticmethod
    def ranked_compositions_from_indeterminate_formula(
        fuzzy_formula: str, lock_if_strict: bool = True
    ) -> list[Composition]:
        """Takes in a formula where capitalization might not be correctly entered,
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
        factor: float = 1,
    ) -> Generator[tuple[Composition, int], None, None]:
        """A recursive helper method for formula parsing that helps in
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

        def _parse_chomp_and_rank(match, formula: str, m_dict: dict[str, float], m_points: int) -> tuple:
            """A helper method for formula parsing that helps in interpreting and
            ranking indeterminate formulas.

            Author: Anubhav Jain.

            Args:
                match: A regex match, with the first group being the element and
                    the second group being the amount
                formula: The formula part containing the match
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
            el = match[1]
            if len(el) > 2 or len(el) < 1:
                raise ValueError("Invalid element symbol entered!")
            amt = float(match[2]) if match[2].strip() != "" else 1

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
                return formula.replace(match.group(), "", 1), m_dict, m_points + points

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
                m_form1, m_dict1, m_points1 = _parse_chomp_and_rank(m1, m_form1, m_dict1, m_points1)
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
                m_form2, m_dict2, m_points2 = _parse_chomp_and_rank(m2, m_form2, m_dict2, m_points2)
                if m_dict2:
                    # there was a real match
                    for match in Composition._comps_from_fuzzy_formula(m_form2, m_dict2, m_points2, factor):
                        yield match


def reduce_formula(
    sym_amt: dict[str, float] | dict[str, int],
    iupac_ordering: bool = False,
) -> tuple[str, float]:
    """Helper function to reduce a sym_amt dict to a reduced formula and factor.

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
        tuple[str, float]: reduced formula and factor.
    """
    syms = sorted(sym_amt, key=lambda x: [get_el_sp(x).X, x])

    syms = list(filter(lambda x: abs(sym_amt[x]) > Composition.amount_tolerance, syms))

    factor = 1
    # Enforce integers for doing gcd.
    if all(int(i) == i for i in sym_amt.values()):
        factor = abs(gcd(*(int(i) for i in sym_amt.values())))

    poly_anions = []
    # if the composition contains a poly anion
    if len(syms) >= 3 and get_el_sp(syms[-1]).X - get_el_sp(syms[-2]).X < 1.65:
        poly_sym_amt = {syms[i]: sym_amt[syms[i]] / factor for i in [-2, -1]}
        poly_form, poly_factor = reduce_formula(poly_sym_amt, iupac_ordering=iupac_ordering)

        if poly_factor != 1:
            poly_anions.append(f"({poly_form}){poly_factor}")

    syms = syms[: len(syms) - 2 if poly_anions else len(syms)]

    if iupac_ordering:
        syms = sorted(syms, key=lambda x: [get_el_sp(x).iupac_ordering, x])

    reduced_form: list[str] = []
    for sym in syms:
        norm_amt = sym_amt[sym] * 1.0 / factor
        reduced_form.extend((sym, str(formula_double_format(norm_amt))))

    return "".join([*reduced_form, *poly_anions]), factor


class ChemicalPotential(dict, MSONable):
    """Represent set of chemical potentials. Can be: multiplied/divided by a Number
    multiplied by a Composition (returns an energy) added/subtracted with other ChemicalPotentials.
    """

    def __init__(self, *args, **kwargs) -> None:
        """
        Args:
            *args: any valid dict init arguments
            **kwargs: any valid dict init arguments.
        """
        dct = dict(*args, **kwargs)
        super().__init__((get_el_sp(key), val) for key, val in dct.items())
        if len(dct) != len(self):
            raise ValueError("Duplicate potential specified")

    def __mul__(self, other: object) -> Self:
        if isinstance(other, (int, float)):
            return type(self)({key: val * other for key, val in self.items()})
        return NotImplemented

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> Self:
        if isinstance(other, (int, float)):
            return type(self)({key: val / other for key, val in self.items()})
        return NotImplemented

    __div__ = __truediv__

    def __sub__(self, other: object) -> Self:
        if isinstance(other, type(self)):
            els = {*self} | {*other}
            return type(self)({e: self.get(e, 0) - other.get(e, 0) for e in els})
        return NotImplemented

    def __add__(self, other: object) -> Self:
        if isinstance(other, type(self)):
            els = {*self} | {*other}
            return type(self)({e: self.get(e, 0) + other.get(e, 0) for e in els})
        return NotImplemented

    def __repr__(self) -> str:
        return f"ChemPots: {super()!r}"

    def get_energy(self, composition: Composition, strict: bool = True) -> float:
        """Calculate the energy of a composition.

        Args:
            composition (Composition): input composition
            strict (bool): Whether all potentials must be specified
        """
        if strict and (missing := set(composition) - set(self)):
            raise ValueError(f"Potentials not specified for {missing}")
        return sum(self.get(key, 0) * val for key, val in composition.items())


class CompositionError(Exception):
    """Exception class for composition errors."""
