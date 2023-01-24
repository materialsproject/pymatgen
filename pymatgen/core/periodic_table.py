# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""Module contains classes presenting Element and Species (Element + oxidation state) and PeriodicTable."""

from __future__ import annotations

import ast
import functools
import json
import re
import warnings
from collections import Counter
from enum import Enum
from itertools import combinations, product
from pathlib import Path
from typing import Any, Callable, Literal

import numpy as np
from monty.json import MSONable

from pymatgen.core.units import SUPPORTED_UNIT_NAMES, FloatWithUnit, Length, Mass, Unit
from pymatgen.util.string import Stringify, formula_double_format

# Loads element data from json file
with open(str(Path(__file__).absolute().parent / "periodic_table.json")) as f:
    _pt_data = json.load(f)

_pt_row_sizes = (2, 8, 8, 18, 18, 32, 32)


@functools.total_ordering
class ElementBase(Enum):
    """Element class defined without any enum values so it can be subclassed."""

    def __init__(self, symbol: str):
        """
        Basic immutable element object with all relevant properties.

        Only one instance of Element for each symbol is stored after creation,
        ensuring that a particular element behaves like a singleton. For all
        attributes, missing data (i.e., data for which is not available) is
        represented by a None unless otherwise stated.

        Args:
            symbol (str): Element symbol, e.g., "H", "Fe"

        .. attribute:: Z

            Atomic number

        .. attribute:: symbol

            Element symbol

        .. attribute:: long_name

           Long name for element. E.g., "Hydrogen".

        .. attribute:: atomic_radius_calculated

            Calculated atomic radius for the element. This is the empirical value.
            Data is obtained from
            http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page).

        .. attribute:: van_der_waals_radius

            Van der Waals radius for the element. This is the empirical
            value determined from critical reviews of X-ray diffraction, gas kinetic
            collision cross-section, and other experimental data by Bondi and later
            workers. The uncertainty in these values is on the order of 0.1 Å.

            Data are obtained from

            "Atomic Radii of the Elements" in CRC Handbook of Chemistry and Physics,
                91st Ed.; Haynes, W.M., Ed.; CRC Press: Boca Raton, FL, 2010.

        .. attribute:: mendeleev_no

            Mendeleev number from definition given by Pettifor, D. G. (1984).
            A chemical scale for crystal-structure maps. Solid State Communications,
            51 (1), 31-34

        .. attribute:: electrical_resistivity

            Electrical resistivity

        .. attribute:: velocity_of_sound

            Velocity of sound

        .. attribute:: reflectivity

            Reflectivity

        .. attribute:: refractive_index

            Refractice index

        .. attribute:: poissons_ratio

            Poisson's ratio

        .. attribute:: molar_volume

            Molar volume

        .. attribute:: electronic_structure

            Electronic structure.
            E.g., The electronic structure for Fe is represented as
            [Ar].3d6.4s2

        .. attribute:: atomic_orbitals

            Atomic Orbitals. Energy of the atomic orbitals as a dict.
            E.g., The orbitals energies in eV are represented as
            {'1s': -1.0, '2s': -0.1}
            Data is obtained from
            https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations
            The LDA values for neutral atoms are used

        .. attribute:: thermal_conductivity

            Thermal conductivity

        .. attribute:: boiling_point

            Boiling point

        .. attribute:: melting_point

            Melting point

        .. attribute:: critical_temperature

            Critical temperature

        .. attribute:: superconduction_temperature

            Superconduction temperature

        .. attribute:: liquid_range

            Liquid range

        .. attribute:: bulk_modulus

            Bulk modulus

        .. attribute:: youngs_modulus

            Young's modulus

        .. attribute:: brinell_hardness

            Brinell hardness

        .. attribute:: rigidity_modulus

            Rigidity modulus

        .. attribute:: mineral_hardness

            Mineral hardness

        .. attribute:: vickers_hardness

            Vicker's hardness

        .. attribute:: density_of_solid

            Density of solid phase

        .. attribute:: coefficient_of_linear_thermal_expansion

            Coefficient of linear thermal expansion

        .. attribute:: ground_level

            Ground level for element

        .. attribute:: ionization_energies

            List of ionization energies. First value is the first ionization energy, second is the second ionization
            energy, etc. Note that this is zero-based indexing! So Element.ionization_energies[0] refer to the 1st
            ionization energy. Values are from the NIST Atomic Spectra Database. Missing values are None.
        """
        self.symbol = symbol
        d = _pt_data[symbol]

        # Store key variables for quick access
        self.Z = d["Atomic no"]

        at_r = d.get("Atomic radius", "no data")
        if str(at_r).startswith("no data"):
            self._atomic_radius = None
        else:
            self._atomic_radius = Length(at_r, "ang")
        self._atomic_mass = Mass(d["Atomic mass"], "amu")
        self.long_name = d["Name"]
        self._data = d

    @property
    def X(self) -> float:
        """
        :return: Electronegativity of element. Note that if an element does not
            have an electronegativity, a NaN float is returned.
        """
        if "X" in self._data:
            return self._data["X"]
        warnings.warn(
            f"No electronegativity for {self.symbol}. Setting to NaN. This has no physical meaning, "
            "and is mainly done to avoid errors caused by the code expecting a float."
        )
        return float("NaN")

    @property
    def atomic_radius(self) -> FloatWithUnit | None:
        """
        Returns:
            float | None: The atomic radius of the element in Ångstroms. Can be None for
                some elements like noble gases.
        """
        return self._atomic_radius

    @property
    def atomic_mass(self) -> FloatWithUnit:
        """
        Returns:
            float: The atomic mass of the element in amu.
        """
        return self._atomic_mass

    def __getattr__(self, item: str) -> Any:
        """Key access to available element data.

        Args:
            item (str): Attribute name.

        Raises:
            AttributeError: If item not in _pt_data.
        """
        if item in [
            "mendeleev_no",
            "electrical_resistivity",
            "velocity_of_sound",
            "reflectivity",
            "refractive_index",
            "poissons_ratio",
            "molar_volume",
            "thermal_conductivity",
            "boiling_point",
            "melting_point",
            "critical_temperature",
            "superconduction_temperature",
            "liquid_range",
            "bulk_modulus",
            "youngs_modulus",
            "brinell_hardness",
            "rigidity_modulus",
            "mineral_hardness",
            "vickers_hardness",
            "density_of_solid",
            "atomic_radius_calculated",
            "van_der_waals_radius",
            "atomic_orbitals",
            "coefficient_of_linear_thermal_expansion",
            "ground_state_term_symbol",
            "valence",
            "ground_level",
            "ionization_energies",
        ]:
            kstr = item.capitalize().replace("_", " ")
            val = self._data.get(kstr, None)
            if str(val).startswith("no data"):
                val = None
            elif isinstance(val, (list, dict)):
                pass
            else:
                try:
                    val = float(val)
                except ValueError:
                    nobracket = re.sub(r"\(.*\)", "", val)
                    toks = nobracket.replace("about", "").strip().split(" ", 1)
                    if len(toks) == 2:
                        try:
                            if "10<sup>" in toks[1]:
                                base_power = re.findall(r"([+-]?\d+)", toks[1])
                                factor = "e" + base_power[1]
                                if toks[0] in ["&gt;", "high"]:
                                    toks[0] = "1"  # return the border value
                                toks[0] += factor
                                if item == "electrical_resistivity":
                                    unit = "ohm m"
                                elif item == "coefficient_of_linear_thermal_expansion":
                                    unit = "K^-1"
                                else:
                                    unit = toks[1]
                                val = FloatWithUnit(toks[0], unit)
                            else:
                                unit = toks[1].replace("<sup>", "^").replace("</sup>", "").replace("&Omega;", "ohm")
                                units = Unit(unit)
                                if set(units).issubset(SUPPORTED_UNIT_NAMES):
                                    val = FloatWithUnit(toks[0], unit)
                        except ValueError:
                            # Ignore error. val will just remain a string.
                            pass
            return val
        raise AttributeError(f"Element has no attribute {item}!")

    @property
    def data(self) -> dict[str, Any]:
        """
        Returns dict of data for element.
        """
        return self._data.copy()

    @property
    def ionization_energy(self) -> float:
        """
        First ionization energy of element.
        """
        return self._data["Ionization energies"][0]

    @property
    def electron_affinity(self) -> float:
        """
        The amount of energy released when an electron is attached to a neutral atom.
        """
        return self._data["Electron affinity"]

    @property
    def electronic_structure(self) -> str:
        """
        Electronic structure as string, with only valence electrons.
        E.g., The electronic structure for Fe is represented as '[Ar].3d6.4s2'
        """
        return re.sub("</*sup>", "", self._data["Electronic structure"])

    @property
    def average_ionic_radius(self) -> FloatWithUnit:
        """
        Average ionic radius for element (with units). The average is taken
        over all oxidation states of the element for which data is present.
        """
        if "Ionic radii" in self._data:
            radii = self._data["Ionic radii"]
            radius = sum(radii.values()) / len(radii)
        else:
            radius = 0.0
        return FloatWithUnit(radius, "ang")

    @property
    def average_cationic_radius(self) -> FloatWithUnit:
        """
        Average cationic radius for element (with units). The average is
        taken over all positive oxidation states of the element for which
        data is present.
        """
        if "Ionic radii" in self._data:
            radii = [v for k, v in self._data["Ionic radii"].items() if int(k) > 0]
            if radii:
                return FloatWithUnit(sum(radii) / len(radii), "ang")
        return FloatWithUnit(0.0, "ang")

    @property
    def average_anionic_radius(self) -> float:
        """
        Average anionic radius for element (with units). The average is
        taken over all negative oxidation states of the element for which
        data is present.
        """
        if "Ionic radii" in self._data:
            radii = [v for k, v in self._data["Ionic radii"].items() if int(k) < 0]
            if radii:
                return FloatWithUnit(sum(radii) / len(radii), "ang")
        return FloatWithUnit(0.0, "ang")

    @property
    def ionic_radii(self) -> dict[int, float]:
        """
        All ionic radii of the element as a dict of
        {oxidation state: ionic radii}. Radii are given in angstrom.
        """
        if "Ionic radii" in self._data:
            return {int(k): FloatWithUnit(v, "ang") for k, v in self._data["Ionic radii"].items()}
        return {}

    @property
    def number(self) -> int:
        """Alternative attribute for atomic number Z"""
        return self.Z

    @property
    def max_oxidation_state(self) -> float:
        """Maximum oxidation state for element"""
        if "Oxidation states" in self._data:
            return max(self._data["Oxidation states"])
        return 0

    @property
    def min_oxidation_state(self) -> float:
        """Minimum oxidation state for element"""
        if "Oxidation states" in self._data:
            return min(self._data["Oxidation states"])
        return 0

    @property
    def oxidation_states(self) -> tuple[int, ...]:
        """Tuple of all known oxidation states"""
        return tuple(int(x) for x in self._data.get("Oxidation states", []))

    @property
    def common_oxidation_states(self) -> tuple[int, ...]:
        """Tuple of common oxidation states"""
        return tuple(self._data.get("Common oxidation states", []))

    @property
    def icsd_oxidation_states(self) -> tuple[int, ...]:
        """Tuple of all oxidation states with at least 10 instances in
        ICSD database AND at least 1% of entries for that element"""
        return tuple(self._data.get("ICSD oxidation states", []))

    @property
    def metallic_radius(self) -> float:
        """
        Metallic radius of the element. Radius is given in ang.
        """
        return FloatWithUnit(self._data["Metallic radius"], "ang")

    @property
    def full_electronic_structure(self) -> list[tuple[int, str, int]]:
        """
        Full electronic structure as tuple.
        E.g., The electronic structure for Fe is represented as:
        [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6),
        (3, "d", 6), (4, "s", 2)]
        """
        estr = self.electronic_structure

        def parse_orbital(orbstr):
            m = re.match(r"(\d+)([spdfg]+)(\d+)", orbstr)
            if m:
                return int(m.group(1)), m.group(2), int(m.group(3))
            return orbstr

        data = [parse_orbital(s) for s in estr.split(".")]
        if data[0][0] == "[":
            sym = data[0].replace("[", "").replace("]", "")
            data = list(Element(sym).full_electronic_structure) + data[1:]
        return data

    @property
    def valence(self):
        """
        From full electron config obtain valence subshell angular moment (L) and number of valence e- (v_e)
        """
        # The number of valence of noble gas is 0
        if self.group == 18:
            return np.nan, 0

        L_symbols = "SPDFGHIKLMNOQRTUVWXYZ"
        valence = []
        full_electron_config = self.full_electronic_structure
        last_orbital = full_electron_config[-1]
        for n, l_symbol, ne in full_electron_config:
            l = L_symbols.lower().index(l_symbol)
            if ne < (2 * l + 1) * 2:
                valence.append((l, ne))
            # check for full last shell (e.g. column 2)
            elif (n, l_symbol, ne) == last_orbital and ne == (2 * l + 1) * 2 and len(valence) == 0:
                valence.append((l, ne))
        if len(valence) > 1:
            raise ValueError("Ambiguous valence")

        return valence[0]

    @property
    def term_symbols(self) -> list[list[str]]:
        """
        All possible  Russell-Saunders term symbol of the Element.
        eg. L = 1, n_e = 2 (s2) returns [['1D2'], ['3P0', '3P1', '3P2'], ['1S0']]
        """
        L_symbols = "SPDFGHIKLMNOQRTUVWXYZ"

        L, v_e = self.valence

        # for one electron in subshell L
        ml = list(range(-L, L + 1))
        ms = [1 / 2, -1 / 2]
        # all possible configurations of ml,ms for one e in subshell L
        ml_ms = list(product(ml, ms))

        # Number of possible configurations for r electrons in subshell L.
        n = (2 * L + 1) * 2
        # the combination of n_e electrons configurations
        # C^{n}_{n_e}
        e_config_combs = list(combinations(range(n), v_e))

        # Total ML = sum(ml1, ml2), Total MS = sum(ms1, ms2)
        TL = [sum(ml_ms[comb[e]][0] for e in range(v_e)) for comb in e_config_combs]
        TS = [sum(ml_ms[comb[e]][1] for e in range(v_e)) for comb in e_config_combs]
        comb_counter = Counter(zip(TL, TS))

        term_symbols = []
        while sum(comb_counter.values()) > 0:
            # Start from the lowest freq combination,
            # which corresponds to largest abs(L) and smallest abs(S)
            L, S = min(comb_counter)

            J = list(np.arange(abs(L - S), abs(L) + abs(S) + 1))
            term_symbols.append([str(int(2 * (abs(S)) + 1)) + L_symbols[abs(L)] + str(j) for j in J])

            # Delete all configurations included in this term
            for ML in range(-L, L - 1, -1):
                for MS in np.arange(S, -S + 1, 1):
                    if (ML, MS) in comb_counter:

                        comb_counter[(ML, MS)] -= 1
                        if comb_counter[(ML, MS)] == 0:
                            del comb_counter[(ML, MS)]
        return term_symbols

    @property
    def ground_state_term_symbol(self):
        """
        Ground state term symbol
        Selected based on Hund's Rule
        """
        L_symbols = "SPDFGHIKLMNOQRTUVWXYZ"

        term_symbols = self.term_symbols
        term_symbol_flat = {
            term: {
                "multiplicity": int(term[0]),
                "L": L_symbols.index(term[1]),
                "J": float(term[2:]),
            }
            for term in sum(term_symbols, [])
        }

        multi = [int(item["multiplicity"]) for terms, item in term_symbol_flat.items()]
        max_multi_terms = {
            symbol: item for symbol, item in term_symbol_flat.items() if item["multiplicity"] == max(multi)
        }

        Ls = [item["L"] for terms, item in max_multi_terms.items()]
        max_L_terms = {symbol: item for symbol, item in term_symbol_flat.items() if item["L"] == max(Ls)}

        J_sorted_terms = sorted(max_L_terms.items(), key=lambda k: k[1]["J"])
        L, v_e = self.valence
        if v_e <= (2 * L + 1):
            return J_sorted_terms[0][0]
        return J_sorted_terms[-1][0]

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Element) and self.Z == other.Z

    def __hash__(self):
        return self.Z

    def __repr__(self):
        return "Element " + self.symbol

    def __str__(self):
        return self.symbol

    def __lt__(self, other):
        """
        Sets a default sort order for atomic species by electronegativity. Very
        useful for getting correct formulas. For example, FeO4PLi is
        automatically sorted into LiFePO4.
        """
        if not hasattr(other, "X") or not hasattr(other, "symbol"):
            return NotImplemented
        x1 = float("inf") if self.X != self.X else self.X
        x2 = float("inf") if other.X != other.X else other.X
        if x1 != x2:
            return x1 < x2

        # There are cases where the electronegativity are exactly equal.
        # We then sort by symbol.
        return self.symbol < other.symbol

    @staticmethod
    def from_Z(Z: int) -> Element:
        """
        Get an element from an atomic number.

        Args:
            Z (int): Atomic number

        Returns:
            Element with atomic number Z.
        """
        for sym, data in _pt_data.items():
            if data["Atomic no"] == Z:
                return Element(sym)
        raise ValueError(f"No element with this atomic number {Z}")

    @staticmethod
    def from_name(name: str) -> Element:
        """
        Get an element from its long name.

        Args:
            name: Long name of the element, e.g. 'Hydrogen' or
                  'Iron'. Not case-sensitive.
        Returns:
            Element with the name 'name'
        """
        for sym, data in _pt_data.items():
            if data["Name"] == name.capitalize():
                return Element(sym)
        raise ValueError(f"No element with the name {name}")

    @staticmethod
    def from_row_and_group(row: int, group: int) -> Element:
        """
        Returns an element from a row and group number.
        Important Note: For lanthanoids and actinoids, the row number must
        be 8 and 9, respectively, and the group number must be
        between 3 (La, Ac) and 17 (Lu, Lr). This is different than the
        value for Element(symbol).row and Element(symbol).group for these
        elements.

        Args:
            row (int): (pseudo) row number. This is the
                standard row number except for the lanthanoids
                and actinoids for which it is 8 or 9, respectively.
            group (int): (pseudo) group number. This is the
                standard group number except for the lanthanoids
                and actinoids for which it is 3 (La, Ac) to 17 (Lu, Lr).

        .. note::
            The 18 group number system is used, i.e., Noble gases are group 18.
        """
        for sym in _pt_data:
            el = Element(sym)
            if 57 <= el.Z <= 71:
                el_pseudorow = 8
                el_pseudogroup = (el.Z - 54) % 32
            elif 89 <= el.Z <= 103:
                el_pseudorow = 9
                el_pseudogroup = (el.Z - 54) % 32
            else:
                el_pseudorow = el.row
                el_pseudogroup = el.group
            if el_pseudorow == row and el_pseudogroup == group:
                return el
        raise ValueError("No element with this row and group!")

    @staticmethod
    def is_valid_symbol(symbol: str) -> bool:
        """
        Returns true if symbol is a valid element symbol.

        Args:
            symbol (str): Element symbol

        Returns:
            True if symbol is a valid element (e.g., "H"). False otherwise
            (e.g., "Zebra").
        """
        return symbol in Element.__members__

    @property
    def row(self) -> int:
        """
        Returns the periodic table row of the element.
        Note: For lanthanoids and actinoids, the row is always 6 or 7,
        respectively.
        """
        z = self.Z
        total = 0
        if 57 <= z <= 71:
            return 6
        if 89 <= z <= 103:
            return 7
        for i, size in enumerate(_pt_row_sizes):
            total += size
            if total >= z:
                return i + 1
        return 8

    @property
    def group(self) -> int:
        """
        Returns the periodic table group of the element.
        Note: For lanthanoids and actinoids, the group is always 3.
        """
        z = self.Z
        if z == 1:
            return 1
        if z == 2:
            return 18

        if 3 <= z <= 18:
            if (z - 2) % 8 == 0:
                return 18
            if (z - 2) % 8 <= 2:
                return (z - 2) % 8
            return 10 + (z - 2) % 8

        if 19 <= z <= 54:
            if (z - 18) % 18 == 0:
                return 18
            return (z - 18) % 18

        if (57 <= z <= 71) or (89 <= z <= 103):
            return 3

        if (z - 54) % 32 == 0:
            return 18
        if (z - 54) % 32 >= 18:
            return (z - 54) % 32 - 14
        return (z - 54) % 32

    @property
    def block(self) -> str:
        """
        Return the block character "s,p,d,f"
        """
        if (self.is_actinoid or self.is_lanthanoid) and self.Z not in [71, 103]:
            return "f"
        if self.is_actinoid or self.is_lanthanoid:
            return "d"
        if self.group in [1, 2]:
            return "s"
        if self.group in range(13, 19):
            return "p"
        if self.group in range(3, 13):
            return "d"
        raise ValueError("unable to determine block")

    @property
    def is_noble_gas(self) -> bool:
        """
        True if element is noble gas.
        """
        return self.Z in (2, 10, 18, 36, 54, 86, 118)

    @property
    def is_transition_metal(self) -> bool:
        """
        True if element is a transition metal.
        """
        ns = list(range(21, 31))
        ns.extend(list(range(39, 49)))
        ns.append(57)
        ns.extend(list(range(72, 81)))
        ns.append(89)
        ns.extend(list(range(104, 113)))
        return self.Z in ns

    @property
    def is_post_transition_metal(self) -> bool:
        """
        True if element is a post-transition or poor metal.
        """
        return self.symbol in ("Al", "Ga", "In", "Tl", "Sn", "Pb", "Bi")

    @property
    def is_rare_earth_metal(self) -> bool:
        """
        True if element is a rare earth metal.
        """
        return self.is_lanthanoid or self.is_actinoid

    @property
    def is_metal(self) -> bool:
        """
        True if is a metal.
        """
        return (
            self.is_alkali
            or self.is_alkaline
            or self.is_post_transition_metal
            or self.is_transition_metal
            or self.is_lanthanoid
            or self.is_actinoid
        )

    @property
    def is_metalloid(self) -> bool:
        """
        True if element is a metalloid.
        """
        return self.symbol in ("B", "Si", "Ge", "As", "Sb", "Te", "Po")

    @property
    def is_alkali(self) -> bool:
        """
        True if element is an alkali metal.
        """
        return self.Z in (3, 11, 19, 37, 55, 87)

    @property
    def is_alkaline(self) -> bool:
        """
        True if element is an alkaline earth metal (group II).
        """
        return self.Z in (4, 12, 20, 38, 56, 88)

    @property
    def is_halogen(self) -> bool:
        """
        True if element is a halogen.
        """
        return self.Z in (9, 17, 35, 53, 85)

    @property
    def is_chalcogen(self) -> bool:
        """
        True if element is a chalcogen.
        """
        return self.Z in (8, 16, 34, 52, 84)

    @property
    def is_lanthanoid(self) -> bool:
        """
        True if element is a lanthanoid.
        """
        return 56 < self.Z < 72

    @property
    def is_actinoid(self) -> bool:
        """
        True if element is a actinoid.
        """
        return 88 < self.Z < 104

    @property
    def is_quadrupolar(self) -> bool:
        """
        Checks if this element can be quadrupolar.
        """
        return len(self.data.get("NMR Quadrupole Moment", {})) > 0

    @property
    def nmr_quadrupole_moment(self) -> dict[str, FloatWithUnit]:
        """
        Get a dictionary the nuclear electric quadrupole moment in units of
        e*millibarns for various isotopes
        """
        return {k: FloatWithUnit(v, "mbarn") for k, v in self.data.get("NMR Quadrupole Moment", {}).items()}

    @property
    def iupac_ordering(self):
        """
        Ordering according to Table VI of "Nomenclature of Inorganic Chemistry
        (IUPAC Recommendations 2005)". This ordering effectively follows the
        groups and rows of the periodic table, except the Lanthanides, Actinides
        and hydrogen.
        """
        return self._data["IUPAC ordering"]

    def __deepcopy__(self, memo):
        return Element(self.symbol)

    @staticmethod
    def from_dict(d) -> Element:
        """
        Makes Element obey the general json interface used in pymatgen for
        easier serialization.
        """
        return Element(d["element"])

    def as_dict(self) -> dict[Literal["element", "@module", "@class"], str]:
        """
        Makes Element obey the general json interface used in pymatgen for
        easier serialization.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "element": self.symbol,
        }

    @staticmethod
    def print_periodic_table(filter_function: Callable | None = None):
        """
        A pretty ASCII printer for the periodic table, based on some
        filter_function.

        Args:
            filter_function: A filtering function taking an Element as input
                and returning a boolean. For example, setting
                filter_function = lambda el: el.X > 2 will print a periodic
                table containing only elements with electronegativity > 2.
        """
        for row in range(1, 10):
            rowstr = []
            for group in range(1, 19):
                try:
                    el = Element.from_row_and_group(row, group)
                except ValueError:
                    el = None
                if el and ((not filter_function) or filter_function(el)):
                    rowstr.append(f"{el.symbol:3s}")
                else:
                    rowstr.append("   ")
            print(" ".join(rowstr))


@functools.total_ordering
class Element(ElementBase):
    """Enum representing an element in the periodic table."""

    # This name = value convention is redundant and dumb, but unfortunately is
    # necessary to preserve backwards compatibility with a time when Element is
    # a regular object that is constructed with Element(symbol).
    H = "H"
    He = "He"
    Li = "Li"
    Be = "Be"
    B = "B"
    C = "C"
    N = "N"
    O = "O"
    F = "F"
    Ne = "Ne"
    Na = "Na"
    Mg = "Mg"
    Al = "Al"
    Si = "Si"
    P = "P"
    S = "S"
    Cl = "Cl"
    Ar = "Ar"
    K = "K"
    Ca = "Ca"
    Sc = "Sc"
    Ti = "Ti"
    V = "V"
    Cr = "Cr"
    Mn = "Mn"
    Fe = "Fe"
    Co = "Co"
    Ni = "Ni"
    Cu = "Cu"
    Zn = "Zn"
    Ga = "Ga"
    Ge = "Ge"
    As = "As"
    Se = "Se"
    Br = "Br"
    Kr = "Kr"
    Rb = "Rb"
    Sr = "Sr"
    Y = "Y"
    Zr = "Zr"
    Nb = "Nb"
    Mo = "Mo"
    Tc = "Tc"
    Ru = "Ru"
    Rh = "Rh"
    Pd = "Pd"
    Ag = "Ag"
    Cd = "Cd"
    In = "In"
    Sn = "Sn"
    Sb = "Sb"
    Te = "Te"
    I = "I"
    Xe = "Xe"
    Cs = "Cs"
    Ba = "Ba"
    La = "La"
    Ce = "Ce"
    Pr = "Pr"
    Nd = "Nd"
    Pm = "Pm"
    Sm = "Sm"
    Eu = "Eu"
    Gd = "Gd"
    Tb = "Tb"
    Dy = "Dy"
    Ho = "Ho"
    Er = "Er"
    Tm = "Tm"
    Yb = "Yb"
    Lu = "Lu"
    Hf = "Hf"
    Ta = "Ta"
    W = "W"
    Re = "Re"
    Os = "Os"
    Ir = "Ir"
    Pt = "Pt"
    Au = "Au"
    Hg = "Hg"
    Tl = "Tl"
    Pb = "Pb"
    Bi = "Bi"
    Po = "Po"
    At = "At"
    Rn = "Rn"
    Fr = "Fr"
    Ra = "Ra"
    Ac = "Ac"
    Th = "Th"
    Pa = "Pa"
    U = "U"
    Np = "Np"
    Pu = "Pu"
    Am = "Am"
    Cm = "Cm"
    Bk = "Bk"
    Cf = "Cf"
    Es = "Es"
    Fm = "Fm"
    Md = "Md"
    No = "No"
    Lr = "Lr"
    Rf = "Rf"
    Db = "Db"
    Sg = "Sg"
    Bh = "Bh"
    Hs = "Hs"
    Mt = "Mt"
    Ds = "Ds"
    Rg = "Rg"
    Cn = "Cn"
    Nh = "Nh"
    Fl = "Fl"
    Mc = "Mc"
    Lv = "Lv"
    Ts = "Ts"
    Og = "Og"


@functools.total_ordering
class Species(MSONable, Stringify):
    """
    An extension of Element with an oxidation state and other optional
    properties. Properties associated with Species should be "idealized"
    values, not calculated values. For example, high-spin Fe2+ may be
    assigned an idealized spin of +5, but an actual Fe2+ site may be
    calculated to have a magmom of +4.5. Calculated properties should be
    assigned to Site objects, and not Species.
    """

    STRING_MODE = "SUPERSCRIPT"
    supported_properties = ("spin",)

    def __init__(
        self,
        symbol: str,
        oxidation_state: float | None = 0.0,
        properties: dict | None = None,
    ):
        """
        Initializes a Species.

        Args:
            symbol (str): Element symbol, e.g., Fe
            oxidation_state (float): Oxidation state of element, e.g., 2 or -2
            properties: Properties associated with the Species, e.g.,
                {"spin": 5}. Defaults to None. Properties must be one of the
                Species supported_properties.

        .. attribute:: oxi_state

            Oxidation state associated with Species

        .. attribute:: ionic_radius

            Ionic radius of Species (with specific oxidation state).

        .. versionchanged:: 2.6.7

            Properties are now checked when comparing two Species for equality.
        """
        self._el = Element(symbol)
        self._oxi_state = oxidation_state
        self._properties = properties or {}
        for k, _ in self._properties.items():
            if k not in Species.supported_properties:
                raise ValueError(f"{k} is not a supported property")

    def __getattr__(self, a):
        # overriding getattr doesn't play nice with pickle, so we
        # can't use self._properties
        p = object.__getattribute__(self, "_properties")
        if a in p:
            return p[a]
        return getattr(self._el, a)

    def __eq__(self, other: object) -> bool:
        """
        Species is equal to other only if element and oxidation states are exactly the same.
        """
        if not hasattr(other, "oxi_state") or not hasattr(other, "symbol") or not hasattr(other, "_properties"):
            return NotImplemented

        return all(getattr(self, attr) == getattr(other, attr) for attr in ["symbol", "oxi_state", "_properties"])

    def __hash__(self):
        """
        Equal Species should have the same str representation, hence
        should hash equally. Unequal Species will have different str
        representations.
        """
        return hash(str(self))

    def __lt__(self, other: object) -> bool:
        """
        Sets a default sort order for atomic species by electronegativity,
        followed by oxidation state, followed by spin.
        """
        if not isinstance(other, type(self)):
            return NotImplemented

        x1 = float("inf") if self.X != self.X else self.X
        x2 = float("inf") if other.X != other.X else other.X
        if x1 != x2:
            return x1 < x2
        if self.symbol != other.symbol:
            # There are cases where the electronegativity are exactly equal.
            # We then sort by symbol.
            return self.symbol < other.symbol
        if self.oxi_state:
            other_oxi = 0 if (isinstance(other, Element) or other.oxi_state is None) else other.oxi_state
            return self.oxi_state < other_oxi
        if getattr(self, "spin", False):
            other_spin = getattr(other, "spin", 0)
            return self.spin < other_spin
        return False

    @property
    def element(self):
        """
        Underlying element object
        """
        return self._el

    @property
    def ionic_radius(self) -> float | None:
        """
        Ionic radius of specie. Returns None if data is not present.
        """
        if self._oxi_state in self.ionic_radii:
            return self.ionic_radii[self._oxi_state]
        if self._oxi_state:
            d = self._el.data
            oxstr = str(int(self._oxi_state))
            if oxstr in d.get("Ionic radii hs", {}):
                warnings.warn(f"No default ionic radius for {self}. Using hs data.")
                return d["Ionic radii hs"][oxstr]
            if oxstr in d.get("Ionic radii ls", {}):
                warnings.warn(f"No default ionic radius for {self}. Using ls data.")
                return d["Ionic radii ls"][oxstr]
        warnings.warn(f"No ionic radius for {self}!")
        return None

    @property
    def oxi_state(self) -> float | None:
        """
        Oxidation state of Species.
        """
        return self._oxi_state

    @staticmethod
    def from_string(species_string: str) -> Species:
        """
        Returns a Species from a string representation.

        Args:
            species_string (str): A typical string representation of a
                species, e.g., "Mn2+", "Fe3+", "O2-".

        Returns:
            A Species object.

        Raises:
            ValueError if species_string cannot be interpreted.
        """
        # e.g. Fe2+,spin=5
        # 1st group: ([A-Z][a-z]*)    --> Fe
        # 2nd group: ([0-9.]*)        --> "2"
        # 3rd group: ([+\-])          --> +
        # 4th group: (.*)             --> everything else, ",spin=5"

        m = re.search(r"([A-Z][a-z]*)([0-9.]*)([+\-]*)(.*)", species_string)
        if m:

            # parse symbol
            sym = m.group(1)

            # parse oxidation state (optional)
            if not m.group(2) and not m.group(3):
                oxi = None
            else:
                oxi = 1 if m.group(2) == "" else float(m.group(2))
                oxi = -oxi if m.group(3) == "-" else oxi

            # parse properties (optional)
            properties = None
            if m.group(4):
                toks = m.group(4).replace(",", "").split("=")
                properties = {toks[0]: ast.literal_eval(toks[1])}

            # but we need either an oxidation state or a property
            if oxi is None and properties is None:
                raise ValueError("Invalid Species String")

            return Species(sym, 0 if oxi is None else oxi, properties)
        raise ValueError("Invalid Species String")

    def __repr__(self):
        return f"Species {self}"

    def __str__(self):
        output = self.symbol
        if self.oxi_state is not None:
            if self.oxi_state >= 0:
                output += formula_double_format(self.oxi_state) + "+"
            else:
                output += formula_double_format(-self.oxi_state) + "-"
        for p, v in self._properties.items():
            output += f",{p}={v}"
        return output

    def to_pretty_string(self) -> str:
        """
        :return: String without properties.
        """
        output = self.symbol
        if self.oxi_state is not None:
            if self.oxi_state >= 0:
                output += formula_double_format(self.oxi_state) + "+"
            else:
                output += formula_double_format(-self.oxi_state) + "-"
        return output

    def get_nmr_quadrupole_moment(self, isotope: str | None = None) -> float:
        """
        Gets the nuclear electric quadrupole moment in units of
        e*millibarns

        Args:
            isotope (str): the isotope to get the quadrupole moment for
                default is None, which gets the lowest mass isotope
        """
        quad_mom = self._el.nmr_quadrupole_moment

        if not quad_mom:
            return 0.0

        if isotope is None:
            isotopes = list(quad_mom)
            isotopes.sort(key=lambda x: int(x.split("-")[1]), reverse=False)
            return quad_mom.get(isotopes[0], 0.0)

        if isotope not in quad_mom:
            raise ValueError(f"No quadrupole moment for isotope {isotope}")
        return quad_mom.get(isotope, 0.0)

    def get_shannon_radius(
        self,
        cn: str,
        spin: Literal["", "Low Spin", "High Spin"] = "",
        radius_type: Literal["ionic", "crystal"] = "ionic",
    ) -> float:
        """
        Get the local environment specific ionic radius for species.

        Args:
            cn (str): Coordination using roman letters. Supported values are
                I-IX, as well as IIIPY, IVPY and IVSQ.
            spin (str): Some species have different radii for different
                spins. You can get specific values using "High Spin" or
                "Low Spin". Leave it as "" if not available. If only one spin
                data is available, it is returned and this spin parameter is
                ignored.
            radius_type (str): Either "crystal" or "ionic" (default).

        Returns:
            Shannon radius for specie in the specified environment.
        """
        radii = self._el.data["Shannon radii"]
        radii = radii[str(int(self._oxi_state))][cn]  # type: ignore
        if len(radii) == 1:
            k, data = list(radii.items())[0]
            if k != spin:
                warnings.warn(
                    f"Specified spin state of {spin} not consistent with database "
                    f"spin of {k}. Only one spin data available, and that value is returned."
                )
        else:
            data = radii[spin]
        return data[f"{radius_type}_radius"]

    def get_crystal_field_spin(
        self, coordination: Literal["oct", "tet"] = "oct", spin_config: Literal["low", "high"] = "high"
    ) -> float:
        """
        Calculate the crystal field spin based on coordination and spin
        configuration. Only works for transition metal species.

        Args:
            coordination ("oct" | "tet"): Tetrahedron or octahedron crystal site coordination
            spin_config ("low" | "high"): Whether the species is in a high or low spin state

        Returns:
            Crystal field spin in Bohr magneton.

        Raises:
            AttributeError if species is not a valid transition metal or has
                an invalid oxidation state.
            ValueError if invalid coordination or spin_config.
        """
        if coordination not in ("oct", "tet") or spin_config not in ("high", "low"):
            raise ValueError("Invalid coordination or spin config.")
        elec = self.full_electronic_structure
        if len(elec) < 4 or elec[-1][1] != "s" or elec[-2][1] != "d":
            raise AttributeError(f"Invalid element {self.symbol} for crystal field calculation.")
        nelectrons = elec[-1][2] + elec[-2][2] - self.oxi_state
        if nelectrons < 0 or nelectrons > 10:
            raise AttributeError(f"Invalid oxidation state {self.oxi_state} for element {self.symbol}")
        if spin_config == "high":
            if nelectrons <= 5:
                return nelectrons
            return 10 - nelectrons
        if spin_config == "low":
            if coordination == "oct":
                if nelectrons <= 3:
                    return nelectrons
                if nelectrons <= 6:
                    return 6 - nelectrons
                if nelectrons <= 8:
                    return nelectrons - 6
                return 10 - nelectrons
            if coordination == "tet":
                if nelectrons <= 2:
                    return nelectrons
                if nelectrons <= 4:
                    return 4 - nelectrons
                if nelectrons <= 7:
                    return nelectrons - 4
                return 10 - nelectrons
        raise RuntimeError()

    def __deepcopy__(self, memo):
        return Species(self.symbol, self.oxi_state, self._properties)

    def as_dict(self) -> dict:
        """
        :return: Json-able dictionary representation.
        """
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "element": self.symbol,
            "oxidation_state": self._oxi_state,
        }
        if self._properties:
            d["properties"] = self._properties
        return d

    @classmethod
    def from_dict(cls, d) -> Species:
        """
        :param d: Dict representation.
        :return: Species.
        """
        return cls(d["element"], d["oxidation_state"], d.get("properties", None))


@functools.total_ordering
class DummySpecies(Species):
    """
    A special specie for representing non-traditional elements or species. For
    example, representation of vacancies (charged or otherwise), or special
    sites, etc.

    .. attribute:: oxi_state

        Oxidation state associated with Species.

    .. attribute:: Z

        DummySpecies is always assigned an atomic number equal to the hash
        number of the symbol. Obviously, it makes no sense whatsoever to use
        the atomic number of a Dummy specie for anything scientific. The purpose
        of this is to ensure that for most use cases, a DummySpecies behaves no
        differently from an Element or Species.

    .. attribute:: X

        DummySpecies is always assigned an electronegativity of 0.
    """

    def __init__(
        self,
        symbol: str = "X",
        oxidation_state: float | None = 0,
        properties: dict | None = None,
    ):
        """
        Args:
            symbol (str): An assigned symbol for the dummy specie. Strict
                rules are applied to the choice of the symbol. The dummy
                symbol cannot have any part of first two letters that will
                constitute an Element symbol. Otherwise, a composition may
                be parsed wrongly. E.g., "X" is fine, but "Vac" is not
                because Vac contains V, a valid Element.
            oxidation_state (float): Oxidation state for dummy specie.
                Defaults to zero.
        """
        # enforce title case to match other elements, reduces confusion
        # when multiple DummySpecies in a "formula" string
        symbol = symbol.title()

        for i in range(1, min(2, len(symbol)) + 1):
            if Element.is_valid_symbol(symbol[:i]):
                raise ValueError(f"{symbol} contains {symbol[:i]}, which is a valid element symbol.")

        # Set required attributes for DummySpecies to function like a Species in
        # most instances.
        self._symbol = symbol
        self._oxi_state = oxidation_state
        self._properties = properties or {}
        for k, _ in self._properties.items():
            if k not in Species.supported_properties:
                raise ValueError(f"{k} is not a supported property")

    def __getattr__(self, a):
        # overriding getattr doesn't play nice with pickle, so we
        # can't use self._properties
        p = object.__getattribute__(self, "_properties")
        if a in p:
            return p[a]
        raise AttributeError(a)

    def __lt__(self, other):
        """
        Sets a default sort order for atomic species by electronegativity,
        followed by oxidation state.
        """
        if self.X != other.X:
            return self.X < other.X
        if self.symbol != other.symbol:
            # There are cases where the electronegativity are exactly equal.
            # We then sort by symbol.
            return self.symbol < other.symbol
        other_oxi = 0 if isinstance(other, Element) else other.oxi_state
        return self.oxi_state < other_oxi

    @property
    def Z(self) -> int:
        """
        DummySpecies is always assigned an atomic number equal to the hash of
        the symbol. The expectation is that someone would be an actual dummy
        to use atomic numbers for a Dummy specie.
        """
        return hash(self.symbol)

    @property
    def oxi_state(self) -> float | None:
        """
        Oxidation state associated with DummySpecies
        """
        return self._oxi_state

    @property
    def X(self) -> float:
        """
        DummySpecies is always assigned an electronegativity of 0. The effect of
        this is that DummySpecies are always sorted in front of actual Species.
        """
        return 0.0

    @property
    def symbol(self) -> str:
        """
        :return: Symbol for DummySpecies.
        """
        return self._symbol

    def __deepcopy__(self, memo):
        return DummySpecies(self.symbol, self._oxi_state)

    @staticmethod
    def from_string(species_string: str) -> DummySpecies:
        """
        Returns a Dummy from a string representation.

        Args:
            species_string (str): A string representation of a dummy
                species, e.g., "X2+", "X3+".

        Returns:
            A DummySpecies object.

        Raises:
            ValueError if species_string cannot be interpreted.
        """
        m = re.search(r"([A-ZAa-z]*)([0-9.]*)([+\-]*)(.*)", species_string)
        if m:
            sym = m.group(1)
            if m.group(2) == "" and m.group(3) == "":
                oxi = 0.0
            else:
                oxi = 1.0 if m.group(2) == "" else float(m.group(2))
                oxi = -oxi if m.group(3) == "-" else oxi
            properties = None
            if m.group(4):
                toks = m.group(4).split("=")
                properties = {toks[0]: float(toks[1])}
            return DummySpecies(sym, oxi, properties)
        raise ValueError("Invalid DummySpecies String")

    def as_dict(self) -> dict:
        """
        :return: MSONAble dict representation.
        """
        d = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "element": self.symbol,
            "oxidation_state": self._oxi_state,
        }
        if self._properties:
            d["properties"] = self._properties  # type: ignore
        return d

    @classmethod
    def from_dict(cls, d) -> DummySpecies:
        """
        :param d: Dict representation
        :return: DummySpecies
        """
        return cls(d["element"], d["oxidation_state"], d.get("properties", None))

    def __repr__(self):
        return f"DummySpecies {self}"

    def __str__(self):
        output = self.symbol
        if self.oxi_state is not None:
            if self.oxi_state >= 0:
                output += formula_double_format(self.oxi_state) + "+"
            else:
                output += formula_double_format(-self.oxi_state) + "-"
        for p, v in self._properties.items():
            output += f",{p}={v}"
        return output


@functools.total_ordering
class Specie(Species):
    """
    This maps the historical grammatically inaccurate Specie to Species
    to maintain backwards compatibility.
    """


@functools.total_ordering
class DummySpecie(DummySpecies):
    """
    This maps the historical grammatically inaccurate DummySpecie to DummySpecies
    to maintain backwards compatibility.
    """


def get_el_sp(obj) -> Element | Species | DummySpecies:
    """
    Utility method to get an Element or Species from an input obj.
    If obj is in itself an element or a specie, it is returned automatically.
    If obj is an int or a string representing an integer, the Element
    with the atomic number obj is returned.
    If obj is a string, Species parsing will be attempted (e.g., Mn2+), failing
    which Element parsing will be attempted (e.g., Mn), failing which
    DummyElement parsing will be attempted.

    Args:
        obj (Element/Species/str/int): An arbitrary object. Supported objects
            are actual Element/Species objects, integers (representing atomic
            numbers) or strings (element symbols or species strings).

    Returns:
        Species or Element, with a bias for the maximum number of properties
        that can be determined.

    Raises:
        ValueError if obj cannot be converted into an Element or Species.
    """
    if isinstance(obj, (Element, Species, DummySpecies)):
        return obj

    try:
        c = float(obj)
        i = int(c)
        i = i if i == c else None  # type: ignore
    except (ValueError, TypeError):
        i = None

    if i is not None:
        return Element.from_Z(i)

    try:
        return Species.from_string(obj)
    except (ValueError, KeyError):
        try:
            return Element(obj)
        except (ValueError, KeyError):
            try:
                return DummySpecies.from_string(obj)
            except Exception:
                raise ValueError(f"Can't parse Element or String from type {type(obj)}: {obj}.")
