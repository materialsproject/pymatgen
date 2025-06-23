"""Classes representing:
- Element: Element in the periodic table.
- Species: Element with optional oxidation state and spin.
- DummySpecies: Non-traditional Elements/Species (vacancies/...).
- ElementType: element types (metal/noble_gas/halogen/s_block/...).
"""

from __future__ import annotations

import ast
import functools
import gzip
import re
import warnings
from collections import Counter
from enum import Enum, EnumMeta, unique
from itertools import combinations, product
from pathlib import Path
from typing import TYPE_CHECKING, overload

import numpy as np
import orjson
from monty.json import MSONable

from pymatgen.core.units import SUPPORTED_UNIT_NAMES, FloatWithUnit, Ha_to_eV, Length, Mass, Unit
from pymatgen.io.core import ParseError
from pymatgen.util.string import Stringify, formula_double_format

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal

    from typing_extensions import Self

    from pymatgen.util.typing import SpeciesLike

# Load element data (periodic table) from JSON file
# NOTE: you should not update the JSON file manually,
# see `dev_scripts/generate_periodic_table_yaml_json.py`
with gzip.open(Path(__file__).absolute().parent / "periodic_table.json.gz", mode="rb") as ptable_json:
    _RAW_PT_DATA: dict = orjson.loads(ptable_json.read())

_PT_UNIT: dict[str, str] = _RAW_PT_DATA.pop("_unit")
_PT_DATA: dict[str, Any] = _RAW_PT_DATA

_PT_ROW_SIZES: tuple[int, ...] = (2, 8, 8, 18, 18, 32, 32)

# Madelung energy ordering rule (lower to higher energy)
_MADELUNG: list[tuple[int, str]] = [
    (1, "s"),
    (2, "s"),
    (2, "p"),
    (3, "s"),
    (3, "p"),
    (4, "s"),
    (3, "d"),
    (4, "p"),
    (5, "s"),
    (4, "d"),
    (5, "p"),
    (6, "s"),
    (4, "f"),
    (5, "d"),
    (6, "p"),
    (7, "s"),
    (5, "f"),
    (6, "d"),
    (7, "p"),
]


@functools.total_ordering
@unique
class ElementBase(Enum):
    """Element class defined without any enum values so it can be subclassed."""

    def __init__(self, symbol: str) -> None:
        """
        This class provides a basic, immutable representation of an element, including
        all relevant chemical and physical properties. It ensures that elements are
        handled as singletons, reducing redundancy and improving efficiency. Missing
        data is represented by `None` unless otherwise specified.

        Attributes:
            Z (int): Atomic number of the element.
            symbol (str): Element symbol (e.g., "H", "Fe").
            long_name (str): Full name of the element (e.g., "Hydrogen").
            A (int, optional): Atomic mass number (sum of protons and neutrons).
            atomic_radius_calculated (float, optional): Calculated atomic radius (Å).
            van_der_waals_radius (float, optional): Van der Waals radius (Å).
            mendeleev_no (int, optional): Mendeleev number based on crystal-structure maps.
            electrical_resistivity (float, optional): Electrical resistivity (Ω·m).
            velocity_of_sound (float, optional): Velocity of sound (m/s).
            reflectivity (float, optional): Reflectivity (%).
            refractive_index (float, optional): Refractive index.
            poissons_ratio (float, optional): Poisson's ratio.
            molar_volume (float, optional): Molar volume (cm³/mol).
            electronic_structure (str): Electronic structure (e.g., "[Ar].3d6.4s2").
            atomic_orbitals (dict): Orbital energies (Hartree units).
            atomic_orbitals_eV (dict): Orbital energies in electron volts (eV).
            thermal_conductivity (float, optional): Thermal conductivity (W/m·K).
            boiling_point (float, optional): Boiling point (K).
            melting_point (float, optional): Melting point (K).
            critical_temperature (float, optional): Critical temperature (K).
            superconduction_temperature (float, optional): Superconducting transition temperature (K).
            liquid_range (float, optional): Temperature range for liquid phase (K).
            bulk_modulus (float, optional): Bulk modulus (GPa).
            youngs_modulus (float, optional): Young's modulus (GPa).
            brinell_hardness (float, optional): Brinell hardness (MPa).
            rigidity_modulus (float, optional): Rigidity modulus (GPa).
            mineral_hardness (float, optional): Mohs hardness.
            vickers_hardness (float, optional): Vickers hardness (MPa).
            density_of_solid (float, optional): Density in solid phase (g/cm³).
            coefficient_of_linear_thermal_expansion (float, optional): Thermal expansion coefficient (K⁻¹).
            ground_level (float, optional): Ground energy level of the element.
            ionization_energies (list[Optional[float]]): Ionization energies (kJ/mol), indexed from 0.

        Examples:
            Create an element instance and access its properties:
                >>> hydrogen = Element("H")
                >>> hydrogen.symbol
                'H'
                >>> hydrogen.Z
                1
                >>> hydrogen.electronic_structure
                '1s1'

            Access additional attributes such as atomic radius:
                >>> hydrogen.atomic_radius_calculated
                0.53

        Notes:
            - This class supports handling of isotopes by incorporating named isotopes
                and their respective properties.
            - Attributes are populated using a JSON file that stores data about all
                known elements.
            - Some attributes are calculated or derived based on predefined constants
                and rules.
        """
        self.symbol = symbol
        data = _PT_DATA[symbol]

        # Store key variables for quick access
        self.Z = data["Atomic no"]

        self._is_named_isotope = data.get("Is named isotope", False)
        if self._is_named_isotope:
            for sym, info in _PT_DATA.items():
                if info["Atomic no"] == self.Z and not info.get("Is named isotope", False):
                    self.symbol = sym
                    break
            # For specified/named isotopes, treat the same as named element
            # (the most common isotope). Then we pad the data block with the
            # entries for the named element.
            data = {**_PT_DATA[self.symbol], **data}

        at_r: float | None = data.get("Atomic radius")
        self._atomic_radius = None if at_r is None else Length(at_r, _PT_UNIT["Atomic radius"])

        self._atomic_mass = Mass(data["Atomic mass"], _PT_UNIT["Atomic mass"])

        self._atomic_mass_number = None
        self.A = data.get("Atomic mass no")
        if self.A is not None:
            self._atomic_mass_number = Mass(self.A, _PT_UNIT["Atomic mass no"])

        self.long_name = data["Name"]
        self._data = data

    def __getattr__(self, item: str) -> Any:
        """Key access to available element data.

        Args:
            item (str): Attribute name.

        Raises:
            AttributeError: If item not in _PT_DATA.
        """
        if item not in {
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
            "metallic_radius",
        }:
            raise AttributeError(f"Element has no attribute {item}!")

        prop_name: str = item.capitalize().replace("_", " ")
        val: Any = self._data.get(prop_name)

        if val is None:
            warnings.warn(f"No data available for {item} for {self.symbol}", stacklevel=2)
            return None

        if isinstance(val, list | dict):
            return val

        unit: str | None = _PT_UNIT.get(prop_name)

        if unit is not None:
            if unit in SUPPORTED_UNIT_NAMES:
                return FloatWithUnit(float(val), Unit(unit))
            return FloatWithUnit(float(val), unit)

        try:
            return float(val)
        except ValueError:
            return val

    def __eq__(self, other: object) -> bool:
        return isinstance(self, Element) and isinstance(other, Element) and self.Z == other.Z and self.A == other.A

    def __hash__(self) -> int:
        # multiply Z by 1000 to avoid hash collisions of element N with isotopes of elements N+/-1,2,3...
        return self.Z * 1000 + self.A if self._is_named_isotope else self.Z * 137 * 100

    def __repr__(self) -> str:
        return f"Element {self.symbol}"

    def __str__(self) -> str:
        return self.symbol

    def __lt__(self, other):
        """Set a default sort order for atomic species by Pauling electronegativity. Very
        useful for getting correct formulas. For example, FeO4PLi is
        automatically sorted into LiFePO4.
        """
        if not hasattr(other, "X") or not hasattr(other, "symbol"):
            return NotImplemented
        x1 = float("inf") if self.X != self.X else self.X
        x2 = float("inf") if other.X != other.X else other.X
        if x1 != x2:
            return x1 < x2

        # There are cases where the Pauling electronegativity are exactly equal.
        # We then sort by symbol.
        return self.symbol < other.symbol

    def __deepcopy__(self, memo) -> Element:
        return Element(self.symbol)

    @property
    def X(self) -> float:
        """Pauling electronegativity of element. Note that if an element does not
        have an Pauling electronegativity, a NaN float is returned.
        """
        if X := self._data.get("X"):
            return X
        warnings.warn(
            f"No Pauling electronegativity for {self.symbol}. Setting to NaN. This has no physical meaning, "
            "and is mainly done to avoid errors caused by the code expecting a float.",
            stacklevel=2,
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

    @property
    def atomic_mass_number(self) -> FloatWithUnit | None:
        """
        Returns:
            float: The atomic mass of the element in amu.
        """
        return self._atomic_mass_number

    @property
    def atomic_orbitals_eV(self) -> dict[str, float]:
        """The LDA energies in eV for neutral atoms, by orbital.

        This property contains the same info as `self.atomic_orbitals`,
        but uses eV for units, per matsci issue https://matsci.org/t/unit-of-atomic-orbitals-energy/54325
        In short, self.atomic_orbitals was meant to be in eV all along but is now kept
        as Hartree for backwards compatibility.
        """
        return {orb_idx: energy * Ha_to_eV for orb_idx, energy in self.atomic_orbitals.items()}

    @property
    def data(self) -> dict[str, Any]:
        """Dict of data for element."""
        return self._data.copy()

    @property
    def ionization_energy(self) -> float | None:
        """First ionization energy of element."""
        if self.ionization_energies is None:
            warnings.warn(f"No data available for ionization_energy for {self.symbol}", stacklevel=2)
            return None
        return self.ionization_energies[0]

    @property
    def electron_affinity(self) -> float:
        """The amount of energy released when an electron is attached to a neutral atom."""
        return self._data["Electron affinity"]

    @property
    def electronic_structure(self) -> str:
        """Electronic structure as string, with only valence electrons. The
        electrons are listed in order of increasing prinicpal quantum number
         (orbital number), irrespective of the actual energy level,
        e.g., The electronic structure for Fe is represented as '[Ar].3d6.4s2'
        even though the 3d electrons are higher in energy than the 4s.
        """
        return self._data["Electronic structure"]["0"]

    @property
    def average_ionic_radius(self) -> FloatWithUnit:
        """Average ionic radius for element (with units). The average is taken
        over all oxidation states of the element for which data is present.
        """
        if "Ionic radii" in self._data:
            radii = self._data["Ionic radii"]
            radius = sum(radii.values()) / len(radii)
        else:
            radius = 0.0
        return FloatWithUnit(radius, _PT_UNIT["Ionic radii"])

    @property
    def average_cationic_radius(self) -> FloatWithUnit:
        """Average cationic radius for element (with units). The average is
        taken over all positive oxidation states of the element for which
        data is present.
        """
        if "Ionic radii" in self._data and (radii := [v for k, v in self._data["Ionic radii"].items() if int(k) > 0]):
            return FloatWithUnit(sum(radii) / len(radii), _PT_UNIT["Ionic radii"])
        return FloatWithUnit(0.0, _PT_UNIT["Ionic radii"])

    @property
    def average_anionic_radius(self) -> FloatWithUnit:
        """Average anionic radius for element (with units). The average is
        taken over all negative oxidation states of the element for which
        data is present.
        """
        if "Ionic radii" in self._data and (radii := [v for k, v in self._data["Ionic radii"].items() if int(k) < 0]):
            return FloatWithUnit(sum(radii) / len(radii), _PT_UNIT["Ionic radii"])
        return FloatWithUnit(0.0, _PT_UNIT["Ionic radii"])

    @property
    def ionic_radii(self) -> dict[int, FloatWithUnit]:
        """All ionic radii of the element as a dict of
        {oxidation state: ionic radii}. Radii are given in angstrom.
        """
        if "Ionic radii" in self._data:
            return {int(k): FloatWithUnit(v, _PT_UNIT["Ionic radii"]) for k, v in self._data["Ionic radii"].items()}
        return {}

    @property
    def number(self) -> int:
        """Alternative attribute for atomic number Z."""
        return self.Z

    @property
    def max_oxidation_state(self) -> float:
        """Maximum oxidation state for element."""
        if "Oxidation states" in self._data:
            return max(self._data["Oxidation states"])
        return 0

    @property
    def min_oxidation_state(self) -> float:
        """Minimum oxidation state for element."""
        if "Oxidation states" in self._data:
            return min(self._data["Oxidation states"])
        return 0

    @property
    def oxidation_states(self) -> tuple[int, ...]:
        """Tuple of all known oxidation states."""
        return tuple(map(int, self._data.get("Oxidation states", ())))

    @property
    def common_oxidation_states(self) -> tuple[int, ...]:
        """Tuple of common oxidation states."""
        return tuple(self._data.get("Common oxidation states", ()))

    @property
    def icsd_oxidation_states(self) -> tuple[int, ...]:
        """Tuple of all oxidation states with at least 10 instances in
        ICSD database AND at least 1% of entries for that element.
        """
        return tuple(self._data.get("ICSD oxidation states", ()))

    @property
    def full_electronic_structure(self) -> list[tuple[int, str, int]]:
        """Full electronic structure in order of increasing
        energy level (according to the Madelung rule). Therefore, the final
        element in the list gives the electronic structure of the valence shell.

        For example, the full electronic structure for Fe is:
            [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6),
            (4, "s", 2), (3, "d", 6)].

        Returns:
            list[tuple[int, str, int]]: A list of tuples representing each subshell,
            where each tuple contains:
                - `n` (int): Principal quantum number.
                - `orbital_type` (str): Orbital type (e.g., "s", "p", "d", "f").
                - `electron_count` (int): Number of electrons in the subshell.
        """
        e_str: str = self.electronic_structure

        def parse_orbital(orb_str: str) -> str | tuple[int, str, int]:
            """Parse orbital information from split electron configuration string."""
            # Parse valence subshell notation (e.g., "3d6" -> (3, "d", 6))
            if match := re.match(r"(\d+)([spdfg]+)(\d+)", orb_str):
                return int(match[1]), match[2], int(match[3])

            # Return core-electron configuration as-is (e.g. "[Ar]")
            return orb_str

        # Split e_str (e.g. for Fe "[Ar].3d6.4s2" into ["[Ar]", "3d6", "4s2"])
        data: list = [parse_orbital(s) for s in e_str.split(".")]

        # Fully expand core-electron configuration (replace noble gas notation string)
        if isinstance(data[0], str):
            sym: str = data[0].replace("[", "").replace("]", "")
            data = list(Element(sym).full_electronic_structure) + data[1:]

        # Sort the final electronic structure by increasing energy level
        return sorted(data, key=lambda x: _MADELUNG.index((x[0], x[1])))

    @property
    def n_electrons(self) -> int:
        """Total number of electrons in the Element."""
        return sum(t[-1] for t in self.full_electronic_structure)

    @property
    def valence(self) -> tuple[int | float, int]:
        """Valence subshell angular moment (L) and number of valence e- (v_e),
        obtained from full electron config, where L=0, 1, 2, or 3 for s, p, d,
        and f orbitals, respectively.
        """
        if self.group == 18:
            return float("nan"), 0  # The number of valence of noble gas is 0

        L_symbols = "SPDFGHIKLMNOQRTUVWXYZ"
        valence: list[tuple[int, int]] = []
        full_electron_config = self.full_electronic_structure
        last_orbital = full_electron_config[-1]
        for n, l_symbol, ne in full_electron_config:
            idx = L_symbols.lower().index(l_symbol)
            if ne < (2 * idx + 1) * 2 or (
                (n, l_symbol, ne) == last_orbital and ne == (2 * idx + 1) * 2 and len(valence) == 0
            ):  # check for full last shell (e.g. column 2)
                valence.append((idx, ne))
        if len(valence) > 1:
            raise ValueError(f"{self} has ambiguous valence")

        return valence[0]

    @property
    def term_symbols(self) -> list[list[str]]:
        """All possible Russell-Saunders term symbol of the Element.
        eg. L = 1, n_e = 2 (s2) returns [['1D2'], ['3P0', '3P1', '3P2'], ['1S0']].
        """
        if self.is_noble_gas:
            return [["1S0"]]

        L, v_e = self.valence

        # for one electron in subshell L
        ml = list(range(-L, L + 1))  # type: ignore[arg-type]
        ms = [1 / 2, -1 / 2]
        # all possible configurations of ml,ms for one e in subshell L
        ml_ms = list(product(ml, ms))

        # Number of possible configurations for r electrons in subshell L.
        n = (2 * L + 1) * 2
        # the combination of n_e electrons configurations
        # C^{n}_{n_e}
        e_config_combs = list(combinations(range(n), v_e))  # type: ignore[arg-type]

        # Total ML = sum(ml1, ml2), Total MS = sum(ms1, ms2)
        TL = [sum(ml_ms[comb[e]][0] for e in range(v_e)) for comb in e_config_combs]
        TS = [sum(ml_ms[comb[e]][1] for e in range(v_e)) for comb in e_config_combs]
        comb_counter = Counter(zip(TL, TS, strict=True))

        term_symbols = []
        L_symbols = "SPDFGHIKLMNOQRTUVWXYZ"
        while sum(comb_counter.values()) > 0:
            # Start from the lowest freq combination,
            # which corresponds to largest abs(L) and smallest abs(S)
            L, S = min(comb_counter)

            J = list(np.arange(abs(L - S), abs(L) + abs(S) + 1))
            term_symbols.append([f"{int(2 * (abs(S)) + 1)}{L_symbols[abs(L)]}{j}" for j in J])

            # Delete all configurations included in this term
            for ML in range(-L, L - 1, -1):
                for MS in np.arange(S, -S + 1, 1):
                    if (ML, MS) in comb_counter:
                        comb_counter[ML, MS] -= 1  # type:ignore[index]
                        if comb_counter[ML, MS] == 0:  # type:ignore[index]
                            del comb_counter[ML, MS]
        return term_symbols

    @property
    def ground_state_term_symbol(self) -> str:
        """Ground state term symbol, selected based on Hund's Rule."""
        L_symbols = "SPDFGHIKLMNOQRTUVWXYZ"

        term_symbols = self.term_symbols
        term_symbol_flat: dict = {
            term: {
                "multiplicity": int(term[0]),
                "L": L_symbols.index(term[1]),
                "J": float(term[2:]),
            }
            for term in [item for sublist in term_symbols for item in sublist]
        }

        multi = [int(item["multiplicity"]) for _terms, item in term_symbol_flat.items()]
        max_multi_terms = {
            symbol: item for symbol, item in term_symbol_flat.items() if item["multiplicity"] == max(multi)
        }

        Ls = [item["L"] for _terms, item in max_multi_terms.items()]
        max_L_terms = {symbol: item for symbol, item in term_symbol_flat.items() if item["L"] == max(Ls)}

        J_sorted_terms = sorted(max_L_terms.items(), key=lambda k: k[1]["J"])
        L, v_e = self.valence
        return J_sorted_terms[0][0] if v_e <= (2 * L + 1) else J_sorted_terms[-1][0]

    @staticmethod
    def from_Z(Z: int, A: int | None = None) -> Element:
        """Get an element from an atomic number.

        Args:
            Z (int): Atomic number (number of protons)
            A (int | None) : Atomic mass number (number of protons + neutrons)

        Returns:
            Element with atomic number Z.
        """
        for sym, data in _PT_DATA.items():
            atomic_mass_num = data.get("Atomic mass no") if A else None
            if data["Atomic no"] == Z and atomic_mass_num == A:
                return Element(sym)

        raise ValueError(f"Unexpected atomic number {Z=}")

    @staticmethod
    def from_name(name: str) -> Element:
        """Get an element from its long name.

        Args:
            name: Long name of the element, e.g. 'Hydrogen' or 'Iron'. Not case-sensitive.

        Returns:
            Element with the name 'name'
        """
        # Accommodate the British-English world
        uk_to_us = {"aluminium": "aluminum", "caesium": "cesium"}
        name = uk_to_us.get(name.lower(), name)

        for sym, data in _PT_DATA.items():
            if data["Name"] == name.capitalize():
                return Element(sym)

        raise ValueError(f"No element with the {name=}")

    @staticmethod
    def from_row_and_group(row: int, group: int) -> Element:
        """Get an element from a row and group number.

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

        Note:
            The 18 group number system is used, i.e. noble gases are group 18.
        """
        for sym in _PT_DATA:
            el = Element(sym)
            if 57 <= el.Z <= 71:
                el_pseudo_row = 8
                el_pseudo_group = (el.Z - 54) % 32
            elif 89 <= el.Z <= 103:
                el_pseudo_row = 9
                el_pseudo_group = (el.Z - 54) % 32
            else:
                el_pseudo_row = el.row
                el_pseudo_group = el.group
            if el_pseudo_row == row and el_pseudo_group == group:
                return el

        raise ValueError("No element with this row and group!")

    @staticmethod
    def is_valid_symbol(symbol: str) -> bool:
        """Check if symbol (e.g., "H") is a valid element symbol.

        Args:
            symbol (str): Element symbol

        Returns:
            bool: True if symbol is a valid element.
        """
        return symbol in Element.__members__

    @property
    def row(self) -> int:
        """The periodic table row of the element.

        Note: For lanthanoids and actinoids, the row is always 6 or 7,
        respectively.
        """
        z = self.Z
        total = 0
        if 57 <= z <= 71:
            return 6
        if 89 <= z <= 103:
            return 7
        for idx, size in enumerate(_PT_ROW_SIZES, start=1):
            total += size
            if total >= z:
                return idx
        return 8

    @property
    def group(self) -> int:
        """The periodic table group of the element.

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
    def block(self) -> Literal["s", "p", "d", "f"]:
        """The block character "s, p, d, f"."""
        if (self.is_actinoid or self.is_lanthanoid) and self.Z not in {71, 103}:
            return "f"
        if self.is_actinoid or self.is_lanthanoid or self.group in range(3, 13):
            return "d"
        if self.group in {1, 2}:
            return "s"
        if self.group in range(13, 19):
            return "p"

        raise ValueError("Unable to determine block.")

    @property
    def is_noble_gas(self) -> bool:
        """True if element is noble gas."""
        return self.Z in {2, 10, 18, 36, 54, 86, 118}

    @property
    def is_transition_metal(self) -> bool:
        """True if element is a transition metal."""
        ns = (*range(21, 31), *range(39, 49), 57, *range(72, 81), 89, *range(104, 113))
        return self.Z in ns

    @property
    def is_post_transition_metal(self) -> bool:
        """True if element is a post-transition or poor metal."""
        return self.symbol in {"Al", "Ga", "In", "Tl", "Sn", "Pb", "Bi"}

    @property
    def is_rare_earth(self) -> bool:
        """True if element is a rare earth element, including Lanthanides (La)
        series, Actinides (Ac) series, Scandium (Sc) and Yttrium (Y).
        """
        return self.is_lanthanoid or self.is_actinoid or self.symbol in {"Sc", "Y"}

    @property
    def is_metal(self) -> bool:
        """True if is a metal."""
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
        """True if element is a metalloid."""
        return self.symbol in {"B", "Si", "Ge", "As", "Sb", "Te", "Po"}

    @property
    def is_alkali(self) -> bool:
        """True if element is an alkali metal."""
        return self.Z in {3, 11, 19, 37, 55, 87}

    @property
    def is_alkaline(self) -> bool:
        """True if element is an alkaline earth metal (group II)."""
        return self.Z in {4, 12, 20, 38, 56, 88}

    @property
    def is_halogen(self) -> bool:
        """True if element is a halogen."""
        return self.Z in {9, 17, 35, 53, 85}

    @property
    def is_chalcogen(self) -> bool:
        """True if element is a chalcogen."""
        return self.Z in {8, 16, 34, 52, 84}

    @property
    def is_lanthanoid(self) -> bool:
        """True if element is a lanthanoid."""
        return 56 < self.Z < 72

    @property
    def is_actinoid(self) -> bool:
        """True if element is a actinoid."""
        return 88 < self.Z < 104

    @property
    def is_radioactive(self) -> bool:
        """True if element is radioactive."""
        return self.Z in {43, 61} or self.Z >= 84

    @property
    def is_quadrupolar(self) -> bool:
        """Check if this element can be quadrupolar."""
        return len(self.data.get("NMR Quadrupole Moment", {})) > 0

    @property
    def nmr_quadrupole_moment(self) -> dict[str, FloatWithUnit]:
        """A dictionary the nuclear electric quadrupole moment in units of
        e*millibarns for various isotopes.
        """
        return {
            k: FloatWithUnit(v, _PT_UNIT["NMR Quadrupole Moment"])
            for k, v in self.data.get("NMR Quadrupole Moment", {}).items()
        }

    @property
    def iupac_ordering(self) -> int:
        """Ordering according to Table VI of "Nomenclature of Inorganic Chemistry
        (IUPAC Recommendations 2005)". This ordering effectively follows the
        groups and rows of the periodic table, except the Lanthanides, Actinides
        and hydrogen.
        """
        return self._data["IUPAC ordering"]

    def as_dict(self) -> dict[str, Any]:
        """Serialize to MSONable dict representation e.g. to write to disk as JSON."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "element": self.symbol,
        }

    @staticmethod
    def from_dict(dct: dict) -> Element:
        """Deserialize from MSONable dict representation."""
        return Element(dct["element"])

    @staticmethod
    def print_periodic_table(filter_function: Callable | None = None) -> None:
        """A pretty ASCII printer for the periodic table, based on some
        filter_function.

        Args:
            filter_function: A filtering function taking an Element as input
                and returning a boolean. For example, setting
                filter_function = lambda el: el.X > 2 will print a periodic
                table containing only elements with Pauling electronegativity > 2.
        """
        for row in range(1, 10):
            row_str = []
            for group in range(1, 19):
                try:
                    el = Element.from_row_and_group(row, group)
                except ValueError:
                    el = None
                if el and ((not filter_function) or filter_function(el)):
                    row_str.append(f"{el.symbol:3s}")
                else:
                    row_str.append("   ")
            print(" ".join(row_str))


class _ElementMeta(EnumMeta):
    """Override Element to handle isotopes."""

    def __iter__(cls):
        """Skip named isotopes when iterating."""
        return (elem for elem in super().__iter__() if not elem._is_named_isotope)

    @property
    def named_isotopes(cls):
        """Get all named isotopes."""
        return tuple(elem for elem in super().__iter__() if elem._is_named_isotope)


class Element(ElementBase, metaclass=_ElementMeta):
    """Enum representing an element in the periodic table."""

    # This name = value convention is redundant and dumb, but unfortunately is
    # necessary to preserve backwards compatibility with a time when Element is
    # a regular object that is constructed with Element(symbol).
    H = "H"
    D = "D"
    T = "T"
    He = "He"
    Li = "Li"
    Be = "Be"
    B = "B"
    C = "C"
    N = "N"
    O = "O"  # noqa: E741
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
    I = "I"  # noqa: E741
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
    """An extension of Element with optional oxidation state and spin. Properties associated
    with Species should be "idealized" values, not calculated values. For example,
    high-spin Fe2+ may be assigned an idealized spin of +5, but an actual Fe2+ site may be
    calculated to have a magmom of +4.5. Calculated properties should be assigned to Site
    objects, and not Species.
    """

    STRING_MODE = "SUPERSCRIPT"

    def __init__(
        self,
        symbol: str,
        oxidation_state: float | None = None,
        spin: float | None = None,
    ) -> None:
        """
        Args:
            symbol (str): Element symbol optionally incl. oxidation state. E.g. Fe, Fe2+, O2-.
            oxidation_state (float): Explicit oxidation state of element, e.g. -2, -1, 0, 1, 2, ...
                If oxidation state is present in symbol, this argument is ignored.
            spin: Spin associated with Species. Defaults to None.

        Raises:
            ValueError: If oxidation state passed both in symbol string and via
                oxidation_state kwarg.
        """
        if oxidation_state is not None and symbol.endswith(("+", "-")):
            raise ValueError(
                f"Oxidation state should be specified either in {symbol=} or as {oxidation_state=}, not both."
            )
        if symbol.endswith(("+", "-")):
            # Extract oxidation state from symbol
            try:
                symbol, oxi = re.match(r"([A-Za-z]+)([0-9\.0-9]*[\+\-])", symbol).groups()  # type: ignore[union-attr]
                self._oxi_state: float | None = (1 if "+" in oxi else -1) * float(oxi[:-1] or 1)
            except AttributeError:
                raise ParseError(f"Failed to parse {symbol=}")
        else:
            self._oxi_state = oxidation_state

        self._el = Element(symbol)

        self._spin = spin

    def __getattr__(self, attr: str) -> Any:
        """Allow Species to inherit properties of underlying Element."""
        return getattr(self._el, attr)

    def __getstate__(self) -> dict:
        return self.__dict__

    def __setstate__(self, dct: dict) -> None:
        self.__dict__.update(dct)

    def __eq__(self, other: object) -> bool:
        """Species is equal to other only if element and oxidation states are exactly the same."""
        attrs = ("oxi_state", "symbol", "spin", "A")
        if not all(hasattr(other, attribute) for attribute in attrs):
            return NotImplemented

        return all(getattr(self, attr) == getattr(other, attr) for attr in attrs)

    def __hash__(self) -> int:
        """Equal Species should have the same str representation, hence
        should hash equally. Unequal Species will have different str
        representations.
        """
        return hash(str(self))

    def __lt__(self, other: object) -> bool:
        """Set a default sort order for atomic species by Pauling electronegativity,
        followed by oxidation state, followed by spin.
        """
        if not isinstance(other, type(self)):
            return NotImplemented

        x1 = float("inf") if self.X != self.X else self.X
        x2 = float("inf") if other.X != other.X else other.X
        if x1 != x2:
            return x1 < x2
        if self.symbol != other.symbol:
            # There are cases where the Pauling electronegativity are exactly equal.
            # We then sort by symbol.
            return self.symbol < other.symbol
        if self.oxi_state:
            other_oxi = 0 if (isinstance(other, Element) or other.oxi_state is None) else other.oxi_state
            return self.oxi_state < other_oxi
        if self.spin is not None:
            return self.spin < other.spin if other.spin is not None else False

        return False

    def __repr__(self) -> str:
        return f"Species {self}"

    def __str__(self) -> str:
        output = self.name if hasattr(self, "name") else self.symbol
        if self.oxi_state is not None:
            abs_charge = formula_double_format(abs(self.oxi_state))
            if isinstance(abs_charge, float):
                abs_charge = f"{abs_charge:.2f}"  # type: ignore[assignment]
            output += f"{abs_charge}{'+' if self.oxi_state >= 0 else '-'}"

        if self._spin is not None:
            spin = self._spin
            output += f",{spin=}"
        return output

    def __deepcopy__(self, memo) -> Self:
        return type(self)(self.symbol, self.oxi_state, spin=self._spin)

    @property
    def element(self) -> Element:
        """Underlying element object."""
        return self._el

    @property
    def oxi_state(self) -> float | None:
        """Oxidation state of Species."""
        return self._oxi_state

    @property
    def spin(self) -> float | None:
        """Spin of Species."""
        return self._spin

    @property
    def electronic_structure(self) -> str:
        """Electronic structure as string, with only valence electrons. The
        electrons are listed in order of increasing prinicpal quantum number
         (orbital number), irrespective of the actual energy level,
        e.g., The electronic structure for Fe is represented as '[Ar].3d6.4s2'
        even though the 3d electrons are higher in energy than the 4s.
        """
        if self._data["Electronic structure"].get(str(self._oxi_state)) is not None:
            return self._data["Electronic structure"][str(self._oxi_state)]

        raise ValueError(f"No electronic structure data for oxidation state {self.oxi_state}")

    # NOTE - copied exactly from Element. Refactoring / inheritance may improve
    # robustness
    @property
    def full_electronic_structure(self) -> list[tuple[int, str, int]]:
        """Full electronic structure in order of increasing
        energy level (according to the Madelung rule). Therefore, the final
        element in the list gives the electronic structure of the valence shell.

        For example, the full electronic structure for Fe is:
            [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6),
            (4, "s", 2), (3, "d", 6)].

        Returns:
            list[tuple[int, str, int]]: A list of tuples representing each subshell,
            where each tuple contains:
                - `n` (int): Principal quantum number.
                - `orbital_type` (str): Orbital type (e.g., "s", "p", "d", "f").
                - `electron_count` (int): Number of electrons in the subshell.
        """
        e_str: str = self.electronic_structure

        def parse_orbital(orb_str: str) -> str | tuple[int, str, int]:
            """Parse orbital information from split electron configuration string."""
            # Parse valence subshell notation (e.g., "3d6" -> (3, "d", 6))
            if match := re.match(r"(\d+)([spdfg]+)(\d+)", orb_str):
                return int(match[1]), match[2], int(match[3])

            # Return core-electron configuration as-is (e.g. "[Ar]")
            return orb_str

        data: list = [parse_orbital(s) for s in e_str.split(".")]
        if isinstance(data[0], str):
            sym = data[0].replace("[", "").replace("]", "")
            data = list(Element(sym).full_electronic_structure) + data[1:]

        # Sort the final electronic structure by increasing energy level
        return sorted(data, key=lambda x: _MADELUNG.index((x[0], x[1])))

    # NOTE - copied exactly from Element. Refactoring / inheritance may improve
    # robustness
    @property
    def n_electrons(self) -> int:
        """Total number of electrons in the Species."""
        return sum(t[-1] for t in self.full_electronic_structure)

    # NOTE - copied exactly from Element. Refactoring / inheritance may improve
    # robustness
    @property
    def valence(self) -> tuple[int | float, int]:
        """Valence subshell angular moment (L) and number of valence e- (v_e),
        obtained from full electron config, where L=0, 1, 2, or 3 for s, p, d,
        and f orbitals, respectively.
        """
        if self.group == 18:
            return float("nan"), 0  # The number of valence of noble gas is 0

        L_symbols = "SPDFGHIKLMNOQRTUVWXYZ"
        valence: list[tuple[int, int]] = []
        full_electron_config = self.full_electronic_structure
        last_orbital = full_electron_config[-1]
        for n, l_symbol, ne in full_electron_config:
            idx = L_symbols.lower().index(l_symbol)
            if ne < (2 * idx + 1) * 2 or (
                (n, l_symbol, ne) == last_orbital and ne == (2 * idx + 1) * 2 and len(valence) == 0
            ):  # check for full last shell (e.g. column 2)
                valence.append((idx, ne))
        if len(valence) > 1:
            raise ValueError(f"{self} has ambiguous valence")

        return valence[0]

    @property
    def ionic_radius(self) -> float | None:
        """Ionic radius of specie. Returns None if data is not present."""
        if self._oxi_state in self.ionic_radii:
            return self.ionic_radii[self._oxi_state]
        if self._oxi_state:
            dct = self._el.data
            oxi_str = str(int(self._oxi_state))
            warn_msg = f"No default ionic radius for {self}."
            if ion_rad := dct.get("Ionic radii hs", {}).get(oxi_str):
                warnings.warn(f"{warn_msg} Using hs data.", stacklevel=2)
                return ion_rad
            if ion_rad := dct.get("Ionic radii ls", {}).get(oxi_str):
                warnings.warn(f"{warn_msg} Using ls data.", stacklevel=2)
                return ion_rad
        warnings.warn(f"No ionic radius for {self}!", stacklevel=2)
        return None

    @classmethod
    def from_str(cls, species_string: str) -> Self:
        """Get a Species from a string representation.

        Args:
            species_string (str): A typical string representation of a
                species, e.g. "Mn2+", "Fe3+", "O2-".

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

        if match := re.search(r"([A-Z][a-z]*)([0-9.]*)([+\-]*)(.*)", species_string):
            sym = match[1]  # parse symbol

            # Parse oxidation state (optional)
            if not match[2] and not match[3]:
                oxi = None
            else:
                oxi = 1 if match[2] == "" else float(match[2])
                oxi = -oxi if match[3] == "-" else oxi

            # Parse properties (optional)
            properties = {}
            if match[4]:  # has Spin properties
                tokens = match[4].replace(",", "").split("=")
                properties = {tokens[0]: ast.literal_eval(tokens[1])}

            # But we need either an oxidation state or a property
            if oxi is None and not properties:
                raise ValueError("Invalid species string")

            return cls(sym, 0 if oxi is None else oxi, **properties)
        raise ValueError("Invalid species string")

    def to_pretty_string(self) -> str:
        """String without properties."""
        output = self.symbol
        if self.oxi_state is not None:
            abs_charge = formula_double_format(abs(self.oxi_state))
            if isinstance(abs_charge, float):
                abs_charge = f"{abs_charge:.2f}"  # type: ignore[assignment]
            output += f"{abs_charge}{'+' if self.oxi_state >= 0 else '-'}"
        return output

    def get_nmr_quadrupole_moment(self, isotope: str | None = None) -> float:
        """Get the nuclear electric quadrupole moment in units of e * millibarns.

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
            raise ValueError(f"No quadrupole moment for {isotope=}")
        return quad_mom.get(isotope, 0.0)

    def get_shannon_radius(
        self,
        cn: str,
        spin: Literal["", "Low Spin", "High Spin"] = "",
        radius_type: Literal["ionic", "crystal"] = "ionic",
    ) -> float:
        """Get the local environment specific ionic radius for species.

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
            Shannon radius for species in the specified environment.
        """
        radii = self._el.data["Shannon radii"]
        if self._oxi_state is None:
            raise ValueError("oxi_state is None.")
        radii = radii[str(int(self._oxi_state))][cn]
        if len(radii) == 1:
            key, data = next(iter(radii.items()))
            if key != spin:
                warnings.warn(
                    f"Specified {spin=} not consistent with database spin of {key}. "
                    "Only one spin data available, and that value is returned.",
                    stacklevel=2,
                )
        else:
            data = radii[spin]
        return data[f"{radius_type}_radius"]

    def get_crystal_field_spin(
        self,
        coordination: Literal["oct", "tet"] = "oct",
        spin_config: Literal["low", "high"] = "high",
    ) -> float:
        """Calculate the crystal field spin based on coordination and spin
        configuration. Only works for transition metal species.

        Args:
            coordination ("oct" | "tet"): Tetrahedron or octahedron crystal site coordination
            spin_config ("low" | "high"): Whether the species is in a high or low spin state

        Returns:
            float: Crystal field spin in Bohr magneton.

        Raises:
            AttributeError if species is not a valid transition metal or has
                an invalid oxidation state.
            ValueError if invalid coordination or spin_config.
        """
        if coordination not in ("oct", "tet") or spin_config not in ("high", "low"):
            raise ValueError("Invalid coordination or spin config")

        elec = self.element.full_electronic_structure
        if len(elec) < 4 or elec[-2][1] != "s" or elec[-1][1] != "d":
            raise AttributeError(f"Invalid element {self.symbol} for crystal field calculation")

        if self.oxi_state is None:
            raise ValueError("oxi_state is None.")
        n_electrons = elec[-1][2] + elec[-2][2] - self.oxi_state
        if n_electrons < 0 or n_electrons > 10:
            raise AttributeError(f"Invalid oxidation state {self.oxi_state} for element {self.symbol}")

        if spin_config == "high":
            return n_electrons if n_electrons <= 5 else 10 - n_electrons

        if spin_config == "low":
            if coordination == "oct":
                if n_electrons <= 3:
                    return n_electrons
                if n_electrons <= 6:
                    return 6 - n_electrons
                if n_electrons <= 8:
                    return n_electrons - 6
                return 10 - n_electrons

            if coordination == "tet":
                if n_electrons <= 2:
                    return n_electrons
                if n_electrons <= 4:
                    return 4 - n_electrons
                if n_electrons <= 7:
                    return n_electrons - 4
                return 10 - n_electrons
            return None
        return None

    def as_dict(self) -> dict:
        """JSON-able dictionary representation."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "element": self.symbol,
            "oxidation_state": self._oxi_state,
            "spin": self._spin,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            Species.
        """
        return cls(dct["element"], dct["oxidation_state"], spin=dct.get("spin"))


@functools.total_ordering
class DummySpecies(Species):
    """A special species for representing non-traditional elements or species. For
    example, representation of vacancies (charged or otherwise), or special
    sites, etc.

    Attributes:
        oxi_state (int): Oxidation state associated with Species.
        Z (int): DummySpecies is always assigned an atomic number equal to the hash
            number of the symbol. Obviously, it makes no sense whatsoever to use
            the atomic number of a Dummy species for anything scientific. The purpose
            of this is to ensure that for most use cases, a DummySpecies behaves no
            differently from an Element or Species.
        A (int): Just as for Z, to get a DummySpecies to behave like an Element,
            it needs atomic mass number A (arbitrarily set to twice Z).
        X (float): DummySpecies is always assigned a Pauling electronegativity of 0.
    """

    def __init__(
        self,
        symbol: str = "X",
        oxidation_state: float | None = 0,
        spin: float | None = None,
    ) -> None:
        """
        Args:
            symbol (str): An assigned symbol for the dummy specie. Strict
                rules are applied to the choice of the symbol. The dummy
                symbol cannot have any part of first two letters that will
                constitute an Element symbol. Otherwise, a composition may
                be parsed wrongly. e.g. "X" is fine, but "Vac" is not
                because Vac contains V, a valid Element.
            oxidation_state (float): Oxidation state for dummy specie. Defaults to 0.
                deprecated and retained purely for backward compatibility.
            spin: Spin associated with Species. Defaults to None.
        """
        # enforce title case to match other elements, reduces confusion
        # when multiple DummySpecies in a "formula" string
        symbol = symbol.title()

        for idx in range(1, min(2, len(symbol)) + 1):
            if Element.is_valid_symbol(symbol[:idx]):
                raise ValueError(f"{symbol} contains {symbol[:idx]}, which is a valid element symbol")

        # Set required attributes for DummySpecies to function like a Species in
        # most instances.
        self._symbol = symbol
        self._oxi_state = oxidation_state
        self._spin = spin

    def __getattr__(self, attr: str) -> None:
        raise AttributeError

    def __lt__(self, other) -> bool:
        """Set a default sort order for atomic species by Pauling electronegativity,
        followed by oxidation state.
        """
        if self.X != other.X:
            return self.X < other.X
        if self.symbol != other.symbol:
            # There are cases where the Pauling electronegativity are exactly equal.
            # We then sort by symbol.
            return self.symbol < other.symbol
        other_oxi = 0 if isinstance(other, Element) else other.oxi_state
        if self.oxi_state is not None and other_oxi is not None:
            return self.oxi_state < other_oxi

        raise RuntimeError("oxi_state for both species are None")

    def __repr__(self) -> str:
        return f"DummySpecies {self}"

    def __str__(self) -> str:
        output = self.symbol
        if self.oxi_state is not None:
            output += f"{formula_double_format(abs(self.oxi_state))}{'+' if self.oxi_state >= 0 else '-'}"
        if self._spin is not None:
            spin = self._spin
            output += f",{spin=}"
        return output

    def __deepcopy__(self, memo) -> Self:
        return type(self)(self.symbol, self._oxi_state)

    @property
    def Z(self) -> int:
        """
        Proton number of DummySpecies.

        DummySpecies is always assigned an atomic number equal to the hash of
        the symbol. This is necessary for the DummySpecies object to behave like
        an ElementBase object (which Species inherits from).
        """
        return hash(self.symbol)

    @property
    def A(self) -> int | None:
        """
        Atomic mass number of a DummySpecies.

        To behave like an ElementBase object (from which Species inherits),
        DummySpecies needs an atomic mass number. Consistent with the
        implementation for ElementBase, this can return an int or None (default).
        """
        return None

    @property
    def oxi_state(self) -> float | None:
        """Oxidation state associated with DummySpecies."""
        return self._oxi_state

    @property
    def X(self) -> float:
        """DummySpecies is always assigned a Pauling electronegativity of 0. The effect of
        this is that DummySpecies are always sorted in front of actual Species.
        """
        return 0.0

    @property
    def symbol(self) -> str:
        """Symbol for DummySpecies."""
        return self._symbol

    @classmethod
    def from_str(cls, species_string: str) -> Self:
        """Get a Dummy from a string representation.

        Args:
            species_string (str): A string representation of a dummy
                species, e.g. "X2+", "X3+".

        Returns:
            A DummySpecies object.

        Raises:
            ValueError if species_string cannot be interpreted.
        """
        if match := re.search(r"([A-ZAa-z]*)([0-9.]*)([+\-]*)(.*)", species_string):
            sym = match[1]
            if match[2] == match[3] == "":
                oxi = 0.0
            else:
                oxi = 1.0 if match[2] == "" else float(match[2])
                oxi = -oxi if match[3] == "-" else oxi
            properties = {}
            if match[4]:  # has Spin property
                tokens = match[4].split("=")
                properties = {tokens[0]: float(tokens[1])}
            return cls(sym, oxi, **properties)
        raise ValueError("Invalid DummySpecies String")

    def as_dict(self) -> dict[str, Any]:
        """MSONable dict representation."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "element": self.symbol,
            "oxidation_state": self._oxi_state,
            "spin": self._spin,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            DummySpecies
        """
        return cls(dct["element"], dct["oxidation_state"], spin=dct.get("spin"))


class Specie(Species):
    """This maps the historical grammatically inaccurate Specie to Species
    to maintain backwards compatibility.
    """


class DummySpecie(DummySpecies):
    """This maps the historical grammatically inaccurate DummySpecie to DummySpecies
    to maintain backwards compatibility.
    """


@overload
def get_el_sp(obj: int) -> Element:
    pass


@overload
def get_el_sp(obj: SpeciesLike) -> Element | Species | DummySpecies:
    pass


@functools.lru_cache(maxsize=1024)
def get_el_sp(obj: int | SpeciesLike) -> Element | Species | DummySpecies:
    """Utility function to get an Element, Species or DummySpecies from any input.

    If obj is an Element or a Species, it is returned as is.
    If obj is an int or a string representing an integer, the Element with the
    atomic number obj is returned.
    If obj is a string, Species parsing will be attempted (e.g. Mn2+). Failing that
    Element parsing will be attempted (e.g. Mn). Failing that DummyElement parsing
    will be attempted.

    Args:
        obj (SpeciesLike): An arbitrary object. Supported objects are actual Element/Species,
            integers (representing atomic numbers) or strings (element symbols or species strings).

    Raises:
        ValueError: if obj cannot be converted into an Element or Species.

    Returns:
        Element | Species | DummySpecies: with a bias for the maximum number
            of properties that can be determined.
    """
    # If obj is already an Element or Species, return as is
    if isinstance(obj, Element | Species | DummySpecies):
        return obj

    # If obj is an integer, return the Element with atomic number obj
    try:
        flt = float(obj)
        assert flt == int(flt)  # noqa: S101
        return Element.from_Z(int(flt))
    except (AssertionError, ValueError, TypeError, KeyError):
        pass

    # If obj is a string, attempt to parse it as a Species
    try:
        return Species.from_str(obj)  # type: ignore[arg-type]
    except (ValueError, TypeError, KeyError):
        pass
    # If Species parsing failed, try Element
    try:
        return Element(obj)  # type: ignore[arg-type]
    except (ValueError, TypeError, KeyError):
        pass

    # If Element parsing failed, try DummySpecies
    try:
        return DummySpecies.from_str(obj)  # type: ignore[arg-type]
    except Exception as exc:
        raise ValueError(f"Can't parse Element or Species from {obj!r}") from exc


@unique
class ElementType(Enum):
    """Enum for element types."""

    noble_gas = "noble_gas"  # He, Ne, Ar, Kr, Xe, Rn
    transition_metal = "transition_metal"  # Sc-Zn, Y-Cd, La-Hg, Ac-Cn
    post_transition_metal = "post_transition_metal"  # Al, Ga, In, Tl, Sn, Pb, Bi, Po
    rare_earth_metal = "rare_earth_metal"  # Ce-Lu, Th-Lr
    metal = "metal"
    metalloid = "metalloid"  # B, Si, Ge, As, Sb, Te, Po
    alkali = "alkali"  # Li, Na, K, Rb, Cs, Fr
    alkaline = "alkaline"  # Be, Mg, Ca, Sr, Ba, Ra
    halogen = "halogen"  # F, Cl, Br, I, At
    chalcogen = "chalcogen"  # O, S, Se, Te, Po
    lanthanoid = "lanthanoid"  # La-Lu
    actinoid = "actinoid"  # Ac-Lr
    radioactive = "radioactive"  # Tc, Pm, Po-Lr
    quadrupolar = "quadrupolar"
    s_block = "s-block"
    p_block = "p-block"
    d_block = "d-block"
    f_block = "f-block"
