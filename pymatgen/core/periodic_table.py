# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""Module contains classes presenting Element and Specie (Element + oxidation state) and PeriodicTable."""

import re
import json
import warnings
from io import open
from pathlib import Path
from enum import Enum
from typing import Optional, Callable
from itertools import product, \
    combinations
from collections import Counter

import numpy as np

from pymatgen.core.units import Mass, Length, FloatWithUnit, Unit, \
    SUPPORTED_UNIT_NAMES
from pymatgen.util.string import formula_double_format
from monty.json import MSONable

# Loads element data from json file
with open(str(Path(__file__).absolute().parent / "periodic_table.json"), "rt") as f:
    _pt_data = json.load(f)

_pt_row_sizes = (2, 8, 8, 18, 18, 32, 32)


class Element(Enum):
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

        .. attribute:: X

            Pauling electronegativity. Elements without an electronegativity
            number are assigned a value of zero by default.

        .. attribute:: number

            Alternative attribute for atomic number

        .. attribute:: max_oxidation_state

            Maximum oxidation state for element

        .. attribute:: min_oxidation_state

            Minimum oxidation state for element

        .. attribute:: oxidation_states

            Tuple of all known oxidation states

        .. attribute:: common_oxidation_states

            Tuple of all common oxidation states

        .. attribute:: full_electronic_structure

            Full electronic structure as tuple.
            E.g., The electronic structure for Fe is represented as:
            [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6),
            (3, "d", 6), (4, "s", 2)]

        .. attribute:: row

            Returns the periodic table row of the element.

        .. attribute:: group

            Returns the periodic table group of the element.

        .. attribute:: block

            Return the block character "s,p,d,f"

        .. attribute:: is_noble_gas

            True if element is noble gas.

        .. attribute:: is_transition_metal

            True if element is a transition metal.

        .. attribute:: is_post_transition_metal

            True if element is a post transition metal.

        .. attribute:: is_rare_earth_metal

            True if element is a rare earth metal.

        .. attribute:: is_metalloid

            True if element is a metalloid.

        .. attribute:: is_alkali

            True if element is an alkali metal.

        .. attribute:: is_alkaline

            True if element is an alkaline earth metal (group II).

        .. attribute:: is_halogen

            True if element is a halogen.

        .. attribute:: is_lanthanoid

            True if element is a lanthanoid.

        .. attribute:: is_actinoid

            True if element is a actinoid.

        .. attribute:: iupac_ordering

            Ordering according to Table VI of "Nomenclature of Inorganic Chemistry
            (IUPAC Recommendations 2005)". This ordering effectively follows the
            groups and rows of the periodic table, except the Lanthanides, Actanides
            and hydrogen.

        .. attribute:: long_name

           Long name for element. E.g., "Hydrogen".

        .. attribute:: atomic_mass

            Atomic mass for the element.

        .. attribute:: atomic_radius

            Atomic radius for the element. This is the empirical value. Data is
            obtained from
            http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page).

        .. attribute:: atomic_radius_calculated

            Calculated atomic radius for the element. This is the empirical value.
            Data is obtained from
            http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page).

        .. attribute:: van_der_waals_radius

            Van der Waals radius for the element. This is the empirical
            value. Data is obtained from
            http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page).

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

        .. attribute:: average_ionic_radius

            Average ionic radius for element in ang. The average is taken over all
            oxidation states of the element for which data is present.

        .. attribute:: average_cationic_radius

            Average cationic radius for element in ang. The average is taken over all
            positive oxidation states of the element for which data is present.

        .. attribute:: average_anionic_radius

            Average ionic radius for element in ang. The average is taken over all
            negative oxidation states of the element for which data is present.

        .. attribute:: ionic_radii

            All ionic radii of the element as a dict of
            {oxidation state: ionic radii}. Radii are given in ang.
        """
        self.symbol = "%s" % symbol
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
    def X(self):
        """
        :return: Electronegativity of element. Note that if an element does not
            have an electronegativity, a NaN float is returned.
        """
        if "X" in self._data:
            return self._data["X"]
        warnings.warn("No electronegativity for %s. Setting to NaN. "
                      "This has no physical meaning, and is mainly done to "
                      "avoid errors caused by the code expecting a float."
                      % self.symbol)
        return float("NaN")

    @property
    def atomic_radius(self):
        """
        Returns: The atomic radius of the element in Ångstroms.
        """
        return self._atomic_radius

    @property
    def atomic_mass(self):
        """
        Returns: The atomic mass of the element in amu.
        """
        return self._atomic_mass

    def __getattr__(self, item):
        if item in ["mendeleev_no", "electrical_resistivity",
                    "velocity_of_sound", "reflectivity",
                    "refractive_index", "poissons_ratio", "molar_volume",
                    "thermal_conductivity", "boiling_point", "melting_point",
                    "critical_temperature", "superconduction_temperature",
                    "liquid_range", "bulk_modulus", "youngs_modulus",
                    "brinell_hardness", "rigidity_modulus",
                    "mineral_hardness", "vickers_hardness",
                    "density_of_solid", "atomic_radius_calculated",
                    "van_der_waals_radius", "atomic_orbitals",
                    "coefficient_of_linear_thermal_expansion",
                    "ground_state_term_symbol", "valence"]:
            kstr = item.capitalize().replace("_", " ")
            val = self._data.get(kstr, None)
            if str(val).startswith("no data"):
                val = None
            elif isinstance(val, dict):
                pass
            else:
                try:
                    val = float(val)
                except ValueError:
                    nobracket = re.sub(r'\(.*\)', "", val)
                    toks = nobracket.replace("about", "").strip().split(" ", 1)
                    if len(toks) == 2:
                        try:
                            if "10<sup>" in toks[1]:
                                base_power = re.findall(r'([+-]?\d+)', toks[1])
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
                                unit = toks[1].replace("<sup>", "^").replace(
                                    "</sup>", "").replace("&Omega;",
                                                          "ohm")
                                units = Unit(unit)
                                if set(units.keys()).issubset(
                                        SUPPORTED_UNIT_NAMES):
                                    val = FloatWithUnit(toks[0], unit)
                        except ValueError:
                            # Ignore error. val will just remain a string.
                            pass
            return val
        raise AttributeError("Element has no attribute %s!" % item)

    @property
    def data(self) -> dict:
        """
        Returns dict of data for element.
        """
        return self._data.copy()

    @property
    def electronic_structure(self) -> str:
        """
        Electronic structure as string, with only valence electrons.
        E.g., The electronic structure for Fe is represented as '[Ar].3d6.4s2'
        """
        return re.sub("</*sup>", "", self._data["Electronic structure"])

    @property
    def average_ionic_radius(self):
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
    def average_cationic_radius(self):
        """
        Average cationic radius for element (with units). The average is
        taken over all positive oxidation states of the element for which
        data is present.
        """
        if "Ionic radii" in self._data:
            radii = [v for k, v in self._data["Ionic radii"].items()
                     if int(k) > 0]
            if radii:
                return FloatWithUnit(sum(radii) / len(radii), "ang")
        return FloatWithUnit(0.0, "ang")

    @property
    def average_anionic_radius(self):
        """
        Average anionic radius for element (with units). The average is
        taken over all negative oxidation states of the element for which
        data is present.
        """
        if "Ionic radii" in self._data:
            radii = [v for k, v in self._data["Ionic radii"].items()
                     if int(k) < 0]
            if radii:
                return FloatWithUnit(sum(radii) / len(radii), "ang")
        return FloatWithUnit(0.0, "ang")

    @property
    def ionic_radii(self):
        """
        All ionic radii of the element as a dict of
        {oxidation state: ionic radii}. Radii are given in ang.
        """
        if "Ionic radii" in self._data:
            return {int(k): FloatWithUnit(v, "ang") for k, v in self._data["Ionic radii"].items()}
        return {}

    @property
    def number(self):
        """Alternative attribute for atomic number"""
        return self.Z

    @property
    def max_oxidation_state(self):
        """Maximum oxidation state for element"""
        if "Oxidation states" in self._data:
            return max(self._data["Oxidation states"])
        return 0

    @property
    def min_oxidation_state(self):
        """Minimum oxidation state for element"""
        if "Oxidation states" in self._data:
            return min(self._data["Oxidation states"])
        return 0

    @property
    def oxidation_states(self):
        """Tuple of all known oxidation states"""
        return tuple(self._data.get("Oxidation states", list()))

    @property
    def common_oxidation_states(self):
        """Tuple of all known oxidation states"""
        return tuple(self._data.get("Common oxidation states", list()))

    @property
    def icsd_oxidation_states(self):
        """Tuple of all oxidation states with at least 10 instances in
        ICSD database AND at least 1% of entries for that element"""
        return tuple(self._data.get("ICSD oxidation states", list()))

    @property
    def metallic_radius(self):
        """
        Metallic radius of the element. Radius is given in ang.
        """
        return FloatWithUnit(self._data["Metallic radius"], "ang")

    @property
    def full_electronic_structure(self):
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
            data = Element(sym).full_electronic_structure + data[1:]
        return data

    @property
    def valence(self):
        """
        # From full electron config obtain valence subshell
        # angular moment (L) and number of valence e- (v_e)

        """
        # the number of valence of noble gas is 0
        if self.group == 18:
            return (np.nan, 0)

        L_symbols = 'SPDFGHIKLMNOQRTUVWXYZ'
        valence = []
        full_electron_config = self.full_electronic_structure
        for _, l_symbol, ne in full_electron_config[::-1]:
            l = L_symbols.lower().index(l_symbol)
            if ne < (2 * l + 1) * 2:
                valence.append((l, ne))
        if len(valence) > 1:
            raise ValueError("Ambiguous valence")

        return valence[0]

    @property
    def term_symbols(self):
        """
        All possible  Russell-Saunders term symbol of the Element
        eg. L = 1, n_e = 2 (s2)
        returns
           [['1D2'], ['3P0', '3P1', '3P2'], ['1S0']]

        """
        L_symbols = 'SPDFGHIKLMNOQRTUVWXYZ'

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
        TL = [sum([ml_ms[comb[e]][0] for e in range(v_e)])
              for comb in e_config_combs]
        TS = [sum([ml_ms[comb[e]][1] for e in range(v_e)])
              for comb in e_config_combs]
        comb_counter = Counter([r for r in zip(TL, TS)])

        term_symbols = []
        while sum(comb_counter.values()) > 0:
            # Start from the lowest freq combination,
            # which corresponds to largest abs(L) and smallest abs(S)
            L, S = min(comb_counter)

            J = list(np.arange(abs(L - S), abs(L) + abs(S) + 1))
            term_symbols.append([str(int(2 * (abs(S)) + 1)) + L_symbols[abs(L)] + str(j) for j in J])
            # Without J
            # term_symbols.append(str(int(2 * (abs(S)) + 1)) \
            #                     + L_symbols[abs(L)])

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
        L_symbols = 'SPDFGHIKLMNOQRTUVWXYZ'

        term_symbols = self.term_symbols
        term_symbol_flat = {term: {"multiplicity": int(term[0]),
                                   "L": L_symbols.index(term[1]),
                                   "J": float(term[2:])}
                            for term in sum(term_symbols, [])}

        multi = [int(item['multiplicity'])
                 for terms, item in term_symbol_flat.items()]
        max_multi_terms = {symbol: item
                           for symbol, item in term_symbol_flat.items()
                           if item['multiplicity'] == max(multi)}

        Ls = [item['L'] for terms, item in max_multi_terms.items()]
        max_L_terms = {symbol: item
                       for symbol, item in term_symbol_flat.items()
                       if item['L'] == max(Ls)}

        J_sorted_terms = sorted(max_L_terms.items(),
                                key=lambda k: k[1]['J'])
        L, v_e = self.valence
        if v_e <= (2 * L + 1):
            return J_sorted_terms[0][0]
        return J_sorted_terms[-1][0]

    def __eq__(self, other):
        return isinstance(other, Element) and self.Z == other.Z

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self.Z

    def __repr__(self):
        return "Element " + self.symbol

    def __str__(self):
        return self.symbol

    def __lt__(self, other):
        """
        Sets a default sort order for atomic species by electronegativity. Very
        useful for getting correct formulas.  For example, FeO4PLi is
        automatically sorted into LiFePO4.
        """
        x1 = float("inf") if self.X != self.X else self.X
        x2 = float("inf") if other.X != other.X else other.X
        if x1 != x2:
            return x1 < x2

        # There are cases where the electronegativity are exactly equal.
        # We then sort by symbol.
        return self.symbol < other.symbol

    @staticmethod
    def from_Z(z: int):
        """
        Get an element from an atomic number.

        Args:
            z (int): Atomic number

        Returns:
            Element with atomic number z.
        """
        for sym, data in _pt_data.items():
            if data["Atomic no"] == z:
                return Element(sym)
        raise ValueError("No element with this atomic number %s" % z)

    @staticmethod
    def from_row_and_group(row: int, group: int):
        """
        Returns an element from a row and group number.

        Args:
            row (int): Row number
            group (int): Group number

        .. note::
            The 18 group number system is used, i.e., Noble gases are group 18.
        """
        for sym in _pt_data.keys():
            el = Element(sym)
            if el.row == row and el.group == group:
                return el
        raise ValueError("No element with this row and group!")

    @staticmethod
    def is_valid_symbol(symbol: str):
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
    def row(self):
        """
        Returns the periodic table row of the element.
        """
        z = self.Z
        total = 0
        if 57 <= z <= 71:
            return 8
        if 89 <= z <= 103:
            return 9
        for i, size in enumerate(_pt_row_sizes):
            total += size
            if total >= z:
                return i + 1
        return 8

    @property
    def group(self):
        """
        Returns the periodic table group of the element.
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

        if (z - 54) % 32 == 0:
            return 18
        if (z - 54) % 32 >= 18:
            return (z - 54) % 32 - 14
        return (z - 54) % 32

    @property
    def block(self):
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
    def is_noble_gas(self):
        """
        True if element is noble gas.
        """
        return self.Z in (2, 10, 18, 36, 54, 86, 118)

    @property
    def is_transition_metal(self):
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
    def is_post_transition_metal(self):
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
        :return: True if is a metal.
        """
        return (self.is_alkali or self.is_alkaline or
                self.is_post_transition_metal or self.is_transition_metal or
                self.is_lanthanoid or self.is_actinoid)

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
    def is_halogen(self):
        """
        True if element is a halogen.
        """
        return self.Z in (9, 17, 35, 53, 85)

    @property
    def is_chalcogen(self):
        """
        True if element is a chalcogen.
        """
        return self.Z in (8, 16, 34, 52, 84)

    @property
    def is_lanthanoid(self):
        """
        True if element is a lanthanoid.
        """
        return 56 < self.Z < 72

    @property
    def is_actinoid(self):
        """
        True if element is a actinoid.
        """
        return 88 < self.Z < 104

    @property
    def is_quadrupolar(self):
        """
        Checks if this element can be quadrupolar
        """
        return len(self.data.get("NMR Quadrupole Moment", {})) > 0

    @property
    def nmr_quadrupole_moment(self):
        """
        Get a dictionary the nuclear electric quadrupole moment in units of
        e*millibarns for various isotopes
        """
        return {k: FloatWithUnit(v, "mbarn")
                for k, v in self.data.get("NMR Quadrupole Moment", {}).items()}

    @property
    def iupac_ordering(self):
        """
        Ordering according to Table VI of "Nomenclature of Inorganic Chemistry
        (IUPAC Recommendations 2005)". This ordering effectively follows the
        groups and rows of the periodic table, except the Lanthanides, Actanides
        and hydrogen.
        """
        return self._data["IUPAC ordering"]

    def __deepcopy__(self, memo):
        return Element(self.symbol)

    @staticmethod
    def from_dict(d):
        """
        Makes Element obey the general json interface used in pymatgen for
        easier serialization.
        """
        return Element(d["element"])

    def as_dict(self):
        """
        Makes Element obey the general json interface used in pymatgen for
        easier serialization.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "element": self.symbol}

    @staticmethod
    def print_periodic_table(filter_function: Optional[Callable] = None):
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
                    rowstr.append("{:3s}".format(el.symbol))
                else:
                    rowstr.append("   ")
            print(" ".join(rowstr))


class Specie(MSONable):
    """
    An extension of Element with an oxidation state and other optional
    properties. Properties associated with Specie should be "idealized"
    values, not calculated values. For example, high-spin Fe2+ may be
    assigned an idealized spin of +5, but an actual Fe2+ site may be
    calculated to have a magmom of +4.5. Calculated properties should be
    assigned to Site objects, and not Specie.
    """

    supported_properties = ("spin",)

    def __init__(self, symbol: str,
                 oxidation_state: float = 0.0,
                 properties: dict = None):
        """
        Initializes a Specie.

        Args:
            symbol (str): Element symbol, e.g., Fe
            oxidation_state (float): Oxidation state of element, e.g., 2 or -2
            properties: Properties associated with the Specie, e.g.,
                {"spin": 5}. Defaults to None. Properties must be one of the
                Specie supported_properties.

        .. attribute:: oxi_state

            Oxidation state associated with Specie

        .. attribute:: ionic_radius

            Ionic radius of Specie (with specific oxidation state).

        .. versionchanged:: 2.6.7

            Properties are now checked when comparing two Species for equality.
        """
        self._el = Element(symbol)
        self._oxi_state = oxidation_state
        self._properties = properties if properties else {}
        for k in self._properties.keys():
            if k not in Specie.supported_properties:
                raise ValueError("{} is not a supported property".format(k))

    def __getattr__(self, a):
        # overriding getattr doesn't play nice with pickle, so we
        # can't use self._properties
        p = object.__getattribute__(self, '_properties')
        if a in p:
            return p[a]
        return getattr(self._el, a)

    def __eq__(self, other):
        """
        Specie is equal to other only if element and oxidation states are
        exactly the same.
        """
        return (isinstance(other, Specie) and self.symbol == other.symbol
                and self.oxi_state == other.oxi_state
                and self._properties == other._properties)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        """
        Equal Specie should have the same str representation, hence
        should hash equally. Unequal Specie will have differnt str
        representations.
        """
        return self.__str__().__hash__()

    def __lt__(self, other):
        """
        Sets a default sort order for atomic species by electronegativity,
        followed by oxidation state, followed by spin.
        """
        x1 = float("inf") if self.X != self.X else self.X
        x2 = float("inf") if other.X != other.X else other.X
        if x1 != x2:
            return x1 < x2
        if self.symbol != other.symbol:
            # There are cases where the electronegativity are exactly equal.
            # We then sort by symbol.
            return self.symbol < other.symbol
        if self.oxi_state:
            other_oxi = 0 if (isinstance(other, Element)
                              or other.oxi_state is None) else other.oxi_state
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
    def ionic_radius(self):
        """
        Ionic radius of specie. Returns None if data is not present.
        """

        if self._oxi_state in self.ionic_radii:
            return self.ionic_radii[self._oxi_state]
        d = self._el.data
        oxstr = str(int(self._oxi_state))
        if oxstr in d.get("Ionic radii hs", {}):
            warnings.warn("No default ionic radius for %s. Using hs data." %
                          self)
            return d["Ionic radii hs"][oxstr]
        if oxstr in d.get("Ionic radii ls", {}):
            warnings.warn("No default ionic radius for %s. Using ls data." %
                          self)
            return d["Ionic radii ls"][oxstr]
        warnings.warn("No ionic radius for {}!".format(self))
        return None

    @property
    def oxi_state(self):
        """
        Oxidation state of Specie.
        """
        return self._oxi_state

    @staticmethod
    def from_string(species_string: str):
        """
        Returns a Specie from a string representation.

        Args:
            species_string (str): A typical string representation of a
                species, e.g., "Mn2+", "Fe3+", "O2-".

        Returns:
            A Specie object.

        Raises:
            ValueError if species_string cannot be intepreted.
        """
        m = re.search(r"([A-Z][a-z]*)([0-9.]*)([+\-])(.*)", species_string)
        if m:
            sym = m.group(1)
            oxi = 1 if m.group(2) == "" else float(m.group(2))
            oxi = -oxi if m.group(3) == "-" else oxi
            properties = None
            if m.group(4):
                toks = m.group(4).replace(",", "").split("=")
                properties = {toks[0]: float(toks[1])}
            return Specie(sym, oxi, properties)
        raise ValueError("Invalid Species String")

    def __repr__(self):
        return "Specie " + self.__str__()

    def __str__(self):
        output = self.symbol
        if self.oxi_state is not None:
            if self.oxi_state >= 0:
                output += formula_double_format(self.oxi_state) + "+"
            else:
                output += formula_double_format(-self.oxi_state) + "-"
        for p, v in self._properties.items():
            output += ",%s=%s" % (p, v)
        return output

    def get_nmr_quadrupole_moment(self, isotope=None):
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
            isotopes = list(quad_mom.keys())
            isotopes.sort(key=lambda x: int(x.split("-")[1]), reverse=False)
            return quad_mom.get(isotopes[0], 0.0)

        if isotope not in quad_mom:
            raise ValueError("No quadrupole moment for isotope {}".format(
                isotope))
        return quad_mom.get(isotope, 0.0)

    def get_shannon_radius(self, cn: str, spin: str = "",
                           radius_type: str = "ionic"):
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
        if len(radii) == 1:  # type: ignore
            k, data = list(radii.items())[0]  # type: ignore
            if k != spin:
                warnings.warn(
                    "Specified spin state of %s not consistent with database "
                    "spin of %s. Only one spin data available, and "
                    "that value is returned." % (spin, k)
                )
        else:
            data = radii[spin]
        return data["%s_radius" % radius_type]

    def get_crystal_field_spin(self, coordination: str = "oct",
                               spin_config: str = "high"):
        """
        Calculate the crystal field spin based on coordination and spin
        configuration. Only works for transition metal species.

        Args:
            coordination (str): Only oct and tet are supported at the moment.
            spin_config (str): Supported keywords are "high" or "low".

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
            raise AttributeError(
                "Invalid element {} for crystal field calculation.".format(self.symbol))
        nelectrons = elec[-1][2] + elec[-2][2] - self.oxi_state
        if nelectrons < 0 or nelectrons > 10:
            raise AttributeError(
                "Invalid oxidation state {} for element {}".format(self.oxi_state, self.symbol))
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
        return Specie(self.symbol, self.oxi_state, self._properties)

    def as_dict(self):
        """
        :return: Json-able dictionary representation.
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "element": self.symbol,
             "oxidation_state": self._oxi_state}
        if self._properties:
            d["properties"] = self._properties
        return d

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation.
        :return: Specie.
        """
        return cls(d["element"], d["oxidation_state"],
                   d.get("properties", None))


class DummySpecie(Specie):
    """
    A special specie for representing non-traditional elements or species. For
    example, representation of vacancies (charged or otherwise), or special
    sites, etc.

    .. attribute:: oxi_state

        Oxidation state associated with Specie.

    .. attribute:: Z

        DummySpecie is always assigned an atomic number equal to the hash
        number of the symbol. Obviously, it makes no sense whatsoever to use
        the atomic number of a Dummy specie for anything scientific. The purpose
        of this is to ensure that for most use cases, a DummySpecie behaves no
        differently from an Element or Specie.

    .. attribute:: X

        DummySpecie is always assigned an electronegativity of 0.
    """

    def __init__(self,
                 symbol: str = "X",
                 oxidation_state: float = 0,
                 properties: dict = None):
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
        # when multiple DummySpecie in a "formula" string
        symbol = symbol.title()

        for i in range(1, min(2, len(symbol)) + 1):
            if Element.is_valid_symbol(symbol[:i]):
                raise ValueError("{} contains {}, which is a valid element "
                                 "symbol.".format(symbol, symbol[:i]))

        # Set required attributes for DummySpecie to function like a Specie in
        # most instances.
        self._symbol = symbol
        self._oxi_state = oxidation_state
        self._properties = properties if properties else {}
        for k in self._properties.keys():
            if k not in Specie.supported_properties:
                raise ValueError("{} is not a supported property".format(k))

    def __getattr__(self, a):
        # overriding getattr doens't play nice with pickle, so we
        # can't use self._properties
        p = object.__getattribute__(self, '_properties')
        if a in p:
            return p[a]
        raise AttributeError(a)

    def __hash__(self):
        return self.symbol.__hash__()

    def __eq__(self, other):
        """
        Specie is equal to other only if element and oxidation states are
        exactly the same.
        """
        if not isinstance(other, DummySpecie):
            return False
        return (isinstance(other, Specie) and
                self.symbol == other.symbol and
                self.oxi_state == other.oxi_state and
                self._properties == other._properties)

    def __ne__(self, other):
        return not self.__eq__(other)

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
        DummySpecie is always assigned an atomic number equal to the hash of
        the symbol. The expectation is that someone would be an actual dummy
        to use atomic numbers for a Dummy specie.
        """
        return self.symbol.__hash__()

    @property
    def oxi_state(self) -> float:
        """
        Oxidation state associated with DummySpecie
        """
        return self._oxi_state

    @property
    def X(self) -> float:
        """
        DummySpecie is always assigned an electronegativity of 0. The effect of
        this is that DummySpecie are always sorted in front of actual Specie.
        """
        return 0.0

    @property
    def symbol(self) -> str:
        """
        :return: Symbol for DummySpecie.
        """
        return self._symbol

    def __deepcopy__(self, memo):
        return DummySpecie(self.symbol, self._oxi_state)

    @staticmethod
    def from_string(species_string: str):
        """
        Returns a Dummy from a string representation.

        Args:
            species_string (str): A string representation of a dummy
                species, e.g., "X2+", "X3+".

        Returns:
            A DummySpecie object.

        Raises:
            ValueError if species_string cannot be intepreted.
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
            return DummySpecie(sym, oxi, properties)
        raise ValueError("Invalid DummySpecies String")

    def as_dict(self):
        """
        :return: MSONAble dict representation.
        """
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "element": self.symbol,
             "oxidation_state": self._oxi_state}
        if self._properties:
            d["properties"] = self._properties
        return d

    @classmethod
    def from_dict(cls, d):
        """
        :param d: Dict representation
        :return: DummySpecie
        """
        return cls(d["element"], d["oxidation_state"],
                   d.get("properties", None))

    def __repr__(self):
        return "DummySpecie " + self.__str__()

    def __str__(self):
        output = self.symbol
        if self.oxi_state is not None:
            if self.oxi_state >= 0:
                output += formula_double_format(self.oxi_state) + "+"
            else:
                output += formula_double_format(-self.oxi_state) + "-"
        for p, v in self._properties.items():
            output += ",%s=%s" % (p, v)
        return output


def get_el_sp(obj):
    """
    Utility method to get an Element or Specie from an input obj.
    If obj is in itself an element or a specie, it is returned automatically.
    If obj is an int or a string representing an integer, the Element
    with the atomic number obj is returned.
    If obj is a string, Specie parsing will be attempted (e.g., Mn2+), failing
    which Element parsing will be attempted (e.g., Mn), failing which
    DummyElement parsing will be attempted.

    Args:
        obj (Element/Specie/str/int): An arbitrary object.  Supported objects
            are actual Element/Specie objects, integers (representing atomic
            numbers) or strings (element symbols or species strings).

    Returns:
        Specie or Element, with a bias for the maximum number of properties
        that can be determined.

    Raises:
        ValueError if obj cannot be converted into an Element or Specie.
    """
    if isinstance(obj, (Element, Specie, DummySpecie)):
        return obj

    if isinstance(obj, (list, tuple)):
        return [get_el_sp(o) for o in obj]

    try:
        c = float(obj)
        i = int(c)
        i = i if i == c else None
    except (ValueError, TypeError):
        i = None

    if i is not None:
        return Element.from_Z(i)

    try:
        return Specie.from_string(obj)
    except (ValueError, KeyError):
        try:
            return Element(obj)
        except (ValueError, KeyError):
            try:
                return DummySpecie.from_string(obj)
            except Exception:
                raise ValueError("Can't parse Element or String from type"
                                 " %s: %s." % (type(obj), obj))
