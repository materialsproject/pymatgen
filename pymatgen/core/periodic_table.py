#!/usr/bin/env python

"""
Module contains classes presenting Element and Specie (Element + oxidation
state) and PeriodicTable.
"""

__author__ = "Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import os
import re
import json

from pymatgen.util.decorators import singleton, cached_class
from pymatgen.util.string_utils import formula_double_format
from pymatgen.serializers.json_coders import MSONable
from functools import total_ordering

#Loads element data from json file
with open(os.path.join(os.path.dirname(__file__), "periodic_table.json")) as f:
    _pt_data = json.load(f)

_pt_row_sizes = (2, 8, 8, 18, 18, 32, 32)

# List with the correspondence Z --> Symbol
# We use a list instead of a mapping so that we can select slices easily.
_z2symbol = 119 * [None]
for (symbol, data) in _pt_data.items():
    _z2symbol[data["Atomic no"]] = symbol


def symbol_from_Z(z):
    """
    Return the symbol of the element from the atomic number.

    Args:
        z:
            Atomic number or slice object
    """
    return _z2symbol[z]


@cached_class
@total_ordering
class Element(object):
    """
    Basic immutable element object with all relevant properties.
    Only one instance of Element for each symbol is stored after creation,
    ensuring that a particular element behaves like a singleton. For all
    attributes, missing data (i.e., data for which is not available) is
    represented by a None unless otherwise stated.

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

    .. attribute:: name

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

        Mendeleev number

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

        Electronic structure. Simplified form with HTML formatting.
        E.g., The electronic structure for Fe is represented as
        [Ar].3d<sup>6</sup>.4s<sup>2</sup>

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

        Average ionic radius for element in pm. The average is taken over all
        oxidation states of the element for which data is present.

    .. attribute:: ionic_radii

        All ionic radii of the element as a dict of
        {oxidation state: ionic radii}. Radii are given in pm.
    """

    def __init__(self, symbol):
        """
        Create immutable element from a symbol.

        Args:
            symbol:
                Element symbol, e.g., "H", "Fe"
        """
        self._data = _pt_data[symbol]

        #Store key variables for quick access
        self._z = self._data["Atomic no"]
        self._symbol = symbol
        self._x = self._data.get("X", 0)

    def __getattr__(self, a):
        if a.lower() in ["name", "atomic_mass", "atomic_radius",
                         "mendeleev_no", "electrical_resistivity",
                         "velocity_of_sound", "reflectivity",
                         "refractive_index", "poissons_ratio", "molar_volume",
                         "electronic_structure", "thermal_conductivity",
                         "boiling_point", "melting_point",
                         "critical_temperature", "superconduction_temperature",
                         "liquid_range", "bulk_modulus", "youngs_modulus",
                         "brinell_hardness", "rigidity_modulus",
                         "mineral_hardness", "vickers_hardness",
                         "density_of_solid", "atomic_radius_calculated",
                         "van_der_waals_radius", "ionic_radii",
                         "coefficient_of_linear_thermal_expansion"]:
            kstr = a.capitalize().replace("_", " ")
            if kstr not in self._data:
                return None
            val = self._data[kstr]
            if str(val).startswith("no data"):
                return None
            else:
                return val
        raise AttributeError(a)

    def __getnewargs__(self):
        #function used by pickle to recreate object
        return self._symbol,

    @property
    def data(self):
        """
        Returns dict of data for element.
        """
        return self._data.copy()

    @property
    def average_ionic_radius(self):
        """
        Average ionic radius for element in pm. The average is taken over all
        oxidation states of the element for which data is present.
        """
        if "Ionic radii" in self._data:
            radii = self._data["Ionic radii"]
            return sum(radii.values()) / len(radii)
        else:
            return 0

    @property
    def ionic_radii(self):
        """
        All ionic radii of the element as a dict of
        {oxidation state: ionic radii}. Radii are given in pm.
        """
        if "Ionic radii" in self._data:
            return {int(k): v for k, v in self._data["Ionic radii"].items()}
        else:
            return {}

    @property
    def Z(self):
        """Atomic number"""
        return self._z

    @property
    def symbol(self):
        """Element symbol"""
        return self._symbol

    @property
    def X(self):
        """Electronegativity"""
        return self._x

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
    def full_electronic_structure(self):
        """
        Full electronic structure as tuple.
        E.g., The electronic structure for Fe is represented as:
        [(1, "s", 2), (2, "s", 2), (2, "p", 6), (3, "s", 2), (3, "p", 6),
        (3, "d", 6), (4, "s", 2)]
        """
        estr = self._data["Electronic structure"]

        def parse_orbital(orbstr):
            m = re.match("(\d+)([spdfg]+)<sup>(\d+)</sup>", orbstr)
            if m:
                return int(m.group(1)), m.group(2), int(m.group(3))
            return orbstr

        data = [parse_orbital(s) for s in estr.split(".")]
        if data[0][0] == "[":
            sym = data[0].replace("[", "").replace("]", "")
            data = Element(sym).full_electronic_structure + data[1:]
        return data

    def __eq__(self, other):
        if not isinstance(other, Element):
            return False
        return self._z == other._z

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._z

    def __repr__(self):
        return "Element " + self._symbol

    def __str__(self):
        return self._symbol

    def __lt__(self, other):
        """
        Sets a default sort order for atomic species by electronegativity. Very
        useful for getting correct formulas.  For example, FeO4PLi is
        automatically sorted into LiFePO4.
        """
        return self._x < other._x

    @staticmethod
    def from_Z(z):
        """
        Get an element from an atomic number.

        Args:
            z:
                Atomic number

        Returns:
            Element with atomic number z.
        """
        for sym, data in _pt_data.items():
            if data["Atomic no"] == z:
                return Element(sym)
        raise ValueError("No element with this atomic number")

    @staticmethod
    def from_row_and_group(row, group):
        """
        Returns an element from a row and group number.

        .. note::
            The 18 group number system is used, i.e., Noble gases are group 18.
        """
        for sym in _pt_data.keys():
            el = Element(sym)
            if el.row == row and el.group == group:
                return el
        return None

    @staticmethod
    def is_valid_symbol(symbol):
        """
        Returns true if symbol is a valid element symbol.

        Args:
            symbol:
                Element symbol

        Returns:
            True if symbol is a valid element (e.g., "H"). False otherwise
            (e.g., "Zebra").
        """
        return symbol in _pt_data

    @property
    def row(self):
        """
        Returns the periodic table row of the element.
        """
        Z = self._z
        total = 0
        if 57 <= Z <= 70:
            return 8
        elif 89 <= Z <= 102:
            return 9

        for i in range(len(_pt_row_sizes)):
            total += _pt_row_sizes[i]
            if total >= Z:
                return i + 1
        return 8

    @property
    def group(self):
        """
        Returns the periodic table group of the element.
        """
        Z = self._z
        if Z == 1:
            return 1
        if Z == 2:
            return 18
        if 3 <= Z <= 18:
            if (Z - 2) % 8 == 0:
                return 18
            elif (Z - 2) % 8 <= 2:
                return (Z - 2) % 8
            else:
                return 10 + (Z - 2) % 8

        if 19 <= Z <= 54:
            if (Z - 18) % 18 == 0:
                return 18
            else:
                return (Z - 18) % 18

        if (Z - 54) % 32 == 0:
            return 18
        elif (Z - 54) % 32 >= 17:
            return (Z - 54) % 32 - 14
        else:
            return (Z - 54) % 32

    @property
    def block(self):
        """
        Return the block character "s,p,d,f"
        """
        block = ""
        if self.group in [1, 2]:
            block = "s"
        elif self.group in range(13, 19):
            block = "p"
        elif self.is_actinoid or self.is_lanthanoid:
            block = "f"
        elif self.group in range(3, 13):
            block = "d"
        else:
            print("unable to determine block")
        return block

    @property
    def is_noble_gas(self):
        """
        True if element is noble gas.
        """
        return self._z in (2, 10, 18, 36, 54, 86, 118)

    @property
    def is_transition_metal(self):
        """
        True if element is a transition metal.
        """
        ns = list(range(21, 31))
        ns.extend(range(39, 49))
        ns.append(57)
        ns.extend(range(72, 81))
        ns.append(89)
        ns.extend(range(104, 113))
        return self._z in ns

    @property
    def is_rare_earth_metal(self):
        """
        True if element is a rare earth metal.
        """
        return self.is_lanthanoid or self.is_actinoid

    @property
    def is_metalloid(self):
        """
        True if element is a metalloid.
        """
        return self._symbol in ("B", "Si", "Ge", "As", "Sb", "Te", "Po")

    @property
    def is_alkali(self):
        """
        True if element is an alkali metal.
        """
        return self._z in (3, 11, 19, 37, 55, 87)

    @property
    def is_alkaline(self):
        """
        True if element is an alkaline earth metal (group II).
        """
        return self._z in (4, 12, 20, 38, 56, 88)

    @property
    def is_halogen(self):
        """
        True if element is a halogen.
        """
        return self._z in (9, 17, 35, 53, 85)

    @property
    def is_chalcogen(self):
        """
        True if element is a chalcogen.
        """
        return self._z in (8, 18, 34, 52, 84)

    @property
    def is_lanthanoid(self):
        """
        True if element is a lanthanoid.
        """
        return 56 < self._z < 72

    @property
    def is_actinoid(self):
        """
        True if element is a actinoid.
        """
        return 88 < self._z < 104

    def __deepcopy__(self, memo):
        return Element(self.symbol)

    @staticmethod
    def from_dict(d):
        """
        Makes Element obey the general json interface used in pymatgen for
        easier serialization.
        """
        return Element(d["element"])

    @property
    def to_dict(self):
        """
        Makes Element obey the general json interface used in pymatgen for
        easier serialization.
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "element": self._symbol}


@cached_class
@total_ordering
class Specie(MSONable):
    """
    An extension of Element with an oxidation state and other optional
    properties. Properties associated with Specie should be "idealized"
    values, not calculated values. For example, high-spin Fe2+ may be
    assigned an idealized spin of +5, but an actual Fe2+ site may be
    calculated to have a magmom of +4.5. Calculated properties should be
    assigned to Site objects, and not Specie.

    .. attribute:: oxi_state

        Oxidation state associated with Specie

    .. attribute:: ionic_radius

        Ionic radius of Specie (with specific oxidation state).

    .. versionchanged:: 2.6.7

        Properties are now checked when comparing two Species for equality.
    """

    supported_properties = ("spin",)

    def __init__(self, symbol, oxidation_state, properties=None):
        """
        Args:
            symbol:
                Element symbol, e.g., Fe
            oxidation_state:
                Oxidation state of element, e.g., 2 or -2
            properties:
                Properties associated with the Specie, e.g.,
                {"spin":5}. Defaults to None. Properties must be one of the
                Specie supported_properties.
        """
        self._el = Element(symbol)
        self._oxi_state = oxidation_state
        self._properties = properties if properties else {}
        for k in self._properties.keys():
            if k not in Specie.supported_properties:
                raise ValueError("{} is not a supported property".format(k))

    def __getattr__(self, a):
        #overriding getattr doens't play nice with pickle, so we
        #can't use self._properties
        p = object.__getattribute__(self, '_properties')
        if a in p:
            return p[a]
        try:
            return getattr(self._el, a)
        except:
            raise AttributeError(a)

    def __eq__(self, other):
        """
        Specie is equal to other only if element and oxidation states are
        exactly the same.
        """
        if not isinstance(other, Specie):
            return False
        return self.symbol == other.symbol \
            and self._oxi_state == other._oxi_state \
            and self._properties == other._properties

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        """
        Given that all oxidation states are below 100 in absolute value, this
        should effectively ensure that no two unequal Specie have the same
        hash.
        """
        return self.Z * 100 + self._oxi_state

    def __lt__(self, other):
        """
        Sets a default sort order for atomic species by electronegativity,
        followed by oxidation state.
        """
        other_oxi = 0 if isinstance(other, Element) else other.oxi_state
        return (self.X - other.X) * 100 + (self.oxi_state - other_oxi)

    @property
    def ionic_radius(self):
        """
        Ionic radius of specie. Returns None if data is not present.
        """
        return self.ionic_radii.get(self._oxi_state, None)

    @property
    def oxi_state(self):
        """
        Oxidation state of Specie.
        """
        return self._oxi_state

    @staticmethod
    def from_string(species_string):
        """
        Returns a Specie from a string representation.

        Args:
            species_string:
                A typical string representation of a species, e.g., "Mn2+",
                "Fe3+", "O2-".

        Returns:
            A Specie object.

        Raises:
            ValueError if species_string cannot be intepreted.
        """
        m = re.search("([A-Z][a-z]*)([0-9\.]*)([\+\-])", species_string)
        if m:
            num = 1 if m.group(2) == "" else float(m.group(2))
            return Specie(m.group(1), num if m.group(3) == "+" else -num)
        else:
            raise ValueError("Invalid Species String")

    def __repr__(self):
        return "Specie " + self.__str__()

    def __str__(self):
        output = self.symbol
        if self._oxi_state >= 0:
            output += formula_double_format(self._oxi_state) + "+"
        else:
            output += formula_double_format(-self._oxi_state) + "-"
        return output

    def get_crystal_field_spin(self, coordination="oct", spin_config="high"):
        """
        Calculate the crystal field spin based on coordination and spin
        configuration. Only works for transition metal species.

        Args:
            coordination:
                Only oct and tet are supported at the moment.
            spin_config:
                Supported keywords are "high" or "low".

        Returns:
            Crystal field spin in Bohr magneton.

        Raises:
            AttributeError if species is not a valid transition metal or has
            an invalid oxidation state.
            ValueError if invalid coordination or spin_config.
        """
        if coordination not in ("oct", "tet") or \
                spin_config not in ("high", "low"):
            raise ValueError("Invalid coordination or spin config.")
        elec = self.full_electronic_structure
        if len(elec) < 4 or elec[-1][1] != "s" or elec[-2][1] != "d":
            raise AttributeError(
                "Invalid element {} for crystal field calculation.".format(
                    self.symbol))
        nelectrons = elec[-1][2] + elec[-2][2] - self.oxi_state
        if nelectrons < 0:
            raise AttributeError(
                "Invalid oxidation state {} for element {}"
                .format(self.oxi_state, self.symbol))
        if spin_config == "high":
            return nelectrons if nelectrons <= 5 else 10 - nelectrons
        elif spin_config == "low":
            if coordination == "oct":
                if nelectrons <= 3:
                    return nelectrons
                elif nelectrons <= 6:
                    return 6 - nelectrons
                elif nelectrons <= 8:
                    return nelectrons - 6
                else:
                    return 10 - nelectrons
            elif coordination == "tet":
                if nelectrons <= 2:
                    return nelectrons
                elif nelectrons <= 4:
                    return 4 - nelectrons
                elif nelectrons <= 7:
                    return nelectrons - 4
                else:
                    return 10 - nelectrons

    def __deepcopy__(self, memo):
        return Specie(self.symbol, self.oxi_state, self._properties)

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "element": self.symbol,
                "oxidation_state": self._oxi_state,
                "properties": self._properties}

    @staticmethod
    def from_dict(d):
        return Specie(d["element"], d["oxidation_state"],
                      d.get("properties", None))


@cached_class
@total_ordering
class DummySpecie(MSONable):
    """
    A special specie for representing non-traditional elements or species. For
    example, representation of vacancies (charged or otherwise), or special
    sites, etc.

    .. attribute:: symbol

        Symbol for the DummySpecie.

    .. attribute:: oxi_state

        Oxidation state associated with Specie.

    .. attribute:: Z

        DummySpecie is always assigned an atomic number of 0.

    .. attribute:: X

        DummySpecie is always assigned an electronegativity of 0.
    """

    def __init__(self, symbol="X", oxidation_state=0, properties=None):
        """
        Args:
            symbol:
                An assigned symbol for the dummy specie. Strict rules are
                applied to the choice of the symbol. The dummy symbol cannot
                have any part of first two letters that will constitute an
                Element symbol. Otherwise, a composition may be parsed wrongly.
                E.g., "X" is fine, but "Vac" is not because Vac contains V, a
                valid Element.
            oxidation_state:
                Oxidation state for dummy specie. Defaults to zero.
        """
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
        #overriding getattr doens't play nice with pickle, so we
        #can't use self._properties
        p = object.__getattribute__(self, '_properties')
        if a in p:
            return p[a]
        try:
            return getattr(self._el, a)
        except:
            raise AttributeError(a)

    def __eq__(self, other):
        """
        Specie is equal to other only if element and oxidation states are
        exactly the same.
        """
        if not isinstance(other, DummySpecie):
            return False
        return self.symbol == other.symbol \
            and self._oxi_state == other._oxi_state

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        """
        Sets a default sort order for atomic species by electronegativity,
        followed by oxidation state.
        """
        other_oxi = 0 if isinstance(other, Element) else other.oxi_state
        return (self.X - other.X) * 100 + (self.oxi_state - other_oxi)

    @property
    def Z(self):
        """
        DummySpecie is always assigned an atomic number of 0.
        """
        return 0

    @property
    def oxi_state(self):
        """
        Oxidation state associated with DummySpecie
        """
        return self._oxi_state

    @property
    def X(self):
        """
        DummySpecie is always assigned an electronegativity of 0.
        """
        return 0

    @property
    def symbol(self):
        return self._symbol

    def __deepcopy__(self, memo):
        return DummySpecie(self._symbol, self._oxi_state)

    @staticmethod
    def from_string(species_string):
        """
        Returns a Dummy from a string representation.

        Args:
            species_string:
                A string representation of a dummy species, e.g., "X2+", "X3+"

        Returns:
            A DummySpecie object.

        Raises:
            ValueError if species_string cannot be intepreted.
        """
        m = re.search("([A-Z][a-z]*)([0-9\.]*)([\+\-]*)", species_string)

        if m:
            if m.group(2) == "" and m.group(3) == "":
                return DummySpecie(m.group(1))
            else:
                num = 1 if m.group(2) == "" else float(m.group(2))
                oxi = num if m.group(3) == "+" else -num
                return DummySpecie(m.group(1), oxidation_state=oxi)
        raise ValueError("Invalid Species String")

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "element": self.symbol,
                "oxidation_state": self._oxi_state,
                "properties": self._properties}

    @staticmethod
    def from_dict(d):
        return DummySpecie(d["element"], d["oxidation_state"],
                           d.get("properties", None))


@singleton
class PeriodicTable(object):
    """
    A Periodic table singleton class. This class contains methods on the
    collection of all known elements. For example, printing all elements, etc.
    """

    def __init__(self):
        """ Implementation of the singleton interface """
        self._all_elements = dict()
        for sym in _pt_data.keys():
            el = Element(sym)
            self._all_elements[sym] = el

    def __getattr__(self, name):
        return self._all_elements[name]

    def __iter__(self):
        for sym in _z2symbol:
            if sym is not None:
                yield self._all_elements[sym]

    def __getitem__(self, Z_or_slice):
        #print Z_or_slice, symbol_from_Z(Z_or_slice)
        try:
            if isinstance(Z_or_slice, slice):
                return [self._all_elements[sym]
                        for sym in symbol_from_Z(Z_or_slice)]
            else:
                return self._all_elements[symbol_from_Z(Z_or_slice)]
        except:
            raise IndexError("Z_or_slice: %s" % str(Z_or_slice))

    @property
    def all_elements(self):
        """
        Returns the list of all known elements as Element objects.
        """
        return self._all_elements.values()

    def print_periodic_table(self, filter_function=None):
        """
        A pretty ASCII printer for the periodic table, based on some
        filter_function.

        Args:
            filter_function:
                A filtering function taking an Element as input and returning
                a boolean. For example, setting
                filter_function = lambda el: el.X > 2 will print
                a periodic table containing only elements with
                electronegativity > 2.
        """
        for row in range(1, 10):
            rowstr = []
            for group in range(1, 19):
                el = Element.from_row_and_group(row, group)
                if el and ((not filter_function) or filter_function(el)):
                    rowstr.append("{:3s}".format(el.symbol))
                else:
                    rowstr.append("   ")
            print(" ".join(rowstr))


def smart_element_or_specie(obj):
    """
    Utility method to get an Element or Specie from an input obj.
    If obj is in itself an element or a specie, it is returned automatically.
    If obj is an int or a string representing an integer, the Element 
    with the atomic number obj is returned.
    If obj is a string, Specie parsing will be attempted (e.g., Mn2+), failing
    which Element parsing will be attempted (e.g., Mn), failing which
    DummyElement parsing will be attempted.

    Args:
        obj:
            An arbitrary object.  Supported objects are actual Element/Specie
            objects, integers (representing atomic numbers) or strings (element
            symbols or species strings).

    Returns:
        Specie or Element, with a bias for the maximum number of properties
        that can be determined.

    Raises:
        ValueError if obj cannot be converted into an Element or Specie.
    """
    if isinstance(obj, (Element, Specie, DummySpecie)):
        return obj

    def string_is_int(s):
        "True is string s represents an integer (with sign)"
        if s[0] in ('-', '+'):
            return s[1:].isdigit()
        return s.isdigit()

    obj = str(obj)

    if string_is_int(obj):
        return Element.from_Z(int(obj))

    try:
        return Specie.from_string(obj)
    except (ValueError, KeyError):
        try:
            return Element(obj)
        except (ValueError, KeyError):
            try:
                return DummySpecie.from_string(obj)
            except:
                raise ValueError("Can't parse Element or String from " +
                                 str(obj))
