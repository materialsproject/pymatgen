#!/usr/bin/env python

"""
Module contains classes presenting Element and Specie (Element + oxidation
state) and PeriodicTable.
"""

__author__ = "Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import os
import re
import json

from pymatgen.util.decorators import singleton, cached_class
from pymatgen.util.string_utils import formula_double_format


def _load__pt_data():
    """Loads element data from json file"""
    module_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(module_dir, "periodic_table.json")) as f:
        return json.load(f)

_pt_data = _load__pt_data()
_pt_row_sizes = (2, 8, 8, 18, 18, 32, 32)

@cached_class
class Element(object):
    '''
    Basic immutable element object with all relevant properties.
    Only one instance of Element for each symbol is stored after creation,
    ensuring that a particular element behaves like a singleton.
    '''

    def __init__(self, symbol):
        '''
        Create immutable element from a symbol.
        Args:
            symbol:
                Element symbol, e.g., "H", "Fe"
        '''
        self._data = _pt_data[symbol]
        #Store key variables for quick access
        self._z = self._data['Atomic no']
        self._symbol = symbol
        self._x = self._data.get('X', 0)

    @property
    def average_ionic_radius(self):
        """
        Average ionic radius for element in pm. The average is taken over all
        oxidation states of the element for which data is present.
        """
        if 'Ionic_radii' in self._data:
            radii = self._data['Ionic_radii']
            return sum(radii.values()) / len(radii)
        else:
            return 0

    @property
    def ionic_radii(self):
        """
        All ionic radii of the element as a dict of
        {oxidation state: ionic radii}. Radii are given in pm.
        """
        if 'Ionic_radii' in self._data:
            return {int(k):v for k, v in self._data['Ionic_radii'].items()}
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
    def name(self):
        """Full name for element"""
        return self._data['Name']

    @property
    def atomic_mass(self):
        """Atomic mass"""
        return self._data['Atomic mass']

    @property
    def atomic_radius(self):
        """Atomic radius"""
        return self._data['Atomic radius']

    @property
    def max_oxidation_state(self):
        """Maximum oxidation state for element"""
        if 'Oxidation states' in self._data:
            return max(self._data['Oxidation states'])
        return 0

    @property
    def min_oxidation_state(self):
        """Minimum oxidation state for element"""
        if 'Oxidation states' in self._data:
            return min(self._data['Oxidation states'])
        return 0

    @property
    def oxidation_states(self):
        """Tuple of all known oxidation states"""
        return tuple(self._data.get('Oxidation states', list()))

    @property
    def common_oxidation_states(self):
        """Tuple of all known oxidation states"""
        return tuple(self._data.get('Common oxidation states', list()))

    @property
    def mendeleev_no(self):
        """Mendeleev number"""
        return self._data['Mendeleev no']

    @property
    def electrical_resistivity(self):
        """Electrical resistivity"""
        return self._data['Electrical resistivity']

    @property
    def velocity_of_sound(self):
        """Velocity of sound"""
        return self._data['Velocity of sound']

    @property
    def reflectivity(self):
        """Reflectivity"""
        return self._data['Reflectivity']

    @property
    def refractive_index(self):
        """Refractice index"""
        return self._data['Refractive index']

    @property
    def poissons_ratio(self):
        """Poisson's ratio"""
        return self._data['Poissons ratio']

    @property
    def molar_volume(self):
        """Molar volume"""
        return self._data['Molar volume']

    @property
    def electronic_structure(self):
        """Electronic structure. Simplified form with HTML formatting.
        E.g., The electronic structure for Fe is represented as [Ar].3d<sup>6</sup>.4s<sup>2</sup>
        """
        return self._data['Electronic structure']

    @property
    def full_electronic_structure(self):
        """
        Full electronic structure as tuple.
        E.g., The electronic structure for Fe is represented as:
        [(1, 's', 2), (2, 's', 2), (2, 'p', 6), (3, 's', 2), (3, 'p', 6), (3, 'd', 6), (4, 's', 2)]
        """
        estr = self._data['Electronic structure']
        def parse_orbital(orbstr):
            m = re.match("(\d+)([spdfg]+)<sup>(\d+)</sup>", orbstr)
            if m:
                return (int(m.group(1)), m.group(2), int(m.group(3)))
            return orbstr
        data = [parse_orbital(s) for s in estr.split(".")]
        if data[0][0] == "[":
            data = Element(data[0].replace("[", "").replace("]", "")).full_electronic_structure + data[1:]
        return data
    
    @property
    def thermal_conductivity(self):
        """Thermal conductivity"""
        return self._data['Thermal conductivity']
    
    @property
    def boiling_point(self):
        """Boiling point"""
        return self._data['Boiling point']
    
    @property
    def melting_point(self):
        """Melting point"""
        return self._data['Melting point']
    
    @property
    def critical_temperature(self):
        """Critical temperature"""
        return self._data['Critical temperature']
    
    @property
    def superconduction_temperature(self):
        """Superconduction temperature"""
        return self._data['Superconduction temperature']
    
    @property
    def liquid_range(self):
        """Liquid range"""
        return self._data['Liquid range']

    @property
    def bulk_modulus(self):
        """Bulk modulus"""
        return self._data['Bulk modulus']
    
    @property
    def youngs_modulus(self):
        """Young's modulus"""
        return self._data['Youngs modulus']
    
    @property
    def brinell_hardness(self):
        """Brinell hardness"""
        return self._data['Brinell hardness']

    @property
    def rigidity_modulus(self):
        """Rigidity modulous"""
        return self._data['Rigidity modulus']
    
    @property
    def mineral_hardness(self):
        """Mineral hardness"""
        return self._data['Mineral hardness']
    
    @property
    def vickers_hardness(self):
        """Vicker's hardness"""
        return self._data['Vickers hardness']
    
    @property
    def density_of_solid(self):
        """Density of solid phase"""
        return self._data['Density of solid']

    @property
    def coefficient_of_linear_thermal_expansion(self):
        """Coefficient of linear thermal expansion"""
        return self._data['Coefficient of linear thermal expansion']

    def __eq__(self,other):
        if other == None:
            return False
        return self.Z == other.Z
    
    def __ne__(self,other):
        if other == None:
            return True
        return self.Z != other.Z

    def __hash__(self):
        return self.Z
    
    def __repr__(self):
        return "Element " + self.symbol

    def __str__(self):
        return self.symbol

    def __cmp__(self, other):
        '''
        Sets a default sort order for atomic species by electronegativity.  Very
        useful for getting correct formulas.  For example, FeO4PLi is automatically
        sorted in LiFePO4.
        '''
        return (self._x - other._x)
    
    def __lt__(self, other):
        '''
        Sets a default sort order for atomic species by electronegativity.  Very
        useful for getting correct formulas.  For example, FeO4PLi is automatically
        sorted in LiFePO4.
        '''
        return (self._x < other._x)

    @staticmethod       
    def from_Z(z):
        '''Get an element from an atomic number'''
        for sym in _pt_data.keys():
            if Element(sym).Z == z:
                return Element(sym)
        raise ValueError("No element with this atomic number")

    @staticmethod
    def from_row_and_group(row,group):
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
        """Returns true if symbol is a valid element symbol.
        
        Args:
            symbol:
                Element symbol
        
        Returns:
            True if symbol is a valid element (e.g., "H"). False otherwise (e.g., "Zebra").
        """
        return symbol in _pt_data

    @property
    def row(self):
        """
        Returns the periodic table row of the element.
        """
        Z = self.Z
        totalEls = 0;
        if (Z >= 57 and Z <= 70):
            return 8
        elif (Z >= 89 and Z <= 102):
            return 9

        for i in xrange(len(_pt_row_sizes)):
            totalEls += _pt_row_sizes[i]
            if totalEls >= Z:
                return i+1
        return 8

    @property
    def group(self):
        """
        Returns the periodic table group of the element.
        """
        Z = self.Z
        if Z == 1:
            return 1
        if Z == 2:
            return 18
        if Z >= 3 and Z <= 18:
            if ((Z - 2) % 8 == 0):
                return 18
            elif ((Z - 2) % 8 <= 2):
                return (Z - 2) % 8;
            else:
                return (10 + (Z - 2) % 8)

        if Z >= 19 and Z <= 54:
            if ((Z - 18) % 18 == 0):
                return 18
            else:
                return (Z - 18) % 18

        if((Z - 54) % 32 == 0):
            return 18
        elif ((Z - 54) % 32 >= 17):
            return (Z - 54) % 32 - 14
        else:
            return (Z - 54) % 32

    @property
    def block(self):
        """
        Return the block character 's,p,d,f'
        """
        block = ''
        if self.group in [1,2]:
            block = 's'
        elif self.group in xrange(13,19):
            block = 'p'
        elif (self.is_actinoid or self.is_lanthanoid):
            block = 'f'
        elif self.group in xrange(3,13):
            block = 'd'
        else:
            print("unable to determine block")
        return block

    @property
    def is_noble_gas(self):
        """
        True if element is noble gas.
        """
        ns = [2, 10, 18, 36, 54, 86, 118]
        return self.Z in ns

    @property
    def is_transition_metal(self):
        """
        True if element is a transition metal.
        """

        ns = range(21, 31)
        ns.extend(xrange(39, 49))
        ns.append(57)
        ns.extend(xrange(72, 81))
        ns.append(89)
        ns.extend(xrange(104, 113))
        return self.Z in ns
            
    @property
    def is_rare_earth_metal(self):
        """
        True if element is a rare earth metal.
        """
        return self.is_lanthanid or self.is_actinoid

    @property
    def is_metalloid(self):
        """
        True if element is a metalloid.
        """
        ns = ['B', 'Si', 'Ge', 'As', 'Sb', 'Te', 'Po']
        return self.symbol in ns
    
    @property
    def is_alkali(self):
        """
        True if element is an alkali metal.
        """
        ns = [3,11, 19,37, 55,87]
        return self.Z in ns
    
    @property
    def is_alkaline(self):
        """
        True if element is an alkaline earth metal (group II).
        """
        ns = [4, 12, 20, 38, 56, 88]
        return self.Z in ns
    
    @property
    def is_halogen(self):
        """
        True if element is a halogen.
        """
        ns = [9, 17, 35, 53, 85]
        return self.Z in ns

    @property
    def is_lanthanoid(self):
        """
        True if element is a lanthanoid.
        """
        return self.Z > 56 and self.Z < 72

    @property
    def is_actinoid(self):
        """
        True if element is a actinoid.
        """
        return self.Z > 88 and self.Z < 104
    
    def __deepcopy__(self, memo):
        return Element(self.symbol)


class Specie(Element):
    """
    An extension of Element with an oxidation state.
    
    .. note::
        While Specie does not directly inherit from Element  (because of certain implementation concerns 
        due to the singleton nature of each element), it does inherit all Element attributes and hence 
        function exactly as an Element would.
    """
           
    def __init__(self, symbol, oxidation_state):
        """
        Args:
            symbol:
                Element symbol, e.g., Fe
            oxidation_state:
                Oxidation state of element, e.g., 2 or -2
        """
        super(Specie, self).__init__(symbol)
        self._oxi_state = oxidation_state
        
    def __eq__(self,other):
        """
        Specie is equal to other only if element and oxidation states are exactly the same.
        """
        if other == None:
            return False
        return self.Z == other.Z and self.oxi_state == other.oxi_state
    
    def __ne__(self,other):
        if other == None:
            return True
        return self.Z != other.Z or self.oxi_state != other.oxi_state

    def __hash__(self):
        """
        Given that all oxidation states are below 100 in absolute value, this should effectively ensure that no two unequal 
        Specie have the same hash.
        """
        return self.Z * 100 + self.oxi_state
    
    def __lt__(self, other):
        '''
        Sets a default sort order for atomic species by electronegativity, followed by oxidation state. 
        '''
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
                A typical string representation of a species, e.g., "Mn2+", "Fe3+", "O2-".
                
        Returns:
            A Specie object.
            
        Raises:
            ValueError if species_string cannot be intepreted.
        """
        m = re.search('([A-Z][a-z]*)([0-9\.]*)([\+\-])', species_string)
        if m:
            num = 1 if m.group(2) == "" else float(m.group(2))
            return Specie(m.group(1), num if m.group(3) == "+" else - num)
        else:
            raise ValueError("Invalid Species String")
    
    def __repr__(self):
        return "Specie " + self.__str__()
    
    def __str__(self):
        output = self.symbol
        if self.oxi_state >= 0:
            output += formula_double_format(self.oxi_state) + "+"
        else:
            output += formula_double_format(- self.oxi_state) + "-"
        return output
    
    def __deepcopy__(self, memo):
        return Specie(self.symbol, self.oxi_state)

@singleton
class PeriodicTable(object):
    '''
    A Periodic table singleton class.
    This class contains methods on the collection of all known elements.  
    For example, printing all elements, etc.
    '''

    def __init__(self):
        """ Implementation of the singleton interface """
        self._all_elements = dict()
        for sym in _pt_data.keys():
            el = Element(sym)
            self._all_elements[sym] = el

    def __getattr__(self, name):
        return self._all_elements[name]
    
    @property
    def all_elements(self):
        """
        Returns the list of all known elements as Element objects.
        """
        return self._all_elements.values()

    @staticmethod
    def print_periodic_table(filter_function = None):
        """
        A pretty ASCII printer for the periodic table, based on some filter_function.
        
        Args:
            filter_function:
                A filtering function taking an Element as input and returning a boolean.
                For example, setting filter_function = lambda el: el.X > 2 will print
                a periodic table containing only elements with electronegativity > 2.
        """
        for row in range(1,10):
            for group in range(1,19):
                el = Element.from_row_and_group(row,group)
                if el != None and ((not filter_function) or filter_function(el)):
                    print "%3s" % (el.symbol),
                else:
                    print "   ",
            print

def smart_element_or_specie(obj):
    """
    Utility method to get an Element or Specie from an input obj.
    If obj is in itself an element or a specie, it is returned automatically.
    If obj is an int, the Element with the atomic number obj is returned.
    If obj is a string, Specie parsing will be attempted (e.g., Mn2+), failing which Element parsing will be attempted (e.g., Mn).
    
    Args:
        obj:
            An arbitrary object.  Supported objects are actual Element/Specie objects, 
            integers (representing atomic numbers) or strings (element symbols or species strings).
            
    Returns:
        Specie or Element, with a bias for the maximum number of properties that can be determined.
        
    Raises:
        ValueError if obj cannot be converted into an Element or Specie.
    """
    if isinstance(obj,(Element, Specie)):
        return obj
    elif isinstance(obj,int):
        return Element.from_Z(obj)
    elif isinstance(obj,basestring):
        try:
            return Specie.from_string(obj)
        except ValueError:
            return Element(obj)
    raise ValueError("Can't parse Element or String from "+ str(obj))
