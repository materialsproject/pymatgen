#!/usr/bin/env python

"""
Module contains classes presenting Element and Specie (Element + oxidation state).
"""

__author__="Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import os
import re
import yaml

from pymatgen.core.design_patterns import singleton
from pymatgen.util.string_utils import formula_double_format

def load_periodic_table_data():
    module_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(module_dir,"periodic_table.yaml")) as f:
        return yaml.load(f)
    
PERIODIC_TABLE_DATA  = load_periodic_table_data()
PERIODIC_TABLE_ROW_SIZES = (2,8,8,18,18,32,32)

class Element(object):
    '''
    Basic immutable element object with all relevant properties.
    Only one instance of Element for each symbol is stored after creation, 
    ensuring that a particular element behaves like a singleton.
    '''
    
    _all_elements = {}

    def __new__(cls, *args):
        '''
        Prevents multiple Element objects from being created.  Faster and reduces memory usage.
        '''
        if tuple(args) in Element._all_elements:
            return Element._all_elements[tuple(args)]
        else:
            return super(Element, cls).__new__(cls)

    def __init__(self,symbol):
        '''
        Create immutable element.
        List of attributes are:
            X, Z, max_oxidation_state, mendeleev_no, min_oxidation_state, name, symbol
        '''
        self._data = PERIODIC_TABLE_DATA[symbol]
        self._z = self._data['atomic no']
        self._symbol = symbol
        
        Element._all_elements[tuple([symbol])] = self
        
    @property
    def Z(self):
        return self._z
    
    @property
    def symbol(self):
        return self._symbol
    
    @property
    def X(self):
        if 'X' in self._data:
            return self._data['X']
        else:
            return 0
        
    @property
    def number(self):
        return self.Z

    @property
    def name(self):
        return self._data['name']
    
    @property
    def atomic_mass(self):
        return self._data['atomic mass']
    
    @property
    def atomic_radius(self):
        return self._data['atomic radius']
    
    @property
    def max_oxidation_state(self):
        return self._data['max oxidation state']
    
    @property
    def min_oxidation_state(self):
        return self._data['min oxidation state']
    
    @property
    def mendeleev_no(self):
        return self._data['mendeleev no']
    
    @property
    def electrical_resistivity(self):
        return self._data['Electrical resistivity']

    @property
    def velocity_of_sound(self):
        return self._data['Velocity of sound']

    @property
    def reflectivity(self):
        return self._data['Reflectivity']
    
    @property
    def refractive_index(self):
        return self._data['Refractive index']

    @property
    def poissons_ratio(self):
        return self._data['Poissons ratio']
    
    @property
    def molar_volume(self):
        return self._data['Molar volume']
    
    @property
    def electronic_structure(self):
        return self._data['electronic']
    
    @property
    def thermal_conductivity(self):
        return self._data['Thermal conductivity']
    
    @property
    def boiling_point(self):
        return self._data['Boiling point']
    
    @property
    def melting_point(self):
        return self._data['Melting point']
    
    @property
    def critical_temperature(self):
        return self._data['Critical temperature']
    
    @property
    def superconduction_temperature(self):
        return self._data['Superconduction temperature']
    
    @property
    def liquid_range(self):
        return self._data['Liquid range']

    @property
    def bulk_modulus(self):
        return self._data['Bulk modulus']
    
    @property
    def youngs_modulus(self):
        return self._data['Youngs modulus']
    
    @property
    def brinell_hardness(self):
        return self._data['Brinell hardness']

    @property
    def rigidity_modulus(self):
        return self._data['Rigidity modulus']
    
    @property
    def mineral_hardness(self):
        return self._data['Mineral hardness']
    
    @property
    def vickers_hardness(self):
        return self._data['Vickers hardness']
    
    @property
    def density_of_solid(self):
        return self._data['Density of solid']

    @property
    def coefficient_of_linear_thermal_expansion(self):
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
        return (self.X - other.X)

    @staticmethod       
    def from_Z(z):
        for sym in PERIODIC_TABLE_DATA.keys():
            if Element(sym).Z == z:
                return Element(sym)
        raise ValueError("No element with this atomic number")

    @staticmethod
    def from_row_and_group(row,group):
        for sym in PERIODIC_TABLE_DATA.keys():
            el = Element(sym)
            if el.row == row and el.group == group:
                return el
        return None

    @staticmethod
    def is_valid_symbol(symbol):
        return symbol in PERIODIC_TABLE_DATA

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

        for i in xrange(len(PERIODIC_TABLE_ROW_SIZES)):
            totalEls += PERIODIC_TABLE_ROW_SIZES[i]
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
        """docstring for block
        return the block character 's,p,d,f'
        """
        block = ''
        if self.group in [1,2]:
            block = 's'
        elif self.group in xrange(13,19):
            block = 'p'
        elif (self.is_actinoid() or self.is_lanthanid()):
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


class Specie(object):
    """
    An extension of Element with an oxidation state.
    Note that while Specie does not directly inherit from Element 
    (because of certain implementation concerns due to the singleton nature of each element),
    it does inherit all Element attributes and hence function exactly as an Element would.
    """
           
    def __init__(self, symbol, oxidation_state):
        """
        Arguments:
            symbol - Element symbol, e.g., Fe
            oxidation_state - e.g., 2 or -2
        """
        self._el = Element(symbol)
        self._oxi_state = oxidation_state
        Element._all_elements[tuple([symbol, oxidation_state])] = self
    
    def __getattr__(self, attr):
        """
        Trick to make Specie inherit all Element properties.
        """
        if hasattr(self._el, attr):
            return getattr(self._el,attr) 
    
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
    
    def __cmp__(self, other):
        '''
        Sets a default sort order for atomic species by electronegativity, followed by oxidation state. 
        '''
        return (self.X - other.X) * 100 + (self.oxi_state - other.oxi_state)
        
    @property
    def oxi_state(self):
        """
        Oxidation state.
        """
        return self._oxi_state
    
    @staticmethod
    def from_string(species_string):
        """
        Returns a Specie from a string representation. 
        For example, Mn2+, Fe3+, O2-.
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

@singleton
class PeriodicTable(object):
    '''
    Periodic table singleton class.
    '''

    def __init__(self):
        """ Implementation of the singleton interface """
        self._all_elements = dict()
        for sym in PERIODIC_TABLE_DATA.keys():
            el = Element(sym)
            self._all_elements[sym] = el

    def __getattr__(self, name):
        return self._all_elements[name]
    
    def all_elements(self):
        return self._all_elements.values()

    @staticmethod
    def print_periodic_table(filter_function = None):
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