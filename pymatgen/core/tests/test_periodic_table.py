#!/usr/bin/python

import unittest
from pymatgen.core.periodic_table import Element, Specie
from copy import deepcopy

class  ElementTestCase(unittest.TestCase):

    def test_init(self):
        self.assertEqual("Fe",Element("Fe").symbol, "Fe test failed")

        fictional_symbols = ["D", "T", "Zebra"]
        
        for sym in fictional_symbols:
            self.assertRaises(KeyError, Element, sym)
            
    def test_block(self):
        testsets = {'O':'p', 'Fe':'d','Li':'s', 'U':'f'}
        for k,v in testsets.items():
            self.assertEqual(Element(k).block, v)
    
    def test_full_electronic_structure(self):
        testsets = {'O':[(1, 's', 2), (2, 's', 2), (2, 'p', 4)], 'Fe':[(1, 's', 2), (2, 's', 2), (2, 'p', 6), (3, 's', 2), (3, 'p', 6), (3, 'd', 6), (4, 's', 2)]
                    ,'Li':[(1, 's', 2), (2, 's', 1)], 
                    'U':[(1, 's', 2), (2, 's', 2), (2, 'p', 6), (3, 's', 2), (3, 'p', 6), (3, 'd', 10), (4, 's', 2), (4, 'p', 6), (4, 'd', 10), (5, 's', 2), (5, 'p', 6), (4, 'f', 14), (5, 'd', 10), (6, 's', 2), (6, 'p', 6), (5, 'f', 3), (6, 'd', 1), (7, 's', 2)]}
        for k,v in testsets.items():
            self.assertEqual(Element(k).full_electronic_structure, v)
    
    def test_attributes(self):
        is_true = {("Xe", "Kr") : "is_noble_gas",
                   ("Fe", "Ni") : 'is_transition_metal',
                   ('Li', 'Cs') : 'is_alkali',
                   ('Ca', 'Mg') : 'is_alkaline',
                   ('F','Br', 'I') : 'is_halogen',
                   ('La',) : 'is_lanthanoid',
                   ('U','Pu') : 'is_actinoid',
                   ('Si','Ge') : 'is_metalloid'
                   } 
        
        for k, v in is_true.items():
            for sym in k:
                self.assertTrue(getattr(Element(sym),v), sym + ' is false')
        
        keys = ["name", "Z", "mendeleev_no", "atomic_mass", "electronic_structure", "X", "atomic_radius", "min_oxidation_state", 
                        "max_oxidation_state", "electrical_resistivity", "velocity_of_sound",
                        "reflectivity", "refractive_index", "poissons_ratio", "molar_volume", "thermal_conductivity", "melting_point", "boiling_point",
                        "liquid_range", "critical_temperature", "superconduction_temperature", 
                        "bulk_modulus", "youngs_modulus", "brinell_hardness", "rigidity_modulus", "mineral_hardness", 
                        "vickers_hardness", "density_of_solid", "coefficient_of_linear_thermal_expansion", "oxidation_states", "common_oxidation_states", 'average_ionic_radius', 'ionic_radii']
        
        #Test all elements up to Uranium
        for i in xrange(1,93):
            for k in keys:
                self.assertIsNotNone(getattr(Element.from_Z(i),k))
            el = Element.from_Z(i)
            if len(el.oxidation_states) > 0:
                self.assertEqual(max(el.oxidation_states), el.max_oxidation_state)
                self.assertEqual(min(el.oxidation_states), el.min_oxidation_state)
            
    
    def test_oxidation_states(self):
        el = Element("Fe")
        self.assertEqual(el.oxidation_states, (-2, -1, 1, 2, 3, 4, 5, 6))
        self.assertEqual(el.common_oxidation_states, (2, 3))
        
    def test_deepcopy(self):
        el1 = Element("Fe")
        el2 = Element("Na")
        list = [el1, el2]
        self.assertEqual(list, deepcopy(list), "Deepcopy operation doesn't produce exact copy of Element list")
    
class  SpecieTestCase(unittest.TestCase):

    def setUp(self):
        self.specie1 = Specie.from_string("Fe2+")
        self.specie2 = Specie("Fe", 3)
        self.specie3 = Specie("Fe", 2)
        
    def test_ionic_radius(self):
        self.assertEqual(self.specie2.ionic_radius, 78.5)
        self.assertEqual(self.specie3.ionic_radius, 92)
    
    def test_eq(self):
        self.assertEqual(self.specie1, self.specie3, "Static and actual constructor for Fe2+_ gives unequal result!")
        self.assertNotEqual(self.specie1, self.specie2, "Fe2+ should not be equal to Fe3+")    
    
    def test_cmp(self):
        self.assertTrue(self.specie1 < self.specie2, "Fe2+ should be < Fe3+")
        
    def test_attr(self):
        self.assertEqual(self.specie1.Z, 26, "Z attribute for Fe2+ should be the same as that for Element Fe.")
        
    def test_deepcopy(self):
        el1 = Specie("Fe", 4)
        el2 = Specie("Na", 1)
        list = [el1, el2]
        self.assertEqual(list, deepcopy(list), "Deepcopy operation doesn't produce exact copy of Specie list")

class  PeriodicTableTestCase(unittest.TestCase):

    def test_element(self):
        symbols = list()
        for i in range(1,102):
            el = Element.from_Z(i)
            self.assertGreater(el.atomic_mass, 0, "Atomic mass cannot be negative!")
            self.assertNotIn(el.symbol, symbols, "Duplicate symbol for " + el.symbol)
            symbols.append('"'+el.symbol+'"')
            self.assertIsNotNone(el.group,"Group cannot be none for Z="+str(i))
            self.assertIsNotNone(el.row, "Row cannot be none for Z="+str(i))
            
            #Test all properties
            all_attr = ['Z', 'symbol', 'X', 'name', 'atomic_mass', 'atomic_radius', 'max_oxidation_state', 'min_oxidation_state', 'mendeleev_no',
            'electrical_resistivity', 'velocity_of_sound','reflectivity', 'refractive_index', 'poissons_ratio', 'molar_volume' , 'electronic_structure',
            'thermal_conductivity', 'boiling_point', 'melting_point', 'critical_temperature', 'superconduction_temperature', 'liquid_range', 'bulk_modulus',
            'youngs_modulus', 'brinell_hardness', 'rigidity_modulus', 'mineral_hardness','vickers_hardness', 'density_of_solid','coefficient_of_linear_thermal_expansion']
            
            for a in all_attr:
                self.assertIsNotNone(el, a)
            
if __name__ == '__main__':
    unittest.main()
    
    

