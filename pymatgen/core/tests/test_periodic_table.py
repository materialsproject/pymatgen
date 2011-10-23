#!/usr/bin/python

import unittest
from pymatgen.core.periodic_table import Element, Specie

class  ElementTestCase(unittest.TestCase):

    def test_init(self):
        self.assertEqual("Fe",Element("Fe").symbol, "Fe test failed")

        fictional_symbols = ["D", "T", "Zebra"]
        
        for sym in fictional_symbols:
            self.assertRaises(KeyError, Element, sym)
    
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


class  SpecieTestCase(unittest.TestCase):

    def setUp(self):
        self.specie1 = Specie.from_string("Fe2+")
        self.specie2 = Specie("Fe", 3)
        self.specie3 = Specie("Fe", 2)      
    
    def test_eq(self):
        self.assertEqual(self.specie1, self.specie3, "Static and actual constructor for Fe2+_ gives unequal result!")
        self.assertNotEqual(self.specie1, self.specie2, "Fe2+ should not be equal to Fe3+")    
    
    def test_cmp(self):
        self.assertTrue(self.specie1 < self.specie2, "Fe2+ should be < Fe3+")
        
    def test_attr(self):
        self.assertEqual(self.specie1.Z, 26, "Z attribute for Fe2+ should be the same as that for Element Fe.")

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
    
    

