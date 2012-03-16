#!/usr/bin/python

import unittest

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor, SupercellMaker, OxidationStateDecorator, OxidationStateRemover
import numpy as np

class StructureEditorTest(unittest.TestCase):
    
    def setUp(self):

        self.si = Element("Si")
        self.fe = Element("Fe")
        coords = list()
        coords.append(np.array([0,0,0]))
        coords.append(np.array([0.75,0.5,0.75]))
        lattice = Lattice.cubic(10)
        s = Structure(lattice,[self.si,self.fe],coords)
        self.modifier = StructureEditor(s)
        
    def test_modified_structure(self):
        self.modifier.translate_sites([0,1], [0.5,0.5,0.5], frac_coords = True)
        self.assertTrue(np.array_equal(self.modifier.modified_structure.frac_coords[0],np.array([ 0.5,  0.5,  0.5])))
        
        self.modifier.translate_sites([0], [0.5,0.5,0.5], frac_coords = False)
        self.assertTrue(np.array_equal(self.modifier.modified_structure.cart_coords[0],np.array([ 5.5,  5.5,  5.5])))
        
        self.modifier.append_site(self.si, [0,0.5,0])
        self.assertEqual(self.modifier.modified_structure.formula, "Fe1 Si2", "Wrong formula!")
        
        self.modifier.insert_site(1, self.si, [0,0.25,0])
        self.assertEqual(self.modifier.modified_structure.formula, "Fe1 Si3", "Wrong formula!")
        
        self.modifier.delete_site(0)
        self.assertEqual(self.modifier.modified_structure.formula, "Fe1 Si2", "Wrong formula!")
        
        self.modifier.replace_single_site(0, Element('Ge'))
        self.assertEqual(self.modifier.modified_structure.formula, "Fe1 Si1 Ge1", "Wrong formula!")
        
        self.assertRaises(ValueError, self.modifier.append_site, self.si, np.array([0,0.5,0]))
        
        d = 0.1
        pre_perturbation_sites = self.modifier.modified_structure.sites
        self.modifier.perturb_structure(distance = d)
        post_perturbation_sites = self.modifier.modified_structure.sites

        for x in pre_perturbation_sites:
            self.assertAlmostEqual(x.distance(post_perturbation_sites.next()), d, 3, "Bad perturbation distance") 
        
        

class SupercellMakerTest(unittest.TestCase):

    def setUp(self):
        si = Element("Si")
        fe = Element("Fe")
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.5,0.75])
        lattice = Lattice.cubic(10)
        s = Structure(lattice,[si,fe],coords)
        self.mod = SupercellMaker(s, [[1,1,0],[-1,1,0],[0,0,2]])
                
    def test_modified_structure(self):
        self.assertEquals(self.mod.modified_structure.formula, "Fe4 Si4", "Wrong formula!")
        
class OxidationStateRemoverTest(unittest.TestCase):

    def setUp(self):
        co_elem = Element("Co")
        o_elem = Element("O")
        co_specie = Specie("Co", 2)
        o_specie = Specie("O", -2)
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.5,0.75])
        lattice = Lattice.cubic(10)
        self.s_elem = Structure(lattice,[co_elem,o_elem],coords)
        self.s_specie = Structure(lattice,[co_specie,o_specie],coords)
        self.mod = OxidationStateRemover(self.s_specie)
        
    def test_modified_structure(self):
        mod_s = self.mod.modified_structure
        self.assertEqual(self.s_elem, mod_s, "Oxidation state remover failed") 
        
class OxidationStateDecoratorTest(unittest.TestCase):

    def setUp(self):
        si = Element("Si")
        fe = Element("Fe")
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.5,0.75])
        lattice = Lattice.cubic(10)
        self.s = Structure(lattice,[si,fe],coords)
        self.oxidation_states = {"Fe":2, "Si":-4}
        self.mod = OxidationStateDecorator(self.s, self.oxidation_states)
    
    def test_init(self):
        oxidation_states = {"Fe":2}
        self.assertRaises(ValueError, OxidationStateDecorator, self.s, oxidation_states)
        
    def test_modified_structure(self):
        mod_s = self.mod.modified_structure
        for site in mod_s:
            for k in site.species_and_occu.keys():
                self.assertEqual(k.oxi_state, self.oxidation_states[k.symbol], "Wrong oxidation state assigned!")

if __name__ == '__main__':
    unittest.main()

