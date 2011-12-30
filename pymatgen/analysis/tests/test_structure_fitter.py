import unittest
import os
import logging
import numpy as np

from pymatgen.analysis.structure_fitter import StructureFitter, shear_invariant, sqrt_matrix
from pymatgen.core.periodic_table import Element
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure_modifier import SupercellMaker
from pymatgen.io.cifio import CifParser

module_dir = os.path.dirname(os.path.abspath(__file__))

class StructureFitterTest(unittest.TestCase):
    
    def setUp(self):
        si = Element("Si")
        fe = Element("Fe")
        logging.basicConfig(level=logging.DEBUG)
        coords = list()
        
        coords.append(np.array([0.75,0.5,0.2]))
        coords.append(np.array([0,0,0]))
        
        lattice = Lattice(np.array([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]]))
        self.a = Structure(lattice,[fe,si],coords)
        coords = list()
        coords.append(np.array([0.75,0.5,0.2]))
        
        coords.append(np.array([0,0,0]))
        lattice = Lattice(np.array([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]]))
        self.b = Structure(lattice,[fe,si],coords)
                
    def test_init(self):
        fitter = StructureFitter(self.b,self.a)
        self.assertTrue(fitter.mapping_op != None, "No fit found!")
        
        #Now to try with rotated structure
        op = SymmOp.from_axis_angle_and_translation([0, 0, 1], 30, False, np.array([0,0,1]))
        editor = StructureEditor(self.a)
        editor.apply_operation(op)
        fitter = StructureFitter(self.b,editor.modified_structure)
        
        self.assertTrue(fitter.mapping_op != None, "No fit found!")
        
        #test with a supercell
        mod = SupercellMaker(self.a,scaling_matrix = [[2,0,0],[0,1,0],[0,0,1]])
        a_super = mod.modified_structure
        fitter = StructureFitter(self.b,a_super)
        self.assertTrue(fitter.mapping_op != None, "No fit found!")
        
        # Test with a structure with a translated point

        editor = StructureEditor(self.a)
        site = self.a[0]
        editor.delete_site(0)
        trans = np.random.randint(0,1000,3)
        editor.insert_site(0, site.species_and_occu, site.frac_coords + trans, False, False)
        fitter = StructureFitter(self.b,editor.modified_structure)
        self.assertTrue(fitter.mapping_op != None, "No fit found for translation {}!".format(trans))
        
        parser = CifParser(os.path.join(module_dir,"FePO4a.cif"))
        a = parser.get_structures()[0]
        parser = CifParser(os.path.join(module_dir,"FePO4b.cif"))
        b = parser.get_structures()[0]
        fitter = StructureFitter(b,a)
        self.assertTrue(fitter.mapping_op != None, "No fit found!")
        
    def test_anonymized_fitting(self):
        parser = CifParser(os.path.join(module_dir,"LiFePO4.cif"))
        a = parser.get_structures()[0]
        parser = CifParser(os.path.join(module_dir,"NaFePO4.cif"))
        b = parser.get_structures()[0]
        fitter = StructureFitter(b,a)
        self.assertTrue(fitter.mapping_op == None, "No fit should be found when NaFePO4 and LiFePo4 are fitted in non-anonymized mode!")
        fitter = StructureFitter(b,a,anonymized = True)
        self.assertTrue(fitter.mapping_op != None, "Fit should be found when NaFePO4 and LiFePo4 are fitted in anonymized mode!")
        self.assertEqual({el1.symbol:el2.symbol for el1, el2 in fitter.el_mapping.items()}, {"O":"O", "Fe":"Fe", "Na":"Li", "P":"P"})

class SupportFunctionTest(unittest.TestCase):
    
    def test_shear_invariant(self):
        mat = np.array([[1.0000002761907518,1.5947378062541873E-7,-2.2550605566218351E-7],[1.5947378062541873E-7,0.9999997238033649,-3.906039923173843E-7],[-2.2550605566218351E-7,-3.906039922896287E-7,1.0000000000061875]])
        self.assertAlmostEqual(0, shear_invariant(mat), 7)
    
    def test_sqrt_matrix(self):
        mat = np.array([[0.1,0, 0],[0.2,0.3,0.1],[0.4,0.7,0.9]])
        expected_ans = np.array([[ 0.09813965, -0.13085287,  0.06542643], [-0.13085287,  0.41807729, -0.17084204], [ 0.06542643, -0.17084204,  1.24722442]])
        self.assertTrue((abs(sqrt_matrix(mat) - expected_ans) < 0.00001).all())
    
if __name__ == '__main__':
    unittest.main()
