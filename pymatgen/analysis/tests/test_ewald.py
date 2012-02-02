import unittest
import os

import pymatgen.io.vaspio 
from pymatgen.core.structure_modifier import OxidationStateDecorator
from pymatgen.analysis.ewald import EwaldSummation, EwaldMinimizer
from pymatgen.io.vaspio import Poscar
import numpy as np

class EwaldSummationTest(unittest.TestCase):

    def test_init(self):
        module_path = os.path.dirname(pymatgen.io.vaspio.__file__)
        filepath = os.path.join(module_path, 'tests','vasp_testfiles', 'POSCAR')
        p = Poscar.from_file(filepath)
        s = p.struct

        modifier = OxidationStateDecorator(s,{"Li":1, "Fe":2, "P":5, "O":-2})
        s = modifier.modified_structure
        ham = EwaldSummation(s)
        self.assertAlmostEqual(ham.real_space_energy, -354.91294268, 4, "Real space energy incorrect!")
        self.assertAlmostEqual(ham.reciprocal_space_energy, 25.475754801, 4, "Reciprocal space energy incorrect!")
        self.assertAlmostEqual(ham.point_energy, -790.463835033, 4, "Point space energy incorrect!")
        self.assertAlmostEqual(ham.total_energy, -1119.90102291, 2, "Total space energy incorrect!")
        self.assertAlmostEqual(sum(sum(abs(ham.forces))), 915.925354346, 4, "Forces incorrect")
        self.assertAlmostEqual(sum(sum(ham.real_space_energy_matrix)), -354.91294268, 4, "Real space energy matrix incorrect!")
        self.assertAlmostEqual(sum(sum(ham.reciprocal_space_energy_matrix)), 25.475754801, 4, "Reciprocal space energy matrix incorrect!")
        self.assertAlmostEqual(sum(ham.point_energy_matrix), -790.463835033, 4, "Point space energy matrix incorrect!")
        self.assertAlmostEqual(sum(sum(ham.total_energy_matrix)), -1119.90102291, 2, "Total space energy matrix incorrect!")
        #note that forces are not individually tested, but should work fine.
        

if __name__ == "__main__":
    unittest.main()
    
class EwaldMinimizerTest(unittest.TestCase):
    
    def test_init(self):
        matrix = np.array([[ -3.,  3.,  4., -0.,  3.,  3.,  1., 14.,  9., -4.],
                           [  1., -3., -3., 12., -4., -1.,  5., 11.,  1., 12.],
                           [ 14.,  7., 13., 15., 13.,  5., -5., 10., 14., -2.],
                           [  9., 13.,  4.,  1.,  3., -4.,  7.,  0.,  6., -4.],
                           [  4., -4.,  6.,  1., 12., -4., -2., 13.,  0.,  6.],
                           [ 13.,  7., -4., 12., -2.,  9.,  8., -5.,  3.,  1.],
                           [  8.,  1., 10., -4., -2.,  4., 13., 12., -3., 13.],
                           [  2., 11.,  8.,  1., -1.,  5., -3.,  4.,  5.,  0.],
                           [ -0., 14.,  4.,  3., -1., -5.,  7., -1., -1.,  3.],
                           [  2., -2., 10.,  1.,  6., -5., -3., 12.,  0., 13.]])
        
        m_list = [[.9,4,[1,2,3,4,6,8],'a'],[-1,2,[5,6,7],'b']]
        
        e_min = EwaldMinimizer(matrix, m_list)
        
        self.assertAlmostEqual(e_min.minimized_sum, 111.63, 3, "Returned wrong minimum value")
        self.assertEqual(len(e_min.best_m_list), 6, "Returned wrong number of permutations")
        
        
        