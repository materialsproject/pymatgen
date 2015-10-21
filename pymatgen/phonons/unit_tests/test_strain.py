import unittest
from pymatgen.phonons.strain import Strain

class StrainTest(unittest.TestCase):

    def setUp(self):
        self.strain = Strain([[ 1., 0., 0.], 
                              [ 0., 1., 0.], 
                              [ 0., 1.000001, 1.]])

    def testreturn(self):
        self.assertEqual(self.matrix1.all(), self.s.deformation_matrix.all())
        
if __name__ == '__main__':
    unittest.main()
