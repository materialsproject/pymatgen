#!/usr/bin/python

import unittest
from pymatgen.core.lattice import Lattice
from numpy import array
     
class  LatticeTestCase(unittest.TestCase):
    
    def setUp(self):
        self.lattice = Lattice.cubic(10.0)
        self.tetragonal = Lattice.tetragonal(10, 20)

    def test_initialization(self):
        a = 9.026
        lattice = Lattice.cubic(a)
        self.assertIsNotNone(lattice, "Initialization from new_cubic failed")
        lattice2 = Lattice([[a, 0, 0], [0, a, 0], [0, 0, a]])
        for i in range(0,3):
            for j in range(0,3):
                self.assertAlmostEqual(lattice.matrix[i][j], lattice2.matrix[i][j], 5, "Inconsistent matrix from two inits!")
        primlatt = lattice.get_primitive_lattice('I')
        (lengths,angles) = primlatt.lengths_and_angles
        for i in range(0,3):
            self.assertAlmostEqual(lengths[i], 7.81674529, 5, "Wrong primitive lattice obtained!")
            self.assertAlmostEqual(angles[i], 109.47122063, 5, "Wrong primitive lattice obtained!")
        
        coord = lattice.get_cartesian_coords(array([0.5,0.5,0.5]))
        prim_frac = primlatt.get_fractional_coords(coord)
        for i in range(0,3):
            self.assertAlmostEqual(coord[i], 4.513, 5, "Wrong coord!")
            self.assertAlmostEqual(prim_frac[i], 1, 5, "Wrong primitive frac coord!")
    
    def test_static_methods(self):
        
        lengths_c=[ 3.840198,    3.84019885,  3.8401976 ]
        angles_c= [ 119.99998575,   90,           60.00000728]
        mat_c = array([[3.840198,0.000000,0.000000],[1.920099,3.325710,0.000000],[0.000000,-2.217138,3.135509]]) #should give the lengths and angles above
        newlatt = Lattice(mat_c)
        (lengths,angles) = newlatt.lengths_and_angles
        for i in range(0,3):
            self.assertAlmostEqual(lengths[i], lengths_c[i], 5, "Lengths incorrect!")
            self.assertAlmostEqual(angles[i], angles_c[i], 5, "Angles incorrect!")
        (lengths,angles) = Lattice.from_lengths_and_angles(lengths,angles).lengths_and_angles
        for i in range(0,3):
            self.assertAlmostEqual(lengths[i], lengths_c[i], 5, "Lengths incorrect!")
            self.assertAlmostEqual(angles[i], angles_c[i], 5, "Angles incorrect!")

    def test_attributes(self):
        """docstring for test_attributes"""
        lattice = Lattice.cubic(10.0)
        self.assertEqual(lattice.a, 10.0)
        self.assertEqual(lattice.b, 10.0)
        self.assertEqual(lattice.c, 10.0)
        self.assertAlmostEqual(lattice.volume, 1000.0)
        xyz = lattice.get_cartesian_coords([0.25, 0.35, 0.45])
        self.assertEqual(xyz[0], 2.5)
        self.assertEqual(xyz[1], 3.5)
        self.assertEqual(xyz[2], 4.5)
        
    def test_consistency(self):
        '''when only lengths and angles are given for constructors, the internal matrix representation is ambiguous since the lattice rotation is not specified.
        This test makes sure that a consistent definition is specified for the lattice rotation when using different constructors from lengths angles
        '''
        l=[ 3.840198,    3.84019885,  3.8401976 ]
        a= [ 119.99998575,   90,           60.00000728]
        mat1=Lattice.from_lengths_and_angles(l, a).matrix
        mat2=Lattice.from_parameters(l[0], l[1], l[2], a[0], a[1], a[2]).matrix
        for i in range(0,3):            
            for j in range(0,3):
                self.assertAlmostEqual(mat1[i][j], mat2[i][j], 5, "Lattice constructors do not define a consistent rotation matrix")

        
        
        
if __name__ == '__main__':
    unittest.main()

