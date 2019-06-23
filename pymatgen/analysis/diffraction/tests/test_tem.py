import unittest
import TEM
from TEM import TEMCalculator
import importlib
import pymatgen as pmg
from pymatgen import Lattice, Structure

from pymatgen.util.testing import PymatgenTest
import matplotlib as mpl

assertions = unittest.TestCase('__init__')

c = TEMCalculator()
a = 4.209
latt = Lattice.cubic(a)
cubic = Structure(latt, ["Cs", "Cl"], [[0,0,0],[0.5,0.5,0.5]])

def test_wavelength_rel():
    # Tests that the relativistic wavelength formula is correct
    assertions.assertAlmostEqual(c.wavelength_rel(), 0.00197*10**-9)

def test_zone_axis_filter():
	# Tests that the appropriate Laue-Zoned points are returned
	empty_points = []
	assertions.assertEqual(c.zone_axis_filter(empty_points), [])

	points = [(-1,-1,-1)]
	assertions.assertEqual(c.zone_axis_filter(points), [])


	laue_1 = [(0,0,1), (0,1,0), (1,0,0), (0,0,-1)]
	assertions.assertEqual(c.zone_axis_filter(laue_1, 1), [(0,0,1)])

def test_get_interplanar_spacings():
	# Tests that the appropriate interplacing spacing is returned
	# TODO: write 6 other tests for the 6 other crystal classes 
	point = [(3,9,0)]
	spacings = c.get_interplanar_spacings(cubic, point)
	for p in point:
		assertions.assertAlmostEqual(spacings[p], 0.4436675557216236)

def test_bragg_angles():
	# Tests that the appropriate bragg angle is returned
	point = [(-10,3,0)]
	spacings = c.get_interplanar_spacings(cubic, point)
	angles = c.bragg_angles(spacings)
	for p in angles:
		assertions.assertAlmostEqual(angles[p], 2.4417132161608178e-12)

def test_get_s2():
	# Tests that the appropriate s2 factor is
	point = [(-10, 3, 0)]
	spacings = c.get_interplanar_spacings(cubic, point)
	angles = c.bragg_angles(spacings)
	s2 = c.get_s2(cubic, angles)
	for p in s2:
		assertions.assertAlmostEqual(s2[p], 1.5381852947115047)

def test_x_ray_factors():
	point = [(-10, 3, 0)]
	spacings = c.get_interplanar_spacings(cubic, point)
	angles = c.bragg_angles(spacings)
	x_ray = c.x_ray_factors(cubic, angles)
	assertions.assertAlmostEqual(x_ray['Cs'][(-10,3,0)], 14.42250869579648)
	assertions.assertAlmostEqual(x_ray['Cl'][(-10,3,0)], 2.7804915737999103)

def test_electron_scattering_factors():
	point = [(-10, 3, 0)]
	spacings = c.get_interplanar_spacings(cubic, point)
	angles = c.bragg_angles(spacings)
	elscatt = c.electron_scattering_factors(cubic, angles)
	assertions.assertAlmostEqual(elscatt['Cs'][(-10,3,0)], 3.0228630359546355e-09)
	assertions.assertAlmostEqual(elscatt['Cl'][(-10,3,0)], 1.0592972859944387e-09)	

def test_cell_scattering_factors():
	point = [(-10, 3, 0)]
	spacings = c.get_interplanar_spacings(cubic, point)
	angles = c.bragg_angles(spacings)
	cellscatt = c.cell_scattering_factors(cubic, angles)
	assertions.assertAlmostEqual(cellscatt[(-10,3,0)], 1.963565749960197e-09+2.0769354180385505e-24j)

def test_cell_intensity():
	point = [(-10, 3, 0)]
	spacings = c.get_interplanar_spacings(cubic, point)
	angles = c.bragg_angles(spacings)
	cellint = c.cell_intensity(cubic, angles)
	assertions.assertAlmostEqual(cellint[(-10,3,0)], 3.85559045441675e-18)	

def test_get_pattern():
    tem = c.get_pattern(cubic, two_theta_range=(0, 90))
    assertions.assertAlmostEqual(tem.x[0], 4.67747420e-13)
    assertions.assertAlmostEqual(tem.y[0], 2.44434019e+01)
    assertions.assertAlmostEqual(tem.d_hkls[0], 3.307473724811734e-12)
    assertions.assertAlmostEqual(tem.x[1], 6.61494745e-13)
    assertions.assertAlmostEqual(tem.y[1], 100)
    assertions.assertAlmostEqual(tem.d_hkls[1], 3.1464489680415293e-12)

def test_is_parallel():
	assertions.assertTrue(c.is_parallel((1,0,0), (3,0,0)))
	assertions.assertFalse(c.is_parallel((1,0,0), (3,0,1)))

if __name__ == '__main__':
    unittest.main()
