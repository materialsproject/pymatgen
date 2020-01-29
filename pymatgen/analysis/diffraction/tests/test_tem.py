# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Unit tests for TEM calculator.
"""

import unittest
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.analysis.diffraction.tem import TEMCalculator
from pymatgen.util.testing import PymatgenTest
import numpy as np
import pandas as pd
import plotly.graph_objs as go

__author__ = "Frank Wan, Jason Liang"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.2"
__maintainer__ = "Jason Liang"
__email__ = "fwan@berkeley.edu, yhljason@berkeley.edu"
__date__ = "1/28/20"


class XRDCalculatorTest(PymatgenTest):
    def test_wavelength_rel(self):
        # Tests that the relativistic wavelength formula is correct
        c = TEMCalculator()
        self.assertAlmostEqual(c.wavelength_rel(), 0.00197 * 10 ** -9)

    def test_generate_points(self):
        # Tests that 3d points are properly generated
        c = TEMCalculator()
        actual = c.generate_points(-1, 1)
        expected = np.array([[-1, -1, -1],
                             [-1, -1, 0],
                             [-1, -1, 1],
                             [0, -1, -1],
                             [0, -1, 0],
                             [0, -1, 1],
                             [1, -1, -1],
                             [1, -1, 0],
                             [1, -1, 1],
                             [-1, 0, -1],
                             [-1, 0, 0],
                             [-1, 0, 1],
                             [0, 0, -1],
                             [0, 0, 0],
                             [0, 0, 1],
                             [1, 0, -1],
                             [1, 0, 0],
                             [1, 0, 1],
                             [-1, 1, -1],
                             [-1, 1, 0],
                             [-1, 1, 1],
                             [0, 1, -1],
                             [0, 1, 0],
                             [0, 1, 1],
                             [1, 1, -1],
                             [1, 1, 0],
                             [1, 1, 1]])
        self.assertArrayEqual(expected, actual)

    def test_zone_axis_filter(self):
        # Tests that the appropriate Laue-Zoned points are returned
        c = TEMCalculator()
        empty_points = np.asarray([])
        self.assertEqual(c.zone_axis_filter(empty_points), [])
        points = np.asarray([[-1, -1, -1]])
        self.assertEqual(c.zone_axis_filter(points), [])
        laue_1 = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, -1]])
        self.assertEqual(c.zone_axis_filter(laue_1, 1), [(0, 0, 1)])

    def test_get_interplanar_spacings(self):
        # Tests that the appropriate interplacing spacing is returned
        c = TEMCalculator()
        point = [(3, 9, 0)]
        latt = Lattice.cubic(4.209)
        cubic = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        tet = self.get_structure("Li10GeP2S12")
        hexa = self.get_structure("Graphite")
        ortho = self.get_structure("K2O2")
        mono = self.get_structure("Li3V2(PO4)3")
        spacings_cubic = c.get_interplanar_spacings(cubic, point)
        spacings_tet = c.get_interplanar_spacings(tet, point)
        spacings_hexa = c.get_interplanar_spacings(hexa, point)
        spacings_ortho = c.get_interplanar_spacings(ortho, point)
        spacings_mono = c.get_interplanar_spacings(mono, point)
        for p in point:
            self.assertAlmostEqual(spacings_cubic[p], 0.4436675557216236)
            self.assertAlmostEqual(spacings_tet[p], 0.9164354445646701)
            self.assertAlmostEqual(spacings_hexa[p], 0.19775826179547752)
            self.assertAlmostEqual(spacings_ortho[p], 0.5072617738916)
            self.assertAlmostEqual(spacings_mono[p], 0.84450786041677972)

    def test_bragg_angles(self):
        # Tests that the appropriate bragg angle is returned.
        c = TEMCalculator()
        latt = Lattice.cubic(4.209)
        cubic = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(-10, 3, 0)]
        spacings = c.get_interplanar_spacings(cubic, point)
        angles = c.bragg_angles(spacings)
        for p in angles:
            self.assertAlmostEqual(angles[p], 2.4417132161608178e-12)

    def test_get_s2(self):
        # Tests that the appropriate s2 factor is returned.
        c = TEMCalculator()
        latt = Lattice.cubic(4.209)
        cubic = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(-10, 3, 0)]
        spacings = c.get_interplanar_spacings(cubic, point)
        angles = c.bragg_angles(spacings)
        s2 = c.get_s2(angles)
        for p in s2:
            self.assertAlmostEqual(s2[p], 1.5381852947115047)

    def test_x_ray_factors(self):
        c = TEMCalculator()
        latt = Lattice.cubic(4.209)
        cubic = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(-10, 3, 0)]
        spacings = c.get_interplanar_spacings(cubic, point)
        angles = c.bragg_angles(spacings)
        x_ray = c.x_ray_factors(cubic, angles)
        self.assertAlmostEqual(x_ray['Cs'][(-10, 3, 0)], 14.42250869579648)
        self.assertAlmostEqual(x_ray['Cl'][(-10, 3, 0)], 2.7804915737999103)

    def test_electron_scattering_factors(self):
        c = TEMCalculator()
        latt = Lattice.cubic(4.209)
        cubic = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(-10, 3, 0)]
        spacings = c.get_interplanar_spacings(cubic, point)
        angles = c.bragg_angles(spacings)
        elscatt = c.electron_scattering_factors(cubic, angles)
        self.assertAlmostEqual(elscatt['Cs'][(-10, 3, 0)], 3.0228630359546355e-09)
        self.assertAlmostEqual(elscatt['Cl'][(-10, 3, 0)], 1.0592972859944387e-09)

    def test_cell_scattering_factors(self):
        c = TEMCalculator()
        latt = Lattice.cubic(4.209)
        cubic = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(-10, 3, 0)]
        spacings = c.get_interplanar_spacings(cubic, point)
        angles = c.bragg_angles(spacings)
        cellscatt = c.cell_scattering_factors(cubic, angles)
        self.assertAlmostEqual(cellscatt[(-10, 3, 0)], 1.963565749960197e-09 + 2.0769354180385505e-24j)

    def test_cell_intensity(self):
        c = TEMCalculator()
        latt = Lattice.cubic(4.209)
        cubic = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(-10, 3, 0)]
        spacings = c.get_interplanar_spacings(cubic, point)
        angles = c.bragg_angles(spacings)
        cellint = c.cell_intensity(cubic, angles)
        self.assertAlmostEqual(cellint[(-10, 3, 0)], 3.85559045441675e-18)

    def test_normalized_cell_intensity(self):
        c = TEMCalculator()
        latt = Lattice.cubic(4.209)
        cubic = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(-10, 3, 0)]
        spacings = c.get_interplanar_spacings(cubic, point)
        angles = c.bragg_angles(spacings)
        cellint = c.normalized_cell_intensity(cubic, angles)
        self.assertAlmostEqual(cellint[(-10, 3, 0)], 1)

    def test_is_parallel(self):
        c = TEMCalculator()
        self.assertTrue(c.is_parallel((1, 0, 0), (3, 0, 0)))
        self.assertFalse(c.is_parallel((1, 0, 0), (3, 0, 1)))

    def test_get_first_point(self):
        c = TEMCalculator()
        latt = Lattice.cubic(4.209)
        points = c.generate_points(-2, 2)
        cubic = Structure(latt, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        first_pt = c.get_first_point(cubic, points)
        self.assertTrue(4.209 in first_pt.values())

    def test_get_plot_coeffs(self):
        # Test if x * p1 + y * p2 yields p3.
        c = TEMCalculator()
        coeffs = c.get_plot_coeffs((1, 1, 0), (1, -1, 0), (2, 0, 0), -2, False)
        self.assertEqual([1, 1], coeffs)

    def test_get_positions(self):
        c = TEMCalculator()
        points = c.generate_points(-2, 2)
        structure = self.get_structure("Si")
        positions = c.get_positions(structure, points)
        self.assertArrayEqual([0, 0], positions[(0, 0, 0)])
        # Test silicon diffraction data spot rough positions:
        # see https://www.doitpoms.ac.uk/tlplib/diffraction-patterns/printall.php
        self.assertArrayAlmostEqual([1, 0], positions[(-1, 0, 0)], 0)

    def test_TEM_dots(self):
        # All dependencies in TEM_dots method are tested. Only make sure each object created is
        # the class desired.
        c = TEMCalculator()
        points = c.generate_points(-2, 2)
        structure = self.get_structure("Si")
        dots = c.tem_dots(structure, points)
        self.assertTrue(all([isinstance(x, tuple) for x in dots]))

    def test_get_pattern(self):
        # All dependencies in get_pattern method are tested.
        # Only make sure result is a pd dataframe.
        c = TEMCalculator()
        structure = self.get_structure("Si")
        self.assertTrue(isinstance(c.get_pattern(structure), pd.DataFrame))

    def test_show_plot_2d(self):
        c = TEMCalculator()
        structure = self.get_structure("Si")
        self.assertTrue(isinstance(c.show_plot_2d(structure), go.Figure))

if __name__ == '__main__':
    unittest.main()
