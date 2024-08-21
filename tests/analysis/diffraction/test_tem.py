"""Unit tests for TEM calculator."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.diffraction.tem import TEMCalculator
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest

__author__ = "Frank Wan, Jason Liang"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.201"
__maintainer__ = "Jason Liang"
__email__ = "fwan@berkeley.edu, yhljason@berkeley.edu"
__date__ = "2/20/20"


class TestTEMCalculator(PymatgenTest):
    def test_wavelength_rel(self):
        # Test that the relativistic wavelength formula (for 200 kV electron beam) is correct
        tem_calc = TEMCalculator()
        assert tem_calc.wavelength_rel() == approx(0.025079, rel=1e-4)

    def test_generate_points(self):
        # Test that 3d points are properly generated
        tem_calc = TEMCalculator()
        actual = tem_calc.generate_points(-1, 1)
        expected = np.array(
            [
                [-1, -1, -1],
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
                [1, 1, 1],
            ]
        )
        assert expected == approx(actual)

    def test_zone_axis_filter(self):
        # Test that the appropriate Laue-Zoned points are returned
        tem_calc = TEMCalculator()
        empty_points = np.asarray([])
        assert tem_calc.zone_axis_filter(empty_points) == []
        points = np.asarray([[-1, -1, -1]])
        assert tem_calc.zone_axis_filter(points) == []
        laue_1 = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, -1]])
        assert tem_calc.zone_axis_filter(laue_1, 1) == [(0, 0, 1)]

    def test_get_interplanar_spacings(self):
        # Test that the appropriate interplanar spacing is returned
        tem_calc = TEMCalculator()
        point = [(3, 9, 0)]
        lattice = Lattice.cubic(4.209)
        cubic = Structure(lattice, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        tet = self.get_structure("Li10GeP2S12")
        hexa = self.get_structure("Graphite")
        ortho = self.get_structure("K2O2")
        mono = self.get_structure("Li3V2(PO4)3")
        spacings_cubic = tem_calc.get_interplanar_spacings(cubic, point)
        spacings_tet = tem_calc.get_interplanar_spacings(tet, point)
        spacings_hexa = tem_calc.get_interplanar_spacings(hexa, point)
        spacings_ortho = tem_calc.get_interplanar_spacings(ortho, point)
        spacings_mono = tem_calc.get_interplanar_spacings(mono, point)
        for p in point:
            assert spacings_cubic[p] == approx(0.4436675557216236)
            assert spacings_tet[p] == approx(0.9164354445646701)
            assert spacings_hexa[p] == approx(0.19775826179547752)
            assert spacings_ortho[p] == approx(0.5072617738916)
            assert spacings_mono[p] == approx(0.84450786041677972)

    def test_bragg_angles(self):
        # Test that the appropriate bragg angle is returned. Testing formula with values of x-ray diffraction in
        # materials project.
        tem_calc = TEMCalculator()
        lattice = Lattice.cubic(4.209)
        cubic = Structure(lattice, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(1, 1, 0)]
        spacings = tem_calc.get_interplanar_spacings(cubic, point)
        bragg_angles_val = np.arcsin(1.5406 / (2 * spacings[point[0]]))
        assert bragg_angles_val == approx(0.26179, rel=1e-4)

    def test_get_s2(self):
        # Test that the appropriate s2 factor is returned.
        tem_calc = TEMCalculator()
        lattice = Lattice.cubic(4.209)
        cubic = Structure(lattice, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(-10, 3, 0)]
        spacings = tem_calc.get_interplanar_spacings(cubic, point)
        angles = tem_calc.bragg_angles(spacings)
        s2 = tem_calc.get_s2(angles)
        for p in s2:
            assert s2[p] == approx(1.5381852947115047)

    def test_x_ray_factors(self):
        tem_calc = TEMCalculator()
        lattice = Lattice.cubic(4.209)
        cubic = Structure(lattice, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(-10, 3, 0)]
        spacings = tem_calc.get_interplanar_spacings(cubic, point)
        angles = tem_calc.bragg_angles(spacings)
        x_ray = tem_calc.x_ray_factors(cubic, angles)
        assert x_ray["Cs"][-10, 3, 0] == approx(14.42250869579648)
        assert x_ray["Cl"][-10, 3, 0] == approx(2.7804915737999103)

    def test_electron_scattering_factors(self):
        # Test the electron atomic scattering factor, values approximate with
        # international table of crystallography volume C. Rounding error when converting hkl to sin(theta)/lambda.
        # Error increases as sin(theta)/lambda is smaller.
        tem_calc = TEMCalculator()
        lattice = Lattice.cubic(4.209)
        cubic = Structure(lattice, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(2, 1, 3)]
        point_nacl = [(4, 2, 0)]
        spacings = tem_calc.get_interplanar_spacings(cubic, point)
        spacings_nacl = tem_calc.get_interplanar_spacings(nacl, point_nacl)
        angles = tem_calc.bragg_angles(spacings)
        angles_nacl = tem_calc.bragg_angles(spacings_nacl)
        el_scatt = tem_calc.electron_scattering_factors(cubic, angles)
        el_scatt_nacl = tem_calc.electron_scattering_factors(nacl, angles_nacl)
        assert el_scatt["Cs"][2, 1, 3] == approx(2.848, rel=1e-3)
        assert el_scatt["Cl"][2, 1, 3] == approx(1.1305, rel=1e-3)
        assert el_scatt_nacl["Na"][4, 2, 0] == approx(0.8352, rel=1e-3)
        assert el_scatt_nacl["Cl"][4, 2, 0] == approx(1.3673, rel=1e-3)

    def test_cell_scattering_factors(self):
        # Test that fcc structure gives 0 intensity for mixed even, odd hkl.
        tem_calc = TEMCalculator()
        nacl = Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.692), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(2, 1, 0)]
        spacings = tem_calc.get_interplanar_spacings(nacl, point)
        angles = tem_calc.bragg_angles(spacings)
        cell_scatt = tem_calc.cell_scattering_factors(nacl, angles)
        assert cell_scatt[2, 1, 0] == approx(0)

    def test_cell_intensity(self):
        # Test that bcc structure gives lower intensity for h + k + l != even.
        tem_calc = TEMCalculator()
        lattice = Lattice.cubic(4.209)
        cubic = Structure(lattice, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(2, 1, 0)]
        point2 = [(2, 2, 0)]
        spacings = tem_calc.get_interplanar_spacings(cubic, point)
        spacings2 = tem_calc.get_interplanar_spacings(cubic, point2)
        angles = tem_calc.bragg_angles(spacings)
        angles2 = tem_calc.bragg_angles(spacings2)
        cell_int = tem_calc.cell_intensity(cubic, angles)
        cell_int2 = tem_calc.cell_intensity(cubic, angles2)
        assert cell_int2[2, 2, 0] > cell_int[2, 1, 0]

    def test_normalized_cell_intensity(self):
        # Test that the method correctly normalizes a value.
        tem_calc = TEMCalculator()
        lattice = Lattice.cubic(4.209)
        cubic = Structure(lattice, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        point = [(2, 0, 0)]
        spacings = tem_calc.get_interplanar_spacings(cubic, point)
        angles = tem_calc.bragg_angles(spacings)
        cell_int = tem_calc.normalized_cell_intensity(cubic, angles)
        assert cell_int[2, 0, 0] == approx(1)

    def test_is_parallel(self):
        tem_calc = TEMCalculator()
        structure = self.get_structure("Si")
        assert tem_calc.is_parallel(structure, (1, 0, 0), (3, 0, 0))
        assert not tem_calc.is_parallel(structure, (1, 0, 0), (3, 0, 1))

    def test_get_first_point(self):
        tem_calc = TEMCalculator()
        lattice = Lattice.cubic(4.209)
        points = tem_calc.generate_points(-2, 2)
        cubic = Structure(lattice, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        first_pt = tem_calc.get_first_point(cubic, points)
        assert 4.209 in first_pt.values()

    def test_interplanar_angle(self):
        # test interplanar angles. Reference values from KW Andrews,
        # Interpretation of Electron Diffraction pp70-90.
        tem_calc = TEMCalculator()
        lattice = Lattice.cubic(4.209)
        cubic = Structure(lattice, ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        phi = tem_calc.get_interplanar_angle(cubic, (0, 0, -1), (0, -1, 0))
        assert phi == approx(90)
        tet = self.get_structure("Li10GeP2S12")
        phi = tem_calc.get_interplanar_angle(tet, (0, 0, 1), (1, 0, 3))
        assert phi == approx(25.7835, rel=1e-4)
        lattice = Lattice.hexagonal(2, 4)
        hexagonal = Structure(lattice, ["Ab"], [[0, 0, 0]])
        phi = tem_calc.get_interplanar_angle(hexagonal, (0, 0, 1), (1, 0, 6))
        assert phi == approx(21.0517, rel=1e-4)

    def test_get_plot_coeffs(self):
        # Test if x * p1 + y * p2 yields p3.
        tem_calc = TEMCalculator()
        coeffs = tem_calc.get_plot_coeffs((1, 1, 0), (1, -1, 0), (2, 0, 0))
        assert_allclose([1, 1], coeffs)

    def test_get_positions(self):
        tem_calc = TEMCalculator()
        points = tem_calc.generate_points(-2, 2)
        structure = self.get_structure("Si")
        positions = tem_calc.get_positions(structure, points)
        assert positions[0, 0, 0].tolist() == [0, 0]
        # Test silicon diffraction data spot rough positions:
        # see https://www.doitpoms.ac.uk/tlplib/diffraction-patterns/printall.php
        assert_allclose([1, 0], positions[-1, 0, 0], atol=1)

    def test_tem_dots(self):
        # All dependencies in TEM_dots method are tested. Only make sure each object created is
        # the class desired.
        tem_calc = TEMCalculator()
        points = tem_calc.generate_points(-2, 2)
        structure = self.get_structure("Si")
        dots = tem_calc.tem_dots(structure, points)
        assert all(isinstance(x, tuple) for x in dots)

    def test_get_pattern(self):
        # All dependencies in get_pattern method are tested.
        # Only make sure result is a pd dataframe.
        tem_calc = TEMCalculator()
        structure = self.get_structure("Si")
        assert isinstance(tem_calc.get_pattern(structure), pd.DataFrame)

    def test_get_plot_2d(self):
        tem_calc = TEMCalculator()
        structure = self.get_structure("Si")
        assert isinstance(tem_calc.get_plot_2d(structure), go.Figure)

    def test_get_plot_2d_concise(self):
        tem_calc = TEMCalculator()
        structure = self.get_structure("Si")
        fig = tem_calc.get_plot_2d_concise(structure)
        width = fig.layout.width
        assert width == 121
        height = fig.layout.height
        assert height == 121
