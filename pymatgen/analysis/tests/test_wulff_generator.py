# coding: utf-8
# !/usr/bin/env python

__author__ = 'Zihan Xu, Richard Tran, Balachandran Radhakrishnan'
__copyright__ = 'Copyright 2013, The Materials Virtual Lab'
__version__ = '0.1'
__maintainer__ = 'Zihan Xu'
__email__ = 'zix009@eng.ucsd.edu'
__date__ = 'May 05 2016'

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import PymatgenTest
from pymatgen.util.coord_utils import in_coord_list
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

from pymatgen.analysis.wulff_generator import wulff_3d


class TestWulff(PymatgenTest):

    def setUp(self):

        # In the case of a high anisotropy material
        # Be: mp-20
        latt_Be = Lattice.cubic(2.5515)
        miller_list_Be = [[1, 1, 0], [2, 1, 1], [1, 1, 1],
                          [2, 2, 1], [2, 1, 0], [1, 0, 0]]
        e_surf_list_Be = [1.0162, 1.6120, 1.8737,
                          1.3692, 1.4282, 1.7367]

        # In the case of an fcc material
        # Ir: mp-101
        latt_Ir = Lattice.cubic(3.8312)
        e_surf_list_Ir = [2.7720, 3.0368, 2.8177,
                          2.7091, 2.5863, 2.2845]
        miller_list_Ir = [[1, 1, 0], [2, 1, 0], [1, 0, 0],
                          [2, 1, 1], [2, 2, 1], [1, 1, 1]]

        # In the case of a hcp material
        # Ti: mp-72
        latt_Ti = Lattice.hexagonal(4.6000, 2.8200)
        e_surf_list_Ti = [2.035417074976807, 2.138851987454173, 2.263788852611296,
                          2.1434230139470034, 2.101264869212902, 2.152215199900508,
                          2.0055604698154053, 2.2190341142135015, 2.0779772084822206,
                          2.155621801227633, 2.251538932874714, 1.9318734349462858]
        miller_list_Ti = [[1, 1, -2, 0], [2, 1, -3, 0], [2, 0, -2, 1], [2, 2, -4, 1],
                          [2, 1, -3, 1], [0, 0, 0, 1], [2, -1, -1, 2], [1, 0, -1, 0],
                          [2, 1, -3, 2], [1, 0, -1, 2], [1, 0, -1, 1], [1, 1, -2, 1]]


        self.ucell_Be = Structure(latt_Be, ["Be", "Be"],
                                  [[0, 0, 0], [0.5, 0.5, 0.5]])
        self.wulff_Be = wulff_3d(latt_Be, miller_list_Be, e_surf_list_Be)

        self.ucell_Ir = Structure(latt_Ir, ["Ir", "Ir", "Ir", "Ir"],
                                  [[0, 0, 0], [0, 0.5, 0.5],
                                   [0.5, 0, 0.5], [0.5, 0.5, 0]])
        self.wulff_Ir = wulff_3d(latt_Ir, miller_list_Ir, e_surf_list_Ir)

        self.ucell_Ti = Structure(latt_Ti, ["Ti", "Ti", "Ti"],
                                  [[0, 0, 0], [0.333333, 0.666667, 0.5],
                                   [0.666667, 0.333333, 0.5]])
        self.wulff_Ti = wulff_3d(latt_Ti, miller_list_Ti, e_surf_list_Ti)

    def symm_check(self, ucell, wulff_vertices):
        """
        # Checks if the point group of the Wulff shape matches
        # the point group of its conventional unit cell

        Args:
            ucell (string): Unit cell that the Wulff shape is based on.
            wulff_vertices (list): List of all vertices on the Wulff
                shape. Use wulff.wulff_pt_list to obtain the list
                (see wulff_generator.py).

        return (bool)
        """

        space_group_analyzer = SpacegroupAnalyzer(ucell)
        symm_ops = space_group_analyzer.get_point_group_operations(cartesian=True)
        for point in wulff_vertices:
            for op in symm_ops:
                symm_point = op.operate(point)
                if in_coord_list(wulff_vertices, symm_point):
                    continue
                else:
                    return False
        return True

    def consistency_tests(self):

        # For a set of given values, these tests will
        # ensure that the general result given by the
        # algorithm does not change as the code is editted

        # For an fcc structure make sure the (111) direction
        # is the most dominant facet on the Wulff shape

        fractional_areas = self.wulff_Ir.area_fraction_dict
        miller_list = [hkl for hkl in fractional_areas.keys()]
        area_list = [fractional_areas[hkl] for hkl in fractional_areas.keys()]
        self.assertEqual(miller_list[area_list.index(max(area_list))], (1,1,1))

        # Overall weighted surface energy of bcc Be should be
        # equal to the energy of the 110 surface, ie. bcc Be
        # is anisotropic, the 110 surface is so low in energy,
        # its the only facet that exists in the Wulff shape

        Be_area_fraction_dict = self.wulff_Be.area_fraction_dict
        for hkl in Be_area_fraction_dict.keys():
            if hkl == (1,1,0):
                self.assertEqual(Be_area_fraction_dict[hkl], 1)
            else:
                self.assertEqual(Be_area_fraction_dict[hkl], 0)

        self.assertEqual(self.wulff_Be.miller_energy_dict[(1,1,0)],
                         self.wulff_Be.weighted_surface_energy)


    def symmetry_test(self):

        # Maintains that all wulff shapes have the same point
        # groups as the conventional unit cell they were
        # derived from. This test should pass for all subsequent
        # updates of the surface_properties collection

        check_symmetry_Be = self.symm_check(self.ucell_Be,
                                         self.wulff_Be.wulff_pt_list)
        check_symmetry_Ir = self.symm_check(self.ucell_Ir,
                                            self.wulff_Ir.wulff_pt_list)
        check_symmetry_Ti = self.symm_check(self.ucell_Ti,
                                            self.wulff_Ti.wulff_pt_list)
        self.assertTrue(check_symmetry_Be)
        self.assertTrue(check_symmetry_Ir)
        self.assertTrue(check_symmetry_Ti)

    def test_get_azimuth_elev(self):
        azim, elev = self.wulff_Ir.get_azimuth_elev((0,0,1))
        self.assertEqual(azim, 0)
        self.assertEqual(elev, 90)
        azim, elev = self.wulff_Ir.get_azimuth_elev((1,1,1))
        self.assertAlmostEqual(azim, 45)

    def test_properties(self):
        anisotropy_Ir = self.wulff_Ir.anisotropy
        anisotropy_Be = self.wulff_Be.anisotropy
        anisotropy_Ti = self.wulff_Ti.anisotropy

        weighted_surf_Ir = self.wulff_Ir.weighted_surface_energy
        weighted_surf_Be = self.wulff_Be.weighted_surface_energy
        weighted_surf_Ti = self.wulff_Ti.weighted_surface_energy

        self.assertEqual(round(anisotropy_Be, 3), 0.032)
        self.assertEqual(round(weighted_surf_Be, 2), 1.90)
        self.assertEqual(round(anisotropy_Ir, 3), 0.084)
        self.assertEqual(round(weighted_surf_Ir, 2), 2.42)
        self.assertEqual(round(anisotropy_Ti, 3), 0.038)
        self.assertEqual(round(weighted_surf_Ti, 2), 2.00)


