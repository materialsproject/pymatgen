from __future__ import division, unicode_literals

__author__ = 'Xiang-Guo Li'
__copyright__ = 'Copyright 2018, The Materials Virtual Lab'
__email__ = 'xil110@eng.ucsd.edu'
__date__ = '07/30/18'

from pymatgen.util.testing import PymatgenTest
import os
import numpy as np
from pymatgen import Structure
from pymatgen.analysis.gb.grain import GrainBoundary, GrainBoundaryGenerator

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "grain_boundary")


class Test_GrainBoundary(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        cls.Cu_conv = Structure.from_file(os.path.join(test_dir,
                                                       "Cu_mp-30_conventional_standard.cif"))
        GB_Cu_conv = GrainBoundaryGenerator(cls.Cu_conv)
        cls.Cu_GB1 = GB_Cu_conv.gb_from_parameters([1, 2, 3], 123.74898859588858,
                                                   expand_times=4, vacuum_thickness=1.5,
                                                   ab_shift=[0.0, 0.0], plane=[1, 3, 1],
                                                   rm_ratio=0.0)
        cls.Cu_GB2 = GB_Cu_conv.gb_from_parameters([1, 2, 3], 123.74898859588858,
                                                   expand_times=4, vacuum_thickness=1.5,
                                                   ab_shift=[0.2, 0.2], rm_ratio=0.0)

    def test_init(self):
        self.assertAlmostEqual(self.Cu_GB1.rotation_angle, 123.74898859588858)
        self.assertAlmostEqual(self.Cu_GB1.vacuum_thickness, 1.5)
        self.assertListEqual(self.Cu_GB2.rotation_axis, [1, 2, 3])
        self.assertArrayAlmostEqual(np.array(self.Cu_GB1.ab_shift),
                                    np.array([0.0, 0.0]))
        self.assertArrayAlmostEqual(np.array(self.Cu_GB2.ab_shift),
                                    np.array([0.2, 0.2]))
        self.assertListEqual(self.Cu_GB1.gb_plane, [1, 3, 1])
        self.assertListEqual(self.Cu_GB2.gb_plane, [1, 2, 3])
        self.assertArrayAlmostEqual(np.array(self.Cu_GB1.init_cell.lattice.matrix),
                                    np.array(self.Cu_conv.lattice.matrix))

    def test_copy(self):
        Cu_GB1_copy = self.Cu_GB1.copy()
        self.assertAlmostEqual(Cu_GB1_copy.sigma, self.Cu_GB1.sigma)
        self.assertAlmostEqual(Cu_GB1_copy.rotation_angle, self.Cu_GB1.rotation_angle)
        self.assertListEqual(Cu_GB1_copy.rotation_axis, self.Cu_GB1.rotation_axis)
        self.assertListEqual(Cu_GB1_copy.gb_plane, self.Cu_GB1.gb_plane)
        self.assertArrayAlmostEqual(Cu_GB1_copy.init_cell.lattice.matrix,
                                    self.Cu_GB1.init_cell.lattice.matrix)
        self.assertArrayAlmostEqual(Cu_GB1_copy.oriented_unit_cell.lattice.matrix,
                                    self.Cu_GB1.oriented_unit_cell.lattice.matrix)
        self.assertArrayAlmostEqual(Cu_GB1_copy.lattice.matrix,
                                    self.Cu_GB1.lattice.matrix)

    def test_sigma(self):
        self.assertAlmostEqual(self.Cu_GB1.sigma, 9)
        self.assertAlmostEqual(self.Cu_GB2.sigma, 9)

    def test_top_grain(self):
        self.assertAlmostEqual(self.Cu_GB1.num_sites, self.Cu_GB1.top_grain.num_sites * 2)
        self.assertArrayAlmostEqual(self.Cu_GB1.lattice.matrix,
                                    self.Cu_GB1.top_grain.lattice.matrix)

    def test_bottom_grain(self):
        self.assertAlmostEqual(self.Cu_GB1.num_sites, self.Cu_GB1.bottom_grain.num_sites * 2)
        self.assertArrayAlmostEqual(self.Cu_GB1.lattice.matrix,
                                    self.Cu_GB1.bottom_grain.lattice.matrix)

    def test_coincidents(self):
        self.assertAlmostEqual(self.Cu_GB1.num_sites / self.Cu_GB1.sigma, len(self.Cu_GB1.coincidents))
        self.assertAlmostEqual(self.Cu_GB2.num_sites / self.Cu_GB2.sigma, len(self.Cu_GB2.coincidents))

    def test_as_dict_and_from_dict(self):
        d1 = self.Cu_GB1.as_dict()
        d2 = self.Cu_GB2.as_dict()
        Cu_GB1_new = GrainBoundary.from_dict(d1)
        Cu_GB2_new = GrainBoundary.from_dict(d2)
        self.assertAlmostEqual(Cu_GB1_new.sigma, self.Cu_GB1.sigma)
        self.assertAlmostEqual(Cu_GB1_new.rotation_angle, self.Cu_GB1.rotation_angle)
        self.assertListEqual(Cu_GB1_new.rotation_axis, self.Cu_GB1.rotation_axis)
        self.assertListEqual(Cu_GB1_new.gb_plane, self.Cu_GB1.gb_plane)
        self.assertArrayAlmostEqual(Cu_GB1_new.init_cell.lattice.matrix,
                                    self.Cu_GB1.init_cell.lattice.matrix)
        self.assertArrayAlmostEqual(Cu_GB1_new.oriented_unit_cell.lattice.matrix,
                                    self.Cu_GB1.oriented_unit_cell.lattice.matrix)
        self.assertArrayAlmostEqual(Cu_GB1_new.lattice.matrix,
                                    self.Cu_GB1.lattice.matrix)
        self.assertAlmostEqual(Cu_GB2_new.sigma, self.Cu_GB2.sigma)
        self.assertAlmostEqual(Cu_GB2_new.rotation_angle, self.Cu_GB2.rotation_angle)
        self.assertListEqual(Cu_GB2_new.rotation_axis, self.Cu_GB2.rotation_axis)
        self.assertListEqual(Cu_GB2_new.gb_plane, self.Cu_GB2.gb_plane)
        self.assertArrayAlmostEqual(Cu_GB2_new.init_cell.lattice.matrix,
                                    self.Cu_GB2.init_cell.lattice.matrix)
        self.assertArrayAlmostEqual(Cu_GB2_new.oriented_unit_cell.lattice.matrix,
                                    self.Cu_GB2.oriented_unit_cell.lattice.matrix)
        self.assertArrayAlmostEqual(Cu_GB2_new.lattice.matrix,
                                    self.Cu_GB2.lattice.matrix)


class Test_GrainBoundaryGenerator(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        cls.Cu_prim = Structure.from_file(os.path.join(test_dir, "Cu_mp-30_primitive.cif"))
        cls.GB_Cu_prim = GrainBoundaryGenerator(cls.Cu_prim)
        cls.Cu_conv = Structure.from_file(os.path.join(test_dir,
                                                       "Cu_mp-30_conventional_standard.cif"))
        cls.GB_Cu_conv = GrainBoundaryGenerator(cls.Cu_conv)
        cls.Be = Structure.from_file(os.path.join(test_dir,
                                                  "Be_mp-87_conventional_standard.cif"))
        cls.GB_Be = GrainBoundaryGenerator(cls.Be)
        cls.Pa = Structure.from_file(os.path.join(test_dir,
                                                  "Pa_mp-62_conventional_standard.cif"))
        cls.GB_Pa = GrainBoundaryGenerator(cls.Pa)
        cls.Br = Structure.from_file(os.path.join(test_dir,
                                                  "Br_mp-23154_conventional_standard.cif"))
        cls.GB_Br = GrainBoundaryGenerator(cls.Br)
        cls.Bi = Structure.from_file(os.path.join(test_dir,
                                                  "Bi_mp-23152_primitive.cif"))
        cls.GB_Bi = GrainBoundaryGenerator(cls.Bi)

    def test_gb_from_parameters(self):
        # from fcc primitive cell,axis[1,2,3],sigma 9.
        gb_cu_123_prim1 = self.GB_Cu_prim.gb_from_parameters([1, 2, 3],
                                                             123.74898859588858,
                                                             expand_times=2)
        lat_mat1 = gb_cu_123_prim1.lattice.matrix
        c_vec1 = np.cross(lat_mat1[0], lat_mat1[1]) / np.linalg.norm(np.cross(lat_mat1[0], lat_mat1[1]))
        c_len1 = np.dot(lat_mat1[2], c_vec1)
        vol_ratio = gb_cu_123_prim1.volume / self.Cu_prim.volume
        self.assertAlmostEqual(vol_ratio, 9 * 2 * 2, 8)
        # test expand_times and vacuum layer
        gb_cu_123_prim2 = self.GB_Cu_prim.gb_from_parameters([1, 2, 3],
                                                             123.74898859588858,
                                                             expand_times=4, vacuum_thickness=1.5)
        lat_mat2 = gb_cu_123_prim2.lattice.matrix
        c_vec2 = np.cross(lat_mat2[0], lat_mat2[1]) / np.linalg.norm(np.cross(lat_mat2[0], lat_mat2[1]))
        c_len2 = np.dot(lat_mat2[2], c_vec2)
        self.assertAlmostEqual((c_len2 - 1.5 * 2) / c_len1, 2)

        # test normal
        gb_cu_123_prim3 = self.GB_Cu_prim.gb_from_parameters([1, 2, 3],
                                                             123.74898859588858,
                                                             expand_times=2, normal=True)
        lat_mat3 = gb_cu_123_prim3.lattice.matrix
        c_vec3 = np.cross(lat_mat3[0], lat_mat3[1]) / np.linalg.norm(np.cross(lat_mat3[0], lat_mat3[1]))
        ab_len3 = np.linalg.norm(np.cross(lat_mat3[2], c_vec3))
        self.assertAlmostEqual(ab_len3, 0)

        # from fcc conventional cell,axis [1,2,3], siamg 9
        gb_cu_123_conv1 = self.GB_Cu_conv.gb_from_parameters([1, 2, 3],
                                                             123.74898859588858,
                                                             expand_times=4,
                                                             vacuum_thickness=1.5)
        lat_mat1 = gb_cu_123_conv1.lattice.matrix
        self.assertAlmostEqual(np.dot(lat_mat1[0], [1, 2, 3]), 0)
        self.assertAlmostEqual(np.dot(lat_mat1[1], [1, 2, 3]), 0)
        # test plane
        gb_cu_123_conv2 = self.GB_Cu_conv.gb_from_parameters([1, 2, 3],
                                                             123.74898859588858,
                                                             expand_times=2,
                                                             vacuum_thickness=1.5,
                                                             normal=False, plane=[1, 3, 1])
        lat_mat2 = gb_cu_123_conv2.lattice.matrix
        self.assertAlmostEqual(np.dot(lat_mat2[0], [1, 3, 1]), 0)
        self.assertAlmostEqual(np.dot(lat_mat2[1], [1, 3, 1]), 0)

        # from hex cell,axis [1,1,1], sigma 21
        gb_Be_111_1 = self.GB_Be.gb_from_parameters([1, 1, 1],
                                                    147.36310249644626,
                                                    ratio=[5, 2],
                                                    expand_times=4,
                                                    vacuum_thickness=1.5,
                                                    plane=[1, 2, 1])
        lat_priv = self.Be.lattice.matrix
        lat_mat1 = np.matmul(gb_Be_111_1.lattice.matrix, np.linalg.inv(lat_priv))
        self.assertAlmostEqual(np.dot(lat_mat1[0], [1, 2, 1]), 0)
        self.assertAlmostEqual(np.dot(lat_mat1[1], [1, 2, 1]), 0)
        # test volume associated with sigma value
        gb_Be_111_2 = self.GB_Be.gb_from_parameters([1, 1, 1], 147.36310249644626, ratio=[5, 2],
                                                    expand_times=4)
        vol_ratio = gb_Be_111_2.volume / self.Be.volume
        self.assertAlmostEqual(vol_ratio, 19 * 2 * 4)
        # test ratio = None, axis [0,0,1], sigma 7
        gb_Be_111_3 = self.GB_Be.gb_from_parameters([0, 0, 1], 21.786789298261812,
                                                    ratio=[5, 2],
                                                    expand_times=4)
        gb_Be_111_4 = self.GB_Be.gb_from_parameters([0, 0, 1], 21.786789298261812, ratio=None,
                                                    expand_times=4)
        self.assertTupleEqual(gb_Be_111_3.lattice.abc, gb_Be_111_4.lattice.abc)
        self.assertTupleEqual(gb_Be_111_3.lattice.angles, gb_Be_111_4.lattice.angles)
        gb_Be_111_5 = self.GB_Be.gb_from_parameters([3, 1, 0], 180.0, ratio=[5, 2],
                                                    expand_times=4)
        gb_Be_111_6 = self.GB_Be.gb_from_parameters([3, 1, 0], 180.0, ratio=None,
                                                    expand_times=4)
        self.assertTupleEqual(gb_Be_111_5.lattice.abc, gb_Be_111_6.lattice.abc)
        self.assertTupleEqual(gb_Be_111_5.lattice.angles, gb_Be_111_6.lattice.angles)

        # gb from tetragonal cell, axis[1,1,1], sigma 15
        gb_Pa_111_1 = self.GB_Pa.gb_from_parameters([1, 1, 1], 151.92751306414706, ratio=[2, 3],
                                                    expand_times=4, max_search=10)
        vol_ratio = gb_Pa_111_1.volume / self.Pa.volume
        self.assertAlmostEqual(vol_ratio, 17 * 2 * 4)

        # gb from orthorhombic cell, axis[1,1,1], sigma 83
        gb_Br_111_1 = self.GB_Br.gb_from_parameters([1, 1, 1], 131.5023374652235,
                                                    ratio=[21, 20, 5],
                                                    expand_times=4, max_search=10)
        vol_ratio = gb_Br_111_1.volume / self.Br.volume
        self.assertAlmostEqual(vol_ratio, 83 * 2 * 4)

        # gb from rhombohedra cell, axis[1,2,0], sigma 63
        gb_Bi_120_1 = self.GB_Bi.gb_from_parameters([1, 2, 0], 63.310675060280246,
                                                    ratio=[19, 5],
                                                    expand_times=4, max_search=5)
        vol_ratio = gb_Bi_120_1.volume / self.Bi.volume
        self.assertAlmostEqual(vol_ratio, 59 * 2 * 4)

    def test_get_ratio(self):
        # hexagnal
        Be_ratio = self.GB_Be.get_ratio(max_denominator=2)
        self.assertListEqual(Be_ratio, [5, 2])
        Be_ratio = self.GB_Be.get_ratio(max_denominator=5)
        self.assertListEqual(Be_ratio, [12, 5])
        # tetragonal
        Pa_ratio = self.GB_Pa.get_ratio(max_denominator=5)
        self.assertListEqual(Pa_ratio, [2, 3])
        # orthorombic
        Br_ratio = self.GB_Br.get_ratio(max_denominator=5)
        self.assertListEqual(Br_ratio, [21, 20, 5])
        # orthorombic
        Bi_ratio = self.GB_Bi.get_ratio(max_denominator=5)
        self.assertListEqual(Bi_ratio, [19, 5])

    def test_enum_sigma_cubic(self):
        true_100 = [5, 13, 17, 25, 29, 37, 41]
        true_110 = [3, 9, 11, 17, 19, 27, 33, 41, 43]
        true_111 = [3, 7, 13, 19, 21, 31, 37, 39, 43, 49]
        sigma_100 = list(GrainBoundaryGenerator.enum_sigma_cubic(50, [1, 0, 0]).keys())
        sigma_110 = list(GrainBoundaryGenerator.enum_sigma_cubic(50, [1, 1, 0]).keys())
        sigma_111 = list(GrainBoundaryGenerator.enum_sigma_cubic(50, [1, 1, 1]).keys())
        sigma_222 = list(GrainBoundaryGenerator.enum_sigma_cubic(50, [2, 2, 2]).keys())
        sigma_888 = list(GrainBoundaryGenerator.enum_sigma_cubic(50, [8, 8, 8]).keys())

        self.assertListEqual(sorted(true_100), sorted(sigma_100))
        self.assertListEqual(sorted(true_110), sorted(sigma_110))
        self.assertListEqual(sorted(true_111), sorted(sigma_111))
        self.assertListEqual(sorted(true_111), sorted(sigma_222))
        self.assertListEqual(sorted(true_111), sorted(sigma_888))

    def test_enum_sigma_hex(self):
        true_100 = [17, 18, 22, 27, 38, 41]
        true_001 = [7, 13, 19, 31, 37, 43, 49]
        true_210 = [10, 11, 14, 25, 35, 49]
        sigma_100 = list(GrainBoundaryGenerator.enum_sigma_hex(50, [1, 0, 0], [8, 3]).keys())
        sigma_001 = list(GrainBoundaryGenerator.enum_sigma_hex(50, [0, 0, 1], [8, 3]).keys())
        sigma_210 = list(GrainBoundaryGenerator.enum_sigma_hex(50, [2, 1, 0], [8, 3]).keys())
        sigma_420 = list(GrainBoundaryGenerator.enum_sigma_hex(50, [4, 2, 0], [8, 3]).keys())
        sigma_840 = list(GrainBoundaryGenerator.enum_sigma_hex(50, [8, 4, 0], [8, 3]).keys())

        self.assertListEqual(sorted(true_100), sorted(sigma_100))
        self.assertListEqual(sorted(true_001), sorted(sigma_001))
        self.assertListEqual(sorted(true_210), sorted(sigma_210))
        self.assertListEqual(sorted(true_210), sorted(sigma_420))
        self.assertListEqual(sorted(true_210), sorted(sigma_840))

    def test_enum_sigma_tet(self):
        true_100 = [5, 37, 41, 13, 3, 15, 39, 25, 17, 29]
        true_331 = [9, 3, 21, 39, 7, 31, 43, 13, 19, 37, 49]
        sigma_100 = list(GrainBoundaryGenerator.enum_sigma_tet(50, [1, 0, 0], [9, 1]).keys())
        sigma_331 = list(GrainBoundaryGenerator.enum_sigma_tet(50, [3, 3, 1], [9, 1]).keys())

        self.assertListEqual(sorted(true_100), sorted(sigma_100))
        self.assertListEqual(sorted(true_331), sorted(sigma_331))

    def test_enum_sigma_ort(self):
        true_100 = [41, 37, 39, 5, 15, 17, 13, 3, 25, 29]
        sigma_100 = list(GrainBoundaryGenerator.enum_sigma_ort(50, [1, 0, 0], [270, 30, 29]).keys())

        self.assertListEqual(sorted(true_100), sorted(sigma_100))

    def test_enum_sigma_rho(self):
        true_100 = [7, 11, 43, 13, 41, 19, 47, 31]
        sigma_100 = list(GrainBoundaryGenerator.enum_sigma_rho(50, [1, 0, 0], [15, 4]).keys())

        self.assertListEqual(sorted(true_100), sorted(sigma_100))

    def test_enum_possible_plane_cubic(self):
        all_plane = GrainBoundaryGenerator.enum_possible_plane_cubic(4, [1, 1, 1], 60)
        self.assertEqual(len(all_plane['Twist']), 1)
        self.assertEqual(len(all_plane['Symmetric tilt']), 6)
        self.assertEqual(len(all_plane['Normal tilt']), 12)

    def test_get_trans_mat(self):
        mat1, mat2 = GrainBoundaryGenerator.get_trans_mat([1, 1, 1], 95.55344419565849,
                                                          lat_type='o', ratio=[10, 20, 21],
                                                          surface=[21, 20, 10], normal=True)
        self.assertAlmostEqual(np.dot(mat1[0], [21, 20, 10]), 0)
        self.assertAlmostEqual(np.dot(mat1[1], [21, 20, 10]), 0)
        self.assertAlmostEqual(np.linalg.det(mat1), np.linalg.det(mat2))
        ab_len1 = np.linalg.norm(np.cross(mat1[2], [1, 1, 1]))
        self.assertAlmostEqual(ab_len1, 0)

    def test_get_rotation_angle_from_sigma(self):
        true_angle = [12.680383491819821, 167.3196165081802]
        angle = GrainBoundaryGenerator.get_rotation_angle_from_sigma(41, [1, 0, 0], lat_type='o', ratio=[270, 30, 29])
        self.assertArrayAlmostEqual(true_angle, angle)
        close_angle = [36.86989764584403, 143.13010235415598]
        angle = GrainBoundaryGenerator.get_rotation_angle_from_sigma(6, [1, 0, 0], lat_type='o', ratio=[270, 30, 29])
        self.assertArrayAlmostEqual(close_angle, angle)
