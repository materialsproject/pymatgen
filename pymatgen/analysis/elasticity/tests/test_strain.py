from __future__ import absolute_import

import unittest2 as unittest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.analysis.elasticity.strain import Strain, Deformation, DeformedStructureSet,\
                                    IndependentStrain
from pymatgen.util.testing import PymatgenTest
import warnings


class DeformationTest(PymatgenTest):
    def setUp(self):
        self.norm_defo = Deformation.from_index_amount((0, 0), 0.02)
        self.ind_defo = Deformation.from_index_amount((0, 1), 0.02)
        self.non_ind_defo = Deformation([[1.0, 0.02, 0.02],
                                         [0.0, 1.0, 0.0],
                                         [0.0, 0.0, 1.0]])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        self.structure = Structure(lattice, ["Si", "Si"], [[0, 0, 0],
                                                           [0.75, 0.5, 0.75]])

    def test_properties(self):
        # green_lagrange_strain
        self.assertArrayAlmostEqual(self.ind_defo.green_lagrange_strain,
                                    [[0., 0.01, 0.],
                                     [0.01, 0.0002, 0.],
                                     [0., 0., 0.]])
        self.assertArrayAlmostEqual(self.non_ind_defo.green_lagrange_strain,
                                    [[0., 0.01, 0.01],
                                     [0.01, 0.0002, 0.0002],
                                     [0.01, 0.0002, 0.0002]])

    def test_check_independent(self):
        self.assertRaises(ValueError, self.non_ind_defo.check_independent)
        self.assertEqual(self.ind_defo.check_independent(), (0, 1))

    def test_apply_to_structure(self):
        strained_norm = self.norm_defo.apply_to_structure(self.structure)
        strained_ind = self.ind_defo.apply_to_structure(self.structure)
        strained_non = self.non_ind_defo.apply_to_structure(self.structure)
        # Check lattices
        self.assertArrayAlmostEqual(strained_norm.lattice.matrix,
                                    [[3.9170018886, 0, 0],
                                     [1.958500946136, 3.32571019, 0],
                                     [0, -2.21713849, 3.13550906]])
        self.assertArrayAlmostEqual(strained_ind.lattice.matrix,
                                    [[3.84019793, 0.07680396, 0],
                                     [1.92009897, 3.36411217, 0],
                                     [0, -2.21713849, 3.13550906]])
        self.assertArrayAlmostEqual(strained_non.lattice.matrix,
                                    [[3.84019793, 0.07680396, 0.07680396],
                                     [1.92009897, 3.36411217, 0.0384019794],
                                     [0, -2.21713849, 3.13550906]])
        # Check coordinates
        self.assertArrayAlmostEqual(strained_norm.sites[1].coords,
                                    [3.91700189, 1.224e-06, 2.3516318])
        self.assertArrayAlmostEqual(strained_ind.sites[1].coords,
                                    [3.84019793, 0.07680518, 2.3516318])
        self.assertArrayAlmostEqual(strained_non.sites[1].coords,
                                    [3.84019793, 0.07680518, 2.42843575])


class StrainTest(PymatgenTest):
    def setUp(self):
        self.norm_str = Strain.from_deformation([[1.02, 0, 0],
                                                 [0, 1, 0],
                                                 [0, 0, 1]])
        self.ind_str = Strain.from_deformation([[1, 0.02, 0],
                                                [0, 1, 0],
                                                [0, 0, 1]])

        self.non_ind_str = Strain.from_deformation([[1, 0.02, 0.02],
                                                    [0, 1, 0],
                                                    [0, 0, 1]])

        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            self.no_dfm = Strain([[0., 0.01, 0.],
                                  [0.01, 0.0002, 0.],
                                  [0., 0., 0.]])

    def test_new(self):
        with warnings.catch_warnings(record=True) as w:
            Strain([[0., 0.01, 0.],
                    [0.01, 0.0002, 0.],
                    [0., 0., 0.]])
            self.assertEqual(len(w), 1)
        self.assertRaises(ValueError, Strain, [[0.1, 0.1, 0],
                                               [0, 0, 0],
                                               [0, 0, 0]])

    def test_from_deformation(self):
        self.assertArrayAlmostEqual(self.norm_str,
                                    [[0.0202, 0, 0],
                                     [0, 0, 0],
                                     [0, 0, 0]])
        self.assertArrayAlmostEqual(self.ind_str,
                                    [[0., 0.01, 0.],
                                     [0.01, 0.0002, 0.],
                                     [0., 0., 0.]])
        self.assertArrayAlmostEqual(self.non_ind_str,
                                    [[0., 0.01, 0.01],
                                     [0.01, 0.0002, 0.0002],
                                     [0.01, 0.0002, 0.0002]])

    def test_properties(self):
        # deformation matrix
        self.assertArrayAlmostEqual(self.ind_str.deformation_matrix,
                                    [[1, 0.02, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]])
        self.assertArrayAlmostEqual(self.no_dfm.deformation_matrix,
                                    [[0.99995,0.0099995, 0],
                                     [0.0099995,1.00015, 0],
                                     [0, 0, 1]])


        # independent deformation
        self.assertArrayEqual(self.ind_str.independent_deformation, (0, 1))
        with self.assertRaises(ValueError):
            self.no_dfm.independent_deformation

        # voigt
        self.assertArrayAlmostEqual(self.non_ind_str.voigt,
                                    [0, 0.0002, 0.0002, 0.0004, 0.02, 0.02])
    def test_convert_strain_to_deformation(self):
        defo = Deformation.from_index_amount((1,2), 0.01)
        pass

class DeformedStructureSetTest(PymatgenTest):
    def setUp(self):
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        self.structure = Structure(lattice, ["Si", "Si"], [[0, 0, 0],
                                                           [0.75, 0.5, 0.75]])
        self.default_dss = DeformedStructureSet(self.structure)

    def test_init(self):
        with self.assertRaises(ValueError):
            DeformedStructureSet(self.structure, num_norm=5)
        with self.assertRaises(ValueError):
            DeformedStructureSet(self.structure, num_shear=5)
        self.assertEqual(self.structure, self.default_dss.undeformed_structure)

    def test_as_strain_dict(self):
        strain_dict = self.default_dss.as_strain_dict()
        for i, def_struct in enumerate(self.default_dss):
            test_strain = IndependentStrain(self.default_dss.deformations[i])
            strain_keys = [strain for strain in list(strain_dict.keys())
                           if (strain == test_strain).all()]
            self.assertEqual(len(strain_keys), 1)
            self.assertEqual(self.default_dss.def_structs[i],
                             strain_dict[strain_keys[0]])


class IndependentStrainTest(PymatgenTest):
    def setUp(self):
        self.ind_strain = IndependentStrain([[1, 0.1, 0],
                                             [0, 1, 0],
                                             [0, 0, 1]])

    def test_new(self):
        with self.assertRaises(ValueError):
            IndependentStrain([[0.1, 0.1, 0],
                               [0, 0, 0],
                               [0, 0, 0]])

    def test_properties(self):
        self.assertEqual(self.ind_strain.i, 0)
        self.assertEqual(self.ind_strain.j, 1)

if __name__ == '__main__':
    unittest.main()
