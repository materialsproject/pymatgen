import unittest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.tensors import Tensor
from pymatgen.analysis.elasticity.strain import Strain, Deformation, \
    convert_strain_to_deformation, DeformedStructureSet
from pymatgen.util.testing import PymatgenTest
import numpy as np
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

    def test_independence(self):
        self.assertFalse(self.non_ind_defo.is_independent())
        self.assertEqual(self.ind_defo.get_perturbed_indices()[0], (0, 1))

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
                                    [[3.84019793, 0, 0],
                                     [1.9866132, 3.32571019, 0],
                                     [-0.04434277, -2.21713849, 3.13550906]])
        self.assertArrayAlmostEqual(strained_non.lattice.matrix,
                                    [[3.84019793, 0, 0],
                                     [1.9866132, 3.3257102, 0],
                                     [0.0183674, -2.21713849, 3.13550906]])
        # Check coordinates
        self.assertArrayAlmostEqual(strained_norm.sites[1].coords,
                                    [3.91700189, 1.224e-06, 2.3516318])
        self.assertArrayAlmostEqual(strained_ind.sites[1].coords,
                                    [3.84019793, 1.224e-6, 2.3516318])
        self.assertArrayAlmostEqual(strained_non.sites[1].coords,
                                    [3.8872306, 1.224e-6, 2.3516318])

        # Check convention for applying transformation
        for vec, defo_vec in zip(self.structure.lattice.matrix,
                                 strained_non.lattice.matrix):
            new_vec = np.dot(self.non_ind_defo, np.transpose(vec))
            self.assertArrayAlmostEqual(new_vec, defo_vec)
        for coord, defo_coord in zip(self.structure.cart_coords, strained_non.cart_coords):
            new_coord = np.dot(self.non_ind_defo, np.transpose(coord))
            self.assertArrayAlmostEqual(new_coord, defo_coord)


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
        test_strain = Strain([[0., 0.01, 0.],
                              [0.01, 0.0002, 0.],
                              [0., 0., 0.]])
        self.assertArrayAlmostEqual(
            test_strain,
            test_strain.get_deformation_matrix().green_lagrange_strain)
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

    def test_from_index_amount(self):
        # From voigt index
        test = Strain.from_index_amount(2, 0.01)
        should_be = np.zeros((3, 3))
        should_be[2, 2] = 0.01
        self.assertArrayAlmostEqual(test, should_be)
        # from full-tensor index
        test = Strain.from_index_amount((1, 2), 0.01)
        should_be = np.zeros((3, 3))
        should_be[1, 2] = should_be[2, 1] = 0.01
        self.assertArrayAlmostEqual(test, should_be)

    def test_properties(self):
        # deformation matrix
        self.assertArrayAlmostEqual(self.ind_str.get_deformation_matrix(),
                                    [[1, 0.02, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]])
        symm_dfm = Strain(self.no_dfm).get_deformation_matrix(shape="symmetric")
        self.assertArrayAlmostEqual(symm_dfm, [[0.99995, 0.0099995, 0],
                                               [0.0099995, 1.00015, 0],
                                               [0, 0, 1]])
        self.assertArrayAlmostEqual(self.no_dfm.get_deformation_matrix(),
                                    [[1, 0.02, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]])

        # voigt
        self.assertArrayAlmostEqual(self.non_ind_str.voigt,
                                    [0, 0.0002, 0.0002, 0.0004, 0.02, 0.02])

    def test_convert_strain_to_deformation(self):
        strain = Tensor(np.random.random((3, 3))).symmetrized
        while not (np.linalg.eigvals(strain) > 0).all():
            strain = Tensor(np.random.random((3, 3))).symmetrized
        upper = convert_strain_to_deformation(strain, shape="upper")
        symm = convert_strain_to_deformation(strain, shape="symmetric")
        self.assertArrayAlmostEqual(np.triu(upper), upper)
        self.assertTrue(Tensor(symm).is_symmetric())
        for defo in upper, symm:
            self.assertArrayAlmostEqual(defo.green_lagrange_strain, strain)


class DeformedStructureSetTest(PymatgenTest):
    def setUp(self):
        self.structure = self.get_structure("Sn")
        self.default_dss = DeformedStructureSet(self.structure)

    def test_init(self):
        self.assertEqual(self.structure, self.default_dss.undeformed_structure)
        # Test symmetry
        dss_symm = DeformedStructureSet(self.structure, symmetry=True)
        # Should be 4 strains for normal, 2 for shear (since +/- shear
        # are symmetrically equivalent)
        self.assertEqual(len(dss_symm), 6)


if __name__ == '__main__':
    unittest.main()
