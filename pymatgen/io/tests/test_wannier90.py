# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Tests for pymatgen.io.wannier90
"""

import numpy as np
from monty.tempfile import ScratchDir
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.wannier90 import Unk


class UnkTest(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        self.data_std = np.random.rand(10, 5, 5, 5).astype(np.complex128)
        self.unk_std = Unk(1, self.data_std)
        self.data_ncl = np.random.rand(10, 2, 5, 5, 5).astype(np.complex128)
        self.unk_ncl = Unk(1, self.data_ncl)

    def test_init(self):
        # standard unk file
        self.assertEqual(self.unk_std.ik, 1)
        self.assertEqual(self.unk_std.nbnd, 10)
        self.assertEqual(self.unk_std.ng[0], 5)
        self.assertEqual(self.unk_std.ng[1], 5)
        self.assertEqual(self.unk_std.ng[2], 5)
        self.assertTrue(np.allclose(self.unk_std.data.shape, (10, 5, 5, 5)))
        self.assertTrue(np.allclose(self.unk_std.data, self.data_std))
        self.assertFalse(self.unk_std.is_noncollinear)

        # too small data
        with self.assertRaises(ValueError):
            data_bad_shape = \
                np.random.rand(2, 2, 2).astype(np.complex128)
            Unk(1, data_bad_shape)

        # too big data
        with self.assertRaises(ValueError):
            data_bad_shape = \
                np.random.rand(2, 2, 2, 2, 2, 2).astype(np.complex128)
            Unk(1, data_bad_shape)

        # noncollinear unk file
        self.assertEqual(self.unk_ncl.ik, 1)
        self.assertEqual(self.unk_ncl.nbnd, 10)
        self.assertEqual(self.unk_ncl.ng[0], 5)
        self.assertEqual(self.unk_ncl.ng[1], 5)
        self.assertEqual(self.unk_ncl.ng[2], 5)
        self.assertTrue(np.allclose(self.unk_ncl.data.shape, (10, 2, 5, 5, 5)))
        self.assertTrue(np.allclose(self.unk_ncl.data, self.data_ncl))
        self.assertTrue(self.unk_ncl.is_noncollinear)

        # too big data
        with self.assertRaises(ValueError):
            data_bad_ncl = \
                np.random.rand(2, 3, 2, 2, 2).astype(np.complex128)
            Unk(1, data_bad_ncl)

    def test_from_file(self):
        unk = Unk.from_file(self.TEST_FILES_DIR / 'UNK.std')
        self.assertEqual(unk.ik, 1)
        self.assertEqual(unk.nbnd, 5)
        self.assertEqual(unk.ng[0], 6)
        self.assertEqual(unk.ng[1], 6)
        self.assertEqual(unk.ng[2], 8)
        self.assertFalse(unk.is_noncollinear)
        self.assertTrue(np.allclose(unk.data.shape, (5, 6, 6, 8)))

        unk = Unk.from_file(self.TEST_FILES_DIR / 'UNK.ncl')
        self.assertEqual(unk.ik, 1)
        self.assertEqual(unk.nbnd, 5)
        self.assertEqual(unk.ng[0], 6)
        self.assertEqual(unk.ng[1], 6)
        self.assertEqual(unk.ng[2], 8)
        self.assertTrue(unk.is_noncollinear)
        self.assertTrue(np.allclose(unk.data.shape, (5, 2, 6, 6, 8)))
        self.assertNotEqual(unk.data[0, 0, 0, 0, 0].real, 0.)
        self.assertAlmostEqual(unk.data[0, 1, 0, 0, 0].real, 0.)

    def test_write_file(self):
        with ScratchDir('.'):
            self.unk_std.write_file('UNK00001.1')
            temp_unk = Unk.from_file('UNK00001.1')
            self.assertEqual(self.unk_std, temp_unk)

        with ScratchDir('.'):
            self.unk_ncl.write_file('UNK00001.NC')
            temp_unk = Unk.from_file('UNK00001.NC')
            self.assertEqual(self.unk_ncl, temp_unk)

    def test_read_write(self):
        unk0 = Unk.from_file(self.TEST_FILES_DIR / 'UNK.std')
        with ScratchDir('.'):
            unk0.write_file('UNK00001.1')
            unk1 = Unk.from_file('UNK00001.1')
            self.assertEqual(unk0, unk1)

        unk0 = Unk.from_file(self.TEST_FILES_DIR / 'UNK.ncl')
        with ScratchDir('.'):
            unk0.write_file('UNK00001.NC')
            unk1 = Unk.from_file('UNK00001.NC')
            self.assertEqual(unk0, unk1)

    def test_repr(self):
        self.assertNotEqual(repr(self.unk_std), '')
        self.assertNotEqual(repr(self.unk_ncl), '')

    def test_eq(self):
        # not implemented
        self.assertFalse(self.unk_std == 'poop')

        # ng
        tmp_unk = Unk(1, np.random.rand(10, 5, 5, 4))
        self.assertFalse(self.unk_std == tmp_unk)

        # ik
        tmp_unk = Unk(2, self.data_std)
        self.assertFalse(self.unk_std == tmp_unk)

        # noncol
        self.assertFalse(self.unk_std == self.unk_ncl)

        # nbnd
        tmp_unk = Unk(1, np.random.rand(9, 5, 5, 5))
        self.assertFalse(self.unk_std == tmp_unk)

        # data
        tmp_unk = Unk(1, np.random.rand(10, 5, 5, 5))
        self.assertFalse(self.unk_std == tmp_unk)
        tmp_unk = Unk(1, np.random.rand(10, 2, 5, 5, 5))
        self.assertFalse(self.unk_ncl == tmp_unk)

        # same
        self.assertTrue(self.unk_std == self.unk_std)
        self.assertTrue(self.unk_std == Unk(1, self.data_std))
        self.assertTrue(self.unk_ncl == self.unk_ncl)
        self.assertTrue(self.unk_ncl == Unk(1, self.data_ncl))
