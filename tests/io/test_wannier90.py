"""Tests for pymatgen.io.wannier90."""

from __future__ import annotations

import numpy as np
import pytest
from pytest import approx

from pymatgen.io.wannier90 import Unk
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestUnk(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        self.data_std = np.random.rand(10, 5, 5, 5).astype(np.complex128)
        self.unk_std = Unk(1, self.data_std)
        self.data_ncl = np.random.rand(10, 2, 5, 5, 5).astype(np.complex128)
        self.unk_ncl = Unk(1, self.data_ncl)

    def test_init(self):
        # standard unk file
        assert self.unk_std.ik == 1
        assert self.unk_std.nbnd == 10
        assert self.unk_std.ng[0] == 5
        assert self.unk_std.ng[1] == 5
        assert self.unk_std.ng[2] == 5
        assert np.allclose(self.unk_std.data.shape, (10, 5, 5, 5))
        assert np.allclose(self.unk_std.data, self.data_std)
        assert not self.unk_std.is_noncollinear

        # too small data
        data_bad_shape = np.random.rand(2, 2, 2).astype(np.complex128)
        with pytest.raises(
            ValueError,
            match=r"invalid data shape, must be \(nbnd, ngx, ngy, ngz\) or \(nbnd, 2, ngx, ngy, ngz\) "
            r"for noncollinear data, given \(2, 2, 2\)",
        ):
            Unk(1, data_bad_shape)

        # too big data
        data_bad_shape = np.random.rand(2, 2, 2, 2, 2, 2).astype(np.complex128)
        with pytest.raises(
            ValueError,
            match=r"invalid data shape, must be \(nbnd, ngx, ngy, ngz\) or \(nbnd, 2, ngx, ngy, ngz\) for noncollinear",
        ):
            Unk(1, data_bad_shape)

        # noncollinear unk file
        assert self.unk_ncl.ik == 1
        assert self.unk_ncl.nbnd == 10
        assert self.unk_ncl.ng[0] == 5
        assert self.unk_ncl.ng[1] == 5
        assert self.unk_ncl.ng[2] == 5
        assert np.allclose(self.unk_ncl.data.shape, (10, 2, 5, 5, 5))
        assert np.allclose(self.unk_ncl.data, self.data_ncl)
        assert self.unk_ncl.is_noncollinear

        # too big data
        data_bad_ncl = np.random.rand(2, 3, 2, 2, 2).astype(np.complex128)
        with pytest.raises(
            ValueError,
            match=r"invalid noncollinear data, shape should be \(nbnd, 2, ngx, ngy, ngz\), given \(2, 3, 2, 2, 2\)",
        ):
            Unk(1, data_bad_ncl)

    def test_from_file(self):
        unk = Unk.from_file(f"{TEST_FILES_DIR}/UNK.std")
        assert unk.ik == 1
        assert unk.nbnd == 5
        assert unk.ng[0] == 6
        assert unk.ng[1] == 6
        assert unk.ng[2] == 8
        assert not unk.is_noncollinear
        assert np.allclose(unk.data.shape, (5, 6, 6, 8))

        unk = Unk.from_file(f"{TEST_FILES_DIR}/UNK.ncl")
        assert unk.ik == 1
        assert unk.nbnd == 5
        assert unk.ng[0] == 6
        assert unk.ng[1] == 6
        assert unk.ng[2] == 8
        assert unk.is_noncollinear
        assert np.allclose(unk.data.shape, (5, 2, 6, 6, 8))
        assert unk.data[0, 0, 0, 0, 0].real != 0.0
        assert unk.data[0, 1, 0, 0, 0].real == approx(0.0)

    def test_write_file(self):
        self.unk_std.write_file("UNK00001.1")
        temp_unk = Unk.from_file("UNK00001.1")
        assert self.unk_std == temp_unk

        self.unk_ncl.write_file("UNK00001.NC")
        temp_unk = Unk.from_file("UNK00001.NC")
        assert self.unk_ncl == temp_unk

    def test_read_write(self):
        unk0 = Unk.from_file(f"{TEST_FILES_DIR}/UNK.std")
        unk0.write_file("UNK00001.1")
        unk1 = Unk.from_file("UNK00001.1")
        assert unk0 == unk1

        unk0 = Unk.from_file(f"{TEST_FILES_DIR}/UNK.ncl")
        unk0.write_file("UNK00001.NC")
        unk1 = Unk.from_file("UNK00001.NC")
        assert unk0 == unk1

    def test_repr(self):
        assert repr(self.unk_std) != ""
        assert repr(self.unk_ncl) != ""

    def test_eq(self):
        # not implemented
        assert self.unk_std != "poop"

        # ng
        tmp_unk = Unk(1, np.random.rand(10, 5, 5, 4))
        assert self.unk_std != tmp_unk

        # ik
        tmp_unk = Unk(2, self.data_std)
        assert self.unk_std != tmp_unk

        # noncol
        assert self.unk_std != self.unk_ncl

        # nbnd
        tmp_unk = Unk(1, np.random.rand(9, 5, 5, 5))
        assert self.unk_std != tmp_unk

        # data
        tmp_unk = Unk(1, np.random.rand(10, 5, 5, 5))
        assert self.unk_std != tmp_unk
        tmp_unk = Unk(1, np.random.rand(10, 2, 5, 5, 5))
        assert self.unk_ncl != tmp_unk

        # same
        assert self.unk_std == self.unk_std
        assert self.unk_std == Unk(1, self.data_std)
        assert self.unk_ncl == self.unk_ncl
        assert self.unk_ncl == Unk(1, self.data_ncl)
