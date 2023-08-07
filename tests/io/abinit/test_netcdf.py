from __future__ import annotations

import os
import unittest

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from pymatgen.core.structure import Structure
from pymatgen.io.abinit import ETSF_Reader
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

try:
    import netCDF4
except ImportError:
    netCDF4 = None

_test_dir = f"{TEST_FILES_DIR}/abinit"


def ref_file(filename):
    return os.path.join(_test_dir, filename)


class ETSF_Reader_TestCase(PymatgenTest):
    def setUp(self):
        formulas = ["Si2"]
        self.GSR_paths = dct = {}
        for formula in formulas:
            dct[formula] = ref_file(formula + "_GSR.nc")

    @unittest.skipIf(netCDF4 is None, "Requires Netcdf4")
    def test_read_Si2(self):
        path = self.GSR_paths["Si2"]

        ref_dims = {"number_of_spins": 1}

        ref_int_values = {
            "space_group": 227,
            "number_of_states": np.reshape([15, 15], (1, 2)),
        }

        ref_float_values = {
            "etotal": -8.85911566912484,
            "primitive_vectors": np.reshape([0, 5.125, 5.125, 5.125, 0, 5.125, 5.125, 5.125, 0], (3, 3)),
        }

        with ETSF_Reader(path) as data:
            assert data.ngroups == 1

            # Test dimensions.
            for dimname, int_ref in ref_dims.items():
                value = data.read_dimvalue(dimname)
                assert_array_equal(value, int_ref)

            # Test int variables
            for varname, int_ref in ref_int_values.items():
                value = data.read_value(varname)
                assert_array_equal(value, int_ref)

            # Test float variables
            for varname, float_ref in ref_float_values.items():
                value = data.read_value(varname)
                assert np.allclose(value, float_ref)
            # assert 0

            # Reading non-existent variables or dims should raise
            # a subclass of NetcdReaderError
            with pytest.raises(data.Error):
                data.read_value("foobar")

            with pytest.raises(data.Error):
                data.read_dimvalue("foobar")

            # Unless default is given
            assert data.read_value("foobar", default=None) is None

            # Initialize pymatgen structure from GSR.
            structure = data.read_structure()
            assert isinstance(structure, Structure)

            # Read ixc.
            # TODO: Upgrade GSR file.
            # xc = data.read_abinit_xcfunc()
            # assert xc == "LDA"


class TestAbinitHeader(PymatgenTest):
    def test_api(self):
        from pymatgen.io.abinit.netcdf import AbinitHeader

        head = AbinitHeader(foo=1, bar=2)
        assert head.foo == 1
        assert str(head)
        assert head.to_str(verbose=2, title="title")
