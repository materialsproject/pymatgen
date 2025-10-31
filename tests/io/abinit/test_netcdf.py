from __future__ import annotations

import os
import sys
import tarfile

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from pymatgen.core.structure import Structure
from pymatgen.io.abinit import EtsfReader
from pymatgen.io.abinit.netcdf import AbinitHeader
from pymatgen.util.testing import TEST_FILES_DIR, MatSciTest

try:
    import netCDF4
except ImportError:
    netCDF4 = None

TEST_DIR = f"{TEST_FILES_DIR}/io/abinit"


class TestEtsfReader(MatSciTest):
    def setup_method(self):
        formulas = ["Si2"]
        self.GSR_paths = dct = {}
        for formula in formulas:
            dct[formula] = f"{TEST_DIR}/{formula}_GSR.nc"

    @pytest.mark.skipif(netCDF4 is None or True, reason="Requires Netcdf4")
    def test_read_si2(self):
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

        with EtsfReader(path) as data:
            assert data.ngroups == 1

            # Test dimensions.
            for dim_name, int_ref in ref_dims.items():
                value = data.read_dimvalue(dim_name)
                assert_array_equal(value, int_ref)

            # Test int variables
            for varname, int_ref in ref_int_values.items():
                value = data.read_value(varname)
                assert_array_equal(value, int_ref)

            # Test float variables
            for varname, float_ref in ref_float_values.items():
                value = data.read_value(varname)
                assert_allclose(value, float_ref)

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
            assert "magmom" not in structure.site_properties

            # Read ixc.
            # TODO: Upgrade GSR file.
            # xc = data.read_abinit_xcfunc()
            # assert xc == "LDA"

    @pytest.mark.skipif(netCDF4 is None, reason="Requires NetCDF4")
    @pytest.mark.xfail(
        sys.platform.startswith("linux") and os.getenv("CI") and netCDF4.__version__ >= "1.6.5",
        reason="Fails with netCDF4 >= 1.6.5 on Ubuntu CI",
    )
    @pytest.mark.parametrize(
        ("tarname", "expected"),
        [
            ("Fe_magmoms_collinear_GSR.tar.xz", [-0.5069359730980665]),
            ("Fe_magmoms_noncollinear_GSR.tar.xz", [[0.357939487, 0.357939487, 0]]),
        ],
    )
    def test_read_fe(self, tarname, expected):
        with tarfile.open(os.path.join(TEST_DIR, tarname), mode="r:xz") as t:
            # TODO: remove attr check after only 3.12+
            if hasattr(tarfile, "data_filter"):
                t.extractall(".", filter="data")
            else:
                t.extractall(".")  # noqa: S202

            path = os.path.basename(tarname).replace(".tar.xz", ".nc")

            with EtsfReader(path) as data:
                structure = data.read_structure()
                assert structure.site_properties["magmom"] == expected


class TestAbinitHeader(MatSciTest):
    def test_api(self):
        head = AbinitHeader(foo=1, bar=2)
        assert head.foo == 1
        assert str(head)
        assert head.to_str(verbose=2, title="title")
        # PLEASE DO NOT REMOVE THIS LINE AS THIS API HAS BEEN AROUND FOR SEVERAL YEARS,
        with pytest.warns(FutureWarning, match="to_string is deprecated"):
            assert head.to_string(verbose=2, title="title")
