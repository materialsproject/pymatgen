# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

import os
import unittest
import numpy as np

from pymatgen import Structure
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinit import ETSF_Reader

try:
    import netCDF4
except ImportError:
    netCDF4 = None


_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", 
                         'test_files', "abinit")

def ref_file(filename):
    return os.path.join(_test_dir, filename)


class ETSF_Reader_TestCase(PymatgenTest):

    def setUp(self):
        formulas = ["Si2",]
        self.GSR_paths = d = {}
        for formula in formulas:
            d[formula] = ref_file(formula + "_GSR.nc")

    @unittest.skipIf(netCDF4 is None, "Requires Netcdf4")
    def test_read_Si2(self):
        path = self.GSR_paths["Si2"]

        ref_dims = {
            "number_of_spins": 1
        }

        ref_int_values = {
            "space_group": 227,
            "number_of_states": np.reshape([15, 15], (1,2)),
        }

        ref_float_values = {
            "etotal": -8.85911566912484,
            "primitive_vectors": np.reshape([0, 5.125, 5.125, 5.125, 0, 5.125,
                                             5.125, 5.125, 0], (3,3)),
        }

        with ETSF_Reader(path) as data:

            self.assertEqual(data.ngroups, 1)

            print(data.read_varnames())

            # Test dimensions.
            for (dimname, int_ref) in ref_dims.items():
                value = data.read_dimvalue(dimname)
                self.assert_equal(value, int_ref)

            # Test int variables
            for (varname, int_ref) in ref_int_values.items():
                value = data.read_value(varname)
                print(varname, value)
                self.assert_equal(value, int_ref)

            # Test float variables
            for (varname, float_ref) in ref_float_values.items():
                value = data.read_value(varname)
                print(varname, value)
                self.assert_almost_equal(value, float_ref)
            #assert 0

            # Reading non-existent variables or dims should raise
            # a subclass of NetcdReaderError
            with self.assertRaises(data.Error):
                data.read_value("foobar")

            with self.assertRaises(data.Error):
                data.read_dimvalue("foobar")

            # Unless default is given
            assert data.read_value("foobar", default=None) is None

            data.print_tree()
            for group in data.walk_tree():
                print("group: " + str(group))

            # Initialize pymatgen structure from GSR.
            structure = data.read_structure()
            self.assertTrue(isinstance(structure, Structure))


if __name__ == "__main__":
    unittest.main()
