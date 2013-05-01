from __future__ import division, print_function

import os
import unittest
import collections
import numpy as np

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinitio import GSR_Reader

import numpy.testing.utils as nptu

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')

def get_netcdf_version():
    """
    Returns the Netcdf version we can read with the available libraries.
    0 if not library is found.
    """
    try:
        import netCDF4
        return 4
    except ImportError:
        try:
            from scipy.io import netcdf
            return 3
        except ImportError:
            return 0

NC_VERSION = get_netcdf_version()

def filepath(basename):
    return os.path.join(test_dir, basename)


class GSR_Reader_TestCase(PymatgenTest):

    def setUp(self):
        formulas = ["Si2",]
        self.GSR_paths = d = {}
        for formula in formulas:
            d[formula] = filepath(formula + "_GSR.nc")

    @unittest.skipIf(NC_VERSION!=4, "Requires Netcdf4")
    def test_read_Si2(self):
        path = self.GSR_paths["Si2"]

        ref_int_values = {
            "space_group" : 227,
            "number_of_states" : np.reshape([15, 15], (1,2)),
        }

        ref_float_values = {
            "etotal" : -8.85911566912484,
            "primitive_vectors" : np.reshape([0, 5.125, 5.125, 5.125, 0, 5.125,
                                              5.125, 5.125, 0], (3,3)),
        }

        with GSR_Reader(path) as data:

            self.assertEqual(data.ngroups, 1)

            print(data.get_varnames())

            # Test int variables
            for (varname, int_ref) in ref_int_values.items():
                value = data.get_value(varname)

                self.assert_equal(value, int_ref)

            # Test float variables
            for (varname, float_ref) in ref_float_values.items():
                value = data.get_value(varname)

                self.assert_almost_equal(value, float_ref)

            data.print_tree()
            for group in data.walk_tree():
                print("group: " + str(group))

            # Init structure from GSR.
            # TODO change name in the ETSF file since site does not accept cartesian_forces
            #site_properties = ["forces",]
            site_properties = []

            structure = data.get_structure(site_properties=site_properties)


if __name__ == "__main__":
    unittest.main()
