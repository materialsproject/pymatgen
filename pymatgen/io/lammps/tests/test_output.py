# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

import numpy as np

from pymatgen.io.lammps.output import LammpsLog, LammpsDump

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        "test_files", "lammps")


class TestLammpsDump(unittest.TestCase):

    def test_init(self):
        # general tests + gzipped
        rdx_10 = LammpsDump(filename=os.path.join(test_dir, "dump.rdx.gz"))
        np.testing.assert_array_equal(rdx_10.timesteps, np.arange(0, 101, 10))
        obox = rdx_10[0]["box"]
        np.testing.assert_array_equal(obox.bounds, np.array([(35, 48)] * 3))
        atom = rdx_10[-1]["atoms_data"][-1]
        np.testing.assert_array_equal(atom,
                                      [19, 2, 0.42369, 0.47347, 0.555425])
        # timestep wildcard
        rdx_25 = LammpsDump(filename=os.path.join(test_dir, "dump.rdx_wc.*"),
                            parse_box=False)
        self.assertEqual(len(rdx_25), 5)
        np.testing.assert_array_equal(rdx_25.timesteps, np.arange(0, 101, 25))
        self.assertNotIn("box", rdx_25[0])
        # tilted box
        tatb = LammpsDump(filename=os.path.join(test_dir, "dump.tatb"))
        tbox = tatb[0]["box"]
        bounds = [[0, 13.624], [0, 17.1149153805], [0, 15.1826391451]]
        tilt = [-5.75315630927, -6.325466, 7.4257288]
        np.testing.assert_array_almost_equal(tbox.bounds, bounds)
        np.testing.assert_array_almost_equal(tbox.tilt, tilt)



if __name__ == "__main__":
    unittest.main()
