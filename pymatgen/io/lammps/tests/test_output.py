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

    @classmethod
    def setUpClass(cls):
        cls.rdx_10 = LammpsDump(filename=os.path.join(test_dir,
                                                      "dump.rdx.gz"))
        cls.rdx_25 = LammpsDump(filename=os.path.join(test_dir,
                                                      "dump.rdx_wc.*"),
                                parse_box=False)
        cls.tatb = LammpsDump(filename=os.path.join(test_dir, "dump.tatb"))

    def test_init(self):
        files = ["dump.rdx_wc.0", "dump.rdx_wc.25", "dump.rdx_wc.50",
                 "dump.rdx_wc.75", "dump.rdx_wc.100"]
        self.assertListEqual([os.path.join(test_dir, f) for f in files],
                             self.rdx_25.all_files)

    def test_read(self):
        # general tests + gzipped
        rdx_10_data = list(self.rdx_10.read())
        timesteps_10 = [d["timestep"] for d in rdx_10_data]
        np.testing.assert_array_equal(timesteps_10, np.arange(0, 101, 10))
        obox = rdx_10_data[0]["box"]
        np.testing.assert_array_equal(obox.bounds, np.array([(35, 48)] * 3))
        atom = rdx_10_data[-1]["data"][-1]
        np.testing.assert_array_equal(atom,
                                      [19, 2, 0.42369, 0.47347, 0.555425])
        # timestep wildcard
        rdx_25_data = list(self.rdx_25.read())
        timesteps_25 = [d["timestep"] for d in rdx_25_data]
        np.testing.assert_array_equal(timesteps_25, np.arange(0, 101, 25))
        self.assertNotIn("box", rdx_25_data[0])
        # tilted box
        tbox = list(self.tatb.read())[0]["box"]
        bounds = [[0, 13.624], [0, 17.1149153805], [0, 15.1826391451]]
        tilt = [-5.75315630927, -6.325466, 7.4257288]
        np.testing.assert_array_almost_equal(tbox.bounds, bounds)
        np.testing.assert_array_almost_equal(tbox.tilt, tilt)


if __name__ == "__main__":
    unittest.main()
