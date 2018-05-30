# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

import numpy as np

from pymatgen.io.lammps.output import LammpsRun, LammpsLog, LammpsDump

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


class TestLammpsRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data_file = os.path.join(test_dir, "nvt.data")
        traj_file = os.path.join(test_dir, "nvt.dump")
        log_file = os.path.join(test_dir, "nvt.log")
        cls.lmps_log = LammpsLog(log_file=log_file)
        cls.lammpsrun = LammpsRun(data_file, traj_file, log_file)

    def test_lammps_log(self):
        fields = "step vol temp press ke pe etotal enthalpy evdwl ecoul epair " \
                 "ebond eangle edihed eimp " \
                 "emol elong etail lx ly lz xy xz yz density"
        fields = fields.split()
        thermo_data_ans = np.loadtxt(
            os.path.join(test_dir, "thermo_data.txt"))
        thermo_data = self.lammpsrun.log.thermo_data
        self.assertEqual(sorted(list(thermo_data.keys())), sorted(fields))
        self.assertEqual(self.lammpsrun.log.nmdsteps + 1,
                         len(thermo_data['step']))
        data = [thermo_data[k] for k in fields]
        np.testing.assert_almost_equal(data, np.transpose(thermo_data_ans), decimal=10)

    def test_lammps_trajectory(self):
        fields = "Atoms_id atom_type x y z vx vy vz mol mass"
        fields = fields.split()
        timestep_ans = 82
        trajectory_ans = np.loadtxt(os.path.join(test_dir,
                                                 "trajectory_timestep_82_sorted.txt"))
        begin = int(timestep_ans / 2) * self.lammpsrun.natoms
        end = (int(timestep_ans / 2) + 1) * self.lammpsrun.natoms
        trajectory = self.lammpsrun.trajectory[begin:end]
        # atom ids in the trajectory starts from 0
        np.testing.assert_almost_equal(trajectory[:][fields[0]],
                                       trajectory_ans[:, 0] - 1, decimal=10)
        for i, fld in enumerate(fields[1:]):
            np.testing.assert_almost_equal(trajectory[:][fld],
                                           trajectory_ans[:, i + 1],
                                           decimal=10)

    def test_get_structures_from_trajectory(self):
        structures = self.lammpsrun.get_structures_from_trajectory()
        self.assertEqual(len(structures), len(self.lammpsrun.timesteps))

    def test_get_displacements(self):
        structure, disp = self.lammpsrun.get_displacements()
        self.assertEqual(disp.shape[0], len(structure))
        self.assertEqual(disp.shape[1], len(self.lammpsrun.timesteps) - 1)
        self.assertEqual(disp.shape[2], 3)
        self.assertAlmostEqual(disp[-1, -1, -1], 0.077079999999999788)

    def test_serialization(self):
        d = self.lammpsrun.as_dict()
        lmps_run = LammpsRun.from_dict(d)
        self.assertDictEqual(d, lmps_run.as_dict())
        d2 = self.lmps_log.as_dict()
        lmps_log = LammpsLog.from_dict(d2)
        self.assertDictEqual(d2, lmps_log.as_dict())


if __name__ == "__main__":
    unittest.main()
