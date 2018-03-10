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

    def setUp(self):
        dump_file_1 = os.path.join(test_dir, "dump_1")
        dump_file_2 = os.path.join(test_dir, "dump_2")
        self.dump1 = LammpsDump.from_file(dump_file_1)
        self.dump2 = LammpsDump.from_file(dump_file_2)
        dump_file_nvt = os.path.join(test_dir, "nvt.dump")
        self.dump_nvt = LammpsDump.from_file(dump_file_nvt)

    def test_non_atoms_data(self):
        self.assertEqual(self.dump1.box_bounds, [[0.0, 25.0],
                                                 [0.0, 25.0],
                                                 [0.0, 25.0]])
        self.assertEqual(self.dump1.natoms, 123)
        self.assertEqual(self.dump1.timesteps, [0.0])
        self.assertEqual(self.dump1.box_bounds, self.dump2.box_bounds)
        self.assertEqual(self.dump1.natoms, self.dump2.natoms)
        self.assertEqual(self.dump1.timesteps, self.dump2.timesteps)

    def test_atoms_data(self):
        self.assertEqual(self.dump1.natoms, len(self.dump1.atoms_data))
        self.assertEqual(self.dump2.natoms, len(self.dump2.atoms_data))

        ans1 = [float(x) for x in "1 2 1 3 3 2 4 1 3 5 4 3 5 2 4 6 1 3 5 7 5 4 " \
                                  "6 3 5 7 2 4 6 1 3 5 7 6 5 7 4 6 3 5 7 2 4 6 " \
                                  "1 3 5 7 7 6 5 7 4 6 3 5 7 2 4 6 1 3 5 7".split()]
        np.testing.assert_almost_equal(self.dump1.atoms_data[115], ans1,
                                       decimal=6)
        np.testing.assert_almost_equal(self.dump2.atoms_data[57],
                                       [99, 2, 0.909816, 0.883438, 0.314853],
                                       decimal=6)

    def test_timesteps_and_atoms_data(self):
        self.assertEqual(self.dump_nvt.natoms * len(self.dump_nvt.timesteps),
                         len(self.dump_nvt.atoms_data))


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
