# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import os
import unittest

import numpy as np
from pymatgen.io.lammps.output import LammpsRun

__author__ = 'Kiran Mathew'
__email__ = 'kmathew@lbl.gov'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestLammpsOutput(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data_file = os.path.join(module_dir, "test_files", "nvt.data")
        traj_file = os.path.join(module_dir, "test_files", "nvt.dump")
        log_file = os.path.join(module_dir, "test_files", "nvt.log")
        cls.lammpsrun = LammpsRun(data_file, traj_file, log_file,
                                  is_forcefield=True)

    def test_lammps_log(self):
        fields = "step vol temp press ke pe etotal enthalpy evdwl ecoul epair ebond eangle edihed eimp " \
                 "emol elong etail lx ly lz xy xz yz density"
        fields = fields.split()
        thermo_data_ans = np.loadtxt(
            os.path.join(module_dir, "test_files", "thermo_data.txt"))
        thermo_data = self.lammpsrun.lammps_log.thermo_data
        self.assertEqual(list(thermo_data.dtype.names), fields)
        self.assertEqual(self.lammpsrun.lammps_log.nmdsteps + 1,
                         thermo_data.shape[0])
        thermo_data = thermo_data[:].view(np.float64).reshape(
            thermo_data.shape + (-1,))
        np.testing.assert_almost_equal(thermo_data, thermo_data_ans,
                                       decimal=10)

    def test_lammps_trajectory(self):
        fields = "Atoms_id atom_type x y z vx vy vz mol mass"
        fields = fields.split()
        timestep_ans = 82
        trajectory_ans = np.loadtxt(os.path.join(module_dir, "test_files",
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


if __name__ == "__main__":
    unittest.main()
