#!/usr/bin/env python
from __future__ import division, print_function
import os

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinitio.tasks import *

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", 'test_files')

class TaskManagerTest(PymatgenTest):

    def test_base(self):
        """
        Simple unit tests for Qadapter subclasses.
        A more complete coverage would require integration testing.
        """
        # Initialize the object from YAML file.
        slurm_manager = TaskManager.from_file(os.path.join(test_dir, "taskmanager.yml"))

        print(slurm_manager)
        self.assertTrue(slurm_manager.tot_ncpus == 2)
        self.assertTrue(slurm_manager.mpi_ncpus == 2)
        self.assertTrue(slurm_manager.omp_ncpus == 1)

        # Make a simple shell manager that will inherit the initial configuration.
        shell_manager = slurm_manager.to_shell_manager(mpi_ncpus=1)
        self.assertTrue(shell_manager.tot_ncpus == 1)
        self.assertTrue(shell_manager.mpi_ncpus == 1)

        # check that the initial slurm_manger has not been modified
        self.assertTrue(slurm_manager.tot_ncpus == 2)

        # Test pickle
        self.serialize_with_pickle(slurm_manager, test_eq=False)


if __name__ == '__main__':
    import unittest
    unittest.main()
