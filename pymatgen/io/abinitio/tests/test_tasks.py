from __future__ import division, print_function
import os

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinitio.tasks import *
from pymatgen.io.abinitio.tasks import TaskPolicy

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


class ParalHintsTest(PymatgenTest):
    def test_base(self):
        """Testing ParalHints."""
        s = \
"""--- !Autoparal
#Autoparal section for Sigma runs.
info:
    autoparal: 1
    max_ncpus: 4
    nkpt: 6
    nsppol: 1
    nspinor: 1
    nbnds: 10
configurations:
    - tot_ncpus: 1
      mpi_ncpus: 1
      efficiency:  1.000000000
      mem_per_cpu:        11.54
    - tot_ncpus: 2
      mpi_ncpus: 2
      efficiency:  1.000000000
      mem_per_cpu:         7.42
    - tot_ncpus: 3
      mpi_ncpus: 3
      efficiency:  0.833333333
      mem_per_cpu:         6.60
    - tot_ncpus: 4
      mpi_ncpus: 4
      efficiency:  0.833333333
      mem_per_cpu:         15.77
...
"""
        tmpfile = self.tmpfile_write(s)
        aequal = self.assertEqual

        # Parse the file with the configurations.
        confs = ParalHintsParser().parse(tmpfile)
        #print("all_confs:\n", confs)

        # When autoparal is 1, max_ncpus must be specified
        with self.assertRaises(ValueError):
            policy = TaskPolicy(autoparal=1)
            optimal = confs.select_optimal_conf(self, policy)

        # Optimize speedup with ncpus <= max_ncpus
        policy = TaskPolicy(autoparal=1, max_ncpus=3)
        optimal = confs.select_optimal_conf(policy)
        aequal(optimal.tot_ncpus, 3)

        # Optimize speedup with ncpus <= max_ncpus and condition on efficiency.
        policy = TaskPolicy(autoparal=1, max_ncpus=4, condition={"efficiency": {"$ge": 0.9}})
        optimal = confs.select_optimal_conf(policy)
        aequal(optimal.tot_ncpus, 2)

        # Optimize speedup with ncpus <= max_ncpus and conditions on efficiency and mem_per_cpu.
        policy = TaskPolicy(autoparal=1, mode="default", max_ncpus=4, 
            condition={"$and": [{"efficiency": {"$ge": 0.8}}, {"mem_per_cpu": {"$le": 7.0}}] })
        optimal = confs.select_optimal_conf(policy)
        aequal(optimal.tot_ncpus, 3)

        # If no configuration satisfies the constraints, we return the conf with the highest speedup.
        policy = TaskPolicy(autoparal=1, max_ncpus=4, condition={"efficiency": {"$ge": 100}})
        optimal = confs.select_optimal_conf(policy)
        aequal(optimal.tot_ncpus, 4)

        # Wrong conditions --> dump a warning and return the conf with the highest speedup.
        policy = TaskPolicy(autoparal=1, max_ncpus=4, condition={"foobar": {"$ge": 100}})
        optimal = confs.select_optimal_conf(policy)
        aequal(optimal.tot_ncpus, 4)


if __name__ == '__main__':
    import unittest
    unittest.main()
