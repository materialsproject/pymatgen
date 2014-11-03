# coding: utf-8
from __future__ import unicode_literals, division, print_function

import yaml
import unittest

from collections import OrderedDict
from monty.collections import AttrDict
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinitio.tasks import ParalConf
from pymatgen.io.abinitio.qadapters import *
from pymatgen.io.abinitio.qadapters import QueueAdapter


class ParseTimestr(PymatgenTest):
    def test_parse_slurm_timestr(self):
        days, hours, minutes, secs = 24*60*60, 60*60, 60, 1
        aequal = self.assertEqual

        # "days-hours",
        aequal(parse_slurm_timestr("2-1"), 2*days + hours)
        # "days-hours:minutes",                                        
        aequal(parse_slurm_timestr("2-1:1"), 2*days + hours + minutes)
        # "days-hours:minutes:seconds".                                
        aequal(parse_slurm_timestr("3-4:2:20"), 3*days + 4*hours + 2*minutes + 20*secs)
        # "minutes",
        aequal(parse_slurm_timestr("10"), 10*minutes)
        # "minutes:seconds",
        aequal(parse_slurm_timestr("3:20"), 3*minutes + 20*secs)
        # "hours:minutes:seconds",
        aequal(parse_slurm_timestr("3:2:5"), 3*hours + 2*minutes + 5*secs)


# TODO
#class PartitionTest(PymatgenTest):
#    def test_partition_methods(self):
#        aequal, atrue, afalse = self.assertEqual, self.assertTrue, self.assertFalse
#
#        # Test mandatory arguments
#        p = AttrDict(qname="test_partition", num_nodes=3, sockets_per_node=2, cores_per_socket=4)
#        with self.assertRaises(ValueError):
#            part = Partition(**p)
#
#        p.update(mem_per_node="8 Mb", timelimit="2:00", priority=1, min_cores=1, max_cores=24)
#        part = Partition(**p)
#        print("partition", str(part))
#        aequal(part.timelimit, 120)
#
#        # Test whether Partition can be serialized with Pickle.
#        deserialized_part = self.serialize_with_pickle(part, test_eq=False)
#
#        # Test properties
#        aequal(part.num_cores, p.num_nodes * p.sockets_per_node * p.cores_per_socket)
#        aequal(part.cores_per_node, p.sockets_per_node * p.cores_per_socket)
#        aequal(part.can_use_omp_threads(p.sockets_per_node * p.cores_per_socket), True)
#        aequal(part.can_use_omp_threads(p.sockets_per_node * p.cores_per_socket + 1), False)
#
#        # Test can_run and distribute
#        #partition has num_nodes=3, sockets_per_node=2, cores_per_socket=4, mem_per_node="8 Mb"
#        afalse(part.can_run(ParalConf(mpi_ncpus=part.num_cores+1, omp_ncpus=1, mem_per_cpu=0.1)))
#        afalse(part.can_run(ParalConf(mpi_ncpus=4, omp_ncpus=9, mem_per_cpu=0.1)))
#        afalse(part.can_run(ParalConf(mpi_ncpus=4, omp_ncpus=1, mem_per_cpu=1024**3)))
#
#        d = part.distribute(mpi_procs=4, omp_threads=1, mem_per_proc=0.1)
#        assert d.num_nodes == 1 and d.mpi_per_node == 4 and d.exact
#
#        d = part.distribute(mpi_procs=16, omp_threads=1, mem_per_proc=1)
#        assert d.num_nodes == 2 and d.mpi_per_node == 8 and d.exact
#
#        # not enough memory per node but can distribute.
#        d = part.distribute(mpi_procs=8, omp_threads=1, mem_per_proc=2)
#        assert d.num_nodes == 2 and d.mpi_per_node == 4 and not d.exact
#
#        # not commensurate with node
#        d = part.distribute(mpi_procs=9, omp_threads=1, mem_per_proc=0.1)
#        assert d.num_nodes == 3 and d.mpi_per_node == 3 and not d.exact
#
#        # mem_per_proc > mem_per_node!
#        with self.assertRaises(part.DistribError):
#            d = part.distribute(mpi_procs=9, omp_threads=1, mem_per_proc=1024)


class QadapterTest(PymatgenTest):
    QDICT = yaml.load("""\
priority: 1
queue:
    qtype: slurm
    qname: Oban
limits:
    timelimit: 0:10:00
    min_cores: 1
    max_cores: 24
    #condition: {"$eq": {omp_threads: 2}}
job:
    modules:
        - intel/compilerpro/13.0.1.117
        - fftw3/intel/3.3
    shell_env:
        PATH: /home/user/tmp_intel13/src/98_main/:/home/user//NAPS/intel13/bin:$PATH
    mpi_runner: mpirun
hardware:
   num_nodes: 100
   sockets_per_node: 2
   cores_per_socket: 4
   mem_per_node: 8 Gb""")

    def test_base(self):
        """
        unit tests for Qadapter subclasses. A more complete coverage would require integration testing.
        """
        aequal, atrue, afalse = self.assertEqual, self.assertTrue, self.assertFalse
        sub_classes = QueueAdapter.__subclasses__()

        # Test if we can instantiate the concrete classes with the abc protocol.
        for subc in sub_classes:
            print("subclass: ", subc)

            # Create the adapter subclass.
            self.QDICT["queue"]["qtype"] = subc.QTYPE
            qad = make_qadapter(**self.QDICT)
            print(qad)

            # Test the programmatic interface used to change job parameters.
            aequal(qad.num_attempts, 0)
            afalse(qad.has_omp)
            atrue(qad.has_mpirun)
            qad.set_mpi_procs(2)
            aequal(qad.mpi_procs, 2)
            atrue(qad.pure_mpi)
            afalse(qad.pure_omp)
            afalse(qad.hybrid_mpi_omp)
            aequal(qad.mem_per_proc, 1024)
            qad.set_mem_per_proc(1024*2)
            aequal(qad.mem_per_proc, 1024*2)

            # Enable OMP
            qad.set_omp_threads(2)
            aequal(qad.omp_threads, 2)
            atrue(qad.has_omp)
            afalse(qad.pure_mpi)
            afalse(qad.pure_omp)
            atrue(qad.hybrid_mpi_omp)

            # Test the creation of the script
            script = qad.get_script_str("job.sh", "/launch/dir", "executable", "qout_path", "qerr_path", 
                                        stdin="STDIN", stdout="STDOUT", stderr="STDERR")

            # Test whether qad can be serialized with Pickle.
            deserialized_qads = self.serialize_with_pickle(qad, test_eq=False)

            for new_qad in deserialized_qads:
                new_script = new_qad.get_script_str("job.sh", "/launch/dir", "executable", "qout_path", "qerr_path", 
                                                    stdin="STDIN", stdout="STDOUT", stderr="STDERR")
                aequal(new_script, script)

class ShellAdapterTest(PymatgenTest):
    """Test suite for Shell adapter."""
    QDICT = yaml.load("""\
priority: 1
queue:
    qname: localhost
    qtype: shell
job:
    mpi_runner: /home/my_mpirun
    pre_run:
        - "source ~/env1.sh"
limits:
    timelimit: 10:00
    min_cores: 1
    max_cores: 1
hardware:
    num_nodes: 1
    sockets_per_node: 1
    cores_per_socket: 1
    mem_per_node: 4 Gb
""")
    def test_methods(self):
        qad = make_qadapter(**self.QDICT)
        print(qad)
        print(qad.mpi_runner)

        assert qad.QTYPE == "shell" and qad.has_mpi and not qad.has_omp
        assert (qad.mpi_procs, qad.omp_threads) == (1, 1)
        assert qad.priority == 1 and qad.num_attempts == 0 and qad.last_attempt is None

        #with self.assertRaises(qad.Error):
        #    qad.set_mpi_procs(2)
        qad.set_omp_threads(1)
        assert qad.has_omp
        #with self.assertRaises(qad.Error):
        #    qad.set_omp_threads(2)

        s = qad.get_script_str("job_name", "/launch_dir", "executable", "qout_path", "qerr_path",
                               stdin="stdin", stdout="stdout", stderr="stderr")
        self.assertMultiLineEqual(s, """\
#!/bin/bash
# OpenMp Environment
export OMP_NUM_THREADS=1

cd /launch_dir
# Commands before execution
source ~/env1.sh

/home/my_mpirun -n 1 executable < stdin > stdout 2> stderr
""")


class SlurmAdapterTest(PymatgenTest):
    """Test suite for Slurm adapter."""
    QDICT = yaml.load("""\
priority: 5
queue:
  qtype: slurm
  qname: Oban
  qparams:
      account: user_account
      mail_user: user@mail.com
limits:
  timelimit: 10:00
  min_cores: 3
  max_cores: 16
job:
  mpi_runner: mpirun
  # pre_run is a string in verbatim mode (note |)
  pre_run: |
      echo ${SLURM_NODELIST}
      ulimit
  modules:
      - intel/compilerpro/13.0.1.117
      - fftw3/intel/3.3
  shell_env:
      PATH: /home/user/bin:$PATH
hardware:
   # Mandatory
   num_nodes: 2
   sockets_per_node: 2
   cores_per_socket: 4
   mem_per_node: 8 Gb
""")

    def test_methods(self):
        self.maxDiff = None

        qad = make_qadapter(**self.QDICT)
        print(qad)
        print(qad.mpi_runner)

        assert qad.QTYPE == "slurm" 
        assert qad.has_mpi and not qad.has_omp
        assert (qad.mpi_procs, qad.omp_threads) == (3, 1)
        assert qad.priority == 5 and qad.num_attempts == 0 and qad.last_attempt is None

        qad.set_mpi_procs(4)

        s = qad.get_script_str("job_name", "/launch_dir", "executable", "qout_path", "qerr_path",
                               stdin="stdin", stdout="stdout", stderr="stderr")
        print(s)
        self.assertMultiLineEqual(s, """\
#!/bin/bash

#SBATCH --partition=Oban
#SBATCH --job-name=job_name
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=1024
#SBATCH --time=0-0:10:0
#SBATCH --account=user_account
#SBATCH --mail-user=user@mail.com
#SBATCH --output=qout_path
#SBATCH --error=qerr_path
# Load Modules
module purge
module load intel/compilerpro/13.0.1.117
module load fftw3/intel/3.3

# Shell Environment
export PATH=/home/user/bin:$PATH

cd /launch_dir
# Commands before execution
echo ${SLURM_NODELIST}
ulimit


mpirun -n 4 executable < stdin > stdout 2> stderr
""")
        #assert 0
        #qad.set_omp_threads(1)
        #assert qad.has_omp


class PbsProadapterTest(PymatgenTest):
    """Test suite for PbsPro adapter."""
    QDICT = yaml.load("""\
priority: 1
queue:
    qtype: pbspro
    qname: fat
    qparams:
        group_list: naps
limits:
    timelimit: 0:0:10
    min_cores: 3
    max_cores: 200
job:
    mpi_runner: mpirun
hardware:
    num_nodes: 100
    sockets_per_node: 2
    cores_per_socket: 4
    mem_per_node: 8 Gb""")

    def test_methods(self):
        self.maxDiff = None
        aequal = self.assertEqual

        qad = make_qadapter(**self.QDICT)
        print(qad)
        print(qad.mpi_runner)

        assert qad.QTYPE == "pbspro" and qad.has_mpi and not qad.has_omp
        assert (qad.mpi_procs, qad.omp_threads) == (3, 1)
        assert qad.priority == 1 and qad.num_attempts == 0 and qad.last_attempt is None

        #qad.set_mpi_procs(4)
        s = qad.get_script_str("job_name", "/launch_dir", "executable", "qout_path", "qerr_path",
                               stdin="stdin", stdout="stdout", stderr="stderr")
        print(s)
        self.assertMultiLineEqual(s, """\
#!/bin/bash

#PBS -q fat
#PBS -N job_name
#PBS -l select=1:ncpus=3:vmem=3072mb:mpiprocs=3
#PBS -l walltime=0:0:10
#PBS -W group_list=naps
# Submission environment
####PBS -V
#PBS -o qout_path
#PBS -e qerr_path
cd /launch_dir
mpirun -n 3 executable < stdin > stdout 2> stderr
""")
        mem = 1024
        qad.set_mem_per_proc(mem)
        print(qad)

        qad.set_mpi_procs(4)
        s, params = qad.get_select(ret_dict=True)
        # IN_CORE PURE MPI: MPI: 4, OMP: 1
        aequal(params, 
          {'ompthreads': 1, 'ncpus': 4, 'chunks': 1, 'mpiprocs': 4, "vmem": mem*4})

        qad.set_omp_threads(2)
        s, params = qad.get_select(ret_dict=True)
        # HYBRID MPI-OPENMP run, perfectly divisible among nodes:  MPI: 4, OMP: 2
        aequal(params, 
            {'vmem': mem*4, 'ncpus': 8, 'chunks': 1, 'ompthreads': 2, 'mpiprocs': 4})

        qad.set_mpi_procs(12)
        s, params = qad.get_select(ret_dict=True)
        # HYBRID MPI-OPENMP run, perfectly divisible among nodes:  MPI: 12, OMP: 2
        aequal(params, 
            {'vmem': mem*4, 'ncpus': 8, 'chunks': 3, 'ompthreads': 2, 'mpiprocs': 4})

        qad.set_omp_threads(5)
        qad.set_mpi_procs(3)
        s, params = qad.get_select(ret_dict=True)
        # HYBRID MPI-OPENMP, NOT commensurate with nodes:  MPI: 3, OMP: 5
        aequal(params, 
            {'vmem': mem, 'ncpus': 5, 'chunks': 3, 'ompthreads': 5, 'mpiprocs': 1})


if __name__ == '__main__':
    import unittest
    unittest.main()
