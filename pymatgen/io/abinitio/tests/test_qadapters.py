# coding: utf-8
from __future__ import unicode_literals, division, print_function

import yaml

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


class PartitionTest(PymatgenTest):
    def test_partition_methods(self):
        aequal, atrue, afalse = self.assertEqual, self.assertTrue, self.assertFalse

        # Test mandatory arguments
        p = AttrDict(qname="test_partition", num_nodes=3, sockets_per_node=2, cores_per_socket=4)
        with self.assertRaises(ValueError):
            part = Partition(**p)

        p.update(mem_per_node="8 Mb", timelimit="2:00", priority=1, min_cores=1, max_cores=24)
        part = Partition(**p)
        print("partition", str(part))
        aequal(part.timelimit, 120)

        # Test whether Partition can be serialized with Pickle.
        deserialized_part = self.serialize_with_pickle(part, test_eq=False)

        # Test properties
        aequal(part.num_cores, p.num_nodes * p.sockets_per_node * p.cores_per_socket)
        aequal(part.cores_per_node, p.sockets_per_node * p.cores_per_socket)
        aequal(part.can_use_omp_threads(p.sockets_per_node * p.cores_per_socket), True)
        aequal(part.can_use_omp_threads(p.sockets_per_node * p.cores_per_socket + 1), False)

        # Test can_run and distribute
        #partition has num_nodes=3, sockets_per_node=2, cores_per_socket=4, mem_per_node="8 Mb"
        afalse(part.can_run(ParalConf(mpi_ncpus=part.num_cores+1, omp_ncpus=1, mem_per_cpu=0.1)))
        afalse(part.can_run(ParalConf(mpi_ncpus=4, omp_ncpus=9, mem_per_cpu=0.1)))
        afalse(part.can_run(ParalConf(mpi_ncpus=4, omp_ncpus=1, mem_per_cpu=1024**3)))

        d = part.distribute(mpi_procs=4, omp_threads=1, mem_per_proc=0.1)
        assert d.num_nodes == 1 and d.mpi_per_node == 4 and d.exact

        d = part.distribute(mpi_procs=16, omp_threads=1, mem_per_proc=1)
        assert d.num_nodes == 2 and d.mpi_per_node == 8 and d.exact

        # not enough memory per node but can distribute.
        d = part.distribute(mpi_procs=8, omp_threads=1, mem_per_proc=2)
        assert d.num_nodes == 2 and d.mpi_per_node == 4 and not d.exact

        # not commensurate with node
        d = part.distribute(mpi_procs=9, omp_threads=1, mem_per_proc=0.1)
        assert d.num_nodes == 3 and d.mpi_per_node == 3 and not d.exact

        # mem_per_proc > mem_per_node!
        with self.assertRaises(part.DistribError):
            d = part.distribute(mpi_procs=9, omp_threads=1, mem_per_proc=1024)

        # TODO: Test OpenMP


class QadapterTest(PymatgenTest):

    def test_base(self):
        """
        unit tests for Qadapter subclasses. A more complete coverage would require integration testing.
        """
        aequal, atrue, afalse = self.assertEqual, self.assertTrue, self.assertFalse
        sub_classes = QueueAdapter.__subclasses__()

        modules = ["intel/compilerpro/13.0.1.117",
                   "fftw3/intel/3.3",
        ]

        shell_env = OrderedDict(
             [("PATH", "/user/abinit-7.4.3-public/tmp_intel13/src/98_main/:/user/bin:$PATH"),
              ("LD_LIBRARY_PATH", "/NAPS/intel13/lib:$LD_LIBRARY_PATH")])

        mpi_runner = MpiRunner("mpirun")

        #omp_env = OmpEnv(OMP_NUM_THREADS=2)
        qname="test_partition",
        partition = dict(num_nodes=100, sockets_per_node=2, cores_per_socket=4, mem_per_node="1 Gb", 
                         timelimit=10, min_cores=1, max_cores=24, priority=1)

        # Test if we can instantiate the concrete classes with the abc protocol.
        for subc in sub_classes:
            print("subclass: ", subc)

            # Create the adapter
            qad = make_qadapter(qtype=subc.QTYPE, qname="test_queue", setup=None, modules=modules, shell_env=shell_env, omp_env=None, 
                                pre_run=None, post_run=None, mpi_runner=mpi_runner, partition=partition)

            # Test the programmatic interface used to change job parameters.
            aequal(qad.num_attempts, 0)
            afalse(qad.has_omp)
            atrue(qad.has_mpirun)
            qad.set_mpi_procs(2)
            aequal(qad.mpi_procs, 2)
            atrue(qad.pure_mpi)
            afalse(qad.pure_omp)
            afalse(qad.hybrid_mpi_omp)
            # TODO
            #qad.set_mem_per_proc(13)
            #aequal(qad.mem_per_proc, 13)

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


class ShelladapterTest(PymatgenTest):
    """Test suite for Shell adapter."""
    def test_methods(self):
        conf = """
          qtype: shell
          qname: localhost
          mpi_runner: /home/my_mpirun
          pre_run:
             - "source ~/env1.sh"
          partition:
             # Mandatory
             num_nodes: 1
             sockets_per_node: 1
             cores_per_socket: 1
             mem_per_node: 4 Gb
             timelimit: 10:00
             priority: 1
             min_cores: 1
             max_cores: 1
        """
        d = yaml.load(conf)
        qtype = d.pop("qtype")
        qad = make_qadapter(qtype, **d)
        print(qad)
        print(qad.mpi_runner)

        assert qad.QTYPE == "shell" and qad.supported_qparams == ["MPI_PROCS"]
        assert qad.has_mpi and not qad.has_omp
        assert (qad.mpi_procs, qad.omp_threads) == (1, 1)
        assert qad.priority == 1 and qad.num_attempts == 0 and qad.last_attempt is None

        qad.set_mpi_procs(2)
        qad.set_omp_threads(4)
        assert qad.has_omp
        s = qad.get_script_str("job_name", "/launch_dir", "executable", "qout_path", "qerr_path",
                               stdin="stdin", stdout="stdout", stderr="stderr")
        self.assertMultiLineEqual(s, """\
#!/bin/bash

export MPI_PROCS=2
# OpenMp Environment
export OMP_NUM_THREADS=4

cd /launch_dir
# Commands before execution
source ~/env1.sh

/home/my_mpirun -n 2 executable < stdin > stdout 2> stderr
""")


class SlurmadapterTest(PymatgenTest):
    """Test suite for Slurm adapter."""
    def test_methods(self):
        conf = """
          qtype: slurm
          qname: Oban
          qparams:
            account: user_account
            mail_user: user@mail.com
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
          partition:
             # Mandatory
             num_nodes: 2
             sockets_per_node: 2
             cores_per_socket: 4
             mem_per_node: 8 Gb
             timelimit: 10:00
             priority: 5
             min_cores: 3
             max_cores: 16
        """
        self.maxDiff = None

        d = yaml.load(conf)
        qtype = d.pop("qtype")
        qad = make_qadapter(qtype, **d)
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
    def test_methods(self):
        self.maxDiff = None
        d = yaml.load("""\
          qtype: pbspro
          qname: fat
          qparams:
             group_list: naps
          mpi_runner: mpirun
          partition:
             # Mandatory
             num_nodes: 2
             sockets_per_node: 2
             cores_per_socket: 4
             mem_per_node: 8 Gb
             timelimit: 1-0:0:0
             priority: 5
             min_cores: 3
             max_cores: 16
        """)

        qtype = d.pop("qtype")
        qad = make_qadapter(qtype, **d)
        print(qad)
        print(qad.mpi_runner)

        assert qad.QTYPE == "pbspro" 
        assert qad.has_mpi and not qad.has_omp
        #assert (qad.mpi_procs, qad.omp_threads) == (3, 1)
        assert qad.priority == 5 and qad.num_attempts == 0 and qad.last_attempt is None

        #qad.set_mpi_procs(4)

        s = qad.get_script_str("job_name", "/launch_dir", "executable", "qout_path", "qerr_path",
                               stdin="stdin", stdout="stdout", stderr="stderr")
        print(s)
        self.assertMultiLineEqual(s, """\
#!/bin/bash

#PBS -q fat
#PBS -N job_name
#PBS -l select=1:ncpus=3:vmem=3072mb:mpiprocs=3
#PBS -l walltime=24:0:0
#PBS -W group_list=naps
# Submission environment
#PBS -V
#PBS -o qout_path
#PBS -e qerr_path
cd /launch_dir
mpirun -n 3 executable < stdin > stdout 2> stderr
""")
        #assert 0

    def test_optimize_params(self):
        aequal = self.assertEqual

        partition = dict(num_nodes=100, sockets_per_node=2, 
                      cores_per_socket=4, mem_per_node="8 Gb", timelimit=10, priority=1, 
                      min_cores=1, max_cores=200)

        qad = make_qadapter(qtype="pbspro", qname="test_pbspro", partition=partition)
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
