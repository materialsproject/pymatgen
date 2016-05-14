# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
from __future__ import unicode_literals, division, print_function

import yaml
import unittest2 as unittest

from collections import OrderedDict
from monty.collections import AttrDict
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinit.tasks import ParalConf
from pymatgen.io.abinit.qadapters import *
from pymatgen.io.abinit.qadapters import QueueAdapter, SlurmAdapter 
from pymatgen.io.abinit import qutils as qu

class ParseTimestr(PymatgenTest):
    def test_slurm_parse_timestr(self):
        days, hours, minutes, secs = 24*60*60, 60*60, 60, 1
        aequal = self.assertEqual

        slurm_parse_timestr = qu.slurm_parse_timestr 

        # "days-hours",
        aequal(slurm_parse_timestr("2-1"), 2*days + hours)
        # "days-hours:minutes",                                        
        aequal(slurm_parse_timestr("2-1:1"), 2*days + hours + minutes)
        # "days-hours:minutes:seconds".                                
        aequal(slurm_parse_timestr("3-4:2:20"), 3*days + 4*hours + 2*minutes + 20*secs)
        # "minutes",
        aequal(slurm_parse_timestr("10"), 10*minutes)
        # "minutes:seconds",
        aequal(slurm_parse_timestr("3:20"), 3*minutes + 20*secs)
        # "hours:minutes:seconds",
        aequal(slurm_parse_timestr("3:2:5"), 3*hours + 2*minutes + 5*secs)


class QadapterTest(PymatgenTest):
    QDICT = yaml.load("""\
priority: 1
queue:
    qtype: slurm
    qname: Oban
limits:
    timelimit: 2:00
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
   num_nodes: 3
   sockets_per_node: 2
   cores_per_socket: 4
   mem_per_node: 8 Gb
""")

    def test_base(self):
        """unit tests for Qadapter subclasses. A more complete coverage would require integration testing."""
        self.maxDiff = None
        aequal, atrue, afalse = self.assertEqual, self.assertTrue, self.assertFalse
        sub_classes = QueueAdapter.__subclasses__()

        # Test if we can instantiate the concrete classes with the abc protocol.
        for subc in sub_classes:
            print("subclass: ", subc)

            # Create the adapter subclass.
            self.QDICT["queue"]["qtype"] = subc.QTYPE
            qad = make_qadapter(**self.QDICT)
            print(qad)
            hw = qad.hw
            giga = 1024

            # Test the programmatic interface used to change job parameters.
            aequal(qad.num_launches, 0)
            afalse(qad.has_omp)
            atrue(qad.has_mpi)
            qad.set_mpi_procs(2)
            aequal(qad.mpi_procs, 2)
            atrue(qad.pure_mpi)
            afalse(qad.pure_omp)
            afalse(qad.hybrid_mpi_omp)
            aequal(qad.mem_per_proc, giga)
            qad.set_mem_per_proc(2 * giga)
            aequal(qad.mem_per_proc, 2 * giga)
            aequal(qad.timelimit, 120)

            # Enable OMP
            qad.set_omp_threads(2)
            aequal(qad.omp_threads, 2)
            atrue(qad.has_omp)
            afalse(qad.pure_mpi)
            afalse(qad.pure_omp)
            atrue(qad.hybrid_mpi_omp)

            atrue(qad.hw.can_use_omp_threads(hw.sockets_per_node * hw.cores_per_socket))
            afalse(qad.hw.can_use_omp_threads(hw.sockets_per_node * hw.cores_per_socket + 1))

            # Test the creation of the script
            script = qad.get_script_str("job.sh", "/launch/dir", "executable", "qout_path", "qerr_path", 
                                        stdin="STDIN", stdout="STDOUT", stderr="STDERR")

            # Test whether qad can be serialized with Pickle.
            deserialized_qads = self.serialize_with_pickle(qad, test_eq=False)

            for new_qad in deserialized_qads:
                new_script = new_qad.get_script_str("job.sh", "/launch/dir", "executable", "qout_path", "qerr_path", 
                                                    stdin="STDIN", stdout="STDOUT", stderr="STDERR")
                aequal(new_script, script)

            # Test can_run and distribute
            # The hardware has num_nodes=3, sockets_per_node=2, cores_per_socket=4, mem_per_node="8 Gb"
            afalse(qad.can_run_pconf(ParalConf(mpi_ncpus=hw.num_cores+1, omp_ncpus=1, mem_per_cpu=0.1)))
            afalse(qad.can_run_pconf(ParalConf(mpi_ncpus=4, omp_ncpus=9, mem_per_cpu=0.1)))
            afalse(qad.can_run_pconf(ParalConf(mpi_ncpus=4, omp_ncpus=1, mem_per_cpu=10 * giga)))

            d = qad.distribute(mpi_procs=4, omp_threads=1, mem_per_proc=giga)
            assert d.num_nodes == 1 and d.mpi_per_node == 4 and d.exact

            d = qad.distribute(mpi_procs=16, omp_threads=1, mem_per_proc=giga)
            assert d.num_nodes == 2 and d.mpi_per_node == 8 and d.exact

            # not enough memory per node but can distribute.
            d = qad.distribute(mpi_procs=8, omp_threads=1, mem_per_proc=2 * giga)
            assert d.num_nodes == 2 and d.mpi_per_node == 4 and not d.exact

            # mem_per_proc > mem_per_node!
            with self.assertRaises(qad.Error):
                d = qad.distribute(mpi_procs=9, omp_threads=1, mem_per_proc=10 * giga)

            # TODO
            # not commensurate with node
            #d = qad.distribute(mpi_procs=9, omp_threads=1, mem_per_proc=giga)
            #assert d.num_nodes == 3 and d.mpi_per_node == 3 and not d.exact

            with self.assertRaises(qad.Error): 
                qad.set_mpi_procs(25)
                qad.validate()
            with self.assertRaises(qad.Error): 
                qad.set_mpi_procs(100)
                qad.validate()
            with self.assertRaises(qad.Error): 
                qad.set_omp_threads(10)
                qad.validate()
            with self.assertRaises(qad.Error): 
                qad.set_mem_per_proc(9 * giga)
                qad.validate()

        # Test if one can register a customized class.
        class MyAdapter(SlurmAdapter):
            QTYPE = "myslurm"

        SlurmAdapter.register(MyAdapter)
        assert issubclass(MyAdapter, QueueAdapter)

        self.QDICT["queue"]["qtype"] = "myslurm"
        qad = make_qadapter(**self.QDICT)
        assert isinstance(qad, MyAdapter)


class ShellAdapterTest(PymatgenTest):
    """Test suite for Shell adapter."""
    QDICT = yaml.load("""\
priority: 1
queue:
    qname: localhost
    qtype: shell
job:
    mpi_runner: /home/local/bin/mpirun
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
        assert qad.priority == 1 and qad.num_launches == 0 and qad.last_launch is None

        qad.set_omp_threads(1)
        assert qad.has_omp

        s = qad.get_script_str("job_name", "/launch_dir", "executable", "qout_path", "qerr_path",
                               stdin="stdin", stdout="stdout", stderr="stderr")
        self.assertMultiLineEqual(s, """\
#!/bin/bash
cd /launch_dir
# OpenMp Environment
export OMP_NUM_THREADS=1

# Commands before execution
source ~/env1.sh

/home/local/bin/mpirun -n 1 executable < stdin > stdout 2> stderr
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
  setup:
      - echo ${SLURM_JOB_NODELIST}
      - ulimit -s unlimited
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
        assert qad.priority == 5 and qad.num_launches == 0 and qad.last_launch is None

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
cd /launch_dir
# Setup section
echo ${SLURM_JOB_NODELIST}
ulimit -s unlimited

# Load Modules
module purge
module load intel/compilerpro/13.0.1.117 2>> mods.err
module load fftw3/intel/3.3 2>> mods.err

# OpenMp Environment
export OMP_NUM_THREADS=1
# Shell Environment
export PATH=/home/user/bin:$PATH

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
    QDICT_SHARED = yaml.load("""\
priority: 1
queue:
    qtype: pbspro
    qname: fat_shared
    qnodes: shared
    qparams:
        group_list: naps
limits:
    timelimit: 0:0:10
    min_cores: 3
    max_cores: 200
    min_mem_per_proc: 2000
master_mem_overhead: 1000
job:
    mpi_runner: mpirun
hardware:
    num_nodes: 100
    sockets_per_node: 2
    cores_per_socket: 12
    mem_per_node: 48000 Mb""")
    QDICT_EXCLUSIVE = yaml.load("""\
priority: 1
queue:
    qtype: pbspro
    qname: fat_exclusive
    qnodes: exclusive
    qparams:
        group_list: naps
limits:
    timelimit: 0:0:10
    min_cores: 3
    max_cores: 200
    min_mem_per_proc: 2000
master_mem_overhead: 1000
job:
    mpi_runner: mpirun
hardware:
    num_nodes: 100
    sockets_per_node: 2
    cores_per_socket: 12
    mem_per_node: 48000 Mb""")

    def test_methods(self):
        self.maxDiff = None
        aequal = self.assertEqual

        qad = make_qadapter(**self.QDICT)
        self.assertMSONable(qad)
        print(qad)
        print(qad.mpi_runner)

        assert qad.QTYPE == "pbspro" and qad.has_mpi and not qad.has_omp
        assert (qad.mpi_procs, qad.omp_threads) == (3, 1)
        assert qad.priority == 1 and qad.num_launches == 0 and qad.last_launch is None

        #qad.set_mpi_procs(4)
        s = qad.get_script_str("job_name", "/launch_dir", "executable", "qout_path", "qerr_path",
                               stdin="stdin", stdout="stdout", stderr="stderr")
        print(s)
        self.assertMultiLineEqual(s, """\
#!/bin/bash

#PBS -q fat
#PBS -N job_name
#PBS -l select=3:ncpus=1:vmem=1024mb:mpiprocs=1
#PBS -l pvmem=1024mb
#PBS -l walltime=0:0:10
#PBS -W group_list=naps
#PBS -o qout_path
#PBS -e qerr_path
cd /launch_dir
# OpenMp Environment
export OMP_NUM_THREADS=1
mpirun -n 3 executable < stdin > stdout 2> stderr
""")
        mem = 1024
        qad.set_mem_per_proc(mem)
        print(qad)

        qad.set_mpi_procs(4)
        s, params = qad.get_select(ret_dict=True)
        # IN_CORE PURE MPI: MPI: 4, OMP: 1
        aequal(params, 
          {'ncpus': 1, 'chunks': 4, 'mpiprocs': 1, "vmem": mem})

        qad.set_omp_threads(2)
        s, params = qad.get_select(ret_dict=True)
        # HYBRID MPI-OPENMP run, perfectly divisible among nodes:  MPI: 4, OMP: 2
        aequal(params, 
            {'vmem': mem, 'ncpus': 2, 'chunks': 4, 'ompthreads': 2, 'mpiprocs': 1})

        qad.set_mpi_procs(12)
        s, params = qad.get_select(ret_dict=True)
        # HYBRID MPI-OPENMP run, perfectly divisible among nodes:  MPI: 12, OMP: 2
        aequal(params, 
            {'vmem': mem, 'ncpus': 2, 'chunks': 12, 'ompthreads': 2, 'mpiprocs': 1})

        qad.set_omp_threads(5)
        qad.set_mpi_procs(3)
        s, params = qad.get_select(ret_dict=True)
        # HYBRID MPI-OPENMP, NOT commensurate with nodes:  MPI: 3, OMP: 5
        aequal(params, 
            {'vmem': mem, 'ncpus': 5, 'chunks': 3, 'ompthreads': 5, 'mpiprocs': 1})

        # Testing the handling of master memory overhead
        # Shared mode (the nodes might be shared amongst different jobs from different users)
        qad_shared = make_qadapter(**self.QDICT_SHARED)
        aequal(qad_shared.hw.mem_per_node, 48000)
        qad_shared.set_mpi_procs(15)
        qad_shared.set_mem_per_proc(6000)
        aequal(qad_shared.get_select(), '1:ncpus=1:vmem=7000mb:mpiprocs=1+'
                                        '14:ncpus=1:vmem=6000mb:mpiprocs=1')
        qad_shared.set_mpi_procs(64)
        qad_shared.set_mem_per_proc(3500)
        qad_shared.set_master_mem_overhead(4000)
        self.assertMSONable(qad_shared)
        aequal(qad_shared.get_select(), '1:ncpus=1:vmem=7500mb:mpiprocs=1+'
                                        '63:ncpus=1:vmem=3500mb:mpiprocs=1')

        # Exclusive mode (the nodes are attributed exclusively to a given user)
        qad_exclusive = make_qadapter(**self.QDICT_EXCLUSIVE)
        aequal(qad_exclusive.hw.mem_per_node, 48000)
        qad_exclusive.set_mpi_procs(47)
        qad_exclusive.set_mem_per_proc(2000)
        qad_exclusive.set_master_mem_overhead(1)
        self.assertMSONable(qad_exclusive)
        aequal(qad_exclusive.get_select(), '1:ncpus=23:vmem=48000mb:mpiprocs=23+'
                                           '1:ncpus=24:vmem=48000mb:mpiprocs=24')
        qad_exclusive.set_mpi_procs(48)
        aequal(qad_exclusive.get_select(), '1:ncpus=1:vmem=48000mb:mpiprocs=1+'
                                           '1:ncpus=24:vmem=48000mb:mpiprocs=24+'
                                           '1:ncpus=23:vmem=48000mb:mpiprocs=23')
        qad_exclusive.set_mpi_procs(50)
        aequal(qad_exclusive.get_select(), '1:ncpus=2:vmem=48000mb:mpiprocs=2+'
                                           '2:ncpus=24:vmem=48000mb:mpiprocs=24')

if __name__ == '__main__':
    import unittest2 as unittest
    unittest.main()
