#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals

import sys
import os
import argparse
import time

from pymatgen.core.units import Memory
from pymatgen.io.abinitio.qadapters import SlurmJob, PbsProJob, make_qadapter


def show_job():
    job_id = sys.argv[1]
    #job = SlurmJob(job_id)
    job = PbsProJob(job_id)

    print("start_time: ", job.estimated_start_time())
    print("get_info:" ,job.get_info())
    print("get_stats:" ,job.get_stats())

    #if job.unknown_status:
    #    job.parse_qout()

    #if job.is_completed:
    #if job.timeout
    #if job.is_running:
    #if job.is_failed:
    #if job.has_node_failures:
    #  nodes = job.save_nodes_in_blacklist()

def show_slurm():
    import yaml
    QDICT = yaml.load("""\
priority: 5
queue:
  qtype: slurm
  qname: Oban
  qparams:
      mail_user: user@mail.com
limits:
  timelimit: 10:00
  min_cores: 1
  max_cores: 16
job:
  shell_env:
      PATH: /home/user/bin:$PATH
hardware:
   # Mandatory
   num_nodes: 2
   sockets_per_node: 2
   cores_per_socket: 4
   mem_per_node: 8 Gb
""")
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-p', '--mpi-procs', default=1, type=int, help='Number of MPI processes')
    parser.add_argument('-o', '--omp-threads', default=1, type=int, help='Number of OpenMP threads')
    parser.add_argument('-m', '--mem_per_proc', type=lambda s: float(Memory.from_string(s).to("Mb")), 
                        default="1Mb", help='Memory per processor in Mb')

    opts = parser.parse_args()
    print(opts.mem_per_proc)

    qad = make_qadapter(**QDICT)
    qad.set_mpi_procs(opts.mpi_procs)
    qad.set_omp_threads(opts.omp_threads)
    qad.set_mem_per_proc(opts.mem_per_proc)

    script = qad.get_script_str("job_name", "launch_dir", "hostname", "qout_path", "qerr_path",
                                stdin=None, stdout=None, stderr=None)

    print(script)
    # Write the script.
    #import tempfile
    #_, script_file = tempfile.mkstemp(text=True)
    #with open(script_file, "w") as fh: 
    #    fh.write(script)
    #process, queue_id = qad.submit_to_queue(script_file)


if __name__ == "__main__":
    #show_job()
    show_slurm()
