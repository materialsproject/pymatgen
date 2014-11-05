#!/usr/bin/env python
import sys

from pymatgen.io.abinitio.qadapters import SlurmJob

if __name__ == "__main__":
    job_id = sys.argv[1]
    job = SlurmJob(job_id)
    print("start_time: ", job.get_start_time())
    print("get_info:" ,job.get_info())
    print("get_stats:" ,job.get_stats())

    #if job.is_completed:
    #if job.timeout
    #if job.is_running:
    #if job.is_failed:
    #if job.has_node_failures:
    #  nodes = job.save_nodes_in_blacklist()
