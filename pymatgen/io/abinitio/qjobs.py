#!/usr/bin/env python
import sys

from pymatgen.io.abinitio.qadapters import SlurmJob, PbsProJob

if __name__ == "__main__":
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
