#!/usr/bin/env python
"""This scripts drives the execution of the pseudo dojo tests."""
from __future__ import division, print_function

import sys

from argparse import ArgumentParser
from pprint import pprint

from pymatgen.io.abinitio.task import RunMode
from pymatgen.io.abinitio.pseudo_dojo import Dojo

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

################################################################################

def main():

    parser = ArgumentParser()

    parser.add_argument('-m', '--max_ncpus', type=int, default=1,
                        help="Maximum number of CPUs used by the DOJO.")

    parser.add_argument('-n', '--mpi-ncpus', type=int, default=1,
                        help="Number of MPI Cpus per run).")

    parser.add_argument('-l', '--max-level', type=int, default=0, help="Maximum DOJO level).")

    #parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                     help='verbose, can be supplied multiple times to increase verbosity')

    parser.add_argument('pseudos', nargs='+', help='List of pseudopotential files')

    options = parser.parse_args()

    max_ncpus = options.max_ncpus
    mpi_ncpus = options.mpi_ncpus

    if mpi_ncpus > max_ncpus:
        raise ValueError("mpi_cpus %(mpi_ncpus)s > max_ncpus %(max_ncpus)s" % locals())

    #runmode = RunMode.sequential()
    #runmode = RunMode.load_user_configuration()
    runmode = RunMode.mpi_parallel(mpi_ncpus=mpi_ncpus)
    pprint(runmode)

    dojo = Dojo(runmode=runmode, max_ncpus=max_ncpus, max_level=options.max_level)

    for pseudo in options.pseudos:
        dojo.challenge_pseudo(pseudo)

################################################################################

if __name__ == "__main__":
    sys.exit(main())
