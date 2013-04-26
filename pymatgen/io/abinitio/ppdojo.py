#!/usr/bin/env python
from __future__ import division, print_function

import sys

from argparse import ArgumentParser 

from pymatgen.io.abinitio.task import RunMode
from pymatgen.io.abinitio.pseudo_dojo import Dojo

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

##########################################################################################

def main():

    parser = ArgumentParser()

    parser.add_argument('-j', '--py-nthreads', type=int, default=1,
                        help="The number of threads used (run PY_NTHREADS calculations simultaneously).")

    parser.add_argument('-l', '--max-level', type=int, default=0,
                        help="Maximum DOJO level).")

    #parser.add_argument('-v', '--verbose', default=0, action='count', # -vv --> verbose=2
    #                     help='verbose, can be supplied multiple times to increase verbosity')  

    parser.add_argument('pseudos', nargs='+', help='List of pseudopotential files')

    options = parser.parse_args()

    runmode = RunMode.sequential()
    dojo = Dojo(max_level=options.max_level, max_ncpus=options.py_nthreads)
    print(dojo)

    for pseudo in options.pseudos:
        dojo_masters = dojo.challenge_pseudo(pseudo, runmode)

##########################################################################################

if __name__ == "__main__":
    sys.exit(main())
