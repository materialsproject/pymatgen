#!/usr/bin/env python

"""
A convenience script engine using VaspObjects to do all manner of simple outputs.
Written by Shyue Ping Ong
"""

import argparse
import os
import logging
import re

from pymatgen.util.string_utils import str_aligned
from pymatgen.apps.borg.hive import GaussianToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
import multiprocessing

save_file = "gaussian_analyzer_data.gz"

def get_energies(rootdir, reanalyze, verbose, pretty):
    if verbose:
        FORMAT = "%(relativeCreated)d msecs : %(message)s"
        logging.basicConfig(level=logging.INFO, format=FORMAT)
    drone = GaussianToComputedEntryDrone(inc_structure=True, parameters=['filename'])
    ncpus = multiprocessing.cpu_count()
    logging.info('Detected {} cpus'.format(ncpus))
    queen = BorgQueen(drone, number_of_drones=ncpus)
    if os.path.exists(save_file) and not reanalyze:
        msg = 'Using previously assimilated data from {}. Use -f to force re-analysis'.format(save_file)
        queen.load_data(save_file)
    else:
        queen.parallel_assimilate(rootdir)
        msg = 'Analysis results saved to {} for faster subsequent loading.'.format(save_file)
        queen.save_data(save_file)

    entries = queen.get_data()
    entries = sorted(entries, key=lambda x:x.parameters['filename'])
    all_data = [(e.parameters['filename'].replace("./", ""),
                 re.sub("\s+", "", e.composition.formula),
                 "{}".format(e.parameters['charge']),
                 "{}".format(e.parameters['spin_mult']),
                 "{:.5f}".format(e.energy), "{:.5f}".format(e.energy_per_atom),
                 ) for e in entries]
    headers = ("Directory", "Formula", "Charge", "Spin Mult.", "Energy",
               "E/Atom")
    if pretty:
        from prettytable import PrettyTable
        t = PrettyTable(headers)
        t.set_field_align("Directory", "l")
        map(t.add_row, all_data)
        print t
    else:
        print str_aligned(all_data, headers)
    print msg


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Convenient Gaussian run analyzer which can recursively go into a directory to search results.
    Author: Shyue Ping Ong
    Version: 1.0
    Last updated: Jul 6 2012''')
    parser.add_argument('directories', metavar='dir', default='.', type=str, nargs='*', help='directory to process')
    parser.add_argument('-v', '--verbose', dest="verbose", action='store_const', const=True, help='verbose mode. Provides detailed output on progress.')
    parser.add_argument('-p', '--pretty', dest="pretty", action='store_const', const=True, help='pretty mode. Uses prettytable to format output. Must have prettytable module installed.')
    parser.add_argument('-f', '--force', dest="reanalyze", action='store_const', const=True, help='force reanalysis. Typically, gaussian_analyzer will just reuse a gaussian_analyzer_data.gz if present. This forces the analyzer to reanalyze.')

    args = parser.parse_args()
    for d in args.directories:
        get_energies(d, args.reanalyze, args.verbose, args.pretty)
