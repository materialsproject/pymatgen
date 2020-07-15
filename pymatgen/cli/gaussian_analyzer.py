#!/usr/bin/env python
#  coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
A convenience script engine to read Gaussian output in a directory tree.
"""


import argparse
import os
import logging
import re
import multiprocessing
from tabulate import tabulate
from pymatgen.apps.borg.hive import GaussianToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen


save_file = "gau_data.gz"


def get_energies(rootdir, reanalyze, verbose):
    """
    :param rootdir:
    :param reanalyze:
    :param verbose:
    :return:
    """
    if verbose:
        FORMAT = "%(relativeCreated)d msecs : %(message)s"
        logging.basicConfig(level=logging.INFO, format=FORMAT)
    drone = GaussianToComputedEntryDrone(inc_structure=True,
                                         parameters=['filename'])
    ncpus = multiprocessing.cpu_count()
    logging.info('Detected {} cpus'.format(ncpus))
    queen = BorgQueen(drone, number_of_drones=ncpus)
    if os.path.exists(save_file) and not reanalyze:
        msg = 'Using previously assimilated data from {}.'.format(save_file) + \
              ' Use -f to force re-analysis.'
        queen.load_data(save_file)
    else:
        queen.parallel_assimilate(rootdir)
        msg = 'Results saved to {} for faster reloading.'.format(save_file)
        queen.save_data(save_file)

    entries = queen.get_data()
    entries = sorted(entries, key=lambda x: x.parameters['filename'])
    all_data = [(e.parameters['filename'].replace("./", ""),
                 re.sub(r"\s+", "", e.composition.formula),
                 "{}".format(e.parameters['charge']),
                 "{}".format(e.parameters['spin_mult']),
                 "{:.5f}".format(e.energy), "{:.5f}".format(e.energy_per_atom),
                 ) for e in entries]
    headers = ("Directory", "Formula", "Charge", "Spin Mult.", "Energy",
               "E/Atom")
    print(tabulate(all_data, headers=headers))
    print("")
    print(msg)


def main():
    """
    Main function
    """
    desc = '''
    Convenient Gaussian run analyzer which can recursively go into a directory
    to search results.
    Author: Shyue Ping Ong
    Version: 1.0
    Last updated: Jul 6 2012'''

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('directories', metavar='dir', default='.', type=str,
                        nargs='*', help='directory to process')
    parser.add_argument('-v', '--verbose', dest="verbose",
                        action='store_const', const=True,
                        help='Verbose mode. Provides detailed output on progress.')
    parser.add_argument('-f', '--force', dest="reanalyze",
                        action='store_const',
                        const=True,
                        help='Force reanalysis, instead of reusing gaussian_analyzer_data.gz.')

    args = parser.parse_args()
    for d in args.directories:
        get_energies(d, args.reanalyze, args.verbose)


if __name__ == "__main__":
    main()
