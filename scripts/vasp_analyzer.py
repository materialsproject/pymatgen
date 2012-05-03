#!/usr/bin/env python

"""
A convenience script engine using VaspObjects to do all manner of simple outputs.
Written by Shyue Ping Ong
"""

import argparse
import os
import re
import logging

from pymatgen.io.vaspio import Outcar
from pymatgen.util.string_utils import str_aligned
from pymatgen.borg.hive import VaspToComputedEntryDrone
from pymatgen.borg.queen import BorgQueen
import multiprocessing


def get_energies(rootdir, reanalyze, verbose, pretty):
    if verbose:
        FORMAT = "%(relativeCreated)d msecs : %(message)s"
        logging.basicConfig(level=logging.INFO, format=FORMAT)
    drone = VaspToComputedEntryDrone(inc_structure=True, data=['filename'])
    ncpus = multiprocessing.cpu_count()
    logging.info('Detected {} cpus'.format(ncpus))
    queen = BorgQueen(drone, number_of_drones=ncpus)
    if os.path.exists('vasp_analyzer_data.gz') and not reanalyze:
        msg = 'Using previously assimilated data from vasp_analyzer_data.gz. Use -f to force re-analysis'
        queen.load_data('vasp_analyzer_data.gz')
    else:
        queen.parallel_assimilate(rootdir)
        msg = 'Analysis results saved to vasp_analyzer_data.gz for faster subsequent loading.'
        queen.save_data('vasp_analyzer_data.gz')

    entries = queen.get_data()
    entries = sorted(entries, key=lambda x:x.data['filename'])
    all_data = [(e.data['filename'].replace("./", ""), e.composition.formula, "{:.5f}".format(e.energy), "{:.5f}".format(e.energy_per_atom), "{:.2f}".format(e.structure.volume)) for e in entries]
    headers = ("Directory", "Formula", "Energy", "E/Atom", "Vol")
    if pretty:
        from prettytable import PrettyTable
        t = PrettyTable(headers)
        t.set_field_align("Directory", "l")
        map(t.add_row, all_data)
        print t
    else:
        print str_aligned(all_data, headers)
    print msg


def get_magnetizations(mydir, ionList):
    print "%10s | %7s" % ("Ion", "Magmoms")
    print "-" * 20
    for (parent, subdirs, files) in os.walk(mydir):
        for f in files:
            if re.match("OUTCAR*", f):
                fullpath = os.path.join(parent, f)
                outcar = Outcar(fullpath)
                mags = outcar.magnetization
                mags = [m['tot'] for m in mags]
                allIons = xrange(len(mags))
                print "%16s" % (fullpath.lstrip("./"))
                if len(ionList) > 0:
                    allIons = ionList
                for ion in allIons:
                    print "%10d | %3.4f" % (ion, mags[ion])
                print "-" * 20


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Convenient vasp run analyzer which can recursively go into a directory to search results.
    Author: Shyue Ping Ong
    Version: 1.0
    Last updated: Nov 15 2011''')
    parser.add_argument('directories', metavar='dir', default='.', type=str, nargs='*', help='directory to process (default to .)')
    parser.add_argument('-e', '--energies', dest='get_energies', action='store_const', const=True, help='print energies')
    parser.add_argument('-m', '--mag', dest="ion_list", type=str, nargs=1, help='print magmoms. ION LIST can be a range (e.g., 1-2) or the string "All" for all ions.')
    parser.add_argument('-f', '--force', dest="reanalyze", action='store_const', const=True, help='force reanalysis. Typically, vasp_analyzer will just reuse a vasp_analyzer_data.gz if present. This forces the analyzer to reanalyzer.')
    parser.add_argument('-v', '--verbose', dest="verbose", action='store_const', const=True, help='verbose mode. Provides detailed output on progress.')
    parser.add_argument('-p', '--pretty', dest="pretty", action='store_const', const=True, help='pretty mode. Uses prettytable to format output. Must have prettytable module installed.')

    args = parser.parse_args()
    if args.get_energies:
        for d in args.directories:
            get_energies(d, args.reanalyze, args.verbose, args.pretty)
    if args.ion_list:
        ion_list = list()
        (start, end) = map(int, re.split("-", args.ion_list[0]))
        ion_list = range(start, end + 1)
        for d in args.directories:
            get_magnetizations(d, ion_list)
