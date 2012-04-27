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


def get_energies(rootdir, reanalyze):
    FORMAT = "%(relativeCreated)d msecs : %(message)s"
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    drone = VaspToComputedEntryDrone(inc_structure=True, data=['filename'])
    ncpus = multiprocessing.cpu_count()
    logging.info('Detected {} cpus'.format(ncpus))
    queen = BorgQueen(drone, number_of_drones=ncpus)
    if os.path.exists('vasp_analyzer_data.gz') and not reanalyze:
        logging.info('Using previously assimilated data file vasp_analyzer_data.gz. Use -f to force re-analysis')
        queen.load_data('vasp_analyzer_data.gz')
    else:
        queen.parallel_assimilate(rootdir)
        queen.save_data('vasp_analyzer_data.gz')
    entries = queen.get_assimilated_data()
    entries = sorted(entries, key=lambda x:x.data['filename'])
    all_data = [(e.data['filename'].replace("./", ""), e.composition.formula, "{:.5f}".format(e.energy), "{:.5f}".format(e.energy_per_atom), "{:.2f}".format(e.structure.volume)) for e in entries]
    print str_aligned(all_data, ("Directory", "Formula", "Energy", "E/Atom", "Vol"))

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

parser = argparse.ArgumentParser(description='''Convenient vasp run analyzer which can recursively go into a directory to search results.
Author: Shyue Ping Ong
Version: 1.0
Last updated: Nov 15 2011''')
parser.add_argument('directories', metavar='dir', default='.', type=str, nargs='*', help='directory to process (default to .)')
parser.add_argument('-e', '--energies', dest='get_energies', action='store_const', const=True, help='print energies')
parser.add_argument('-m', '--mag', dest="ion_list", type=str, nargs=1, help='print magmoms. ION LIST can be a range (e.g., 1-2) or the string "All" for all ions.')
parser.add_argument('-f', '--force', dest="reanalyze", action='store_const', const=True, help='force reanalysis. Typically, vasp_analyzer will just reuse a vasp_analyzer_data.gz if present. This forces the analyzer to reanalyzer.')

args = parser.parse_args()
if args.get_energies:
    for d in args.directories:
        get_energies(d, args.reanalyze)
if args.ion_list:
    ion_list = list()
    (start, end) = map(int, re.split("-", args.ion_list[0]))
    ion_list = range(start, end + 1)
    get_magnetizations(args.directories, ion_list)
