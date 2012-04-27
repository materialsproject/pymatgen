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


def detect_ncpus():
    """
    Detects the number of CPUs on a system. Cribbed from pp.
    """
    # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
    else: # OSX:
        return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    # Windows:
    if os.environ.has_key("NUMBER_OF_PROCESSORS"):
        ncpus = int(os.environ["NUMBER_OF_PROCESSORS"]);
        if ncpus > 0:
            return ncpus
    return 1 # Default

def get_energies(rootdir, reanalyze):
    FORMAT = "%(relativeCreated)d msecs : %(message)s"
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    drone = VaspToComputedEntryDrone(data=['filename'])
    ncpus = detect_ncpus()
    logging.info('Detect {} cpus'.format(ncpus))
    queen = BorgQueen(drone, number_of_drones=ncpus)
    if os.path.exists('vasp_analyzer_data.gz') and not reanalyze:
        logging.info('Using previously assimilated data file vasp_analyzer_data.gz. Use -f to force re-analysis')
        queen.load_data('vasp_analyzer_data.gz')
    else:
        queen.parallel_assimilate(rootdir)
        queen.save_data('vasp_analyzer_data.gz')
    entries = queen.get_assimilated_data()
    entries = sorted(entries, key=lambda x:x.data['filename'])
    all_data = [(e.data['filename'], e.energy, e.energy_per_atom) for e in entries]
    print str_aligned(all_data, ("Directory", "Energy", "Energy/Atom"))

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
