#!/usr/bin/env python

"""
A convenience script engine to read vasp output in a directory tree.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 9, 2012"


import argparse
import os
import re
import logging

from pymatgen.io.vaspio import Outcar
from pymatgen.util.string_utils import str_aligned
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone, \
    VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
import multiprocessing

save_file = 'vasp_data.gz'


def get_energies(rootdir, reanalyze, verbose, pretty, detailed, sort):
    if verbose:
        FORMAT = "%(relativeCreated)d msecs : %(message)s"
        logging.basicConfig(level=logging.INFO, format=FORMAT)

    if not detailed:
        drone = SimpleVaspToComputedEntryDrone(inc_structure=True)
    else:
        drone = VaspToComputedEntryDrone(inc_structure=True,
                                         data=['filename',
                                               'initial_structure'])

    ncpus = multiprocessing.cpu_count()
    logging.info('Detected {} cpus'.format(ncpus))
    queen = BorgQueen(drone, number_of_drones=ncpus)
    if os.path.exists(save_file) and not reanalyze:
        msg = 'Using previously assimilated data from {}.'.format(save_file) \
            + ' Use -f to force re-analysis.'
        queen.load_data(save_file)
    else:
        queen.parallel_assimilate(rootdir)
        msg = 'Analysis results saved to {} for faster '.format(save_file) + \
              'subsequent loading.'
        queen.save_data(save_file)

    entries = queen.get_data()
    if sort == "energy_per_atom":
        entries = sorted(entries, key=lambda x: x.energy_per_atom)
    elif sort == "filename":
        entries = sorted(entries, key=lambda x: x.data['filename'])

    all_data = []
    for e in entries:
        if not detailed:
            delta_vol = "{:.2f}".format(e.data['delta_volume'] * 100)
        else:
            delta_vol = e.structure.volume / \
                e.data['initial_structure'].volume - 1
            delta_vol = "{:.2f}".format(delta_vol * 100)
        all_data.append((e.data['filename'].replace("./", ""),
                     re.sub("\s+", "", e.composition.formula),
                     "{:.5f}".format(e.energy),
                     "{:.5f}".format(e.energy_per_atom),
                     delta_vol))
    if len(all_data) > 0:
        headers = ("Directory", "Formula", "Energy", "E/Atom", "% vol chg")
        if pretty:
            from prettytable import PrettyTable
            t = PrettyTable(headers)
            t.set_field_align("Directory", "l")
            map(t.add_row, all_data)
            print t
        else:
            print str_aligned(all_data, headers)
        print msg
    else:
        print "No valid vasp run found."


def get_magnetizations(mydir, ion_list):
    data = []
    max_row = 0
    for (parent, subdirs, files) in os.walk(mydir):
        for f in files:
            if re.match("OUTCAR*", f):
                try:
                    row = []
                    fullpath = os.path.join(parent, f)
                    outcar = Outcar(fullpath)
                    mags = outcar.magnetization
                    mags = [m['tot'] for m in mags]
                    all_ions = xrange(len(mags))
                    row.append(fullpath.lstrip("./"))
                    if ion_list:
                        all_ions = ion_list
                    for ion in all_ions:
                        row.append(str(mags[ion]))
                    data.append(row)
                    if len(all_ions) > max_row:
                        max_row = len(all_ions)
                except:
                    pass

    for d in data:
        if len(d) < max_row + 1:
            d.extend([""] * (max_row + 1 - len(d)))
    headers = ["Filename"]
    for i in xrange(max_row):
        headers.append(str(i))
    print str_aligned(data, headers)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''
    Convenient vasp run analyzer which can recursively go into a directory to
    search results.
    Author: Shyue Ping Ong
    Version: 1.0
    Last updated: Jul 14 2012''')
    parser.add_argument('directories', metavar='dir', default='.',
                        type=str, nargs='*',
                        help='directory to process (default to .)')
    parser.add_argument('-e', '--energies', dest='get_energies',
                        action='store_true', help='Print energies')
    parser.add_argument('-m', '--mag', dest="ion_list", type=str, nargs=1,
                        help='Print magmoms. ION LIST can be a range '
                        '(e.g., 1-2) or the string "All" for all ions.')
    parser.add_argument('-f', '--force', dest="reanalyze", action='store_true',
                        help='Force reanalysis. Typically, vasp_analyzer'
                        ' will just reuse a vasp_analyzer_data.gz if present. '
                        'This forces the analyzer to reanalyze the data.')
    parser.add_argument('-v', '--verbose', dest="verbose", action='store_true',
                        help='verbose mode. Provides detailed output on '
                        'progress.')
    parser.add_argument('-p', '--pretty', dest="pretty", action='store_const',
                        const=True, help='pretty mode. Uses prettytable to '
                        'format output. Must have prettytable module '
                        'installed.')
    parser.add_argument('-d', '--detailed', dest="detailed",
                        action='store_true',
                        help='Detailed mode. Parses vasprun.xml instead of '
                        'separate vasp input. Slower.')
    parser.add_argument('-s', '--sort', dest="sort", type=str, nargs=1,
                        default=['energy_per_atom'],
                        help='Sort criteria. Defaults to energy / atom.')

    args = parser.parse_args()

    default_energies = not (args.get_energies or args.ion_list)

    if args.get_energies or default_energies:
        for d in args.directories:
            get_energies(d, args.reanalyze, args.verbose, args.pretty,
                         args.detailed, args.sort[0])
    if args.ion_list:
        ion_list = list()
        if args.ion_list[0] == "All":
            ion_list = None
        else:
            (start, end) = map(int, re.split("-", args.ion_list[0]))
            ion_list = range(start, end + 1)
        for d in args.directories:
            get_magnetizations(d, ion_list)
