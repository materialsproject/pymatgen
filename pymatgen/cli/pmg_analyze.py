# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
Implementation for `pmg analyze` CLI.
"""

import os
import re
import logging
import multiprocessing

from tabulate import tabulate

from pymatgen.io.vasp import Outcar
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone, \
    VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "4.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "Aug 13 2016"

SAVE_FILE = "vasp_data.gz"


def get_energies(rootdir, reanalyze, verbose, quick, sort, fmt):
    """
    Get energies of all vaspruns in directory (nested).
    Args:
        rootdir (str): Root directory.
        reanalyze (bool): Whether to ignore saved results and reanalyze
        verbose (bool): Verbose mode or not.
        quick (bool): Whether to perform a quick analysis (using OSZICAR instead
            of vasprun.xml
        sort (bool): Whether to sort the results in ascending order.
        fmt (str): tablefmt passed to tabulate.
    """
    if verbose:
        logformat = "%(relativeCreated)d msecs : %(message)s"
        logging.basicConfig(level=logging.INFO, format=logformat)

    if quick:
        drone = SimpleVaspToComputedEntryDrone(inc_structure=True)
    else:
        drone = VaspToComputedEntryDrone(inc_structure=True,
                                         data=["filename",
                                               "initial_structure"])

    ncpus = multiprocessing.cpu_count()
    logging.info("Detected {} cpus".format(ncpus))
    queen = BorgQueen(drone, number_of_drones=ncpus)
    if os.path.exists(SAVE_FILE) and not reanalyze:
        msg = "Using previously assimilated data from {}.".format(SAVE_FILE) \
              + " Use -r to force re-analysis."
        queen.load_data(SAVE_FILE)
    else:
        if ncpus > 1:
            queen.parallel_assimilate(rootdir)
        else:
            queen.serial_assimilate(rootdir)
        msg = "Analysis results saved to {} for faster ".format(SAVE_FILE) + \
              "subsequent loading."
        queen.save_data(SAVE_FILE)

    entries = queen.get_data()
    if sort == "energy_per_atom":
        entries = sorted(entries, key=lambda x: x.energy_per_atom)
    elif sort == "filename":
        entries = sorted(entries, key=lambda x: x.data["filename"])

    all_data = []
    for e in entries:
        if quick:
            delta_vol = "NA"
        else:
            delta_vol = e.structure.volume / \
                        e.data["initial_structure"].volume - 1
            delta_vol = "{:.2f}".format(delta_vol * 100)
        all_data.append((e.data["filename"].replace("./", ""),
                         re.sub(r"\s+", "", e.composition.formula),
                         "{:.5f}".format(e.energy),
                         "{:.5f}".format(e.energy_per_atom),
                         delta_vol))
    if len(all_data) > 0:
        headers = ("Directory", "Formula", "Energy", "E/Atom", "% vol chg")
        print(tabulate(all_data, headers=headers, tablefmt=fmt))
        print("")
        print(msg)
    else:
        print("No valid vasp run found.")
        os.unlink(SAVE_FILE)
    return 0


def get_magnetizations(mydir, ion_list):
    """
    Get magnetization info from OUTCARs.

    Args:
        mydir (str): Directory name
        ion_list (List): List of ions to obtain magnetization information for.

    Returns:

    """
    data = []
    max_row = 0
    for (parent, subdirs, files) in os.walk(mydir):
        for f in files:
            if re.match(r"OUTCAR*", f):
                try:
                    row = []
                    fullpath = os.path.join(parent, f)
                    outcar = Outcar(fullpath)
                    mags = outcar.magnetization
                    mags = [m["tot"] for m in mags]
                    all_ions = list(range(len(mags)))
                    row.append(fullpath.lstrip("./"))
                    if ion_list:
                        all_ions = ion_list
                    for ion in all_ions:
                        row.append(str(mags[ion]))
                    data.append(row)
                    if len(all_ions) > max_row:
                        max_row = len(all_ions)
                except Exception:
                    pass

    for d in data:
        if len(d) < max_row + 1:
            d.extend([""] * (max_row + 1 - len(d)))
    headers = ["Filename"]
    for i in range(max_row):
        headers.append(str(i))
    print(tabulate(data, headers))
    return 0


def analyze(args):
    """
    Master function controlling which analysis to call.

    Args:
        args (dict): args from argparse.
    """
    default_energies = not (args.get_energies or args.ion_list)

    if args.get_energies or default_energies:
        for d in args.directories:
            return get_energies(d, args.reanalyze, args.verbose,
                                args.quick, args.sort, args.format)
    if args.ion_list:
        if args.ion_list[0] == "All":
            ion_list = None
        else:
            (start, end) = [int(i) for i in re.split(r"-", args.ion_list[0])]
            ion_list = list(range(start, end + 1))
        for d in args.directories:
            return get_magnetizations(d, ion_list)

    return -1
