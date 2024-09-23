"""Implementation for `pmg analyze` CLI."""

from __future__ import annotations

import logging
import multiprocessing
import os
import re

from tabulate import tabulate

from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone, VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.io.vasp import Outcar

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "4.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "Aug 13 2016"

SAVE_FILE = "vasp_data.gz"


def get_energies(rootdir, reanalyze, verbose, quick, sort, fmt):
    """Get energies of all vaspruns in directory (nested).

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
        log_fmt = "%(relativeCreated)d msecs : %(message)s"
        logging.basicConfig(level=logging.INFO, format=log_fmt)

    if quick:
        drone = SimpleVaspToComputedEntryDrone(inc_structure=True)
    else:
        drone = VaspToComputedEntryDrone(inc_structure=True, data=["filename", "initial_structure"])

    n_cpus = multiprocessing.cpu_count()
    logging.info(f"Detected {n_cpus} cpus")
    queen = BorgQueen(drone, number_of_drones=n_cpus)
    if os.path.isfile(SAVE_FILE) and not reanalyze:
        msg = f"Using previously assimilated data from {SAVE_FILE}. Use -r to force re-analysis."
        queen.load_data(SAVE_FILE)
    else:
        if n_cpus > 1:
            queen.parallel_assimilate(rootdir)
        else:
            queen.serial_assimilate(rootdir)
        msg = f"Analysis results saved to {SAVE_FILE} for faster subsequent loading."
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
            delta_vol = e.structure.volume / e.data["initial_structure"].volume - 1
            delta_vol = f"{delta_vol * 100:.2f}"
        all_data.append(
            (
                e.data["filename"].replace("./", ""),
                re.sub(r"\s+", "", e.formula),
                f"{e.energy:.5f}",
                f"{e.energy_per_atom:.5f}",
                delta_vol,
            )
        )
    if len(all_data) > 0:
        headers = ("Directory", "Formula", "Energy", "E/Atom", "% vol chg")
        print(tabulate(all_data, headers=headers, tablefmt=fmt))
        print()
        print(msg)
    else:
        print("No valid vasp run found.")
        os.unlink(SAVE_FILE)
    return 0


def get_magnetizations(dirc: str, ion_list: list[int]):
    """Get magnetization info from OUTCARs.

    Args:
        dirc (str): Directory name
        ion_list (list[int]): List of ions to obtain magnetization information for.

    Returns:
        int: 0 if successful.
    """
    data = []
    max_row = 0
    for parent, _subdirs, files in os.walk(dirc):
        for file in files:
            if re.match(r"OUTCAR*", file):
                try:
                    row = []
                    fullpath = os.path.join(parent, file)
                    outcar = Outcar(fullpath)
                    mags = outcar.magnetization
                    _mags: list = [m["tot"] for m in mags]
                    all_ions = list(range(len(_mags)))
                    row.append(fullpath.lstrip("./"))
                    if ion_list:
                        all_ions = ion_list
                    for ion in all_ions:
                        row.append(str(_mags[ion]))
                    data.append(row)
                    max_row = max(len(all_ions), max_row)
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
    """Master function controlling which analysis to call.

    Args:
        args (dict): args from argparse.
    """
    default_energies = not (args.get_energies or args.ion_list)

    if args.get_energies or default_energies:
        for folder in args.directories:
            return get_energies(folder, args.reanalyze, args.verbose, args.quick, args.sort, args.format)
    if args.ion_list:
        if args.ion_list[0] == "All":
            ion_list = None
        else:
            start, end = (int(i) for i in re.split(r"-", args.ion_list[0]))
            ion_list = list(range(start, end + 1))
        for folder in args.directories:
            return get_magnetizations(folder, ion_list)

    return -1
