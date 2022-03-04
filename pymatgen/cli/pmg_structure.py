#!/usr/bin/env python
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Implementation for `pmg structure` CLI.
"""

import sys

from tabulate import tabulate

from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "4.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "Aug 13 2016"


def convert_fmt(args):
    """
    Convert files from one format to another

    Args:
        args (dict): Args from argparse.
    """
    if len(args.filenames) != 2:
        print("File format conversion takes in only two filenames.")
    s = Structure.from_file(args.filenames[0], primitive="prim" in args.filenames[1].lower())
    s.to(filename=args.filenames[1])


def analyze_symmetry(args):
    """
    Analyze symmetry of structures in files.

    Args:
        args (dict): Args from argparse.
    """
    tolerance = args.symmetry
    t = []
    for filename in args.filenames:
        s = Structure.from_file(filename, primitive=False)
        finder = SpacegroupAnalyzer(s, tolerance)
        dataset = finder.get_symmetry_dataset()
        t.append([filename, dataset["international"], dataset["number"], dataset["hall"]])
    print(tabulate(t, headers=["Filename", "Int Symbol", "Int number", "Hall"]))


def analyze_localenv(args):
    """
    Analyze local env of structures in files.

    Args:
        args (dict): Args for argparse.
    """
    bonds = {}
    for bond in args.localenv:
        toks = bond.split("=")
        species = toks[0].split("-")
        bonds[(species[0], species[1])] = float(toks[1])
    for filename in args.filenames:
        print(f"Analyzing {filename}...")
        data = []
        s = Structure.from_file(filename)
        for i, site in enumerate(s):
            for species, dist in bonds.items():
                if species[0] in [sp.symbol for sp in site.species.keys()]:
                    dists = [
                        d
                        for nn, d in s.get_neighbors(site, dist)
                        if species[1] in [sp.symbol for sp in nn.species.keys()]
                    ]
                    dists = ", ".join([f"{d:.3f}" for d in sorted(dists)])
                    data.append([i, species[0], species[1], dists])
        print(tabulate(data, headers=["#", "Center", "Ligand", "Dists"]))


def compare_structures(args):
    """
    Compare structures in files for similarity using structure matcher.

    Args:
        args (dict): Args from argparse.
    """
    filenames = args.filenames
    if len(filenames) < 2:
        print("You need more than one structure to compare!")
        sys.exit(-1)
    try:
        structures = [Structure.from_file(fn) for fn in filenames]
    except Exception as ex:
        print("Error converting file. Are they in the right format?")
        print(str(ex))
        sys.exit(-1)

    m = StructureMatcher() if args.group == "species" else StructureMatcher(comparator=ElementComparator())
    for i, grp in enumerate(m.group_structures(structures)):
        print(f"Group {i}: ")
        for s in grp:
            print(f"- {filenames[structures.index(s)]} ({s.formula})")
        print()


def analyze_structures(args):
    """
    Master function to handle which operation to perform.

    Args:
        args (dict): Args from argparse.
    """
    if args.convert:
        convert_fmt(args)
    elif args.symmetry:
        analyze_symmetry(args)
    elif args.group:
        compare_structures(args)
    elif args.localenv:
        analyze_localenv(args)
