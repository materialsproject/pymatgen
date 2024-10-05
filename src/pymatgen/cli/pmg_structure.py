#!/usr/bin/env python

"""Implementation for `pmg structure` CLI."""

from __future__ import annotations

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
    """Convert files from one format to another.

    Args:
        args (dict): Args from argparse.
    """
    if len(args.filenames) != 2:
        print("File format conversion takes in only two filenames.")
    struct = Structure.from_file(args.filenames[0], primitive="prim" in args.filenames[1].lower())
    struct.to(filename=args.filenames[1])


def analyze_symmetry(args):
    """Analyze symmetry of structures in files.

    Args:
        args (dict): Args from argparse.
    """
    tolerance = args.symmetry
    table_rows = []
    for filename in args.filenames:
        struct = Structure.from_file(filename, primitive=False)
        finder = SpacegroupAnalyzer(struct, tolerance)
        dataset = finder.get_symmetry_dataset()
        table_rows.append([filename, dataset.international, dataset.number, dataset.hall])
    print(tabulate(table_rows, headers=["Filename", "Int Symbol", "Int number", "Hall"]))


def analyze_localenv(args):
    """Analyze local env of structures in files.

    Args:
        args (dict): Args for argparse.
    """
    bonds = {}
    for bond in args.localenv:
        tokens = bond.split("=")
        species = tokens[0].split("-")
        bonds[species[0], species[1]] = float(tokens[1])
    for filename in args.filenames:
        print(f"Analyzing {filename}...")
        data = []
        struct = Structure.from_file(filename)
        for idx, site in enumerate(struct):
            for species, dist in bonds.items():
                if species[0] in [sp.symbol for sp in site.species]:
                    dists = [
                        nn.nn_distance
                        for nn in struct.get_neighbors(site, dist)
                        if species[1] in [sp.symbol for sp in nn.species]
                    ]
                    dists = ", ".join(f"{d:.3f}" for d in sorted(dists))
                    data.append([idx, species[0], species[1], dists])
        print(tabulate(data, headers=["#", "Center", "Ligand", "Dists"]))


def compare_structures(args):
    """Compare structures in files for similarity using structure matcher.

    Args:
        args (dict): Args from argparse.
    """
    filenames = args.filenames
    if len(filenames) < 2:
        raise SystemExit("You need more than one structure to compare!")
    try:
        structures = [Structure.from_file(fn) for fn in filenames]
    except Exception as exc:
        print("Error converting file. Are they in the right format?")
        raise SystemExit(exc)

    matcher = StructureMatcher() if args.group == "species" else StructureMatcher(comparator=ElementComparator())
    for idx, grp in enumerate(matcher.group_structures(structures)):
        print(f"Group {idx}: ")
        for s in grp:
            print(f"- {filenames[structures.index(s)]} ({s.formula})")
        print()


def analyze_structures(args):
    """Master function to handle which operation to perform.

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
