#!/usr/bin/env python
# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import sys

from tabulate import tabulate

from pymatgen import Structure

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher, \
        ElementComparator

"""
A master convenience script with many tools for vasp and structure analysis.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "4.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "Aug 13 2016"


def convert_fmt(args):
    if len(args.filenames) != 2:
        print("File format conversion takes in only two filenames.")
    s = Structure.from_file(args.filenames[0],
                            primitive="prim" in args.filenames[1].lower())
    s.to(filename=args.filenames[1])


def analyze_symmetry(args):
    tolerance = args.symmetry
    t = []
    for filename in args.filenames:
        s = Structure.from_file(filename, primitive=False)
        finder = SpacegroupAnalyzer(s, tolerance)
        dataset = finder.get_symmetry_dataset()
        t.append([filename, dataset["international"], dataset["number"],
                  dataset["hall"]])
    print(tabulate(t, headers=["Filename", "Int Symbol", "Int number", "Hall"]))


def analyze_localenv(args):
    bonds = {}
    for bond in args.localenv:
        toks = bond.split("=")
        species = toks[0].split("-")
        bonds[(species[0], species[1])] = float(toks[1])
    for filename in args.filenames:
        print("Analyzing %s..." % filename)
        data = []
        s = Structure.from_file(filename)
        for i, site in enumerate(s):
            for species, dist in bonds.items():
                if species[0] in [sp.symbol
                                  for sp in site.species_and_occu.keys()]:
                    dists = [d for nn, d in s.get_neighbors(site, dist)
                             if species[1] in
                             [sp.symbol for sp in nn.species_and_occu.keys()]]
                    dists = ", ".join(["%.3f" % d for d in sorted(dists)])
                    data.append([i, species[0], species[1], dists])
        print(tabulate(data, headers=["#", "Center", "Ligand", "Dists"]))


def compare_structures(args):
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

    m = StructureMatcher() if args.group == "species" \
        else StructureMatcher(comparator=ElementComparator())
    for i, grp in enumerate(m.group_structures(structures)):
        print("Group {}: ".format(i))
        for s in grp:
            print("- {} ({})".format(filenames[structures.index(s)],
                                     s.formula))
        print()


def analyze_structures(args):
    if args.convert:
        convert_fmt(args)
    elif args.symmetry:
        analyze_symmetry(args)
    elif args.group:
        compare_structures(args)
    elif args.localenv:
        analyze_localenv(args)