#!/usr/bin/env python
# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import argparse
import sys
import itertools

from tabulate import tabulate, tabulate_formats
from pymatgen import Structure, SETTINGS
from pymatgen.io.vasp import Incar, Potcar

from pymatgen.cli.pmg_analyze import analyze
from pymatgen.cli.pmg_config import configure_pmg
from pymatgen.cli.pmg_potcar import generate_potcar
from pymatgen.cli.pmg_plot import plot
from pymatgen.cli.pmg_structure import analyze_structures
from pymatgen.cli.pmg_query import do_query

"""
A master convenience script with many tools for vasp and structure analysis.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "5.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "Sep 19 2016"


def parse_view(args):
    from pymatgen.vis.structure_vtk import StructureVis
    excluded_bonding_elements = args.exclude_bonding[0].split(",") \
        if args.exclude_bonding else []
    s = Structure.from_file(args.filename[0])
    vis = StructureVis(excluded_bonding_elements=excluded_bonding_elements)
    vis.set_structure(s)
    vis.show()


def diff_incar(args):
    filepath1 = args.incars[0]
    filepath2 = args.incars[1]
    incar1 = Incar.from_file(filepath1)
    incar2 = Incar.from_file(filepath2)

    def format_lists(v):
        if isinstance(v, (tuple, list)):
            return " ".join(["%d*%.2f" % (len(tuple(group)), i)
                             for (i, group) in itertools.groupby(v)])
        return v

    d = incar1.diff(incar2)
    output = [['SAME PARAMS', '', ''], ['---------------', '', ''],
              ['', '', ''], ['DIFFERENT PARAMS', '', ''],
              ['----------------', '', '']]
    output.extend([(k, format_lists(d['Same'][k]), format_lists(d['Same'][k]))
                   for k in sorted(d['Same'].keys()) if k != "SYSTEM"])
    output.extend([(k, format_lists(d['Different'][k]['INCAR1']),
                    format_lists(d['Different'][k]['INCAR2']))
                   for k in sorted(d['Different'].keys()) if k != "SYSTEM"])
    print(tabulate(output, headers=['', filepath1, filepath2]))


def main():
    parser = argparse.ArgumentParser(description="""
    pmg is a convenient script that uses pymatgen to perform many
    analyses, plotting and format conversions. This script works based on
    several sub-commands with their own options. To see the options for the
    sub-commands, type "pmg sub-command -h".""",
                                     epilog="""
    Author: Shyue Ping Ong
    Version: {}
    Last updated: {}""".format(__version__, __date__))

    subparsers = parser.add_subparsers()

    parser_config = subparsers.add_parser(
        "config", help="Tools for configuring pymatgen, e.g., "
                       "potcar setup, modifying .pmgrc.yaml "
                       "configuration file.")
    groups = parser_config.add_mutually_exclusive_group(required=True)
    groups.add_argument("-p", "--potcar", dest="potcar_dirs",
                        metavar="dir_name",
                        nargs=2,
                        help="Initial directory where downloaded VASP "
                             "POTCARs are extracted to, and the "
                             "output directory where the reorganized "
                             "potcars will be stored. The input "
                             "directory should be "
                             "the parent directory that contains the "
                             "POT_GGA_PAW_PBE or potpaw_PBE type "
                             "subdirectories.")
    groups.add_argument("-i", "--install", dest="install",
                        metavar="package_name",
                        choices=["enumlib", "bader"],
                        help="Install various optional command line "
                             "tools needed for full functionality.")

    groups.add_argument("-a", "--add", dest="var_spec", nargs="+",
                        help="Variables to add in the form of space "
                             "separated key value pairs. E.g., "
                             "PMG_VASP_PSP_DIR ~/psps")
    parser_config.set_defaults(func=configure_pmg)

    parser_analyze = subparsers.add_parser(
        "analyze", help="Vasp calculation analysis tools.")
    parser_analyze.add_argument("directories", metavar="dir", default=".",
                                type=str, nargs="*",
                                help="directory to process (default to .)")
    parser_analyze.add_argument("-e", "--energies", dest="get_energies",
                                action="store_true", help="Print energies")
    parser_analyze.add_argument(
        "-m", "--mag", dest="ion_list", type=str, nargs=1,
        help="Print magmoms. ION LIST can be a range "
             "(e.g., 1-2) or the string 'All' for all ions.")
    parser_analyze.add_argument(
        "-r", "--reanalyze", dest="reanalyze", action="store_true",
        help="Force reanalysis. Typically, vasp_analyzer"
             " will just reuse a vasp_analyzer_data.gz if "
             "present. This forces the analyzer to reanalyze "
             "the data.")
    parser_analyze.add_argument(
        "-f", "--format", dest="format", choices=tabulate_formats,
        default="simple",
        help="Format for table. Supports all options in tabulate package.")
    parser_analyze.add_argument(
        "-v", "--verbose", dest="verbose", action="store_true",
        help="Verbose mode. Provides detailed output on progress.")
    parser_analyze.add_argument(
        "-d", "--detailed", dest="detailed", action="store_true",
        help="Detailed, but slower mode. Parses vasprun.xml instead of "
             "separate vasp outputs.")
    parser_analyze.add_argument(
        "-s", "--sort", dest="sort", choices=["energy_per_atom", "filename"],
        default="energy_per_atom",
        help="Sort criteria. Defaults to energy / atom.")
    parser_analyze.set_defaults(func=analyze)

    parser_query = subparsers.add_parser(
        "query",
        help="Search for structures and data from the Materials Project.")
    parser_query.add_argument(
        "criteria", metavar="criteria",
        help="Search criteria. Supported formats in formulas, chemical "
             "systems, Materials Project ids, etc.")
    group = parser_query.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-s", "--structure", dest="structure", metavar="format",
        choices=["poscar", "cif", "cssr"], type=str.lower,
        help="Get structures from Materials Project and write them to a "
             "specified format.")
    group.add_argument(
        "-e", "--entries", dest="entries", metavar="filename",
        help="Get entries from Materials Project and write them to "
             "serialization file. JSON and YAML supported.")
    group.add_argument(
        "-d", "--data", dest="data", metavar="fields", nargs="*",
        help="Print a summary of entries in the Materials Project satisfying "
             "the criteria. Supply field names to include additional data. "
             "By default, the Materials Project id, formula, spacegroup, "
             "energy per atom, energy above hull are shown.")
    parser_query.set_defaults(func=do_query)

    parser_plot = subparsers.add_parser("plot", help="Plotting tool for "
                                                     "DOS, CHGCAR, XRD, etc.")
    group = parser_plot.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--dos', dest="dos_file", metavar="vasprun.xml",
                       help="Plot DOS from a vasprun.xml")
    group.add_argument('-c', '--chgint', dest="chgcar_file", metavar="CHGCAR",
                       help="Generate charge integration plots from any "
                            "CHGCAR")
    group.add_argument('-x', '--xrd', dest="xrd_structure_file",
                       metavar="structure_file",
                       help="Generate XRD plots from any supported structure "
                            "file, e.g., CIF, POSCAR, vasprun.xml, etc.")

    parser_plot.add_argument("-s", "--site", dest="site", action="store_const",
                             const=True, help="Plot site projected DOS")
    parser_plot.add_argument("-e", "--element", dest="element", type=str,
                             nargs=1,
                             help="List of elements to plot as comma-separated"
                                  " values e.g., Fe,Mn")
    parser_plot.add_argument("-o", "--orbital", dest="orbital",
                             action="store_const", const=True,
                             help="Plot orbital projected DOS")

    parser_plot.add_argument("-i", "--indices", dest="inds", type=str,
                             nargs=1,
                             help="Comma-separated list of indices to plot "
                                  "charge integration, e.g., 1,2,3,4. If not "
                                  "provided, the code will plot the chgint "
                                  "for all symmetrically distinct atoms "
                                  "detected.")
    parser_plot.add_argument("-r", "--radius", dest="radius", type=float,
                             default=3,
                             help="Radius of integration for charge "
                                  "integration plot.")
    parser_plot.add_argument("--out_file", dest="out_file", type=str,
                             help="Save plot to file instead of displaying.")
    parser_plot.set_defaults(func=plot)

    parser_structure = subparsers.add_parser(
        "structure",
        help="Structure conversion and analysis tools.")

    parser_structure.add_argument(
        "-f", "--filenames", dest="filenames",
        metavar="filename", nargs="+",
        help="List of structure files.")

    groups = parser_structure.add_mutually_exclusive_group(required=True)
    groups.add_argument("-c", "--convert", dest="convert", action="store_true",
                        help="Convert from structure file 1 to structure "
                             "file 2. Format determined from filename. "
                             "Supported formats include POSCAR/CONTCAR, "
                             "CIF, CSSR, etc. If the keyword'prim' is within "
                             "the filename, the code will automatically attempt "
                             "to find a primitive cell.")
    groups.add_argument("-s", "--symmetry", dest="symmetry",
                        metavar="tolerance", type=float,
                        help="Determine the spacegroup using the "
                             "specified tolerance. 0.1 is usually a good "
                             "value for DFT calculations.")
    groups.add_argument("-g", "--group", dest="group",
                        choices=["element", "species"],
                        metavar="mode",
                        help="Compare a set of structures for similarity. "
                             "Element mode does not compare oxidation states. "
                             "Species mode will take into account oxidations "
                             "states.")
    groups.add_argument(
        "-l", "--localenv", dest="localenv", nargs="+",
        help="Local environment analysis. Provide bonds in the format of"
             "Center Species-Ligand Species=max_dist, e.g., H-O=0.5.")

    parser_structure.set_defaults(func=analyze_structures)

    parser_view = subparsers.add_parser("view", help="Visualize structures")
    parser_view.add_argument("filename", metavar="filename", type=str,
                             nargs=1, help="Filename")
    parser_view.add_argument("-e", "--exclude_bonding", dest="exclude_bonding",
                             type=str, nargs=1,
                             help="List of elements to exclude from bonding "
                                  "analysis. E.g., Li,Na")
    parser_view.set_defaults(func=parse_view)

    parser_diff = subparsers.add_parser(
        "diff", help="Diffing tool. For now, only INCAR supported.")
    parser_diff.add_argument("-i", "--incar", dest="incars",
                             metavar="INCAR", required=True,
                             nargs=2, help="List of INCARs to compare.")
    parser_diff.set_defaults(func=diff_incar)

    parser_potcar = subparsers.add_parser("potcar",
                                          help="Generate POTCARs")
    parser_potcar.add_argument("-f", "--functional", dest="functional",
                               type=str,
                               choices=sorted(Potcar.FUNCTIONAL_CHOICES),
                               default=SETTINGS.get("PMG_DEFAULT_FUNCTIONAL",
                                                    "PBE"),
                               help="Functional to use. Unless otherwise "
                                    "stated (e.g., US), "
                                    "refers to PAW psuedopotential.")
    group = parser_potcar.add_mutually_exclusive_group(required=True)

    group.add_argument("-s", "--symbols", dest="symbols",
                       type=str, nargs="+",
                       help="List of POTCAR symbols. Use -f to set "
                            "functional. Defaults to PBE.")
    group.add_argument("-r", "--recursive", dest="recursive",
                       type=str, nargs="+",
                       help="Dirname to find and generate from POTCAR.spec.")
    parser_potcar.set_defaults(func=generate_potcar)


    try:
        import argcomplete
        argcomplete.autocomplete(parser)
    except ImportError:
        # argcomplete not present.
        pass

    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)


if __name__ == "__main__":
    main()
