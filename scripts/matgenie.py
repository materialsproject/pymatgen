#!/usr/bin/env python

"""
A master convenience script with many tools for vasp and structure analysis.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Dec 1, 2012"

import argparse
import os
import re
import logging
import multiprocessing
import sys
import datetime

from collections import OrderedDict

from pymatgen.io.vaspio import Outcar, Vasprun, Chgcar
from pymatgen.util.string_utils import str_aligned
from pymatgen.apps.borg.hive import SimpleVaspToComputedEntryDrone, \
    VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser, CifWriter
from pymatgen.io.vaspio_set import MPVaspInputSet, MITVaspInputSet
from pymatgen.io.smartio import read_structure, write_structure
from pymatgen.io.cssrio import Cssr
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.alchemy.materials import TransformedStructure

save_file = "vasp_data.gz"


def get_energies(rootdir, reanalyze, verbose, pretty, detailed, sort):
    if verbose:
        FORMAT = "%(relativeCreated)d msecs : %(message)s"
        logging.basicConfig(level=logging.INFO, format=FORMAT)

    if not detailed:
        drone = SimpleVaspToComputedEntryDrone(inc_structure=True)
    else:
        drone = VaspToComputedEntryDrone(inc_structure=True,
                                         data=["filename",
                                               "initial_structure"])

    ncpus = multiprocessing.cpu_count()
    logging.info("Detected {} cpus".format(ncpus))
    queen = BorgQueen(drone, number_of_drones=ncpus)
    if os.path.exists(save_file) and not reanalyze:
        msg = "Using previously assimilated data from {}.".format(save_file) \
            + " Use -f to force re-analysis."
        queen.load_data(save_file)
    else:
        if ncpus > 1:
            queen.parallel_assimilate(rootdir)
        else:
            queen.serial_assimilate(rootdir)
        msg = "Analysis results saved to {} for faster ".format(save_file) + \
              "subsequent loading."
        queen.save_data(save_file)

    entries = queen.get_data()
    if sort == "energy_per_atom":
        entries = sorted(entries, key=lambda x: x.energy_per_atom)
    elif sort == "filename":
        entries = sorted(entries, key=lambda x: x.data["filename"])

    all_data = []
    for e in entries:
        if not detailed:
            delta_vol = "{:.2f}".format(e.data["delta_volume"] * 100)
        else:
            delta_vol = e.structure.volume / \
                e.data["initial_structure"].volume - 1
            delta_vol = "{:.2f}".format(delta_vol * 100)
        all_data.append((e.data["filename"].replace("./", ""),
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
                    mags = [m["tot"] for m in mags]
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


def plot_dos(args):
    v = Vasprun(args.filename[0])
    dos = v.complete_dos

    all_dos = OrderedDict()
    all_dos["Total"] = dos

    structure = v.final_structure

    if args.site:
        for i in xrange(len(structure)):
            site = structure[i]
            all_dos["Site " + str(i) + " " + site.specie.symbol] = \
                dos.get_site_dos(site)

    if args.element:
        syms = [tok.strip() for tok in args.element[0].split(",")]
        all_dos = {}
        for el, dos in dos.get_element_dos().items():
            if el.symbol in syms:
                all_dos[el] = dos
    if args.orbital:
        all_dos = dos.get_spd_dos()

    plotter = DosPlotter()
    plotter.add_dos_dict(all_dos)
    if args.file:
        plotter.get_plot().savefig(args.file[0])
    else:
        plotter.show()


def plot_chgint(args):
    chgcar = Chgcar.from_file(args.filename[0])
    s = chgcar.structure

    if args.inds:
        atom_ind = map(int, args.inds[0].split(","))
    else:
        finder = SymmetryFinder(s, symprec=0.1)
        sites = [sites[0] for sites in
                 finder.get_symmetrized_structure().equivalent_sites]
        atom_ind = [s.sites.index(site) for site in sites]

    from pymatgen.util.plotting_utils import get_publication_quality_plot
    plt = get_publication_quality_plot(12, 8)
    for i in atom_ind:
        d = chgcar.get_integrated_diff(i, args.radius, 30)
        plt.plot(d[:, 0], d[:, 1],
                 label="Atom {} - {}".format(i, s[i].species_string))
    plt.legend(loc="upper left")
    plt.xlabel("Radius (A)")
    plt.ylabel("Integrated charge (e)")
    plt.tight_layout()
    plt.show()


def parse_vasp(args):

    default_energies = not (args.get_energies or args.ion_list)

    if args.get_energies or default_energies:
        for d in args.directories:
            get_energies(d, args.reanalyze, args.verbose, args.pretty,
                         args.detailed, args.sort[0])
    if args.ion_list:
        if args.ion_list[0] == "All":
            ion_list = None
        else:
            (start, end) = map(int, re.split("-", args.ion_list[0]))
            ion_list = range(start, end + 1)
        for d in args.directories:
            get_magnetizations(d, ion_list)


def convert_fmt(args):
    iformat = args.input_format[0]
    oformat = args.output_format[0]
    filename = args.input_filename[0]
    out_filename = args.output_filename[0]

    try:
        if iformat == "smart":
            structure = read_structure(filename)
        if iformat == "POSCAR":
            p = Poscar.from_file(filename)
            structure = p.structure
        elif iformat == "CIF":
            r = CifParser(filename)
            structure = r.get_structures()[0]
        elif iformat == "CSSR":
            structure = Cssr.from_file(filename).structure

        if oformat == "smart":
            write_structure(structure, out_filename)
        elif oformat == "POSCAR":
            p = Poscar(structure)
            p.write_file(out_filename)
        elif oformat == "CIF":
            w = CifWriter(structure)
            w.write_file(out_filename)
        elif oformat == "CSSR":
            c = Cssr(structure)
            c.write_file(out_filename)
        elif oformat == "VASP":
            input_set = MPVaspInputSet()
            ts = TransformedStructure(
                structure, [],
                history=[{"source": "file",
                          "datetime": str(datetime.datetime.now()),
                          "original_file": open(filename).read()}])
            ts.write_vasp_input(input_set, output_dir=out_filename)
        elif oformat == "MITVASP":
            input_set = MITVaspInputSet()
            ts = TransformedStructure(
                structure, [],
                history=[{"source": "file",
                          "datetime": str(datetime.datetime.now()),
                          "original_file": open(filename).read()}])
            ts.write_vasp_input(input_set, output_dir=out_filename)

    except Exception as ex:
        print "Error converting file. Are they in the right format?"
        print str(ex)


def parse_symmetry(args):

    tolerance = float(args.tolerance[0])

    for filename in args.filenames:
        s = read_structure(filename)
        if args.spacegroup:
            finder = SymmetryFinder(s, tolerance)
            dataset = finder.get_symmetry_dataset()
            print filename
            print "Spacegroup  : {}".format(dataset["international"])
            print "Int number  : {}".format(dataset["number"])
            print "Hall symbol : {}".format(dataset["hall"])
            print


def parse_view(args):
    from pymatgen.vis.structure_vtk import StructureVis
    excluded_bonding_elements = args.exclude_bonding[0].split(",") \
        if args.exclude_bonding else []
    s = read_structure(args.filename[0])
    vis = StructureVis(excluded_bonding_elements=excluded_bonding_elements)
    vis.set_structure(s)
    vis.show()


def compare_structures(args):
    filenames = args.filenames
    if len(filenames) < 2:
        print "You need more than one structure to compare!"
        sys.exit(-1)
    try:
        structures = map(read_structure, filenames)
    except Exception as ex:
        print "Error converting file. Are they in the right format?"
        print str(ex)
        sys.exit(-1)

    from pymatgen.analysis.structure_matcher import StructureMatcher, \
        ElementComparator
    m = StructureMatcher() if args.oxi \
        else StructureMatcher(comparator=ElementComparator())
    for i, grp in enumerate(m.group_structures(structures)):
        print "Group {}: ".format(i)
        for s in grp:
            print "- {} ({})".format(filenames[structures.index(s)],
                                     s.formula)
        print


def generate_files(args):
    from pymatgen.io.vaspio.vasp_input import Potcar
    if args.symbols:
        try:
            p = Potcar(args.symbols, functional=args.functional)
            p.write_file("POTCAR")
        except Exception as ex:
            print "An error has occurred: {}".format(str(ex))

    else:
        print "No valid options selected."


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    matgenie is a convenient script that uses pymatgen to perform many
    analyses, plotting and format conversions. This script works based on
    several sub-commands with their own options. To see the options for the
    sub-commands, type "matgenie.py sub-command -h".""",
                                     epilog="""
    Author: Shyue Ping Ong
    Version: {}
    Last updated: {}""".format(__version__, __date__))

    subparsers = parser.add_subparsers()

    parser_vasp = subparsers.add_parser("analyze", help="Vasp run analysis.")
    parser_vasp.add_argument("directories", metavar="dir", default=".",
                             type=str, nargs="*",
                             help="directory to process (default to .)")
    parser_vasp.add_argument("-e", "--energies", dest="get_energies",
                             action="store_true", help="Print energies")
    parser_vasp.add_argument("-m", "--mag", dest="ion_list", type=str, nargs=1,
                             help="Print magmoms. ION LIST can be a range "
                             "(e.g., 1-2) or the string 'All' for all ions.")
    parser_vasp.add_argument("-f", "--force", dest="reanalyze",
                             action="store_true",
                             help="Force reanalysis. Typically, vasp_analyzer"
                             " will just reuse a vasp_analyzer_data.gz if "
                             "present. This forces the analyzer to reanalyze "
                             "the data.")
    parser_vasp.add_argument("-v", "--verbose", dest="verbose",
                             action="store_true",
                             help="verbose mode. Provides detailed output on "
                             "progress.")
    parser_vasp.add_argument("-p", "--pretty", dest="pretty",
                             action="store_const",
                             const=True, help="pretty mode. Uses prettytable "
                             "to format output. Must have prettytable module "
                             "installed.")
    parser_vasp.add_argument("-d", "--detailed", dest="detailed",
                             action="store_true",
                             help="Detailed mode. Parses vasprun.xml instead "
                             "of separate vasp input. Slower.")
    parser_vasp.add_argument("-s", "--sort", dest="sort", type=str, nargs=1,
                             default=["energy_per_atom"],
                             help="Sort criteria. Defaults to energy / atom.")
    parser_vasp.set_defaults(func=parse_vasp)

    parser_plot = subparsers.add_parser("plotdos", help="Plotting for dos.")
    parser_plot.add_argument("filename", metavar="filename", type=str, nargs=1,
                             help="vasprun.xml file to plot")
    parser_plot.add_argument("-s", "--site", dest="site", action="store_const",
                             const=True, help="Plot site projected DOS")
    parser_plot.add_argument("-e", "--element", dest="element", type=str,
                             nargs=1,
                             help="List of elements to plot as comma-separated"
                             " values e.g., Fe,Mn")
    parser_plot.add_argument("-o", "--orbital", dest="orbital",
                             action="store_const", const=True,
                             help="Plot orbital projected DOS")
    parser_plot.add_argument("-f", "--file", dest="file", type=str, nargs=1,
                             help="Save to file.")
    parser_plot.set_defaults(func=plot_dos)

    parser_plotchg = subparsers.add_parser("plotchgint",
                                           help="Plotting for the charge "
                                                "integration.")
    parser_plotchg.add_argument("filename", metavar="filename", type=str,
                                nargs=1, help="CHGCAR file to plot")
    parser_plotchg.add_argument("-i", "--indices", dest="inds", type=str,
                                nargs=1,
                                help="Comma-separated list of indices to plot"
                                     ", e.g., 1,2,3,4. If not provided, "
                                     "the code will plot the chgint for all "
                                     "symmetrically distinct atoms detected.")
    parser_plotchg.add_argument("-r", "--radius", dest="radius", type=float,
                                default=3,
                                help="Radius of integration.")
    parser_plotchg.set_defaults(func=plot_chgint)

    parser_convert = subparsers.add_parser("convert",
                                           help="File format conversion tools."
                                           )
    parser_convert.add_argument("input_filename", metavar="input_filename",
                                type=str, nargs=1, help="Input filename.")
    parser_convert.add_argument("output_filename", metavar="output_filename",
                                type=str, nargs=1,
                                help="Output filename (for POSCAR/CIF/CSSR "
                                "output) / dirname (VASP output)")
    parser_convert.add_argument("-i", "--input", dest="input_format", type=str,
                                nargs=1,
                                choices=["POSCAR", "CIF", "CSSR", "smart"],
                                default=["smart"],
                                help="Input file format. By default, smart is "
                                "selected, which guesses the format from the "
                                "filename. Other formats can be enforced as "
                                "needed.")
    parser_convert.add_argument("-o", "--output", dest="output_format",
                                type=str, nargs=1,
                                choices=["POSCAR", "CIF", "CSSR", "VASP",
                                         "MITVASP"
                                         "smart"],
                                default=["smart"],
                                help="Output file format. By default, smart is"
                                " selected, which guesses the format from the "
                                "filename. Other formats can be enforced as "
                                "needed. VASP is a special output form, which "
                                "outputs a set of VASP input files to a "
                                "directory. MITVASP uses the MIT input set "
                                "instead of the default Materials project "
                                "input set.")
    parser_convert.set_defaults(func=convert_fmt)

    parser_symm = subparsers.add_parser("symm", help="Symmetry tools.")
    parser_symm.add_argument("filenames", metavar="filenames", type=str,
                             nargs="+",
                             help="Filenames to determine symmetry.")
    parser_symm.add_argument("-t", "--tolerance", dest="tolerance", type=float,
                             nargs=1, default=[0.1],
                             help="Tolerance for symmetry determination")
    parser_symm.add_argument("-s", "--spacegroup", dest="spacegroup",
                             action="store_true",
                             help="Determine symmetry")
    parser_symm.set_defaults(func=parse_symmetry)

    parser_view = subparsers.add_parser("view", help="Visualize structures")
    parser_view.add_argument("filename", metavar="filename", type=str,
                             nargs=1, help="Filename")
    parser_view.add_argument("-e", "--exclude_bonding", dest="exclude_bonding",
                             type=str, nargs=1,
                             help="List of elements to exclude from bonding "
                             "analysis. E.g., Li,Na")
    parser_view.set_defaults(func=parse_view)

    parser_cmp = subparsers.add_parser("compare", help="Compare structures")
    parser_cmp.add_argument("filenames", metavar="filenames", type=str,
                            nargs="*", help="List of filenames to compare.")
    parser_cmp.add_argument("-o", "--oxi", dest="oxi",
                            action="store_true",
                            help="Oxi mode means that different oxidation "
                                 "states will not match to each other, i.e.,"
                                 " Fe2+ amd Fe3+ will be treated as "
                                 "different species for the purposes of "
                                 "matching.")
    parser_cmp.set_defaults(func=compare_structures)

    parser_generate = subparsers.add_parser("generate",
                                            help="Generate input files")
    parser_generate.add_argument("-f", "--functional", dest="functional",
                                 type=str,
                                 choices=["LDA", "PBE", "PW91", "LDA_US"],
                                 default="PBE",
                                 help="Functional to use. Unless otherwise "
                                      "stated (e.g., US), "
                                      "refers to PAW psuedopotential.")
    parser_generate.add_argument("-p", "--potcar", dest="symbols",
                                 type=str, nargs="+",
                                 help="List of POTCAR symbols. Use -f to set "
                                      "functional. Defaults to PBE.")
    parser_generate.set_defaults(func=generate_files)


    args = parser.parse_args()
    args.func(args)
