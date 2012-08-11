#!/usr/bin/env python

from __future__ import division

"""
Created on Nov 14, 2011
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 14, 2011"

import argparse

from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser, CifWriter
from pymatgen.io.vaspio_set import MaterialsProjectVaspInputSet
from pymatgen.io.smartio import read_structure, write_structure
from pymatgen.io.cssrio import Cssr

parser = argparse.ArgumentParser(description="""
Convenient file format convertor.
Author: Shyue Ping Ong
Version: 2.0
Last updated: Oct 26 2011""")
parser.add_argument("input_file", metavar="input file", type=str, nargs=1,
                    help="input file")
parser.add_argument("output_file", metavar="output file", type=str, nargs=1,
                    help="Output file (for POSCAR or CIF output)/dirname (VASP output)")

parser.add_argument("-i", "--input", dest="input_format", type=str,
                    nargs=1,
                    choices=["POSCAR", "CIF", "CSSR", "smart"],
                    default=["smart"],
                    help="Input file format. By default, smart is selected, which guesses the format from the filename. Other formats can be enforced as needed.")

parser.add_argument("-o", "--output", dest="output_format", type=str, nargs=1,
                    choices=["POSCAR", "CIF", "CSSR", "VASP", "smart"],
                    default=["smart"],
                    help="Output file format. By default, smart is selected, which guesses the format from the filename. Other formats can be enforced as needed.")

args = parser.parse_args()
iformat = args.input_format[0]
oformat = args.output_format[0]
filename = args.input_file[0]
out_filename = args.output_file[0]

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
        input_set = MaterialsProjectVaspInputSet()
        input_set.write_input(structure, output_dir=out_filename)

except Exception as ex:
    print "Error converting file. Are they in the right format?"
    print str(ex)
