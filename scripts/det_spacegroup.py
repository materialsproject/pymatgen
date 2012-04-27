#!/usr/bin/env python

'''
Simple structure visualizer script.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 29, 2011"

import argparse

from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser
from pymatgen.symmetry.spglib_adaptor import SymmetryFinder

parser = argparse.ArgumentParser(description='''Convenient structure spacegroup determination.  Currently only cif and POSCAR formats supported.
Author: Shyue Ping Ong
Version: 1.0
Last updated: Apr 24 2012''')
parser.add_argument('input_file', metavar='input file', type=str, nargs=1, help='input file')

parser.add_argument('-f', '--format', dest='format', type=str, nargs=1, choices=['cif', 'poscar'], default='poscar', help='Format of input file. Defaults to POSCAR')
parser.add_argument('-t', '--tolerance', dest='tolerance', type=float, nargs=1, help='Tolerance for symmetry determination')

args = parser.parse_args()
tolerance = float(args.tolerance[0]) if args.tolerance else 0.1

file_format = args.format

if args.input_file[0].lower().endswith(".cif"):
    file_format = "cif"
elif args.input_file[0].upper().startswith("POSCAR"):
    file_format = "poscar"

s = None

if file_format == 'poscar':
    p = Poscar.from_file(args.input_file[0])
    s = p.struct
else:
    r = CifParser(args.input_file[0])
    s = r.get_structures(False)[0]

if s:
    finder = SymmetryFinder(s, tolerance)
    dataset = finder.get_symmetry_dataset()
    print "Spacegroup  : {}".format(dataset['international'])
    print "Int number  : {}".format(dataset['number'])
    print "Hall symbol : {}".format(dataset["hall"])

