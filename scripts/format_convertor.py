#!/usr/bin/env python

from __future__ import division

'''
Created on Nov 14, 2011
'''

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 14, 2011"

import argparse

from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser, CifWriter

parser = argparse.ArgumentParser(description='''Convenient file format convertor. 
Author: Shyue Ping Ong
Version: 1.0
Last updated: Oct 26 2011''')
parser.add_argument('input_file', metavar='input file', type=str, nargs = 1, help='input file')
parser.add_argument('output_file', metavar='output file', type=str, nargs = 1, help='output file')

parser.add_argument('-c', '--conversion', dest='conversion', type=str, nargs = 1, choices=['poscar2cif','cif2poscar'], default='poscar2cif', help='Format conversion desired. ')

args = parser.parse_args()
try:
    if args.conversion[0] == 'poscar2cif':
        p = Poscar.from_file(args.input_file[0])
        w = CifWriter(p.struct)
        w.write_file(args.output_file[0])
    else:
        r = CifParser(args.input_file[0])
        p = Poscar(r.get_structures()[0])
        p.write_file(args.output_file[0])
except Exception as ex:
    print "Error converting file. Are they in the right format?"
    print str(ex)
    