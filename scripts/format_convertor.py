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
import os
import re

from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser, CifWriter
from pymatgen.io.cssrio import Cssr

parser = argparse.ArgumentParser(description='''Convenient file format convertor. 
Author: Shyue Ping Ong
Version: 1.0
Last updated: Oct 26 2011''')
parser.add_argument('input_file', metavar='input file', type=str, nargs = 1, help='input file')
parser.add_argument('output_file', metavar='output file', type=str, nargs = 1, help='output file')
parser.add_argument('-i', '--input-format', dest='input_format', type=str, nargs = 1, choices=['poscar','cif'], default='', help='Input file format.')
parser.add_argument('-o', '--output-format', dest='output_format', type=str, nargs = 1, choices=['poscar', 'cif', 'cssr'], default='', help='Output file format.')

args = parser.parse_args()

def guess_file_format(filename, specified):
    filename, ext = os.path.splitext(filename)
    if specified != '':
        return specified[0]
    else:
        if re.match('POSCAR|CONTCAR', filename):
            return 'poscar'
        else:
            return ext.lower().strip(".")

input_format = guess_file_format(args.input_file[0], args.input_format)
output_format = guess_file_format(args.output_file[0], args.output_format)

try:
    if input_format == 'poscar':
        p = Poscar.from_file(args.input_file[0])
        s = p.struct
    else:
        r = CifParser(args.input_file[0])
        s = r.get_structures()[0]
    
    
    if output_format == 'cif':
        w = CifWriter(s)
    elif output_format == 'poscar':
        w = Poscar(s)
    elif output_format == 'cssr':
        w = Cssr(s)
    
    w.write_file(args.output_file[0])
except:
    print "Error converting file. Are they in the right format?"
    
