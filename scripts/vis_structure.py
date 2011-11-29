#!/usr/bin/env python

'''
Simple structure visualizer script.
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 29, 2011"

import argparse

from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser
from pymatgen.vis.structure_vtk import StructureVis

parser = argparse.ArgumentParser(description='''Convenient structure file viewer.  Currently only cif and POSCAR formats supported.
Author: Shyue Ping Ong
Version: 1.0
Last updated: Nov 29 2011''')
parser.add_argument('input_file', metavar='input file', type=str, nargs = 1, help='input file')

parser.add_argument('-f', '--format', dest='format', type=str, nargs = 1, choices=['cif','poscar'], default='poscar', help='Format of input file. Defaults to POSCAR')
parser.add_argument('-e', '--exclude_bonding', dest='exclude_bonding', type=str, nargs = 1, help='List of elements to exclude from bonding analysis. E.g., Li,Na')

args = parser.parse_args()
excluded_bonding_elements = args.exclude_bonding[0].split(',') if args.exclude_bonding else []

if args.format == 'poscar':
    p = Poscar.from_file(args.input_file[0])
    s = p.struct
    vis = StructureVis(excluded_bonding_elements=excluded_bonding_elements)
    vis.set_structure(s)
    vis.show()
else:
    r = CifParser(args.input_file[0])
    s = r.get_structures(False)[0]
    vis = StructureVis(excluded_bonding_elements=excluded_bonding_elements)
    vis.set_structure(s)
    vis.show()
    
    