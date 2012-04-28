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
import re
import json

from pymatgen.io.vaspio import Poscar
from pymatgen.io.cifio import CifParser
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.vis.structure_vtk import StructureVis
from pymatgen.util.io_utils import file_open_zip_aware

parser = argparse.ArgumentParser(description='''Convenient structure file viewer.  Currently only cif and POSCAR formats supported.
Author: Shyue Ping Ong
Version: 1.0
Last updated: Nov 29 2011''')
parser.add_argument('input_file', metavar='input file', type=str, nargs=1, help='input file')

parser.add_argument('-f', '--format', dest='format', type=str, nargs=1, choices=['cif', 'poscar', 'mpjson'], default='poscar', help='Format of input file. Defaults to POSCAR')
parser.add_argument('-e', '--exclude_bonding', dest='exclude_bonding', type=str, nargs=1, help='List of elements to exclude from bonding analysis. E.g., Li,Na')

args = parser.parse_args()
excluded_bonding_elements = args.exclude_bonding[0].split(',') if args.exclude_bonding else []

file_format = args.format
filename = args.input_file[0]

s = None

if filename.endswith(".cif"):
    file_format = "cif"
elif filename.startswith("POSCAR"):
    file_format = "poscar"
elif re.search('\.json', filename):
    file_format = 'mpjson'


if file_format == 'poscar':
    p = Poscar.from_file(filename)
    s = p.struct
elif file_format == 'cif':
    r = CifParser(filename)
    s = r.get_structures(False)[0]
else:
    d = json.load(file_open_zip_aware(filename))
    ts = TransformedStructure.from_dict(d)
    s = ts.final_structure

if s:
    vis = StructureVis(excluded_bonding_elements=excluded_bonding_elements)
    vis.set_structure(s)
    vis.show()
