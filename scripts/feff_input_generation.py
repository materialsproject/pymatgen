#!/usr/bin/env python
#  coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
Example script that generates FEFF input files from a cif file
Remove comment # on write line to actually write files to disk
"""

from __future__ import division

__author__ = "Alan Dozier"
__credits__= "Anubhav Jain, Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0.2"
__maintainer__ = "Alan Dozier"
__email__ = "adozier@uky.edu"
__date__ = "April 7, 2013"

import argparse

from pymatgen.io.feff.sets import *
from pymatgen.io.vasp import *
from pymatgen.io.cif import CifParser

parser = argparse.ArgumentParser(description='''
Example script to generate FEFF input files from a cif file
Author: Alan Dozier
Version: 1.0
Last updated: August, 2012''')


parser.add_argument('cif_file', metavar='cif_file', type=str, nargs=1, help='cif_file to use')
parser.add_argument('central_atom', metavar='central_atom', type=str, nargs=1, help='symbol of absorbing atom')
parser.add_argument('calc_type', metavar='calc_type', type=str, nargs=1, help='type of calc, currently XANES or EXAFS')

args = parser.parse_args()
cif_file = args.cif_file[0]
central_atom = args.central_atom[0]
calc_type = args.calc_type[0]

r = CifParser(cif_file)
structure = r.get_structures()[0]
x = FeffInputSet("MaterialsProject")

source = os.path.basename(cif_file)
comment = 'From cif file'

header = FeffInputSet.get_header(x, structure, source, comment)
print "\n\nHEADER\n"
print header

tags = FeffInputSet.get_feff_tags(x, calc_type)
print "\n\nPARAMETERS\n"
print tags

POT = FeffInputSet.get_feff_pot(x, structure, central_atom)
print "\n\nPOTENTIALS\n"
print POT

ATOMS = FeffInputSet.get_feff_atoms(x, structure, central_atom)
print"\n\nATOMS\n"
print ATOMS

#x.write_input(structure, calc_type, source, central_atom, comment, "./feffinput")