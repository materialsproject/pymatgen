#!/usr/bin/env python

'''
Created on Aug 12, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Aug 12, 2012"

import argparse
import itertools

from pymatgen.io.smartio import read_structure
from pymatgen.analysis.structure_fitter import StructureFitter

parser = argparse.ArgumentParser(description="""
Convenient structure comparison.""", epilog="""
Author: {}
Version: {}
Last updated: {}""".format(__author__, __version__, __date__))
parser.add_argument("structure_files", metavar="structure_files", type=str,
                    nargs='+',
                    help="Names of files to compare")

args = parser.parse_args()
structure_files = args.structure_files
try:
    structures = [(filename, read_structure(filename))
                  for filename in structure_files]
except Exception as ex:
    print "Error reading files. Are they in the right format?"
    print str(ex)

try:
    for (s1, s2) in itertools.combinations(structures, 2):
        fitter = StructureFitter(s1[1], s2[1])
        print "{} ({}) fits {} ({}) = {}".format(s1[0], s1[1].composition,
                                                 s2[0], s2[1].composition,
                                                 fitter.fit_found)
except Exception as ex:
    print "Error fitting."
    print str(ex)
