#!/usr/bin/env python
#  coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
Created on Nov 12, 2011
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 12, 2011"

import itertools
import argparse

from tabulate import tabulate

from pymatgen.io.vasp import Incar


parser = argparse.ArgumentParser(description='''Convenient INCAR diff.
Author: Shyue Ping Ong
Version: 1.0
Last updated: Oct 26 2011''')
parser.add_argument('incar_file', metavar='filename', type=str, nargs=2,
                    help='files to process')

args = parser.parse_args()

filepath1 = args.incar_file[0]
filepath2 = args.incar_file[1]
incar1 = Incar.from_file(filepath1)
incar2 = Incar.from_file(filepath2)


def format_lists(v):
    if isinstance(v, (tuple, list)):
        return " ".join(["%d*%.2f" % (len(tuple(group)), i)
                         for (i, group) in itertools.groupby(v)])
    return v

d = incar1.diff(incar2)
output = [['SAME PARAMS', '', '']]
output.append(['---------------', '', ''])
output.extend([(k, format_lists(d['Same'][k]), format_lists(d['Same'][k]))
               for k in sorted(d['Same'].keys()) if k != "SYSTEM"])
output.append(['', '', ''])
output.append(['DIFFERENT PARAMS', '', ''])
output.append(['----------------', '', ''])
output.extend([(k, format_lists(d['Different'][k]['INCAR1']),
                format_lists(d['Different'][k]['INCAR2']))
               for k in sorted(d['Different'].keys()) if k != "SYSTEM"])
print(tabulate(output, headers=['', filepath1, filepath2]))
