#!/usr/bin/env python

'''
Created on Nov 12, 2011
'''

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 12, 2011"

import itertools
import argparse

from pymatgen.io.vaspio import Incar
from pymatgen.util.string_utils import str_aligned


parser = argparse.ArgumentParser(description='''Convenient INCAR diff. 
Author: Shyue Ping Ong
Version: 1.0
Last updated: Oct 26 2011''')
parser.add_argument('incar_file', metavar='filename', type=str, nargs = 2, help='files to process')

args = parser.parse_args()

filepath1 = args.incar_file[0]
filepath2 = args.incar_file[1]
incar1 = Incar.from_file(filepath1)
incar2 = Incar.from_file(filepath2)

print type(incar1['EDIFF'])
print type(incar2['EDIFF'])

def format_lists(v):
    if isinstance(v, (tuple, list)):
        return " ".join(["%d*%.2f" % (len(tuple(group)), i) for (i,group) in itertools.groupby(v)])
    return v

d = incar1.diff(incar2)
output = [['SAME PARAMS','', '']]
output.append(['---------------','', ''])
output.extend([(k,format_lists(v),format_lists(v)) for k,v in d['Same'].items() if k != "SYSTEM"])
output.append(['','', ''])
output.append(['DIFFERENT PARAM','', ''])
output.append(['---------------','', ''])
output.extend([(k,format_lists(v['INCAR1']),format_lists(v['INCAR2'])) for k, v in d['Different'].items() if k != "SYSTEM"])
print str_aligned(output, ['', filepath1, filepath2])