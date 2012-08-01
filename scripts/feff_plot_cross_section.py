#!/usr/bin/env python

from __future__ import division

'''
Created on Nov 8, 2011
'''

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Nov 8, 2011"

import argparse
from collections import OrderedDict

from pymatgen.io.feffio import *
from pymatgen.util.plotting_utils import get_publication_quality_plot
from pymatgen.electronic_structure import plotter

parser = argparse.ArgumentParser(description='''Convenient DOS Plotter for Feff runs.
Author: Shyue Ping Ong, Alan Dozier
Version: 1.0
Last updated: May 11, 2012''')
parser.add_argument('filename', metavar='filename', type=str, nargs=1, help='xmu file to plot')
parser.add_argument('filename1', metavar='filename1', type=str, nargs=1, help='feff.inp filename to import')

plt = get_publication_quality_plot(12, 8)
color_order = ['r', 'b', 'g', 'c', 'k', 'm', 'y']

args = parser.parse_args()
xmu = Xmu(args.filename[0], args.filename1[0])

data=xmu.to_dict

plt.title(data['calc'] +' Feff9.1 Calculation for ' + data['atom'] + ' in ' + data['formula'] + ' unit cell')
plt.xlabel('Energies (eV)')
plt.ylabel('Absorbtion Cross-section')

x=data['energies']
y=data['scross']
tle='Single ' + data['atom'] + ' ' + data['edge'] + ' edge'
plt.plot(x,y,color_order[1 % 7], label=tle)

y=data['across']
tle=data['atom'] + ' '  + data['edge'] + ' edge in ' + data['formula']
plt.plot(x,y,color_order[2 % 7], label=tle)

plt.legend()
leg = plt.gca().get_legend()
ltext = leg.get_texts()  # all the text.Text instance in the legend
plt.setp(ltext, fontsize = 15)
plt.tight_layout()
plt.show()
