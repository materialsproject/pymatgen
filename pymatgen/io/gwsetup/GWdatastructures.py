#!/usr/bin/env python

"""
Script to test writing GW Input for VASP.
Reads the POSCAR_name in the the current folder and outputs GW input to
subfolders name
"""

from __future__ import division

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "Oct 23, 2013"

import os
import os.path

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class GWConvergenceData():
    def __init__(self, structure, spec):
        self.structure = structure
        self.spec = spec
        self.data = {}
        self.name = structure.composition.reduced_formula

    def read(self):
        if self.spec['code'] == 'ABINIT':
            read_next = True
            n = 3
            while read_next:
                output = os.path.join(self.name,  'work_0', 'task_' + str(n), 'run.abo')
                if os.path.isfile(output):
                    n += 1
                else:
                    read_next = False

    def print_plot_data(self):
        for data in self.data:
            print data['nbands'], data['encuteps'], data['gap']