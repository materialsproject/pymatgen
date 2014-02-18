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

from pymatgen.io.gwsetup.GWsetup import GWSpecs

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


if __name__ == "__main__":
    spec = GWSpecs()
    spec.read_from_file('spec.in')
    print 'Found setup for ', spec.data['code']
    spec.loop_structures('o')