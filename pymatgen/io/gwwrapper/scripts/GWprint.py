# coding: utf-8

from __future__ import unicode_literals, division, print_function

#!/usr/bin/env python

"""
Script to print GW results for VASP and ABINIT.
"""

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "May 2014"

import os
import os.path

from pymatgen.io.gwwrapper.datastructures import get_spec
from pymatgen.io.gwwrapper.helpers import load_ps

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


if __name__ == "__main__":
    load_ps()
    counter = 0
    spec_in = get_spec('GW')
    spec_in.read_from_file('spec.in')
    spec_in.loop_structures('w')
