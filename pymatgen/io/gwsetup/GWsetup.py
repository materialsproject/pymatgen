#!/usr/bin/env python

"""
Script to write GW Input for VASP and ABINIT.
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

from pymatgen.io.gwsetup.GWdatastructures import GWSpecs

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


if __name__ == "__main__":
    spec_in = GWSpecs()
    spec_in.update_interactive()
    spec_in.test()
    spec_in.write_to_file('spec.in')
    spec_in.loop_structures('i')