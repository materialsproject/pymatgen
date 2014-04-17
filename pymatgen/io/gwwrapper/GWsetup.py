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

from pymatgen.io.gwwrapper.GWdatastructures import get_spec

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


if __name__ == "__main__":
    import this
    spec_in = get_spec('GW')
    try:
        spec_in.read_from_file('spec.in')
    except (IOError, OSError):
        try:
            spec_in.read_from_file('$HOME/.abinit/abipy/spec.in')
        except (IOError, OSError):
            pass
        pass
    spec_in.update_interactive()
    spec_in.test()
    spec_in.write_to_file('spec.in')
    spec_in.loop_structures('i')