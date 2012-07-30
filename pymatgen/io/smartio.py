#!/usr/bin/env python

"""
This class implements smart io classes that performs intelligent io based on
file extensions.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 29, 2012"

from pymatgen.io.vaspio import Vasprun, Poscar
from pymatgen.io.cifio import CifParser


def read_structure(filename):
    """
    Reads a structure based on file extension. For example, anything ending in
    a "cif" is assumed to be a Crystallographic Information Format file.
    """
    lower_filename = filename.lower()
    if lower_filename.endswith(".cif"):
        parser = CifParser(filename)
        return parser.get_structures(True)[0]
    elif lower_filename.startswith("poscar") \
        or lower_filename.startswith("contcar"):
        return Poscar(filename).structure
    elif lower_filename.startswith("vasprun"):
        return Vasprun(filename).final_structure
    raise ValueError("Unrecognized file extension!")
