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

import re
import os

from pymatgen.io.vaspio import Vasprun, Poscar
from pymatgen.io.cifio import CifParser, CifWriter
from pymatgen.io.cssrio import Cssr


def read_structure(filename):
    """
    Reads a structure based on file extension. For example, anything ending in
    a "cif" is assumed to be a Crystallographic Information Format file.

    Args:
        filename:
            A filename to read from.

    Returns:
        A Structure object.
    """
    lower_filename = os.path.basename(filename).lower()
    if re.search("\.cif", lower_filename):
        parser = CifParser(filename)
        return parser.get_structures(True)[0]
    elif lower_filename.startswith("poscar") \
            or lower_filename.startswith("contcar"):
        return Poscar.from_file(filename, False).structure
    elif re.search("vasprun", lower_filename) \
            and re.search("xml", lower_filename):
        return Vasprun(filename).final_structure
    elif re.search("\.cssr", lower_filename):
        cssr = Cssr.from_file(filename)
        return cssr.structure

    raise ValueError("Unrecognized file extension!")


def write_structure(structure, filename):
    """
    Write a structure to a file based on file extension. For example, anything
    ending in a "cif" is assumed to be a Crystallographic Information Format
    file.

    Args:
        structure:
            Structure to write
        filename:
            A filename to write to.
    """
    lower_filename = os.path.basename(filename).lower()
    if re.search("\.cif", lower_filename):
        writer = CifWriter(structure)
    elif lower_filename.startswith("poscar") \
            or lower_filename.startswith("contcar"):
        writer = Poscar(structure)
    elif re.search("\.cssr", lower_filename):
        writer = Cssr(structure)
    else:
        raise ValueError("Unrecognized file extension!")

    writer.write_file(filename)
