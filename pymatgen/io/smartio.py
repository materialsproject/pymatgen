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
__email__ = "shyuep@gmail.com"
__date__ = "Jul 29, 2012"

import re
import os
import json

from pymatgen.core.structure import Structure, Molecule
from pymatgen.io.vaspio import Vasprun, Poscar, Chgcar
from pymatgen.io.cifio import CifParser, CifWriter
from pymatgen.io.cssrio import Cssr
from pymatgen.io.xyzio import XYZ
from pymatgen.io.gaussianio import GaussianInput, GaussianOutput
from pymatgen.util.io_utils import zopen
from pymatgen.serializers.json_coders import PMGJSONDecoder, PMGJSONEncoder
from pymatgen.io.babelio import BabelMolAdaptor


def read_structure(filename):
    """
    Reads a structure based on file extension. For example, anything ending in
    a "cif" is assumed to be a Crystallographic Information Format file.
    Supported formats include CIF, POSCAR/CONTCAR, CHGCAR, LOCPOT,
    vasprun.xml, CSSR and pymatgen's JSON serialized structures.

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
    elif lower_filename.startswith("chgcar") \
            or lower_filename.startswith("locpot"):
        return Chgcar.from_file(filename).structure
    elif re.search("vasprun", lower_filename) \
            and re.search("xml", lower_filename):
        return Vasprun(filename).final_structure
    elif re.search("\.cssr", lower_filename):
        cssr = Cssr.from_file(filename)
        return cssr.structure
    elif re.search("\.[mj]son", lower_filename):
        with zopen(lower_filename) as f:
            s = json.load(f, cls=PMGJSONDecoder)
            if type(s) != Structure:
                raise IOError("File does not contain a valid serialized "
                              "structure")
            return s

    raise ValueError("Unrecognized file extension!")


def write_structure(structure, filename):
    """
    Write a structure to a file based on file extension. For example, anything
    ending in a "cif" is assumed to be a Crystallographic Information Format
    file. Supported formats include CIF, POSCAR, CSSR and pymatgen's JSON
    serialized structures.

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
    elif re.search("\.[mj]son", lower_filename):
        with zopen(lower_filename, "w") as f:
            json.dump(structure, f, cls=PMGJSONEncoder)
            return
    else:
        raise ValueError("Unrecognized file extension!")

    writer.write_file(filename)


def read_mol(filename):
    """
    Reads a molecule based on file extension. For example, anything ending in
    a "xyz" is assumed to be a XYZ file. Supported formats include xyz,
    gaussian input (gjf|g03|g09|com|inp), Gaussian output (.out|and
    pymatgen's JSON serialized molecules. Using openbabel,
    many more extensions are supported but requires openbabel to be installed.

    Args:
        filename:
            A filename to read from.

    Returns:
        A Molecule object.
    """
    lower_filename = os.path.basename(filename).lower()
    if re.search("\.xyz", lower_filename):
        return XYZ.from_file(filename).molecule
    elif re.search("\.(gjf|g03|g09|com|inp)", lower_filename):
        return GaussianInput.from_file(filename).molecule
    elif re.search("\.(out|lis|log)", lower_filename):
        return GaussianOutput(filename).final_structure
    elif re.search("\.[mj]son", lower_filename):
        with zopen(lower_filename) as f:
            s = json.load(f, cls=PMGJSONDecoder)
            if type(s) != Molecule:
                raise IOError("File does not contain a valid serialized "
                              "molecule")
            return s
    else:
        m = re.search("\.(pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv)",
                      lower_filename)
        if m:
            return BabelMolAdaptor.from_file(filename,
                                             m.group(1)).pymatgen_mol

    raise ValueError("Unrecognized file extension!")


def write_mol(mol, filename):
    """
    Write a molecule to a file based on file extension. For example, anything
    ending in a "xyz" is assumed to be a XYZ file. Supported formats include
    xyz, Gaussian input (gjf|g03|g09|com|inp), and pymatgen's JSON serialized
    molecules.

    Args:
        mol:
            Molecule to write
        filename:
            A filename to write to.
    """
    lower_filename = os.path.basename(filename).lower()
    if re.search("\.xyz", lower_filename):
        return XYZ(mol).write_file(filename)
    elif re.search("\.(gjf|g03|g09|com|inp)", lower_filename):
        return GaussianInput(mol).write_file(filename)
    elif re.search("\.[mj]son", lower_filename):
        with zopen(lower_filename, "w") as f:
            return json.dump(mol, f, cls=PMGJSONEncoder)
    else:
        m = re.search("\.(pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv)",
                      lower_filename)
        if m:
            return BabelMolAdaptor(mol).write_file(filename, m.group(1))

    raise ValueError("Unrecognized file extension!")

