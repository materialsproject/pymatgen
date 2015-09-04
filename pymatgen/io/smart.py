# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This class implements smart io classes that performs intelligent io based on
file extensions.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 29, 2012"

import os
import json
import re
from fnmatch import fnmatch

from monty.dev import deprecated

from pymatgen.core.structure import Structure, Molecule
from pymatgen.io.vasp import Vasprun, Poscar, Chgcar
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.io.cssr import Cssr
from pymatgen.io.xyz import XYZ
from pymatgen.io.gaussian import GaussianInput, GaussianOutput
from monty.io import zopen
from monty.json import MontyDecoder, MontyEncoder
from monty.string import str2unicode
from pymatgen.io.babel import BabelMolAdaptor


@deprecated(Structure.from_file, message="Will be removed in pymatgen 4.0.")
def read_structure(filename, primitive=True, sort=False):
    """
    Reads a structure based on file extension. For example, anything ending in
    a "cif" is assumed to be a Crystallographic Information Format file.
    Supported formats include CIF, POSCAR/CONTCAR, CHGCAR, LOCPOT,
    vasprun.xml, CSSR and pymatgen's JSON serialized structures.

    Args:
        filename (str): A filename to read from.
        primitive (bool): Whether to convert to a primitive cell for cifs.
            Defaults to True.
        sort (bool): Whether to sort sites. Default to False.

    Returns:
        A Structure object.
    """
    fname = os.path.basename(filename)
    if fnmatch(fname.lower(), "*.cif*"):
        parser = CifParser(filename)
        s = parser.get_structures(primitive=primitive)[0]
    elif fnmatch(fname, "POSCAR*") or fnmatch(fname, "CONTCAR*"):
        s = Poscar.from_file(filename, False).structure
    elif fnmatch(fname, "CHGCAR*") or fnmatch(fname, "LOCPOT*"):
        s = Chgcar.from_file(filename).structure
    elif fnmatch(fname, "vasprun*.xml*"):
        s = Vasprun(filename).final_structure
    elif fnmatch(fname.lower(), "*.cssr*"):
        cssr = Cssr.from_file(filename)
        s = cssr.structure
    elif fnmatch(fname, "*.json*") or fnmatch(fname, "*.mson*"):
        with zopen(filename) as f:
            s = json.load(f, cls=MontyDecoder)
            if type(s) != Structure:
                raise IOError("File does not contain a valid serialized "
                              "structure")
    else:
        raise ValueError("Unrecognized file extension!")
    if sort:
        s = s.get_sorted_structure()
    return s


@deprecated(replacement=Structure.to, message="Will be removed in pymatgen 4.0")
def write_structure(structure, filename):
    """
    Write a structure to a file based on file extension. For example, anything
    ending in a "cif" is assumed to be a Crystallographic Information Format
    file. Supported formats include CIF, POSCAR, CSSR and pymatgen's JSON
    serialized structures.

    Args:
        structure (Structure/IStructure): Structure to write
        filename (str): A filename to write to.
    """
    fname = os.path.basename(filename)
    if fnmatch(fname, "*.cif*"):
        writer = CifWriter(structure)
    elif fnmatch(fname, "POSCAR*") or fnmatch(fname, "CONTCAR*"):
        writer = Poscar(structure)
    elif fnmatch(fname.lower(), "*.cssr*"):
        writer = Cssr(structure)
    elif fnmatch(fname, "*.json*") or fnmatch(fname, "*.mson*"):
        with zopen(filename, "wt") as f:
            f.write(str2unicode(json.dumps(structure, cls=MontyEncoder)))
            return
    else:
        raise ValueError("Unrecognized file extension!")

    writer.write_file(filename)


@deprecated(Molecule.from_file, message="Will be removed in pymatgen 4.0.")
def read_mol(filename):
    """
    Reads a molecule based on file extension. For example, anything ending in
    a "xyz" is assumed to be a XYZ file. Supported formats include xyz,
    gaussian input (gjf|g03|g09|com|inp), Gaussian output (.out|and
    pymatgen's JSON serialized molecules. Using openbabel,
    many more extensions are supported but requires openbabel to be installed.

    Args:
        filename (str): A filename to read from.

    Returns:
        A Molecule object.
    """
    fname = os.path.basename(filename)
    if fnmatch(fname.lower(), "*.xyz*"):
        return XYZ.from_file(filename).molecule
    elif any([fnmatch(fname.lower(), "*.{}*".format(r))
              for r in ["gjf", "g03", "g09", "com", "inp"]]):
        return GaussianInput.from_file(filename).molecule
    elif any([fnmatch(fname.lower(), "*.{}*".format(r))
              for r in ["out", "lis", "log"]]):
        return GaussianOutput(filename).final_structure
    elif fnmatch(fname, "*.json*") or fnmatch(fname, "*.mson*"):
        with zopen(filename) as f:
            s = json.load(f, cls=MontyDecoder)
            if type(s) != Molecule:
                raise IOError("File does not contain a valid serialized "
                              "molecule")
            return s
    else:
        m = re.search("\.(pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv)",
                      filename.lower())
        if m:
            return BabelMolAdaptor.from_file(filename,
                                             m.group(1)).pymatgen_mol

    raise ValueError("Unrecognized file extension!")


@deprecated(replacement=Molecule.to, message="Will be removed in pymatgen 4.0.")
def write_mol(mol, filename):
    """
    Write a molecule to a file based on file extension. For example, anything
    ending in a "xyz" is assumed to be a XYZ file. Supported formats include
    xyz, Gaussian input (gjf|g03|g09|com|inp), and pymatgen's JSON serialized
    molecules.

    Args:
        mol (Molecule/IMolecule): Molecule to write
        filename (str): A filename to write to.
    """
    fname = os.path.basename(filename)
    if fnmatch(fname.lower(), "*.xyz*"):
        return XYZ(mol).write_file(filename)
    elif any([fnmatch(fname.lower(), "*.{}*".format(r))
              for r in ["gjf", "g03", "g09", "com", "inp"]]):
        return GaussianInput(mol).write_file(filename)
    elif fnmatch(fname, "*.json*") or fnmatch(fname, "*.mson*"):
        with zopen(filename, "wt") as f:
            return f.write(str2unicode(json.dumps(mol, cls=MontyEncoder)))
    else:
        m = re.search("\.(pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv)",
                      filename.lower())
        if m:
            return BabelMolAdaptor(mol).write_file(filename, m.group(1))

    raise ValueError("Unrecognized file extension!")
