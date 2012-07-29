#!/usr/bin/env python

"""
OpenBabel interface module, which opens up access to the hundreds of file
formats supported by OpenBabel. Requires openbabel with python bindings to be
installed. Please consult the
`openbabel documentation <http://openbabel.org/wiki/Main_Page>`_.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 28, 2012"

from pymatgen.core.structure import Molecule

try:
    import openbabel as ob
    import pybel as pb
    babel_loaded = True
except ImportError:
    babel_loaded = False


class BabelMolAdaptor(object):
    """
    Adaptor serves as a bridge between OpenBabel"s Molecule and pymatgen"s
    Molecule.
    """

    def __init__(self, mol):
        """
        Initializes with pymatgen Molecule or OpenBabel"s OBMol.

        Args:
            mol:
                pymatgen"s Molecule or OpenBabel Molecule
        """
        if isinstance(mol, Molecule):
            if not mol.is_ordered:
                raise ValueError("OpenBabel Molecule only supports ordered "
                                 "molecules.")
            obmol = ob.OBMol()
            for site in mol:
                a = obmol.NewAtom()
                a.SetAtomicNum(site.specie.Z)
                a.SetVector(site.x, site.y, site.z)
            self._mol = mol
            self._obmol = obmol
        elif isinstance(mol, ob.OBMol):
            sp = []
            coords = []
            for atom in ob.OBMolAtomIter(mol):
                sp.append(atom.GetAtomicNum())
                coords.append([atom.GetX(), atom.GetY(), atom.GetZ()])
            self._mol = Molecule(sp, coords)
            self._obmol = mol

    @property
    def pymatgen_mol(self):
        """
        Returns pymatgen"s Molecule representation.
        """
        return self._mol

    @property
    def openbabel_mol(self):
        """
        Returns OpenBabel's OBMol representation.
        """
        return self._obmol

    def write_file(self, filename, file_format="xyz"):
        """
        Uses OpenBabel to output all supported formats.

        Args:
            filename:
                Filename of file to output
            file_format:
                String specifying any OpenBabel supported formats.
        """
        mol = pb.Molecule(self._obmol)
        mol.write(file_format, filename, overwrite=True)

    @staticmethod
    def from_file(filename, file_format="xyz"):
        """
        Uses OpenBabel to read a molecule from all supported formats.

        Args:
            filename:
                Filename of file to output
            file_format:
                String specifying any OpenBabel supported formats.

        Returns:
            BabelMolAdaptor object
        """
        mols = list(pb.readfile(file_format, filename))
        return BabelMolAdaptor(mols[0].OBMol)
