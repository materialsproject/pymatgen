#!/usr/bin/env python

"""
This module defines the FeffInputSet abstract base class and a concrete
implementation for the Materials Project.  The basic concept behind an input
set is to specify a scheme to generate a consistent set of Feff inputs from a
structure without further user intervention. This ensures comparability across
runs.
"""

from __future__ import division

__author__ = "Alan Dozier"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0.1"
__maintainer__ = "Alan Dozier"
__email__ = "adozier@uky.edu"
__date__ = "September 17, 2012"

import os
import abc
import ConfigParser

from pymatgen.io.feffio import FeffAtoms, FeffTags, FeffPot, Header


class AbstractFeffInputSet(object):
    """
    Abstract base class representing a set of Feff input parameters.
    The idea is that using a FeffInputSet, a complete set of input files
    (feffPOT,feffXANES, feffEXAFS, ATOMS, feff.inp)
    can be generated in an automated fashion for any structure.
    """
    __metaclass__ = abc.ABCMeta

    def get_feff_atoms(self, structure, central_atom):
        """
        Returns Atoms string from a structure that goes in feff.inp file.
        """
        return FeffAtoms(structure)

    @abc.abstractmethod
    def get_feff_tags(self, calctype):
        """
        Returns standard calculation paramters for either an FEFF XANES or
        EXAFS input.
        """
        return

    @abc.abstractmethod
    def get_feff_pot(self, structure, central_atom):
        """
        Returns POTENTIAL section used in feff.inp from a structure.

        Args:
            structure:
                Structure object
        """
        return

    @abc.abstractmethod
    def get_header(self, structure, cif_file):
        """
        Returns header to be used in feff.inp file from a structure and
        cif_file

        Args:
            structure:
                A structure
            cif_file:
                Name of cif_file used to creat structure
        """
        return

    def get_all_feff_input(self, structure, calc_type, cif_file, central_atom):
        """
        Returns all input files as a dict of {filename: feffio object}

        Args:
            structure:
                Structure object
            calc_type:
                XANES or EXAFS
            cif_file:
                path and name of cif file of material
            central:
                symbol of absorbing atom

        Returns:
            dict of {filename: file_as_string}, e.g., {"INCAR":"EDIFF=1e-4..."}
        """
        feff = {"HEADER": self.get_header(structure, cif_file),
                "PARAMETERS": self.get_feff_tags(calc_type),
                "POTENTIALS": self.get_feff_pot(structure, central_atom),
                "ATOMS": self.get_feff_atoms(structure, central_atom)}

        return feff

    def write_input(self, structure, calc_type, cif_file, output_dir,
                    central_atom, make_dir_if_not_present=True):
        """
        Writes a set of FEFF input to a directory.

        Args:
            structure:
                Structure object
            calc_type:
                XANES or EXAFS
            cif_file:
                path and name of cif file of material
            central:
                symbol of absorbing atom
            output_dir:
                Directory to output the FEFF input files
            make_dir_if_not_present:
                Set to True if you want the directory (and the whole path) to
                be created if it is not present.
        """
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        feff = self.get_all_feff_input(structure, calc_type, cif_file,
                                       central_atom)
        feff_input = str(feff["HEADER"]) + "\n\n" + str(feff["PARAMETERS"]) + \
            "\n\n" + str(feff["POTENTIALS"]) + "\n\n" + str(feff["ATOMS"])
        for k, v in self.get_all_feff_input(structure, calc_type, cif_file,
                                            central_atom).items():
            with open(os.path.join(output_dir, k), "w") as f:
                f.write(str(v))
        with open(os.path.join(output_dir, "feff.inp"), "w") as f:
            f.write(feff_input)


class FeffInputSet(AbstractFeffInputSet):
    """
    Standard implementation of FeffInputSet, which can be extended by specific
    implementations.
    """

    def __init__(self, name):
        """
        Args:
            name:
                The name of a grouping of input parameter sets such as
                "MaterialsProject".
        """
        self.name = name
        module_dir = os.path.dirname(os.path.abspath(__file__))
        self._config = ConfigParser.SafeConfigParser()
        self._config.optionxform = str
        self._config.readfp(open(os.path.join(module_dir,
                                              "FeffInputSets.cfg")))
        self.xanes_settings = dict(self._config.items(self.name + "feffXANES"))
        self.exafs_settings = dict(self._config.items(self.name + "feffEXAFS"))

    def get_header(self, structure, cif_file):
        """
        Creates header string from cif file and structure object

        Args:
            structure:
                pymatgen structure object
            cif_file:
                path and cif_file name
        Returns:
            HEADER string
        """
        space_number, space_group = Header.structure_symmetry(structure)
        header = Header(structure, cif_file, space_number, space_group)
        return header.get_string()

    def get_feff_tags(self, calc_type):
        """
        Reads standard parameters for XANES or EXAFS calculation
        from FeffInputSets.cfg file.

        Args:
            calc_type:
                XANES or EXAFS

        Returns:
            FeffTags object
        """
        if calc_type == "XANES":
            fefftags = FeffTags(self.xanes_settings)
        else:
            fefftags = FeffTags(self.exafs_settings)
        return fefftags

    def get_feff_pot(self, structure, central_atom):
        """
        Creates string representation of potentials used in POTENTIAL file and
        feff.inp.

        Args:
            structure:
                structure object
            central_atom:
                symbol for absorbing atom

        Returns:
            string representation of potential indicies, etc. used in POTENTIAL
            file.
        """
        pot = FeffPot(structure, central_atom)
        return pot.get_string()

    def get_feff_atoms(self, structure, central_atom):
        """
        Creates string representation of atomic shell coordinates using in
        ATOMS file and feff.inp.

        Args:
            structure:
                structure object
            central_atom:
                symbol for absorbing atom

        Returns:
            String representation of atoms file.
        """
        atoms = FeffAtoms(structure, central_atom)
        return atoms.get_string()

    def __str__(self):
        output = [self.name]
        output.append("")
        section_names = ["XANES", "EXAFS"]
        count = 0
        for d in [self.xanes_settings, self.exafs_settings]:
            output.append(section_names[count])
            for k, v in d.items():
                output.append("%s = %s" % (k, str(v)))
            output.append("")
            count += 1

        return "\n".join(output)


class MaterialsProjectFeffInputSet(FeffInputSet):
    """
    Implementation of FeffInputSet utilizing parameters in the public
    Materials Project.
    """
    def __init__(self):
        super(MaterialsProjectFeffInputSet, self).__init__("MaterialsProject")
