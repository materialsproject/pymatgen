# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module defines the FeffInputSet abstract base class and a concrete
implementation for the Materials Project.  The basic concept behind an input
set is to specify a scheme to generate a consistent set of Feff inputs from a
structure without further user intervention. This ensures comparability across
runs.
"""

import six

__author__ = "Alan Dozier"
__credits__ = "Anubhav Jain, Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0.3"
__maintainer__ = "Alan Dozier"
__email__ = "adozier@uky.edu"
__date__ = "April 7, 2013"

import os
import abc

from monty.serialization import loadfn

from pymatgen.io.feff import FeffAtoms, FeffTags, FeffPot, Header


class AbstractFeffInputSet(six.with_metaclass(abc.ABCMeta, object)):
    """
    Abstract base class representing a set of Feff input parameters.
    The idea is that using a FeffInputSet, a complete set of input files
    (feffPOT,feffXANES, feffEXAFS, ATOMS, feff.inp)set_
    can be generated in an automated fashion for any structure.
    """

    @abc.abstractmethod
    def get_feff_atoms(self, structure, central_atom):
        """
        Returns Atoms string from a structure that goes in feff.inp file.

        Args:
            structure: pymatgen structure object
            central_atom: atom symbol string for absorbing atom

        Returns:
            FeffAtoms object.
        """
        return

    @abc.abstractmethod
    def get_feff_tags(self, calc_type):
        """
        Returns standard calculation paramters for either an FEFF XANES or
        EXAFS input.

        Args:
            calc_type: At this time either 'XANES' or 'EXAFS' string is
                supported for K shell excitation. In the future this will be
                expanded to include other shells and material class
                differentiation.
        """
        return

    @abc.abstractmethod
    def get_feff_pot(self, structure, central_atom):
        """
        Returns POTENTIAL section used in feff.inp from a structure.

        Args:
            structure: pymatgen structure object
            central_atom: atom symbol string for absorbing atom
        """
        return

    @abc.abstractmethod
    def get_header(self, structure, source, comment):
        """
        Returns header to be used in feff.inp file from a pymatgen structure

        Args:
            structure: A pymatgen structure object
            source: Source identifier used to create structure, can be defined
                however user wants to organize structures, calculations, etc.
                example would be Materials Project material ID number.
        """
        return

    def get_all_feff_input(self, structure, calc_type, source, central_atom,
                           comment=''):
        """
        Returns all input files as a dict of {filename: feffio object}

        Args:
            structure: Structure object
            calc_type: At this time either 'XANES' or 'EXAFS' string is
                supported for K shell excitation. In the future this will be
                expanded to inlude other shells and material class
                differentiation.
            source: Source identifier used to create structure, can be defined
                however user wants to organize structures, calculations, etc.
                example would be Materials Project material ID number.
            central_atom: Atom symbol string for absorbing atom
            comment: Comment to appear in Header.

        Returns:
            dict of objects used to create feff.inp file i.e. Header, FeffTags,
            FeffPot, FeffAtoms
        """

        feff = {"HEADER": self.get_header(structure, source, comment),
                "PARAMETERS": self.get_feff_tags(calc_type),
                "POTENTIALS": self.get_feff_pot(structure, central_atom),
                "ATOMS": self.get_feff_atoms(structure, central_atom)}

        return feff

    def write_input(self, structure, calc_type, source, central_atom,
                    comment='', output_dir=".", make_dir_if_not_present=True):
        """
        Writes a set of FEFF input to a directory.

        Args:
            structure: Structure object
            calc_type: At this time either 'XANES' or 'EXAFS' string is
                supported for K shell excitation. In the future this will be
                expanded to include other shells and material class
                differentiation.
            source: Source identifier used to create structure, can be defined
                however user wants to organize structures, calculations, etc.
                example would be Materials Project material ID number.
            central_atom: Atom symbol string for absorbing atom
            output_dir: Directory to output the FEFF input files
            comment: comment for Header
            make_dir_if_not_present: Set to True if you want the directory (
                and the whole path) to be created if it is not present.
        """

        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        feff = self.get_all_feff_input(structure, calc_type, source,
                                       central_atom, comment)

        feff_input = "\n\n".join(str(feff[f]) for f in ["HEADER", "PARAMETERS",
                                 "POTENTIALS", "ATOMS"])

        for k, v in six.iteritems(feff):
            with open(os.path.join(output_dir, k), "w") as f:
                f.write(str(v))

        with open(os.path.join(output_dir, "feff.inp"), "w") as f:
            f.write(feff_input)
        f.close()

    def as_dict(self, structure, calc_type, source, central_atom,
                comment=''):
        """Creates a feff.inp dictionary as a string"""

        feff = self.get_all_feff_input(structure, calc_type, source,
                                       central_atom, comment)
        feff_input = "\n\n".join(str(feff[f]) for f in ["HEADER", "PARAMETERS",
                                 "POTENTIALS", "ATOMS"])
        return {'@module': self.__class__.__module__,
                '@class': self.__class__.__name__,
                'feff.inp': feff_input}

    @staticmethod
    def from_dict(d):
        """Return feff.inp from a dictionary string representation"""
        return d['feff.inp']


class FeffInputSet(AbstractFeffInputSet):
    """
    Standard implementation of FeffInputSet, which can be extended by specific
    implementations.

    Args:
        name: The name of a grouping of input parameter sets such as
            "MaterialsProject".
    """

    def __init__(self, name):
        self.name = name
        module_dir = os.path.dirname(os.path.abspath(__file__))
        config = loadfn(os.path.join(module_dir, "FeffInputSets.yaml"))
        self.xanes_settings = config[self.name + "feffXANES"]
        self.exafs_settings = config[self.name + "feffEXAFS"]

    def get_header(self, structure, source='', comment=''):
        """
        Creates header string from structure object

        Args:
            structure: A pymatgen structure object
            source: Source identifier used to create structure, can be defined
                however user wants to organize structures, calculations, etc.
                example would be Materials Project material ID number.
            comment: comment to include in header

        Returns:
            Header object to be used in feff.inp file from a pymatgen structure
        """
        return Header(structure, source, comment)

    def get_feff_tags(self, calc_type):
        """
        Reads standard parameters for XANES or EXAFS calculation
        from FeffInputSets.yaml file.

        Args:
            calc_type: At this time either 'XANES' or 'EXAFS' string is
                supported for K shell excitation. In the future this will be
                expanded to include other shells and material class
                differentiation.

        Returns:
            FeffTags object
        """

        if calc_type.upper() == "XANES":
            fefftags = FeffTags(self.xanes_settings)
        elif calc_type.upper() == "EXAFS":
            fefftags = FeffTags(self.exafs_settings)
        else:
            raise ValueError("{} is not a valid calculation type"
                             .format(calc_type))

        return fefftags

    def get_feff_pot(self, structure, central_atom):
        """
        Creates string representation of potentials used in POTENTIAL file and
        feff.inp.

        Args:
            structure: pymatgen structure object
            central_atom: atom symbol string for absorbing atom

        Returns:
            FeffPot object
        """
        return FeffPot(structure, central_atom)

    def get_feff_atoms(self, structure, central_atom):
        """
        Creates string representation of atomic shell coordinates using in
        ATOMS file and feff.inp.

        Args:
            structure: pymatgen structure object
            central_atom: atom symbol string for absorbing atom

        Returns:
            FeffAtoms object
        """
        return FeffAtoms(structure, central_atom)

    def __str__(self):
        output = [self.name]
        section_names = ["XANES", "EXAFS"]
        for ns in section_names:
            for d in [self.xanes_settings, self.exafs_settings]:
                output.append(ns)
                for k, v in six.iteritems(d):
                    output.append("%s = %s" % (k, str(v)))
                output.append("")

        return "\n".join(output)


class MaterialsProjectFeffInputSet(FeffInputSet):
    """
    Implementation of FeffInputSet utilizing parameters in the public
    Materials Project.
    """
    def __init__(self):
        super(MaterialsProjectFeffInputSet, self).__init__("MaterialsProject")
