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

__author__ = "Alan Dozier, Kiran Mathew"
__credits__ = "Anubhav Jain, Shyue Ping Ong"
__version__ = "1.1"
__maintainer__ = "Kiran Mathew"
__email__ = "kmathew@lbl.gov"
__date__ = "Sept 10, 2016"

import os
import abc
import six
from copy import deepcopy

from monty.serialization import loadfn
from monty.json import MSONable

from pymatgen.io.feff.inputs import Atoms, Tags, Potential, Header


MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class AbstractFeffInputSet(six.with_metaclass(abc.ABCMeta, MSONable)):
    """
    Abstract base class representing a set of Feff input parameters.
    The idea is that using a FeffInputSet, a complete set of input files
    (feffPOT, feffXANES, feffEXAFS, ATOMS, feff.inp)set_
    can be generated in an automated fashion for any structure.
    """

    @abc.abstractmethod
    def header(self):
        """
        Returns header to be used in feff.inp file from a pymatgen structure
        """
        pass

    @abc.abstractproperty
    def atoms(self):
        """
        Returns Atoms string from a structure that goes in feff.inp file.

        Returns:
            Atoms object.
        """
        pass

    @abc.abstractproperty
    def tags(self):
        """
        Returns standard calculation parameters.
        """
        return

    @abc.abstractproperty
    def potential(self):
        """
        Returns POTENTIAL section used in feff.inp from a structure.
        """
        pass

    def all_input(self):
        """
        Returns all input files as a dict of {filename: feffio object}
        """
        return {"HEADER": self.header(),
                "PARAMETERS": self.tags,
                "POTENTIALS": self.potential,
                "ATOMS": self.atoms}

    def write_input(self, output_dir=".", make_dir_if_not_present=True):
        """
        Writes a set of FEFF input to a directory.

        Args:
            output_dir: Directory to output the FEFF input files
            make_dir_if_not_present: Set to True if you want the directory (
                and the whole path) to be created if it is not present.
        """
        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        feff = self.all_input()

        feff_input = "\n\n".join(str(feff[f]) for f in ["HEADER", "PARAMETERS",
                                 "POTENTIALS", "ATOMS"])

        for k, v in six.iteritems(feff):
            with open(os.path.join(output_dir, k), "w") as f:
                f.write(str(v))

        with open(os.path.join(output_dir, "feff.inp"), "w") as f:
            f.write(feff_input)
        f.close()


class FEFFDictSet(AbstractFeffInputSet):
    """
    Standard implementation of FeffInputSet, which can be extended by specific
    implementations.
    """

    def __init__(self, absorbing_atom, structure, radius, config_dict,
                  spectrum="EXAFS", user_tag_settings=None):
        """

        Args:
            absorbing_atom (str): absorbing atom symbol
            structure (Structure): input structure
            radius (float): cluster radius
            config_dict (dict): control tag settings dict
            spectrum (str): type of spectrum to calculate, available options :
                EXAFS, XANES, DANES, XMCD, ELNES, EXELFS, FPRIME, NRIXS, XES.
                The default is EXAFS.
            user_tag_settings (dict): override default tag settings. To delete
                tags, set the key '_del' in the user_tag_settings.
                eg: user_tag_settings={"_del": ["COREHOLE", "EXCHANGE"]}
        """
        self.absorbing_atom = absorbing_atom
        self.structure = structure
        self.radius = radius
        self.config_dict = deepcopy(config_dict)
        self.spectrum = spectrum
        self.user_tag_settings = user_tag_settings or {}
        self.config_dict.update(self.user_tag_settings)
        if "_del" in self.user_tag_settings:
            for tag in self.user_tag_settings["_del"]:
                if tag in self.config_dict:
                    del self.config_dict[tag]
            del self.config_dict["_del"]

    def header(self, source='', comment=''):
        """
        Creates header string from structure object

        Args:
            source: Source identifier used to create structure, can be defined
                however user wants to organize structures, calculations, etc.
                example would be Materials Project material ID number.
            comment: comment to include in header

        Returns:
            Header
        """
        return Header(self.structure, source, comment)

    @property
    def tags(self):
        """
        FEFF job parameters.

        Returns:
            Tags
        """
        return Tags(self.config_dict)

    @property
    def potential(self):
        """
        FEFF potential

        Returns:
            Potential
        """
        return Potential(self.structure, self.absorbing_atom)

    @property
    def atoms(self):
        """
        absorber + the rest

        Returns:
            Atoms
        """
        return Atoms(self.structure, self.absorbing_atom, self.radius)

    def __str__(self):
        output = [self.spectrum]
        output.extend(["%s = %s" % (k, str(v))
                       for k, v in six.iteritems(self.config_dict)])
        output.append("")
        return "\n".join(output)


class MPXANESSet(FEFFDictSet):
    """
    FeffDictSet for XANES spectroscopy.
    """

    CONFIG = loadfn(os.path.join(MODULE_DIR, "MPXANESSet.yaml"))

    def __init__(self, absorbing_atom, structure, radius=10., **kwargs):
        """
        Args:
            absorbing_atom (str): absorbing atom symbol
            structure (Structure): input structure
            radius (float): cluster radius in Angstroms.
            **kwargs
        """
        super(MPXANESSet, self).__init__(absorbing_atom, structure, radius,
                                         MPXANESSet.CONFIG, spectrum="XANES",
                                         **kwargs)
        self.kwargs = kwargs


class MPEXAFSSet(FEFFDictSet):
    """
    FeffDictSet for EXAFS spectroscopy.
    """

    CONFIG = loadfn(os.path.join(MODULE_DIR, "MPEXAFSSet.yaml"))

    def __init__(self, absorbing_atom, structure, radius=10., **kwargs):
        """
        Args:
            absorbing_atom (str): absorbing atom symbol
            structure (Structure): input structure
            radius (float): cluster radius in Angstroms.
            **kwargs
        """
        super(MPEXAFSSet, self).__init__(absorbing_atom, structure, radius,
                                         MPEXAFSSet.CONFIG, spectrum="EXAFS",
                                         **kwargs)
        self.kwargs = kwargs
