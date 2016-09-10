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
        Returns standard calculation paramters.
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
                "PARAMETERS": self.tags(),
                "POTENTIALS": self.potential(),
                "ATOMS": self.atoms()}

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

    def __init__(self, absorbing_atom, structure, config_dict, name="MPFEFF",
                 user_tag_settings=None):
        """

        Args:
            absorbing_atom (str): absorbing atom symbol
            structure (Structure): input structure
            config_dict (dict): control tag settings dict
            name (str)
            user_tag_settings (dict): override default tag settings
        """
        self.absorbing_atom = absorbing_atom
        self.structure = structure
        self.name = name
        self.config_dict = deepcopy(config_dict)
        self.user_tag_settings = user_tag_settings or {}

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
        settings = dict(self.config_dict)
        settings.update(self.user_tag_settings)

        return Tags(settings)

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
        return Atoms(self.structure, self.absorbing_atom)

    def __str__(self):
        d = self.config_dict
        output = [self.name]
        output.extend(["%s = %s" % (k, str(v)) for k, v in six.iteritems(d)])
        output.append("")
        return "\n".join(output)


class MPXANESSet(FEFFDictSet):
    """
    FeffDictSet for XANES
    """

    CONFIG = loadfn(os.path.join(MODULE_DIR, "MPXANESSet.yaml"))

    def __init__(self, absorbing_atom, structure, name="MPXANES", **kwargs):
        super(MPXANESSet, self).__init__(absorbing_atom, structure,
                                         MPXANESSet.CONFIG, name, **kwargs)
        self.kwargs = kwargs


class MPEXAFSSet(FEFFDictSet):
    """
    FeffDictSet for EXAFS.
    """

    CONFIG = loadfn(os.path.join(MODULE_DIR, "MPEXAFSSet.yaml"))

    def __init__(self, absorbing_atom, structure, name="MPEXAFS", **kwargs):
        super(MPEXAFSSet, self).__init__(absorbing_atom, structure,
                                         MPEXAFSSet.CONFIG, name, **kwargs)
        self.kwargs = kwargs
