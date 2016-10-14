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

import sys
import os
import abc
import six
from copy import deepcopy
import logging

from monty.serialization import loadfn
from monty.json import MSONable

from pymatgen.io.feff.inputs import Atoms, Tags, Potential, Header


__author__ = "Kiran Mathew"
__credits__ = "Alan Dozier, Anubhav Jain, Shyue Ping Ong"
__version__ = "1.1"
__maintainer__ = "Kiran Mathew"
__email__ = "kmathew@lbl.gov"
__date__ = "Sept 10, 2016"


MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(name)s: %(message)s')
sh = logging.StreamHandler(stream=sys.stdout)
sh.setFormatter(formatter)
logger.addHandler(sh)


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
        d = {"HEADER": self.header(), "PARAMETERS": self.tags}

        if "RECIPROCAL" not in self.tags:
            d.update({"POTENTIALS": self.potential, "ATOMS": self.atoms})

        return d

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

        feff_input = "\n\n".join(str(feff[k]) for k in
                                 ["HEADER", "PARAMETERS", "POTENTIALS", "ATOMS"]
                                 if k in feff)

        for k, v in six.iteritems(feff):
            with open(os.path.join(output_dir, k), "w") as f:
                f.write(str(v))

        with open(os.path.join(output_dir, "feff.inp"), "w") as f:
            f.write(feff_input)

        # write the structure to cif file
        if "ATOMS" not in feff:
            self.atoms.struct.to(fmt="cif",
                                 filename=os.path.join(
                                     output_dir, feff["PARAMETERS"]["CIF"]))


class FEFFDictSet(AbstractFeffInputSet):
    """
    Standard implementation of FeffInputSet, which can be extended by specific
    implementations.
    """

    def __init__(self, absorbing_atom, structure, radius, config_dict,
                 edge="K", spectrum="EXAFS", nkpts=1000, user_tag_settings=None):
        """

        Args:
            absorbing_atom (str/int): absorbing atom symbol or site index
            structure (Structure): input structure
            radius (float): cluster radius
            config_dict (dict): control tag settings dict
            edge (str): absorption edge
            spectrum (str): type of spectrum to calculate, available options :
                EXAFS, XANES, DANES, XMCD, ELNES, EXELFS, FPRIME, NRIXS, XES.
                The default is EXAFS.
            nkpts (int): Total number of kpoints in the brillouin zone. Used
                only when feff is run in the reciprocal space mode.
            user_tag_settings (dict): override default tag settings. To delete
                tags, set the key '_del' in the user_tag_settings.
                eg: user_tag_settings={"_del": ["COREHOLE", "EXCHANGE"]}
        """
        self.absorbing_atom = absorbing_atom
        self.structure = structure
        self.radius = radius
        self.config_dict = deepcopy(config_dict)
        self.edge = edge
        self.spectrum = spectrum
        self.nkpts = nkpts
        self.user_tag_settings = user_tag_settings or {}
        self.config_dict["EDGE"] = self.edge
        self.config_dict.update(self.user_tag_settings)
        if "_del" in self.user_tag_settings:
            for tag in self.user_tag_settings["_del"]:
                if tag in self.config_dict:
                    del self.config_dict[tag]
            del self.config_dict["_del"]
        # k-space feff only for small systems. The hardcoded system size in
        # feff is around 14 atoms.
        self.small_system = True if len(self.structure) < 14 else False

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
        if "RECIPROCAL" in self.config_dict:
            if self.small_system:
                self.config_dict["CIF"] = "{}.cif".format(
                    self.structure.formula.replace(" ", ""))
                self.config_dict["TARGET"] = self.atoms.center_index + 1
                self.config_dict["COREHOLE"] = "RPA"
                logger.warn("Setting COREHOLE = RPA for K-space calculation")
                if not self.config_dict.get("KMESH", None):
                    abc = self.structure.lattice.abc
                    mult = (self.nkpts * abc[0] * abc[1] * abc[2]) ** (1 / 3)
                    self.config_dict["KMESH"] = [int(round(mult / l)) for l in abc]
            else:
                logger.warn("Large system(>=14 atoms), removing K-space settings")
                del self.config_dict["RECIPROCAL"]
                self.config_dict.pop("CIF", None)
                self.config_dict.pop("TARGET", None)
                self.config_dict.pop("KMESH", None)
                self.config_dict.pop("STRFAC", None)

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

    def __init__(self, absorbing_atom, structure, edge="K", radius=10.,
                 nkpts=1000, **kwargs):
        """
        Args:
            absorbing_atom (str/int): absorbing atom symbol or site index
            structure (Structure): input
            edge (str): absorption edge
            radius (float): cluster radius in Angstroms.
            nkpts (int): Total number of kpoints in the brillouin zone. Used
                only when feff is run in the reciprocal space mode.
            **kwargs
        """
        super(MPXANESSet, self).__init__(absorbing_atom, structure, radius,
                                         MPXANESSet.CONFIG, edge=edge,
                                         spectrum="XANES", nkpts=nkpts, **kwargs)
        self.kwargs = kwargs


class MPEXAFSSet(FEFFDictSet):
    """
    FeffDictSet for EXAFS spectroscopy.
    """

    CONFIG = loadfn(os.path.join(MODULE_DIR, "MPEXAFSSet.yaml"))

    def __init__(self, absorbing_atom, structure, edge="K", radius=10.,
                 nkpts=1000, **kwargs):
        """
        Args:
            absorbing_atom (str/int): absorbing atom symbol or site index
            structure (Structure): input structure
            edge (str): absorption edge
            radius (float): cluster radius in Angstroms.
            nkpts (int): Total number of kpoints in the brillouin zone. Used
                only when feff is run in the reciprocal space mode.
            **kwargs
        """
        super(MPEXAFSSet, self).__init__(absorbing_atom, structure, radius,
                                         MPEXAFSSet.CONFIG, edge=edge,
                                         spectrum="EXAFS", nkpts=nkpts, **kwargs)
        self.kwargs = kwargs


class MPEELSDictSet(FEFFDictSet):
    """
    FeffDictSet for ELNES spectroscopy.
    """

    def __init__(self, absorbing_atom, structure, edge, spectrum, radius,
                 beam_energy, beam_direction, collection_angle,
                 convergence_angle, config_dict, user_eels_settings=None,
                 nkpts=1000, **kwargs):
        """
        Args:
            absorbing_atom (str/int): absorbing atom symbol or site index
            structure (Structure): input structure
            edge (str): absorption edge
            spectrum (str): ELNES or EXELFS
            radius (float): cluster radius in Angstroms.
            beam_energy (float): Incident beam energy in keV
            beam_direction (list): Incident beam direction. If None, the
                cross section will be averaged.
            collection_angle (float): Detector collection angle in mrad.
            convergence_angle (float): Beam convergence angle in mrad.
            user_eels_settings (dict): override default EELS config.
                See MPELNESSet.yaml for supported keys.
            nkpts (int): Total number of kpoints in the brillouin zone. Used
                only when feff is run in the reciprocal space mode.
            **kwargs
        """
        self.beam_energy = beam_energy
        self.beam_direction = beam_direction
        self.collection_angle = collection_angle
        self.convergence_angle = convergence_angle
        self.user_eels_settings = user_eels_settings
        eels_config_dict = deepcopy(config_dict)

        if beam_direction:
            beam_energy_list = [beam_energy, 0, 1, 1]
            eels_config_dict[spectrum]["BEAM_DIRECTION"] = beam_direction
        else:
            beam_energy_list = [beam_energy, 1, 0, 1]
            del eels_config_dict[spectrum]["BEAM_DIRECTION"]
        eels_config_dict[spectrum]["BEAM_ENERGY"] = beam_energy_list
        eels_config_dict[spectrum]["ANGLES"] = [collection_angle,
                                                convergence_angle]

        if user_eels_settings:
            eels_config_dict[spectrum].update(user_eels_settings)

        super(MPEELSDictSet, self).__init__(absorbing_atom, structure, radius,
                                            eels_config_dict, edge=edge,
                                            spectrum=spectrum, nkpts=nkpts,
                                            **kwargs)
        self.kwargs = kwargs


class MPELNESSet(MPEELSDictSet):
    """
    FeffDictSet for ELNES spectroscopy.
    """

    CONFIG = loadfn(os.path.join(MODULE_DIR, "MPELNESSet.yaml"))

    def __init__(self, absorbing_atom, structure, edge="K", radius=10.,
                 beam_energy=100, beam_direction=None, collection_angle=1,
                 convergence_angle=1, user_eels_settings=None, nkpts=1000,
                 **kwargs):
        """
        Args:
            absorbing_atom (str/int): absorbing atom symbol or site index
            structure (Structure): input structure
            edge (str): absorption edge
            radius (float): cluster radius in Angstroms.
            beam_energy (float): Incident beam energy in keV
            beam_direction (list): Incident beam direction. If None, the
                cross section will be averaged.
            collection_angle (float): Detector collection angle in mrad.
            convergence_angle (float): Beam convergence angle in mrad.
            user_eels_settings (dict): override default EELS config.
                See MPELNESSet.yaml for supported keys.
            nkpts (int): Total number of kpoints in the brillouin zone. Used
                only when feff is run in the reciprocal space mode.
            **kwargs
        """

        super(MPELNESSet, self).__init__(absorbing_atom, structure, edge,
                                         "ELNES", radius, beam_energy,
                                         beam_direction, collection_angle,
                                         convergence_angle, MPELNESSet.CONFIG,
                                         user_eels_settings=user_eels_settings,
                                         nkpts=nkpts, **kwargs)
        self.kwargs = kwargs


class MPEXELFSSet(MPEELSDictSet):
    """
    FeffDictSet for EXELFS spectroscopy.
    """

    CONFIG = loadfn(os.path.join(MODULE_DIR, "MPEXELFSSet.yaml"))

    def __init__(self, absorbing_atom, structure, edge="K", radius=10.,
                 beam_energy=100, beam_direction=None, collection_angle=1,
                 convergence_angle=1, user_eels_settings=None, nkpts=1000,
                 **kwargs):
        """
        Args:
            absorbing_atom (str/int): absorbing atom symbol or site index
            structure (Structure): input structure
            edge (str): absorption edge
            radius (float): cluster radius in Angstroms.
            beam_energy (float): Incident beam energy in keV
            beam_direction (list): Incident beam direction. If None, the
                cross section will be averaged.
            collection_angle (float): Detector collection angle in mrad.
            convergence_angle (float): Beam convergence angle in mrad.
            user_eels_settings (dict): override default EELS config.
                See MPEXELFSSet.yaml for supported keys.
            nkpts (int): Total number of kpoints in the brillouin zone. Used
                only when feff is run in the reciprocal space mode.
            **kwargs
        """

        super(MPEXELFSSet, self).__init__(absorbing_atom, structure, edge,
                                          "EXELFS", radius, beam_energy,
                                          beam_direction, collection_angle,
                                          convergence_angle, MPEXELFSSet.CONFIG,
                                          user_eels_settings=user_eels_settings,
                                          nkpts=nkpts, **kwargs)
        self.kwargs = kwargs
