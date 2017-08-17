# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from monty.json import MSONable
import abc
from scipy.interpolate import interp1d
import json

"""
This module defines classes to represent all spectra
"""

__author__ = "Chen Zheng"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Chen Zheng"
__email__ = "chz022@ucsd.edu"
__date__ = "Aug 9, 2017"


class Spectrum(abc.ABC):
    """
    Base class for Xanespectrum and Exafspectrum.
    Not meant to be instantiated directly
    """

    @property
    def structure(self):
        """
        Returns the structure associated with the spectrum
        """
        return self.structure

    @property
    def spectrum_type(self):
        """
        Returns the spectrum type
        """
        return self.type

    @property
    def spectrum(self):
        """
        Returns the spectrum
        """
        return self.spectrum

    @classmethod
    @abc.abstractmethod
    def from_file(cls, filename):
        """
        Initiate Spectrum from a filename
        """
        pass


class Xanespectrum(Spectrum, MSONable):
    def __init__(self, formula, absorption_specie, edge, structure, spectrum):
        """
        Create a XANES spectrum object
        :param formula: The reduced formula of the structure associated with the spectrum
        :param absorption_specie: absorption specie of the spectrum
        :param edge: absorption spectrum edge
        :param structure: materials structure associated with the spectrum, need to be a pymatgen structure object
        :param spectrum: Spectrum data
        """

        self.formula = formula
        self.absorption_specie = absorption_specie
        self.edge = edge
        self.structure = structure
        self.spectrum = spectrum


class XanesFEFF(Xanespectrum):
    """
    This is a standard class for FEFF calculated XANES spectrum
    """

    def __init__(self, dir_name, formula, absorption_specie, edge, structure,
                 spectrum, input_parameters, mp_id):
        """
        Create an XanesFEFF spectrum object
        :param dir_name: XANES calculation directory name
        :param formula: The formula of the structure
        :param absorption_specie: absorption species
        :param edge: spectrum edge
        :param corehole: corehole setting in FEFF
        :param structure: pymatgen structure object used for FEFF XANES calculation
        :param spectrum: spectrum data
        :param scf: scf value of FEFF calculation
        :param fms: fms value of FEFF calculation
        :param exchange: exchange value of FEFF calculation
        :param s02: s02 value of FEFF calculation
        :param xanes: xanes value of FEFF calculation
        :param mp_id: mp_id of the structure
        """
        self.dir_name = dir_name
        self.formula = formula
        self.absorption_specie = absorption_specie
        self.edge = edge
        self.structure = structure
        self.spectrum = spectrum
        self.input_parameters = input_parameters
        self.mp_id = mp_id

    def e0_interpolate(self):
        """
        Use scipy interp1d function and relationship between 'relative_energies' and 'energies', calculate e0 of spectrum
        :return: e0
        """

        f = interp1d(self.spectrum['relative_energies'], self.spectrum['energies'])
        self.e0 = f(0).item()

    @classmethod
    def from_file(cls, filename):
        """
        Initiate XANES entry object from file, currently support json file or dictionary object query from database directly
        """

        if filename.endswith('.json'):
            data_entry = json.load(filename)

            return cls(data_entry['dir_name'], data_entry['pretty_formula'], data_entry['absorbing_atom_specie'],
                       data_entry['edge'], data_entry['structure'], data_entry['spectrum'],
                       data_entry['input_parameters'],
                       data_entry['mp_id'])
        elif isinstance(filename, dict):
            data_entry = filename
            return cls(data_entry['dir_name'], data_entry['pretty_formula'], data_entry['absorbing_atom_specie'],
                       data_entry['edge'], data_entry['structure'], data_entry['spectrum'],
                       data_entry['input_parameters'],
                       data_entry['mp_id'])
