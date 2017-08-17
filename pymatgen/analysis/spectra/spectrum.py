# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from monty.json import MSONable
import abc
from scipy.interpolate import interp1d
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen import Structure
import json
import numpy as np

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
    def spectrum_structure(self):
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
    def spectrum_data(self):
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

    def __init__(self, dir_name, formula, absorbing_specie, edge, structure,
                 spectrum, input_parameters, absorbing_atom_index, mp_id, type = "XANES"):
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
        self.absorbing_specie = absorbing_specie
        self.edge = edge
        try:
            self.structure = Structure.from_dict(structure)
        except:
            self.structure = structure
        self.spectrum = np.array(spectrum)
        self.input_parameters = input_parameters
        self.absorbing_atom_index = absorbing_atom_index
        self.mp_id = mp_id
        self.type = type

    def e0_interpolate(self):
        """
        Use scipy interp1d function and relationship between 'relative_energies' and 'energies', calculate e0 of spectrum
        :return: e0
        """

        f = interp1d(self.spectrum[:,1], self.spectrum[:,0])
        self.e0 = f(0).item()

    def site_multiplicity(self):
        """
        Use SpacegroupAnalysis and SymmetrizedStructure to find multiplicity number of absorbing site in Structure,
            i.e. number of equivalent sites in Structure w.r.t absorbing site
        """
        absorbing_structure = Structure.from_dict(self.structure)

        spaceg_analysis = SpacegroupAnalyzer(absorbing_structure)
        sym_structure = spaceg_analysis.get_symmetrized_structure()
        equivalent_sites = sym_structure.find_equivalent_sites(absorbing_structure[self.absorbing_atom_index])
        self.equivalent_sites = equivalent_sites
        self.absorber_multiplicity = len(equivalent_sites)

    @classmethod
    def from_file(cls, filename):
        """
        Initiate XANES entry object from file, currently support json file or dictionary object query from database directly
        """

        if isinstance(filename, dict):
            data_entry = filename
            return cls(data_entry['dir_name'], data_entry['pretty_formula'], data_entry['absorbing_atom_specie'],
                       data_entry['edge'], data_entry['structure'], data_entry['spectrum'],
                       data_entry['input_parameters'], data_entry['metadata']['absorbing_atom_index'],
                       data_entry['mp_id'])

        elif filename.endswith('.json'):
            with open(filename, 'r') as f:
                data_entry = json.load(f)

            return cls(data_entry['dir_name'], data_entry['pretty_formula'], data_entry['absorbing_atom_specie'],
                       data_entry['edge'], data_entry['structure'], data_entry['spectrum'],
                       data_entry['input_parameters'], data_entry['metadata']['absorbing_atom_index'],
                       data_entry['mp_id'])

        else:
            raise Exception('Unknown data type')
