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


class XANES(Spectrum, MSONable):
    """
    Basic XANES object.

    Args:
        energies: A sequence of x-ray energies in eV
        mu: A sequence of mu(E)

    .. attribute: energies
        The sequence of energies
    .. attribute: mu
        The sequence of mu(E)
    .. attribute: absorption_specie
        The absorption_species of the spectrum
    """

    def __init__(self, structure, energies, mu, absorption_specie, edge):
        try:
            self.structure = Structure.from_dict(structure)
        except:
            self.structure = structure
        self.energies = np.array(energies)
        self.mu = np.array(mu)
        self.absorption_specie = absorption_specie
        self.edge = edge
        self.type = 'XANES'
        self.spectrum = np.column_stack((self.energies, self.mu))

    def find_e0(self):
        """
        Use the maximum gradient to find e0
        """
        self.e0 = self.energies[np.argmax(np.gradient(self.mu) / np.gradient(self.energies))]

    def spectrum_norm(self):
        """
        Normalize the peak intentsity according to the cumulative sum of the intensity.
        Therefore, peaks retain properties required from probability mass function
        """
        return np.column_stack((self.energies, (self.mu / self.mu.sum())))

    @classmethod
    def from_dict(cls, d):
        """
        Return XANES object from dict representation of XANES
        """
        return XANES(d['structure'], d['energies'], d['mu'],
                     d['absorption_specie'], d['edge'])

    def as_dict(self):
        """
        Json-serializable dict representation of XANES
        """
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "energeis": list(self.energies),
                "mu": list(self.mu),
                "absorption_specie": list(self.absorption_specie),
                "edge": self.edge
                }
