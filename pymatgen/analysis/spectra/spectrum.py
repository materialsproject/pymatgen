# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from monty.json import MSONable
import abc

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



