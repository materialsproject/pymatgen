# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
from pymatgen.analysis.elasticity.tensors import SquareTensor
from collections import namedtuple
import numpy as np
"""
A module for NMR analysis
"""


__author__ = "Shyam Dwaraknath"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.2"
__maintainer__ = "Xiaohui Qu"
__credits__ = "Xiaohui Qu"
__email__ = "shyamd@lbl.gov"
__date__ = "Mar 1, 2018"


class ChemicalShift(SquareTensor):
    """
    This class extends the SquareTensor to perform extra analysis unique to
    NMR Chemical shift tensors
    
    Three notations to describe chemical shift tensor (RK Harris; Magn. Reson.
    Chem. 2008, 46, 582–598; DOI: 10.1002/mrc.2225) are supported.

    Args:
        sigma_1 (float): chemical shift tensor principle component 1
        sigma_2 (float): chemical shift tensor principle component 2
        sigma_3 (float): chemical shift tensor principle component 3

    .. attribute:: sigma_11, simga_22, sigma33
        principle components in Mehring notation

    Authors: Xiaohui Qu
    """

    HaeberlenNotation = namedtuple(typename="HaeberlenNotion",
                                   field_names="sigma_iso, delta_sigma_iso, zeta, eta")
    MehringNotation = namedtuple(typename="MehringNotation",
                                 field_names="sigma_iso, sigma_11, sigma_22, sigma_33")
    MarylandNotation = namedtuple(typename="MarylandNotation",
                                  field_names="sigma_iso, omega, kappa")


    def __new__(cls, cs_matrix):
        """
        Create a Chemical Shift tensor.
        Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.

        Args:
            cs_matrix (1x3 or 3x3 array-like): the 3x3 array-like
                representing the chemical shift tensor
                or a 1x3 array of the primary sigma values corresponding to the principal axis system
        """
        t_array = np.array(cs_matrix)

        if t_array.shape == (3,):
            return super(ChemicalShift, cls).__new__(cls, np.diag(cs_matrix))
        elif t_array.shape == (3,3):
            return super(ChemicalShift, cls).__new__(cls, cs_matrix)

    @property
    def principal_axis_system(self):
        """
        Returns a chemical shift tensor aligned to the principle axis system so that only the 3 diagnol components are non-zero
        """
        return ChemicalShift(np.diag(np.sort(np.linalg.eigvals(self))))

    @property
    def haeberlen_values(self):
        """
        Returns: the Chemical shift tensor in Haeberlen Notation
        """
        pas = self.principal_axis_system
        sigma_iso = pas.trace()/3
        sigma_indexes = np.argsort(np.abs(np.diag(pas) - sigma_iso))
        sigmas = np.diag(pas)[sigma_indexes]
        sigma_yy = sigmas[0]
        sigma_xx = sigmas[1]
        sigma_zz = sigmas[2]
        delta_sigma = sigma_zz - 0.5 * (sigma_xx + sigma_yy)
        zeta = sigma_zz - sigma_iso
        eta = (sigma_yy - sigma_xx) / zeta
        return self.HaeberlenNotation(sigma_iso, delta_sigma, zeta, eta)

    @property
    def mehring_values(self):
        """
        Returns: the Chemical shift tensor in Mehring Notation
        """
        pas = self.principal_axis_system
        sigma_iso = pas.trace()/3
        sigma_11 = np.diag(pas)[0]
        sigma_22 = np.diag(pas)[1]
        sigma_33 = np.diag(pas)[2]
        return self.MehringNotation(sigma_iso, sigma_11, sigma_22, sigma_33)

    @property
    def maryland_values(self):
        """
        Returns: the Chemical shift tensor in Maryland Notation
        """
        pas = self.principal_axis_system
        sigma_iso = pas.trace()/3
        omega = np.diag(pas)[2] - np.diag(pas)[0]
        # There is a typo in equation 20 from Magn. Reson. Chem. 2008, 46, 582–598, the sign is wrong.
        # There correct order is presented in Solid State Nucl. Magn. Reson. 1993, 2, 285-288.
        kappa = 3.0 * (np.diag(pas)[1] - sigma_iso) / omega
        return self.MarylandNotation(sigma_iso, omega, kappa)

    @classmethod
    def from_maryland_notation(cls, sigma_iso, omega, kappa):
        sigma_22 = sigma_iso + kappa * omega / 3.0
        sigma_11 = (3.0 * sigma_iso - omega - sigma_22) / 2.0
        sigma_33 = 3.0 * sigma_iso - sigma_22 - sigma_11
        return cls(np.diag([sigma_11, sigma_22, sigma_33]))
