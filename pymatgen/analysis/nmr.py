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
__maintainer__ = "Shyam Dwaraknath"
__credits__ = "Xiaohui Qu"
__email__ = "shyamd@lbl.gov"
__date__ = "Mar 1, 2018"


class ChemicalShielding(SquareTensor):
    """
    This class extends the SquareTensor to perform extra analysis unique to
    NMR Chemical shielding tensors

    Three notations to describe chemical shielding tensor (RK Harris; Magn. Reson.
    Chem. 2008, 46, 582–598; DOI: 10.1002/mrc.2225) are supported.

    Authors: Shyam Dwaraknath, Xiaohui Qu
    """

    HaeberlenNotation = namedtuple(typename="HaeberlenNotion", field_names="sigma_iso, delta_sigma_iso, zeta, eta")
    MehringNotation = namedtuple(typename="MehringNotation", field_names="sigma_iso, sigma_11, sigma_22, sigma_33")
    MarylandNotation = namedtuple(typename="MarylandNotation", field_names="sigma_iso, omega, kappa")

    def __new__(cls, cs_matrix):
        """
        Create a Chemical Shielding tensor.
        Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.

        Args:
            cs_matrix (1x3 or 3x3 array-like): the 3x3 array-like
                representing the chemical shielding tensor
                or a 1x3 array of the primary sigma values corresponding to the principal axis system
        """
        t_array = np.array(cs_matrix)

        if t_array.shape == (3, ):
            return super(ChemicalShielding, cls).__new__(cls, np.diag(cs_matrix))
        elif t_array.shape == (3, 3):
            return super(ChemicalShielding, cls).__new__(cls, cs_matrix)

    @property
    def principal_axis_system(self):
        """
        Returns a chemical shielding tensor aligned to the principle axis system so that only the 3 diagnol components are non-zero
        """
        return ChemicalShielding(np.diag(np.sort(np.linalg.eigvals(self))))

    @property
    def haeberlen_values(self):
        """
        Returns: the Chemical shielding tensor in Haeberlen Notation
        """
        pas = self.principal_axis_system
        sigma_iso = pas.trace() / 3
        sigmas = np.diag(pas)
        sigmas = sorted(sigmas, key=lambda x: np.abs(x-sigma_iso))
        sigma_yy, sigma_xx, sigma_zz = sigmas
        delta_sigma = sigma_zz - 0.5 * (sigma_xx + sigma_yy)
        zeta = sigma_zz - sigma_iso
        eta = (sigma_yy - sigma_xx) / zeta
        return self.HaeberlenNotation(sigma_iso, delta_sigma, zeta, eta)

    @property
    def mehring_values(self):
        """
        Returns: the Chemical shielding tensor in Mehring Notation
        """
        pas = self.principal_axis_system
        sigma_iso = pas.trace() / 3
        sigma_11, sigma_22, sigma_33 = np.diag(pas)
        return self.MehringNotation(sigma_iso, sigma_11, sigma_22, sigma_33)

    @property
    def maryland_values(self):
        """
        Returns: the Chemical shielding tensor in Maryland Notation
        """
        pas = self.principal_axis_system
        sigma_iso = pas.trace() / 3
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
