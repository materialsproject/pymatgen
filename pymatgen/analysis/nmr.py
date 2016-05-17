# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

from collections import namedtuple
from monty.json import MSONable

"""
A module to perform NMR data analysis/processing.
"""


__author__ = "Xiaohui Qu"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Xiaohui Qu"
__email__ = "xhqu1981@gmail.com"
__date__ = "Apr 17, 2016"


class NMRChemicalShiftNotation(MSONable):
    """
    Helper class to convert between different chemical shift conventions
    internally using the Mehring notation. Note that this is different than the
    default notion adopted by VASP which is the Maryland notation.

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
                                   field_names="sigma_iso, delta_sigma, zeta, eta")
    MehringNotation = namedtuple(typename="MehringNotation",
                                 field_names="sigma_iso, sigma_11, sigma_22, sigma_33")
    MarylandNotation = namedtuple(typename="MarylandNotation",
                                  field_names="sigma_iso, omega, kappa")

    def __init__(self, sigma_1, sigma_2, sigma_3):
        sigmas = sorted([sigma_1, sigma_2, sigma_3])
        self.sigma_11, self.sigma_22, self.sigma_33 = sigmas

    @property
    def haeberlen_values(self):
        """
        Returns: the Chemical shift tensor in Haeberlen Notation
        """
        sigma_iso = (self.sigma_11 + self.sigma_22 + self.sigma_33) / 3.0
        h_order_sigmas = sorted([self.sigma_11, self.sigma_22, self.sigma_33],
                                key=lambda x: abs(x - sigma_iso),
                                reverse=True)
        sigma_zz, sigma_xx, sigma_yy = h_order_sigmas
        delta_sigma = sigma_zz - 0.5 * (sigma_xx + sigma_yy)
        zeta = sigma_zz - sigma_iso
        assert abs(delta_sigma - 1.5 * zeta) < 1.0E-5
        eta = (sigma_yy - sigma_xx) / zeta
        return self.HaeberlenNotation(sigma_iso, delta_sigma, zeta, eta)

    @property
    def mehring_values(self):
        """
        Returns: the Chemical shift tensor in Mehring Notation
        """
        sigma_iso = (self.sigma_11 + self.sigma_22 + self.sigma_33) / 3.0
        return self.MehringNotation(sigma_iso, self.sigma_11,
                                    self.sigma_22, self.sigma_33)

    @property
    def maryland_values(self):
        """
        Returns: the Chemical shift tensor in Maryland Notation
        """
        sigma_iso = (self.sigma_11 + self.sigma_22 + self.sigma_33) / 3.0
        omega = self.sigma_33 - self.sigma_11
        # There is a typo in equation 20 from Magn. Reson. Chem. 2008, 46, 582–598, the sign is wrong.
        # There correct order is presented in Solid State Nucl. Magn. Reson. 1993, 2, 285-288.
        kappa = 3.0 * (self.sigma_22 - sigma_iso) / omega
        return self.MarylandNotation(sigma_iso, omega, kappa)

    @classmethod
    def from_maryland_notation(cls, sigma_iso, omega, kappa):
        sigma_22 = sigma_iso + kappa * omega / 3.0
        sigma_11 = (3.0 * sigma_iso - omega - sigma_22) / 2.0
        sigma_33 = 3.0 * sigma_iso - sigma_22 - sigma_11
        return cls(sigma_11, sigma_22, sigma_33)

    def as_dict(self):
        d = {"sigma_11": self.sigma_11,
             "sigma_22": self.sigma_22,
             "sigma_33": self.sigma_33}
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d["sigma_11"], d["sigma_22"], d["sigma_33"])
