# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
A module for NMR analysis
"""

from pymatgen.core.tensors import SquareTensor
from collections import namedtuple

from pymatgen.core.units import FloatWithUnit

from pymatgen.core.periodic_table import Specie
from pymatgen.core.structure import Site

import numpy as np


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

    HaeberlenNotation = namedtuple("HaeberlenNotion", "sigma_iso, delta_sigma_iso, zeta, eta")
    MehringNotation = namedtuple("MehringNotation", "sigma_iso, sigma_11, sigma_22, sigma_33")
    MarylandNotation = namedtuple("MarylandNotation", "sigma_iso, omega, kappa")

    def __new__(cls, cs_matrix, vscale=None):
        """
        Create a Chemical Shielding tensor.
        Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.

        Args:
            cs_matrix (1x3 or 3x3 array-like): the 3x3 array-like
                representing the chemical shielding tensor
                or a 1x3 array of the primary sigma values corresponding
                to the principal axis system
            vscale (6x1 array-like): 6x1 array-like scaling the
                voigt-notation vector with the tensor entries
        """
        t_array = np.array(cs_matrix)

        if t_array.shape == (3,):
            return super().__new__(cls, np.diag(cs_matrix), vscale)
        elif t_array.shape == (3, 3):
            return super().__new__(cls, cs_matrix, vscale)

    @property
    def principal_axis_system(self):
        """
        Returns a chemical shielding tensor aligned to the principle axis system
        so that only the 3 diagnol components are non-zero
        """
        return ChemicalShielding(np.diag(np.sort(np.linalg.eigvals(self.symmetrized))))

    @property
    def haeberlen_values(self):
        """
        Returns: the Chemical shielding tensor in Haeberlen Notation
        """
        pas = self.principal_axis_system
        sigma_iso = pas.trace() / 3
        sigmas = np.diag(pas)
        sigmas = sorted(sigmas, key=lambda x: np.abs(x - sigma_iso))
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
        """
        Initialize from Maryland notation.

        Args:
            sigma_iso ():
            omega ():
            kappa ():

        Returns:
            ChemicalShielding
        """
        sigma_22 = sigma_iso + kappa * omega / 3.0
        sigma_11 = (3.0 * sigma_iso - omega - sigma_22) / 2.0
        sigma_33 = 3.0 * sigma_iso - sigma_22 - sigma_11
        return cls(np.diag([sigma_11, sigma_22, sigma_33]))


class ElectricFieldGradient(SquareTensor):
    """
    This class extends the SquareTensor to perform extra analysis unique to
    NMR Electric Field Gradient tensors in units of V/Angstrom^2

    Authors: Shyam Dwaraknath, Xiaohui Qu
    """

    def __new__(cls, efg_matrix, vscale=None):
        """
        Create a Chemical Shielding tensor.
        Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.

        Args:
            efg_matrix (1x3 or 3x3 array-like): the 3x3 array-like
                representing the electric field tensor
                or a 1x3 array of the primary values corresponding
                to the principal axis system
            vscale (6x1 array-like): 6x1 array-like scaling the
                voigt-notation vector with the tensor entries
        """
        t_array = np.array(efg_matrix)

        if t_array.shape == (3,):
            return super().__new__(cls, np.diag(efg_matrix), vscale)
        elif t_array.shape == (3, 3):
            return super().__new__(cls, efg_matrix, vscale)

    @property
    def principal_axis_system(self):
        """
        Returns a electric field gradient tensor aligned to the principle axis system so that only the 3 diagnol
        components are non-zero
        """
        return ElectricFieldGradient(np.diag(np.sort(np.linalg.eigvals(self))))

    @property
    def V_xx(self):
        """
        Returns: First diagonal element
        """
        diags = np.diag(self.principal_axis_system)
        return sorted(diags, key=np.abs)[0]

    @property
    def V_yy(self):
        """
        Returns: Second diagonal element
        """
        diags = np.diag(self.principal_axis_system)
        return sorted(diags, key=np.abs)[1]

    @property
    def V_zz(self):
        """
        Returns: Third diagonal element
        """
        diags = np.diag(self.principal_axis_system)
        return sorted(diags, key=np.abs)[2]

    @property
    def asymmetry(self):
        """
        Asymmetry of the electric field tensor defined as:
            (V_yy - V_xx)/V_zz
        """
        diags = np.diag(self.principal_axis_system)
        V = sorted(diags, key=np.abs)
        return np.abs((V[1] - V[0]) / V[2])

    def coupling_constant(self, specie):
        """
        Computes the couplling constant C_q as defined in:
            Wasylishen R E, Ashbrook S E, Wimperis S. NMR of quadrupolar nuclei
            in solid materials[M]. John Wiley & Sons, 2012. (Chapter 3.2)

        C_q for a specific atom type for this electric field tensor:
                C_q=e*Q*V_zz/h
            h: planck's constant
            Q: nuclear electric quadrupole moment in mb (millibarn
            e: elementary proton charge

        Args:
            specie: flexible input to specify the species at this site.
                    Can take a isotope or element string, Specie object,
                    or Site object

        Return:

            the coupling constant as a FloatWithUnit in MHz
        """
        planks_constant = FloatWithUnit(6.62607004E-34, "m^2 kg s^-1")
        Vzz = FloatWithUnit(self.V_zz, "V ang^-2")
        e = FloatWithUnit(-1.60217662E-19, "C")

        # Convert from string to Specie object
        if isinstance(specie, str):
            # isotope was provided in string format
            if len(specie.split("-")) > 1:
                isotope = str(specie)
                specie = Specie(specie.split("-")[0])
                Q = specie.get_nmr_quadrupole_moment(isotope)
            else:
                specie = Specie(specie)
                Q = specie.get_nmr_quadrupole_moment()
        elif isinstance(specie, Site):
            specie = specie.specie
            Q = specie.get_nmr_quadrupole_moment()
        elif isinstance(specie, Specie):
            Q = specie.get_nmr_quadrupole_moment()
        else:
            raise ValueError("Invalid speciie provided for quadrupolar coupling constant calcuations")

        return (e * Q * Vzz / planks_constant).to("MHz")
