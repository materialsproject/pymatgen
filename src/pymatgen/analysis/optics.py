"""
Perform optical property calculations.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
import scipy.constants as const

if TYPE_CHECKING:
    from pymatgen.io.vasp.outputs import Vasprun


class DielectricAnalysis:
    """
    Class to compute optical properties of materials based on provided energy levels, real and imaginary
    components of the dielectric function. The resulting properties include wavelength, refractive index
    (`n`), extinction coefficient (`k`), reflectivity (`R`), absorptivity (`L`), and transmittance (`T`)
    for up to six configurations.

    This class provides capabilities to compute these properties from supplied dielectric data or parsed
    data from VASP calculation results.

    Attributes:
        energies (numpy.ndarray): Array of energy levels in electron volts (eV).
        eps_real (numpy.ndarray): Array of real parts of the dielectric function for different configurations.
            The structure aligns with the input list provided at initialization.
        eps_imag (numpy.ndarray): Array of imaginary parts of the dielectric function
            matching the structure of eps_real.
        wavelengths (numpy.ndarray): Array of wavelengths (in nanometers) corresponding to the energy levels.
            Computed using Planck's constant and the speed of light.
        n (numpy.ndarray): Array of refractive indices for each energy level and corresponding configuration.
        k (numpy.ndarray): Array of extinction coefficients for each energy level and corresponding configuration.
        R (numpy.ndarray): Array of reflectivity values for each energy level and corresponding configuration.
        L (numpy.ndarray): Array of absorptivity values for each energy level and corresponding configuration.
        T (numpy.ndarray): Array of transmittance values for each energy level and corresponding configuration.
    """

    def __init__(self, energies: list[float], eps_real: list[list[float]], eps_imag: list[list[float]]) -> None:
        """
        Class to compute optical properties of materials based on provided energy levels, real and imaginary components
        of the dielectric function. The resulting properties include wavelength, refractive index (`n`), extinction
        coefficient (`k`), reflectivity (`R`), absorptivity (`L`), and transmittance (`T`).

        Args:
            energies (list[float]): List of energy levels in electron volts (eV).
            eps_real (list[list[float]]): 2D list of real parts of the dielectric function for
                different configurations. The outer list corresponds to different energy levels,
                and the inner list corresponds to configurations.
            eps_imag (list[list[float]]): 2D list of imaginary parts of the dielectric function
                matching the structure of eps_real.
        """
        self.energies = np.asarray(energies)
        self.eps_real = np.asarray(eps_real)
        self.eps_imag = np.asarray(eps_imag)

        # Wavelengths are in nm.
        self.wavelengths = const.h / const.e * const.c / self.energies * 1e9

        self.n = (1 / math.sqrt(2)) * ((self.eps_real**2 + self.eps_imag**2) ** 0.5 + self.eps_real) ** 0.5

        self.k = (1 / math.sqrt(2)) * ((self.eps_real**2 + self.eps_imag**2) ** 0.5 - self.eps_real) ** 0.5

        self.R = ((self.n - 1) ** 2 + self.k**2) / ((self.n + 1) ** 2 + self.k**2)  # reflectivity
        self.L = self.eps_imag / (self.eps_real**2 + self.eps_imag**2)  # absorptivity
        self.T = 1 - self.R - self.L

    @classmethod
    def from_vasprun(cls, vasprun: Vasprun) -> DielectricAnalysis:
        """
        Creates an instance of the DielectricAnalysis class using the dielectric data
        extracted from a VASP calculation parsed by the vasprun object. The dielectric
        data typically includes the energies, real part of the dielectric tensor, and
        imaginary part of the dielectric tensor.

        Args:
            vasprun (Vasprun): Object containing parsed VASP calculation results, specifically
            the dielectric data extracted from the calculation.

        Returns:
            DielectricAnalysis: An instance of the DielectricAnalysis class initialized with the
            dielectric data (energies, real dielectric tensor, and imaginary dielectric
            tensor).
        """
        _energies, _real_diel, _imag_diel = vasprun.dielectric
        return DielectricAnalysis(
            *vasprun.dielectric,
        )
