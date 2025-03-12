"""
Perform optical property calculations.
"""

from __future__ import annotations

import math

import numpy as np
import scipy.constants as const


class DielectricAnalysis:
    """
    Class to compute optical properties of materials based on provided energy levels, real and imaginary
    components of the dielectric function. The resulting properties include wavelength, refractive index
    (`n`), extinction coefficient (`k`), reflectivity (`R`), absorptivity (`L`), and transmittance (`T`)
    for up to six configurations.

    This class provides capabilities to compute these properties from supplied dielectric data or parsed
    data from VASP calculation results.

    :ivar energies: Array of energy levels in electron volts (eV).
    :type energies: numpy.ndarray
    :ivar eps_real: Array of real parts of the dielectric function for different configurations.
                    The structure aligns with the input list provided at initialization.
    :type eps_real: numpy.ndarray
    :ivar eps_imag: Array of imaginary parts of the dielectric function matching the structure of eps_real.
    :type eps_imag: numpy.ndarray
    :ivar wavelengths: Array of wavelengths (in nanometers) corresponding to the energy levels.
                       Computed using Planck's constant and the speed of light.
    :type wavelengths: numpy.ndarray
    :ivar n: Array of refractive indices for each energy level and corresponding configuration.
    :type n: numpy.ndarray
    :ivar k: Array of extinction coefficients for each energy level and corresponding configuration.
    :type k: numpy.ndarray
    :ivar R: Array of reflectivity values for each energy level and corresponding configuration.
    :type R: numpy.ndarray
    :ivar L: Array of absorptivity values for each energy level and corresponding configuration.
    :type L: numpy.ndarray
    :ivar T: Array of transmittance values for each energy level and corresponding configuration.
    :type T: numpy.ndarray
    """

    def __init__(self, energies: list[float], eps_real: list[list[float]], eps_imag: list[list[float]]) -> None:
        """
        Class to compute optical properties of materials based on provided energy levels, real and imaginary components
        of the dielectric function. The resulting properties include wavelength, refractive index (`n`), extinction
        coefficient (`k`), reflectivity (`R`), absorptivity (`L`), and transmittance (`T`).

        :param energies: List of energy levels in electron volts (eV).
        :type energies: list[float]
        :param eps_real: 2D list of real parts of the dielectric function for different configurations.
                         The outer list corresponds to different energy levels, and the inner list corresponds to
                         configurations.
        :type eps_real: list[list[float]]
        :param eps_imag: 2D list of imaginary parts of the dielectric function matching the structure of eps_real.
        :type eps_imag: list[list[float]]
        """
        self.energies = np.array(energies)
        self.eps_real = np.array(eps_real)
        self.eps_imag = np.array(eps_imag)

        # Wavelengths are in nm.
        self.wavelengths = const.h / const.e * const.c / self.energies * 1e9

        self.n = (1 / math.sqrt(2)) * ((self.eps_real**2 + self.eps_imag**2) ** 0.5 + self.eps_real) ** 0.5

        self.k = (1 / math.sqrt(2)) * ((self.eps_real**2 + self.eps_imag**2) ** 0.5 - self.eps_real) ** 0.5

        self.R = ((self.n - 1) ** 2 + self.k**2) / ((self.n + 1) ** 2 + self.k**2)  # reflectivity
        self.L = self.eps_imag / (self.eps_real**2 + self.eps_imag**2)  # absorptivity
        self.T = 1 - self.R - self.L

    @classmethod
    def from_vasprun(cls, vasprun):
        """
        Creates an instance of the DielectricAnalysis class using the dielectric data
        extracted from a VASP calculation parsed by the vasprun object. The dielectric
        data typically includes the energies, real part of the dielectric tensor, and
        imaginary part of the dielectric tensor.

        :param vasprun: Object containing parsed VASP calculation results, specifically
            the dielectric data extracted from the calculation.
        :type vasprun: Vasprun
        :return: An instance of the DielectricAnalysis class initialized with the
            dielectric data (energies, real dielectric tensor, and imaginary dielectric
            tensor).
        :rtype: DielectricAnalysis
        """
        energies, real_diel, imag_diel = vasprun.dielectric
        return DielectricAnalysis(
            *vasprun.dielectric,
        )
