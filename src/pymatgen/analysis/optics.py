"""
Perform optical property calculations.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import scipy.constants as const


class DielectricAnalysis:
    """
    Represents the analysis of dielectric properties derived from energy-dependent
    dielectric functions. This includes calculations for wavelength, refractive index,
    extinction coefficient, reflectivity, absorptivity, and transmittance based on
    real and imaginary parts of the dielectric tensors.

    This class is designed to process dielectric tensors (both real and imaginary)
    across different energy levels and provide detailed tabular data. It is useful
    for understanding optical and electronic properties of materials.

    :ivar data: DataFrame containing calculated properties such as wavelength,
        refractive index, extinction coefficient, reflectivity, absorptivity, and
        transmittance for each tensor component.
    :type data: pandas.DataFrame
    """

    def __init__(self, energies: list[float], eps_real: list[list[float]], eps_imag: list[list[float]]) -> None:
        """
        Class to compute optical properties of materials based on provided energy levels, real and imaginary components
        of the dielectric function. The resulting properties include wavelength, refractive index (`n`), extinction
        coefficient (`k`), reflectivity (`R`), absorptivity (`L`), and transmittance (`T`) for up to six configurations.

        :param energies: List of energy levels in electron volts (eV).
        :type energies: list[float]
        :param eps_real: 2D list of real parts of the dielectric function for different configurations.
                         The outer list corresponds to different energy levels, and the inner list corresponds to
                         configurations.
        :type eps_real: list[list[float]]
        :param eps_imag: 2D list of imaginary parts of the dielectric function matching the structure of eps_real.
        :type eps_imag: list[list[float]]
        """

        data = pd.DataFrame(
            {
                "energy": energies,
            }
        )
        data["wavelength"] = (
            const.physical_constants["Planck constant in eV s"][0]
            * const.physical_constants["speed of light in vacuum"][0]
            / data["energy"]
            * 1e9
        )
        eps_real = np.array(eps_real)
        eps_imag = np.array(eps_imag)

        for i in range(6):
            data[f"eps_{i}_real"] = eps_real[:, i]  # type: ignore[call-overload]
            data[f"eps_{i}_imag"] = eps_imag[:, i]  # type: ignore[call-overload]

            data[f"n_{i}"] = (1 / math.sqrt(2)) * (
                (data[f"eps_{i}_real"] ** 2 + data[f"eps_{i}_imag"] ** 2) ** 0.5 + data[f"eps_{i}_real"]
            ) ** 0.5
            data[f"k_{i}"] = (1 / math.sqrt(2)) * (
                (data[f"eps_{i}_real"] ** 2 + data[f"eps_{i}_imag"] ** 2) ** 0.5 - data[f"eps_{i}_real"]
            ) ** 0.5

            # reflectivity
            data[f"R_{i}"] = ((data[f"n_{i}"] - 1) ** 2 + data[f"k_{i}"] ** 2) / (
                (data[f"n_{i}"] + 1) ** 2 + data[f"k_{i}"] ** 2
            )
            # absorptivity
            data[f"L_{i}"] = data[f"eps_{i}_imag"] / (data[f"eps_{i}_real"] ** 2 + data[f"eps_{i}_imag"] ** 2)
            data[f"T_{i}"] = 1 - data[f"R_{i}"] - data[f"L_{i}"]

        self.data = data

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
