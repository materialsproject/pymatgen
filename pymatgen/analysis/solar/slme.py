"""
Calculate spectroscopy limited maximum efficiency (SLME) given dielectric function data.

Forked and adjusted from :
https://github.com/usnistgov/jarvis

References: 1) https://doi.org/10.1021/acs.chemmater.9b02166  &
            2) https://doi.org/10.1103/PhysRevLett.108.068701
"""

from __future__ import annotations

import os
from math import pi

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from scipy.constants import physical_constants, speed_of_light
from scipy.integrate import simps
from scipy.interpolate import interp1d

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.util.due import Doi, due

due.cite(
    Doi("10.1021/acs.chemmater.9b02166"),
    description="Accelerated Discovery of Efficient Solar Cell Materials Using Quantum and Machine-Learning Methods",
)
due.cite(
    Doi("10.1103/PhysRevLett.108.068701"),
    description="Identification of Potential Photovoltaic Absorbers Based on First-Principles "
    "Spectroscopic Screening of Materials",
)


eV_to_recip_cm = 1.0 / (physical_constants["Planck constant in eV s"][0] * speed_of_light * 1e2)


def get_dir_indir_gap(run=""):
    """Get direct and indirect bandgaps for a vasprun.xml."""
    vasp_run = Vasprun(run)
    bandstructure = vasp_run.get_band_structure()
    dir_gap = bandstructure.get_direct_band_gap()
    indir_gap = bandstructure.get_band_gap()["energy"]
    return dir_gap, indir_gap


def matrix_eigvals(matrix):
    """
    Calculate the eigenvalues of a matrix.

    Args:
        matrix (np.array): The matrix to diagonalise.

    Returns:
        (np.array): Array of the matrix eigenvalues.
    """
    eigvals, eigvecs = np.linalg.eig(matrix)
    return eigvals


def to_matrix(xx, yy, zz, xy, yz, xz):
    """
    Convert a list of matrix components to a symmetric 3x3 matrix.
    Inputs should be in the order xx, yy, zz, xy, yz, xz.

    Args:
        xx (float): xx component of the matrix.
        yy (float): yy component of the matrix.
        zz (float): zz component of the matrix.
        xy (float): xy component of the matrix.
        yz (float): yz component of the matrix.
        xz (float): xz component of the matrix.

    Returns:
        (np.array): The matrix, as a 3x3 numpy array.
    """
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])


def parse_dielectric_data(data):
    """
    Convert a set of 2D vasprun formatted dielectric data to
    the eigenvalues of each corresponding 3x3 symmetric numpy matrices.

    Args:
        data (list): length N list of dielectric data. Each entry should be
                     a list of ``[xx, yy, zz, xy, xz, yz ]`` dielectric
                     tensor elements.

    Returns:
        (np.array):  a Nx3 numpy array. Each row contains the eigenvalues
                     for the corresponding row in `data`.
    """
    return np.array([matrix_eigvals(to_matrix(*e)) for e in data])


def absorption_coefficient(dielectric):
    """
    Calculate the optical absorption coefficient from an input set of
    pymatgen vasprun dielectric constant data.

    Args:
        dielectric (list): A list containing the dielectric response function
                           in the pymatgen vasprun format.
                           | element 0: list of energies
                           | element 1: real dielectric tensors, in ``[xx, yy, zz, xy, xz, yz]`` format.
                           | element 2: imaginary dielectric tensors, in ``[xx, yy, zz, xy, xz, yz]`` format.

    Returns:
        (np.array): absorption coefficient using eV as frequency units (cm^-1).
    """
    energies_in_eV = np.array(dielectric[0])
    real_dielectric = parse_dielectric_data(dielectric[1])
    imag_dielectric = parse_dielectric_data(dielectric[2])
    epsilon_1 = np.mean(real_dielectric, axis=1)
    epsilon_2 = np.mean(imag_dielectric, axis=1)
    return (
        energies_in_eV,
        (
            2.0
            * np.sqrt(2.0)
            * pi
            * eV_to_recip_cm
            * energies_in_eV
            * np.sqrt(-epsilon_1 + np.sqrt(epsilon_1**2 + epsilon_2**2))
        ),
    )


def optics(path=""):
    """Helper function to calculate optical absorption coefficient."""
    dirgap, indirgap = get_dir_indir_gap(path)

    run = Vasprun(path, occu_tol=1e-2)
    new_en, new_abs = absorption_coefficient(run.dielectric)
    return (
        np.array(new_en, dtype=np.float64),
        np.array(new_abs, dtype=np.float64),
        dirgap,
        indirgap,
    )


def slme(
    material_energy_for_absorbance_data,
    material_absorbance_data,
    material_direct_allowed_gap,
    material_indirect_gap,
    thickness=50e-6,
    temperature=293.15,
    absorbance_in_inverse_centimeters=False,
    cut_off_absorbance_below_direct_allowed_gap=True,
    plot_current_voltage=False,
):
    """
    Calculate the SLME.

    Args:
        material_energy_for_absorbance_data: energy grid for absorbance data
        material_absorbance_data: absorption coefficient in m^-1
        material_direct_allowed_gap: direct bandgap in eV
        material_indirect_gap: indirect bandgap in eV
        thickness: thickness of the material in m
        temperature: working temperature in K
        absorbance_in_inverse_centimeters: whether the absorbance data is in the unit of cm^-1
        cut_off_absorbance_below_direct_allowed_gap: whether to discard all absorption below bandgap
        plot_current_voltage: whether to plot the current-voltage curve

    Returns:
        The calculated maximum efficiency.

    """
    # Defining constants for tidy equations
    c = constants.c  # speed of light, m/s
    h = constants.h  # Planck's constant J*s (W)
    h_e = constants.h / constants.e  # Planck's constant eV*s
    k = constants.k  # Boltzmann's constant J/K
    k_e = constants.k / constants.e  # Boltzmann's constant eV/K
    e = constants.e  # Coulomb

    # Make sure the absorption coefficient has the right units (m^{-1})
    if absorbance_in_inverse_centimeters:
        material_absorbance_data = material_absorbance_data * 100

    # Load the Air Mass 1.5 Global tilt solar spectrum
    solar_spectrum_data_file = str(os.path.join(os.path.dirname(__file__), "am1.5G.dat"))

    solar_spectra_wavelength, solar_spectra_irradiance = np.loadtxt(
        solar_spectrum_data_file, usecols=[0, 1], unpack=True, skiprows=2
    )

    solar_spectra_wavelength_meters = solar_spectra_wavelength * 1e-9

    delta = material_direct_allowed_gap - material_indirect_gap
    fr = np.exp(-delta / (k_e * temperature))

    # need to convert solar irradiance from Power/m**2(nm) into
    # photon#/s*m**2(nm) power is Watt, which is Joule / s
    # E = hc/wavelength
    # at each wavelength, Power * (wavelength(m)/(h(Js)*c(m/s))) = ph#/s
    solar_spectra_photon_flux = solar_spectra_irradiance * (solar_spectra_wavelength_meters / (h * c))

    # Calculation of total solar power incoming
    power_in = simps(solar_spectra_irradiance, solar_spectra_wavelength)

    # calculation of blackbody irradiance spectra
    # units of W/(m**3), different than solar_spectra_irradiance!!! (This
    # is intentional, it is for convenience)
    blackbody_irradiance = (2.0 * h * c**2 / (solar_spectra_wavelength_meters**5)) * (
        1.0 / ((np.exp(h * c / (solar_spectra_wavelength_meters * k * temperature))) - 1.0)
    )

    # now to convert the irradiance from Power/m**2(m) into photon#/s*m**2(m)
    blackbody_photon_flux = blackbody_irradiance * (solar_spectra_wavelength_meters / (h * c))

    # units of nm
    material_wavelength_for_absorbance_data = ((c * h_e) / (material_energy_for_absorbance_data + 0.00000001)) * 10**9

    # absorbance interpolation onto each solar spectrum wavelength

    # creates cubic spline interpolating function, set up to use end values
    #  as the guesses if leaving the region where data exists
    material_absorbance_data_function = interp1d(
        material_wavelength_for_absorbance_data,
        material_absorbance_data,
        kind="cubic",
        fill_value=(material_absorbance_data[0], material_absorbance_data[-1]),
        bounds_error=False,
    )

    material_interpolated_absorbance = np.zeros(len(solar_spectra_wavelength_meters))
    for i in range(0, len(solar_spectra_wavelength_meters)):
        # Cutting off absorption data below the gap. This is done to deal
        # with VASPs broadening of the calculated absorption data

        if (
            solar_spectra_wavelength[i] < 1e9 * ((c * h_e) / material_direct_allowed_gap)
            or cut_off_absorbance_below_direct_allowed_gap is False
        ):
            material_interpolated_absorbance[i] = material_absorbance_data_function(solar_spectra_wavelength[i])

    absorbed_by_wavelength = 1.0 - np.exp(-2.0 * material_interpolated_absorbance * thickness)

    #  Numerically integrating irradiance over wavelength array
    #  Note: elementary charge, not math e!  # units of A/m**2   W/(V*m**2)
    J_0_r = (
        e
        * np.pi
        * simps(
            blackbody_photon_flux * absorbed_by_wavelength,
            solar_spectra_wavelength_meters,
        )
    )

    J_0 = J_0_r / fr

    J_sc = e * simps(solar_spectra_photon_flux * absorbed_by_wavelength, solar_spectra_wavelength)

    def J(V):
        return J_sc - J_0 * (np.exp(e * V / (k * temperature)) - 1.0)

    def power(V):
        return J(V) * V

    test_voltage = 0
    voltage_step = 0.001
    while power(test_voltage + voltage_step) > power(test_voltage):
        test_voltage += voltage_step

    max_power = power(test_voltage)

    # Calculate the maximized efficiency
    efficiency = max_power / power_in

    if plot_current_voltage:
        V = np.linspace(0, 2, 200)
        plt.plot(V, J(V))
        plt.plot(V, power(V), linestyle="--")
        plt.savefig("pp.png")
        plt.close()

    return 100.0 * efficiency
