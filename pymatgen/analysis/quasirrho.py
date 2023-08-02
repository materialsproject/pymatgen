# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
A module to calculate free energies using the Quasi-Rigid Rotor Harmonic
 Oscillator approximation. Modified from a script by Steven Wheeler.
See: Grimme, S. Chem. Eur. J. 2012, 18, 9955
"""

from __future__ import annotations

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Alex Epstein"
__email__ = "aepstein@lbl.gov"
__date__ = "August 1, 2023"
__credits__ = "Ryan Kingsbury, Steven Wheeler, Trevor Seguin, Evan Spotte-Smith"

from math import isclose
from typing import TYPE_CHECKING

import numpy as np
import scipy.constants as const

from pymatgen.core.units import kb as kb_ev

if TYPE_CHECKING:
    from pymatgen.core import Molecule
    from pymatgen.io.gaussian import GaussianOutput
    from pymatgen.io.qchem.outputs import QCOutput

# Define useful constants
kb = kb_ev * const.eV  # Pymatgen kb [J/K]
c = const.speed_of_light * 100  # [cm/s]
h = const.h  # Planck's constant [J.s]
R = const.R / const.calorie  # Ideal gas constant [cal/mol/K]
R_ha = const.R / const.value("Hartree energy") / const.Avogadro  # Ideal gas
# constant [Ha/K]
R_volume = const.R / const.atm * 1000  # Ideal gas constant [L.atm.K^-1.mol^-1]

# Define useful conversion factors
amu_to_kg = const.value("atomic mass unit-kilogram relationship")  # AMU to kg
# kcal2hartree = 0.0015936  # kcal/mol to hartree/mol
kcal2hartree = 1000 * const.calorie / const.value("Hartree energy") / const.Avogadro


def get_avg_mom_inertia(mol):
    """
    Calculate the average moment of inertia of a molecule

    Args:
        mol (Molecule): Pymatgen Molecule

    Returns:
        int, list: average moment of inertia, eigenvalues of the inertia tensor
    """
    centered_mol = mol.get_centered_molecule()
    inertia_tensor = np.zeros((3, 3))
    for site in centered_mol:
        c = site.coords
        wt = site.specie.atomic_mass
        for i in range(3):
            inertia_tensor[i, i] += wt * (c[(i + 1) % 3] ** 2 + c[(i + 2) % 3] ** 2)
        for i, j in [(0, 1), (1, 2), (0, 2)]:
            inertia_tensor[i, j] += -wt * c[i] * c[j]
            inertia_tensor[j, i] += -wt * c[j] * c[i]

    inertia_eigenvals = np.multiply(np.linalg.eig(inertia_tensor)[0], amu_to_kg * 1e-20).tolist()  # amuangs^2 to kg m^2

    iav = np.average(inertia_eigenvals)

    return iav, inertia_eigenvals


class QuasiRRHO:
    """
    Class to calculate thermochemistry using Grimme's Quasi-RRHO approximation.
    All outputs are in atomic units, e.g. energy outputs are in Hartrees.
    Citation: Grimme, S. Chemistry - A European Journal 18, 9955-9964 (2012).

    Attributes:
        temp (float): Temperature [K]
        press (float): Pressure [Pa]
        conc (float): Solvent concentration. Assumes 1M unless specified [M]
        v0 (float): Cutoff frequency for Quasi-RRHO method [1/cm]
        entropy_quasiRRHO (float): Quasi-RRHO entropy [Ha/K]
        entropy_ho (float): Total entropy calculated with a harmonic
            oscillator approximation for the vibrational entropy [Ha/K]
        h_corrected (float): Thermal correction to the enthalpy [Ha]
        free_energy_quasiRRHO (float): Quasi-RRHO free energy [Ha]
        concentration_corrected_g_quasiRRHO (float): Quasi-RRHO free energy
            with a standard state correction for the solvent concentration [Ha]
        free_energy_ho (float): Free energy calculated without the Quasi-RRHO
            method, i.e. with a harmonic oscillator approximation for the
            vibrational entropy [Ha]

    """

    def __init__(
        self,
        mol: Molecule,
        frequencies: list,
        energy: float,
        mult: int,
        sigma_r=1,
        temp=298.15,
        press=101317,
        conc=1,
        v0=100,
    ):
        """
        Args:
            mol (Molecule): Pymatgen molecule
            frequencies (list): List of frequencies (float) [cm^-1]
            energy (float): Electronic energy [Ha]
            mult (int): Spin muliplicity
            sigma_r (int): Rotational symmetry number
            temp (float): Temperature [K]
            press (float): Pressure [Pa]
            conc (float): Solvent concentration [M]
            v0 (float): Cutoff frequency for Quasi-RRHO method [cm^-1]
        """
        # TO-DO: calculate sigma_r with PointGroupAnalyzer
        # and/or edit Gaussian and QChem io to parse for sigma_r

        self.temp = temp
        self.press = press
        self.conc = conc
        self.v0 = v0

        self.entropy_quasiRRHO = None  # Ha/K
        self.free_energy_quasiRRHO = None  # Ha
        self.concentration_corrected_g_quasiRRHO = None  # Ha
        self.h_corrected = None  # Ha

        self.entropy_ho = None  # Ha/K
        self.free_energy_ho = None  # Ha

        self._get_quasirrho_thermo(mol=mol, mult=mult, frequencies=frequencies, elec_energy=energy, sigma_r=sigma_r)

    @classmethod
    def from_GaussianOutput(cls, output: GaussianOutput, **kwargs):
        """

        Args:
            output (GaussianOutput): Pymatgen GaussianOutput object

        Returns:
            QuasiRRHO: QuasiRRHO class instantiated from a Gaussian Output
        """
        mult = output.spin_multiplicity
        elec_e = output.final_energy
        mol = output.final_structure
        vib_freqs = [f["frequency"] for f in output.frequencies[-1]]
        return cls(mol=mol, frequencies=vib_freqs, energy=elec_e, mult=mult, **kwargs)

    @classmethod
    def from_QCOutput(cls, output: QCOutput, **kwargs):
        """

        Args:
            output (QCOutput): Pymatgen QCOutput object

        Returns:
            QuasiRRHO: QuasiRRHO class instantiated from a QChem Output

        """
        mult = output.data["multiplicity"]
        elec_e = output.data["SCF_energy_in_the_final_basis_set"]
        if output.data["optimization"]:
            mol = output.data["molecule_from_last_geometry"]
        else:
            mol = output.data["initial_molecule"]
        frequencies = output.data["frequencies"]

        return cls(mol=mol, frequencies=frequencies, energy=elec_e, mult=mult, **kwargs)

    @due.dcite(
        Doi("10.1002/chem.201200497"),
        description="Supramolecular Binding Thermodynamics by Dispersion-Corrected Density Functional Theory",
    )
    def _get_quasirrho_thermo(self, mol, mult, sigma_r, frequencies, elec_energy):
        """
        Calculate Quasi-RRHO thermochemistry
        Args:
            mol (Molecule): Pymatgen molecule
            mult (int): Spin multiplicity
            sigma_r (int): Rotational symmetry number
            frequencies (list): List of frequencies [cm^-1]
            elec_energy (float): Electronic energy [Ha]
        """
        # Calculate mass in kg
        mass = 0
        for site in mol.sites:
            mass += site.specie.atomic_mass
        mass *= amu_to_kg

        # Calculate vibrational temperatures
        vib_temps = []
        for f in frequencies:
            if f > 0:
                vib_temps.append(f * c * h / kb)

        # Translational component of entropy and energy
        qt = (2 * np.pi * mass * kb * self.temp / (h * h)) ** (3 / 2) * kb * self.temp / self.press
        st = R * (np.log(qt) + 5 / 2)
        et = 3 * R * self.temp / 2

        # Electronic component of Entropy
        se = R * np.log(mult)

        # Get properties related to rotational symmetry. Bav is average moment of inertia
        Bav, i_eigen = get_avg_mom_inertia(mol)

        # Check if linear
        coords = mol.cart_coords
        v0 = coords[1] - coords[0]
        linear = True
        for coord in coords[1:]:
            theta = abs(np.dot(coord - coords[0], v0) / np.linalg.norm(coord - coords[0]) / np.linalg.norm(v0))
            if not isclose(theta, 1, abs_tol=1e-4):
                linear = False

        # Rotational component of Entropy and Energy
        if linear:
            i = np.amax(i_eigen)
            qr = 8 * np.pi**2 * i * kb * self.temp / (sigma_r * (h * h))
            sr = R * (np.log(qr) + 1)
            er = R * self.temp
        else:
            rot_temps = [h**2 / (np.pi**2 * kb * 8 * i) for i in i_eigen]
            qr = np.sqrt(np.pi) / sigma_r * self.temp ** (3 / 2) / np.sqrt(rot_temps[0] * rot_temps[1] * rot_temps[2])
            sr = R * (np.log(qr) + 3 / 2)
            er = 3 * R * self.temp / 2

        # Vibrational component of Entropy and Energy
        ev = 0
        sv_quasiRRHO = 0
        sv = 0

        for vt in vib_temps:
            ev += vt * (1 / 2 + 1 / (np.exp(vt / self.temp) - 1))
            sv_temp = vt / (self.temp * (np.exp(vt / self.temp) - 1)) - np.log(1 - np.exp(-vt / self.temp))
            sv += sv_temp

            mu = h / (8 * np.pi**2 * vt * c)
            mu_prime = mu * Bav / (mu + Bav)
            srotor = 1 / 2 + np.log(np.sqrt(8 * np.pi**3 * mu_prime * kb * self.temp / h**2))
            weight = 1 / (1 + (self.v0 / vt) ** 4)
            sv_quasiRRHO += weight * sv_temp + (1 - weight) * srotor

        sv_quasiRRHO *= R
        sv *= R
        ev *= R
        etot = (et + er + ev) * kcal2hartree / 1000
        self.h_corrected = etot + R * self.temp * kcal2hartree / 1000
        molarity_corr = R_ha * self.temp * np.log(R_volume * self.temp * self.conc)
        self.entropy_ho = st + sr + sv + se
        self.free_energy_ho = elec_energy + self.h_corrected - (self.temp * self.entropy_ho * kcal2hartree / 1000)
        self.entropy_quasiRRHO = st + sr + sv_quasiRRHO + se
        self.free_energy_quasiRRHO = (
            elec_energy + self.h_corrected - (self.temp * self.entropy_quasiRRHO * kcal2hartree / 1000)
        )
        self.concentration_corrected_g_quasiRRHO = self.free_energy_quasiRRHO + molarity_corr
