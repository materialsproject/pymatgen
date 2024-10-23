# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
A module to calculate free energies using the Quasi-Rigid Rotor Harmonic
Oscillator approximation. Modified from a script by Steven Wheeler.
See: Grimme, S. Chem. Eur. J. 2012, 18, 9955.
"""

from __future__ import annotations

from math import isclose
from typing import TYPE_CHECKING

import numpy as np
import scipy.constants as const

from pymatgen.core.units import kb as kb_ev
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing_extensions import Self

    from pymatgen.core import Molecule
    from pymatgen.io.gaussian import GaussianOutput
    from pymatgen.io.qchem.outputs import QCOutput

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Alex Epstein"
__email__ = "aepstein@lbl.gov"
__date__ = "August 1, 2023"
__credits__ = "Ryan Kingsbury, Steven Wheeler, Trevor Seguin, Evan Spotte-Smith"

# Define useful constants
kb = kb_ev * const.eV  # Pymatgen kb [J/K]
light_speed = const.speed_of_light * 100  # [cm/s]
h_plank = const.h  # Planck's constant [J.s]
ideal_gas_const = const.R / const.calorie  # Ideal gas constant [cal/mol/K]
R_ha = const.R / const.value("Hartree energy") / const.Avogadro  # Ideal gas

# Define useful conversion factors
amu_to_kg = const.value("atomic mass unit-kilogram relationship")  # AMU to kg
# kcal2hartree = 0.0015936  # kcal/mol to hartree/mol
kcal_to_hartree = 1000 * const.calorie / const.value("Hartree energy") / const.Avogadro


def get_avg_mom_inertia(mol):
    """
    Calculate the average moment of inertia of a molecule.

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
        for dim in range(3):
            inertia_tensor[dim, dim] += wt * (c[(dim + 1) % 3] ** 2 + c[(dim + 2) % 3] ** 2)
        for ii, jj in [(0, 1), (1, 2), (0, 2)]:
            inertia_tensor[ii, jj] += -wt * c[ii] * c[jj]
            inertia_tensor[jj, ii] += -wt * c[jj] * c[ii]

    # atomic mass unit * angstrom^2 to kg m^2
    inertia_eigen_vals = np.linalg.eig(inertia_tensor)[0] * amu_to_kg * 1e-20

    iav = np.mean(inertia_eigen_vals)

    return iav, inertia_eigen_vals


class QuasiRRHO:
    """Calculate thermochemistry using Grimme's Quasi-RRHO approximation.
    All outputs are in atomic units, e.g. energy outputs are in Hartrees.
    Citation: Grimme, S. Chemistry - A European Journal 18, 9955-9964 (2012).

    Attributes:
        temp (float): Temperature [K]
        press (float): Pressure [Pa]
        v0 (float): Cutoff frequency for Quasi-RRHO method [1/cm]
        entropy_quasiRRHO (float): Quasi-RRHO entropy [Ha/K]
        entropy_ho (float): Total entropy calculated with a harmonic
            oscillator approximation for the vibrational entropy [Ha/K]
        h_corrected (float): Thermal correction to the enthalpy [Ha]
        free_energy_quasiRRHO (float): Quasi-RRHO free energy [Ha]
        free_energy_ho (float): Free energy calculated without the Quasi-RRHO
            method, i.e. with a harmonic oscillator approximation for the
            vibrational entropy [Ha]
    """

    def __init__(
        self,
        mol: Molecule,
        frequencies: list[float],
        energy: float,
        mult: int,
        sigma_r: float = 1,
        temp: float = 298.15,
        press: float = 101_317,
        v0: float = 100,
    ) -> None:
        """
        Args:
            mol (Molecule): Pymatgen molecule
            frequencies (list): List of frequencies (float) [cm^-1]
            energy (float): Electronic energy [Ha]
            mult (int): Spin multiplicity
            sigma_r (int): Rotational symmetry number. Defaults to 1.
            temp (float): Temperature [K]. Defaults to 298.15.
            press (float): Pressure [Pa]. Defaults to 101_317.
            v0 (float): Cutoff frequency for Quasi-RRHO method [cm^-1]. Defaults to 100.
        """
        # TO-DO: calculate sigma_r with PointGroupAnalyzer
        # and/or edit Gaussian and QChem io to parse for sigma_r

        self.temp = temp
        self.press = press
        self.v0 = v0

        self.entropy_quasiRRHO = None  # Ha/K
        self.free_energy_quasiRRHO = None  # Ha
        self.h_corrected = None  # Ha

        self.entropy_ho = None  # Ha/K
        self.free_energy_ho = None  # Ha

        self._get_quasirrho_thermo(
            mol=mol,
            mult=mult,
            frequencies=frequencies,
            elec_energy=energy,
            sigma_r=sigma_r,
        )

    @classmethod
    def from_gaussian_output(cls, output: GaussianOutput, **kwargs) -> Self:
        """
        Args:
            output (GaussianOutput): Pymatgen GaussianOutput object.

        Returns:
            QuasiRRHO: QuasiRRHO class instantiated from a Gaussian Output
        """
        mult = output.spin_multiplicity
        elec_e = output.final_energy
        mol = output.final_structure
        vib_freqs = [freq["frequency"] for freq in output.frequencies[-1]]
        return cls(mol=mol, frequencies=vib_freqs, energy=elec_e, mult=mult, **kwargs)

    @classmethod
    def from_qc_output(cls, output: QCOutput, **kwargs) -> Self:
        """
        Args:
            output (QCOutput): Pymatgen QCOutput object.

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
    def _get_quasirrho_thermo(
        self,
        mol: Molecule,
        mult: int,
        sigma_r: int,
        frequencies: list[float],
        elec_energy: float,
    ) -> None:
        """
        Calculate Quasi-RRHO thermochemistry.

        Args:
            mol (Molecule): Pymatgen molecule
            mult (int): Spin multiplicity
            sigma_r (int): Rotational symmetry number
            frequencies (list): List of frequencies [cm^-1]
            elec_energy (float): Electronic energy [Ha]
        """
        # Calculate mass in kg
        mass: float = 0
        for site in mol:
            mass += site.specie.atomic_mass
        mass *= amu_to_kg

        # Calculate vibrational temperatures
        vib_temps = [freq * light_speed * h_plank / kb for freq in frequencies if freq > 0]

        # Translational component of entropy and energy
        qt = (2 * np.pi * mass * kb * self.temp / (h_plank * h_plank)) ** (3 / 2) * kb * self.temp / self.press
        st = ideal_gas_const * (np.log(qt) + 5 / 2)
        et = 3 * ideal_gas_const * self.temp / 2

        # Electronic component of Entropy
        se = ideal_gas_const * np.log(mult)

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
            qr = 8 * np.pi**2 * i * kb * self.temp / (sigma_r * (h_plank * h_plank))
            sr = ideal_gas_const * (np.log(qr) + 1)
            er = ideal_gas_const * self.temp
        else:
            rot_temps = [h_plank**2 / (np.pi**2 * kb * 8 * i) for i in i_eigen]
            qr = np.sqrt(np.pi) / sigma_r * self.temp ** (3 / 2) / np.sqrt(rot_temps[0] * rot_temps[1] * rot_temps[2])
            sr = ideal_gas_const * (np.log(qr) + 3 / 2)
            er = 3 * ideal_gas_const * self.temp / 2

        # Vibrational component of Entropy and Energy
        ev = sv_quasiRRHO = sv = 0

        for vt in vib_temps:
            ev += vt * (1 / 2 + 1 / (np.exp(vt / self.temp) - 1))
            sv_temp = vt / (self.temp * (np.exp(vt / self.temp) - 1)) - np.log(1 - np.exp(-vt / self.temp))
            sv += sv_temp

            mu = h_plank / (8 * np.pi**2 * vt * light_speed)
            mu_prime = mu * Bav / (mu + Bav)
            s_rotor = 1 / 2 + np.log(np.sqrt(8 * np.pi**3 * mu_prime * kb * self.temp / h_plank**2))
            weight = 1 / (1 + (self.v0 / vt) ** 4)
            sv_quasiRRHO += weight * sv_temp + (1 - weight) * s_rotor

        sv_quasiRRHO *= ideal_gas_const
        sv *= ideal_gas_const
        ev *= ideal_gas_const
        e_tot = (et + er + ev) * kcal_to_hartree / 1000
        self.h_corrected = e_tot + ideal_gas_const * self.temp * kcal_to_hartree / 1000
        self.entropy_ho = st + sr + sv + se
        self.free_energy_ho = elec_energy + self.h_corrected - (self.temp * self.entropy_ho * kcal_to_hartree / 1000)
        self.entropy_quasiRRHO = st + sr + sv_quasiRRHO + se
        self.free_energy_quasiRRHO = (
            elec_energy + self.h_corrected - (self.temp * self.entropy_quasiRRHO * kcal_to_hartree / 1000)
        )
