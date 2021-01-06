# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
A module to calculate free energies using the Quasi-Rigid Rotor Harmonic
 Oscillator approximation.
See: Grimme, S. Chem. Eur. J. 2012, 18, 9955
"""

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Alex Epstein"
__email__ = "aepstein@lbl.gov"
__date__ = "January 5, 2021"
__credits__ = "Trevor Seguin"

from typing import Union

import numpy as np
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.qchem.outputs import QCOutput

# Define useful constants
kb = 1.380662E-23  # J/K
c = 29979245800  # cm/s
h = 6.62606957E-34  # J.s
R = 1.987204  # kcal/mol

# Define useful conversion factors
cm2kcal = 0.0028591459
cm2hartree = 4.5563352812122E-06
kcal2hartree = 0.0015936
amu2kg = 1.66053886E-27
bohr2angs = 0.52917721092


class QuasiRRHO:
    """
    Class to calculate thermochemistry using Grimme's Quasi-RRHO approximation
    """

    def __init__(self, output: Union[GaussianOutput, QCOutput, dict],
                 sigma_r=1, temp=298.15, press=101317, conc=1,
                 v0=100):
        """

        :param output: Requires input of a Gaussian output file,
                        QChem output file, or dictionary of necessary inputs:
                        {"mol": Molecule, "mult": spin multiplicity (int),
                        "frequencies": list of vibrational frequencies [a.u.],
                        elec_energy": electronic energy [a.u.]}
        :param sigma_r: Rotational symmetry number
        :param temp: Temperature [K]
        :param press: Pressure [Pa]
        :param conc: Solvent concentration [M]
        :param v0: Cutoff frequency for Quasi-RRHO method [cm^1]
        """

        # TO-DO: calculate sigma_r with PointGroupAnalyzer
        # and/or edit Gaussian and QChem io to parse for sigma_r

        self.sigma_r = sigma_r
        self.temp = temp
        self.press = press
        self.conc = conc
        self.v0 = v0

        if isinstance(output, GaussianOutput):
            mult = output.spin_multiplicity
            sigma_r = self.sigma_r
            elec_e = output.final_energy
            mol = output.final_structure
            vib_freqs = [f['frequency'] for f in output.frequencies[-1]]

            self._get_quasirrho_thermo(mol=mol, mult=mult, sigma_r=sigma_r,
                                       frequencies=vib_freqs,
                                       elec_energy=elec_e)

        if isinstance(output, QCOutput):
            mult = output.data['multiplicity']
            elec_e = output.data["SCF_energy_in_the_final_basis_set"]

            if output.data["optimization"]:
                mol = output.data["molecule_from_last_geometry"]
            else:
                mol = output.data["initial_molecule"]
            frequencies = output.data['frequencies']

            self._get_quasirrho_thermo(mol=mol, mult=mult,
                                       sigma_r=self.sigma_r,
                                       frequencies=frequencies,
                                       elec_energy=elec_e)

        if isinstance(output, dict):
            self._get_quasirrho_thermo(mol=output["mol"],
                                       mult=output["mult"],
                                       sigma_r=self.sigma_r,
                                       frequencies=output["frequencies"],
                                       elec_energy=output["elec_energy"])

    def get_avg_mom_inertia(self, mol):
        """
        Caclulate the average moment of inertia of a molecule
        :param mol: Molecule
        :return: average moment of inertia, eigenvalues of inertia tensor
        """
        centered_mol = mol.get_centered_molecule()
        inertia_tensor = np.zeros((3, 3))
        for site in centered_mol:
            c = site.coords
            wt = site.specie.atomic_mass
            for i in range(3):
                inertia_tensor[i, i] += wt * (
                        c[(i + 1) % 3] ** 2 + c[(i + 2) % 3] ** 2
                )
            for i, j in [(0, 1), (1, 2), (0, 2)]:
                inertia_tensor[i, j] += -wt * c[i] * c[j]
                inertia_tensor[j, i] += -wt * c[j] * c[i]

        inertia_eigenvals = np.multiply(
            np.linalg.eig(inertia_tensor)[0],
            amu2kg * 1E-20).tolist()  # amuangs^2 to kg m^2

        iav = np.average(inertia_eigenvals)

        return iav, inertia_eigenvals

    def _get_quasirrho_thermo(self, mol, mult, sigma_r, frequencies,
                              elec_energy):
        """
        Caclulate Quasi-RRHO thermochemistry
        :param mol: Molecule
        :param mult: Spin multiplicity
        :param sigma_r: Rotational symmetry number
        :param frequencies: List of frequencies [a.u.]
        :param elec_energy: Electornic energy [a.u.]
        """

        # Calculate mass in kg
        mass = 0
        for site in mol.sites:
            mass += site.specie.atomic_mass
        mass *= amu2kg

        # Calculate vibrational temperatures
        vib_temps = []
        for f in frequencies:
            if f > 0:
                vib_temps.append(f * c * h / kb)

        # Get properties related to rotational symmetry. Bav is average moment of inertia
        Bav, i_eigen = self.get_avg_mom_inertia(mol)
        rot_temps = [h ** 2 / (np.pi ** 2 * kb * 8 * i) for i in i_eigen]

        # Translational component of entropy and energy
        qt = (2 * np.pi * mass * kb * self.temp / (h * h)) ** (
                3 / 2) * kb * self.temp / self.press
        st = R * (np.log(qt) + 5 / 2)
        et = 3 * R * self.temp / 2

        # Electronic component of Entropy
        se = R * np.log(mult)

        # Rotational component of Entropy and Energy
        qr = np.sqrt(np.pi) / sigma_r * self.temp ** (3 / 2) / np.sqrt(
            rot_temps[0] * rot_temps[1] * rot_temps[2])
        sr = R * (np.log(qr) + 3 / 2)
        er = 3 * R * self.temp / 2

        # Vibrational component of Entropy and Energy
        ev = 0
        sv_quasiRRHO = 0

        for vt in vib_temps:
            ev += vt * (1 / 2 + 1 / (np.exp(vt / self.temp) - 1))
            sv_temp = vt / (self.temp * (np.exp(vt / self.temp) - 1)) - np.log(
                1 - np.exp(-vt / self.temp))
            mu = h / (8 * np.pi ** 2 * vt * c)
            mu_prime = mu * Bav / (mu + Bav)
            srotor = 1 / 2 + np.log(
                np.sqrt(8 * np.pi ** 3 * mu_prime * kb * self.temp / h ** 2))
            weight = 1 / (1 + (self.v0 / vt) ** 4)
            sv_quasiRRHO += weight * sv_temp + (1 - weight) * srotor

        sv_quasiRRHO *= R
        ev *= R
        etot = (et + er + ev) * kcal2hartree / 1000
        h_corrected = etot + R * self.temp * kcal2hartree / 1000

        molarity_corr = 0.000003166488448771253 * self.temp * np.log(
            0.082057338 * self.temp * self.conc)
        self.entropy_quasiRRHO = st + sr + sv_quasiRRHO + se
        self.free_energy_quasiRRHO = (elec_energy + h_corrected -
                                      (self.temp * self.entropy_quasiRRHO
                                       * kcal2hartree / 1000))

        self.concentration_corrected_g_quasiRRHO = (
                self.free_energy_quasiRRHO + molarity_corr)
