# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
A module to calculate free energeis using a Quasi-Rigid Rotor Harmonic Oscillator approximation
See: Grimme, S. Chem. Eur. J. 2012, 18, 9955
"""

__author__ = "Alex Epstein"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Alex Epstein"
__email__ = "aepstein@lbl.gov"
__date__ = ""
__credits__ = "Trevor Seguin, Steven Wheeler"

import numpy as np

from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.qchem.outputs import QCOutput

from typing import Union

# Some useful constants
kb = 1.380662E-23  # J/K
c = 29979245800  # cm/s
h = 6.62606957E-34  # J.s
R = 1.987204  # kcal/mol
cm2kcal = 0.0028591459
cm2hartree = 4.5563352812122E-06
kcal2hartree = 0.0015936
amu2kg = 1.66053886E-27


class QuasiRRHO():
    """
    Grimme's Quasi-RRHO approximation
    """

    def __init__(self, output: Union[GaussianOutput, QCOutput, dict], temp=298.15, press=101317, conc=1, v0=100):
        self.temp = temp
        self.press = press
        self.conc = conc
        self.v0 = v0

        if isinstance(output, GaussianOutput):
            self._get_thermochemistry_from_GaussianOutput(output)

        if isinstance(output, QCOutput):
            self._get_thermochemistry_from_QCOutput(output)

        if isinstance(output, dict):
            self._get_thermochemistry_from_dict(output)

    def _get_thermochemistry_from_GaussianOutput(self, GaussianOutput):
        mult = GaussianOutput.spin_multiplicity
        mass = GaussianOutput.mass * amu2kg
        sigma_r = GaussianOutput.rot_sym_num
        rot_temps = GaussianOutput.rot_temps
        h_corrected = GaussianOutput.corrections['Enthalpy']
        final_energy = GaussianOutput.final_energy

        # Get frequencies for last frequency run
        vib_freqs = []
        vib_temps = []
        for freq_dict in GaussianOutput.frequencies[-1]:
            f = freq_dict['frequency']
        if f > 0:
            vib_freqs.append(f)
            vib_temps.append(f * c * h / kb)

        Bav = (h ** 2 / (24 * np.pi ** 2 * kb)) * (
                1 / rot_temps[0] + 1 / rot_temps[1] + 1 / rot_temps[2])
        print('bav', Bav)

        # Translational component of entropy
        qt = (2 * np.pi * mass * kb * self.temp / (h * h)) ** (3 / 2) * kb * self.temp / self.press
        st = R * (np.log(qt) + 5 / 2)

        # Electronic component of Entropy
        se = R * np.log(mult)

        # Rotational component of Entropy
        qr = np.sqrt(np.pi) / sigma_r * self.temp ** (3 / 2) / np.sqrt(
            rot_temps[0] * rot_temps[1] * rot_temps[2])
        sr = R * (np.log(qr) + 3 / 2)

        # Vibrational component of Entropy and Energy
        sv = 0
        sv_quasiRRHO = 0

        for vt in vib_temps:
            sv_temp = vt / (self.temp * (np.exp(vt / self.temp) - 1)) - np.log(1 - np.exp(-vt / self.temp))
            sv += sv_temp

            # quasi-RRHO
            mu = h / (8 * np.pi ** 2 * vt * c)
            mu_prime = mu * Bav / (mu + Bav)
            srotor = 1 / 2 + np.log(np.sqrt(8 * np.pi ** 3 * mu_prime * kb * self.temp / h ** 2))
            weight = 1 / (1 + (self.v0 / vt) ** 4)
            sv_quasiRRHO += weight * sv_temp + (1 - weight) * srotor

        sv *= R
        sv_quasiRRHO *= R

        molarity_corr = 0.000003166488448771253 * self.temp * np.log(0.082057338 * self.temp * self.conc)

        self.entropy_quasiRRHO = st + sr + sv_quasiRRHO + se
        self.free_energy_quasiRRHO = final_energy + h_corrected - self.temp * self.entropy_quasiRRHO * kcal2hartree / 1000
        self.concentration_corrected_g_quasiRRHO = self.free_energy_quasiRRHO + molarity_corr

    def _get_thermochemistry_from_QCOutput(self, QCOutput):
        mult = QCOutput.data['multiplicity']
        mol = QCOutput.data['initial_molecule']

        mass = 0
        for site in mol.sites:
            mass += site.specie._atomic_mass
        mass *= amu2kg

        # TO-DO fix QCOutput to get rotational symmetry number
        sigma_r = 1
        rot_temps = GaussianOutput.rot_temps
        h_corrected = GaussianOutput.corrections['Enthalpy']
        final_energy = GaussianOutput.final_energy

        # Get frequencies for last frequency run
        vib_freqs = []
        vib_temps = []
        for freq_dict in GaussianOutput.frequencies[-1]:
            f = freq_dict['frequency']
        if f > 0:
            vib_freqs.append(f)
            vib_temps.append(f * c * h / kb)

        Bav = (h ** 2 / (24 * np.pi ** 2 * kb)) * (
                1 / rot_temps[0] + 1 / rot_temps[1] + 1 / rot_temps[2])

        # Translational component of entropy
        qt = (2 * np.pi * mass * kb * self.temp / (h * h)) ** (3 / 2) * kb * self.temp / self.press
        st = R * (np.log(qt) + 5 / 2)

        # Electronic component of Entropy
        se = R * np.log(mult)

        # Rotational component of Entropy
        qr = np.sqrt(np.pi) / sigma_r * self.temp ** (3 / 2) / np.sqrt(
            rot_temps[0] * rot_temps[1] * rot_temps[2])
        sr = R * (np.log(qr) + 3 / 2)

        # Vibrational component of Entropy and Energy
        sv = 0
        sv_quasiRRHO = 0

        for vt in vib_temps:
            sv_temp = vt / (self.temp * (np.exp(vt / self.temp) - 1)) - np.log(1 - np.exp(-vt / self.temp))
            sv += sv_temp

            # quasi-RRHO
            mu = h / (8 * np.pi ** 2 * vt * c)
            mu_prime = mu * Bav / (mu + Bav)
            srotor = 1 / 2 + np.log(np.sqrt(8 * np.pi ** 3 * mu_prime * kb * self.temp / h ** 2))
            weight = 1 / (1 + (self.v0 / vt) ** 4)
            sv_quasiRRHO += weight * sv_temp + (1 - weight) * srotor

        sv *= R
        sv_quasiRRHO *= R

        molarity_corr = 0.000003166488448771253 * self.temp * np.log(0.082057338 * self.temp * self.conc)

        self.entropy_quasiRRHO = st + sr + sv_quasiRRHO + se
        self.free_energy_quasiRRHO = final_energy + h_corrected - self.temp * self.entropy_quasiRRHO * kcal2hartree / 1000
        self.concentration_corrected_g_quasiRRHO = self.free_energy_quasiRRHO + molarity_corr

    @classmethod
    def _get_thermochemistry_from_params(self, param_dict):
        self.data = {'final_thermochemistry': {}, 'thermal_corrections': {}}

        zpve = param_dict['zpve']
        elec_zvpe = param_dict['elec_zvpe']
        rot_temps = param_dict['rot_temps']
        vib_temps = param_dict['vib_temps']
        vib_freqs = param_dict['vib_freqs']
        mult = param_dict['mult']
        mass = param_dict['mass']
        sigma_r = param_dict['sigma_r']

        Bav = (h ** 2 / (24 * np.pi ** 2 * kb)) * (
                1 / rot_temps[0] + 1 / rot_temps[1] + 1 / rot_temps[2])

        # Translational entropy
        qt = (2 * np.pi * self.mass * kb * self.temp / (h * h)) ** (3 / 2) * kb * self.temp / self.press
        st = R * (np.log(qt) + 5 / 2)

        # Translation component of Energy
        et = 3 * R * self.temp / 2;

        # Electronic component of Entropy
        se = R * np.log(self.mult)

        # Rotational component of Entropy
        qr = np.sqrt(np.pi) / self.sigma_r * self.temp ** (3 / 2) / np.sqrt(
            rot_temps[0] * rot_temps[1] * rot_temps[2])
        sr = R * (np.log(qr) + 3 / 2)

        # Rotational component of Energy
        er = 3 * R * self.temp / 2

        # Vibrational component of Entropy and Energy
        ev = 0
        sv = 0
        sv_quasiRRHO = 0

        for vt in vib_temps:
            sv_temp = vt / (self.temp * (np.exp(vt / self.temp) - 1)) - np.log(1 - np.exp(-vt / self.temp))
            sv += sv_temp
            ev += vt * (1 / 2 + 1 / (np.exp(vt / self.temp) - 1))

            # quasi-RRHO
            mu = h / (8 * np.pi ** 2 * vt * c)
            mu_prime = mu * Bav / (mu + Bav)
            sr = 1 / 2 + np.log(np.sqrt(8 * np.pi ** 3 * mu_prime * kb * self.temp / h ** 2))
            weight = 1 / (1 + (self.v0 / vt) ** 4)
            sv_quasiRRHO += weight * sv_temp + (1 - weight) * sr

        sv *= R
        ev *= R
        sv_quasiRRHO *= R

        self.data['final_thermochemistry']['electronic_energy'] = (self.elec_zvpe - self.zpve) * kcal2hartree / 1000
        self.data['thermal_corrections']['e_corrected'] = et + er + ev
        self.data['thermal_corrections']['h_corrected'] = (self.data['thermal_e_correction'] + R * self.temp) \
                                                          * kcal2hartree / 1000
        self.data['thermal_corrections']['quasiRRHO_entropy'] = st + sr + sv_quasiRRHO + se
        self.data['thermal_corrections']['entropy'] = st + sr + sv + se
        self.data['thermal_corrections']['harmonic_g_corrected'] = (self.data['thermal_corrections'][
                                                                        'h_corrected'] - self.temp *
                                                                    self.data['thermal_corrections'][
                                                                        'entropy']) * kcal2hartree / 1000
        self.data['thermal_corrections']['quasiRRHO_g_corrected'] = (self.data['thermal_corrections'][
                                                                         'h_corrected'] - self.temp *
                                                                     self.data['thermal_corrections'][
                                                                         'quasiRRHO_entropy']) * kcal2hartree / 1000

        molarity_corr = 0.000003166488448771253 * self.temp * np.log(0.082057338 * self.T * self.conc)
        self.data['final_thermochemistry']['energies'] = self.data['thermal_corrections']['e_corrected'] + \
                                                         self.data['final_thermochemistry']['electronic_energy']
        self.data['final_thermochemistry']['enthalpy'] = self.data['final_thermochemistry']['electronic_energy'] + \
                                                         self.data['thermal_corrections']['h_corrected']

        self.data['final_thermochemistry']['harmonic_g'] = self.data['final_thermochemistry']['electronic_energy'] + \
                                                           self.data['thermal_corrections'][
                                                               'harmonic_g_corrected']

        self.data['final_thermochemistry']['quasiRRHO_g'] = self.data['final_thermochemistry'][
                                                                'electronic_energy'] + \
                                                            self.data['thermal_corrections'][
                                                                'quasiRRHO_g_corrected']
        self.data['final_thermochemistry']['concentration_corrected_quasiRRHO_g'] = \
            self.data['final_thermochemistry']['electronic_energy'] + \
            self.data['thermal_corrections'][
                'quasiRRHO_g_corrected'] + molarity_corr

        self.quasiRRHO_g = self.data['final_thermochemistry']['quasiRRHO_g']
