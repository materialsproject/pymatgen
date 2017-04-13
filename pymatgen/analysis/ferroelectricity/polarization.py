# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals
from __future__ import absolute_import

import os
from math import *
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.cif import CifWriter
from pymatgen.core.lattice import Lattice
import numpy as np
import yaml

"""
This module provides the classes needed to analyze the change in polarization
from a nonpolar reference phase to a polar ferroelectric phase.
"""

__author__ = "Tess Smidt"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "1.0"
__email__ = "tsmidt@berkeley.edu"
__status__ = "Development"
__date__ = "April 8, 2017"


# TO DO
# Figure out calc_ionic
# Give species_potcar_dict example and create YAML for that example
# Allow Polarization to accept arrays or Outcars
# LOTS O' COMMENTS
# unit tests

# Load ZVAL dictionaries
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(MODULE_DIR, "ZVAL.yaml"),'r') as f:
    ZVAL = yaml.load(f)
VASP_PBE_ZVAL = ZVAL['VASP_PBE_ZVAL']

# Load examples species_potcar_dict

def calc_ionic(pos, s, zval, center=None, tiny=0.001):
    """
    Function for calculating the ionic dipole moment for a site.

    pos (pymatgen.core.site.Site) : pymatgen Site
    s (pymatgen.core.structure.Structure) : pymatgen Structure
    zval (float) : number of core electrons of pseudopotential
    center (np.array with shape [3,1]) : dipole center used by VASP
    tiny (float) : tolerance for determining boundary of calculation.
    """
    # lattice vector lenghts
    norms = s.lattice.lengths_and_angles[0]
    # Define center of dipole moment. If not set default is used.
    center = np.array([-100.0,-100.0,-100.0]) if center == None else center
    # Following VASP dipol.F. SUBROUTINE POINT_CHARGE_DIPOL
    temp = (pos.frac_coords - center + 10.5) % 1 - 0.5
    for i in range(3):
        if abs(abs(temp[i]) - 0.5) < tiny/norms[i]:
            temp[i] = 0.0
    # Convert to Angstroms from fractional coords before returning.
    return np.dot(np.transpose(s.lattice.matrix), -temp*zval)

def calc_ionic_2(pos, s, zval, center=None, tiny=0.001):
    """
    Function for calculating the ionic dipole moment for a site.

    pos (pymatgen.core.site.Site) : pymatgen Site
    s (pymatgen.core.structure.Structure) : pymatgen Structure
    zval (float) : number of core electrons of pseudopotential
    center (np.array with shape [3,1]) : dipole center used by VASP
    tiny (float) : tolerance for determining boundary of calculation.
    """
    # lattice vector lenghts
    norms = s.lattice.lengths_and_angles[0]
    # Define center of dipole moment. If not set default is used.
    center = np.array([-100.0,-100.0,-100.0]) if center == None else center
    # Following VASP dipol.F. SUBROUTINE POINT_CHARGE_DIPOL
    temp = (pos.frac_coords - center + 10.5) % 1 - 0.5
    for i in range(3):
        if abs(abs(temp[i]) - 0.5) < tiny/norms[i]:
            temp[i] = 0.0
    # Convert to Angstroms from fractional coords before returning.
    return np.multiply(norms, -temp*zval)

def calc_ionic_naive(pos, s, zval):
    return np.dot(np.transpose(s.lattice.matrix), -pos.frac_coords * zval)

def calc_ionic_naive_2(pos, s, zval):
    norms = s.lattice.lengths_and_angles[0]
    return np.multiply(norms,-pos.frac_coords * zval)

def get_total_ionic_dipole(structure, species_potcar_dict,
                           pseudo_dict = None, center = None, tiny= 0.001):
    """
    Get the total ionic dipole moment for a structure.

    structure: pymatgen Structure
    species_potcar_dict: dictionary for species and pseudopotential name.
        Example {‘Li’ : ‘Li' , ’Nb’: ‘Nb_sv', ‘O’: ‘O’}
    pseudo_dict: default is PBE_ZVAL
    center (np.array with shape [3,1]) : dipole center used by VASP
    tiny (float) : tolerance for determining boundary of calculation.
    """

    pseudo_dict = VASP_PBE_ZVAL if pseudo_dict == None else pseudo_dict
    center = np.array([-100.0, -100.0, -100.0]) if center == None else center

    tot_ionic = []
    for site in structure:
        zval = pseudo_dict[species_potcar_dict[str(site.specie)]]
        tot_ionic.append(calc_ionic(site, structure, zval, center, tiny))
    return np.sum(tot_ionic,axis=0)

def get_total_ionic_dipole_2(structure, species_potcar_dict,
                               pseudo_dict=None, center=None, tiny=0.001):
    """
    Get the total ionic dipole moment for a structure.

    structure: pymatgen Structure
    species_potcar_dict: dictionary for species and pseudopotential name.
        Example {‘Li’ : ‘Li' , ’Nb’: ‘Nb_sv', ‘O’: ‘O’}
    pseudo_dict: default is PBE_ZVAL
    center (np.array with shape [3,1]) : dipole center used by VASP
    tiny (float) : tolerance for determining boundary of calculation.
    """

    pseudo_dict = VASP_PBE_ZVAL if pseudo_dict == None else pseudo_dict
    center = np.array([-100.0, -100.0, -100.0]) if center == None else center

    tot_ionic = []
    for site in structure:
        zval = pseudo_dict[species_potcar_dict[str(site.specie)]]
        tot_ionic.append(calc_ionic_2(site, structure, zval, center, tiny))
    return np.sum(tot_ionic, axis=0)

def get_total_ionic_dipole_naive(structure, species_potcar_dict, pseudo_dict=None):

    """
    Get the total ionic dipole moment for a structure.

    structure: pymatgen Structure
    species_potcar_dict: dictionary for species and pseudopotential name.
        Example {‘Li’ : ‘Li' , ’Nb’: ‘Nb_sv', ‘O’: ‘O’}
    pseudo_dict: default is PBE_ZVAL
    center (np.array with shape [3,1]) : dipole center used by VASP
    tiny (float) : tolerance for determining boundary of calculation.
    """

    pseudo_dict = VASP_PBE_ZVAL if pseudo_dict == None else pseudo_dict

    tot_ionic = []
    for site in structure:
        zval = pseudo_dict[species_potcar_dict[str(site.specie)]]
        tot_ionic.append(calc_ionic_naive(site, structure, zval))
    return np.sum(tot_ionic, axis=0)

def get_total_ionic_dipole_naive_2(structure, species_potcar_dict, pseudo_dict=None):

    """
    Get the total ionic dipole moment for a structure.

    structure: pymatgen Structure
    species_potcar_dict: dictionary for species and pseudopotential name.
        Example {‘Li’ : ‘Li' , ’Nb’: ‘Nb_sv', ‘O’: ‘O’}
    pseudo_dict: default is PBE_ZVAL
    center (np.array with shape [3,1]) : dipole center used by VASP
    tiny (float) : tolerance for determining boundary of calculation.
    """

    pseudo_dict = VASP_PBE_ZVAL if pseudo_dict == None else pseudo_dict

    tot_ionic = []
    for site in structure:
        zval = pseudo_dict[species_potcar_dict[str(site.specie)]]
        tot_ionic.append(calc_ionic_naive_2(site, structure, zval))
    return np.sum(tot_ionic, axis=0)

def get_total_ionic_dipole_consistent(structures, species_potcar_dict, pseudo_dict=None):
    """
    Attempts to use consistent images across distortion

    Parameters
    ----------
    structures
    species_potcar_dict
    pseudo_dict

    Returns
    -------

    """

    pseudo_dict = VASP_PBE_ZVAL if pseudo_dict == None else pseudo_dict

    S = len(structures)
    L = len(structures[0])
    p_ion = np.zeros((S,3))
    for i in range(L):
        orig_sites = [s[i] for s in structures]
        adjust_sites = []
        for j in range(S):
            if j == 0:
                adjust_sites.append(orig_sites[j])
            else:
                new_site, dist = structures[j].get_nearest_site(adjust_sites[-1].coords,
                                                                orig_sites[j])
                adjust_sites.append(new_site)

            site = adjust_sites[-1]
            zval = pseudo_dict[species_potcar_dict[str(site.specie)]]
            p_ion[j,:] += calc_ionic_naive(site, structures[j], zval)
    return p_ion

class Polarization(object):
    """
    Revised object for getting polarization

    list should be given in order of nonpolar to polar
    """
    def __init__(self, p_elecs, p_ions, structures):
        if len(p_elecs) != len(p_ions) or len(p_elecs) != len(structures):
            raise ValueError("The number of electronic polarization and ionic polarization values must be equal.")
        self.p_elecs = np.matrix(p_elecs)
        self.p_ions = np.matrix(p_ions)
        self.structures = structures

    def get_pelecs_and_pions(self, convert_to_muC_per_cm2=False):

        if not convert_to_muC_per_cm2:
            return self.p_elecs, self.p_ions

        if convert_to_muC_per_cm2:
            p_elecs = np.matrix(self.p_elecs).T
            p_ions = np.matrix(self.p_ions).T

            volumes = [s.lattice.volume for s in self.structures]
            e_to_muC = -1.6021766e-13
            cm2_to_A2 = 1e16
            units = 1.0 / np.matrix(volumes)
            units *= e_to_muC * cm2_to_A2

            p_elecs = np.multiply(units, p_elecs)
            p_ions = np.multiply(units, p_ions)

            p_elecs, p_ions = p_elecs.T, p_ions.T

            return p_elecs, p_ions

    def get_same_branch_polarization_data(self, convert_to_muC_per_cm2=False):
        """
        Get same branch polarization for given polarization data.

        convert_to_muC_per_cm2: convert polarization from electron * Angstroms to microCoulomb per centimeter**2
        abc: return polarization in coordinates of a,b,c (versus x,y,z)
        """

        p_elec, p_ion = self.get_pelecs_and_pions()
        p_tot = p_elec + p_ion
        p_tot = np.matrix(p_tot)

        lattices = [s.lattice for s in self.structures]
        volumes = np.matrix([s.lattice.volume for s in self.structures])

        L = len(p_elec)

        # convert polarizations and lattice lengths prior to adjustment
        if convert_to_muC_per_cm2:
            e_to_muC = -1.6021766e-13
            cm2_to_A2 = 1e16
            units = 1.0 / np.matrix(volumes)
            units *= e_to_muC * cm2_to_A2
            # Convert the total polarization
            p_tot = np.multiply(units.T, p_tot)
            # adjust lattices
            for i in range(L):
                lattice = lattices[i]
                l,a = lattice.lengths_and_angles
                lattices[i] = Lattice.from_lengths_and_angles(np.array(l)*units.A1[i],a)

        d_structs = []
        sites = []

        for i in range(L):
            l = lattices[i]
            frac_coord = np.divide(np.matrix(p_tot[i]), np.matrix([l.a, l.b, l.c]))
            d = Structure(l, ["C"], [np.matrix(frac_coord).A1])
            d_structs.append(d)
            site = d[0]
            if i == 0:
                # Adjust nonpolar polarization to be closest to zero.
                # This is compatible with both a polarization of zero or a half quantum.
                prev_site = [0, 0, 0]
                # An alternative method which leaves the nonpolar polarization where it is.
                #sites.append(site)
                #continue
            else:
                prev_site = sites[-1].coords
            new_site = d.get_nearest_site(prev_site, site)
            sites.append(new_site[0])

        adjust_pol = []
        for s, d in zip(sites, d_structs):
            l = d.lattice
            adjust_pol.append(np.multiply(s.frac_coords, np.matrix([l.a, l.b, l.c])).A1)
        adjust_pol = np.matrix(adjust_pol)

        return adjust_pol

    def get_lattice_quanta(self, convert_to_muC_per_cm2 = True):
        """
        Returns the quanta along a, b, and c for all structures.
        """
        lattices = [s.lattice for s in self.structures]
        volumes = np.matrix([s.lattice.volume for s in self.structures])

        L = len(self.structures)

        # convert polarizations and lattice lengths prior to adjustment
        if convert_to_muC_per_cm2:
            e_to_muC = -1.6021766e-13
            cm2_to_A2 = 1e16
            units = 1.0 / np.matrix(volumes)
            units *= e_to_muC * cm2_to_A2
            # adjust lattices
            for i in range(L):
                lattice = lattices[i]
                l,a = lattice.lengths_and_angles
                lattices[i] = Lattice.from_lengths_and_angles(np.array(l)*units.A1[i],a)

        quanta = np.matrix([np.array(l.lengths_and_angles[0]) for l in lattices])

        return quanta

    def get_polarization_change(self):
        tot = self.get_same_branch_polarization_data(convert_to_muC_per_cm2=True)
        return (tot[-1] - tot[0])

    def get_polarization_change_norm(self):
        polar = self.structures[-1]
        a, b, c = polar.lattice.matrix
        a, b, c = a/np.linalg.norm(a), b/np.linalg.norm(b), c/np.linalg.norm(c)
        P = self.get_polarization_change().A1
        P_norm = np.linalg.norm(a * P[0] + b * P[1] + c * P[2])
        return P_norm

    def same_branch_splines(self):
        from scipy.interpolate import UnivariateSpline
        tot = self.get_same_branch_polarization_data(convert_to_muC_per_cm2=True)
        L = tot.shape[0]
        try:
            sp_a = UnivariateSpline(range(L),tot[:,0].A1)
        except:
            sp_a = None
        try:
            sp_b = UnivariateSpline(range(L),tot[:,1].A1)
        except:
            sp_b = None
        try:
            sp_c = UnivariateSpline(range(L),tot[:,2].A1)
        except:
            sp_c = None
        return sp_a, sp_b, sp_c

    def max_spline_jumps(self):
        tot = self.get_same_branch_polarization_data(convert_to_muC_per_cm2=True)
        sps = self.same_branch_splines()
        max_jumps = [None,None,None]
        for i,sp in enumerate(sps):
            if sp != None:
                max_jumps[i] = max(tot[:,i].A1 - sp(range(len(tot[:,i].A1))))
        return max_jumps

    def smoothness(self):
        tot = self.get_same_branch_polarization_data(convert_to_muC_per_cm2=True)
        L = tot.shape[0]
        try:
            sp = self.same_branch_splines()
        except:
            print("Something went wrong.")
            return None
        sp_latt = [sp[i](range(L)) for i in range(3)]
        diff = [sp_latt[i] - tot[:,i].A1 for i in range(3)]
        rms = [np.sqrt(np.sum(np.square(diff[i])) / L) for i in range(3)]
        #rms_mag_norm = [rms[i] / (max(tot[:,i].A1) - min(tot[:,i].A1)) for i in range(3)]
        return rms

    def is_smooth(self, rms_mag_norm_tol = 1e-2):
        """
        Returns whether spline fitted to adjusted a, b, c polarizations are smooth relative to a tolerance.
        """
        tot = self.get_same_branch_polarization_data(convert_to_muC_per_cm2=True)
        L = tot.shape[0]
        try:
            sp = self.same_branch_splines()
        except:
            print("Something went wrong.")
            return None
        sp_latt = [sp[i](range(L)) for i in range(3)]
        diff = [sp_latt[i] - tot[:,i].A1 for i in range(3)]
        rms = [np.sqrt(np.sum(np.square(diff[i])) / L) for i in range(3)]
        rms_mag_norm = [rms[i] / (max(tot[:,i].A1) - min(tot[:,i].A1)) for i in range(3)]
        return [rms_mag_norm[i] <= rms_mag_norm_tol for i in range(3)]


class EnergyTrend(object):
    def __init__(self,energies):
        self.energies = energies

    def spline(self):
        from scipy.interpolate import UnivariateSpline
        sp = UnivariateSpline(range(len(self.energies)),self.energies, k=4)
        return sp

    def smoothness(self):
        energies = self.energies
        try:
            sp = self.spline()
        except:
            print("Energy spline failed.")
            return None
        spline_energies = sp(range(len(energies)))
        diff = spline_energies - energies
        rms = np.sqrt(np.sum(np.square(diff))/len(energies))
        #rms_mag_norm = rms / (max(energies) - min(energies))
        return rms

    def max_spline_jump(self):
        sp = self.spline()
        return max(self.energies - sp(range(len(self.energies))))

    def is_smooth(self, rms_mag_norm_tol = 1e-2):
        energies = self.energies
        try:
            sp = self.spline()
        except:
            print("Energy spline failed.")
            return None
        spline_energies = sp(range(len(energies)))
        diff = spline_energies - energies
        rms = np.sqrt(np.sum(np.square(diff))/len(energies))
        rms_mag_norm = rms / (max(energies) - min(energies))
        return rms_mag_norm <= rms_mag_norm_tol

    def endpoints_minima(self, slope_cutoff = 5e-3):
        energies = self.energies
        try:
            sp = self.spline()
        except:
            print("Energy spline failed.")
            return None
        der = sp.derivative()
        spline_energies = sp(range(len(energies)))
        der_energies = der(range(len(energies)))
        return {"polar" : abs(der_energies[-1]) <= slope_cutoff,
                "nonpolar" : abs(der_energies[0]) <= slope_cutoff}
