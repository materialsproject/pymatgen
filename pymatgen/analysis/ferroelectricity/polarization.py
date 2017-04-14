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

"""
This module contains classes for recovering the spontaneous
polarization from multiple calculations along a nonpolar to polar
ferroelectric distortion.

We recommend using our calc_ionic function for calculating the ionic
polarization rather than the values from OUTCAR.

We find that the ionic dipole moment reported in OUTCAR differs from
the naive calculation of \sum_i Z_i r_i where i is the index of the
atom, r is the distance in Angstroms along the lattice vectors.
Compare to VASP dipol.F. SUBROUTINE POINT_CHARGE_DIPOL.

We are able to recover a smooth same branch polarization more frequently
using the naive calculation in calc_ionic than using the ionic dipole
moment reported in the OUTCAR.
"""

# Load ZVAL dictionaries
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(MODULE_DIR, "ZVAL.yaml"),'r') as f:
    ZVAL = yaml.load(f)
VASP_PBE_ZVAL = ZVAL['VASP_PBE_ZVAL']


def zval_dict_from_potcar(potcar):
    """
    Creates zval_dictionary for calculating the ionic polarization from
    Potcar Object

    potcar: Potcar object
    """
    zval_dict = {}
    for p in potcar:
        zval_dict.update({p.element: p.ZVAL})

def create_zval_dict_from_pseudo_dict(species_potcar_dict, pseudo_dict = None):
    """
    Return a dictionary of pseudopotential ZVALs given a dictionary of
    which pseudopotentials

    species_potcar_dict: dictionary for species and pseudopotential name.
        Example {‘Li’ : ‘Li' , ’Nb’: ‘Nb_sv', ‘O’: ‘O’}
    pseudo_dict: default is VASP_PBE_ZVAL
        Default matches pseudopotentials used for MPRelaxSet and MPStaticSet.

    """
    if pseudo_dict == None:
        pseudo_dict = VASP_PBE_ZVAL
    zval_dict = {}
    for key,value in species_potcar_dict.iteritems():
        zval_dict.update({key: pseudo_dict[value]})
    return zval_dict

def calc_ionic(site, structure, zval):
    """
    Calculate the ionic dipole moment using ZVAL from pseudopotential

    site: PeriodicSite
    structure: Structure
    zval: Charge value for ion (ZVAL for VASP pseudopotential)

    Returns polarization in electron Angstroms.
    """
    norms = structure.lattice.lengths_and_angles[0]
    return np.multiply(norms,-site.frac_coords * zval)

def get_total_ionic_dipole(structure, zval_dict):
    """
    Get the total ionic dipole moment for a structure.

    structure: pymatgen Structure
    zval_dict: specie, zval dictionary pairs
    center (np.array with shape [3,1]) : dipole center used by VASP
    tiny (float) : tolerance for determining boundary of calculation.
    """

    tot_ionic = []
    for site in structure:
        zval = zval_dict[str(site.specie)]
        tot_ionic.append(calc_ionic(site, structure, zval))
    return np.sum(tot_ionic, axis=0)

class Polarization(object):
    """
    Class for recovering the same branch polarization for a set of
    polarization calculations along the nonpolar - polar distortion
    path of a ferroelectric.

    p_elecs, p_ions, and structures lists should be given in order
    of nonpolar to polar!

    It is assumed that the electronic and ionic dipole moment values
    are given in electron Angstroms along the three lattice directions
    (a,b,c) rather than (x,y,z).

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
        Returns the dipole / polarization quanta along a, b, and c for
        all structures.
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
        return rms


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
