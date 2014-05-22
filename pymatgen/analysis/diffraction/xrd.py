#!/usr/bin/env python

"""
This module implements an XRD pattern calculator.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "5/22/14"


import math, cmath
import os

import numpy as np
import json


#XRD wavelengths in angstroms
WAVELENGTHS = {
    "CuKa": 0.1542 * 10
}

with open(os.path.join(os.path.dirname(__file__),
                       "atomic_scattering_factors.json")) as f:
    ATOMIC_SCATTERING_FACTORS = json.load(f)


class XRDCalculator(object):

    def __init__(self, radiation="CuKa"):
        self.radiation = radiation
        self.wavelength = WAVELENGTHS[radiation]

    def get_xrd_data(self, structure):
        """
        Calculates the XRD data for a structure. This takes into account the
        atomic scattering factors and the Lorentz polarization factor and
        multiplicity factors, but not the Debye-Waller (temperature) factor.

        Args:
            structure: Input structure

        Returns:
            {XRD data} in the form of [
                [two_theta, [scaled_intensity, [h, k, l]]
            ]
        """
        wavelength = self.wavelength
        latt = structure.lattice
        # Obtain reciprocal lattice and points within limiting sphere
        recip_latt = latt.reciprocal_lattice_crystallographic
        recip_pts = recip_latt.get_points_in_sphere([[0, 0, 0]], [0, 0, 0],
                                                    2 / wavelength)

        intensities = {}
        for hkl, g_hkl, ind in sorted(recip_pts, key=lambda d: (d[1], -d[0][0],
                                                            -d[0][1], -d[0][2])):
            if g_hkl != 0:
                theta = math.asin(wavelength / (2 / g_hkl))
                s = g_hkl / 2
                s_2 = s ** 2
                F_hkl = 0
                for site in structure:
                    el = site.specie
                    asf = ATOMIC_SCATTERING_FACTORS[el.symbol]
                    fs = el.Z - 41.78214 * s_2 * sum([d[0] * math.exp(-d[1] *
                                                                      s_2)
                                                      for d in asf])
                    F_hkl += fs * cmath.exp(2j * math.pi
                                            * np.dot(hkl, site.frac_coords))

                I_hkl = (F_hkl * F_hkl.conjugate()).real

                lorentz_factor = (1 + math.cos(2 * theta) ** 2) / (math.sin(theta) ** 2 * math.cos(theta))
                twotheta = 2 * theta / math.pi * 180
                if twotheta in intensities:
                    intensities[twotheta][0] += I_hkl * lorentz_factor
                else:
                    intensities[twotheta] = [I_hkl * lorentz_factor,
                                             tuple(hkl)]

        cleaned_intensities = {}
        prev = 0
        for k in sorted(intensities.keys()):
            if abs(k - prev) < 1e-8:
                cleaned_intensities[prev][0] += intensities[k][0]
            else:
                cleaned_intensities[k] = intensities[k]
                prev = k
        max_intensity = max([v[0] for v in cleaned_intensities.values()])
        data = []
        for k in sorted(cleaned_intensities.keys()):
            v = cleaned_intensities[k]
            data.append([k, [v[0] / max_intensity * 100, v[1]]])
        return data

    def show_xrd_plot(self, structure, two_theta_range=None,
                      annotate_peaks=True):
        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(16, 10)
        two_theta_range = [-1, float("inf")] if two_theta_range is None else two_theta_range
        for d in self.get_xrd_data(structure):
            if two_theta_range[0] <= d[0] <= two_theta_range[1]:
                plt.plot([d[0], d[0]], [0, d[1][0]], color='k',
                         linewidth=3, label=str(d[1][1]))
                if annotate_peaks:
                    plt.annotate(str(d[1][1]), xy=[d[0], d[1][0]],
                                 xytext=[d[0], d[1][0]], fontsize=16)
        plt.show()
