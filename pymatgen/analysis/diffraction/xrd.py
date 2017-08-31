# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

from math import sin, cos, asin, pi, degrees, radians
import os
import collections

import numpy as np
import json

from pymatgen.core.spectrum import Spectrum
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.plotting import add_fig_kwargs

"""
This module implements an XRD pattern calculator.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "5/22/14"


# XRD wavelengths in angstroms
WAVELENGTHS = {
    "CuKa": 1.54184,
    "CuKa2": 1.54439,
    "CuKa1": 1.54056,
    "CuKb1": 1.39222,
    "MoKa": 0.71073,
    "MoKa2": 0.71359,
    "MoKa1": 0.70930,
    "MoKb1": 0.63229,
    "CrKa": 2.29100,
    "CrKa2": 2.29361,
    "CrKa1": 2.28970,
    "CrKb1": 2.08487,
    "FeKa": 1.93735,
    "FeKa2": 1.93998,
    "FeKa1": 1.93604,
    "FeKb1": 1.75661,
    "CoKa": 1.79026,
    "CoKa2": 1.79285,
    "CoKa1": 1.78896,
    "CoKb1": 1.63079,
    "AgKa": 0.560885,
    "AgKa2": 0.563813,
    "AgKa1": 0.559421,
    "AgKb1": 0.497082,
}

with open(os.path.join(os.path.dirname(__file__),
                       "atomic_scattering_params.json")) as f:
    ATOMIC_SCATTERING_PARAMS = json.load(f)


class XRDPattern(Spectrum):
    """
    A representation of an XRDPattern
    """

    XLABEL = "$2\\Theta$"
    YLABEL = "Intensity"

    def __init__(self, x, y, hkls, d_hkls):
        """
        Args:
            x: Two theta angles.
            y: Intensities
            hkls: [{(h, k, l): mult}] {(h, k, l): mult} is a dict of Miller
                indices for all diffracted lattice facets contributing to each
                intensity.
            d_hkls: List of interplanar spacings.
        """
        super(XRDPattern, self).__init__(x, y, hkls,d_hkls)
        self.hkls = hkls
        self.d_hkls = d_hkls


class XRDCalculator(object):
    """
    Computes the XRD pattern of a crystal structure.

    This code is implemented by Shyue Ping Ong as part of UCSD's NANO106 -
    Crystallography of Materials. The formalism for this code is based on
    that given in Chapters 11 and 12 of Structure of Materials by Marc De
    Graef and Michael E. McHenry. This takes into account the atomic
    scattering factors and the Lorentz polarization factor, but not
    the Debye-Waller (temperature) factor (for which data is typically not
    available). Note that the multiplicity correction is not needed since
    this code simply goes through all reciprocal points within the limiting
    sphere, which includes all symmetrically equivalent facets. The algorithm
    is as follows

    1. Calculate reciprocal lattice of structure. Find all reciprocal points
       within the limiting sphere given by :math:`\\frac{2}{\\lambda}`.

    2. For each reciprocal point :math:`\\mathbf{g_{hkl}}` corresponding to
       lattice plane :math:`(hkl)`, compute the Bragg condition
       :math:`\\sin(\\theta) = \\frac{\\lambda}{2d_{hkl}}`

    3. Compute the structure factor as the sum of the atomic scattering
       factors. The atomic scattering factors are given by

       .. math::

           f(s) = Z - 41.78214 \\times s^2 \\times \\sum\\limits_{i=1}^n a_i \
           \\exp(-b_is^2)

       where :math:`s = \\frac{\\sin(\\theta)}{\\lambda}` and :math:`a_i`
       and :math:`b_i` are the fitted parameters for each element. The
       structure factor is then given by

       .. math::

           F_{hkl} = \\sum\\limits_{j=1}^N f_j \\exp(2\\pi i \\mathbf{g_{hkl}}
           \\cdot \\mathbf{r})

    4. The intensity is then given by the modulus square of the structure
       factor.

       .. math::

           I_{hkl} = F_{hkl}F_{hkl}^*

    5. Finally, the Lorentz polarization correction factor is applied. This
       factor is given by:

       .. math::

           P(\\theta) = \\frac{1 + \\cos^2(2\\theta)}
           {\\sin^2(\\theta)\\cos(\\theta)}
    """

    # Tuple of available radiation keywords.
    AVAILABLE_RADIATION = tuple(WAVELENGTHS.keys())

    # Tolerance in which to treat two peaks as having the same two theta.
    TWO_THETA_TOL = 1e-5

    # Tolerance in which to treat a peak as effectively 0 if the scaled
    # intensity is less than this number. Since the max intensity is 100,
    # this means the peak must be less than 1e-5 of the peak intensity to be
    # considered as zero. This deals with numerical issues where systematic
    # absences do not cancel exactly to zero.
    SCALED_INTENSITY_TOL = 1e-3

    def __init__(self, wavelength="CuKa", symprec=0, debye_waller_factors=None):
        """
        Initializes the XRD calculator with a given radiation.

        Args:
            wavelength (str/float): The wavelength can be specified as either a
                float or a string. If it is a string, it must be one of the
                supported definitions in the AVAILABLE_RADIATION class
                variable, which provides useful commonly used wavelengths.
                If it is a float, it is interpreted as a wavelength in
                angstroms. Defaults to "CuKa", i.e, Cu K_alpha radiation.
            symprec (float): Symmetry precision for structure refinement. If
                set to 0, no refinement is done. Otherwise, refinement is
                performed using spglib with provided precision.
            debye_waller_factors ({element symbol: float}): Allows the
                specification of Debye-Waller factors. Note that these
                factors are temperature dependent.
        """
        if isinstance(wavelength, float):
            self.wavelength = wavelength
        else:
            self.radiation = wavelength
            self.wavelength = WAVELENGTHS[wavelength]
        self.symprec = symprec
        self.debye_waller_factors = debye_waller_factors or {}

    def get_xrd_pattern(self, structure, scaled=True, two_theta_range=(0, 90)):
        """
        Calculates the XRD pattern for a structure.

        Args:
            structure (Structure): Input structure
            scaled (bool): Whether to return scaled intensities. The maximum
                peak is set to a value of 100. Defaults to True. Use False if
                you need the absolute values to combine XRD plots.
            two_theta_range ([float of length 2]): Tuple for range of
                two_thetas to calculate in degrees. Defaults to (0, 90). Set to
                None if you want all diffracted beams within the limiting
                sphere of radius 2 / wavelength.

        Returns:
            (XRDPattern)
        """
        if self.symprec:
            finder = SpacegroupAnalyzer(structure, symprec=self.symprec)
            structure = finder.get_refined_structure()

        wavelength = self.wavelength
        latt = structure.lattice
        is_hex = latt.is_hexagonal()

        # Obtained from Bragg condition. Note that reciprocal lattice
        # vector length is 1 / d_hkl.
        min_r, max_r = (0, 2 / wavelength) if two_theta_range is None else \
            [2 * sin(radians(t / 2)) / wavelength for t in two_theta_range]

        # Obtain crystallographic reciprocal lattice points within range
        recip_latt = latt.reciprocal_lattice_crystallographic
        recip_pts = recip_latt.get_points_in_sphere(
            [[0, 0, 0]], [0, 0, 0], max_r)
        if min_r:
            recip_pts = [pt for pt in recip_pts if pt[1] >= min_r]

        # Create a flattened array of zs, coeffs, fcoords and occus. This is
        # used to perform vectorized computation of atomic scattering factors
        # later. Note that these are not necessarily the same size as the
        # structure as each partially occupied specie occupies its own
        # position in the flattened array.
        zs = []
        coeffs = []
        fcoords = []
        occus = []
        dwfactors = []

        for site in structure:
            for sp, occu in site.species_and_occu.items():
                zs.append(sp.Z)
                try:
                    c = ATOMIC_SCATTERING_PARAMS[sp.symbol]
                except KeyError:
                    raise ValueError("Unable to calculate XRD pattern as "
                                     "there is no scattering coefficients for"
                                     " %s." % sp.symbol)
                coeffs.append(c)
                dwfactors.append(self.debye_waller_factors.get(sp.symbol, 0))
                fcoords.append(site.frac_coords)
                occus.append(occu)

        zs = np.array(zs)
        coeffs = np.array(coeffs)
        fcoords = np.array(fcoords)
        occus = np.array(occus)
        dwfactors = np.array(dwfactors)
        peaks = {}
        two_thetas = []

        for hkl, g_hkl, ind in sorted(
                recip_pts, key=lambda i: (i[1], -i[0][0], -i[0][1], -i[0][2])):
            # Force miller indices to be integers.
            hkl = [int(round(i)) for i in hkl]
            if g_hkl != 0:

                d_hkl = 1 / g_hkl

                # Bragg condition
                theta = asin(wavelength * g_hkl / 2)

                # s = sin(theta) / wavelength = 1 / 2d = |ghkl| / 2 (d =
                # 1/|ghkl|)
                s = g_hkl / 2

                # Store s^2 since we are using it a few times.
                s2 = s ** 2

                # Vectorized computation of g.r for all fractional coords and
                # hkl.
                g_dot_r = np.dot(fcoords, np.transpose([hkl])).T[0]

                # Highly vectorized computation of atomic scattering factors.
                # Equivalent non-vectorized code is::
                #
                #   for site in structure:
                #      el = site.specie
                #      coeff = ATOMIC_SCATTERING_PARAMS[el.symbol]
                #      fs = el.Z - 41.78214 * s2 * sum(
                #          [d[0] * exp(-d[1] * s2) for d in coeff])
                fs = zs - 41.78214 * s2 * np.sum(
                    coeffs[:, :, 0] * np.exp(-coeffs[:, :, 1] * s2), axis=1)

                dw_correction = np.exp(-dwfactors * s2)

                # Structure factor = sum of atomic scattering factors (with
                # position factor exp(2j * pi * g.r and occupancies).
                # Vectorized computation.
                f_hkl = np.sum(fs * occus * np.exp(2j * pi * g_dot_r)
                               * dw_correction)

                # Lorentz polarization correction for hkl
                lorentz_factor = (1 + cos(2 * theta) ** 2) / \
                    (sin(theta) ** 2 * cos(theta))

                # Intensity for hkl is modulus square of structure factor.
                i_hkl = (f_hkl * f_hkl.conjugate()).real

                two_theta = degrees(2 * theta)

                if is_hex:
                    # Use Miller-Bravais indices for hexagonal lattices.
                    hkl = (hkl[0], hkl[1], - hkl[0] - hkl[1], hkl[2])
                # Deal with floating point precision issues.
                ind = np.where(np.abs(np.subtract(two_thetas, two_theta)) <
                               XRDCalculator.TWO_THETA_TOL)
                if len(ind[0]) > 0:
                    peaks[two_thetas[ind[0][0]]][0] += i_hkl * lorentz_factor
                    peaks[two_thetas[ind[0][0]]][1].append(tuple(hkl))
                else:
                    peaks[two_theta] = [i_hkl * lorentz_factor, [tuple(hkl)],
                                        d_hkl]
                    two_thetas.append(two_theta)

        # Scale intensities so that the max intensity is 100.
        max_intensity = max([v[0] for v in peaks.values()])
        x = []
        y = []
        hkls = []
        d_hkls = []
        for k in sorted(peaks.keys()):
            v = peaks[k]
            fam = get_unique_families(v[1])
            if v[0] / max_intensity * 100 > XRDCalculator.SCALED_INTENSITY_TOL:
                x.append(k)
                y.append(v[0])
                hkls.append(fam)
                d_hkls.append(v[2])
        xrd = XRDPattern(x, y, hkls, d_hkls)
        if scaled:
            xrd.normalize(mode="max", value=100)
        return xrd

    def get_xrd_plot(self, structure, two_theta_range=(0, 90),
                     annotate_peaks=True, ax=None, with_labels=True,
                     fontsize=16):
        """
        Returns the XRD plot as a matplotlib.pyplot.

        Args:
            structure: Input structure
            two_theta_range ([float of length 2]): Tuple for range of
                two_thetas to calculate in degrees. Defaults to (0, 90). Set to
                None if you want all diffracted beams within the limiting
                sphere of radius 2 / wavelength.
            annotate_peaks: Whether to annotate the peaks with plane
                information.
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            with_labels: True to add xlabels and ylabels to the plot.
            fontsize: (int) fontsize for peak labels.

        Returns:
            (matplotlib.pyplot)
        """
        if ax is None:
            from pymatgen.util.plotting import pretty_plot
            plt = pretty_plot(16, 10)
            ax = plt.gca()
        else:
            # This to maintain the type of the return value.
            import matplotlib.pyplot as plt

        xrd = self.get_xrd_pattern(structure, two_theta_range=two_theta_range)

        for two_theta, i, hkls, d_hkl in zip(xrd.x, xrd.y, xrd.hkls, xrd.d_hkls):
            if two_theta_range[0] <= two_theta <= two_theta_range[1]:
                label = ", ".join([str(hkl) for hkl in hkls.keys()])
                ax.plot([two_theta, two_theta], [0, i], color='k',
                         linewidth=3, label=label)
                if annotate_peaks:
                    ax.annotate(label, xy=[two_theta, i],
                                xytext=[two_theta, i], fontsize=fontsize)

        if with_labels:
            ax.set_xlabel(r"$2\theta$ ($^\circ$)")
            ax.set_ylabel("Intensities (scaled)")

        if hasattr(ax, "tight_layout"):
            ax.tight_layout()

        return plt

    def show_xrd_plot(self, structure, two_theta_range=(0, 90),
                      annotate_peaks=True):
        """
        Shows the XRD plot.

        Args:
            structure (Structure): Input structure
            two_theta_range ([float of length 2]): Tuple for range of
                two_thetas to calculate in degrees. Defaults to (0, 90). Set to
                None if you want all diffracted beams within the limiting
                sphere of radius 2 / wavelength.
            annotate_peaks (bool): Whether to annotate the peaks with plane
                information.
        """
        self.get_xrd_plot(structure, two_theta_range=two_theta_range,
                          annotate_peaks=annotate_peaks).show()

    @add_fig_kwargs
    def plot_structures(self, structures, two_theta_range=(0, 90),
                       annotate_peaks=True, fontsize=6, **kwargs):
        """
        Plot XRD for multiple structures on the same figure.

        Args:
            structures (Structure): List of structures
            two_theta_range ([float of length 2]): Tuple for range of
                two_thetas to calculate in degrees. Defaults to (0, 90). Set to
                None if you want all diffracted beams within the limiting
                sphere of radius 2 / wavelength.
            annotate_peaks (bool): Whether to annotate the peaks with plane
                information.
            fontsize: (int) fontsize for peak labels.
        """
        import matplotlib.pyplot as plt
        nrows = len(structures)
        fig, axes = plt.subplots(nrows=nrows, ncols=1, sharex=True, squeeze=False)

        for i, (ax, structure) in enumerate(zip(axes.ravel(), structures)):
            self.get_xrd_plot(structure, two_theta_range=two_theta_range,
                              annotate_peaks=annotate_peaks,
                              fontsize=fontsize, ax=ax, with_labels=i == nrows - 1)
            spg_symbol, spg_number = structure.get_space_group_info()
            ax.set_title("{} {} ({}) ".format(structure.formula, spg_symbol, spg_number))

        return fig


def get_unique_families(hkls):
    """
    Returns unique families of Miller indices. Families must be permutations
    of each other.

    Args:
        hkls ([h, k, l]): List of Miller indices.

    Returns:
        {hkl: multiplicity}: A dict with unique hkl and multiplicity.
    """
    # TODO: Definitely can be sped up.
    def is_perm(hkl1, hkl2):
        h1 = np.abs(hkl1)
        h2 = np.abs(hkl2)
        return all([i == j for i, j in zip(sorted(h1), sorted(h2))])

    unique = collections.defaultdict(list)
    for hkl1 in hkls:
        found = False
        for hkl2 in unique.keys():
            if is_perm(hkl1, hkl2):
                found = True
                unique[hkl2].append(hkl1)
                break
        if not found:
            unique[hkl1].append(hkl1)

    pretty_unique = {}
    for k, v in unique.items():
        pretty_unique[sorted(v)[-1]] = len(v)

    return pretty_unique
