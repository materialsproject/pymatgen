# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import os
import collections
import abc

import numpy as np
import json

from pymatgen.core.spectrum import Spectrum
from pymatgen.util.plotting import add_fig_kwargs

"""
This module implements an diffraction pattern calculator.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "5/22/14"


with open(os.path.join(os.path.dirname(__file__),
                       "atomic_scattering_params.json")) as f:
    ATOMIC_SCATTERING_PARAMS = json.load(f)


class DiffractionPattern(Spectrum):
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
        super(DiffractionPattern, self).__init__(x, y, hkls,d_hkls)
        self.hkls = hkls
        self.d_hkls = d_hkls


class DiffractionPatternCalculator(abc.ABC):
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
    # Tolerance in which to treat two peaks as having the same two theta.
    TWO_THETA_TOL = 1e-5

    # Tolerance in which to treat a peak as effectively 0 if the scaled
    # intensity is less than this number. Since the max intensity is 100,
    # this means the peak must be less than 1e-5 of the peak intensity to be
    # considered as zero. This deals with numerical issues where systematic
    # absences do not cancel exactly to zero.
    SCALED_INTENSITY_TOL = 1e-3

    @abc.abstractmethod
    def get_pattern(self, structure, scaled=True, two_theta_range=(0, 90)):
        """
        Calculates the diffraction pattern for a structure.

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
            (DiffractionPattern)
        """
        pass

    def get_plot(self, structure, two_theta_range=(0, 90),
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

        xrd = self.get_pattern(structure, two_theta_range=two_theta_range)

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

    def show_plot(self, structure, two_theta_range=(0, 90),
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
        self.get_plot(structure, two_theta_range=two_theta_range,
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
            self.get_plot(structure, two_theta_range=two_theta_range,
                          annotate_peaks=annotate_peaks,
                          fontsize=fontsize, ax=ax, with_labels=i == nrows - 1)
            spg_symbol, spg_number = structure.get_space_group_info()
            ax.set_title("{} {} ({}) ".format(structure.formula, spg_symbol,
                                              spg_number))

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
