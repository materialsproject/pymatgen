# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines classes to represent all xas and stitching methods
"""
import math
import sys
import warnings
from typing import List

import numpy as np
from scipy.interpolate import interp1d

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.spectrum import Spectrum
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if sys.version_info >= (3, 8):
    from typing import Literal
else:
    from typing_extensions import Literal

__author__ = "Chen Zheng, Yiming Chen"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "3.0"
__maintainer__ = "Yiming Chen"
__email__ = "chz022@ucsd.edu, yic111@ucsd.edu"
__date__ = "July 17, 2020"


class XAS(Spectrum):
    """
    Basic XAS object.

    Args:
        x: A sequence of x-ray energies in eV
        y: A sequence of mu(E)
        structure (Structure): Structure associated with the spectrum
        absorbing_element (Element): Element associated with the spectrum
        edge (str): Absorption edge associated with the spectrum
        spectrum_type (str): 'XANES' or 'EXAFS'
        absorbing_index (None or int): If None, the spectrum is assumed to be a
         site-weighted spectrum, which is comparable to experimental one.
         Otherwise, it indicates that the absorbing_index for a site-wise spectrum.


    .. attribute: x
        The sequence of energies

    .. attribute: y
        The sequence of mu(E)

    .. attribute: absorbing_element
        The absorbing_element of the spectrum

    .. attribute: edge
        The edge of the spectrum

    .. attribute: spectrum_type
        XANES or EXAFS spectrum

    .. attribute: absorbing_index
        The absorbing_index of the spectrum



    """

    XLABEL = "Energy"
    YLABEL = "Intensity"

    def __init__(
        self,
        x,
        y,
        structure,
        absorbing_element,
        edge="K",
        spectrum_type="XANES",
        absorbing_index=None,
    ):
        """
        Initializes a spectrum object.
        """

        super().__init__(x, y, structure, absorbing_element, edge)
        self.structure = structure
        self.absorbing_element = absorbing_element
        self.edge = edge
        self.spectrum_type = spectrum_type
        self.e0 = self.x[np.argmax(np.gradient(self.y) / np.gradient(self.x))]
        # Wavenumber, k is calculated based on equation
        #             k^2=2*m/(hbar)^2*(E-E0)
        self.k = [np.sqrt((i - self.e0) / 3.8537) if i > self.e0 else -np.sqrt((self.e0 - i) / 3.8537) for i in self.x]
        self.absorbing_index = absorbing_index
        # check for empty spectra and negative intensities
        if sum(1 for i in self.y if i <= 0) / len(self.y) > 0.05:
            raise ValueError("Please double check the intensities. Most of them are non-positive values. ")

    def __str__(self):
        return "%s %s Edge %s for %s: %s" % (
            self.absorbing_element,
            self.edge,
            self.spectrum_type,
            self.structure.composition.reduced_formula,
            super().__str__(),
        )

    def stitch(self, other: "XAS", num_samples: int = 500, mode: Literal["XAFS", "L23"] = "XAFS") -> "XAS":
        """
        Stitch XAS objects to get the full XAFS spectrum or L23 edge XANES
        spectrum depending on the mode.

        1. Use XAFS mode for stitching XANES and EXAFS with same absorption edge.
            The stitching will be performed based on wavenumber, k.
            for k <= 3, XAS(k) = XAS[XANES(k)]
            for 3 < k < max(xanes_k), will interpolate according to
                XAS(k)=f(k)*mu[XANES(k)]+(1-f(k))*mu[EXAFS(k)]
                where f(k)=cos^2((pi/2) (k-3)/(max(xanes_k)-3)
            for k > max(xanes_k), XAS(k) = XAS[EXAFS(k)]
        2. Use L23 mode for stitching L2 and L3 edge XANES for elements with
            atomic number <=30.

        Args:
            other: Another XAS object.
            num_samples(int): Number of samples for interpolation.
            mode("XAFS" | "L23"): Either XAFS mode for stitching XANES and EXAFS
                or L23 mode for stitching L2 and L3.

        Returns:
            XAS object: The stitched spectrum.
        """
        m = StructureMatcher()
        if not m.fit(self.structure, other.structure):
            raise ValueError("The input structures for spectra mismatch")
        if not self.absorbing_element == other.absorbing_element:
            raise ValueError("The absorbing elements for spectra are different")
        if not self.absorbing_index == other.absorbing_index:
            raise ValueError("The absorbing indexes for spectra are different")

        if mode == "XAFS":
            if not self.edge == other.edge:
                raise ValueError("Only spectrum with the same absorption edge can be stitched in XAFS mode.")
            if self.spectrum_type == other.spectrum_type:
                raise ValueError("Need one XANES and one EXAFS spectrum to stitch in XAFS mode")

            xanes = self if self.spectrum_type == "XANES" else other
            exafs = self if self.spectrum_type == "EXAFS" else other
            if max(xanes.x) < min(exafs.x):
                raise ValueError("Energy overlap between XANES and EXAFS is needed for stitching")

            # for k <= 3
            wavenumber, mu = [], []  # type: List[float],  List[float]
            idx = xanes.k.index(min(self.k, key=lambda x: (abs(x - 3), x)))
            mu.extend(xanes.y[:idx])
            wavenumber.extend(xanes.k[:idx])

            # for 3 < k < max(xanes.k)
            fs = []  # type: List[float]
            ks = np.linspace(3, max(xanes.k), 50)
            for k in ks:
                f = np.cos((math.pi / 2) * (k - 3) / (max(xanes.k) - 3)) ** 2
                fs.append(f)
            f_xanes = interp1d(
                np.asarray(xanes.k),
                np.asarray(xanes.y),
                bounds_error=False,
                fill_value=0,
            )
            f_exafs = interp1d(
                np.asarray(exafs.k),
                np.asarray(exafs.y),
                bounds_error=False,
                fill_value=0,
            )
            mu_xanes = f_xanes(ks)
            mu_exafs = f_exafs(ks)
            mus = [fs[i] * mu_xanes[i] + (1 - fs[i]) * mu_exafs[i] for i in np.arange(len(ks))]
            mu.extend(mus)
            wavenumber.extend(ks)

            # for k > max(xanes.k)
            idx = exafs.k.index(min(exafs.k, key=lambda x: (abs(x - max(xanes.k)))))
            mu.extend(exafs.y[idx:])
            wavenumber.extend(exafs.k[idx:])

            # interpolation
            f_final = interp1d(np.asarray(wavenumber), np.asarray(mu), bounds_error=False, fill_value=0)
            wavenumber_final = np.linspace(min(wavenumber), max(wavenumber), num=num_samples)
            mu_final = f_final(wavenumber_final)
            energy_final = [
                3.8537 * i ** 2 + xanes.e0 if i > 0 else -3.8537 * i ** 2 + xanes.e0 for i in wavenumber_final
            ]

            return XAS(
                energy_final,
                mu_final,
                self.structure,
                self.absorbing_element,
                xanes.edge,
                "XAFS",
            )

        if mode == "L23":
            if self.spectrum_type != "XANES" or other.spectrum_type != "XANES":
                raise ValueError("Only XANES spectrum can be stitched in L23 mode.")
            if self.edge not in ["L2", "L3"] or other.edge not in ["L2", "L3"] or self.edge == other.edge:
                raise ValueError("Need one L2 and one L3 edge spectrum to stitch in L23 mode.")
            l2_xanes = self if self.edge == "L2" else other
            l3_xanes = self if self.edge == "L3" else other
            if l2_xanes.absorbing_element.number > 30:
                raise ValueError("Does not support L2,3-edge XANES for {} element".format(l2_xanes.absorbing_element))

            l2_f = interp1d(
                l2_xanes.x,
                l2_xanes.y,
                bounds_error=False,
                fill_value="extrapolate",
                kind="cubic",
            )
            l3_f = interp1d(l3_xanes.x, l3_xanes.y, bounds_error=True, fill_value=0, kind="cubic")
            energy = list(np.linspace(min(l3_xanes.x), max(l3_xanes.x), num=num_samples))
            mu = [i + j for i, j in zip([0 if i < 0 else i for i in l2_f(energy)], l3_f(energy))]
            # check for jumps at the onset of L2-edge XANES
            idx = energy.index(min(energy, key=lambda x: (abs(x - l2_xanes.x[0]))))
            if abs(mu[idx] - mu[idx - 1]) / (mu[idx - 1]) > 0.1:
                warnings.warn(
                    "There might exist a jump at the L2 and L3-edge junction.",
                    UserWarning,
                )

            return XAS(energy, mu, self.structure, self.absorbing_element, "L23", "XANES")

        raise ValueError("Invalid mode. Only XAFS and L23 are supported.")


def site_weighted_spectrum(xas_list: List["XAS"], num_samples: int = 500) -> "XAS":
    """
    Obtain site-weighted XAS object based on site multiplicity for each
    absorbing index and its corresponding site-wise spectrum.

    Args:
        xas_list([XAS]): List of XAS object to be weighted
        num_samples(int): Number of samples for interpolation

    Returns:
        XAS object: The site-weighted spectrum
    """
    m = StructureMatcher()
    groups = m.group_structures([i.structure for i in xas_list])
    if len(groups) > 1:
        raise ValueError("The input structures mismatch")
    if not len({i.absorbing_element for i in xas_list}) == len({i.edge for i in xas_list}) == 1:
        raise ValueError(
            "Can only perform site-weighting for spectra with same absorbing element and same absorbing edge."
        )
    if len({i.absorbing_index for i in xas_list}) == 1 or None in {i.absorbing_index for i in xas_list}:
        raise ValueError("Need at least two site-wise spectra to perform site-weighting")

    sa = SpacegroupAnalyzer(groups[0][0])
    ss = sa.get_symmetrized_structure()
    maxes, mines = [], []
    fs = []
    multiplicities = []

    for xas in xas_list:
        multiplicity = len(ss.find_equivalent_sites(ss[xas.absorbing_index]))
        multiplicities.append(multiplicity)
        maxes.append(max(xas.x))
        mines.append(min(xas.x))
        # use 3rd-order spline interpolation for mu (idx 3) vs energy (idx 0).
        f = interp1d(
            np.asarray(xas.x),
            np.asarray(xas.y),
            bounds_error=False,
            fill_value=0,
            kind="cubic",
        )
        fs.append(f)
    # Interpolation within the intersection of x-axis ranges.
    x_axis = np.linspace(max(mines), min(maxes), num=num_samples)
    weighted_spectrum = np.zeros(num_samples)
    sum_multiplicities = sum(multiplicities)

    for i, j in enumerate(multiplicities):
        weighted_spectrum += (j * fs[i](x_axis)) / sum_multiplicities

    return XAS(
        x_axis,
        weighted_spectrum,
        ss,
        xas.absorbing_element,
        xas.edge,
        xas.spectrum_type,
    )
