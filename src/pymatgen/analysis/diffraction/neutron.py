"""This module implements a neutron diffraction (ND) pattern calculator."""

from __future__ import annotations

import math
import os
from typing import TYPE_CHECKING

import numpy as np
import orjson

from pymatgen.analysis.diffraction.core import (
    AbstractDiffractionPatternCalculator,
    DiffractionPattern,
    get_unique_families,
)
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from pymatgen.core import Structure

__author__ = "Yuta Suzuki"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Yuta Suzuki"
__email__ = "resnant@outlook.jp"
__date__ = "4/19/18"

# This table was cited from "Neutron Data Booklet" 2nd ed (Old City 2003).
with open(
    os.path.join(os.path.dirname(__file__), "neutron_scattering_length.json"),
    "rb",
) as file:
    ATOMIC_SCATTERING_LEN = orjson.loads(file.read())


class NDCalculator(AbstractDiffractionPatternCalculator):
    """
    Computes the powder neutron diffraction pattern of a crystal structure.
    This code is a slight modification of XRDCalculator in
    pymatgen.analysis.diffraction.xrd. See it for details of the algorithm.
    Main changes by using neutron instead of X-ray are as follows:

    1. Atomic scattering length is a constant.
    2. Polarization correction term of Lorentz factor is unnecessary.

    Reference:
    Marc De Graef and Michael E. McHenry, Structure of Materials 2nd ed,
    Chapter13, Cambridge University Press 2003.
    """

    def __init__(self, wavelength=1.54184, symprec: float = 0, debye_waller_factors=None):
        """Initialize the ND calculator with a given radiation.

        Args:
            wavelength (float): The wavelength of neutron in angstroms.
                Defaults to 1.54, corresponds to Cu K_alpha x-ray radiation.
            symprec (float): Symmetry precision for structure refinement. If
                set to 0, no refinement is done. Otherwise, refinement is
                performed using spglib with provided precision.
            debye_waller_factors ({element symbol: float}): Allows the
                specification of Debye-Waller factors. Note that these
                factors are temperature dependent.
        """
        self.wavelength = wavelength
        self.symprec = symprec
        self.debye_waller_factors = debye_waller_factors or {}

    def get_pattern(self, structure: Structure, scaled=True, two_theta_range=(0, 90)):
        """
        Calculates the powder neutron diffraction pattern for a structure.

        Args:
            structure (Structure): Input structure
            scaled (bool): Whether to return scaled intensities. The maximum
                peak is set to a value of 100. Defaults to True. Use False if
                you need the absolute values to combine ND plots.
            two_theta_range ([float of length 2]): Tuple for range of
                two_thetas to calculate in degrees. Defaults to (0, 90). Set to
                None if you want all diffracted beams within the limiting
                sphere of radius 2 / wavelength.

        Returns:
            DiffractionPattern: ND pattern
        """
        if self.symprec:
            finder = SpacegroupAnalyzer(structure, symprec=self.symprec)
            structure = finder.get_refined_structure()

        wavelength = self.wavelength
        lattice = structure.lattice
        is_hex = lattice.is_hexagonal()

        # Obtained from Bragg condition. Note that reciprocal lattice
        # vector length is 1 / d_hkl.
        min_r, max_r = (
            (0, 2 / wavelength)
            if two_theta_range is None
            else [2 * math.sin(math.radians(t / 2)) / wavelength for t in two_theta_range]
        )

        # Obtain crystallographic reciprocal lattice points within range
        recip_lattice = lattice.reciprocal_lattice_crystallographic
        recip_pts = recip_lattice.get_points_in_sphere([[0, 0, 0]], [0, 0, 0], max_r)
        if min_r:
            recip_pts = [pt for pt in recip_pts if pt[1] >= min_r]

        # --- Build per-site arrays ---
        _coeffs, _frac_coords, _occus, _dw_factors = [], [], [], []
        for site in structure:
            for sp, occu in site.species.items():
                try:
                    c = ATOMIC_SCATTERING_LEN[sp.symbol]
                except KeyError:
                    raise ValueError(
                        f"Unable to calculate ND pattern as there is no scattering coefficients for {sp.symbol}."
                    )
                _coeffs.append(c)
                _dw_factors.append(self.debye_waller_factors.get(sp.symbol, 0))
                _frac_coords.append(site.frac_coords)
                _occus.append(occu)

        coeffs = np.asarray(_coeffs)  # (N,)
        fcoords = np.asarray(_frac_coords)  # (N,3)
        occus = np.asarray(_occus)  # (N,)
        dwfactors = np.asarray(_dw_factors)  # (N,)

        # --- Unpack reciprocal points ---
        recip_pts_sorted = sorted(
            recip_pts,
            key=lambda i: (i[1], -i[0][0], -i[0][1], -i[0][2]),
        )

        hkls_raw = np.array([pt[0] for pt in recip_pts_sorted])
        g_hkls = np.array([pt[1] for pt in recip_pts_sorted])

        nonzero = g_hkls != 0
        hkls_raw = hkls_raw[nonzero]
        g_hkls = g_hkls[nonzero]

        hkls_int = np.round(hkls_raw).astype(int)  # (M, 3)

        # --- Vectorized diffraction calculation ---
        theta = np.arcsin(np.clip(wavelength * g_hkls / 2, -1, 1))
        s2 = (g_hkls / 2) ** 2

        dw = np.exp(-dwfactors[None, :] * s2[:, None])

        g_dot_r = hkls_int.astype(float) @ fcoords.T

        f_hkl = np.sum(
            coeffs[None, :] * occus[None, :] * np.exp(2j * np.pi * g_dot_r) * dw,
            axis=1,
        )

        i_hkl = (f_hkl * f_hkl.conjugate()).real

        sint = np.sin(theta)
        cost = np.cos(theta)
        lorentz = 1.0 / (sint**2 * cost)

        intensities = i_hkl * lorentz
        two_thetas_arr = np.degrees(2 * theta)

        # --- Merge peaks within tolerance ---
        tol = AbstractDiffractionPatternCalculator.TWO_THETA_TOL
        bin_keys = np.round(two_thetas_arr / tol).astype(int)

        peaks: dict[int, list] = {}

        for m in range(len(g_hkls)):
            hkl = tuple(int(v) for v in hkls_int[m])

            if is_hex:
                hkl = (hkl[0], hkl[1], -hkl[0] - hkl[1], hkl[2])

            key = bin_keys[m]
            d_hkl = 1.0 / g_hkls[m]

            if key in peaks:
                peaks[key][0] += float(intensities[m])
                peaks[key][1].append(hkl)
            else:
                peaks[key] = [
                    float(intensities[m]),
                    [hkl],
                    float(two_thetas_arr[m]),
                    float(d_hkl),
                ]

        max_intensity = max(v[0] for v in peaks.values())

        x, y, hkls_out, d_hkls_out = [], [], [], []
        tol_scaled = AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL

        for key in sorted(peaks):
            intensity, hkls, two_theta, d_hkl = peaks[key]
            fam = get_unique_families(hkls)

            if intensity / max_intensity * 100 > tol_scaled:
                x.append(two_theta)
                y.append(intensity)
                hkls_out.append([{"hkl": hkl, "multiplicity": mult} for hkl, mult in fam.items()])
                d_hkls_out.append(d_hkl)

        nd = DiffractionPattern(x, y, hkls_out, d_hkls_out)

        if scaled:
            nd.normalize(mode="max", value=100)

        return nd
