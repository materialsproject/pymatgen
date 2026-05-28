"""This module implements an XRD pattern calculator."""

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

with open(
    os.path.join(os.path.dirname(__file__), "atomic_scattering_params.json"),
    "rb",
) as file:
    ATOMIC_SCATTERING_PARAMS = orjson.loads(file.read())


class XRDCalculator(AbstractDiffractionPatternCalculator):
    r"""
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
       within the limiting sphere given by \frac{2}{\lambda}.

    2. For each reciprocal point \mathbf{g_{hkl}} corresponding to
       lattice plane (hkl), compute the Bragg condition
       \sin(\theta) = \frac{ \lambda}{2d_{hkl}}

    3. Compute the structure factor as the sum of the atomic scattering
       factors. The atomic scattering factors are given by

           f(s) = Z - 41.78214 \times s^2 \times \sum \limits_{i=1}^n a_i \exp(-b_is^2)

       where s = \ frac{\ sin(\ theta)}{\ lambda} and a_i
       and b_i are the fitted parameters for each element. The
       structure factor is then given by

           F_{hkl} = \sum \limits_{j=1}^N f_j  \exp(2 \pi i  \mathbf{g_{hkl}} \cdot  \mathbf{r})

    4. The intensity is then given by the modulus square of the structure factor.

           I_{hkl} = F_{hkl}F_{hkl}^*

    5. Finally, the Lorentz polarization correction factor is applied. This
       factor is given by:

           P(\theta) = \frac{1 + \cos^2(2 \theta)}{\sin^2(\theta) \cos(\theta)}
    """

    # Tuple of available radiation keywords.
    AVAILABLE_RADIATION = tuple(WAVELENGTHS)

    def __init__(self, wavelength="CuKa", symprec: float = 0, debye_waller_factors=None):
        """Initialize the XRD calculator with a given radiation.

        Args:
            wavelength (str | float): The wavelength can be specified as either a
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
        if isinstance(wavelength, float | int):
            self.wavelength = wavelength
        elif isinstance(wavelength, str):
            self.radiation = wavelength
            self.wavelength = WAVELENGTHS[wavelength]
        else:
            wavelength_type = type(wavelength).__name__
            raise TypeError(f"{wavelength_type=} must be either float, int or str")
        self.symprec = symprec
        self.debye_waller_factors = debye_waller_factors or {}

    def get_pattern(self, structure: Structure, scaled=True, two_theta_range=(0, 90)):
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
            DiffractionPattern: XRD pattern
        """
        if self.symprec:
            finder = SpacegroupAnalyzer(structure, symprec=self.symprec)
            structure = finder.get_refined_structure()

        wavelength = self.wavelength
        latt = structure.lattice
        is_hex = latt.is_hexagonal()

        min_r, max_r = (
            (0, 2 / wavelength)
            if two_theta_range is None
            else [2 * math.sin(math.radians(t / 2)) / wavelength for t in two_theta_range]
        )

        recip_latt = latt.reciprocal_lattice_crystallographic
        recip_pts = recip_latt.get_points_in_sphere([[0, 0, 0]], [0, 0, 0], max_r)
        if min_r:
            recip_pts = [pt for pt in recip_pts if pt[1] >= min_r]

        # --- Build per-site arrays ---
        _zs, _coeffs, _fcoords, _occus, _dwfactors = [], [], [], [], []
        for site in structure:
            for sp, occu in site.species.items():
                _zs.append(sp.Z)
                try:
                    c = ATOMIC_SCATTERING_PARAMS[sp.symbol]
                except KeyError:
                    raise ValueError(
                        f"Unable to calculate XRD pattern as there is no scattering coefficients for {sp.symbol}."
                    )
                _coeffs.append(c)
                _dwfactors.append(self.debye_waller_factors.get(sp.symbol, 0))
                _fcoords.append(site.frac_coords)
                _occus.append(occu)

        zs = np.array(_zs)  # (N,)
        coeffs = np.array(_coeffs)  # (N, 4, 2)
        fcoords = np.array(_fcoords)  # (N, 3)
        occus = np.array(_occus)  # (N,)
        dwfactors = np.array(_dwfactors)  # (N,)

        # --- Unpack reciprocal points & filter g_hkl == 0 ---
        recip_pts_sorted = sorted(recip_pts, key=lambda i: (i[1], -i[0][0], -i[0][1], -i[0][2]))

        hkls_raw = np.array([pt[0] for pt in recip_pts_sorted])  # (M, 3)
        g_hkls = np.array([pt[1] for pt in recip_pts_sorted])  # (M,)

        nonzero = g_hkls != 0
        hkls_raw = hkls_raw[nonzero]
        g_hkls = g_hkls[nonzero]

        hkls_int = np.round(hkls_raw).astype(int)  # (M, 3)

        # --- Fully vectorized computation over all M hkl points ---
        # shapes: (M,)
        theta = np.arcsin(np.clip(wavelength * g_hkls / 2, -1, 1))
        s2 = (g_hkls / 2) ** 2  # (M,)

        # Atomic scattering factors: (M, N)
        #   fs[m, n] = zs[n] - 41.78214 * s2[m] * sum_k(coeffs[n,k,0] * exp(-coeffs[n,k,1]*s2[m]))
        # coeffs: (N, 4, 2)  →  broadcast s2: (M, 1, 1)
        s2_mnk = s2[:, None, None]  # (M, 1, 1)
        gauss = np.sum(
            coeffs[None, :, :, 0] * np.exp(-coeffs[None, :, :, 1] * s2_mnk),
            axis=2,
        )  # (M, N)
        fs = zs[None, :] - 41.78214 * s2[:, None] * gauss  # (M, N)

        # Debye-Waller per atom, per hkl: (M, N)
        dw = np.exp(-dwfactors[None, :] * s2[:, None])

        # g·r for all hkl and all atoms: (M, N)
        g_dot_r = hkls_int.astype(float) @ fcoords.T  # (M, N)

        # Structure factors: (M,)
        f_hkl = np.sum(
            fs * occus[None, :] * np.exp(2j * math.pi * g_dot_r) * dw,
            axis=1,
        )
        i_hkl = (f_hkl * f_hkl.conjugate()).real  # (M,)

        # Lorentz-polarization factor: (M,)
        cos2t = np.cos(2 * theta)
        sint = np.sin(theta)
        cost = np.cos(theta)
        lorentz = (1 + cos2t**2) / (sint**2 * cost)

        intensities = i_hkl * lorentz  # (M,)
        two_thetas_arr = np.degrees(2 * theta)  # (M,)

        # --- Merge peaks within TWO_THETA_TOL using rounding-based binning ---
        tol = AbstractDiffractionPatternCalculator.TWO_THETA_TOL
        bin_keys = np.round(two_thetas_arr / tol).astype(int)

        peaks: dict[int, list] = {}
        for m in range(len(g_hkls)):
            hkl = tuple(int(v) for v in hkls_int[m])  # np.int64 -> int
            if is_hex:
                hkl = (hkl[0], hkl[1], -hkl[0] - hkl[1], hkl[2])
            key = bin_keys[m]
            d_hkl = 1.0 / g_hkls[m]
            if key in peaks:
                peaks[key][0] += float(intensities[m])  # np.float64 -> float
                peaks[key][1].append(hkl)
            else:
                peaks[key] = [float(intensities[m]), [hkl], float(two_thetas_arr[m]), float(d_hkl)]

        # --- Build output ---
        max_intensity = max(v[0] for v in peaks.values())
        tol_scaled = AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL
        x, y, hkls_out, d_hkls_out = [], [], [], []

        for key in sorted(peaks):
            v = peaks[key]
            fam = get_unique_families(v[1])
            if v[0] / max_intensity * 100 > tol_scaled:
                x.append(float(v[2]))
                y.append(float(v[0]))
                hkls_out.append([{"hkl": hkl, "multiplicity": mult} for hkl, mult in fam.items()])
                d_hkls_out.append(float(v[3]))

        xrd = DiffractionPattern(x, y, hkls_out, d_hkls_out)
        if scaled:
            xrd.normalize(mode="max", value=100)
        return xrd
