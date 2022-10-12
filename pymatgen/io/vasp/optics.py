# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""Classes for parsing and manipulating VASP optical properties calculations."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
from monty.json import MSONable
from scipy import constants, special
from tqdm import tqdm

from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun, Waveder

au2ang = constants.physical_constants["atomic unit of length"][0] / 1e-10
ryd2ev = constants.physical_constants["Rydberg constant times hc in eV"][0]
edeps = 4 * np.pi * 2 * ryd2ev * au2ang  # from constant.inc in VASP

KB = constants.physical_constants["Boltzmann constant in eV/K"][0]


@dataclass
class DielectricFunctionCalculator(MSONable):
    """Class for postprocessing VASP optical properties calculations.

    This objects helps load the different parameters from the vasprun.xml file but allows users to override
    them as needed.

    The standard vasprun.xml from an ``LOPTICS=.True.`` calculation already contains
    the complex frequency dependent dielectric functions.  However you have no way to decompose
    the different contributions.  Since the ``WAVEDER`` file is also written during an optical calculation,
    you can reconstruct the dielectric functions purely in Python and have full control over contribution
    from different bands and k-points.

    VASP's linear optics follow these steps:
        - Calculate the imaginary part
        - Perform symmetry operations (this is not implemented here)
        - Calculate the real part

    Currently, this Calculator only works for ``ISYM=0`` calculations since we cannot gauranttee that our
    externally defined symmetry operations are the same as VASP's.  This can be fixed by printing the
    symmetry operators into the vasprun.xml file.  If this happens in future versions of VASP,
    we can dramatically speed up the calculations here by considering only the irreducible kpoints.
    """

    cder: npt.NDArray
    eigs: npt.NDArray
    kweights: npt.NDArray
    nedos: int
    deltae: float
    ismear: int
    sigma: float
    efermi: float
    cshift: float
    ispin: int
    volume: float

    @classmethod
    def from_vasp_objects(cls, vrun: Vasprun, waveder: Waveder):
        """Construct a DielectricFunction from Vasprun, Kpoint, and Waveder objects.

        Args:
            vrun: Vasprun object
            kpoint: Kpoint object
            waveder: Waveder object
        """
        bands = vrun.eigenvalues
        sspins = [Spin.up, Spin.down]
        eigs = np.stack([bands[spin] for spin in sspins[: vrun.parameters["ISPIN"]]], axis=2)[..., 0]
        eigs = np.swapaxes(eigs, 0, 1)
        cder = waveder.cder
        kweights = vrun.actual_kpoints_weights
        nedos = vrun.parameters["NEDOS"]
        deltae = vrun.dielectric[0][1]
        ismear = vrun.parameters["ISMEAR"]
        sigma = vrun.parameters["SIGMA"]
        cshift = vrun.parameters["CSHIFT"]
        efermi = vrun.efermi
        ispin = vrun.parameters["ISPIN"]
        volume = vrun.final_structure.volume
        if vrun.parameters["ISYM"] != 0:
            raise NotImplementedError("ISYM != 0 is not implemented yet")

        return DielectricFunctionCalculator(
            cder=cder,
            eigs=eigs,
            kweights=kweights,
            nedos=nedos,
            deltae=deltae,
            ismear=ismear,
            sigma=sigma,
            efermi=efermi,
            cshift=cshift,
            ispin=ispin,
            volume=volume,
        )

    @classmethod
    def from_directory(cls, directory: Path | str):
        """Construct a DielectricFunction from a directory containing vasprun.xml and WAVEDER files."""
        d_ = Path(directory)
        vrun = Vasprun(d_ / "vasprun.xml")
        waveder = Waveder.from_binary(d_ / "WAVEDER")
        return cls.from_vasp_objects(vrun, waveder)

    def get_epsilon(
        self,
        idir: int,
        jdir: int,
        efermi: float = None,
        nedos: int = None,
        deltae: float = None,
        ismear: int = None,
        sigma: float = None,
        cshift: float = None,
    ) -> npt.NDArray:
        """Compute the frequency dependent dielectric function.

        Args:
            idir: First direction of the dielectric tensor
            jdir: Second direction of the dielectric tensor
            efermi: Fermi energy
            nedos: Number of points in the DOS
            deltae: Energy step in the DOS
            ismear: Smearing method (only has 0:gaussian, >0:Methfessel-Paxton)
            sigma: Smearing width
            cshift: Complex shift used for Kramer-Kronig transformation
        """

        def _use_default(param, default):
            return param if param is not None else default

        efermi = _use_default(efermi, self.efermi)
        nedos = _use_default(nedos, self.nedos)
        deltae = _use_default(deltae, self.deltae)
        ismear = _use_default(ismear, self.ismear)
        sigma = _use_default(sigma, self.sigma)
        cshift = _use_default(cshift, self.cshift)

        egrid, eps_imag = epsilon_imag(  # type: ignore
            cder=self.cder,
            eigs=self.eigs,
            kweights=self.kweights,
            efermi=efermi,  # type: ignore
            nedos=nedos,  # type: ignore
            deltae=deltae,  # type: ignore
            ismear=ismear,  # type: ignore
            sigma=sigma,  # type: ignore
            idir=idir,
            jdir=jdir,
        )
        # scaling constant: edeps * np.pi / structure.volume
        eps_in = eps_imag * edeps * np.pi / self.volume
        eps = kramers_kronig(eps_in, nedos=nedos, deltae=deltae, cshift=cshift)  # type: ignore
        if idir == jdir:
            eps += 1.0 + 0.0j
        return egrid, eps


def delta_methfessel_paxton(x, n):
    """
    D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
    where H is a Hermite polynomial and
    A_i = (-1)^i / ( i! 4^i sqrt(pi) )
    """
    ii = np.arange(0, n + 1)
    A = (-1) ** ii / (special.factorial(ii) * 4**ii * np.sqrt(np.pi))
    H = special.eval_hermite(ii * 2, np.tile(x, (len(ii), 1)).T)
    return np.exp(-(x * x)) * np.dot(A, H.T)


def step_methfessel_paxton(x, n):
    """
    S_n (x) = (1 + erf x)/2 - exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
    where H is a Hermite polynomial and
    A_i = (-1)^i / ( i! 4^i sqrt(pi) )
    """
    ii = np.arange(1, n + 1)
    A = (-1) ** ii / (special.factorial(ii) * 4**ii * np.sqrt(np.pi))
    H = special.eval_hermite(ii * 2 - 1, np.tile(x, (len(ii), 1)).T)
    return (1.0 + special.erf(x)) / 2.0 - np.exp(-(x * x)) * np.dot(A, H.T)


def delta_func(x, ismear):
    """Replication of VASP's delta function"""
    if ismear < -1:
        raise ValueError("Delta function not implemented for ismear < -1")
    elif ismear == -1:
        return step_func(x, -1) * (1 - step_func(x, -1))
    elif ismear < 0:
        return np.exp(-(x * x)) / np.sqrt(np.pi)
    return delta_methfessel_paxton(x, ismear)


def step_func(x, ismear):
    """Replication of VASP's step function"""
    if ismear < -1:
        raise ValueError("Delta function not implemented for ismear < -1")
    elif ismear == -1:
        return 1 / (1.0 + np.exp(-x))
    elif ismear < 0:
        return 0.5 + 0.5 * special.erf(x)
    return step_methfessel_paxton(x, ismear)


def get_delta(x0: float, sigma: float, nx: int, dx: float, ismear: int = 3):
    """Get the smeared delta function to be added to form the spectrum.

    This replaces the `SLOT` function from VASP. Uses finite differences instead of
    evaluating the delta function since the step function is more likely to have analytic form.

    Args:
        x0: The center of the dielectric function.
        sigma: The width of the smearing
        nx: The number of grid points in the output grid.
        dx: The gridspacing of the output grid.
        ismear: The smearing parameter used by the ``step_func``.

    Return:
        np.array: Array of size `nx` with delta function on the desired outputgrid.

    """
    xgrid = np.arange(0, nx * dx, dx)
    xgrid -= x0
    x_scaled = (xgrid + (dx / 2)) / sigma
    sfun = step_func(x_scaled, ismear)
    dfun = np.zeros_like(xgrid)
    dfun[1:] = (sfun[1:] - sfun[:-1]) / dx
    return dfun


def get_step(x0, sigma, nx, dx, ismear):
    """Get the smeared step function to be added to form the spectrum.

    This replaces the `SLOT` function from VASP.

    Args:
        x0: The center of the dielectric function.
        sigma: The width of the smearing
        nx: The number of grid points in the output grid.
        dx: The gridspacing of the output grid.
        ismear: The smearing parameter used by the ``step_func``.

    Return:
        np.array: Array of size `nx` with step function on the desired outputgrid.
    """
    xgrid = np.arange(0, nx * dx, dx)
    xgrid -= x0
    x_scaled = (xgrid + (dx / 2)) / sigma
    return step_func(x_scaled, ismear)


def epsilon_imag(
    cder: npt.NDArray,
    eigs: npt.NDArray,
    kweights: npt.ArrayLike,
    efermi: float,
    nedos: int,
    deltae: float,
    ismear: int,
    sigma: float,
    idir: int,
    jdir: int,
):
    """Replicate the EPSILON_IMAG function of VASP.

    Args:
        cder: The data written to the WAVEDER (nbands, nbands, nkpoints, nspin, diri, dirj)
        eigs: The eigenvalues (nbands, nkpoints, nspin)
        kweights: The kpoint weights (nkpoints)
        efermi: The fermi energy
        nedos: The sampling of the energy values
        deltae: The energy grid spacing
        ismear: The smearing parameter used by the ``step_func``.
        sigma: The width of the smearing

    Return:
        np.array: Array of size `nedos` with the imaginary part of the dielectric function.

    """
    norm_kweights = np.array(kweights) / np.sum(kweights)
    egrid = np.arange(0, nedos * deltae, deltae)
    eigs_shifted = eigs - efermi
    # np.subtract.outer results in a matrix of shape (nband, nband)
    rspin = 3 - cder.shape[3]

    # for the transition between two bands at one kpoint the contributions is:
    #  (fermi[band_i] - fermi[band_j]) * rspin * normalized_kpoint_weight
    epsdd = np.zeros_like(egrid, dtype=np.complex128)
    for ib, jb, ik, ispin in tqdm(np.ndindex(cder.shape[:4]), total=np.prod(cder.shape[:4])):
        # print(f"ib={ib}, jb={jb}, ik={ik}, ispin={ispin}")
        fermi_w_i = step_func((eigs_shifted[ib, ik, ispin]) / sigma, ismear)
        fermi_w_j = step_func((eigs_shifted[jb, ik, ispin]) / sigma, ismear)
        weight = (fermi_w_j - fermi_w_i) * rspin * norm_kweights[ik]
        decel = eigs[jb, ik, ispin] - eigs[ib, ik, ispin]
        A = cder[ib, jb, ik, ispin, idir] * np.conjugate(cder[ib, jb, ik, ispin, jdir])
        # Reproduce the `SLOT` function calls in VASP:
        # CALL SLOT( REAL(DECEL,q), ISMEAR, SIGMA, NEDOS, DELTAE,  WEIGHT*A*CONST, EPSDD)
        # The conjugate part is not needed since we are running over all pairs of ib, jb
        # vasp just does the conjugate trick to save loop time
        smeared = get_delta(x0=decel, sigma=sigma, nx=nedos, dx=deltae, ismear=ismear) * weight * A
        epsdd += smeared
    return egrid, epsdd


def kramers_kronig(
    eps: npt.ArrayLike,
    nedos: int,
    deltae: float,
    cshift: float = 0.1,
) -> npt.NDArray:
    """Perform the Kramers-Kronig transformation.

    Perform the Kramers-Kronig transformation exactly as VASP does it.
    The input eps should be complex and the imaginary part of the dielectric function
    should be stored as the real part of the complex input array.
    The output should be the complex dielectric function.

    Args:
        eps: The dielectric function with the imaginary part stored as the real part and nothing in the imaginary part.
        nedos: The sampling of the energy values
        deltae: The energy grid spacing
        cshift: The shift of the imaginary part of the dielectric function.

    Return:
        np.array: Array of size `nedos` with the complex dielectric function.
    """
    egrid = np.linspace(0, deltae * nedos, nedos)
    csfhit = cshift * 1.0j
    cdiff = np.subtract.outer(egrid, egrid) + csfhit
    csum = np.add.outer(egrid, egrid) + csfhit
    vals = -0.5 * ((eps / cdiff) - (np.conj(eps) / csum))
    return np.sum(vals, axis=1) * 2 / np.pi * deltae
