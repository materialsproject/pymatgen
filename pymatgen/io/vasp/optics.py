"""Classes for parsing and manipulating VASP optical properties calculations."""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import scipy.constants
import scipy.special
from monty.json import MSONable
from tqdm import tqdm

from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun, Waveder

if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import ArrayLike, NDArray

__author__ = "Jimmy-Xuan Shen"
__copyright__ = "Copyright 2022, The Materials Project"
__maintainer__ = "Jimmy-Xuan Shen"
__email__ = "jmmshn@gmail.com"

au2ang = scipy.constants.physical_constants["atomic unit of length"][0] / 1e-10
ryd2ev = scipy.constants.physical_constants["Rydberg constant times hc in eV"][0]
edeps = 4 * np.pi * 2 * ryd2ev * au2ang  # from constant.inc in VASP

KB = scipy.constants.physical_constants["Boltzmann constant in eV/K"][0]


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

    Currently, this Calculator only works for ``ISYM=0`` calculations since we cannot guarantee that our
    externally defined symmetry operations are the same as VASP's. This can be fixed by printing the
    symmetry operators into the vasprun.xml file. If this happens in future versions of VASP,
    we can dramatically speed up the calculations here by considering only the irreducible kpoints.
    """

    cder_real: NDArray
    cder_imag: NDArray
    eigs: NDArray
    kweights: NDArray
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
            cder_real=waveder.cder_real,
            cder_imag=waveder.cder_imag,
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

        def _try_reading(dtypes):
            """Return None if failed."""
            for dtype in dtypes:
                try:
                    return Waveder.from_binary(f"{directory}/WAVEDER", data_type=dtype)
                except ValueError as e:
                    if "reshape" in str(e):
                        continue
                    raise e
            return None

        vrun = Vasprun(f"{directory}/vasprun.xml")
        if "gamma" in vrun.generator["subversion"].lower():
            waveder = _try_reading(["float64", "float32"])  # large one first should give value error
        else:
            waveder = _try_reading(["complex128", "complex64"])
        return cls.from_vasp_objects(vrun, waveder)

    @property
    def cder(self):
        """Complex CDER from WAVEDER."""
        return self.cder_real + self.cder_imag * 1.0j

    def get_epsilon(
        self,
        idir: int,
        jdir: int,
        efermi: float | None = None,
        nedos: int | None = None,
        deltae: float | None = None,
        ismear: int | None = None,
        sigma: float | None = None,
        cshift: float | None = None,
        mask: NDArray | None = None,
    ) -> tuple[NDArray, NDArray]:
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
            mask: Mask for the bands/kpoint/spin index to include in the calculation
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
            mask=mask,
        )
        # scaling constant: edeps * np.pi / structure.volume
        eps_in = eps_imag * edeps * np.pi / self.volume
        eps = kramers_kronig(eps_in, nedos=nedos, deltae=deltae, cshift=cshift)  # type: ignore
        if idir == jdir:
            eps += 1.0 + 0.0j
        return egrid, eps

    def plot_weighted_transition_data(self, idir: int, jdir: int, mask: NDArray | None = None, min_val: float = 0.0):
        """Data for plotting the weight matrix elements as a scatter plot.

        Since the computation of the final spectrum (especially the smearing part)
        is still fairly expensive.  This function can be used to check the values
        of some portion of the spectrum (defined by the mask).
        In a sense, we are lookin at the imaginary part of the dielectric function
        before the smearing is applied.

        Args:
            idir: First direction of the dielectric tensor.
            jdir: Second direction of the dielectric tensor.
            mask: Mask to apply to the CDER for the bands/kpoint/spin
                index to include in the calculation
            min_val: Minimum value below this value the matrix element will not be shown.
        """
        cderm = self.cder * mask if mask is not None else self.cder

        norm_kweights = np.array(self.kweights) / np.sum(self.kweights)
        eigs_shifted = self.eigs - self.efermi
        rspin = 3 - cderm.shape[3]
        # limit the first two indices based on the mask
        try:
            min_band0, max_band0 = np.min(np.where(cderm)[0]), np.max(np.where(cderm)[0])
            min_band1, max_band1 = np.min(np.where(cderm)[1]), np.max(np.where(cderm)[1])
        except ValueError as exc:
            if "zero-size array" in str(exc):
                raise ValueError("No matrix elements found. Check the mask.")
            raise

        x_val = []
        y_val = []
        text = []
        _, _, nk, nspin = cderm.shape[:4]
        iter_idx = [
            range(min_band0, max_band0 + 1),
            range(min_band1, max_band1 + 1),
            range(nk),
            range(nspin),
        ]
        num_ = (max_band0 - min_band0) * (max_band1 - min_band1) * nk * nspin
        for ib, jb, ik, ispin in tqdm(itertools.product(*iter_idx), total=num_):
            fermi_w_i = step_func((eigs_shifted[ib, ik, ispin]) / self.sigma, self.ismear)
            fermi_w_j = step_func((eigs_shifted[jb, ik, ispin]) / self.sigma, self.ismear)
            weight = (fermi_w_j - fermi_w_i) * rspin * norm_kweights[ik]
            A = cderm[ib, jb, ik, ispin, idir] * np.conjugate(cderm[ib, jb, ik, ispin, jdir])
            decel = self.eigs[jb, ik, ispin] - self.eigs[ib, ik, ispin]
            matrix_el = np.abs(A) * float(weight)  # can have negative weight due to fermi function
            if matrix_el > min_val:
                x_val.append(decel)
                y_val.append(matrix_el)
                text.append(f"s:{ispin}, k:{ik}, {ib} -> {jb} ({decel:.2f})")
        return x_val, y_val, text


def delta_methfessel_paxton(x, n):
    """
    D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
    where H is a Hermite polynomial and
    A_i = (-1)^i / ( i! 4^i sqrt(pi) ).
    """
    ii = np.arange(0, n + 1)
    A = (-1) ** ii / (scipy.special.factorial(ii) * 4**ii * np.sqrt(np.pi))
    H = scipy.special.eval_hermite(ii * 2, np.tile(x, (len(ii), 1)).T)
    return np.exp(-(x * x)) * np.dot(A, H.T)


def step_methfessel_paxton(x, n):
    """
    S_n (x) = (1 + erf x)/2 - exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
    where H is a Hermite polynomial and
    A_i = (-1)^i / ( i! 4^i sqrt(pi) ).
    """
    ii = np.arange(1, n + 1)
    A = (-1) ** ii / (scipy.special.factorial(ii) * 4**ii * np.sqrt(np.pi))
    H = scipy.special.eval_hermite(ii * 2 - 1, np.tile(x, (len(ii), 1)).T)
    return (1.0 + scipy.special.erf(x)) / 2.0 - np.exp(-(x * x)) * np.dot(A, H.T)


def delta_func(x, ismear):
    """Replication of VASP's delta function."""
    if ismear < -1:
        raise ValueError("Delta function not implemented for ismear < -1")
    if ismear == -1:
        return step_func(x, -1) * (1 - step_func(x, -1))
    if ismear == 0:
        return np.exp(-(x * x)) / np.sqrt(np.pi)
    return delta_methfessel_paxton(x, ismear)


def step_func(x, ismear):
    """Replication of VASP's step function."""
    if ismear < -1:
        raise ValueError("Delta function not implemented for ismear < -1")
    if ismear == -1:
        return 1 / (1.0 + np.exp(-x))
    if ismear == 0:
        return 0.5 + 0.5 * scipy.special.erf(x)
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
    xgrid = np.linspace(0, nx * dx, nx, endpoint=False)
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
    xgrid = np.linspace(0, nx * dx, nx, endpoint=False)
    xgrid -= x0
    x_scaled = (xgrid + (dx / 2)) / sigma
    return step_func(x_scaled, ismear)


def epsilon_imag(
    cder: NDArray,
    eigs: NDArray,
    kweights: ArrayLike,
    efermi: float,
    nedos: int,
    deltae: float,
    ismear: int,
    sigma: float,
    idir: int,
    jdir: int,
    mask: NDArray | None = None,
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
        idir: The first direction of the dielectric tensor
        jdir: The second direction of the dielectric tensor
        mask: Mask for the bands/kpoint/spin index to include in the calculation

    Return:
        np.array: Array of size `nedos` with the imaginary part of the dielectric function.

    """
    norm_kweights = np.array(kweights) / np.sum(kweights)
    egrid = np.linspace(0, nedos * deltae, nedos, endpoint=False)
    eigs_shifted = eigs - efermi
    # np.subtract.outer results in a matrix of shape (nband, nband)
    rspin = 3 - cder.shape[3]

    # for the transition between two bands at one kpoint the contributions is:
    #  (fermi[band_i] - fermi[band_j]) * rspin * normalized_kpoint_weight
    cderm = cder * mask if mask is not None else cder

    # min_band0, max_band0 = np.min(np.where(cderm)[0]), np.max(np.where(cderm)[0])
    # min_band1, max_band1 = np.min(np.where(cderm)[1]), np.max(np.where(cderm)[1])
    # limit the first two indices based on the mask
    try:
        min_band0, max_band0 = np.min(np.where(cderm)[0]), np.max(np.where(cderm)[0])
        min_band1, max_band1 = np.min(np.where(cderm)[1]), np.max(np.where(cderm)[1])
    except ValueError as e:
        if "zero-size array" in str(e):
            return egrid, np.zeros_like(egrid, dtype=np.complex_)
        raise e
    _, _, nk, nspin = cderm.shape[:4]
    iter_idx = [
        range(min_band0, max_band0 + 1),
        range(min_band1, max_band1 + 1),
        range(nk),
        range(nspin),
    ]
    num_ = (max_band0 - min_band0) * (max_band1 - min_band1) * nk * nspin
    epsdd = np.zeros_like(egrid, dtype=np.complex128)
    for ib, jb, ik, ispin in tqdm(itertools.product(*iter_idx), total=num_):
        fermi_w_i = step_func((eigs_shifted[ib, ik, ispin]) / sigma, ismear)
        fermi_w_j = step_func((eigs_shifted[jb, ik, ispin]) / sigma, ismear)
        weight = (fermi_w_j - fermi_w_i) * rspin * norm_kweights[ik]
        decel = eigs[jb, ik, ispin] - eigs[ib, ik, ispin]
        A = cderm[ib, jb, ik, ispin, idir] * np.conjugate(cderm[ib, jb, ik, ispin, jdir])
        # Reproduce the `SLOT` function calls in VASP:
        # CALL SLOT( REAL(DECEL,q), ISMEAR, SIGMA, NEDOS, DELTAE,  WEIGHT*A*CONST, EPSDD)
        # The conjugate part is not needed since we are running over all pairs of ib, jb
        # vasp just does the conjugate trick to save loop time
        smeared = get_delta(x0=decel, sigma=sigma, nx=nedos, dx=deltae, ismear=ismear) * weight * A
        epsdd += smeared
    return egrid, epsdd


def kramers_kronig(
    eps: np.ndarray,
    nedos: int,
    deltae: float,
    cshift: float = 0.1,
) -> NDArray:
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

    # loop over that
