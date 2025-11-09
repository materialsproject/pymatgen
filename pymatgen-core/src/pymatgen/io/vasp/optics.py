"""Classes for parsing and manipulating VASP optical properties calculations."""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import TYPE_CHECKING, overload

import numpy as np
import scipy.constants
import scipy.special
from monty.json import MSONable
from tqdm import tqdm

from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun, Waveder

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray
    from typing_extensions import Self

    from pymatgen.util.typing import PathLike

__author__ = "Jimmy-Xuan Shen"
__copyright__ = "Copyright 2022, The Materials Project"
__maintainer__ = "Jimmy-Xuan Shen"
__email__ = "jmmshn@gmail.com"

au2ang: float = scipy.constants.physical_constants["atomic unit of length"][0] / 1e-10
ryd2ev: float = scipy.constants.physical_constants["Rydberg constant times hc in eV"][0]
edeps: float = 4 * np.pi * 2 * ryd2ev * au2ang  # from constant.inc in VASP

KB: float = scipy.constants.physical_constants["Boltzmann constant in eV/K"][0]


@dataclass
class DielectricFunctionCalculator(MSONable):
    """Post-process VASP optical properties calculations.

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
    def from_vasp_objects(cls, vrun: Vasprun, waveder: Waveder) -> Self:
        """Construct a DielectricFunction from Vasprun, Kpoint, and Waveder.

        Args:
            vrun: Vasprun object
            kpoint: Kpoint object
            waveder: Waveder object
        """
        bands = vrun.eigenvalues
        if bands is None:
            raise RuntimeError("eigenvalues cannot be None.")

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

        if efermi is None:
            raise ValueError("efermi cannot be None.")

        return cls(
            cder_real=waveder.cder_real,
            cder_imag=waveder.cder_imag,
            eigs=eigs,
            kweights=kweights,  # type:ignore[arg-type]
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
    def from_directory(cls, directory: PathLike) -> Self:
        """Construct a DielectricFunction from a directory containing
        vasprun.xml and WAVEDER files.
        """

        def _try_reading(dtypes):
            """Return None if failed."""
            for dtype in dtypes:
                try:
                    return Waveder.from_binary(f"{directory}/WAVEDER", data_type=dtype)
                except ValueError as exc:
                    if "reshape" in str(exc):
                        continue
                    raise
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

        @overload
        def _use_default(param: int | None, default: int) -> int:
            pass

        @overload
        def _use_default(param: float | None, default: float) -> float:
            pass

        def _use_default(param: float | None, default: float) -> float:
            return param if param is not None else default

        _efermi = _use_default(efermi, self.efermi)
        _nedos = _use_default(nedos, self.nedos)
        _deltae = _use_default(deltae, self.deltae)
        _ismear = _use_default(ismear, self.ismear)
        _sigma = _use_default(sigma, self.sigma)
        _cshift = _use_default(cshift, self.cshift)

        egrid, eps_imag = epsilon_imag(
            cder=self.cder,
            eigs=self.eigs,
            kweights=self.kweights,
            efermi=_efermi,
            nedos=_nedos,
            deltae=_deltae,
            ismear=_ismear,
            sigma=_sigma,
            idir=idir,
            jdir=jdir,
            mask=mask,
        )
        # scaling constant: edeps * np.pi / structure.volume
        eps_in = eps_imag * edeps * np.pi / self.volume
        eps = kramers_kronig(eps_in, nedos=_nedos, deltae=_deltae, cshift=_cshift)
        if idir == jdir:
            eps += 1.0 + 0.0j
        return egrid, eps

    def plot_weighted_transition_data(
        self,
        idir: int,
        jdir: int,
        mask: NDArray | None = None,
        min_val: float = 0.0,
    ):
        """Data for plotting the weight matrix elements as a scatter plot.

        Since the computation of the final spectrum (especially the smearing part)
        is still fairly expensive.  This function can be used to check the values
        of some portion of the spectrum (defined by the mask).
        In a sense, we are looking at the imaginary part of the dielectric function
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

        # Limit the first two indices based on the mask
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


def delta_methfessel_paxton(x: NDArray, n: int) -> NDArray:
    """
    D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
    where H is a Hermite polynomial and
    A_i = (-1)^i / ( i! 4^i sqrt(pi) ).
    """
    ii = np.arange(0, n + 1)
    A = (-1) ** ii / (scipy.special.factorial(ii) * 4**ii * np.sqrt(np.pi))
    H = scipy.special.eval_hermite(ii * 2, np.tile(x, (len(ii), 1)).T)
    return np.exp(-(x * x)) * np.dot(A, H.T)


def step_methfessel_paxton(x: NDArray, n: int) -> NDArray:
    """
    S_n (x) = (1 + erf x)/2 - exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
    where H is a Hermite polynomial and
    A_i = (-1)^i / ( i! 4^i sqrt(pi) ).
    """
    ii = np.arange(1, n + 1)
    A = (-1) ** ii / (scipy.special.factorial(ii) * 4**ii * np.sqrt(np.pi))
    H = scipy.special.eval_hermite(ii * 2 - 1, np.tile(x, (len(ii), 1)).T)
    return (1.0 + scipy.special.erf(x)) / 2.0 - np.exp(-(x * x)) * np.dot(A, H.T)


def delta_func(x: NDArray, ismear: int) -> NDArray:
    """Replication of VASP's delta function."""
    if ismear < -1:
        raise ValueError("Delta function not implemented for ismear < -1")
    if ismear == -1:
        return step_func(x, -1) * (1 - step_func(x, -1))
    if ismear == 0:
        return np.exp(-(x * x)) / np.sqrt(np.pi)
    return delta_methfessel_paxton(x, ismear)


def step_func(x: NDArray, ismear: int) -> NDArray:
    """Replication of VASP's step function."""
    if ismear < -1:
        raise ValueError("Delta function not implemented for ismear < -1")
    if ismear == -1:
        return 1 / (1.0 + np.exp(-x))
    if ismear == 0:
        return 0.5 + 0.5 * scipy.special.erf(x)
    return step_methfessel_paxton(x, ismear)


def get_delta(x0: float, sigma: float, nx: int, dx: float, ismear: int = 3) -> NDArray:
    """Get the smeared delta function to be added to form the spectrum.

    This replaces the `SLOT` function from VASP. Uses finite differences instead of
    evaluating the delta function since the step function is more likely to have analytic form.

    Args:
        x0: The center of the dielectric function.
        sigma: The width of the smearing
        nx: The number of grid points in the output grid.
        dx: The gridspacing of the output grid.
        ismear: The smearing parameter used by the ``step_func``.

    Returns:
        np.ndarray: Array of size `nx` with delta function on the desired outputgrid.
    """
    xgrid = np.linspace(0, nx * dx, nx, endpoint=False)
    xgrid -= x0
    x_scaled = (xgrid + (dx / 2)) / sigma
    sfun = step_func(x_scaled, ismear)
    dfun = np.zeros_like(xgrid)
    dfun[1:] = (sfun[1:] - sfun[:-1]) / dx
    return dfun


def get_step(x0: float, sigma: float, nx: int, dx: float, ismear: int) -> NDArray:
    """Get the smeared step function to be added to form the spectrum.

    This replaces the `SLOT` function from VASP.

    Args:
        x0: The center of the dielectric function.
        sigma: The width of the smearing
        nx: The number of grid points in the output grid.
        dx: The gridspacing of the output grid.
        ismear: The smearing parameter used by the ``step_func``.

    Returns:
        np.ndarray: Array of size `nx` with step function on the desired outputgrid.
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
) -> tuple[NDArray, NDArray]:
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

    Returns:
        np.ndarray: Array of size `nedos` with the imaginary part of the dielectric function.
    """
    norm_kweights = np.array(kweights) / np.sum(kweights)
    egrid = np.linspace(0, nedos * deltae, nedos, endpoint=False)
    eigs_shifted = eigs - efermi
    # np.subtract.outer results in a matrix of shape (nband, nband)
    rspin = 3 - cder.shape[3]

    # For the transition between two bands at one kpoint the contributions is:
    #  (fermi[band_i] - fermi[band_j]) * rspin * normalized_kpoint_weight
    cderm = cder * mask if mask is not None else cder

    # min_band0, max_band0 = np.min(np.where(cderm)[0]), np.max(np.where(cderm)[0])
    # min_band1, max_band1 = np.min(np.where(cderm)[1]), np.max(np.where(cderm)[1])
    # limit the first two indices based on the mask
    try:
        min_band0, max_band0 = np.min(np.where(cderm)[0]), np.max(np.where(cderm)[0])
        min_band1, max_band1 = np.min(np.where(cderm)[1]), np.max(np.where(cderm)[1])
    except ValueError as exc:
        if "zero-size array" in str(exc):
            return egrid, np.zeros_like(egrid, dtype=np.complex128)
        raise
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
    eps: NDArray,
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
        eps: The dielectric function with the imaginary part stored as the real part
            and nothing in the imaginary part.
        nedos: The sampling of the energy values
        deltae: The energy grid spacing
        cshift: The shift of the imaginary part of the dielectric function.

    Returns:
        np.ndarray: Array of size `nedos` with the complex dielectric function.
    """
    egrid = np.linspace(0, deltae * nedos, nedos)
    csfhit = cshift * 1.0j
    cdiff = np.subtract.outer(egrid, egrid) + csfhit
    csum = np.add.outer(egrid, egrid) + csfhit
    vals = -0.5 * ((eps / cdiff) - (np.conj(eps) / csum))
    return np.sum(vals, axis=1) * 2 / np.pi * deltae
