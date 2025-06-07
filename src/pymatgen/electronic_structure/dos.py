"""This module defines classes to represent the density of states (DOS), etc."""

from __future__ import annotations

import functools
import warnings
from typing import TYPE_CHECKING, NamedTuple, cast

import numpy as np
from monty.json import MSONable
from packaging import version
from scipy.constants import value as _constant
from scipy.ndimage import gaussian_filter1d
from scipy.signal import hilbert
from scipy.special import expit
from scipy.stats import wasserstein_distance

from pymatgen.core import Structure, get_el_sp
from pymatgen.core.spectrum import Spectrum
from pymatgen.electronic_structure.core import Orbital, OrbitalType, Spin
from pymatgen.util.coord import get_linear_interpolated_value

if version.parse(np.__version__) < version.parse("2.0.0"):
    np.trapezoid = np.trapz  # type:ignore[assignment] # noqa: NPY201

if TYPE_CHECKING:
    from collections.abc import Mapping
    from typing import Any, Literal

    from numpy.typing import ArrayLike, NDArray
    from typing_extensions import Self

    from pymatgen.core.sites import PeriodicSite
    from pymatgen.util.typing import SpeciesLike


class DOS(Spectrum):
    """(Work in progress) Replacement of basic DOS object.
        All other DOS objects are extended versions of this.

    Attributes:
        energies (Sequence[float]): Energies.
        densities (dict[Spin, NDArray]): Spin densities,
            e.g. {Spin.up: DOS_up, Spin.down: DOS_down}.
        efermi (float): The Fermi level.
    """

    XLABEL = "Energy"
    YLABEL = "Density"

    def __init__(self, energies: ArrayLike, densities: ArrayLike, efermi: float) -> None:
        """
        Args:
            energies (Sequence[float]): The Energies.
            densities (NDArray): A Nx1 or Nx2 array. If former, it is
                interpreted as a Spin.up only density. Otherwise, the first column
                is interpreted as Spin.up and the other Spin.down.
            efermi (float): The Fermi level.
        """
        super().__init__(energies, densities, efermi)
        self.efermi = efermi

    def __str__(self) -> str:
        """A string which can be easily plotted."""
        if Spin.down in self.densities:
            str_arr = [f"#{'Energy':30s} {'DensityUp':30s} {'DensityDown':30s}"]
            for idx, energy in enumerate(self.energies):
                str_arr.append(f"{energy:.5f} {self.densities[Spin.up][idx]:.5f} {self.densities[Spin.down][idx]:.5f}")

        else:
            str_arr = [f"#{'Energy':30s} {'DensityUp':30s}"]
            for idx, energy in enumerate(self.energies):
                str_arr.append(f"{energy:.5f} {self.densities[Spin.up][idx]:.5f}")

        return "\n".join(str_arr)

    def get_interpolated_gap(
        self,
        tol: float = 1e-4,
        abs_tol: bool = False,
        spin: Spin | None = None,
    ) -> tuple[float, float, float]:
        """Find the interpolated band gap.

        Args:
            tol (float): Tolerance in occupations for determining the gap.
            abs_tol (bool): Use absolute (True) or relative (False) tolerance.
            spin (Spin | None): Find the gap:
                - None: In the summed DOS.
                - Up: In the spin up channel.
                - Down: In the spin down channel.

        Returns:
            tuple[float, float, float]: Energies in eV corresponding to the
                band gap, CBM and VBM.
        """
        if spin is None:
            tdos = self.y if len(self.ydim) == 1 else np.sum(self.y, axis=1)
        elif spin == Spin.up:
            tdos = self.y[:, 0]
        else:
            tdos = self.y[:, 1]

        if not abs_tol:
            tol = tol * tdos.sum() / tdos.shape[0]

        energies = self.x
        below_fermi = [i for i in range(len(energies)) if energies[i] < self.efermi and tdos[i] > tol]
        above_fermi = [i for i in range(len(energies)) if energies[i] > self.efermi and tdos[i] > tol]
        if not below_fermi or not above_fermi:
            return 0.0, self.efermi, self.efermi

        vbm_start = max(below_fermi)
        cbm_start = min(above_fermi)
        if vbm_start in [cbm_start, cbm_start - 1]:
            return 0.0, self.efermi, self.efermi

        # Interpolate between adjacent values
        terminal_dens = tdos[vbm_start : vbm_start + 2][::-1]
        terminal_energies = energies[vbm_start : vbm_start + 2][::-1]
        vbm = get_linear_interpolated_value(terminal_dens, terminal_energies, tol)
        terminal_dens = tdos[cbm_start - 1 : cbm_start + 1]
        terminal_energies = energies[cbm_start - 1 : cbm_start + 1]
        cbm = get_linear_interpolated_value(terminal_dens, terminal_energies, tol)
        return cbm - vbm, cbm, vbm

    def get_cbm_vbm(self, tol: float = 1e-4, abs_tol: bool = False, spin: Spin | None = None) -> tuple[float, float]:
        """
        Expects a DOS object and finds the CBM and VBM eigenvalues,
        using interpolation to determine the points at which the
        DOS crosses the threshold `tol`.

        `tol` may need to be increased for systems with noise/disorder.

        Args:
            tol (float): Tolerance in occupations for determining the gap.
            abs_tol (bool): Use absolute (True) or relative (False) tolerance.
            spin (Spin | None): Find the gap:
                - None: In the summed DOS.
                - Up: In the spin up channel.
                - Down: In the spin down channel.

        Returns:
            tuple[float, float]: Energies in eV corresponding to the CBM and VBM.
        """
        _gap, cbm, vbm = self.get_interpolated_gap(tol, abs_tol, spin)
        return cbm, vbm

    def get_gap(self, tol: float = 1e-4, abs_tol: bool = False, spin: Spin | None = None) -> float:
        """
        Expects a DOS object and finds the band gap, using the determined
        VBM and CBM eigenvalues.

        `tol` may need to be increased for systems with noise/disorder.

        Args:
            tol (float): Tolerance in occupations for determining the gap.
            abs_tol (bool): Use absolute (True) or relative (False) tolerance.
            spin (Spin | None): Find the gap:
                - None: In the summed DOS.
                - Up: In the spin up channel.
                - Down: In the spin down channel.

        Returns:
            float: Gap in eV.
        """
        cbm, vbm = self.get_cbm_vbm(tol, abs_tol, spin)
        return max(cbm - vbm, 0.0)


class Dos(MSONable):
    """Basic DOS object. All other DOS objects are extended versions of this.

    Attributes:
        energies (Sequence[float]): Energies.
        densities (dict[Spin, NDArray): Spin densities,
            e.g. {Spin.up: DOS_up, Spin.down: DOS_down}.
        efermi (float): The Fermi level.
    """

    def __init__(
        self,
        efermi: float,
        energies: ArrayLike,
        densities: Mapping[Spin, ArrayLike],
        norm_vol: float | None = None,
    ) -> None:
        """
        Args:
            efermi (float): The Fermi level.
            energies (Sequence[float]): Energies.
            densities (dict[Spin, NDArray]): The density of states for each Spin.
            norm_vol (float | None): The volume used to normalize the DOS.
                Defaults to 1 if None which will not perform any normalization.
                If None, the result will be in unit of states/eV,
                otherwise will be in states/eV/Angstrom^3.
        """
        self.efermi = efermi
        self.energies = np.asarray(energies)
        self.norm_vol = norm_vol
        vol = norm_vol or 1
        self.densities = {k: np.asarray(d) / vol for k, d in densities.items()}

    def __add__(self, other):
        """Add two Dos.

        Args:
            other (Dos): Another Dos object.

        Raises:
            ValueError: If energy scales are different.

        Returns:
            Sum of the two Dos.
        """
        if not all(np.equal(self.energies, other.energies)):
            raise ValueError("Energies of both DOS are not compatible!")

        densities = {spin: self.densities[spin] + other.densities[spin] for spin in self.densities}
        return type(self)(self.efermi, self.energies, densities)

    def __str__(self) -> str:
        """A string which can be easily plotted."""
        if Spin.down in self.densities:
            str_arr = [f"#{'Energy':30s} {'DensityUp':30s} {'DensityDown':30s}"]
            for idx, energy in enumerate(self.energies):
                str_arr.append(f"{energy:.5f} {self.densities[Spin.up][idx]:.5f} {self.densities[Spin.down][idx]:.5f}")

        else:
            str_arr = [f"#{'Energy':30s} {'DensityUp':30s}"]
            for idx, energy in enumerate(self.energies):
                str_arr.append(f"{energy:.5f} {self.densities[Spin.up][idx]:.5f}")

        return "\n".join(str_arr)

    def get_densities(self, spin: Spin | None = None) -> None | NDArray:
        """Get the DOS for a particular spin.

        Args:
            spin (Spin): Spin.

        Returns:
            NDArray: The DOS for the particular spin. Or the sum of both spins
                if Spin is None.
        """
        if self.densities is None:
            return None

        if spin is not None:
            return self.densities[spin]

        if Spin.down in self.densities:
            return self.densities[Spin.up] + self.densities[Spin.down]

        return self.densities[Spin.up]

    def get_smeared_densities(self, sigma: float) -> dict[Spin, NDArray]:
        """Get the the DOS with a Gaussian smearing.

        Args:
            sigma (float): Standard deviation of Gaussian smearing.

        Returns:
            {Spin: NDArray}: Gaussian-smeared DOS by spin.
        """
        diff = [self.energies[idx + 1] - self.energies[idx] for idx in range(len(self.energies) - 1)]
        avg_diff = sum(diff) / len(diff)
        return {spin: gaussian_filter1d(dens, sigma / avg_diff) for spin, dens in self.densities.items()}

    def get_interpolated_value(self, energy: float) -> dict[Spin, float]:
        """Get interpolated density for a particular energy.

        Args:
            energy (float): Energy to return the density for.

        Returns:
            dict[Spin, float]: Density for energy for each spin.
        """
        return {
            spin: get_linear_interpolated_value(self.energies, self.densities[spin], energy) for spin in self.densities
        }

    def get_interpolated_gap(
        self,
        tol: float = 1e-4,
        abs_tol: bool = False,
        spin: Spin | None = None,
    ) -> tuple[float, float, float]:
        """Find the interpolated band gap.

        Args:
            tol (float): Tolerance in occupations for determining the band gap.
            abs_tol (bool): Use absolute (True) or relative (False) tolerance.
            spin (Spin | None): Find the gap:
                None - In the summed DOS.
                Up - In the spin up channel.
                Down - In the spin down channel.

        Returns:
            tuple[float, float, float]: Energies in eV corresponding to the
                band gap, CBM and VBM.
        """
        tdos = self.get_densities(spin)
        if tdos is None:
            raise ValueError("tdos is None")
        if not abs_tol:
            tol = tol * tdos.sum() / tdos.shape[0]

        energies = self.energies
        below_fermi = [i for i in range(len(energies)) if energies[i] < self.efermi and tdos[i] > tol]
        above_fermi = [i for i in range(len(energies)) if energies[i] > self.efermi and tdos[i] > tol]
        if not below_fermi or not above_fermi:
            return 0.0, self.efermi, self.efermi

        vbm_start = max(below_fermi)
        cbm_start = min(above_fermi)
        if vbm_start in [cbm_start, cbm_start - 1]:
            return 0.0, self.efermi, self.efermi

        # Interpolate between adjacent values
        terminal_dens = tdos[vbm_start : vbm_start + 2][::-1]
        terminal_energies = energies[vbm_start : vbm_start + 2][::-1]
        vbm = get_linear_interpolated_value(terminal_dens, terminal_energies, tol)

        terminal_dens = tdos[cbm_start - 1 : cbm_start + 1]
        terminal_energies = energies[cbm_start - 1 : cbm_start + 1]
        cbm = get_linear_interpolated_value(terminal_dens, terminal_energies, tol)

        return cbm - vbm, cbm, vbm

    def get_cbm_vbm(self, tol: float = 1e-4, abs_tol: bool = False, spin: Spin | None = None) -> tuple[float, float]:
        """
        Expects a DOS object and finds the CBM and VBM eigenvalues,
        using interpolation to determine the points at which the
        DOS crosses the threshold `tol`.

        `tol` may need to be increased for systems with noise/disorder.

        Args:
            tol (float): Tolerance in occupations for determining the gap.
            abs_tol (bool): Use absolute (True) or relative (False) tolerance.
            spin (Spin | None): Find the gap:
                None - In the summed DOS.
                Up - In the spin up channel.
                Down - In the spin down channel.

        Returns:
            tuple[float, float]: Energies in eV corresponding to the CBM and VBM.
        """
        _gap, cbm, vbm = self.get_interpolated_gap(tol, abs_tol, spin)
        return cbm, vbm

    def get_gap(
        self,
        tol: float = 1e-4,
        abs_tol: bool = False,
        spin: Spin | None = None,
    ) -> float:
        """Find the band gap.

        Args:
            tol (float): Tolerance in occupations for determining the band gap.
            abs_tol (bool): Use absolute (True) or relative (False) tolerance.
            spin (Spin | None): Find the band gap:
                None - In the summed DOS.
                Up - In the spin up channel.
                Down - In the spin down channel.

        Returns:
            float: Band gap in eV.
        """
        cbm, vbm = self.get_cbm_vbm(tol, abs_tol, spin)
        return max(cbm - vbm, 0.0)

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Get Dos from a dict representation."""
        return cls(
            dct["efermi"],
            dct["energies"],
            {Spin(int(k)): v for k, v in dct["densities"].items()},
        )

    def as_dict(self) -> dict[str, Any]:
        """JSON-serializable dict representation of Dos."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "efermi": self.efermi,
            "energies": self.energies.tolist(),
            "densities": {str(spin): list(dens) for spin, dens in self.densities.items()},
        }


class FermiDos(Dos, MSONable):
    """Relate the density of states, doping levels
    (i.e. carrier concentrations) and corresponding Fermi levels.

    A negative doping concentration indicates the majority carriers are
    electrons (N-type); a positive doping concentration indicates holes
    are the majority carriers (P-type).
    """

    def __init__(
        self,
        dos: Dos,
        structure: Structure | None = None,
        nelecs: float | None = None,
        bandgap: float | None = None,
    ) -> None:
        """
        Args:
            dos (Dos): Pymatgen Dos object.
            structure (Structure): A structure. If None, the Structure
                of the Dos will be used. If the Dos does not have an
                associated Structure, an ValueError will be raised.
            nelecs (float): The number of electrons included in the energy range of
                Dos. It is used for normalizing the DOS. Default None to
                the total number of electrons in the structure.
            bandgap (float): If not None, the energy values are scissored so that
                the electronic band gap matches this value.
        """
        super().__init__(
            dos.efermi,
            energies=dos.energies,
            densities={k: np.array(d) for k, d in dos.densities.items()},
        )

        if structure is None:
            if hasattr(dos, "structure"):
                structure = dos.structure
            else:
                raise ValueError("Structure object is not provided and not present in dos")

        self.structure = structure
        self.nelecs = nelecs or self.structure.composition.total_electrons

        self.volume = self.structure.volume
        self.energies = np.array(dos.energies)
        self.de = np.hstack((self.energies[1:], self.energies[-1])) - self.energies

        # Normalize total density of states based on integral at 0 K
        tdos = np.array(self.get_densities())
        self.tdos = tdos * self.nelecs / (tdos * self.de)[self.energies <= self.efermi].sum()

        ecbm, evbm = self.get_cbm_vbm()
        self.idx_vbm = int(np.argmin(abs(self.energies - evbm)))
        self.idx_cbm = int(np.argmin(abs(self.energies - ecbm)))
        self.idx_mid_gap = int(self.idx_vbm + (self.idx_cbm - self.idx_vbm) / 2)
        self.A_to_cm = 1e-8

        if bandgap:
            eref = self.efermi if evbm < self.efermi < ecbm else (evbm + ecbm) / 2.0

            idx_fermi = int(np.argmin(abs(self.energies - eref)))

            if idx_fermi == self.idx_vbm:
                # Fermi level and VBM should have different indices
                idx_fermi += 1

            self.energies[:idx_fermi] -= (bandgap - (ecbm - evbm)) / 2.0
            self.energies[idx_fermi:] += (bandgap - (ecbm - evbm)) / 2.0

    def get_doping(self, fermi_level: float, temperature: float) -> float:
        """
        Calculate the doping (majority carrier concentration) at a given
        Fermi level and temperature. A simple Left Riemann sum is used for
        integrating the density of states over energy & equilibrium Fermi-Dirac
        distribution.

        Args:
            fermi_level (float): The Fermi level in eV.
            temperature (float): The temperature in Kelvin.

        Returns:
            float: The doping concentration in units of 1/cm^3. Negative values
                indicate that the majority carriers are electrons (N-type),
                whereas positive values indicates the majority carriers are holes
                (P-type).
        """
        cb_integral = np.sum(
            self.tdos[max(self.idx_mid_gap, self.idx_vbm + 1) :]
            * f0(self.energies[max(self.idx_mid_gap, self.idx_vbm + 1) :], fermi_level, temperature)
            * self.de[max(self.idx_mid_gap, self.idx_vbm + 1) :],
            axis=0,
        )
        vb_integral = np.sum(
            self.tdos[: min(self.idx_mid_gap, self.idx_cbm - 1) + 1]
            * f0(-self.energies[: min(self.idx_mid_gap, self.idx_cbm - 1) + 1], -fermi_level, temperature)
            * self.de[: min(self.idx_mid_gap, self.idx_cbm - 1) + 1],
            axis=0,
        )
        return (vb_integral - cb_integral) / (self.volume * self.A_to_cm**3)

    def get_fermi(
        self,
        concentration: float,
        temperature: float,
        rtol: float = 0.01,
        nstep: int = 50,
        step: float = 0.1,
        precision: int = 8,
    ) -> float:
        """Find the Fermi level at which the doping concentration at the given
        temperature (T) is equal to concentration. An algorithm is used
        where the relative error is minimized by calculating the doping at a
        grid which continually becomes finer.

        Args:
            concentration (float): The doping concentration in 1/cm^3. Negative
                values represent N-type doping and positive values represent P-type.
            temperature (float): The temperature in Kelvin.
            rtol (float): The maximum acceptable relative error.
            nstep (int): The number of steps checked around a given Fermi level.
            step (float): The initial Energy step length when searching.
            precision (int): The decimal places of calculated Fermi level.

        Raises:
            ValueError: If the Fermi level cannot be found.

        Returns:
            float: The Fermi level in eV. Note that this is different from
                the default Dos.efermi.
        """
        fermi = self.efermi  # initialize target Fermi
        relative_error: list | NDArray = [float("inf")]
        for _ in range(precision):
            fermi_range = np.arange(-nstep, nstep + 1) * step + fermi
            calc_doping = np.array([self.get_doping(fermi_lvl, temperature) for fermi_lvl in fermi_range])
            relative_error = np.abs(calc_doping / concentration - 1.0)
            fermi = fermi_range[np.argmin(relative_error)]
            step /= 10.0

        if min(relative_error) > rtol:
            raise ValueError(f"Could not find fermi within {rtol:.1%} of {concentration=}")
        return fermi

    def get_fermi_interextrapolated(
        self,
        concentration: float,
        temperature: float,
        warn: bool = True,
        c_ref: float = 1e10,
        **kwargs,
    ) -> float:
        """Similar to get_fermi method except that when it fails to converge, an
        interpolated or extrapolated Fermi level is returned, with the assumption
        that the Fermi level changes linearly with log(abs(concentration)),
        and therefore must be used with caution.

        Args:
            concentration (float): The doping concentration in 1/cm^3. Negative
                value represents N-type doping and positive value represents P-type.
            temperature (float): The temperature in Kelvin.
            warn (bool): Whether to give a warning the first time the Fermi level
                cannot be found.
            c_ref (float): A doping concentration where get_fermi returns a
                value without error for both c_ref and -c_ref.
            **kwargs: Keyword arguments passed to the get_fermi function.

        Returns:
            float: The possibly interpolated or extrapolated Fermi level.
        """
        try:
            return self.get_fermi(concentration, temperature, **kwargs)
        except ValueError as exc:
            if warn:
                warnings.warn(str(exc), stacklevel=2)

            if abs(concentration) < c_ref:
                if abs(concentration) < 1e-10:
                    concentration = 1e-10

                # max(10, ...) is to avoid log(0<x<1) and log(1+x) where both are slow
                f2 = self.get_fermi_interextrapolated(
                    max(10, abs(concentration) * 10.0),
                    temperature,
                    warn=False,
                    **kwargs,
                )
                f1 = self.get_fermi_interextrapolated(
                    -max(10, abs(concentration) * 10.0),
                    temperature,
                    warn=False,
                    **kwargs,
                )
                c2 = np.log(abs(1 + self.get_doping(f2, temperature)))
                c1 = -np.log(abs(1 + self.get_doping(f1, temperature)))
                slope = (f2 - f1) / (c2 - c1)
                return f2 + slope * (np.sign(concentration) * np.log(abs(1 + concentration)) - c2)

            f_ref = self.get_fermi_interextrapolated(np.sign(concentration) * c_ref, temperature, warn=False, **kwargs)
            f_new = self.get_fermi_interextrapolated(concentration / 10.0, temperature, warn=False, **kwargs)
            clog = np.sign(concentration) * np.log(abs(concentration))
            c_new_log = np.sign(concentration) * np.log(abs(self.get_doping(f_new, temperature)))
            slope = (f_new - f_ref) / (c_new_log - np.sign(concentration) * 10.0)
            return f_new + slope * (clog - c_new_log)

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Get FermiDos object from a dict representation."""
        dos = Dos(
            dct["efermi"],
            dct["energies"],
            {Spin(int(k)): v for k, v in dct["densities"].items()},
        )
        return cls(dos, structure=Structure.from_dict(dct["structure"]), nelecs=dct["nelecs"])

    def as_dict(self) -> dict[str, Any]:
        """JSON-serializable dict representation of FermiDos."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "efermi": self.efermi,
            "energies": self.energies.tolist(),
            "densities": {str(spin): dens.tolist() for spin, dens in self.densities.items()},
            "structure": self.structure,
            "nelecs": self.nelecs,
        }


class DosFingerprint(NamedTuple):
    """
    Represents a Density of States (DOS) fingerprint.

    This named tuple is used to store information related to the Density of States (DOS)
    in a material. It includes the energies, densities, type, number of bins, and bin width.

    Args:
        energies: The energy values associated with the DOS.
        densities: The corresponding density values for each energy.
        fp_type: The type of DOS fingerprint.
        n_bins: The number of bins used in the fingerprint.
        bin_width: The width of each bin in the DOS fingerprint.
    """

    energies: np.ndarray
    densities: np.ndarray
    fp_type: str
    n_bins: int
    bin_width: float


class CompleteDos(Dos):
    """Define total DOS, and projected DOS (PDOS).

    Mainly used by pymatgen.io.vasp.Vasprun to create a complete DOS from
    a vasprun.xml file. You are unlikely to generate this object manually.

    Attributes:
        structure (Structure): Structure associated with the CompleteDos.
        pdos (dict[PeriodicSite, dict[Orbital, dict[Spin, NDArray]]]): PDOS.
    """

    def __init__(
        self,
        structure: Structure,
        total_dos: Dos,
        pdoss: Mapping[PeriodicSite, Mapping[Orbital, Mapping[Spin, ArrayLike]]],
        normalize: bool = False,
    ) -> None:
        """
        Args:
            structure (Structure): Structure associated with this DOS.
            total_dos (Dos): Total DOS for the structure.
            pdoss (dict): The PDOSs supplied as {Site: {Orbital: {Spin: Densities}}}.
            normalize (bool): Whether to normalize the DOS by the volume of
                the structure. If True, the units of the DOS are states/eV/Angstrom^3.
                Otherwise, the units are states/eV.
        """
        vol = structure.volume if normalize else None
        super().__init__(
            total_dos.efermi,
            energies=total_dos.energies,
            densities={k: np.array(d) for k, d in total_dos.densities.items()},
            norm_vol=vol,
        )
        self.pdos = pdoss
        self.structure = structure

    def __str__(self) -> str:
        return f"Complete DOS for {self.structure}"

    def get_normalized(self) -> Self:
        """Get normalized CompleteDos."""
        if self.norm_vol is not None:
            return self

        return type(self)(
            structure=self.structure,
            total_dos=self,
            pdoss=self.pdos,
            normalize=True,
        )

    def get_site_orbital_dos(self, site: PeriodicSite, orbital: Orbital) -> Dos:
        """Get the Dos for a particular orbital of a particular site.

        Args:
            site: Site in Structure associated with CompleteDos.
            orbital: Orbital in the site.

        Returns:
            Dos: for a particular orbital of a particular site.
        """
        return Dos(self.efermi, self.energies, self.pdos[site][orbital])

    def get_site_dos(self, site: PeriodicSite) -> Dos:
        """Get the total DOS for a site with all orbitals.

        Args:
            site (PeriodicSite): Site in Structure associated with CompleteDos.

        Returns:
            Dos: Total DOS for a site with all orbitals.
        """
        site_dos = functools.reduce(add_densities, self.pdos[site].values())
        return Dos(self.efermi, self.energies, site_dos)

    def get_site_spd_dos(self, site: PeriodicSite) -> dict[OrbitalType, Dos]:
        """Get orbital projected DOS of a particular site.

        Args:
            site (PeriodicSite): Site in Structure associated with CompleteDos.

        Returns:
            dict[OrbitalType, Dos]
        """
        spd_dos: dict[OrbitalType, dict[Spin, np.ndarray]] = {}
        for orb, pdos in self.pdos[site].items():
            orbital_type = _get_orb_type(orb)
            if orbital_type in spd_dos:
                spd_dos[orbital_type] = add_densities(spd_dos[orbital_type], pdos)
            else:
                spd_dos[orbital_type] = pdos  # type: ignore[assignment]
        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in spd_dos.items()}

    def get_site_t2g_eg_resolved_dos(
        self,
        site: PeriodicSite,
    ) -> dict[Literal["e_g", "t2g"], Dos]:
        """Get the t2g/e_g projected DOS for a particular site.

        Args:
            site (PeriodicSite): Site in Structure associated with CompleteDos.

        Returns:
            dict[Literal["e_g", "t2g"], Dos]: Summed e_g and t2g DOS for the site.
        """
        t2g_dos = []
        eg_dos = []
        for s, atom_dos in self.pdos.items():
            if s == site:
                for orb, pdos in atom_dos.items():
                    if orb in (Orbital.dxy, Orbital.dxz, Orbital.dyz):
                        t2g_dos.append(pdos)
                    elif orb in (Orbital.dx2, Orbital.dz2):
                        eg_dos.append(pdos)
        return {
            "t2g": Dos(self.efermi, self.energies, functools.reduce(add_densities, t2g_dos)),
            "e_g": Dos(self.efermi, self.energies, functools.reduce(add_densities, eg_dos)),
        }

    def get_spd_dos(self) -> dict[OrbitalType, Dos]:
        """Get orbital projected DOS.

        Returns:
            dict[OrbitalType, Dos]
        """
        spd_dos = {}
        for atom_dos in self.pdos.values():
            for orb, pdos in atom_dos.items():
                orbital_type = _get_orb_type(orb)
                if orbital_type not in spd_dos:
                    spd_dos[orbital_type] = pdos
                else:
                    spd_dos[orbital_type] = add_densities(spd_dos[orbital_type], pdos)
        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in spd_dos.items()}

    def get_element_dos(self) -> dict[SpeciesLike, Dos]:
        """Get element projected DOS.

        Returns:
            dict[Element, Dos]
        """
        el_dos: dict[SpeciesLike, dict[Spin, ArrayLike]] = {}
        for site, atom_dos in self.pdos.items():
            el = site.specie
            for pdos in atom_dos.values():
                el_dos[el] = add_densities(el_dos[el], pdos) if el in el_dos else pdos  # type: ignore[assignment]

        return {el: Dos(self.efermi, self.energies, densities) for el, densities in el_dos.items()}

    def get_element_spd_dos(self, el: SpeciesLike) -> dict[OrbitalType, Dos]:
        """Get element and orbital (spd) projected DOS.

        Args:
            el (SpeciesLike): Element associated with CompleteDos.

        Returns:
            dict[OrbitalType, Dos]
        """
        el = get_el_sp(el)
        el_dos = {}
        for site, atom_dos in self.pdos.items():
            if site.specie == el:
                for orb, pdos in atom_dos.items():
                    orbital_type = _get_orb_type(orb)
                    if orbital_type not in el_dos:
                        el_dos[orbital_type] = pdos
                    else:
                        el_dos[orbital_type] = add_densities(el_dos[orbital_type], pdos)

        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in el_dos.items()}

    @property
    def spin_polarization(self) -> float | None:
        """Spin polarization at Fermi level.

        Examples:
            See Sanvito et al., DOI: 10.1126/sciadv.1602241 for an example usage.

        Returns:
            float | None: Spin polarization in range [0, 1], will return NaN
                if spin polarization is ill-defined (e.g. for insulator).
                None if the calculation is not spin-polarized.
        """
        n_F = self.get_interpolated_value(self.efermi)

        n_F_up = n_F[Spin.up]
        if Spin.down not in n_F:
            return None

        n_F_down = n_F[Spin.down]
        # Only well defined for metals or half-metals
        if (n_F_up + n_F_down) == 0:
            return float("NaN")

        spin_polarization = (n_F_up - n_F_down) / (n_F_up + n_F_down)

        return abs(spin_polarization)

    def get_band_filling(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
    ) -> float:
        """Compute the orbital-projected band filling, defined as the zeroth moment
        up to the Fermi level.

        "elements" and "sites" cannot be used together.

        Args:
            band (OrbitalType): Orbital to get the band center of (default is d-band).
            elements (list[SpeciesLike]): Elements to get the band center of.
            sites (list[PeriodicSite]): Sites to get the band center of.
            spin (Spin): Spin channel to use. If None, both spin channels will be combined.

        Returns:
            float: Band filling in eV, often denoted f_d for the d-band.
        """
        # Get the projected DOS
        if elements and sites:
            raise ValueError("Both element and site cannot be specified.")

        densities: dict[Spin, NDArray] = {}
        if elements:
            for idx, el in enumerate(elements):
                spd_dos = self.get_element_spd_dos(el)[band]
                densities = spd_dos.densities if idx == 0 else add_densities(densities, spd_dos.densities)
            dos = Dos(self.efermi, self.energies, densities)

        elif sites:
            for idx, site in enumerate(sites):
                spd_dos = self.get_site_spd_dos(site)[band]
                densities = spd_dos.densities if idx == 0 else add_densities(densities, spd_dos.densities)
            dos = Dos(self.efermi, self.energies, densities)

        else:
            dos = self.get_spd_dos()[band]

        energies = dos.energies - dos.efermi
        dos_densities = dos.get_densities(spin=spin)
        if dos_densities is None:
            raise ValueError("dos_densities is None")

        # Only integrate up to Fermi level
        energies = dos.energies - dos.efermi
        return np.trapezoid(dos_densities[energies < 0], x=energies[energies < 0]) / np.trapezoid(
            dos_densities, x=energies
        )

    def get_band_center(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: tuple[float, float] | None = None,
    ) -> float:
        """Compute the orbital-projected band center, defined as the first moment
        relative to the Fermi level as:
            int_{-inf}^{+inf} rho(E)*E dE/int_{-inf}^{+inf} rho(E) dE.

        Note that the band center is often highly sensitive to the selected energy range.

        "elements" and "sites" cannot be used together.

        References:
            Hammer and Norskov, Surf. Sci., 343 (1995).

        Args:
            band (OrbitalType): Orbital to get the band center of (default is d-band)
            elements (list[SpeciesLike]): Elements to get the band center of.
            sites (list[PeriodicSite]): Sites to get the band center of.
            spin (Spin): Spin channel to use. If None, both spin channels will be combined.
            erange (tuple(min, max)): The energy range to consider, with respect to the
                Fermi level. Default to None for all energies.

        Returns:
            float: The band center in eV, often denoted epsilon_d for the d-band center.
        """
        return self.get_n_moment(
            1,
            elements=elements,
            sites=sites,
            band=band,
            spin=spin,
            erange=erange,
            center=False,
        )

    def get_band_width(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: tuple[float, float] | None = None,
    ) -> float:
        """Get the orbital-projected band width, defined as the square root
        of the second moment:
            sqrt(int_{-inf}^{+inf} rho(E)*(E-E_center)^2 dE/int_{-inf}^{+inf} rho(E) dE)
        where E_center is the orbital-projected band center.

        Note that the band width is often highly sensitive to the selected energy range.

        "elements" and "sites" cannot be used together.

        Args:
            band (OrbitalType): Orbital to get the band center of (default is d-band).
            elements (list[SpeciesLike]): Elements to get the band center of.
            sites (list[PeriodicSite]): Sites to get the band center of.
            spin (Spin): Spin channel to use. By default, both spin channels will be combined.
            erange (tuple(min, max)): The energy range to consider, with respect to the
                Fermi level. Default to None for all energies.

        Returns:
            float: Orbital-projected band width in eV.
        """
        return np.sqrt(self.get_n_moment(2, elements=elements, sites=sites, band=band, spin=spin, erange=erange))

    def get_band_skewness(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: tuple[float, float] | None = None,
    ) -> float:
        """Get the orbital-projected skewness, defined as the third standardized moment:
            int_{-inf}^{+inf} rho(E)*(E-E_center)^3 dE/int_{-inf}^{+inf} rho(E) dE)
            /
            (int_{-inf}^{+inf} rho(E)*(E-E_center)^2 dE/int_{-inf}^{+inf} rho(E) dE))^(3/2)
        where E_center is the orbital-projected band center.

        Note that the skewness is often highly sensitive to the selected energy range.

        "elements" and "sites" cannot be used together.

        Args:
            band (OrbitalType): Orbitals to get the band center of (default is d-band).
            elements (list[SpeciesLike]): Elements to get the band center of.
            sites (list[PeriodicSite]): Sites to get the band center of.
            spin (Spin): Spin channel to use. By default, both spin channels will be combined.
            erange (tuple(min, max)): The energy range to consider, with respect to the
                Fermi level. Default to None for all energies.

        Returns:
            float: The orbital-projected skewness (dimensionless).
        """
        kwds: dict = dict(elements=elements, sites=sites, band=band, spin=spin, erange=erange)
        return self.get_n_moment(3, **kwds) / self.get_n_moment(2, **kwds) ** (3 / 2)

    def get_band_kurtosis(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: tuple[float, float] | None = None,
    ) -> float:
        """Get the orbital-projected kurtosis, defined as the fourth standardized moment
            int_{-inf}^{+inf} rho(E)*(E-E_center)^4 dE/int_{-inf}^{+inf} rho(E) dE)
            /
            (int_{-inf}^{+inf} rho(E)*(E-E_center)^2 dE/int_{-inf}^{+inf} rho(E) dE))^2
        where E_center is the orbital-projected band center.

        Note that the kurtosis is often highly sensitive to the selected energy range.

        "elements" and "sites" cannot be used together.

        Args:
            band (OrbitalType): Orbitals to get the band center of (default is d-band).
            elements (list[SpeciesLike]): Elements to get the band center of.
            sites (list[PeriodicSite]): Sites to get the band center of.
            spin (Spin): Spin channel to use. By default, both spin channels will be combined.
            erange (tuple(min, max)): The energy range to consider, with respect to the
                Fermi level. Default to None for all energies.

        Returns:
            float: The orbital-projected kurtosis (dimensionless).
        """
        kwds: dict = dict(elements=elements, sites=sites, band=band, spin=spin, erange=erange)
        return self.get_n_moment(4, **kwds) / self.get_n_moment(2, **kwds) ** 2

    def get_n_moment(
        self,
        n: int,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: tuple[float, float] | None = None,
        center: bool = True,
    ) -> float:
        """Get the nth moment of the DOS centered around the orbital-projected
        band center, defined as:
            int_{-inf}^{+inf} rho(E)*(E-E_center)^n dE/int_{-inf}^{+inf} rho(E) dE
        where n is the order, E_center is the orbital-projected band center, and
        E is the set of energies taken with respect to the Fermi level.

        "elements" and "sites" cannot be used together.

        Args:
            n (int): The order for the moment.
            band (OrbitalType): Orbital to get the band center of (default is d-band).
            elements (list[PeriodicSite]): Elements to get the band center of.
            sites (list[PeriodicSite]): Sites to get the band center of.
            spin (Spin): Spin channel to use. By default, both spin channels will be combined.
            erange (tuple(min, max)): The energy range to consider, with respect to the
                Fermi level. Default to None for all energies.
            center (bool): Take moments with respect to the band center.

        Returns:
            Orbital-projected nth moment in eV
        """
        # Get the projected DOS
        if elements and sites:
            raise ValueError("Both element and site cannot be specified.")

        densities: dict[Spin, NDArray] = {}
        if elements:
            for idx, el in enumerate(elements):
                spd_dos = self.get_element_spd_dos(el)[band]
                densities = spd_dos.densities if idx == 0 else add_densities(densities, spd_dos.densities)
            dos = Dos(self.efermi, self.energies, densities)

        elif sites:
            for idx, site in enumerate(sites):
                spd_dos = self.get_site_spd_dos(site)[band]
                densities = spd_dos.densities if idx == 0 else add_densities(densities, spd_dos.densities)
            dos = Dos(self.efermi, self.energies, densities)

        else:
            dos = self.get_spd_dos()[band]

        energies = dos.energies - dos.efermi
        dos_densities = dos.get_densities(spin=spin)
        if dos_densities is None:
            raise ValueError("dos_densities is None")

        # Only consider a given energy range
        if erange:
            dos_densities = dos_densities[(energies >= erange[0]) & (energies <= erange[1])]
            energies = energies[(energies >= erange[0]) & (energies <= erange[1])]

        # Center the energies about the band center if requested
        if center:
            band_center = self.get_band_center(elements=elements, sites=sites, band=band, spin=spin, erange=erange)
            p = energies - band_center
        else:
            p = energies

        # Take the nth moment
        return np.trapezoid(p**n * dos_densities, x=energies) / np.trapezoid(dos_densities, x=energies)

    def get_hilbert_transform(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
    ) -> Dos:
        """Get the Hilbert transform of the orbital-projected DOS,
        often plotted for a Newns-Anderson analysis.

        "elements" and "sites" cannot be used together.

        Args:
            band (OrbitalType): Orbital to get the band center of (default is d-band).
            elements (list[SpeciesLike]): Elements to get the band center of.
            sites (list[PeriodicSite]): Sites to get the band center of.

        Returns:
            Dos: Hilbert transformation of the projected DOS.
        """
        # Get the projected DOS
        if elements and sites:
            raise ValueError("Both element and site cannot be specified.")

        densities: dict[Spin, NDArray] = {}
        if elements:
            for idx, el in enumerate(elements):
                spd_dos = self.get_element_spd_dos(el)[band]
                densities = spd_dos.densities if idx == 0 else add_densities(densities, spd_dos.densities)
            dos = Dos(self.efermi, self.energies, densities)

        elif sites:
            for idx, site in enumerate(sites):
                spd_dos = self.get_site_spd_dos(site)[band]
                densities = spd_dos.densities if idx == 0 else add_densities(densities, spd_dos.densities)
            dos = Dos(self.efermi, self.energies, densities)

        else:
            dos = self.get_spd_dos()[band]

        # Get Hilbert-transformed densities
        densities_transformed = {Spin.up: np.imag(hilbert(dos.get_densities(spin=Spin.up)))}
        if Spin.down in self.densities:
            densities_transformed[Spin.down] = np.imag(hilbert(dos.get_densities(spin=Spin.down)))

        return Dos(self.efermi, self.energies, densities_transformed)

    def get_upper_band_edge(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: tuple[float, float] | None = None,
    ) -> float:
        """Get the orbital-projected upper band edge.

        The definition by Xin et al. Phys. Rev. B, 89, 115114 (2014) is used,
        which is the highest peak position of the Hilbert transform of
        the orbital-projected DOS.

        "elements" and "sites" cannot be used together.

        Args:
            band (OrbitalType): Orbital to get the band center of (default is d-band).
            elements (list[SpeciesLike]): Elements to get the band center of.
            sites (list[PeriodicSite]): Sites to get the band center of.
            spin (Spin): Spin channel to use. Both spin channels will be combined by default.
            erange (tuple(min, max)): The energy range to consider, with respect to the
                Fermi level. Default to None for all energies.

        Returns:
            float: Upper band edge in eV, often denoted epsilon_u.
        """
        # Get the Hilbert-transformed DOS
        transformed_dos = self.get_hilbert_transform(elements=elements, sites=sites, band=band)

        energies = transformed_dos.energies - transformed_dos.efermi
        densities = transformed_dos.get_densities(spin=spin)
        if densities is None:
            raise ValueError("densities is None")

        # Only consider a given energy range, if specified
        if erange:
            densities = densities[(energies >= erange[0]) & (energies <= erange[1])]
            energies = energies[(energies >= erange[0]) & (energies <= erange[1])]

        # Calculate the upper band edge
        return energies[np.argmax(densities)]

    def get_dos_fp(
        self,
        fp_type: str = "summed_pdos",
        binning: bool = True,
        min_e: float | None = None,
        max_e: float | None = None,
        n_bins: int = 256,
        normalize: bool = True,
    ) -> DosFingerprint:
        """Generate the DOS fingerprint.

        Based on the work of:
            F. Knoop, T. A. r Purcell, M. Scheffler, C. Carbogno, J. Open Source Softw. 2020, 5, 2671.
            Source - https://gitlab.com/vibes-developers/vibes/-/tree/master/vibes/materials_fp
            Copyright (c) 2020 Florian Knoop, Thomas A.R.Purcell, Matthias Scheffler, Christian Carbogno.
            Please also see and cite related work by:
            M. Kuban, S. Rigamonti, C. Draxl, Digital Discovery 2024, 3, 2448.

        Args:
            fp_type (str): The FingerPrint type, can be "{s/p/d/f/summed}_{pdos/tdos}"
                (default is summed_pdos).
            binning (bool): Whether to bin the DOS FingerPrint using np.linspace and n_bins.
                Default is True.
            min_e (float): The minimum energy to include (default is None).
            max_e (float): The maximum energy to include (default is None).
            n_bins (int): Number of bins to be used, if binning (default is 256).
            normalize (bool): Whether to normalize the integrated DOS to 1. Default is True.

        Raises:
            ValueError: If "fp_type" is not one of the accepted values.

        Returns:
            DosFingerprint(energies, densities, type, n_bins): The DOS fingerprint.
        """
        energies = self.energies - self.efermi

        if max_e is None:
            max_e = np.max(energies)

        if min_e is None:
            min_e = np.min(energies)

        pdos_obj = self.get_spd_dos()

        pdos = {key.name: pdos_obj[key].get_densities() for key in pdos_obj}

        pdos["summed_pdos"] = np.sum(list(pdos.values()), axis=0)  # type:ignore[arg-type]
        pdos["tdos"] = self.get_densities()

        try:
            densities = pdos[fp_type]
            if densities is None:
                raise ValueError("densities is None")
            if len(energies) < n_bins:
                inds = np.where((energies >= min_e) & (energies <= max_e))
                return DosFingerprint(
                    energies[inds],
                    densities[inds],
                    fp_type,
                    len(energies),
                    np.diff(energies)[0],
                )

            if binning:
                ener_bounds = np.linspace(min_e, max_e, n_bins + 1)
                ener = ener_bounds[:-1] + (ener_bounds[1] - ener_bounds[0]) / 2.0
                bin_width = np.diff(ener)[0]
            else:
                ener_bounds = np.array(energies)
                ener = np.append(energies, [energies[-1] + np.abs(energies[-1]) / 10])
                n_bins = len(energies)
                bin_width = np.diff(energies)[0]

            dos_rebin = np.zeros(ener.shape)

            for ii, e1, e2 in zip(range(len(ener)), ener_bounds[:-1], ener_bounds[1:], strict=False):
                inds = np.where((energies >= e1) & (energies < e2))
                dos_rebin[ii] = np.sum(densities[inds])

            # Scale DOS bins to make area under histogram equal 1
            if normalize:
                area = np.sum(dos_rebin * bin_width)
                dos_rebin_sc = dos_rebin / area
            else:
                dos_rebin_sc = dos_rebin

            return DosFingerprint(np.array([ener]), dos_rebin_sc, fp_type, n_bins, bin_width)

        except KeyError as exc:
            raise ValueError(
                "Please recheck fp_type requested, either the orbital projections unavailable in input DOS or "
                "there's a typo in type."
            ) from exc

    @staticmethod
    def fp_to_dict(fp: DosFingerprint) -> dict[str, NDArray]:
        """Convert a DOS FingerPrint into a dict.

        Args:
            fp (DosFingerprint): The DOS FingerPrint to convert.

        Returns:
            dict(Keys=type, Values=np.array(energies, densities)): The FingerPrint as dict.
        """
        return {fp[2]: np.array([fp[0], fp[1]], dtype="object").T}

    @staticmethod
    def get_dos_fp_similarity(
        fp1: DosFingerprint,
        fp2: DosFingerprint,
        col: int = 1,
        pt: int | Literal["All"] = "All",
        normalize: bool = False,
        metric: Literal["tanimoto", "wasserstein", "cosine-sim"] = "tanimoto",
    ) -> float:
        """Calculate the similarity index (dot product) of two DOS FingerPrints.

        Args:
            fp1 (DosFingerprint): The 1st dos fingerprint object
            fp2 (DosFingerprint): The 2nd dos fingerprint object
            col (int): The item in the fingerprints (0:energies,1: densities) to compute
                the similarity index of (default is 1)
            pt (int or str) : The index of the point that the dot product is to be taken (default is All)
            normalize (bool): If True normalize the scalar product to 1 (default is False)
            metric (Literal): Metric used to compute similarity default is "tanimoto".

        Raises:
            ValueError: If metric other than tanimoto, wasserstein and "cosine-sim" is requested.
            ValueError:  If normalize is set to True along with the metric.

        Returns:
            float: Similarity index given by the dot product.
        """
        valid_metrics = ("tanimoto", "wasserstein", "cosine-sim")
        if metric not in valid_metrics:
            raise ValueError(f"Invalid {metric=}, choose from {valid_metrics}.")

        fp1_dict = CompleteDos.fp_to_dict(fp1) if not isinstance(fp1, dict) else fp1
        fp2_dict = CompleteDos.fp_to_dict(fp2) if not isinstance(fp2, dict) else fp2

        if pt == "All":
            vec1 = np.array([pt[col] for pt in fp1_dict.values()]).flatten()
            vec2 = np.array([pt[col] for pt in fp2_dict.values()]).flatten()
        else:
            vec1 = fp1_dict[fp1[2][pt]][col]
            vec2 = fp2_dict[fp2[2][pt]][col]

        if not normalize and metric == "tanimoto":
            rescale = np.linalg.norm(vec1) ** 2 + np.linalg.norm(vec2) ** 2 - np.dot(vec1, vec2)
            return np.dot(vec1, vec2) / rescale

        if not normalize and metric == "wasserstein":
            return wasserstein_distance(
                u_values=np.cumsum(vec1 * fp1.bin_width),
                v_values=np.cumsum(vec2 * fp2.bin_width),
            )

        if normalize and metric == "cosine-sim":
            rescale = np.linalg.norm(vec1) * np.linalg.norm(vec2)
            return np.dot(vec1, vec2) / rescale

        if not normalize and metric == "cosine-sim":
            rescale = 1.0
            return np.dot(vec1, vec2) / rescale

        raise ValueError("Cannot compute similarity index. When normalize=True, then please set metric=cosine-sim")

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Get CompleteDos object from dict representation."""
        tdos = Dos.from_dict(dct)
        struct = Structure.from_dict(dct["structure"])
        pdoss = {}
        for idx in range(len(dct["pdos"])):
            at = struct[idx]
            orb_dos = {}
            for orb_str, odos in dct["pdos"][idx].items():
                orb = Orbital[orb_str]
                orb_dos[orb] = {Spin(int(k)): v for k, v in odos["densities"].items()}
            pdoss[at] = orb_dos
        return cls(struct, tdos, pdoss)

    def as_dict(self) -> dict[str, Any]:
        """JSON-serializable dict representation of CompleteDos."""
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "efermi": self.efermi,
            "structure": self.structure.as_dict(),
            "energies": self.energies.tolist(),
            "densities": {str(spin): dens.tolist() for spin, dens in self.densities.items()},
            "pdos": [],
        }
        if len(self.pdos) > 0:
            for at in self.structure:
                dd = {}
                for orb, pdos in self.pdos[at].items():
                    dd[str(orb)] = {"densities": {str(int(spin)): list(dens) for spin, dens in pdos.items()}}  # type:ignore[arg-type]
                dct["pdos"].append(dd)
            dct["atom_dos"] = {str(at): dos.as_dict() for at, dos in self.get_element_dos().items()}
            dct["spd_dos"] = {str(orb): dos.as_dict() for orb, dos in self.get_spd_dos().items()}
        return dct


_lobster_orb_labs = (
    "s",
    "p_y",
    "p_z",
    "p_x",
    "d_xy",
    "d_yz",
    "d_z^2",
    "d_xz",
    "d_x^2-y^2",
    "f_y(3x^2-y^2)",
    "f_xyz",
    "f_yz^2",
    "f_z^3",
    "f_xz^2",
    "f_z(x^2-y^2)",
    "f_x(x^2-3y^2)",
)


class LobsterCompleteDos(CompleteDos):
    """Extended CompleteDos for LOBSTER."""

    def get_site_orbital_dos(self, site: PeriodicSite, orbital: str) -> Dos:
        """Get the DOS for a particular orbital of a particular site.

        Args:
            site (PeriodicSite): Site in Structure associated with LobsterCompleteDos.
            orbital (str): Principal quantum number and orbital, e.g. "4s".
                    Possible orbitals are: "s", "p_y", "p_z", "p_x", "d_xy", "d_yz", "d_z^2",
                        "d_xz", "d_x^2-y^2", "f_y(3x^2-y^2)", "f_xyz",
                        "f_yz^2", "f_z^3", "f_xz^2", "f_z(x^2-y^2)", "f_x(x^2-3y^2)".
                    In contrast to the Cohpcar and the Cohplist objects,
                        the strings from the LOBSTER files are used.

        Returns:
            Dos: DOS of an orbital of a specific site.
        """
        if orbital[1:] not in _lobster_orb_labs:
            raise ValueError("orbital is not correct")

        return Dos(self.efermi, self.energies, self.pdos[site][orbital])  # type: ignore[index]

    def get_site_t2g_eg_resolved_dos(
        self,
        site: PeriodicSite,
    ) -> dict[Literal["e_g", "t2g"], Dos]:
        """Get the t2g/e_g projected DOS for a particular site.

        Args:
            site (PeriodicSite): Site in Structure associated with LobsterCompleteDos.

        Returns:
            dict[Literal["e_g", "t2g"], Dos]: Summed e_g and t2g DOS for the site.
        """
        warnings.warn("Are the orbitals correctly oriented? Are you sure?", stacklevel=2)

        t2g_dos = []
        eg_dos = []
        for s, atom_dos in self.pdos.items():
            if s == site:
                for orb, pdos in atom_dos.items():
                    orbital = _get_orb_lobster(str(orb))
                    if orbital is None:
                        raise ValueError("orbital is None")

                    if orbital in (Orbital.dxy, Orbital.dxz, Orbital.dyz):
                        t2g_dos.append(pdos)
                    elif orbital in (Orbital.dx2, Orbital.dz2):
                        eg_dos.append(pdos)
        return {
            "t2g": Dos(self.efermi, self.energies, functools.reduce(add_densities, t2g_dos)),
            "e_g": Dos(self.efermi, self.energies, functools.reduce(add_densities, eg_dos)),
        }

    def get_spd_dos(self) -> dict[str, Dos]:
        """Get orbital projected DOS.

        For example, if 3s and 4s are included in the basis of some element,
        they will be both summed in the orbital projected DOS.

        Returns:
            {orbital: Dos}
        """
        spd_dos = {}
        orb = None
        for atom_dos in self.pdos.values():
            for orb, pdos in atom_dos.items():
                orbital_type = _get_orb_type_lobster(str(orb))
                if orbital_type not in spd_dos:
                    spd_dos[orbital_type] = pdos
                else:
                    spd_dos[orbital_type] = add_densities(spd_dos[orbital_type], pdos)

        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in spd_dos.items()}  # type: ignore[misc]

    def get_element_spd_dos(self, el: SpeciesLike) -> dict[str, Dos]:
        """Get element and s/p/d projected DOS.

        Args:
            el (SpeciesLike): Element associated with LobsterCompleteDos.

        Returns:
            dict of {OrbitalType.s: Dos, OrbitalType.p: Dos, OrbitalType.d: Dos}
        """
        el = get_el_sp(el)
        el_dos = {}
        for site, atom_dos in self.pdos.items():
            if site.specie == el:
                for orb, pdos in atom_dos.items():
                    orbital_type = _get_orb_type_lobster(str(orb))
                    if orbital_type not in el_dos:
                        el_dos[orbital_type] = pdos
                    else:
                        el_dos[orbital_type] = add_densities(el_dos[orbital_type], pdos)

        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in el_dos.items()}  # type: ignore[misc]

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Get LobsterCompleteDos from a dict representation."""
        tdos = Dos.from_dict(dct)
        struct = Structure.from_dict(dct["structure"])
        pdos = {}
        for idx in range(len(dct["pdos"])):
            pdos[struct[idx]] = {
                orb_str: {Spin(int(k)): v for k, v in odos["densities"].items()}
                for orb_str, odos in dct["pdos"][idx].items()
            }
        return cls(struct, tdos, pdos)


def add_densities(
    density1: Mapping[Spin, ArrayLike],
    density2: Mapping[Spin, ArrayLike],
) -> dict[Spin, NDArray]:
    """Sum two DOS along each spin channel.

    Args:
        density1 (dict[Spin, NDArray]): First DOS.
        density2 (dict[Spin, NDArray]): Second DOS.

    Returns:
        dict[Spin, NDArray]
    """
    return {spin: np.array(density1[spin]) + np.array(density2[spin]) for spin in density1}


def _get_orb_type(orb: Orbital | OrbitalType) -> OrbitalType:
    """Get OrbitalType."""
    try:
        return cast("Orbital", orb).orbital_type
    except AttributeError:
        return cast("OrbitalType", orb)


def f0(E: float | NDArray, fermi: float, T: float) -> float:
    """Fermi-Dirac distribution function.

    Args:
        E (float): Energy in eV.
        fermi (float): The Fermi level in eV.
        T (float): The temperature in kelvin.

    Returns:
        float: The Fermi-Dirac occupation probability at energy E.
    """
    exponent = (E - fermi) / (_constant("Boltzmann constant in eV/K") * T)
    return expit(-exponent)  # scipy logistic sigmoid function; expit(x) = 1/(1+exp(-x))


def _get_orb_type_lobster(orb: str) -> OrbitalType | None:
    """Get OrbitalType from str representation of the orbital.

    Args:
        orb (str): String representation of the orbital.

    Returns:
        OrbitalType
    """
    try:
        orbital = Orbital(_lobster_orb_labs.index(orb[1:]))
        return orbital.orbital_type

    except AttributeError:
        print("Orb not in list")
    return None


def _get_orb_lobster(orb: str) -> Orbital | None:
    """Get Orbital from str representation of the orbital.

    Args:
        orb (str): String representation of the orbital.

    Returns:
        pymatgen.electronic_structure.core.Orbital
    """
    try:
        return Orbital(_lobster_orb_labs.index(orb[1:]))

    except AttributeError:
        print("Orb not in list")
        return None
