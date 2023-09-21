"""This module defines classes to represent the density of states, etc."""

from __future__ import annotations

import functools
import warnings
from collections import namedtuple
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
from monty.json import MSONable
from scipy.constants import value as _cd
from scipy.signal import hilbert

from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.spectrum import Spectrum
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Orbital, OrbitalType, Spin
from pymatgen.util.coord import get_linear_interpolated_value

if TYPE_CHECKING:
    from collections.abc import Mapping

    from numpy.typing import ArrayLike

    from pymatgen.core.sites import PeriodicSite
    from pymatgen.util.typing import SpeciesLike


class DOS(Spectrum):
    """Replacement basic DOS object. All other DOS objects are extended versions
    of this object. Work in progress.

    Attributes:
        energies (Sequence[float]): The sequence of energies.
        densities (dict[Spin, Sequence[float]]): A dict of spin densities, e.g., {Spin.up: [...], Spin.down: [...]}.
        efermi (float): Fermi level.
    """

    XLABEL = "Energy"
    YLABEL = "Density"

    def __init__(self, energies: ArrayLike, densities: ArrayLike, efermi: float):
        """
        Args:
            energies: A sequence of energies
            densities (ndarray): Either a Nx1 or a Nx2 array. If former, it is
                interpreted as a Spin.up only density. Otherwise, the first column
                is interpreted as Spin.up and the other is Spin.down.
            efermi: Fermi level energy.
        """
        super().__init__(energies, densities, efermi)
        self.efermi = efermi

    def get_interpolated_gap(self, tol: float = 0.001, abs_tol: bool = False, spin: Spin | None = None):
        """Expects a DOS object and finds the gap.

        Args:
            tol: tolerance in occupations for determining the gap
            abs_tol: Set to True for an absolute tolerance and False for a
                relative one.
            spin: Possible values are None - finds the gap in the summed
                densities, Up - finds the gap in the up spin channel,
                Down - finds the gap in the down spin channel.

        Returns:
            (gap, cbm, vbm):
                Tuple of floats in eV corresponding to the gap, cbm and vbm.
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
        vbm_start = max(below_fermi)
        cbm_start = min(above_fermi)
        if vbm_start == cbm_start:
            return 0.0, self.efermi, self.efermi
        # Interpolate between adjacent values
        terminal_dens = tdos[vbm_start : vbm_start + 2][::-1]
        terminal_energies = energies[vbm_start : vbm_start + 2][::-1]
        start = get_linear_interpolated_value(terminal_dens, terminal_energies, tol)
        terminal_dens = tdos[cbm_start - 1 : cbm_start + 1]
        terminal_energies = energies[cbm_start - 1 : cbm_start + 1]
        end = get_linear_interpolated_value(terminal_dens, terminal_energies, tol)
        return end - start, end, start

    def get_cbm_vbm(self, tol: float = 0.001, abs_tol: bool = False, spin=None):
        """Expects a DOS object and finds the cbm and vbm.

        Args:
            tol: tolerance in occupations for determining the gap
            abs_tol: An absolute tolerance (True) and a relative one (False)
            spin: Possible values are None - finds the gap in the summed
                densities, Up - finds the gap in the up spin channel,
                Down - finds the gap in the down spin channel.

        Returns:
            (cbm, vbm): float in eV corresponding to the gap
        """
        # determine tolerance
        if spin is None:
            tdos = self.y if len(self.ydim) == 1 else np.sum(self.y, axis=1)
        elif spin == Spin.up:
            tdos = self.y[:, 0]
        else:
            tdos = self.y[:, 1]

        if not abs_tol:
            tol = tol * tdos.sum() / tdos.shape[0]

        # find index of fermi energy
        i_fermi = 0
        while self.x[i_fermi] <= self.efermi:
            i_fermi += 1

        # work backwards until tolerance is reached
        i_gap_start = i_fermi
        while i_gap_start - 1 >= 0 and tdos[i_gap_start - 1] <= tol:
            i_gap_start -= 1

        # work forwards until tolerance is reached
        i_gap_end = i_gap_start
        while i_gap_end < tdos.shape[0] and tdos[i_gap_end] <= tol:
            i_gap_end += 1
        i_gap_end -= 1
        return self.x[i_gap_end], self.x[i_gap_start]

    def get_gap(self, tol: float = 0.001, abs_tol: bool = False, spin: Spin | None = None):
        """Expects a DOS object and finds the gap.

        Args:
            tol: tolerance in occupations for determining the gap
            abs_tol: An absolute tolerance (True) and a relative one (False)
            spin: Possible values are None - finds the gap in the summed
                densities, Up - finds the gap in the up spin channel,
                Down - finds the gap in the down spin channel.

        Returns:
            gap in eV
        """
        (cbm, vbm) = self.get_cbm_vbm(tol, abs_tol, spin)
        return max(cbm - vbm, 0.0)

    def __str__(self):
        """Returns a string which can be easily plotted (using gnuplot)."""
        if Spin.down in self.densities:
            stringarray = [f"#{'Energy':30s} {'DensityUp':30s} {'DensityDown':30s}"]
            for i, energy in enumerate(self.energies):
                stringarray.append(f"{energy:.5f} {self.densities[Spin.up][i]:.5f} {self.densities[Spin.down][i]:.5f}")
        else:
            stringarray = [f"#{'Energy':30s} {'DensityUp':30s}"]
            for i, energy in enumerate(self.energies):
                stringarray.append(f"{energy:.5f} {self.densities[Spin.up][i]:.5f}")
        return "\n".join(stringarray)


class Dos(MSONable):
    """Basic DOS object. All other DOS objects are extended versions of this
    object.

    Attributes:
        energies (Sequence[float]): The sequence of energies.
        densities (dict[Spin, Sequence[float]]): A dict of spin densities, e.g., {Spin.up: [...], Spin.down: [...]}.
        efermi (float): Fermi level.
    """

    def __init__(
        self, efermi: float, energies: ArrayLike, densities: Mapping[Spin, ArrayLike], norm_vol: float | None = None
    ) -> None:
        """
        Args:
            efermi: Fermi level energy
            energies: A sequences of energies
            densities (dict[Spin: np.array]): representing the density of states for each Spin.
            norm_vol: The volume used to normalize the densities. Defaults to 1 if None which will not perform any
                normalization. If not None, the resulting density will have units of states/eV/Angstrom^3, otherwise
                the density will be in states/eV.
        """
        self.efermi = efermi
        self.energies = np.array(energies)
        self.norm_vol = norm_vol
        vol = norm_vol or 1
        self.densities = {k: np.array(d) / vol for k, d in densities.items()}

    def get_densities(self, spin: Spin | None = None):
        """Returns the density of states for a particular spin.

        Args:
            spin: Spin

        Returns:
            Returns the density of states for a particular spin. If Spin is
            None, the sum of all spins is returned.
        """
        if self.densities is None:
            result = None
        elif spin is None:
            if Spin.down in self.densities:
                result = self.densities[Spin.up] + self.densities[Spin.down]
            else:
                result = self.densities[Spin.up]
        else:
            result = self.densities[spin]
        return result

    def get_smeared_densities(self, sigma: float):
        """Returns the Dict representation of the densities, {Spin: densities},
        but with a Gaussian smearing of std dev sigma.

        Args:
            sigma: Std dev of Gaussian smearing function.

        Returns:
            Dict of Gaussian-smeared densities.
        """
        from scipy.ndimage import gaussian_filter1d

        smeared_dens = {}
        diff = [self.energies[i + 1] - self.energies[i] for i in range(len(self.energies) - 1)]
        avgdiff = sum(diff) / len(diff)
        for spin, dens in self.densities.items():
            smeared_dens[spin] = gaussian_filter1d(dens, sigma / avgdiff)
        return smeared_dens

    def __add__(self, other):
        """Adds two DOS together. Checks that energy scales are the same.
        Otherwise, a ValueError is thrown.

        Args:
            other: Another DOS object.

        Returns:
            Sum of the two DOSs.
        """
        if not all(np.equal(self.energies, other.energies)):
            raise ValueError("Energies of both DOS are not compatible!")
        densities = {spin: self.densities[spin] + other.densities[spin] for spin in self.densities}
        return Dos(self.efermi, self.energies, densities)

    def get_interpolated_value(self, energy: float):
        """Returns interpolated density for a particular energy.

        Args:
            energy: Energy to return the density for.
        """
        f = {}
        for spin in self.densities:
            f[spin] = get_linear_interpolated_value(self.energies, self.densities[spin], energy)
        return f

    def get_interpolated_gap(self, tol: float = 0.001, abs_tol: bool = False, spin: Spin | None = None):
        """Expects a DOS object and finds the gap.

        Args:
            tol: tolerance in occupations for determining the gap
            abs_tol: Set to True for an absolute tolerance and False for a
                relative one.
            spin: Possible values are None - finds the gap in the summed
                densities, Up - finds the gap in the up spin channel,
                Down - finds the gap in the down spin channel.

        Returns:
            (gap, cbm, vbm):
                Tuple of floats in eV corresponding to the gap, cbm and vbm.
        """
        tdos = self.get_densities(spin)
        if not abs_tol:
            tol = tol * tdos.sum() / tdos.shape[0]
        energies = self.energies
        below_fermi = [i for i in range(len(energies)) if energies[i] < self.efermi and tdos[i] > tol]
        above_fermi = [i for i in range(len(energies)) if energies[i] > self.efermi and tdos[i] > tol]
        vbm_start = max(below_fermi)
        cbm_start = min(above_fermi)
        if vbm_start == cbm_start:
            return 0.0, self.efermi, self.efermi

        # Interpolate between adjacent values
        terminal_dens = tdos[vbm_start : vbm_start + 2][::-1]
        terminal_energies = energies[vbm_start : vbm_start + 2][::-1]
        start = get_linear_interpolated_value(terminal_dens, terminal_energies, tol)
        terminal_dens = tdos[cbm_start - 1 : cbm_start + 1]
        terminal_energies = energies[cbm_start - 1 : cbm_start + 1]
        end = get_linear_interpolated_value(terminal_dens, terminal_energies, tol)
        return end - start, end, start

    def get_cbm_vbm(self, tol: float = 0.001, abs_tol: bool = False, spin: Spin | None = None):
        """Expects a DOS object and finds the cbm and vbm.

        Args:
            tol: tolerance in occupations for determining the gap
            abs_tol: An absolute tolerance (True) and a relative one (False)
            spin: Possible values are None - finds the gap in the summed
                densities, Up - finds the gap in the up spin channel,
                Down - finds the gap in the down spin channel.

        Returns:
            (cbm, vbm): float in eV corresponding to the gap
        """
        # determine tolerance
        tdos = self.get_densities(spin)
        if not abs_tol:
            tol = tol * tdos.sum() / tdos.shape[0]

        # find index of fermi energy
        i_fermi = 0
        while self.energies[i_fermi] <= self.efermi:
            i_fermi += 1

        # work backwards until tolerance is reached
        i_gap_start = i_fermi
        while i_gap_start - 1 >= 0 and tdos[i_gap_start - 1] <= tol:
            i_gap_start -= 1

        # work forwards until tolerance is reached
        i_gap_end = i_gap_start
        while i_gap_end < tdos.shape[0] and tdos[i_gap_end] <= tol:
            i_gap_end += 1
        i_gap_end -= 1
        return self.energies[i_gap_end], self.energies[i_gap_start]

    def get_gap(self, tol: float = 0.001, abs_tol: bool = False, spin: Spin | None = None):
        """Expects a DOS object and finds the gap.

        Args:
            tol: tolerance in occupations for determining the gap
            abs_tol: An absolute tolerance (True) and a relative one (False)
            spin: Possible values are None - finds the gap in the summed
                densities, Up - finds the gap in the up spin channel,
                Down - finds the gap in the down spin channel.

        Returns:
            gap in eV
        """
        (cbm, vbm) = self.get_cbm_vbm(tol, abs_tol, spin)
        return max(cbm - vbm, 0.0)

    def __str__(self):
        """Returns a string which can be easily plotted (using gnuplot)."""
        if Spin.down in self.densities:
            stringarray = [f"#{'Energy':30s} {'DensityUp':30s} {'DensityDown':30s}"]
            for i, energy in enumerate(self.energies):
                stringarray.append(f"{energy:.5f} {self.densities[Spin.up][i]:.5f} {self.densities[Spin.down][i]:.5f}")
        else:
            stringarray = [f"#{'Energy':30s} {'DensityUp':30s}"]
            for i, energy in enumerate(self.energies):
                stringarray.append(f"{energy:.5f} {self.densities[Spin.up][i]:.5f}")
        return "\n".join(stringarray)

    @classmethod
    def from_dict(cls, d) -> Dos:
        """Returns Dos object from dict representation of Dos."""
        return Dos(
            d["efermi"],
            d["energies"],
            {Spin(int(k)): v for k, v in d["densities"].items()},
        )

    def as_dict(self) -> dict:
        """JSON-serializable dict representation of Dos."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "efermi": self.efermi,
            "energies": self.energies.tolist(),
            "densities": {str(spin): dens.tolist() for spin, dens in self.densities.items()},
        }


class FermiDos(Dos, MSONable):
    """This wrapper class helps relate the density of states, doping levels
    (i.e. carrier concentrations) and corresponding fermi levels. A negative
    doping concentration indicates the majority carriers are electrons
    (n-type doping); a positive doping concentration indicates holes are the
    majority carriers (p-type doping).
    """

    def __init__(
        self,
        dos: Dos,
        structure: Structure | None = None,
        nelecs: float | None = None,
        bandgap: float | None = None,
    ):
        """
        Args:
            dos: Pymatgen Dos object.
            structure: A structure. If not provided, the structure
                of the dos object will be used. If the dos does not have an
                associated structure object, an error will be thrown.
            nelecs: The number of electrons included in the energy range of
                dos. It is used for normalizing the densities. Default is the total
                number of electrons in the structure.
            bandgap: If set, the energy values are scissored so that the electronic
                band gap matches this value.
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

        # normalize total density of states based on integral at 0K
        tdos = np.array(self.get_densities())
        self.tdos = tdos * self.nelecs / (tdos * self.de)[self.energies <= self.efermi].sum()

        ecbm, evbm = self.get_cbm_vbm()
        self.idx_vbm = int(np.argmin(abs(self.energies - evbm)))
        self.idx_cbm = int(np.argmin(abs(self.energies - ecbm)))
        self.A_to_cm = 1e-8

        if bandgap:
            eref = self.efermi if evbm < self.efermi < ecbm else (evbm + ecbm) / 2.0

            idx_fermi = int(np.argmin(abs(self.energies - eref)))

            if idx_fermi == self.idx_vbm:
                # Fermi level and vbm should be different indices
                idx_fermi += 1

            self.energies[:idx_fermi] -= (bandgap - (ecbm - evbm)) / 2.0
            self.energies[idx_fermi:] += (bandgap - (ecbm - evbm)) / 2.0

    def get_doping(self, fermi_level: float, temperature: float) -> float:
        """Calculate the doping (majority carrier concentration) at a given
        Fermi level  and temperature. A simple Left Riemann sum is used for
        integrating the density of states over energy & equilibrium Fermi-Dirac
        distribution.

        Args:
            fermi_level: The fermi_level level in eV.
            temperature: The temperature in Kelvin.

        Returns:
            The doping concentration in units of 1/cm^3. Negative values
            indicate that the majority carriers are electrons (n-type doping)
            whereas positive values indicates the majority carriers are holes
            (p-type doping).
        """
        cb_integral = np.sum(
            self.tdos[self.idx_cbm :]
            * f0(self.energies[self.idx_cbm :], fermi_level, temperature)
            * self.de[self.idx_cbm :],
            axis=0,
        )
        vb_integral = np.sum(
            self.tdos[: self.idx_vbm + 1]
            * f0(-self.energies[: self.idx_vbm + 1], -fermi_level, temperature)
            * self.de[: self.idx_vbm + 1],
            axis=0,
        )
        return (vb_integral - cb_integral) / (self.volume * self.A_to_cm**3)

    def get_fermi_interextrapolated(
        self, concentration: float, temperature: float, warn: bool = True, c_ref: float = 1e10, **kwargs
    ) -> float:
        """Similar to get_fermi except that when get_fermi fails to converge,
        an interpolated or extrapolated fermi is returned with the assumption
        that the Fermi level changes linearly with log(abs(concentration)).

        Args:
            concentration: The doping concentration in 1/cm^3. Negative values
                represent n-type doping and positive values represent p-type
                doping.
            temperature: The temperature in Kelvin.
            warn: Whether to give a warning the first time the fermi cannot be
                found.
            c_ref: A doping concentration where get_fermi returns a
                value without error for both c_ref and -c_ref.
            **kwargs: Keyword arguments passed to the get_fermi function.

        Returns:
            The Fermi level. Note, the value is possibly interpolated or
            extrapolated and must be used with caution.
        """
        try:
            return self.get_fermi(concentration, temperature, **kwargs)
        except ValueError as e:
            if warn:
                warnings.warn(str(e))

            if abs(concentration) < c_ref:
                if abs(concentration) < 1e-10:
                    concentration = 1e-10

                # max(10, ) is to avoid log(0<x<1) and log(1+x) both of which
                # are slow
                f2 = self.get_fermi_interextrapolated(
                    max(10, abs(concentration) * 10.0), temperature, warn=False, **kwargs
                )
                f1 = self.get_fermi_interextrapolated(
                    -max(10, abs(concentration) * 10.0), temperature, warn=False, **kwargs
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

    def get_fermi(
        self,
        concentration: float,
        temperature: float,
        rtol: float = 0.01,
        nstep: int = 50,
        step: float = 0.1,
        precision: int = 8,
    ) -> float:
        """Finds the Fermi level at which the doping concentration at the given
        temperature (T) is equal to concentration. A greedy algorithm is used
        where the relative error is minimized by calculating the doping at a
        grid which continually becomes finer.

        Args:
            concentration: The doping concentration in 1/cm^3. Negative values
                represent n-type doping and positive values represent p-type
                doping.
            temperature: The temperature in Kelvin.
            rtol: The maximum acceptable relative error.
            nstep: The number of steps checked around a given Fermi level.
            step: Initial step in energy when searching for the Fermi level.
            precision: Essentially the decimal places of calculated Fermi level.

        Raises:
            ValueError: If the Fermi level cannot be found.

        Returns:
            The Fermi level in eV. Note that this is different from the default
            dos.efermi.
        """
        fermi = self.efermi  # initialize target fermi
        relative_error = [float("inf")]
        for _ in range(precision):
            f_range = np.arange(-nstep, nstep + 1) * step + fermi
            calc_doping = np.array([self.get_doping(f, temperature) for f in f_range])
            relative_error = np.abs(calc_doping / concentration - 1.0)  # type: ignore
            fermi = f_range[np.argmin(relative_error)]
            step /= 10.0

        if min(relative_error) > rtol:
            raise ValueError(f"Could not find fermi within {rtol:.1%} of {concentration=}")
        return fermi

    @classmethod
    def from_dict(cls, d) -> FermiDos:
        """Returns Dos object from dict representation of Dos."""
        dos = Dos(
            d["efermi"],
            d["energies"],
            {Spin(int(k)): v for k, v in d["densities"].items()},
        )
        return FermiDos(dos, structure=Structure.from_dict(d["structure"]), nelecs=d["nelecs"])

    def as_dict(self) -> dict:
        """JSON-serializable dict representation of Dos."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "efermi": self.efermi,
            "energies": self.energies.tolist(),
            "densities": {str(spin): dens.tolist() for spin, dens in self.densities.items()},
            "structure": self.structure,
            "nelecs": self.nelecs,
        }


class CompleteDos(Dos):
    """This wrapper class defines a total dos, and also provides a list of PDos.
    Mainly used by pymatgen.io.vasp.Vasprun to create a complete Dos from
    a vasprun.xml file. You are unlikely to try to generate this object
    manually.

    Attributes:
        structure (Structure): Structure associated with the CompleteDos.
        pdos (dict): Dict of partial densities of the form {Site:{Orbital:{Spin:Densities}}}.
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
            structure: Structure associated with this particular DOS.
            total_dos: total Dos for structure
            pdoss: The pdoss are supplied as an {Site: {Orbital: {Spin:Densities}}}
            normalize: Whether to normalize the densities by the volume of the structure.
                If True, the units of the densities are states/eV/Angstrom^3. Otherwise,
                the units are states/eV.
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

    def get_normalized(self) -> CompleteDos:
        """Returns a normalized version of the CompleteDos."""
        if self.norm_vol is not None:
            return self
        return CompleteDos(
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
            Dos containing densities for orbital of site.
        """
        return Dos(self.efermi, self.energies, self.pdos[site][orbital])

    def get_site_dos(self, site: PeriodicSite) -> Dos:
        """Get the total Dos for a site (all orbitals).

        Args:
            site: Site in Structure associated with CompleteDos.

        Returns:
            Dos containing summed orbital densities for site.
        """
        site_dos = functools.reduce(add_densities, self.pdos[site].values())
        return Dos(self.efermi, self.energies, site_dos)

    def get_site_spd_dos(self, site: PeriodicSite) -> dict[OrbitalType, Dos]:
        """Get orbital projected Dos of a particular site.

        Args:
            site: Site in Structure associated with CompleteDos.

        Returns:
            dict of {OrbitalType: Dos}, e.g. {OrbitalType.s: Dos object, ...}
        """
        spd_dos: dict[OrbitalType, dict[Spin, np.ndarray]] = {}
        for orb, pdos in self.pdos[site].items():
            orbital_type = _get_orb_type(orb)
            if orbital_type in spd_dos:
                spd_dos[orbital_type] = add_densities(spd_dos[orbital_type], pdos)
            else:
                spd_dos[orbital_type] = pdos  # type: ignore[assignment]
        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in spd_dos.items()}

    def get_site_t2g_eg_resolved_dos(self, site: PeriodicSite) -> dict[str, Dos]:
        """Get the t2g, eg projected DOS for a particular site.

        Args:
            site: Site in Structure associated with CompleteDos.

        Returns:
            dict[str, Dos]: A dict {"e_g": Dos, "t2g": Dos} containing summed e_g and t2g DOS for the site.
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
        """Get orbital projected Dos.

        Returns:
            dict of {OrbitalType: Dos}, e.g. {OrbitalType.s: Dos object, ...}
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
        """Get element projected Dos.

        Returns:
            dict of {Element: Dos}
        """
        el_dos = {}
        for site, atom_dos in self.pdos.items():
            el = site.specie
            for pdos in atom_dos.values():
                if el not in el_dos:
                    el_dos[el] = pdos
                else:
                    el_dos[el] = add_densities(el_dos[el], pdos)
        return {el: Dos(self.efermi, self.energies, densities) for el, densities in el_dos.items()}

    def get_element_spd_dos(self, el: SpeciesLike) -> dict[OrbitalType, Dos]:
        """Get element and spd projected Dos.

        Args:
            el: Element in Structure.composition associated with CompleteDos

        Returns:
            dict of {OrbitalType: Dos}, e.g. {OrbitalType.s: Dos object, ...}
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
        """Calculates spin polarization at Fermi level. If the
        calculation is not spin-polarized, None will be
        returned.

        See Sanvito et al., doi: 10.1126/sciadv.1602241 for
        an example usage.

        :return (float): spin polarization in range [0, 1],
        will also return NaN if spin polarization ill-defined
        (e.g. for insulator)
        """
        n_F = self.get_interpolated_value(self.efermi)

        n_F_up = n_F[Spin.up]
        if Spin.down not in n_F:
            return None
        n_F_down = n_F[Spin.down]

        if (n_F_up + n_F_down) == 0:
            # only well defined for metals or half-metals
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

        Args:
            band: Orbital type to get the band center of (default is d-band)
            elements: Elements to get the band center of (cannot be used in conjunction with site)
            sites: Sites to get the band center of (cannot be used in conjunction with el)
            spin: Spin channel to use. By default, the spin channels will be combined.

        Returns:
            band filling in eV, often denoted f_d for the d-band
        """
        # Get the projected DOS
        if elements and sites:
            raise ValueError("Both element and site cannot be specified.")

        densities: dict[Spin, ArrayLike] = {}
        if elements:
            for idx, el in enumerate(elements):
                spd_dos = self.get_element_spd_dos(el)[band]
                densities = (
                    spd_dos.densities if idx == 0 else add_densities(densities, spd_dos.densities)  # type: ignore
                )
            dos = Dos(self.efermi, self.energies, densities)
        elif sites:
            for idx, site in enumerate(sites):
                spd_dos = self.get_site_spd_dos(site)[band]
                densities = (
                    spd_dos.densities if idx == 0 else add_densities(densities, spd_dos.densities)  # type: ignore
                )
            dos = Dos(self.efermi, self.energies, densities)
        else:
            dos = self.get_spd_dos()[band]

        energies = dos.energies - dos.efermi
        dos_densities = dos.get_densities(spin=spin)

        # Only consider up to Fermi level in numerator
        energies = dos.energies - dos.efermi
        return np.trapz(dos_densities[energies < 0], x=energies[energies < 0]) / np.trapz(dos_densities, x=energies)

    def get_band_center(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: list[float] | None = None,
    ) -> float:
        """Compute the orbital-projected band center, defined as the first moment
        relative to the Fermi level
            int_{-inf}^{+inf} rho(E)*E dE/int_{-inf}^{+inf} rho(E) dE
        based on the work of Hammer and Norskov, Surf. Sci., 343 (1995) where the
        limits of the integration can be modified by erange and E is the set
        of energies taken with respect to the Fermi level. Note that the band center
        is often highly sensitive to the selected erange.

        Args:
            band: Orbital type to get the band center of (default is d-band)
            elements: Elements to get the band center of (cannot be used in conjunction with site)
            sites: Sites to get the band center of (cannot be used in conjunction with el)
            spin: Spin channel to use. By default, the spin channels will be combined.
            erange: [min, max] energy range to consider, with respect to the Fermi level.
                Default is None, which means all energies are considered.

        Returns:
            band center in eV, often denoted epsilon_d for the d-band center
        """
        return self.get_n_moment(1, elements=elements, sites=sites, band=band, spin=spin, erange=erange, center=False)

    def get_band_width(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: list[float] | None = None,
    ) -> float:
        """Get the orbital-projected band width, defined as the square root of the second moment
            sqrt(int_{-inf}^{+inf} rho(E)*(E-E_center)^2 dE/int_{-inf}^{+inf} rho(E) dE)
        where E_center is the orbital-projected band center, the limits of the integration can be
        modified by erange, and E is the set of energies taken with respect to the Fermi level.
        Note that the band width is often highly sensitive to the selected erange.

        Args:
            band: Orbital type to get the band center of (default is d-band)
            elements: Elements to get the band center of (cannot be used in conjunction with site)
            sites: Sites to get the band center of (cannot be used in conjunction with el)
            spin: Spin channel to use. By default, the spin channels will be combined.
            erange: [min, max] energy range to consider, with respect to the Fermi level.
                Default is None, which means all energies are considered.

        Returns:
            Orbital-projected band width in eV
        """
        return np.sqrt(self.get_n_moment(2, elements=elements, sites=sites, band=band, spin=spin, erange=erange))

    def get_band_skewness(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: list[float] | None = None,
    ) -> float:
        """Get the orbital-projected skewness, defined as the third standardized moment
            int_{-inf}^{+inf} rho(E)*(E-E_center)^3 dE/int_{-inf}^{+inf} rho(E) dE)
            /
            (int_{-inf}^{+inf} rho(E)*(E-E_center)^2 dE/int_{-inf}^{+inf} rho(E) dE))^(3/2)
        where E_center is the orbital-projected band center, the limits of the integration can be
        modified by erange, and E is the set of energies taken with respect to the Fermi level.
        Note that the skewness is often highly sensitive to the selected erange.

        Args:
            band: Orbitals to get the band center of (default is d-band)
            elements: Elements to get the band center of (cannot be used in conjunction with site)
            sites: Sites to get the band center of (cannot be used in conjunction with el)
            spin: Spin channel to use. By default, the spin channels will be combined.
            erange: [min, max] energy range to consider, with respect to the Fermi level.
                Default is None, which means all energies are considered.

        Returns:
            Orbital-projected skewness in eV
        """
        return self.get_n_moment(
            3, elements=elements, sites=sites, band=band, spin=spin, erange=erange
        ) / self.get_n_moment(2, elements=elements, sites=sites, band=band, spin=spin, erange=erange) ** (3 / 2)

    def get_band_kurtosis(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: list[float] | None = None,
    ) -> float:
        """Get the orbital-projected kurtosis, defined as the fourth standardized moment
            int_{-inf}^{+inf} rho(E)*(E-E_center)^4 dE/int_{-inf}^{+inf} rho(E) dE)
            /
            (int_{-inf}^{+inf} rho(E)*(E-E_center)^2 dE/int_{-inf}^{+inf} rho(E) dE))^2
        where E_center is the orbital-projected band center, the limits of the integration can be
        modified by erange, and E is the set of energies taken with respect to the Fermi level.
        Note that the skewness is often highly sensitive to the selected erange.

        Args:
            band: Orbital type to get the band center of (default is d-band)
            elements: Elements to get the band center of (cannot be used in conjunction with site)
            sites: Sites to get the band center of (cannot be used in conjunction with el)
            spin: Spin channel to use. By default, the spin channels will be combined.
            erange: [min, max] energy range to consider, with respect to the Fermi level.
                Default is None, which means all energies are considered.

        Returns:
            Orbital-projected kurtosis in eV
        """
        return (
            self.get_n_moment(4, elements=elements, sites=sites, band=band, spin=spin, erange=erange)
            / self.get_n_moment(2, elements=elements, sites=sites, band=band, spin=spin, erange=erange) ** 2
        )

    def get_n_moment(
        self,
        n: int,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
        spin: Spin | None = None,
        erange: list[float] | None = None,
        center: bool = True,
    ) -> float:
        """Get the nth moment of the DOS centered around the orbital-projected band center, defined as
            int_{-inf}^{+inf} rho(E)*(E-E_center)^n dE/int_{-inf}^{+inf} rho(E) dE
        where n is the order, E_center is the orbital-projected band center, the limits of the integration can be
        modified by erange, and E is the set of energies taken with respect to the Fermi level. If center is False,
        then the E_center reference is not used.

        Args:
            n: The order for the moment
            band: Orbital type to get the band center of (default is d-band)
            elements: Elements to get the band center of (cannot be used in conjunction with site)
            sites: Sites to get the band center of (cannot be used in conjunction with el)
            spin: Spin channel to use. By default, the spin channels will be combined.
            erange: [min, max] energy range to consider, with respect to the Fermi level.
                Default is None, which means all energies are considered.
            center: Take moments with respect to the band center

        Returns:
            Orbital-projected nth moment in eV
        """
        # Get the projected DOS
        if elements and sites:
            raise ValueError("Both element and site cannot be specified.")

        densities: Mapping[Spin, ArrayLike] = {}
        if elements:
            for i, el in enumerate(elements):
                spd_dos = self.get_element_spd_dos(el)[band]
                densities = spd_dos.densities if i == 0 else add_densities(densities, spd_dos.densities)
            dos = Dos(self.efermi, self.energies, densities)
        elif sites:
            for i, site in enumerate(sites):
                spd_dos = self.get_site_spd_dos(site)[band]
                densities = spd_dos.densities if i == 0 else add_densities(densities, spd_dos.densities)
            dos = Dos(self.efermi, self.energies, densities)
        else:
            dos = self.get_spd_dos()[band]

        energies = dos.energies - dos.efermi
        dos_densities = dos.get_densities(spin=spin)

        # Only consider a given erange, if desired
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
        return np.trapz(p**n * dos_densities, x=energies) / np.trapz(dos_densities, x=energies)

    def get_hilbert_transform(
        self,
        band: OrbitalType = OrbitalType.d,
        elements: list[SpeciesLike] | None = None,
        sites: list[PeriodicSite] | None = None,
    ) -> Dos:
        """Return the Hilbert transform of the orbital-projected density of states,
        often plotted for a Newns-Anderson analysis.

        Args:
            elements: Elements to get the band center of (cannot be used in conjunction with site)
            sites: Sites to get the band center of (cannot be used in conjunction with el)
            band: Orbitals to get the band center of (default is d-band)

        Returns:
            Hilbert transformation of the projected DOS.
        """
        # Get the projected DOS
        if elements and sites:
            raise ValueError("Both element and site cannot be specified.")

        densities: Mapping[Spin, ArrayLike] = {}
        if elements:
            for i, el in enumerate(elements):
                spd_dos = self.get_element_spd_dos(el)[band]
                densities = spd_dos.densities if i == 0 else add_densities(densities, spd_dos.densities)
            dos = Dos(self.efermi, self.energies, densities)
        elif sites:
            for i, site in enumerate(sites):
                spd_dos = self.get_site_spd_dos(site)[band]
                densities = spd_dos.densities if i == 0 else add_densities(densities, spd_dos.densities)
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
        erange: list[float] | None = None,
    ) -> float:
        """Get the orbital-projected upper band edge. The definition by Xin et al.
        Phys. Rev. B, 89, 115114 (2014) is used, which is the highest peak position of the
        Hilbert transform of the orbital-projected DOS.

        Args:
            band: Orbital type to get the band center of (default is d-band)
            elements: Elements to get the band center of (cannot be used in conjunction with site)
            sites: Sites to get the band center of (cannot be used in conjunction with el)
            spin: Spin channel to use. By default, the spin channels will be combined.
            erange: [min, max] energy range to consider, with respect to the Fermi level.
                Default is None, which means all energies are considered.

        Returns:
            Upper band edge in eV, often denoted epsilon_u
        """
        # Get the Hilbert-transformed DOS
        transformed_dos = self.get_hilbert_transform(elements=elements, sites=sites, band=band)

        energies = transformed_dos.energies - transformed_dos.efermi
        densities = transformed_dos.get_densities(spin=spin)

        # Only consider a given erange, if specified
        if erange:
            densities = densities[(energies >= erange[0]) & (energies <= erange[1])]
            energies = energies[(energies >= erange[0]) & (energies <= erange[1])]

        # Calculate the upper band edge
        return energies[np.argmax(densities)]

    def get_dos_fp(
        self,
        type: str = "summed_pdos",
        binning: bool = True,
        min_e: float | None = None,
        max_e: float | None = None,
        n_bins: int = 256,
        normalize: bool = True,
    ) -> NamedTuple:
        """Generates the DOS fingerprint based on work of
        F. Knoop, T. A. r Purcell, M. Scheffler, C. Carbogno, J. Open Source Softw. 2020, 5, 2671.
        Source - https://gitlab.com/vibes-developers/vibes/-/tree/master/vibes/materials_fp
        Copyright (c) 2020 Florian Knoop, Thomas A.R.Purcell, Matthias Scheffler, Christian Carbogno.

        Args:
            type (str): Specify fingerprint type needed can accept '{s/p/d/f/}summed_{pdos/tdos}'
            (default is summed_pdos)
            binning (bool): If true, the DOS fingerprint is binned using np.linspace and n_bins.
                Default is True.
            min_e (float): The minimum mode energy to include in the fingerprint (default is None)
            max_e (float): The maximum mode energy to include in the fingerprint (default is None)
            n_bins (int): Number of bins to be used in the fingerprint (default is 256)
            normalize (bool): If true, normalizes the area under fp to equal to 1. Default is True.

        Raises:
            ValueError: If type is not one of the accepted values {s/p/d/f/}summed_{pdos/tdos}.

        Returns:
            Fingerprint(namedtuple) : The electronic density of states fingerprint
            of format (energies, densities, type, n_bins)
        """
        fingerprint = namedtuple("fingerprint", "energies densities type n_bins bin_width")
        energies = self.energies - self.efermi

        if max_e is None:
            max_e = np.max(energies)

        if min_e is None:
            min_e = np.min(energies)

        pdos_obj = self.get_spd_dos()

        pdos = {}
        for key in pdos_obj:
            dens = pdos_obj[key].get_densities()

            pdos[key.name] = dens

        pdos["summed_pdos"] = np.sum(list(pdos.values()), axis=0)
        pdos["tdos"] = self.get_densities()

        try:
            densities = pdos[type]
            if len(energies) < n_bins:
                inds = np.where((energies >= min_e) & (energies <= max_e))
                return fingerprint(energies[inds], densities[inds], type, len(energies), np.diff(energies)[0])

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

            for ii, e1, e2 in zip(range(len(ener)), ener_bounds[0:-1], ener_bounds[1:]):
                inds = np.where((energies >= e1) & (energies < e2))
                dos_rebin[ii] = np.sum(densities[inds])
            if normalize:  # scale DOS bins to make area under histogram equal 1
                area = np.sum(dos_rebin * bin_width)
                dos_rebin_sc = dos_rebin / area
            else:
                dos_rebin_sc = dos_rebin

            return fingerprint(np.array([ener]), dos_rebin_sc, type, n_bins, bin_width)

        except KeyError:
            raise ValueError(
                "Please recheck type requested, either the orbital projections unavailable in input DOS or "
                "there's a typo in type."
            )

    @staticmethod
    def fp_to_dict(fp: NamedTuple) -> dict:
        """Converts a fingerprint into a dictionary.

        Args:
            fp: The DOS fingerprint to be converted into a dictionary

        Returns:
            dict: A dict of the fingerprint Keys=type, Values=np.ndarray(energies, densities)
        """
        fp_dict = {}
        fp_dict[fp[2]] = np.array([fp[0], fp[1]], dtype="object").T

        return fp_dict

    @staticmethod
    def get_dos_fp_similarity(
        fp1: NamedTuple,
        fp2: NamedTuple,
        col: int = 1,
        pt: int | str = "All",
        normalize: bool = False,
        tanimoto: bool = False,
    ) -> float:
        """Calculates the similarity index (dot product) of two fingerprints.

        Args:
            fp1 (NamedTuple): The 1st dos fingerprint object
            fp2 (NamedTuple): The 2nd dos fingerprint object
            col (int): The item in the fingerprints (0:energies,1: densities) to take the dot product of (default is 1)
            pt (int or str) : The index of the point that the dot product is to be taken (default is All)
            normalize (bool): If True normalize the scalar product to 1 (default is False)
            tanimoto (bool): If True will compute Tanimoto index (default is False)

        Raises:
            ValueError: If both tanimoto and normalize are set to True.

        Returns:
            float: Similarity index given by the dot product
        """
        fp1_dict = CompleteDos.fp_to_dict(fp1) if not isinstance(fp1, dict) else fp1

        fp2_dict = CompleteDos.fp_to_dict(fp2) if not isinstance(fp2, dict) else fp2

        if pt == "All":
            vec1 = np.array([pt[col] for pt in fp1_dict.values()]).flatten()
            vec2 = np.array([pt[col] for pt in fp2_dict.values()]).flatten()
        else:
            vec1 = fp1_dict[fp1[2][pt]][col]
            vec2 = fp2_dict[fp2[2][pt]][col]

        if not normalize and tanimoto:
            rescale = np.linalg.norm(vec1) ** 2 + np.linalg.norm(vec2) ** 2 - np.dot(vec1, vec2)
            return np.dot(vec1, vec2) / rescale

        if not tanimoto and normalize:
            rescale = np.linalg.norm(vec1) * np.linalg.norm(vec2)
            return np.dot(vec1, vec2) / rescale

        if not tanimoto and not normalize:
            rescale = 1.0
            return np.dot(vec1, vec2) / rescale

        raise ValueError(
            "Cannot compute similarity index. Please set either normalize=True or tanimoto=True or both to False."
        )

    @classmethod
    def from_dict(cls, d) -> CompleteDos:
        """Returns CompleteDos object from dict representation."""
        tdos = Dos.from_dict(d)
        struct = Structure.from_dict(d["structure"])
        pdoss = {}
        for i in range(len(d["pdos"])):
            at = struct[i]
            orb_dos = {}
            for orb_str, odos in d["pdos"][i].items():
                orb = Orbital[orb_str]
                orb_dos[orb] = {Spin(int(k)): v for k, v in odos["densities"].items()}
            pdoss[at] = orb_dos
        return CompleteDos(struct, tdos, pdoss)

    def as_dict(self) -> dict:
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
                    dd[str(orb)] = {
                        "densities": {str(int(spin)): list(dens) for spin, dens in pdos.items()}  # type: ignore
                    }
                dct["pdos"].append(dd)
            dct["atom_dos"] = {str(at): dos.as_dict() for at, dos in self.get_element_dos().items()}
            dct["spd_dos"] = {str(orb): dos.as_dict() for orb, dos in self.get_spd_dos().items()}
        return dct

    def __str__(self):
        return f"Complete DOS for {self.structure}"


class LobsterCompleteDos(CompleteDos):
    """Extended CompleteDOS for Lobster."""

    def get_site_orbital_dos(self, site: PeriodicSite, orbital: str) -> Dos:  # type: ignore
        """Get the Dos for a particular orbital of a particular site.

        Args:
            site: Site in Structure associated with CompleteDos.
            orbital: principal quantum number and orbital in string format, e.g. "4s".
                    possible orbitals are: "s", "p_y", "p_z", "p_x", "d_xy", "d_yz", "d_z^2",
                    "d_xz", "d_x^2-y^2", "f_y(3x^2-y^2)", "f_xyz",
                    "f_yz^2", "f_z^3", "f_xz^2", "f_z(x^2-y^2)", "f_x(x^2-3y^2)"
                    In contrast to the Cohpcar and the Cohplist objects, the strings from the Lobster files are used

        Returns:
            Dos containing densities of an orbital of a specific site.
        """
        if orbital[1:] not in [
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
        ]:
            raise ValueError("orbital is not correct")
        return Dos(self.efermi, self.energies, self.pdos[site][orbital])  # type: ignore

    def get_site_t2g_eg_resolved_dos(self, site: PeriodicSite) -> dict[str, Dos]:
        """Get the t2g, eg projected DOS for a particular site.

        Args:
            site: Site in Structure associated with CompleteDos.

        Returns:
            A dict {"e_g": Dos, "t2g": Dos} containing summed e_g and t2g DOS
            for the site.
        """
        warnings.warn("Are the orbitals correctly oriented? Are you sure?")
        t2g_dos = []
        eg_dos = []
        for s, atom_dos in self.pdos.items():
            if s == site:
                for orb, pdos in atom_dos.items():
                    if _get_orb_lobster(orb) in (Orbital.dxy, Orbital.dxz, Orbital.dyz):
                        t2g_dos.append(pdos)
                    elif _get_orb_lobster(orb) in (Orbital.dx2, Orbital.dz2):
                        eg_dos.append(pdos)
        return {
            "t2g": Dos(self.efermi, self.energies, functools.reduce(add_densities, t2g_dos)),
            "e_g": Dos(self.efermi, self.energies, functools.reduce(add_densities, eg_dos)),
        }

    def get_spd_dos(self) -> dict[str, Dos]:  # type: ignore
        """Get orbital projected Dos.
        For example, if 3s and 4s are included in the basis of some element, they will be both summed in the orbital
        projected DOS.

        Returns:
            dict of {orbital: Dos}, e.g. {"s": Dos object, ...}
        """
        spd_dos = {}
        for atom_dos in self.pdos.values():
            for orb, pdos in atom_dos.items():
                orbital_type = _get_orb_type_lobster(orb)
                if orbital_type not in spd_dos:
                    spd_dos[orbital_type] = pdos
                else:
                    spd_dos[orbital_type] = add_densities(spd_dos[orbital_type], pdos)

        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in spd_dos.items()}  # type: ignore

    def get_element_spd_dos(self, el: SpeciesLike) -> dict[str, Dos]:  # type: ignore
        """Get element and spd projected Dos.

        Args:
            el: Element in Structure.composition associated with LobsterCompleteDos

        Returns:
            dict of {OrbitalType.s: densities, OrbitalType.p: densities, OrbitalType.d: densities}
        """
        el = get_el_sp(el)
        el_dos = {}
        for site, atom_dos in self.pdos.items():
            if site.specie == el:
                for orb, pdos in atom_dos.items():
                    orbital_type = _get_orb_type_lobster(orb)
                    if orbital_type not in el_dos:
                        el_dos[orbital_type] = pdos
                    else:
                        el_dos[orbital_type] = add_densities(el_dos[orbital_type], pdos)

        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in el_dos.items()}  # type: ignore

    @classmethod
    def from_dict(cls, d) -> LobsterCompleteDos:
        """Hydrate CompleteDos object from dict representation."""
        tdos = Dos.from_dict(d)
        struct = Structure.from_dict(d["structure"])
        pdoss = {}
        for i in range(len(d["pdos"])):
            at = struct[i]
            orb_dos = {}
            for orb_str, odos in d["pdos"][i].items():
                orb = orb_str
                orb_dos[orb] = {Spin(int(k)): v for k, v in odos["densities"].items()}
            pdoss[at] = orb_dos
        return LobsterCompleteDos(struct, tdos, pdoss)


def add_densities(density1: Mapping[Spin, ArrayLike], density2: Mapping[Spin, ArrayLike]) -> dict[Spin, np.ndarray]:
    """Sum two densities.

    Args:
        density1: First density.
        density2: Second density.

    Returns:
        dict[Spin, np.ndarray]
    """
    return {spin: np.array(density1[spin]) + np.array(density2[spin]) for spin in density1}


def _get_orb_type(orb) -> OrbitalType:
    try:
        return orb.orbital_type
    except AttributeError:
        return orb


def f0(E, fermi, T) -> float:
    """Return the equilibrium fermi-dirac.

    Args:
        E (float): energy in eV
        fermi (float): the Fermi level in eV
        T (float): the temperature in kelvin

    Returns:
        float
    """
    return 1.0 / (1.0 + np.exp((E - fermi) / (_cd("Boltzmann constant in eV/K") * T)))


def _get_orb_type_lobster(orb) -> OrbitalType | None:
    """
    Args:
        orb: string representation of orbital.

    Returns:
        OrbitalType
    """
    orb_labs = ["s", "p_y", "p_z", "p_x", "d_xy", "d_yz", "d_z^2", "d_xz", "d_x^2-y^2"]
    orb_labs += ["f_y(3x^2-y^2)", "f_xyz", "f_yz^2", "f_z^3", "f_xz^2", "f_z(x^2-y^2)", "f_x(x^2-3y^2)"]

    try:
        orbital = Orbital(orb_labs.index(orb[1:]))
        return orbital.orbital_type
    except AttributeError:
        print("Orb not in list")
    return None


def _get_orb_lobster(orb):
    """
    Args:
        orb: string representation of orbital.

    Returns:
        Orbital.
    """
    orb_labs = [
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
    ]

    try:
        return Orbital(orb_labs.index(orb[1:]))
    except AttributeError:
        print("Orb not in list")
    return None
