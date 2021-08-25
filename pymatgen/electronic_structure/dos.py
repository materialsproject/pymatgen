# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines classes to represent the density of states, etc.
"""

import functools
import warnings
from typing import Dict

import numpy as np
from monty.json import MSONable
from scipy.constants.codata import value as _cd

from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.spectrum import Spectrum
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.util.coord import get_linear_interpolated_value
from pymatgen.util.typing import ArrayLike, SpeciesLike


class DOS(Spectrum):
    """
    Replacement basic DOS object. All other DOS objects are extended versions
    of this object. Work in progress.

    .. attribute: energies

        The sequence of energies

    .. attribute: densities

        A dict of spin densities, e.g., {Spin.up: [...], Spin.down: [...]}

    .. attribute: efermi

        Fermi level
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

    def get_interpolated_gap(self, tol: float = 0.001, abs_tol: bool = False, spin: Spin = None):
        """
        Expects a DOS object and finds the gap

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
            tol = tol * tdos.sum() / tdos.shape[0]  # type: ignore
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
        """
        Expects a DOS object and finds the cbm and vbm.

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
            tol = tol * tdos.sum() / tdos.shape[0]  # type: ignore

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

    def get_gap(self, tol: float = 0.001, abs_tol: bool = False, spin: Spin = None):
        """
        Expects a DOS object and finds the gap.

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
        """
        Returns a string which can be easily plotted (using gnuplot).
        """
        if Spin.down in self.densities:
            stringarray = ["#{:30s} {:30s} {:30s}".format("Energy", "DensityUp", "DensityDown")]
            for i, energy in enumerate(self.energies):
                stringarray.append(
                    "{:.5f} {:.5f} {:.5f}".format(energy, self.densities[Spin.up][i], self.densities[Spin.down][i])
                )
        else:
            stringarray = ["#{:30s} {:30s}".format("Energy", "DensityUp")]
            for i, energy in enumerate(self.energies):
                stringarray.append("{:.5f} {:.5f}".format(energy, self.densities[Spin.up][i]))
        return "\n".join(stringarray)


class Dos(MSONable):
    """
    Basic DOS object. All other DOS objects are extended versions of this
    object.

    .. attribute: energies

        The sequence of energies

    .. attribute: densities

        A dict of spin densities, e.g., {Spin.up: [...], Spin.down: [...]}

    .. attribute: efermi

        Fermi level
    """

    def __init__(self, efermi: float, energies: ArrayLike, densities: Dict[Spin, ArrayLike]):
        """
        Args:
            efermi: Fermi level energy
            energies: A sequences of energies
            densities ({Spin: np.array}): representing the density of states
                for each Spin.
        """
        self.efermi = efermi
        self.energies = np.array(energies)
        self.densities = {k: np.array(d) for k, d in densities.items()}

    def get_densities(self, spin: Spin = None):
        """
        Returns the density of states for a particular spin.

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
        """
        Returns the Dict representation of the densities, {Spin: densities},
        but with a Gaussian smearing of std dev sigma applied about the fermi
        level.

        Args:
            sigma: Std dev of Gaussian smearing function.

        Returns:
            Dict of Gaussian-smeared densities.
        """
        from scipy.ndimage.filters import gaussian_filter1d

        smeared_dens = {}
        diff = [self.energies[i + 1] - self.energies[i] for i in range(len(self.energies) - 1)]
        avgdiff = sum(diff) / len(diff)
        for spin, dens in self.densities.items():
            smeared_dens[spin] = gaussian_filter1d(dens, sigma / avgdiff)
        return smeared_dens

    def __add__(self, other):
        """
        Adds two DOS together. Checks that energy scales are the same.
        Otherwise, a ValueError is thrown.

        Args:
            other: Another DOS object.

        Returns:
            Sum of the two DOSs.
        """
        if not all(np.equal(self.energies, other.energies)):
            raise ValueError("Energies of both DOS are not compatible!")
        densities = {spin: self.densities[spin] + other.densities[spin] for spin in self.densities.keys()}
        return Dos(self.efermi, self.energies, densities)

    def get_interpolated_value(self, energy: float):
        """
        Returns interpolated density for a particular energy.

        Args:
            energy: Energy to return the density for.
        """
        f = {}
        for spin in self.densities.keys():
            f[spin] = get_linear_interpolated_value(self.energies, self.densities[spin], energy)
        return f

    def get_interpolated_gap(self, tol: float = 0.001, abs_tol: bool = False, spin: Spin = None):
        """
        Expects a DOS object and finds the gap

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

    def get_cbm_vbm(self, tol: float = 0.001, abs_tol: bool = False, spin: Spin = None):
        """
        Expects a DOS object and finds the cbm and vbm.

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

    def get_gap(self, tol: float = 0.001, abs_tol: bool = False, spin: Spin = None):
        """
        Expects a DOS object and finds the gap.

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
        """
        Returns a string which can be easily plotted (using gnuplot).
        """
        if Spin.down in self.densities:
            stringarray = ["#{:30s} {:30s} {:30s}".format("Energy", "DensityUp", "DensityDown")]
            for i, energy in enumerate(self.energies):
                stringarray.append(
                    "{:.5f} {:.5f} {:.5f}".format(energy, self.densities[Spin.up][i], self.densities[Spin.down][i])
                )
        else:
            stringarray = ["#{:30s} {:30s}".format("Energy", "DensityUp")]
            for i, energy in enumerate(self.energies):
                stringarray.append("{:.5f} {:.5f}".format(energy, self.densities[Spin.up][i]))
        return "\n".join(stringarray)

    @classmethod
    def from_dict(cls, d) -> "Dos":
        """
        Returns Dos object from dict representation of Dos.
        """
        return Dos(
            d["efermi"],
            d["energies"],
            {Spin(int(k)): v for k, v in d["densities"].items()},
        )

    def as_dict(self) -> dict:
        """
        Json-serializable dict representation of Dos.
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "efermi": self.efermi,
            "energies": self.energies.tolist(),
            "densities": {str(spin): dens.tolist() for spin, dens in self.densities.items()},
        }


class FermiDos(Dos, MSONable):
    """
    This wrapper class helps relate the density of states, doping levels
    (i.e. carrier concentrations) and corresponding fermi levels. A negative
    doping concentration indicates the majority carriers are electrons
    (n-type doping); a positive doping concentration indicates holes are the
    majority carriers (p-type doping).
    """

    def __init__(
        self,
        dos: Dos,
        structure: Structure = None,
        nelecs: float = None,
        bandgap: float = None,
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
                raise ValueError("Structure object is not provided and not " "present in dos")

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
            if evbm < self.efermi < ecbm:
                eref = self.efermi
            else:
                eref = (evbm + ecbm) / 2.0

            idx_fermi = int(np.argmin(abs(self.energies - eref)))

            if idx_fermi == self.idx_vbm:
                # Fermi level and vbm should be different indices
                idx_fermi += 1

            self.energies[:idx_fermi] -= (bandgap - (ecbm - evbm)) / 2.0
            self.energies[idx_fermi:] += (bandgap - (ecbm - evbm)) / 2.0

    def get_doping(self, fermi_level: float, temperature: float) -> float:
        """
        Calculate the doping (majority carrier concentration) at a given
        fermi level  and temperature. A simple Left Riemann sum is used for
        integrating the density of states over energy & equilibrium Fermi-Dirac
        distribution.

        Args:
            fermi_level: The fermi_level level in eV.
            temperature: The temperature in Kelvin.

        Returns:
            The doping concentration in units of 1/cm^3. Negative values
            indicate that the majority carriers are electrons (n-type doping)
            whereas positivie values indicates the majority carriers are holes
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
        return (vb_integral - cb_integral) / (self.volume * self.A_to_cm ** 3)

    def get_fermi_interextrapolated(
        self, concentration: float, temperature: float, warn: bool = True, c_ref: float = 1e10, **kwargs
    ) -> float:
        """
        Similar to get_fermi except that when get_fermi fails to converge,
        an interpolated or extrapolated fermi is returned with the assumption
        that the fermi level changes linearly with log(abs(concentration)).

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
            c_newlog = np.sign(concentration) * np.log(abs(self.get_doping(f_new, temperature)))
            slope = (f_new - f_ref) / (c_newlog - np.sign(concentration) * 10.0)
            return f_new + slope * (clog - c_newlog)

    def get_fermi(
        self,
        concentration: float,
        temperature: float,
        rtol: float = 0.01,
        nstep: int = 50,
        step: float = 0.1,
        precision: int = 8,
    ):
        """
        Finds the fermi level at which the doping concentration at the given
        temperature (T) is equal to concentration. A greedy algorithm is used
        where the relative error is minimized by calculating the doping at a
        grid which continually becomes finer.

        Args:
            concentration: The doping concentration in 1/cm^3. Negative values
                represent n-type doping and positive values represent p-type
                doping.
            temperature: The temperature in Kelvin.
            rtol: The maximum acceptable relative error.
            nstep: THe number of steps checked around a given fermi level.
            step: Initial step in energy when searching for the Fermi level.
            precision: Essentially the decimal places of calculated Fermi level.

        Returns:
            The fermi level in eV.. Note that this is different from the default
            dos.efermi.
        """
        fermi = self.efermi  # initialize target fermi
        relative_error = [float("inf")]
        for _ in range(precision):
            frange = np.arange(-nstep, nstep + 1) * step + fermi
            calc_doping = np.array([self.get_doping(f, temperature) for f in frange])
            relative_error = np.abs(calc_doping / concentration - 1.0)
            fermi = frange[np.argmin(relative_error)]
            step /= 10.0

        if min(relative_error) > rtol:
            raise ValueError("Could not find fermi within {}% of concentration={}".format(rtol * 100, concentration))
        return fermi

    @classmethod
    def from_dict(cls, d) -> "FermiDos":
        """
        Returns Dos object from dict representation of Dos.
        """
        dos = Dos(
            d["efermi"],
            d["energies"],
            {Spin(int(k)): v for k, v in d["densities"].items()},
        )
        return FermiDos(dos, structure=Structure.from_dict(d["structure"]), nelecs=d["nelecs"])

    def as_dict(self) -> dict:
        """
        Json-serializable dict representation of Dos.
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "efermi": self.efermi,
            "energies": self.energies.tolist(),
            "densities": {str(spin): dens.tolist() for spin, dens in self.densities.items()},
            "structure": self.structure,
            "nelecs": self.nelecs,
        }


class CompleteDos(Dos):
    """
    This wrapper class defines a total dos, and also provides a list of PDos.
    Mainly used by pymatgen.io.vasp.Vasprun to create a complete Dos from
    a vasprun.xml file. You are unlikely to try to generate this object
    manually.

    .. attribute:: structure

        Structure associated with the CompleteDos.

    .. attribute:: pdos

        Dict of partial densities of the form {Site:{Orbital:{Spin:Densities}}}
    """

    def __init__(
        self, structure: Structure, total_dos: Dos, pdoss: Dict[PeriodicSite, Dict[Orbital, Dict[Spin, ArrayLike]]]
    ):
        """
        Args:
            structure: Structure associated with this particular DOS.
            total_dos: total Dos for structure
            pdoss: The pdoss are supplied as an {Site:{Orbital:{
                Spin:Densities}}}
        """
        super().__init__(
            total_dos.efermi,
            energies=total_dos.energies,
            densities={k: np.array(d) for k, d in total_dos.densities.items()},
        )
        self.pdos = pdoss
        self.structure = structure

    def get_site_orbital_dos(self, site: PeriodicSite, orbital: Orbital) -> Dos:
        """
        Get the Dos for a particular orbital of a particular site.

        Args:
            site: Site in Structure associated with CompleteDos.
            orbital: Orbital in the site.

        Returns:
            Dos containing densities for orbital of site.
        """
        return Dos(self.efermi, self.energies, self.pdos[site][orbital])

    def get_site_dos(self, site: PeriodicSite) -> Dos:
        """
        Get the total Dos for a site (all orbitals).

        Args:
            site: Site in Structure associated with CompleteDos.

        Returns:
            Dos containing summed orbital densities for site.
        """
        site_dos = functools.reduce(add_densities, self.pdos[site].values())
        return Dos(self.efermi, self.energies, site_dos)

    def get_site_spd_dos(self, site: PeriodicSite) -> Dict[Orbital, Dos]:
        """
        Get orbital projected Dos of a particular site

        Args:
            site: Site in Structure associated with CompleteDos.

        Returns:
            dict of {orbital: Dos}, e.g. {"s": Dos object, ...}
        """
        spd_dos: Dict[Orbital, Dict[Spin, ArrayLike]] = dict()
        for orb, pdos in self.pdos[site].items():
            orbital_type = _get_orb_type(orb)
            if orbital_type in spd_dos:
                spd_dos[orbital_type] = add_densities(spd_dos[orbital_type], pdos)
            else:
                spd_dos[orbital_type] = pdos
        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in spd_dos.items()}

    def get_site_t2g_eg_resolved_dos(self, site: PeriodicSite) -> Dict[str, Dos]:
        """
        Get the t2g, eg projected DOS for a particular site.

        Args:
            site: Site in Structure associated with CompleteDos.

        Returns:
            A dict {"e_g": Dos, "t2g": Dos} containing summed e_g and t2g DOS
            for the site.
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

    def get_spd_dos(self) -> Dict[Orbital, Dos]:
        """
        Get orbital projected Dos.

        Returns:
            dict of {orbital: Dos}, e.g. {"s": Dos object, ...}
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

    def get_element_dos(self) -> Dict[SpeciesLike, Dos]:
        """
        Get element projected Dos.

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

    def get_element_spd_dos(self, el: SpeciesLike) -> Dict[Orbital, Dos]:
        """
        Get element and spd projected Dos

        Args:
            el: Element in Structure.composition associated with CompleteDos

        Returns:
            dict of {Element: {"S": densities, "P": densities, "D": densities}}
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
    def spin_polarization(self) -> float:
        """
        Calculates spin polarization at Fermi level.

        See Sanvito et al., doi: 10.1126/sciadv.1602241 for
        an example usage.

        :return (float): spin polarization in range [0, 1],
        will also return NaN if spin polarization ill-defined
        (e.g. for insulator)
        """
        n_F = self.get_interpolated_value(self.efermi)

        n_F_up = n_F[Spin.up]
        n_F_down = n_F[Spin.down]

        if (n_F_up + n_F_down) == 0:
            # only well defined for metals or half-mteals
            return float("NaN")

        spin_polarization = (n_F_up - n_F_down) / (n_F_up + n_F_down)

        return abs(spin_polarization)

    @classmethod
    def from_dict(cls, d) -> "CompleteDos":
        """
        Returns CompleteDos object from dict representation.
        """
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
        """
        Json-serializable dict representation of CompleteDos.
        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
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
                d["pdos"].append(dd)
            d["atom_dos"] = {str(at): dos.as_dict() for at, dos in self.get_element_dos().items()}
            d["spd_dos"] = {str(orb): dos.as_dict() for orb, dos in self.get_spd_dos().items()}
        return d

    def __str__(self):
        return "Complete DOS for " + str(self.structure)


class LobsterCompleteDos(CompleteDos):
    """
    Extended CompleteDOS for Lobster
    """

    def get_site_orbital_dos(self, site: PeriodicSite, orbital: str) -> Dos:  # type: ignore
        """
        Get the Dos for a particular orbital of a particular site.

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

    def get_site_t2g_eg_resolved_dos(self, site: PeriodicSite) -> Dict[str, Dos]:
        """
        Get the t2g, eg projected DOS for a particular site.
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

    def get_spd_dos(self) -> Dict[str, Dos]:  # type: ignore
        """
        Get orbital projected Dos.
        For example, if 3s and 4s are included in the basis of some element, they will be both summed in the orbital
        projected DOS

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

        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in spd_dos.items()}

    def get_element_spd_dos(self, el: SpeciesLike) -> Dict[str, Dos]:  # type: ignore
        """
        Get element and spd projected Dos


        Args:
            el: Element in Structure.composition associated with LobsterCompleteDos

        Returns:
            dict of {"S": densities, "P": densities, "D": densities}
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

        return {orb: Dos(self.efermi, self.energies, densities) for orb, densities in el_dos.items()}

    @classmethod
    def from_dict(cls, d) -> "LobsterCompleteDos":
        """
        Returns: CompleteDos object from dict representation.
        """
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


def add_densities(density1: Dict[Spin, ArrayLike], density2: Dict[Spin, ArrayLike]) -> Dict[Spin, ArrayLike]:
    """
    Method to sum two densities.

    Args:
        density1: First density.
        density2: Second density.

    Returns:
        Dict of {spin: density}.
    """
    return {spin: np.array(density1[spin]) + np.array(density2[spin]) for spin in density1.keys()}


def _get_orb_type(orb):
    try:
        return orb.orbital_type
    except AttributeError:
        return orb


def f0(E, fermi, T):
    """
    Returns the equilibrium fermi-dirac.
    Args:
        E (float): energy in eV
        fermi (float): the fermi level in eV
        T (float): the temperature in kelvin
    """
    return 1.0 / (1.0 + np.exp((E - fermi) / (_cd("Boltzmann constant in eV/K") * T)))


def _get_orb_type_lobster(orb):
    """
    Args:
     orb: string representation of orbital
    Returns:
     OrbitalType
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
        orbital = Orbital(orb_labs.index(orb[1:]))
        return orbital.orbital_type
    except AttributeError:
        print("Orb not in list")
    return None


def _get_orb_lobster(orb):
    """
    Args:
        orb: string representation of orbital
    Returns:
         Orbital
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
        orbital = Orbital(orb_labs.index(orb[1:]))
        return orbital
    except AttributeError:
        print("Orb not in list")
    return None
