"""This module defines classes to represent:
    - Crystal orbital Hamilton population (COHP) and integrated COHP (ICOHP).
    - Crystal orbital overlap population (COOP).
    - Crystal orbital bond index (COBI).

If you use this module, please cite:
J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier,
"Automated Bonding Analysis with Crystal Orbital Hamilton Populations",
ChemPlusChem 2022, e202200123,
DOI: 10.1002/cplu.202200123.
"""

from __future__ import annotations

import re
import sys
import warnings
from typing import TYPE_CHECKING

import numpy as np
from monty.json import MSONable
from scipy.interpolate import InterpolatedUnivariateSpline

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.io.lmto import LMTOCopl
from pymatgen.io.lobster import Cohpcar
from pymatgen.util.coord import get_linear_interpolated_value
from pymatgen.util.due import Doi, due
from pymatgen.util.num import round_to_sigfigs

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, Literal

    from numpy.typing import NDArray
    from typing_extensions import Self

    from pymatgen.util.typing import PathLike, SpinLike, Vector3D

__author__ = "Marco Esters, Janine George"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.2"
__maintainer__ = "Janine George"
__email__ = "janinegeorge.ulfen@gmail.com"
__date__ = "Dec 13, 2017"

due.cite(
    Doi("10.1002/cplu.202200123"),
    description="Automated Bonding Analysis with Crystal Orbital Hamilton Populations",
)


class Cohp(MSONable):
    """Basic COHP object."""

    def __init__(
        self,
        efermi: float,
        energies: Sequence[float],
        cohp: dict[Spin, NDArray],
        are_coops: bool = False,
        are_cobis: bool = False,
        are_multi_center_cobis: bool = False,
        icohp: dict[Spin, NDArray] | None = None,
    ) -> None:
        """
        Args:
            efermi (float): The Fermi level.
            energies (Sequence[float]): Energies.
            cohp ({Spin: NDArrary}): The COHP for each spin.
            are_coops (bool): Whether this object describes COOPs.
            are_cobis (bool): Whether this object describes COBIs.
            are_multi_center_cobis (bool): Whether this object describes multi-center COBIs.
            icohp ({Spin: NDArrary}): The ICOHP for each spin.
        """
        self.are_coops = are_coops
        self.are_cobis = are_cobis
        self.are_multi_center_cobis = are_multi_center_cobis
        self.efermi = efermi
        self.energies = np.array(energies)
        self.cohp = cohp
        self.icohp = icohp

    def __repr__(self) -> str:
        """A string that can be easily plotted (e.g. using gnuplot)."""
        if self.are_coops:
            cohp_str = "COOP"
        elif self.are_cobis or self.are_multi_center_cobis:
            cohp_str = "COBI"
        else:
            cohp_str = "COHP"

        header = ["Energy", f"{cohp_str}Up"]
        data = [self.energies, self.cohp[Spin.up]]
        if Spin.down in self.cohp:
            header.append(f"{cohp_str}Down")
            data.append(self.cohp[Spin.down])
        if self.icohp:
            header.append(f"I{cohp_str}Up")
            data.append(self.icohp[Spin.up])
            if Spin.down in self.cohp:
                header.append(f"I{cohp_str}Down")
                data.append(self.icohp[Spin.down])
        format_header = "#" + " ".join("{:15s}" for __ in header)
        format_data = " ".join("{:.5f}" for __ in header)
        str_arr = [format_header.format(*header)]
        for idx in range(len(self.energies)):
            str_arr.append(format_data.format(*(d[idx] for d in data)))
        return "\n".join(str_arr)

    def as_dict(self) -> dict[str, Any]:
        """JSON-serializable dict representation of COHP."""
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "are_coops": self.are_coops,
            "are_cobis": self.are_cobis,
            "are_multi_center_cobis": self.are_multi_center_cobis,
            "efermi": self.efermi,
            "energies": self.energies.tolist(),
            "COHP": {str(spin): pops.tolist() for spin, pops in self.cohp.items()},
        }
        if self.icohp:
            dct["ICOHP"] = {str(spin): pops.tolist() for spin, pops in self.icohp.items()}
        return dct

    def get_cohp(
        self,
        spin: SpinLike | None = None,
        integrated: bool = False,
    ) -> dict[Spin, NDArray] | None:
        """Get the COHP or ICOHP for a particular spin.

        Args:
            spin (SpinLike): Selected spin. If is None and both
                spins are present, both will be returned.
            integrated: Return ICOHP (True) or COHP (False).

        Returns:
            dict: The COHP or ICOHP for the selected spin.
        """
        populations = self.icohp if integrated else self.cohp

        if populations is None:
            return None
        if spin is None:
            return populations
        if isinstance(spin, int):
            spin = Spin(spin)
        elif isinstance(spin, str):
            spin = Spin({"up": 1, "down": -1}[spin.lower()])
        return {spin: populations[spin]}

    def get_icohp(
        self,
        spin: SpinLike | None = None,
    ) -> dict[Spin, NDArray] | None:
        """Convenient wrapper to get the ICOHP for a particular spin."""
        return self.get_cohp(spin=spin, integrated=True)

    def get_interpolated_value(
        self,
        energy: float,
        integrated: bool = False,
    ) -> dict[Spin, float]:
        """Get the interpolated COHP for a particular energy.

        Args:
            energy (float): Energy to get the COHP value for.
            integrated (bool): Return ICOHP (True) or COHP (False).
        """
        inters = {}
        for spin in self.cohp:
            if not integrated:
                inters[spin] = get_linear_interpolated_value(self.energies, self.cohp[spin], energy)
            elif self.icohp is not None:
                inters[spin] = get_linear_interpolated_value(self.energies, self.icohp[spin], energy)
            else:
                raise ValueError("ICOHP is empty.")
        return inters

    def has_antibnd_states_below_efermi(
        self,
        spin: SpinLike | None = None,
        limit: float = 0.01,
    ) -> dict[Spin, bool] | None:
        """Get dict of antibonding states below the Fermi level for the spin.

        Args:
            spin (SpinLike): Selected spin.
            limit (float): Only COHP higher than this value will be considered.
        """
        populations = self.cohp
        n_energies_below_efermi = sum(energy <= self.efermi for energy in self.energies)

        if populations is None:
            return None

        dict_to_return = {}
        if spin is None:
            for sp, cohp_vals in populations.items():
                # NOTE: Casting to bool is necessary, otherwise ended up
                # getting "bool_" instead of "bool" from NumPy
                dict_to_return[sp] = bool((max(cohp_vals[:n_energies_below_efermi])) > limit)

        else:
            if isinstance(spin, int):
                spin = Spin(spin)
            elif isinstance(spin, str):
                spin = Spin({"up": 1, "down": -1}[spin.lower()])
            dict_to_return[spin] = bool((max(populations[spin][:n_energies_below_efermi])) > limit)

        return dict_to_return

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """Generate Cohp from a dict representation."""
        icohp = {Spin(int(key)): np.array(val) for key, val in dct["ICOHP"].items()} if "ICOHP" in dct else None

        return cls(
            dct["efermi"],
            dct["energies"],
            {Spin(int(key)): np.array(val) for key, val in dct["COHP"].items()},
            icohp=icohp,
            are_coops=dct["are_coops"],
            are_cobis=dct.get("are_cobis", False),
            are_multi_center_cobis=dct.get("are_multi_center_cobis", False),
        )


class CompleteCohp(Cohp):
    """A wrapper that defines an average COHP, and individual COHPs.

    Attributes:
        are_coops (bool): Whether the object consists of COOPs.
        are_cobis (bool): Whether the object consists of COBIs.
        efermi (float): The Fermi level.
        energies (Sequence[float]): Sequence of energies.
        structure (Structure): Structure associated with the COHPs.
        cohp (Sequence[float]): The average COHP.
        icohp (Sequence[float]): The average ICOHP.
        all_cohps (dict[str, Sequence[float]]): COHPs for individual bonds of the form {label: COHP}.
        orb_res_cohp (dict[str, Dict[str, Sequence[float]]]): Orbital-resolved COHPs.
    """

    def __init__(
        self,
        structure: Structure,
        avg_cohp: Cohp,
        cohp_dict: dict[str, Cohp],
        bonds: dict[str, Any] | None = None,
        are_coops: bool = False,
        are_cobis: bool = False,
        are_multi_center_cobis: bool = False,
        orb_res_cohp: dict[str, dict] | None = None,
    ) -> None:
        """
        Args:
            structure (Structure): Structure associated with this COHP.
            avg_cohp (Cohp): The average COHP.
            cohp_dict (dict[str, Cohp]): COHP for individual bonds of the form
                {label: COHP}.
            bonds (dict[str, Any]): Information on the bonds of the form
                {label: {key: val}}. The value can be any information,
                but typically contains the sites, the bond length,
                and the number of bonds. If nothing is
                supplied, it will default to an empty dict.
            are_coops (bool): Whether the Cohp objects are COOPs.
                Defaults to False for COHPs.
            are_cobis (bool): Whether the Cohp objects are COBIs.
                Defaults to False for COHPs.
            are_multi_center_cobis (bool): Whether the Cohp objects are multi-center COBIs.
                Defaults to False for COHPs.
            orb_res_cohp (dict): Orbital-resolved COHPs.
        """
        if (
            (are_coops and are_cobis)
            or (are_coops and are_multi_center_cobis)
            or (are_cobis and are_multi_center_cobis)
        ):
            raise ValueError("You cannot have info about COOPs, COBIs and/or multi-center COBIS in the same file.")

        super().__init__(
            avg_cohp.efermi,
            avg_cohp.energies,
            avg_cohp.cohp,
            are_coops=are_coops,
            are_cobis=are_cobis,
            are_multi_center_cobis=are_multi_center_cobis,
            icohp=avg_cohp.icohp,
        )
        self.structure = structure
        self.are_coops = are_coops
        self.are_cobis = are_cobis
        self.are_multi_center_cobis = are_multi_center_cobis
        self.all_cohps = cohp_dict
        self.orb_res_cohp = orb_res_cohp
        self.bonds = bonds or {label: {} for label in self.all_cohps}

    def __str__(self) -> str:
        if self.are_coops:
            header = "COOPs"
        elif self.are_cobis:
            header = "COBIs"
        else:
            header = "COHPs"

        return f"Complete {header} for {self.structure}"

    def as_dict(self) -> dict[str, Any]:
        """JSON-serializable dict representation of CompleteCohp."""
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "are_coops": self.are_coops,
            "are_cobis": self.are_cobis,
            "are_multi_center_cobis": self.are_multi_center_cobis,
            "efermi": self.efermi,
            "structure": self.structure.as_dict(),
            "energies": self.energies.tolist(),
            "COHP": {"average": {str(spin): pops.tolist() for spin, pops in self.cohp.items()}},
        }

        if self.icohp is not None:
            dct["ICOHP"] = {"average": {str(spin): pops.tolist() for spin, pops in self.icohp.items()}}

        for label in self.all_cohps:
            dct["COHP"] |= {label: {str(spin): pops.tolist() for spin, pops in self.all_cohps[label].cohp.items()}}
            icohp = self.all_cohps[label].icohp
            if icohp is not None:
                if "ICOHP" not in dct:
                    dct["ICOHP"] = {label: {str(spin): pops.tolist() for spin, pops in icohp.items()}}
                else:
                    dct["ICOHP"] |= {label: {str(spin): pops.tolist() for spin, pops in icohp.items()}}

        if False in [bond_dict == {} for bond_dict in self.bonds.values()]:
            dct["bonds"] = {
                bond: {
                    "length": self.bonds[bond]["length"],
                    "sites": [site.as_dict() for site in self.bonds[bond]["sites"]],
                }
                for bond in self.bonds
            }

        if self.orb_res_cohp:
            orb_dict: dict[str, Any] = {}
            for label in self.orb_res_cohp:
                orb_dict[label] = {}
                for orbs in self.orb_res_cohp[label]:
                    orb_dict[label][orbs] = {
                        "COHP": {
                            str(spin): pops.tolist() for spin, pops in self.orb_res_cohp[label][orbs]["COHP"].items()
                        },
                        "ICOHP": {
                            str(spin): pops.tolist() for spin, pops in self.orb_res_cohp[label][orbs]["ICOHP"].items()
                        },
                        "orbitals": [[orb[0], orb[1].name] for orb in self.orb_res_cohp[label][orbs]["orbitals"]],
                    }

            dct["orb_res_cohp"] = orb_dict

        return dct

    def get_cohp_by_label(
        self,
        label: str,
        summed_spin_channels: bool = False,
    ) -> Cohp:
        """Get specific Cohp by the label, to simplify plotting.

        Args:
            label (str): Label for the interaction.
            summed_spin_channels (bool): Sum the spin channels and return the sum as Spin.up.

        Returns:
            The Cohp.
        """
        if label.lower() == "average":
            divided_cohp: dict[Spin, Any] | None = self.cohp
            divided_icohp: dict[Spin, Any] | None = self.icohp
        else:
            divided_cohp = self.all_cohps[label].get_cohp(spin=None, integrated=False)
            divided_icohp = self.all_cohps[label].get_icohp(spin=None)

        if divided_cohp is None:
            raise ValueError("divided_cohp is None")

        if summed_spin_channels and Spin.down in self.cohp:
            if divided_icohp is None:
                raise ValueError("divided_icohp is None")
            final_cohp: dict[Spin, Any] = {Spin.up: np.sum([divided_cohp[Spin.up], divided_cohp[Spin.down]], axis=0)}
            final_icohp: dict[Spin, Any] | None = {
                Spin.up: np.sum([divided_icohp[Spin.up], divided_icohp[Spin.down]], axis=0)
            }
        else:
            final_cohp = divided_cohp
            final_icohp = divided_icohp

        return Cohp(
            efermi=self.efermi,
            energies=self.energies,
            cohp=final_cohp,
            are_coops=self.are_coops,
            are_cobis=self.are_cobis,
            icohp=final_icohp,
        )

    def get_summed_cohp_by_label_list(
        self,
        label_list: list[str],
        divisor: float = 1,
        summed_spin_channels: bool = False,
    ) -> Cohp:
        """Get a Cohp object that includes a summed COHP divided by divisor.

        Args:
            label_list (list[str]): Labels for the COHP to include.
            divisor (float): The summed COHP will be divided by this divisor.
            summed_spin_channels (bool): Sum the spin channels and return the sum in Spin.up.

        Returns:
            A Cohp object for the summed COHP.
        """
        # Check if COHPs are spin polarized
        first_cohpobject = self.get_cohp_by_label(label_list[0])
        summed_cohp = first_cohpobject.cohp.copy()
        if first_cohpobject.icohp is None:
            raise ValueError("icohp of first_cohpobject is None")
        summed_icohp = first_cohpobject.icohp.copy()
        for label in label_list[1:]:
            cohp = self.get_cohp_by_label(label)
            icohp = cohp.icohp
            if icohp is None:
                raise ValueError("icohp is None")

            summed_cohp[Spin.up] = np.sum([summed_cohp[Spin.up], cohp.cohp[Spin.up]], axis=0)

            if Spin.down in summed_cohp:
                summed_cohp[Spin.down] = np.sum([summed_cohp[Spin.down], cohp.cohp[Spin.down]], axis=0)

            summed_icohp[Spin.up] = np.sum([summed_icohp[Spin.up], icohp[Spin.up]], axis=0)

            if Spin.down in summed_icohp:
                summed_icohp[Spin.down] = np.sum([summed_icohp[Spin.down], icohp[Spin.down]], axis=0)

        divided_cohp = {Spin.up: np.divide(summed_cohp[Spin.up], divisor)}
        divided_icohp = {Spin.up: np.divide(summed_icohp[Spin.up], divisor)}
        if Spin.down in summed_cohp:
            divided_cohp[Spin.down] = np.divide(summed_cohp[Spin.down], divisor)
            divided_icohp[Spin.down] = np.divide(summed_icohp[Spin.down], divisor)

        if summed_spin_channels and Spin.down in summed_cohp:
            final_cohp = {Spin.up: np.sum([divided_cohp[Spin.up], divided_cohp[Spin.down]], axis=0)}
            final_icohp = {Spin.up: np.sum([divided_icohp[Spin.up], divided_icohp[Spin.down]], axis=0)}
        else:
            final_cohp = divided_cohp
            final_icohp = divided_icohp

        return Cohp(
            efermi=first_cohpobject.efermi,
            energies=first_cohpobject.energies,
            cohp=final_cohp,
            are_coops=first_cohpobject.are_coops,
            are_cobis=first_cohpobject.are_coops,
            icohp=final_icohp,
        )

    def get_summed_cohp_by_label_and_orbital_list(
        self,
        label_list: list[str],
        orbital_list: list[str],
        divisor: float = 1,
        summed_spin_channels: bool = False,
    ) -> Cohp:
        """Get a Cohp object that includes a summed COHP divided by divisor.

        Args:
            label_list (list[str]): Labels for the COHP that should be included.
            orbital_list (list[str]): Orbitals for the COHPs that should be included
                (same order as label_list).
            divisor (float): The summed COHP will be divided by this divisor.
            summed_spin_channels (bool): Sum the spin channels and return the sum in Spin.up.

        Returns:
            A Cohp object including the summed COHP.
        """
        # Check length of label_list and orbital_list
        if not len(label_list) == len(orbital_list):
            raise ValueError("label_list and orbital_list don't have the same length!")

        # Check if COHPs are spin polarized
        first_cohpobject = self.get_orbital_resolved_cohp(label_list[0], orbital_list[0])
        if first_cohpobject is None:
            raise ValueError("first_cohpobject is None")
        if first_cohpobject.icohp is None:
            raise ValueError("icohp of first_cohpobject is None")
        summed_cohp = first_cohpobject.cohp.copy()
        summed_icohp = first_cohpobject.icohp.copy()

        for idx, label in enumerate(label_list[1:], start=1):
            cohp = self.get_orbital_resolved_cohp(label, orbital_list[idx])
            if cohp is None:
                raise ValueError("cohp is None.")
            if cohp.icohp is None:
                raise ValueError("icohp of cohp is None.")
            summed_cohp[Spin.up] = np.sum([summed_cohp[Spin.up], cohp.cohp.copy()[Spin.up]], axis=0)
            if Spin.down in summed_cohp:
                summed_cohp[Spin.down] = np.sum([summed_cohp[Spin.down], cohp.cohp.copy()[Spin.down]], axis=0)

            summed_icohp[Spin.up] = np.sum([summed_icohp[Spin.up], cohp.icohp.copy()[Spin.up]], axis=0)
            if Spin.down in summed_icohp:
                summed_icohp[Spin.down] = np.sum([summed_icohp[Spin.down], cohp.icohp.copy()[Spin.down]], axis=0)

        divided_cohp = {Spin.up: np.divide(summed_cohp[Spin.up], divisor)}
        divided_icohp = {Spin.up: np.divide(summed_icohp[Spin.up], divisor)}
        if Spin.down in summed_cohp:
            divided_cohp[Spin.down] = np.divide(summed_cohp[Spin.down], divisor)
            divided_icohp[Spin.down] = np.divide(summed_icohp[Spin.down], divisor)

        if summed_spin_channels and Spin.down in divided_cohp:
            final_cohp = {Spin.up: np.sum([divided_cohp[Spin.up], divided_cohp[Spin.down]], axis=0)}
            final_icohp = {Spin.up: np.sum([divided_icohp[Spin.up], divided_icohp[Spin.down]], axis=0)}
        else:
            final_cohp = divided_cohp
            final_icohp = divided_icohp

        return Cohp(
            efermi=first_cohpobject.efermi,
            energies=first_cohpobject.energies,
            cohp=final_cohp,
            are_coops=first_cohpobject.are_coops,
            are_cobis=first_cohpobject.are_cobis,
            icohp=final_icohp,
        )

    def get_orbital_resolved_cohp(
        self,
        label: str,
        orbitals: str | list[tuple[str, Orbital]] | tuple[tuple[str, Orbital], ...],
        summed_spin_channels: bool = False,
    ) -> Cohp | None:
        """Get orbital-resolved COHP.

        Args:
            label (str): Bond labels as in ICOHPLIST/ICOOPLIST.lobster.
            orbitals: The orbitals as a label, or list/tuple of
                [(n1, orbital1), (n2, orbital2), ...].
                Where each orbital can either be str, int, or Orbital.
            summed_spin_channels (bool): Sum the spin channels and return the sum as Spin.up.

        Returns:
            A Cohp object if CompleteCohp contains orbital-resolved COHP,
                or None if it doesn't.

        Note: It currently assumes that orbitals are str if they aren't the
            other valid types. This is not ideal, but is the easiest way to
            avoid unicode issues between Python 2 and Python 3.
        """
        if self.orb_res_cohp is None:
            return None

        if isinstance(orbitals, list | tuple):
            cohp_orbs = [val["orbitals"] for val in self.orb_res_cohp[label].values()]
            orbs = []
            for orbital in orbitals:
                if isinstance(orbital[1], int):
                    orbs.append((orbital[0], Orbital(orbital[1])))
                elif isinstance(orbital[1], Orbital):
                    orbs.append((orbital[0], orbital[1]))
                elif isinstance(orbital[1], str):
                    orbs.append((orbital[0], Orbital[orbital[1]]))
                else:
                    raise TypeError("Orbital must be str, int, or Orbital.")
            orb_index = cohp_orbs.index(orbs)
            orb_label = list(self.orb_res_cohp[label])[orb_index]

        elif isinstance(orbitals, str):
            orb_label = orbitals
        else:
            raise TypeError("Orbitals must be str, list, or tuple.")

        try:
            icohp = self.orb_res_cohp[label][orb_label]["ICOHP"]
        except KeyError:
            icohp = None

        start_cohp = self.orb_res_cohp[label][orb_label]["COHP"]
        start_icohp = icohp

        if summed_spin_channels and Spin.down in start_cohp:
            final_cohp = {}
            final_icohp = {}
            final_cohp[Spin.up] = np.sum([start_cohp[Spin.up], start_cohp[Spin.down]], axis=0)
            if start_icohp is not None:
                final_icohp[Spin.up] = np.sum([start_icohp[Spin.up], start_icohp[Spin.down]], axis=0)
        else:
            final_cohp = start_cohp
            final_icohp = start_icohp

        return Cohp(
            self.efermi,
            self.energies,
            final_cohp,
            icohp=final_icohp,
            are_coops=self.are_coops,
            are_cobis=self.are_cobis,
        )

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """Get CompleteCohp object from a dict representation.

        TODO: Clean this up.
        """
        cohp_dict = {}
        efermi = dct["efermi"]
        energies = dct["energies"]
        structure = Structure.from_dict(dct["structure"])
        are_cobis = dct.get("are_cobis", False)
        are_multi_center_cobis = dct.get("are_multi_center_cobis", False)
        are_coops = dct["are_coops"]
        avg_cohp = None

        if "bonds" in dct:
            bonds = {
                bond: {
                    "length": dct["bonds"][bond]["length"],
                    "sites": tuple(PeriodicSite.from_dict(site) for site in dct["bonds"][bond]["sites"]),
                    "cells": dct["bonds"][bond].get("cells", None),
                }
                for bond in dct["bonds"]
            }
        else:
            bonds = None

        for label in dct["COHP"]:
            cohp = {Spin(int(spin)): np.array(dct["COHP"][label][spin]) for spin in dct["COHP"][label]}
            try:
                icohp = {Spin(int(spin)): np.array(dct["ICOHP"][label][spin]) for spin in dct["ICOHP"][label]}
            except KeyError:
                icohp = None
            if label == "average":
                avg_cohp = Cohp(
                    efermi,
                    energies,
                    cohp,
                    icohp=icohp,
                    are_coops=are_coops,
                    are_cobis=are_cobis,
                    are_multi_center_cobis=are_multi_center_cobis,
                )
            else:
                cohp_dict[label] = Cohp(efermi, energies, cohp, icohp=icohp)

        if "orb_res_cohp" in dct:
            orb_cohp: dict[str, dict[Orbital, dict[str, Any]]] = {}
            for label in dct["orb_res_cohp"]:
                orb_cohp[label] = {}
                for orb in dct["orb_res_cohp"][label]:
                    cohp = {
                        Spin(int(s)): np.array(dct["orb_res_cohp"][label][orb]["COHP"][s], dtype=float)
                        for s in dct["orb_res_cohp"][label][orb]["COHP"]
                    }
                    try:
                        icohp = {
                            Spin(int(s)): np.array(dct["orb_res_cohp"][label][orb]["ICOHP"][s], dtype=float)
                            for s in dct["orb_res_cohp"][label][orb]["ICOHP"]
                        }
                    except KeyError:
                        icohp = None
                    orbitals = [(int(o[0]), Orbital[o[1]]) for o in dct["orb_res_cohp"][label][orb]["orbitals"]]
                    orb_cohp[label][orb] = {
                        "COHP": cohp,
                        "ICOHP": icohp,
                        "orbitals": orbitals,
                    }
                # If no total COHPs are present, calculate the total
                # COHPs from the single-orbital populations.
                # Total COHPs may not be present when the COHP generator keyword
                # is used in LOBSTER versions 2.2.0 and earlier.
                if label not in dct["COHP"] or dct["COHP"][label] is None:
                    cohp = {
                        Spin.up: np.sum(
                            np.array([orb_cohp[label][orb]["COHP"][Spin.up] for orb in orb_cohp[label]]),
                            axis=0,
                        )
                    }
                    try:
                        cohp[Spin.down] = np.sum(
                            np.array([orb_cohp[label][orb]["COHP"][Spin.down] for orb in orb_cohp[label]]),
                            axis=0,
                        )
                    except KeyError:
                        pass

                orb_res_icohp = None in [orb_cohp[label][orb]["ICOHP"] for orb in orb_cohp[label]]
                if (label not in dct["ICOHP"] or dct["ICOHP"][label] is None) and orb_res_icohp:
                    icohp = {
                        Spin.up: np.sum(
                            np.array([orb_cohp[label][orb]["ICOHP"][Spin.up] for orb in orb_cohp[label]]),
                            axis=0,
                        )
                    }
                    try:
                        icohp[Spin.down] = np.sum(
                            np.array([orb_cohp[label][orb]["ICOHP"][Spin.down] for orb in orb_cohp[label]]),
                            axis=0,
                        )
                    except KeyError:
                        pass
        else:
            orb_cohp = {}

        if avg_cohp is None:
            raise ValueError("avg_cohp is None")
        return cls(
            structure,
            avg_cohp,
            cohp_dict,
            bonds=bonds,
            are_coops=dct["are_coops"],
            are_cobis=dct.get("are_cobis", False),
            are_multi_center_cobis=are_multi_center_cobis,
            orb_res_cohp=orb_cohp,
        )

    @classmethod
    def from_file(
        cls,
        fmt: Literal["LMTO", "LOBSTER"],
        filename: PathLike | None = None,
        structure_file: PathLike | None = None,
        are_coops: bool = False,
        are_cobis: bool = False,
        are_multi_center_cobis: bool = False,
    ) -> Self:
        """Create CompleteCohp from an output file of a COHP calculation.

        Args:
            fmt (Literal["LMTO", "LOBSTER"]): The code used to calculate COHPs.
            filename (PathLike): The COHP output file. Defaults to "COPL"
                for LMTO and "COHPCAR.lobster/COOPCAR.lobster" for LOBSTER.
            structure_file (PathLike): The file containing the structure.
                If None, use "CTRL" for LMTO and "POSCAR" for LOBSTER.
            are_coops (bool): Whether the populations are COOPs or COHPs.
                Defaults to False for COHPs.
            are_cobis (bool): Whether the populations are COBIs or COHPs.
                Defaults to False for COHPs.
            are_multi_center_cobis (bool): Whether this file
                includes information on multi-center COBIs.

        Returns:
            A CompleteCohp object.
        """
        if are_coops and are_cobis:
            raise ValueError("You cannot have info about COOPs and COBIs in the same file.")

        fmt = fmt.upper()  # type: ignore[assignment]
        if fmt == "LMTO":
            # TODO: LMTO COOPs and orbital-resolved COHP cannot be handled yet
            are_coops = are_cobis = False
            orb_res_cohp = None
            if structure_file is None:
                structure_file = "CTRL"
            if filename is None:
                filename = "COPL"

            cohp_file: LMTOCopl | Cohpcar = LMTOCopl(filename=filename, to_eV=True)

        elif fmt == "LOBSTER":
            if (
                (are_coops and are_cobis)
                or (are_coops and are_multi_center_cobis)
                or (are_cobis and are_multi_center_cobis)
            ):
                raise ValueError("You cannot have info about COOPs, COBIs and/or multi-center COBIS in the same file.")
            if structure_file is None:
                structure_file = "POSCAR"
            if filename is None and filename is None:
                if are_coops:
                    filename = "COOPCAR.lobster"
                elif are_cobis or are_multi_center_cobis:
                    filename = "COBICAR.lobster"
                else:
                    filename = "COHPCAR.lobster"
            cohp_file = Cohpcar(
                filename=filename,
                are_coops=are_coops,
                are_cobis=are_cobis,
                are_multi_center_cobis=are_multi_center_cobis,
            )
            orb_res_cohp = cohp_file.orb_res_cohp

        else:
            raise ValueError(f"Unknown format {fmt}. Valid formats are LMTO and LOBSTER.")

        structure = Structure.from_file(structure_file)
        efermi = cohp_file.efermi
        cohp_data = cohp_file.cohp_data
        energies = cohp_file.energies

        # LOBSTER shifts the energies so that the Fermi level is at zero.
        # Shifting should be done by the plotter object though.
        spins = [Spin.up, Spin.down] if cohp_file.is_spin_polarized else [Spin.up]
        if fmt == "LOBSTER":
            energies += efermi

        if orb_res_cohp is not None:
            # If no total COHPs are present, calculate the total
            # COHPs from the single-orbital populations. Total COHPs
            # may not be present when the cohpgenerator keyword is used
            # in LOBSTER versions 2.2.0 and earlier.

            # TODO: Test this more extensively

            for label in orb_res_cohp:
                if cohp_file.cohp_data[label]["COHP"] is None:
                    cohp_data[label]["COHP"] = {
                        sp: np.sum(
                            [orb_res_cohp[label][orbs]["COHP"][sp] for orbs in orb_res_cohp[label]],
                            axis=0,
                        )
                        for sp in spins
                    }
                if cohp_file.cohp_data[label]["ICOHP"] is None:
                    cohp_data[label]["ICOHP"] = {
                        sp: np.sum(
                            [orb_res_cohp[label][orbs]["ICOHP"][sp] for orbs in orb_res_cohp[label]],
                            axis=0,
                        )
                        for sp in spins
                    }

        if fmt == "LMTO":
            # Calculate the average COHP for the LMTO file to be consistent with LOBSTER
            avg_data: dict[Literal["COHP", "ICOHP"], dict] = {"COHP": {}, "ICOHP": {}}
            for dtype in avg_data:
                for spin in spins:
                    rows = np.array([v[dtype][spin] for v in cohp_data.values()])
                    avg = np.mean(rows, axis=0)
                    # LMTO COHPs have 5 significant digits
                    avg_data[dtype] |= {spin: np.array([round_to_sigfigs(a, 5) for a in avg], dtype=float)}
            avg_cohp = Cohp(efermi, energies, avg_data["COHP"], icohp=avg_data["ICOHP"])

        elif not are_multi_center_cobis:
            avg_cohp = Cohp(
                efermi,
                energies,
                cohp_data["average"]["COHP"],
                icohp=cohp_data["average"]["ICOHP"],
                are_coops=are_coops,
                are_cobis=are_cobis,
                are_multi_center_cobis=are_multi_center_cobis,
            )
            del cohp_data["average"]

        else:
            # Only include two-center COBIs in average for both spin channels
            cohp = {}
            cohp[Spin.up] = np.array(
                [np.array(c["COHP"][Spin.up]) for c in cohp_file.cohp_data.values() if len(c["sites"]) <= 2]
            ).mean(axis=0)
            try:
                cohp[Spin.down] = np.array(
                    [np.array(c["COHP"][Spin.down]) for c in cohp_file.cohp_data.values() if len(c["sites"]) <= 2]
                ).mean(axis=0)
            except KeyError:
                pass

            try:
                icohp = {}
                icohp[Spin.up] = np.array(
                    [np.array(c["ICOHP"][Spin.up]) for c in cohp_file.cohp_data.values() if len(c["sites"]) <= 2]
                ).mean(axis=0)
                try:
                    icohp[Spin.down] = np.array(
                        [np.array(c["ICOHP"][Spin.down]) for c in cohp_file.cohp_data.values() if len(c["sites"]) <= 2]
                    ).mean(axis=0)
                except KeyError:
                    pass
            except KeyError:
                icohp = None

            avg_cohp = Cohp(
                efermi,
                energies,
                cohp,
                icohp=icohp,
                are_coops=are_coops,
                are_cobis=are_cobis,
                are_multi_center_cobis=are_multi_center_cobis,
            )

        cohp_dict = {
            key: Cohp(
                efermi,
                energies,
                dct["COHP"],
                icohp=dct["ICOHP"],
                are_coops=are_coops,
                are_cobis=are_cobis,
                are_multi_center_cobis=are_multi_center_cobis,
            )
            for key, dct in cohp_data.items()
        }

        bond_dict = {
            key: {
                "length": dct["length"],
                "sites": [structure[site] for site in dct["sites"]],
            }
            for key, dct in cohp_data.items()
        }

        return cls(
            structure,
            avg_cohp,
            cohp_dict,
            bonds=bond_dict,
            are_coops=are_coops,
            are_cobis=are_cobis,
            orb_res_cohp=orb_res_cohp,
        )


class IcohpValue(MSONable):
    """Information for an ICOHP or ICOOP value.

    Attributes:
        energies (NDArray): Energy values for the COHP/ICOHP/COOP/ICOOP.
        densities (NDArray): Density of states for the COHP/ICOHP/COOP/ICOOP.
        energies_are_cartesian (bool): Whether the energies are cartesian.
        are_coops (bool): Whether the object is COOP/ICOOP.
        are_cobis (bool): Whether the object is COBIS/ICOBIS.
        icohp (dict): The ICOHP/COHP values, whose keys are Spin.up and Spin.down.
        summed_icohp (float): The summed ICOHP/COHP values.
        num_bonds (int): The number of bonds used for the average COHP (for LOBSTER versions <3.0).
    """

    def __init__(
        self,
        label: str,
        atom1: str,
        atom2: str,
        length: float,
        translation: Vector3D,
        num: int,
        icohp: dict[Spin, float],
        are_coops: bool = False,
        are_cobis: bool = False,
        orbitals: dict[str, dict[Literal["icohp", "orbitals"], Any]] | None = None,
    ) -> None:
        """
        Args:
            label (str): Label for the ICOHP.
            atom1 (str): The first atom that contributes to the bond.
            atom2 (str): The second atom that contributes to the bond.
            length (float): Bond length.
            translation (Vector3D): cell translation vector, e.g. (0, 0, 0).
            num (int): The number of equivalent bonds.
            icohp (dict[Spin, float]): {Spin.up: ICOHP_up, Spin.down: ICOHP_down}
            are_coops (bool): Whether these are COOPs.
            are_cobis (bool): Whether these are COBIs.
            orbitals (dict): {[str(Orbital1)-str(Orbital2)]: {
                "icohp": {
                            Spin.up: IcohpValue for spin.up,
                            Spin.down: IcohpValue for spin.down
                        },
                "orbitals": [Orbital1, Orbital2, ...]}.
        """
        if are_coops and are_cobis:
            raise ValueError("You cannot have info about COOPs and COBIs in the same file.")

        self._are_coops = are_coops
        self._are_cobis = are_cobis
        self._label = label
        self._atom1 = atom1
        self._atom2 = atom2
        self._length = length
        self._translation = translation
        self._num = num
        self._icohp = icohp
        self._orbitals = orbitals
        self._is_spin_polarized = Spin.down in self._icohp

    def __str__(self) -> str:
        """String representation of the ICOHP/ICOOP."""
        # (are_coops and are_cobis) is never True
        if self._are_coops:
            header = "ICOOP"
        elif self._are_cobis:
            header = "ICOBI"
        else:
            header = "ICOHP"

        if self._is_spin_polarized:
            return (
                f"{header} {self._label} between {self._atom1} and {self._atom2} ({self._translation}): "
                f"{self._icohp[Spin.up]} eV (Spin up) and {self._icohp[Spin.down]} eV (Spin down)"
            )

        return (
            f"{header} {self._label} between {self._atom1} and {self._atom2} ({self._translation}): "
            f"{self._icohp[Spin.up]} eV (Spin up)"
        )

    @property
    def num_bonds(self) -> int:
        """The number of bonds for which the ICOHP value is an average.

        Returns:
            int
        """
        return self._num

    @property
    def are_coops(self) -> bool:
        """Whether these are ICOOPs.

        Returns:
            bool
        """
        return self._are_coops

    @property
    def are_cobis(self) -> bool:
        """Whether these are ICOBIs.

        Returns:
            bool
        """
        return self._are_cobis

    @property
    def is_spin_polarized(self) -> bool:
        """Whether this is a spin polarized calculation.

        Returns:
            bool
        """
        return self._is_spin_polarized

    def icohpvalue(self, spin: Spin = Spin.up) -> float:
        """
        Args:
            spin: Spin.up or Spin.down.

        Returns:
            float: ICOHP value corresponding to chosen spin.
        """
        if not self.is_spin_polarized and spin == Spin.down:
            raise ValueError("The calculation was not performed with spin polarization")

        return self._icohp[spin]

    def icohpvalue_orbital(
        self,
        orbitals: tuple[Orbital, Orbital] | str,
        spin: Spin = Spin.up,
    ) -> float:
        """
        Args:
            orbitals: tuple[Orbital, Orbital] or "str(Orbital0)-str(Orbital1)".
            spin (Spin): Spin.up or Spin.down.

        Returns:
            float: ICOHP value corresponding to chosen spin.
        """
        if not self.is_spin_polarized and spin == Spin.down:
            raise ValueError("The calculation was not performed with spin polarization")

        if isinstance(orbitals, tuple | list):
            orbitals = f"{orbitals[0]}-{orbitals[1]}"

        if self._orbitals is None:
            raise ValueError("self._orbitals is None")
        return self._orbitals[orbitals]["icohp"][spin]

    @property
    def icohp(self) -> dict[Spin, float]:
        """Dict with ICOHPs for spin up and spin down.

        Returns:
            dict[Spin, float]: {Spin.up: ICOHP_up, Spin.down: ICOHP_down}.
        """
        return self._icohp

    @property
    def summed_icohp(self) -> float:
        """Summed ICOHPs of both spin channels if spin polarized.

        Returns:
            float: ICOHP value in eV.
        """
        return self._icohp[Spin.down] + self._icohp[Spin.up] if self._is_spin_polarized else self._icohp[Spin.up]

    @property
    def summed_orbital_icohp(self) -> dict[str, float]:
        """Summed orbital-resolved ICOHPs of both spin channels if spin-polarized.

        Returns:
            dict[str, float]: "str(Orbital1)-str(Ortibal2)": ICOHP value in eV.
        """
        if self._orbitals is None:
            raise ValueError("_orbitals attrib is None.")

        orbital_icohp = {}
        for orb, item in self._orbitals.items():
            orbital_icohp[orb] = (
                item["icohp"][Spin.up] + item["icohp"][Spin.down] if self._is_spin_polarized else item["icohp"][Spin.up]
            )
        return orbital_icohp


class IcohpCollection(MSONable):
    """Collection of IcohpValues.

    Attributes:
        are_coops (bool): Whether these are ICOOPs.
        are_cobis (bool): Whether these are ICOOPs.
        is_spin_polarized (bool): Whether the calculation is spin polarized.
    """

    def __init__(
        self,
        list_labels: list[str],
        list_atom1: list[str],
        list_atom2: list[str],
        list_length: list[float],
        list_translation: list[Vector3D],
        list_num: list[int],
        list_icohp: list[dict[Spin, float]],
        is_spin_polarized: bool,
        list_orb_icohp: list[dict[str, dict[Literal["icohp", "orbitals"], Any]]] | None = None,
        are_coops: bool = False,
        are_cobis: bool = False,
    ) -> None:
        """
        Args:
            list_labels (list[str]): Labels for ICOHP/ICOOP values.
            list_atom1 (list[str]): Atom names, e.g. "O1".
            list_atom2 (list[str]): Atom names, e.g. "O1".
            list_length (list[float]): Bond lengths in Angstrom.
            list_translation (list[Vector3D]): Cell translation vectors.
            list_num (list[int]): Numbers of equivalent bonds, usually 1 starting from LOBSTER 3.0.0.
            list_icohp (list[dict]): Dicts as {Spin.up: ICOHP_up, Spin.down: ICOHP_down}.
            is_spin_polarized (bool): Whether the calculation is spin polarized.
            list_orb_icohp (list[dict]): Dicts as {[str(Orbital1)-str(Orbital2)]: {
                "icohp": {Spin.up: IcohpValue for spin.up, Spin.down: IcohpValue for spin.down},
                "orbitals": [Orbital1, Orbital2]}.
            are_coops (bool): Whether ICOOPs are stored.
            are_cobis (bool): Whether ICOBIs are stored.
        """
        if are_coops and are_cobis:
            raise ValueError("You cannot have info about COOPs and COBIs in the same file.")

        self._are_coops = are_coops
        self._are_cobis = are_cobis
        self._is_spin_polarized = is_spin_polarized
        self._list_labels = list_labels
        self._list_atom1 = list_atom1
        self._list_atom2 = list_atom2
        self._list_length = list_length
        self._list_translation = list_translation
        self._list_num = list_num
        self._list_icohp = list_icohp
        self._list_orb_icohp = list_orb_icohp

        # TODO: DanielYang: self._icohplist name is misleading
        # (not list), and confuses with self._list_icohp
        self._icohplist: dict[str, IcohpValue] = {}
        for idx, label in enumerate(list_labels):
            self._icohplist[label] = IcohpValue(
                label=label,
                atom1=list_atom1[idx],
                atom2=list_atom2[idx],
                length=list_length[idx],
                translation=list_translation[idx],
                num=list_num[idx],
                icohp=list_icohp[idx],
                are_coops=are_coops,
                are_cobis=are_cobis,
                orbitals=None if list_orb_icohp is None else list_orb_icohp[idx],
            )

    def __str__(self) -> str:
        return "\n".join([str(value) for value in self._icohplist.values()])

    def get_icohp_by_label(
        self,
        label: str,
        summed_spin_channels: bool = True,
        spin: Spin = Spin.up,
        orbitals: str | tuple[Orbital, Orbital] | None = None,
    ) -> float:
        """Get an ICOHP value for a certain bond indicated by the label.

        Args:
            label (str): The bond number in Icohplist.lobster/Icooplist.lobster,
                starting from "1".
            summed_spin_channels (bool): Whether the ICOHPs/ICOOPs of both
                spin channels should be summed.
            spin (Spin): If not summed_spin_channels, indicate
                which spin channel should be returned.
            orbitals: List of Orbital or "str(Orbital1)-str(Orbital2)".

        Returns:
            float: ICOHP/ICOOP value.
        """
        icohp: IcohpValue = self._icohplist[label]

        if orbitals is None:
            return icohp.summed_icohp if summed_spin_channels else icohp.icohpvalue(spin)

        if isinstance(orbitals, tuple | list):
            orbitals = f"{orbitals[0]}-{orbitals[1]}"

        if summed_spin_channels:
            return icohp.summed_orbital_icohp[orbitals]

        return icohp.icohpvalue_orbital(spin=spin, orbitals=orbitals)

    def get_summed_icohp_by_label_list(
        self,
        label_list: list[str],
        divisor: float = 1.0,
        summed_spin_channels: bool = True,
        spin: Spin = Spin.up,
    ) -> float:
        """Get the sum of ICOHP values.

        Args:
            label_list (list[str]): Labels of the ICOHPs/ICOOPs that should be summed,
                the same as in ICOHPLIST/ICOOPLIST.
            divisor (float): Divisor used to divide the sum.
            summed_spin_channels (bool): Whether the ICOHPs/ICOOPs of both
                spin channels should be summed.
            spin (Spin): If not summed_spin_channels, indicate
                which spin channel should be returned.

        Returns:
            float: Sum of ICOHPs selected with label_list.
        """
        sum_icohp: float = 0
        for label in label_list:
            icohp = self._icohplist[label]
            if icohp.num_bonds != 1:
                warnings.warn("One of the ICOHP values is an average over bonds. This is currently not considered.")

            if icohp._is_spin_polarized and summed_spin_channels:
                sum_icohp += icohp.summed_icohp
            else:
                sum_icohp += icohp.icohpvalue(spin)

        return sum_icohp / divisor

    def get_icohp_dict_by_bondlengths(
        self,
        minbondlength: float = 0.0,
        maxbondlength: float = 8.0,
    ) -> dict[str, IcohpValue]:
        """Get IcohpValues within certain bond length range.

        Args:
            minbondlength (float): The minimum bond length.
            maxbondlength (float): The maximum bond length.

        Returns:
            dict[str, IcohpValue]: Keys are the labels from the initial list_labels.
        """
        new_icohp_dict = {}
        for value in self._icohplist.values():
            if minbondlength <= value._length <= maxbondlength:
                new_icohp_dict[value._label] = value
        return new_icohp_dict

    def get_icohp_dict_of_site(
        self,
        site: int,
        minsummedicohp: float | None = None,
        maxsummedicohp: float | None = None,
        minbondlength: float = 0.0,
        maxbondlength: float = 8.0,
        only_bonds_to: list[str] | None = None,
    ) -> dict[str, IcohpValue]:
        """Get IcohpValues for a certain site.

        Args:
            site (int): The site of interest, ordered as in Icohplist.lobster/Icooplist.lobster,
                starts from 0.
            minsummedicohp (float): Minimal ICOHP/ICOOP of the bonds that are considered.
                It is the summed ICOHP value from both spin channels for spin polarized cases
            maxsummedicohp (float): Maximal ICOHP/ICOOP of the bonds that are considered.
                It is the summed ICOHP value from both spin channels for spin polarized cases
            minbondlength (float): The minimum bond length.
            maxbondlength (float): The maximum bond length.
            only_bonds_to (list[str]): The bonding partners that are allowed, e.g. ["O"].

        Returns:
            Dict of IcohpValues, the keys correspond to the values from the initial list_labels.
        """
        new_icohp_dict = {}
        for key, value in self._icohplist.items():
            atomnumber1 = int(re.split(r"(\d+)", value._atom1)[1]) - 1
            atomnumber2 = int(re.split(r"(\d+)", value._atom2)[1]) - 1
            if site in (atomnumber1, atomnumber2):
                # Swap order of atoms so that searched one is always atom1
                if site == atomnumber2:
                    save = value._atom1
                    value._atom1 = value._atom2
                    value._atom2 = save

                second_test = True if only_bonds_to is None else re.split("(\\d+)", value._atom2)[0] in only_bonds_to
                if minbondlength <= value._length <= maxbondlength and second_test:
                    # TODO: DanielYang: merge the following condition blocks
                    if minsummedicohp is not None:
                        if value.summed_icohp >= minsummedicohp:
                            if maxsummedicohp is not None:
                                if value.summed_icohp <= maxsummedicohp:
                                    new_icohp_dict[key] = value
                            else:
                                new_icohp_dict[key] = value
                    elif maxsummedicohp is not None:
                        if value.summed_icohp <= maxsummedicohp:
                            new_icohp_dict[key] = value
                    else:
                        new_icohp_dict[key] = value

        return new_icohp_dict

    def extremum_icohpvalue(
        self,
        summed_spin_channels: bool = True,
        spin: Spin = Spin.up,
    ) -> float:
        """Get ICOHP/ICOOP of the strongest bond.

        Args:
            summed_spin_channels (bool): Whether the ICOHPs/ICOOPs of both
                spin channels should be summed.
            spin (Spin): If not summed_spin_channels, this indicates which
                spin channel should be returned.

        Returns:
            Lowest ICOHP/largest ICOOP value (i.e. ICOHP/ICOOP value of strongest bond).
        """
        extremum = -sys.float_info.max if self._are_coops or self._are_cobis else sys.float_info.max

        if not self._is_spin_polarized:
            if spin == Spin.down:
                warnings.warn("This spin channel does not exist. I am switching to Spin.up")
            spin = Spin.up

        for value in self._icohplist.values():
            if not value.is_spin_polarized or not summed_spin_channels:
                if not self._are_coops and not self._are_cobis:
                    extremum = min(value.icohpvalue(spin), extremum)
                elif value.icohpvalue(spin) > extremum:
                    extremum = value.icohpvalue(spin)

            elif not self._are_coops and not self._are_cobis:
                extremum = min(value.summed_icohp, extremum)

            elif value.summed_icohp > extremum:
                extremum = value.summed_icohp

        return extremum

    @property
    def is_spin_polarized(self) -> bool:
        """Whether this is spin polarized."""
        return self._is_spin_polarized

    @property
    def are_coops(self) -> bool:
        """Whether this is COOP."""
        return self._are_coops

    @property
    def are_cobis(self) -> bool:
        """Whether this is COBI."""
        return self._are_cobis


def get_integrated_cohp_in_energy_range(
    cohp: CompleteCohp,
    label: str,
    orbital: str | None = None,
    energy_range: float | tuple[float, float] | None = None,
    relative_E_Fermi: bool = True,
    summed_spin_channels: bool = True,
) -> float | dict[Spin, float]:
    """Integrate CompleteCohps which include data of integrated COHPs (ICOHPs).

    Args:
        cohp (CompleteCohp): CompleteCohp object.
        label (str): Label of the COHP data.
        orbital (str): If not None, a orbital resolved integrated COHP will be returned.
        energy_range: If None, return the ICOHP value at Fermi level.
            If float, integrate from this value up to Fermi level.
            If (float, float), integrate in between.
        relative_E_Fermi (bool): Whether energy scale with Fermi level at 0 eV is chosen.
        summed_spin_channels (bool): Whether Spin channels will be summed.

    Returns:
        If summed_spin_channels:
            float: the ICOHP.
        else:
            dict: {Spin.up: float, Spin.down: float}
    """
    if orbital is None:
        icohps = cohp.all_cohps[label].get_icohp(spin=None)
    else:
        _icohps = cohp.get_orbital_resolved_cohp(label=label, orbitals=orbital)
        if _icohps is None:
            raise ValueError("_icohps is None")
        icohps = _icohps.icohp

    if icohps is None:
        raise ValueError("ichops is None")

    summedicohp = {}
    if summed_spin_channels and Spin.down in icohps:
        summedicohp[Spin.up] = icohps[Spin.up] + icohps[Spin.down]
    else:
        summedicohp = icohps

    if energy_range is None:
        energies_corrected = cohp.energies - cohp.efermi
        spl_spinup = InterpolatedUnivariateSpline(energies_corrected, summedicohp[Spin.up], ext=0)

        if not summed_spin_channels and Spin.down in icohps:
            spl_spindown = InterpolatedUnivariateSpline(energies_corrected, summedicohp[Spin.down], ext=0)
            return {Spin.up: spl_spinup(0.0), Spin.down: spl_spindown(0.0)}

        return spl_spinup(0.0) if summed_spin_channels else {Spin.up: spl_spinup(0.0)}

    # Return ICOHP value at the Fermi level
    if isinstance(energy_range, float):
        if relative_E_Fermi:
            energies_corrected = cohp.energies - cohp.efermi
            spl_spinup = InterpolatedUnivariateSpline(energies_corrected, summedicohp[Spin.up], ext=0)

            if not summed_spin_channels and Spin.down in icohps:
                spl_spindown = InterpolatedUnivariateSpline(energies_corrected, summedicohp[Spin.down], ext=0)
                return {
                    Spin.up: spl_spinup(0) - spl_spinup(energy_range),
                    Spin.down: spl_spindown(0) - spl_spindown(energy_range),
                }
            if summed_spin_channels:
                return spl_spinup(0) - spl_spinup(energy_range)
            return {Spin.up: spl_spinup(0) - spl_spinup(energy_range)}

        energies_corrected = cohp.energies
        spl_spinup = InterpolatedUnivariateSpline(energies_corrected, summedicohp[Spin.up], ext=0)

        if not summed_spin_channels and Spin.down in icohps:
            spl_spindown = InterpolatedUnivariateSpline(energies_corrected, summedicohp[Spin.down], ext=0)
            return {
                Spin.up: spl_spinup(cohp.efermi) - spl_spinup(energy_range),
                Spin.down: spl_spindown(cohp.efermi) - spl_spindown(energy_range),
            }
        if summed_spin_channels:
            return spl_spinup(cohp.efermi) - spl_spinup(energy_range)
        return {Spin.up: spl_spinup(cohp.efermi) - spl_spinup(energy_range)}

    energies_corrected = cohp.energies - cohp.efermi if relative_E_Fermi else cohp.energies

    spl_spinup = InterpolatedUnivariateSpline(energies_corrected, summedicohp[Spin.up], ext=0)

    if not summed_spin_channels and Spin.down in icohps:
        spl_spindown = InterpolatedUnivariateSpline(energies_corrected, summedicohp[Spin.down], ext=0)
        return {
            Spin.up: spl_spinup(energy_range[1]) - spl_spinup(energy_range[0]),
            Spin.down: spl_spindown(energy_range[1]) - spl_spindown(energy_range[0]),
        }
    if summed_spin_channels:
        return spl_spinup(energy_range[1]) - spl_spinup(energy_range[0])

    return {Spin.up: spl_spinup(energy_range[1]) - spl_spinup(energy_range[0])}
