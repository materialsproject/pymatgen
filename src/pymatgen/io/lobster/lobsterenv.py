"""This module provides classes to perform analyses of the local
environments (e.g., finding near neighbors) of single sites in molecules
and structures based on bonding analysis with LOBSTER.

If you use this module, please cite:
J. George, G. Petretto, A. Naik, M. Esters, A. J. Jackson, R. Nelson, R. Dronskowski, G.-M. Rignanese, G. Hautier,
"Automated Bonding Analysis with Crystal Orbital Hamilton Populations",
ChemPlusChem 2022, e202200123,
DOI: 10.1002/cplu.202200123.
"""

from __future__ import annotations

import collections
import copy
import math
import tempfile
from typing import TYPE_CHECKING, NamedTuple

import matplotlib as mpl
import numpy as np
from monty.dev import deprecated

from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.analysis.local_env import NearNeighbors
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import CohpPlotter
from pymatgen.io.lobster import Charge, Icohplist
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing import Any, Literal

    import matplotlib as mpl
    from numpy.typing import NDArray
    from typing_extensions import Self

    from pymatgen.core import PeriodicNeighbor, PeriodicSite, Structure
    from pymatgen.core.periodic_table import Element
    from pymatgen.electronic_structure.cohp import IcohpCollection, IcohpValue
    from pymatgen.util.typing import PathLike

__author__ = "Janine George"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "1.0"
__maintainer__ = "J. George"
__email__ = "janinegeorge.ulfen@gmail.com"
__status__ = "Production"
__date__ = "February 2, 2021"

due.cite(
    Doi("10.1002/cplu.202200123"),
    description="Automated Bonding Analysis with Crystal Orbital Hamilton Populations",
)


class LobsterNeighbors(NearNeighbors):
    """
    This class combines capabilities from LocalEnv and ChemEnv to determine
    coordination environments based on bonding analysis.
    """

    def __init__(
        self,
        structure: Structure,
        filename_icohp: PathLike | None = "ICOHPLIST.lobster",
        obj_icohp: Icohplist | None = None,
        are_coops: bool = False,
        are_cobis: bool = False,
        valences: list[float] | None = None,
        limits: tuple[float, float] | None = None,
        additional_condition: Literal[0, 1, 2, 3, 4, 5, 6] = 0,
        only_bonds_to: list[str] | None = None,
        perc_strength_icohp: float = 0.15,
        noise_cutoff: float = 0.1,
        valences_from_charges: bool = False,
        filename_charge: PathLike | None = None,
        obj_charge: Charge | None = None,
        which_charge: Literal["Mulliken", "Loewdin"] = "Mulliken",
        adapt_extremum_to_add_cond: bool = False,
        add_additional_data_sg: bool = False,
        filename_blist_sg1: PathLike | None = None,
        filename_blist_sg2: PathLike | None = None,
        id_blist_sg1: Literal["icoop", "icobi"] = "icoop",
        id_blist_sg2: Literal["icoop", "icobi"] = "icobi",
    ) -> None:
        """
        Args:
            filename_icohp (PathLike): Path to ICOHPLIST.lobster or
                ICOOPLIST.lobster or ICOBILIST.lobster.
            obj_icohp (Icohplist): Icohplist object.
            structure (Structure): Typically constructed by Structure.from_file("POSCAR").
            are_coops (bool): Whether the file is a ICOOPLIST.lobster (True) or a
                ICOHPLIST.lobster (False). Only tested for ICOHPLIST.lobster so far.
            are_cobis (bool): Whether the file is a ICOBILIST.lobster (True) or
                a ICOHPLIST.lobster (False).
            valences (list[float]): Valence/charge for each element.
            limits (tuple[float, float]): Range to decide which ICOHPs (ICOOP
                or ICOBI) should be considered.
            additional_condition (int): Additional condition that decides
                which kind of bonds will be considered:
                    0 - NO_ADDITIONAL_CONDITION
                    1 - ONLY_ANION_CATION_BONDS
                    2 - NO_ELEMENT_TO_SAME_ELEMENT_BONDS
                    3 - ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS
                    4 - ONLY_ELEMENT_TO_OXYGEN_BONDS
                    5 - DO_NOT_CONSIDER_ANION_CATION_BONDS
                    6 - ONLY_CATION_CATION_BONDS
            only_bonds_to (list[str]): Only consider bonds to certain elements (e.g. ["O"] for oxygen).
            perc_strength_icohp (float): If no "limits" are given, this will decide
                which ICOHPs will be considered (relative to the strongest ICOHP/ICOOP/ICOBI).
            noise_cutoff (float): The lower limit of ICOHPs considered.
            valences_from_charges (bool): If True and path to CHARGE.lobster is provided,
                will use LOBSTER charges (Mulliken) instead of valences.
            filename_charge (PathLike): Path to Charge.lobster.
            obj_charge (Charge): Charge object.
            which_charge ("Mulliken" | "Loewdin"): Source of charge.
            adapt_extremum_to_add_cond (bool): Whether to adapt the limits to only
                focus on the bonds determined by the additional condition.
            add_additional_data_sg (bool): Add the information from filename_add_bondinglist_sg1.
            filename_blist_sg1 (PathLike): Path to additional ICOOP, ICOBI data for structure graphs.
            filename_blist_sg2 (PathLike): Path to additional ICOOP, ICOBI data for structure graphs.
            id_blist_sg1 ("icoop" | "icobi"): Identity of data in filename_blist_sg1.
            id_blist_sg2 ("icoop" | "icobi"): Identity of data in filename_blist_sg2.
        """
        if filename_icohp is not None:
            self.ICOHP = Icohplist(are_coops=are_coops, are_cobis=are_cobis, filename=filename_icohp)
        elif obj_icohp is not None:
            self.ICOHP = obj_icohp
        else:
            raise ValueError("Please provide either filename_icohp or obj_icohp")

        self.Icohpcollection = self.ICOHP.icohpcollection
        self.structure = structure
        self.limits = limits
        self.only_bonds_to = only_bonds_to
        self.adapt_extremum_to_add_cond = adapt_extremum_to_add_cond
        self.are_coops = are_coops
        self.are_cobis = are_cobis
        self.add_additional_data_sg = add_additional_data_sg
        self.filename_blist_sg1 = filename_blist_sg1
        self.filename_blist_sg2 = filename_blist_sg2
        self.noise_cutoff = noise_cutoff

        self.id_blist_sg1 = id_blist_sg1.lower()
        self.id_blist_sg2 = id_blist_sg2.lower()

        allowed_arguments = {"icoop", "icobi"}
        if self.id_blist_sg1 not in allowed_arguments or self.id_blist_sg2 not in allowed_arguments:
            raise ValueError("Algorithm can only work with ICOOPs, ICOBIs")

        if add_additional_data_sg:
            if self.id_blist_sg1 == "icoop":
                are_coops_id1 = True
                are_cobis_id1 = False
            else:
                are_coops_id1 = False
                are_cobis_id1 = True

            self.bonding_list_1 = Icohplist(
                filename=self.filename_blist_sg1,
                are_coops=are_coops_id1,
                are_cobis=are_cobis_id1,
            )

            if self.id_blist_sg2 == "icoop":
                are_coops_id2 = True
                are_cobis_id2 = False
            else:
                are_coops_id2 = False
                are_cobis_id2 = True

            self.bonding_list_2 = Icohplist(
                filename=self.filename_blist_sg2,
                are_coops=are_coops_id2,
                are_cobis=are_cobis_id2,
            )

        # Check the additional condition
        if additional_condition not in range(7):
            raise ValueError(f"Unexpected {additional_condition=}, must be one of {list(range(7))}")
        self.additional_condition = additional_condition

        # Read in valences, will prefer manual setting of valences
        self.valences: list[float] | None
        if valences is None:
            if valences_from_charges and filename_charge is not None:
                chg = Charge(filename=filename_charge)
                if which_charge == "Mulliken":
                    self.valences = chg.Mulliken
                elif which_charge == "Loewdin":
                    self.valences = chg.Loewdin

            elif valences_from_charges and obj_charge is not None:
                chg = obj_charge
                if which_charge == "Mulliken":
                    self.valences = chg.Mulliken
                elif which_charge == "Loewdin":
                    self.valences = chg.Loewdin

            else:
                bv_analyzer = BVAnalyzer()
                try:
                    self.valences = bv_analyzer.get_valences(structure=self.structure)
                except ValueError as exc:
                    self.valences = None
                    if additional_condition in {1, 3, 5, 6}:
                        raise ValueError(
                            "Valences cannot be assigned, additional_conditions 1, 3, 5 and 6 will not work"
                        ) from exc
        else:
            self.valences = valences

        if np.allclose(self.valences or [], np.zeros_like(self.valences)) and additional_condition in {1, 3, 5, 6}:
            raise ValueError("All valences are equal to 0, additional_conditions 1, 3, 5 and 6 will not work")

        if limits is None:
            self.lowerlimit = self.upperlimit = None
        else:
            self.lowerlimit, self.upperlimit = limits

        # Evaluate coordination environments
        self._evaluate_ce(
            lowerlimit=self.lowerlimit,
            upperlimit=self.upperlimit,
            only_bonds_to=only_bonds_to,
            additional_condition=self.additional_condition,
            perc_strength_icohp=perc_strength_icohp,
            adapt_extremum_to_add_cond=adapt_extremum_to_add_cond,
        )

    @property
    def structures_allowed(self) -> Literal[True]:
        """Whether this LobsterNeighbors class can be used with Structure objects."""
        return True

    @property
    def molecules_allowed(self) -> Literal[False]:
        """Whether this LobsterNeighbors class can be used with Molecule objects."""
        return False

    @property
    def anion_types(self) -> set[Element]:
        """The set of anion types in crystal structure.

        Returns:
            set[Element]: Anions in the crystal structure.
        """
        if self.valences is None:
            raise ValueError("No cations and anions defined")

        anion_species = []
        for site, val in zip(self.structure, self.valences, strict=True):
            if val < 0.0:
                anion_species.append(site.specie)

        return set(anion_species)

    @deprecated(anion_types)
    def get_anion_types(self) -> set[Element]:
        return self.anion_types

    def get_nn_info(
        self,
        structure: Structure,
        n: int,
        use_weights: bool = False,
    ) -> dict[str, Any]:
        """Get coordination number (CN) of site by index.

        Args:
            structure (Structure): Input structure.
            n (int): Index of site for which to determine CN.
            use_weights (bool): Whether to use weights for computing
                the CN (True), or each coordinated site has equal weight (False).
                The former is not implemented yet.

        Raises:
            ValueError: If use_weights is True, or if arg "structure" and structure
                used to initialize LobsterNeighbors have different lengths.

        Returns:
            dict[str, Any]: coordination number and a list of nearest neighbors.
        """
        if use_weights:
            raise ValueError("LobsterEnv cannot use weights")

        if len(structure) != len(self.structure):
            raise ValueError(
                f"Length of structure ({len(structure)}) and LobsterNeighbors ({len(self.structure)}) differ"
            )

        return self.sg_list[n]  # type: ignore[return-value]

    def get_light_structure_environment(
        self,
        only_cation_environments: bool = False,
        only_indices: list[int] | None = None,
    ) -> LobsterLightStructureEnvironments:
        """Get a LobsterLightStructureEnvironments object if the structure
        only contains coordination environments smaller 13.

        Args:
            only_cation_environments (bool): Only return data for cations.
            only_indices (list[int]): Only evaluate indexes in this list.

        Returns:
            LobsterLightStructureEnvironments
        """
        lgf = LocalGeometryFinder()
        lgf.setup_structure(structure=self.structure)
        list_ce_symbols = []
        list_csm = []
        list_permut = []
        for idx, _neigh_coords in enumerate(self.list_coords):
            if (len(_neigh_coords)) > 13:
                raise ValueError("Environment cannot be determined. Number of neighbors is larger than 13.")
            # Avoid problems if _neigh_coords is empty
            if _neigh_coords != []:
                lgf.setup_local_geometry(isite=idx, coords=_neigh_coords, optimization=2)
                cncgsm = lgf.get_coordination_symmetry_measures(optimization=2)
                list_ce_symbols.append(min(cncgsm.items(), key=lambda t: t[1]["csm_wcs_ctwcc"])[0])
                list_csm.append(min(cncgsm.items(), key=lambda t: t[1]["csm_wcs_ctwcc"])[1]["csm_wcs_ctwcc"])
                list_permut.append(min(cncgsm.items(), key=lambda t: t[1]["csm_wcs_ctwcc"])[1]["indices"])
            else:
                list_ce_symbols.append(None)
                list_csm.append(None)
                list_permut.append(None)

        new_list_ce_symbols = []
        new_list_csm = []
        new_list_permut = []
        new_list_neighsite = []
        new_list_neighisite = []

        if only_indices is None:
            if not only_cation_environments:
                return LobsterLightStructureEnvironments.from_Lobster(
                    list_ce_symbol=list_ce_symbols,
                    list_csm=list_csm,
                    list_permutation=list_permut,
                    list_neighsite=self.list_neighsite,
                    list_neighisite=self.list_neighisite,
                    structure=self.structure,
                    valences=self.valences,
                )

            if self.valences is None:
                raise ValueError(f"{self.valences=}")

            for idx, val in enumerate(self.valences):
                if val >= 0.0:
                    new_list_ce_symbols.append(list_ce_symbols[idx])
                    new_list_csm.append(list_csm[idx])
                    new_list_permut.append(list_permut[idx])
                    new_list_neighisite.append(self.list_neighisite[idx])
                    new_list_neighsite.append(self.list_neighsite[idx])
                else:
                    new_list_ce_symbols.append(None)
                    new_list_csm.append(None)
                    new_list_permut.append([])
                    new_list_neighisite.append([])
                    new_list_neighsite.append([])

        else:
            for site_idx, _site in enumerate(self.structure):
                if site_idx in only_indices:
                    new_list_ce_symbols.append(list_ce_symbols[site_idx])
                    new_list_csm.append(list_csm[site_idx])
                    new_list_permut.append(list_permut[site_idx])
                    new_list_neighisite.append(self.list_neighisite[site_idx])
                    new_list_neighsite.append(self.list_neighsite[site_idx])
                else:
                    new_list_ce_symbols.append(None)
                    new_list_csm.append(None)
                    new_list_permut.append([])
                    new_list_neighisite.append([])
                    new_list_neighsite.append([])

        return LobsterLightStructureEnvironments.from_Lobster(
            list_ce_symbol=new_list_ce_symbols,
            list_csm=new_list_csm,
            list_permutation=new_list_permut,
            list_neighsite=new_list_neighsite,
            list_neighisite=new_list_neighisite,
            structure=self.structure,
            valences=self.valences,
        )

    def get_info_icohps_to_neighbors(
        self,
        isites: list[int] | None = None,
        onlycation_isites: bool = True,
    ) -> ICOHPNeighborsInfo:
        """Get information on the ICOHPs of neighbors for certain sites
        as identified by their site id.

        This is useful for plotting the COHPs (ICOOPLIST.lobster/
        ICOHPLIST.lobster/ICOBILIST.lobster) of a site in the structure.


        Args:
            isites (list[int]): Site IDs. If is None, all isites will be used
                to add the ICOHPs of the neighbors.
            onlycation_isites (bool): If True and if isite is None, will
                only analyse the cations sites.

        Returns:
            ICOHPNeighborsInfo
        """
        if self.valences is None and onlycation_isites:
            raise ValueError("No valences are provided")

        if isites is None:
            if onlycation_isites:
                if self.valences is None:
                    raise ValueError(f"{self.valences}=")

                isites = [idx for idx in range(len(self.structure)) if self.valences[idx] >= 0.0]
            else:
                isites = list(range(len(self.structure)))

        if self.Icohpcollection is None:
            raise ValueError(f"{self.Icohpcollection=}")

        summed_icohps: float = 0.0
        list_icohps: list[float] = []
        number_bonds: int = 0
        labels: list[str] = []
        atoms: list[list[str]] = []
        final_isites: list[int] = []
        for idx, _site in enumerate(self.structure):
            if idx in isites:
                for key, icohpsum in zip(self.list_keys[idx], self.list_icohps[idx], strict=True):
                    summed_icohps += icohpsum
                    list_icohps.append(icohpsum)
                    labels.append(key)
                    atoms.append(
                        [
                            self.Icohpcollection._list_atom1[int(key) - 1],
                            self.Icohpcollection._list_atom2[int(key) - 1],
                        ]
                    )
                    number_bonds += 1
                    final_isites.append(idx)
        return ICOHPNeighborsInfo(summed_icohps, list_icohps, number_bonds, labels, atoms, final_isites)

    def plot_cohps_of_neighbors(
        self,
        path_to_cohpcar: PathLike | None = "COHPCAR.lobster",
        obj_cohpcar: CompleteCohp | None = None,
        isites: list[int] | None = None,
        onlycation_isites: bool = True,
        only_bonds_to: list[str] | None = None,
        per_bond: bool = False,
        summed_spin_channels: bool = False,
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] = (-10, 6),
        integrated: bool = False,
    ) -> mpl.axes.Axes:
        """Plot summed COHPs or COBIs or COOPs.

        Please be careful in the spin polarized case (plots might overlap).

        Args:
            path_to_cohpcar (PathLike): Path to COHPCAR or COOPCAR or COBICAR.
            obj_cohpcar (CompleteCohp): CompleteCohp object
            isites (list[int]): Site IDs. If empty, all sites will be used to add the ICOHPs of the neighbors.
            onlycation_isites (bool): Only use cations, if isite is empty.
            only_bonds_to (list[str]): Only anions in this list will be considered.
            per_bond (bool): Whether to plot a normalization of the plotted COHP
                per number of bond (True), or the sum (False).
            xlim (tuple[float, float]): Limits of x values.
            ylim (tuple[float, float]): Limits of y values.
            integrated (bool): Whether to show integrated COHP instead of COHP.

        Returns:
            plt of the COHPs or COBIs or COOPs.
        """
        cp = CohpPlotter(are_cobis=self.are_cobis, are_coops=self.are_coops)

        plotlabel, summed_cohp = self.get_info_cohps_to_neighbors(
            path_to_cohpcar,
            obj_cohpcar,
            isites,
            only_bonds_to,
            onlycation_isites,
            per_bond,
            summed_spin_channels=summed_spin_channels,
        )

        cp.add_cohp(plotlabel, summed_cohp)
        ax = cp.get_plot(integrated=integrated)
        if xlim is not None:
            ax.set_xlim(xlim)

        if ylim is not None:
            ax.set_ylim(ylim)

        return ax

    def get_info_cohps_to_neighbors(
        self,
        path_to_cohpcar: PathLike | None = "COHPCAR.lobster",
        obj_cohpcar: CompleteCohp | None = None,
        isites: list[int] | None = None,
        only_bonds_to: list[str] | None = None,
        onlycation_isites: bool = True,
        per_bond: bool = True,
        summed_spin_channels: bool = False,
    ) -> tuple[str | None, CompleteCohp | None]:
        """Get the COHPs (COOPs or COBIs) as a summed Cohp object
        and a label from all sites mentioned in isites with neighbors.

        Args:
            path_to_cohpcar (PathLike): Path to COHPCAR/COOPCAR/COBICAR.
            obj_cohpcar (CompleteCohp): CompleteCohp object.
            isites (list[int]): The indexes of the sites.
            only_bonds_to (list[str]): Only show COHPs to selected element, e.g. ["O"].
            onlycation_isites (bool): If isites is None, only cation sites will be returned.
            per_bond (bool): Whether to normalize per bond.
            summed_spin_channels (bool): Whether to sum both spin channels.

        Returns:
            str: Label for COHP.
            CompleteCohp: Describe all COHPs/COOPs/COBIs of the sites
                as given by isites and the other arguments.
        """
        # TODO: add options for orbital-resolved COHPs
        _summed_icohps, _list_icohps, _number_bonds, labels, atoms, final_isites = self.get_info_icohps_to_neighbors(
            isites=isites, onlycation_isites=onlycation_isites
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            path = f"{tmp_dir}/POSCAR.vasp"

            self.structure.to(filename=path, fmt="poscar")

            if not hasattr(self, "completecohp"):
                if path_to_cohpcar is not None and obj_cohpcar is None:
                    self.completecohp = CompleteCohp.from_file(
                        fmt="LOBSTER",
                        filename=path_to_cohpcar,
                        structure_file=path,
                        are_coops=self.are_coops,
                        are_cobis=self.are_cobis,
                    )
                elif obj_cohpcar is not None:
                    self.completecohp = obj_cohpcar
                else:
                    raise ValueError("Please provide either path_to_cohpcar or obj_cohpcar")

        # Check that the number of bonds in ICOHPLIST and COHPCAR are identical
        # TODO: Further checks could be implemented
        if self.Icohpcollection is None:
            raise ValueError(f"{self.Icohpcollection=}")

        if len(self.Icohpcollection._list_atom1) != len(self.completecohp.bonds):
            raise ValueError("COHPCAR and ICOHPLIST do not fit together")

        is_spin_completecohp = Spin.down in self.completecohp.get_cohp_by_label("1").cohp
        if self.Icohpcollection.is_spin_polarized != is_spin_completecohp:
            raise ValueError("COHPCAR and ICOHPLIST do not fit together")

        if only_bonds_to is None:
            # Sort by anion type
            divisor = len(labels) if per_bond else 1

            plot_label = self._get_plot_label(atoms, per_bond)
            summed_cohp = self.completecohp.get_summed_cohp_by_label_list(
                label_list=labels,
                divisor=divisor,
                summed_spin_channels=summed_spin_channels,
            )

        else:
            # Labels of the COHPs that will be summed
            # Iterate through labels and atoms and check which bonds can be included
            new_labels = []
            new_atoms = []
            if final_isites is None:
                raise ValueError(f"{final_isites=}")

            for key, atompair, isite in zip(labels, atoms, final_isites, strict=True):
                present = False
                for atomtype in only_bonds_to:
                    # This is necessary to identify also bonds between the same elements correctly
                    if str(self.structure[isite].species.elements[0]) != atomtype:
                        if atomtype in {
                            self._split_string(atompair[0])[0],
                            self._split_string(atompair[1])[0],
                        }:
                            present = True
                    elif (
                        atomtype == self._split_string(atompair[0])[0]
                        and atomtype == self._split_string(atompair[1])[0]
                    ):
                        present = True

                if present:
                    new_labels.append(key)
                    new_atoms.append(atompair)
            if new_labels:
                divisor = len(new_labels) if per_bond else 1

                plot_label = self._get_plot_label(new_atoms, per_bond)
                summed_cohp = self.completecohp.get_summed_cohp_by_label_list(
                    label_list=new_labels,
                    divisor=divisor,
                    summed_spin_channels=summed_spin_channels,
                )
            else:
                plot_label = None

                summed_cohp = None

        return plot_label, summed_cohp

    def _get_plot_label(self, atoms: list[list[str]], per_bond: bool) -> str:
        """Count the types of bonds and append a label."""
        all_labels = []
        for atoms_names in atoms:
            new = [self._split_string(atoms_names[0])[0], self._split_string(atoms_names[1])[0]]
            new.sort()
            string_here = f"{new[0]}-{new[1]}"
            all_labels.append(string_here)

        counter = collections.Counter(all_labels)
        plotlabels = [f"{item} x {key}" for key, item in counter.items()]
        label = ", ".join(plotlabels)
        if per_bond:
            label += " (per bond)"
        return label

    def get_info_icohps_between_neighbors(
        self,
        isites: list[int] | None = None,
        onlycation_isites: bool = True,
    ) -> ICOHPNeighborsInfo:
        """Get interactions between neighbors of certain sites.

        Args:
            isites (list[int]): Site IDs. If is None, all sites will be used.
            onlycation_isites (bool): Only use cations, if isite is None.

        Returns:
            ICOHPNeighborsInfo
        """
        lowerlimit = self.lowerlimit
        upperlimit = self.upperlimit

        if self.valences is None and onlycation_isites:
            raise ValueError("No valences are provided")

        if isites is None:
            if onlycation_isites:
                if self.valences is None:
                    raise ValueError(f"{self.valences=}")

                isites = [idx for idx in range(len(self.structure)) if self.valences[idx] >= 0.0]
            else:
                isites = list(range(len(self.structure)))

        summed_icohps: float = 0.0
        list_icohps: list[float] = []
        number_bonds: int = 0
        labels: list[str] = []
        atoms: list[list[str]] = []
        if self.Icohpcollection is None:
            raise ValueError(f"{self.Icohpcollection=}")

        for isite in isites:
            for site1_idx, n_site in enumerate(self.list_neighsite[isite]):
                for site2_idx, n_site2 in enumerate(self.list_neighsite[isite]):
                    if site1_idx < site2_idx:
                        unitcell1 = self._determine_unit_cell(n_site)
                        unitcell2 = self._determine_unit_cell(n_site2)

                        index_n_site = self._get_original_site(self.structure, n_site)
                        index_n_site2 = self._get_original_site(self.structure, n_site2)

                        if index_n_site < index_n_site2:
                            translation = list(np.array(unitcell1) - np.array(unitcell2))
                        elif index_n_site2 < index_n_site:
                            translation = list(np.array(unitcell2) - np.array(unitcell1))
                        else:
                            translation = list(np.array(unitcell1) - np.array(unitcell2))

                        icohps = self._get_icohps(
                            icohpcollection=self.Icohpcollection,
                            site_idx=index_n_site,
                            lowerlimit=lowerlimit,
                            upperlimit=upperlimit,
                            only_bonds_to=self.only_bonds_to,
                        )

                        done = False
                        for icohp in icohps.values():
                            atomnr1 = self._get_atomnumber(icohp._atom1)
                            atomnr2 = self._get_atomnumber(icohp._atom2)
                            label = icohp._label

                            if (index_n_site == atomnr1 and index_n_site2 == atomnr2) or (
                                index_n_site == atomnr2 and index_n_site2 == atomnr1
                            ):
                                if atomnr1 != atomnr2:
                                    if np.all(np.asarray(translation) == np.asarray(icohp._translation)):
                                        summed_icohps += icohp.summed_icohp
                                        list_icohps.append(icohp.summed_icohp)
                                        number_bonds += 1
                                        labels.append(label)
                                        atoms.append(
                                            [
                                                self.Icohpcollection._list_atom1[int(label) - 1],
                                                self.Icohpcollection._list_atom2[int(label) - 1],
                                            ]
                                        )

                                elif not done:
                                    icohp_trans = -np.asarray(
                                        [icohp._translation[0], icohp._translation[1], icohp._translation[2]]
                                    )

                                    if (np.all(np.asarray(translation) == np.asarray(icohp._translation))) or (
                                        np.all(np.asarray(translation) == icohp_trans)
                                    ):
                                        summed_icohps += icohp.summed_icohp
                                        list_icohps.append(icohp.summed_icohp)
                                        number_bonds += 1
                                        labels.append(label)
                                        atoms.append(
                                            [
                                                self.Icohpcollection._list_atom1[int(label) - 1],
                                                self.Icohpcollection._list_atom2[int(label) - 1],
                                            ]
                                        )
                                        done = True

        return ICOHPNeighborsInfo(summed_icohps, list_icohps, number_bonds, labels, atoms, None)

    def _evaluate_ce(
        self,
        lowerlimit: float | None,
        upperlimit: float | None,
        only_bonds_to: list[str] | None = None,
        additional_condition: Literal[0, 1, 2, 3, 4, 5, 6] = 0,
        perc_strength_icohp: float = 0.15,
        adapt_extremum_to_add_cond: bool = False,
    ) -> None:
        """
        Args:
            lowerlimit (float): Lower limit which determines the ICOHPs
                that are considered for the determination of the neighbors.
            upperlimit (float): Upper limit which determines the ICOHPs
                that are considered for the determination of the neighbors.
            only_bonds_to (list[str]): Restrict the types of bonds that will be considered.
            additional_condition (int): Additional condition for the evaluation.
            perc_strength_icohp (float): Determine how strong the ICOHPs
                (percentage * strongest_ICOHP) will be that are still considered.
            adapt_extremum_to_add_cond (bool): Whether to recalculate the limit
                based on the bonding type and not on the overall extremum.
        """
        # Get extremum
        if lowerlimit is None and upperlimit is None:
            if self.Icohpcollection is None:
                raise ValueError(f"{self.Icohpcollection=}")

            limits = self._get_limit_from_extremum(
                self.Icohpcollection,
                percentage=perc_strength_icohp,
                adapt_extremum_to_add_cond=adapt_extremum_to_add_cond,
                additional_condition=additional_condition,
            )

            if limits is None:
                raise ValueError(f"{limits=}")
            lowerlimit, upperlimit = limits

        elif upperlimit is None or lowerlimit is None:
            raise ValueError("Please give two limits or leave them both at None")

        # Find environments based on ICOHP values
        list_icohps, list_keys, list_lengths, list_neighisite, list_neighsite, list_coords = self._find_environments(
            additional_condition, lowerlimit, upperlimit, only_bonds_to
        )

        self.list_icohps = list_icohps
        self.list_lengths = list_lengths
        self.list_keys = list_keys
        self.list_neighsite = list_neighsite
        self.list_neighisite = list_neighisite
        self.list_coords = list_coords

        # Make a structure graph
        # Make sure everything is relative to the given Structure and
        # not just the atoms in the unit cell
        if self.add_additional_data_sg:
            if self.bonding_list_1.icohpcollection is None:
                raise ValueError(f"{self.bonding_list_1.icohpcollection=}")
            if self.bonding_list_2.icohpcollection is None:
                raise ValueError(f"{self.bonding_list_2.icohpcollection=}")

            self.sg_list = [
                [
                    {
                        "site": neighbor,
                        "image": tuple(
                            int(round(idx))
                            for idx in (
                                neighbor.frac_coords
                                - self.structure[
                                    next(
                                        site_idx
                                        for site_idx, site in enumerate(self.structure)
                                        if neighbor.is_periodic_image(site)
                                    )
                                ].frac_coords
                            )
                        ),
                        "weight": 1,
                        # Here, the ICOBIs and ICOOPs are added based on the bond
                        # strength cutoff of the ICOHP
                        # More changes are necessary here if we use ICOBIs for cutoffs
                        "edge_properties": {
                            "ICOHP": self.list_icohps[neighbors_idx][nbr_idx],
                            "bond_length": self.list_lengths[neighbors_idx][nbr_idx],
                            "bond_label": self.list_keys[neighbors_idx][nbr_idx],
                            self.id_blist_sg1.upper(): self.bonding_list_1.icohpcollection.get_icohp_by_label(
                                self.list_keys[neighbors_idx][nbr_idx]
                            ),
                            self.id_blist_sg2.upper(): self.bonding_list_2.icohpcollection.get_icohp_by_label(
                                self.list_keys[neighbors_idx][nbr_idx]
                            ),
                        },
                        "site_index": next(
                            site_idx for site_idx, site in enumerate(self.structure) if neighbor.is_periodic_image(site)
                        ),
                    }
                    for nbr_idx, neighbor in enumerate(neighbors)
                ]
                for neighbors_idx, neighbors in enumerate(self.list_neighsite)
            ]
        else:
            self.sg_list = [
                [
                    {
                        "site": neighbor,
                        "image": tuple(
                            int(round(idx))
                            for idx in (
                                neighbor.frac_coords
                                - self.structure[
                                    next(
                                        site_idx
                                        for site_idx, site in enumerate(self.structure)
                                        if neighbor.is_periodic_image(site)
                                    )
                                ].frac_coords
                            )
                        ),
                        "weight": 1,
                        "edge_properties": {
                            "ICOHP": self.list_icohps[neighbors_idx][nbr_idx],
                            "bond_length": self.list_lengths[neighbors_idx][nbr_idx],
                            "bond_label": self.list_keys[neighbors_idx][nbr_idx],
                        },
                        "site_index": next(
                            site_idx for site_idx, site in enumerate(self.structure) if neighbor.is_periodic_image(site)
                        ),
                    }
                    for nbr_idx, neighbor in enumerate(neighbors)
                ]
                for neighbors_idx, neighbors in enumerate(self.list_neighsite)
            ]

    def _find_environments(
        self,
        additional_condition: Literal[0, 1, 2, 3, 4, 5, 6],
        lowerlimit: float,
        upperlimit: float,
        only_bonds_to: list[str] | None,
    ) -> tuple[
        list[list[IcohpValue]],
        list[list[str]],
        list[list[float]],
        list[list[int]],
        list[list[PeriodicNeighbor]],
        list[list[NDArray]],
    ]:
        """Find all relevant neighbors based on certain restrictions.

        Args:
            additional_condition (int): Additional condition.
            lowerlimit (float): Lower limit that ICOHPs are considered.
            upperlimit (float): Upper limit that ICOHPs are considered.
            only_bonds_to (list[str]): Only bonds to these elements will be considered.

        Returns:
            Tuple of ICOHPs, keys, lengths, neighisite, neighsite, coords.
        """
        list_icohps: list[list[IcohpValue]] = []
        list_keys: list[list[str]] = []
        list_lengths: list[list[float]] = []
        list_neighisite: list[list[int]] = []
        list_neighsite: list[list[PeriodicNeighbor]] = []
        list_coords: list[list[NDArray]] = []

        # Run over structure
        if self.Icohpcollection is None:
            raise ValueError(f"{self.Icohpcollection=}")

        for idx, site in enumerate(self.structure):
            icohps = self._get_icohps(
                icohpcollection=self.Icohpcollection,
                site_idx=idx,
                lowerlimit=lowerlimit,
                upperlimit=upperlimit,
                only_bonds_to=only_bonds_to,
            )

            additional_conds = self._find_relevant_atoms_additional_condition(idx, icohps, additional_condition)
            keys_from_ICOHPs, lengths_from_ICOHPs, neighbors_from_ICOHPs, selected_ICOHPs = additional_conds

            if len(neighbors_from_ICOHPs) > 0:
                centralsite = site

                neighbors_by_distance_start = self.structure.get_sites_in_sphere(
                    pt=centralsite.coords,
                    r=np.max(lengths_from_ICOHPs) + 0.5,
                    include_image=True,
                    include_index=True,
                )

                neighbors_by_distance = []
                list_distances = []
                index_here_list = []
                coords = []
                for neigh_new in sorted(neighbors_by_distance_start, key=lambda x: x[1]):
                    site_here = neigh_new[0].to_unit_cell()
                    index_here = neigh_new[2]
                    index_here_list.append(index_here)
                    cell_here = neigh_new[3]
                    new_coords = [
                        site_here.frac_coords[0] + float(cell_here[0]),
                        site_here.frac_coords[1] + float(cell_here[1]),
                        site_here.frac_coords[2] + float(cell_here[2]),
                    ]
                    coords.append(site_here.lattice.get_cartesian_coords(new_coords))

                    # new_site = PeriodicSite(
                    #     species=site_here.species_string,
                    #     coords=site_here.lattice.get_cartesian_coords(new_coords),
                    #     lattice=site_here.lattice,
                    #     to_unit_cell=False,
                    #     coords_are_cartesian=True,
                    # )
                    neighbors_by_distance.append(neigh_new[0])
                    list_distances.append(neigh_new[1])
                _list_neighsite = []
                _list_neighisite = []
                copied_neighbors_from_ICOHPs = copy.copy(neighbors_from_ICOHPs)
                copied_distances_from_ICOHPs = copy.copy(lengths_from_ICOHPs)
                _neigh_coords = []
                _neigh_frac_coords = []

                for neigh_idx, neigh in enumerate(neighbors_by_distance):
                    index_here2 = index_here_list[neigh_idx]

                    for dist_idx, dist in enumerate(copied_distances_from_ICOHPs):
                        if (
                            np.isclose(dist, list_distances[neigh_idx], rtol=1e-4)
                            and copied_neighbors_from_ICOHPs[dist_idx] == index_here2
                        ):
                            _list_neighsite.append(neigh)
                            _list_neighisite.append(index_here2)
                            _neigh_coords.append(coords[neigh_idx])
                            _neigh_frac_coords.append(neigh.frac_coords)
                            del copied_distances_from_ICOHPs[dist_idx]
                            del copied_neighbors_from_ICOHPs[dist_idx]
                            break

                list_neighisite.append(_list_neighisite)
                list_neighsite.append(_list_neighsite)
                list_lengths.append(lengths_from_ICOHPs)
                list_keys.append(keys_from_ICOHPs)
                list_coords.append(_neigh_coords)
                list_icohps.append(selected_ICOHPs)

            else:
                list_neighsite.append([])
                list_neighisite.append([])
                list_icohps.append([])
                list_lengths.append([])
                list_keys.append([])
                list_coords.append([])
        return (
            list_icohps,
            list_keys,
            list_lengths,
            list_neighisite,
            list_neighsite,
            list_coords,
        )

    def _find_relevant_atoms_additional_condition(
        self,
        site_idx: int,
        icohps: dict[str, IcohpValue],
        additional_condition: Literal[0, 1, 2, 3, 4, 5, 6],
    ) -> tuple[list[str], list[float], list[int], list[IcohpValue]]:
        """Find all relevant atoms that fulfill the additional condition.

        Args:
            site_idx (int): Site index in structure (start from 0).
            icohps (dict[str, IcohpValue]): ICOHP values.
            additional_condition (int): Additional condition.

        Returns:
            tuple: keys, lengths, neighbors from selected ICOHPs and selected ICOHPs.
        """
        keys_from_ICOHPs: list[str] = []
        lengths_from_ICOHPs: list[float] = []
        neighbors_from_ICOHPs: list[int] = []
        icohps_from_ICOHPs: list[IcohpValue] = []

        for key, icohp in icohps.items():
            atomnr1 = self._get_atomnumber(icohp._atom1)
            atomnr2 = self._get_atomnumber(icohp._atom2)

            # Check additional conditions
            val1 = val2 = None
            if self.valences is None:
                raise ValueError(f"{self.valences=}")
            if additional_condition in {1, 3, 5, 6}:
                val1 = self.valences[atomnr1]
                val2 = self.valences[atomnr2]

            # NO_ADDITIONAL_CONDITION
            if additional_condition == 0:
                if atomnr1 == site_idx:
                    neighbors_from_ICOHPs.append(atomnr2)
                    lengths_from_ICOHPs.append(icohp._length)
                    icohps_from_ICOHPs.append(icohp.summed_icohp)
                    keys_from_ICOHPs.append(key)
                elif atomnr2 == site_idx:
                    neighbors_from_ICOHPs.append(atomnr1)
                    lengths_from_ICOHPs.append(icohp._length)
                    icohps_from_ICOHPs.append(icohp.summed_icohp)
                    keys_from_ICOHPs.append(key)

            # ONLY_ANION_CATION_BONDS
            elif additional_condition == 1:
                if (val1 < 0.0 < val2) or (val2 < 0.0 < val1):  # type: ignore[operator]
                    if atomnr1 == site_idx:
                        neighbors_from_ICOHPs.append(atomnr2)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

                    elif atomnr2 == site_idx:
                        neighbors_from_ICOHPs.append(atomnr1)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

            # NO_ELEMENT_TO_SAME_ELEMENT_BONDS
            elif additional_condition == 2:
                if icohp._atom1.rstrip("0123456789") != icohp._atom2.rstrip("0123456789"):
                    if atomnr1 == site_idx:
                        neighbors_from_ICOHPs.append(atomnr2)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

                    elif atomnr2 == site_idx:
                        neighbors_from_ICOHPs.append(atomnr1)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

            # ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS
            elif additional_condition == 3:
                if ((val1 < 0.0 < val2) or (val2 < 0.0 < val1)) and icohp._atom1.rstrip(  # type: ignore[operator]
                    "0123456789"
                ) != icohp._atom2.rstrip("0123456789"):
                    if atomnr1 == site_idx:
                        neighbors_from_ICOHPs.append(atomnr2)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

                    elif atomnr2 == site_idx:
                        neighbors_from_ICOHPs.append(atomnr1)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

            # ONLY_ELEMENT_TO_OXYGEN_BONDS
            elif additional_condition == 4:
                if icohp._atom1.rstrip("0123456789") == "O" or icohp._atom2.rstrip("0123456789") == "O":
                    if atomnr1 == site_idx:
                        neighbors_from_ICOHPs.append(atomnr2)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

                    elif atomnr2 == site_idx:
                        neighbors_from_ICOHPs.append(atomnr1)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

            # DO_NOT_CONSIDER_ANION_CATION_BONDS
            elif additional_condition == 5:
                if (val1 > 0.0 and val2 > 0.0) or (val1 < 0.0 and val2 < 0.0):  # type: ignore[operator]
                    if atomnr1 == site_idx:
                        neighbors_from_ICOHPs.append(atomnr2)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

                    elif atomnr2 == site_idx:
                        neighbors_from_ICOHPs.append(atomnr1)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

            # ONLY_CATION_CATION_BONDS
            elif additional_condition == 6 and val1 > 0.0 and val2 > 0.0:  # type: ignore[operator]
                if atomnr1 == site_idx:
                    neighbors_from_ICOHPs.append(atomnr2)
                    lengths_from_ICOHPs.append(icohp._length)
                    icohps_from_ICOHPs.append(icohp.summed_icohp)
                    keys_from_ICOHPs.append(key)

                elif atomnr2 == site_idx:
                    neighbors_from_ICOHPs.append(atomnr1)
                    lengths_from_ICOHPs.append(icohp._length)
                    icohps_from_ICOHPs.append(icohp.summed_icohp)
                    keys_from_ICOHPs.append(key)

        return keys_from_ICOHPs, lengths_from_ICOHPs, neighbors_from_ICOHPs, icohps_from_ICOHPs

    @staticmethod
    def _get_icohps(
        icohpcollection: IcohpCollection,
        site_idx: int,
        lowerlimit: float | None,
        upperlimit: float | None,
        only_bonds_to: list[str] | None,
    ) -> dict[str, IcohpValue]:
        """Get ICOHP dict for certain site.

        Args:
            icohpcollection (IcohpCollection): IcohpCollection object.
            site_idx (int): Site index.
            lowerlimit (float): Lower limit that ICOHPs are considered.
            upperlimit (float): Upper limit that ICOHPs are considered.
            only_bonds_to (list[str]): Only bonds to these elements will be considered, e.g. ["O"].

        Returns:
            dict of IcohpValues. The keys correspond to the initial list_labels.
        """
        return icohpcollection.get_icohp_dict_of_site(
            site=site_idx,
            maxbondlength=6.0,
            minsummedicohp=lowerlimit,
            maxsummedicohp=upperlimit,
            only_bonds_to=only_bonds_to,
        )

    @staticmethod
    def _get_atomnumber(atomstring: str) -> int:
        """Get the index of the atom within the POSCAR (e.g., Return 0 for "Na1").

        Args:
            atomstring (str): Atom as str, such as "Na1".

        Returns:
            int: Index of the atom in the POSCAR.
        """
        return int(LobsterNeighbors._split_string(atomstring)[1]) - 1

    @staticmethod
    def _split_string(s) -> tuple[str, str]:
        """Split strings such as "Na1" into ["Na", "1"] and return "1".

        Args:
            s (str): String to split.
        """
        head = s.rstrip("0123456789")
        tail = s[len(head) :]
        return head, tail

    @staticmethod
    def _determine_unit_cell(site: PeriodicSite) -> list[int]:
        """Determine the unit cell based on the site.

        Args:
            site (PeriodicSite): The site.
        """
        unitcell = []
        for coord in site.frac_coords:
            value = math.floor(round(coord, 4))
            unitcell.append(value)

        return unitcell

    def _adapt_extremum_to_add_cond(
        self,
        list_icohps: list[float],
        percentage: float,
    ) -> float:
        """Get the extremum from the given ICOHPs or ICOOPs or ICOBIs.

        Args:
            list_icohps (list): ICOHPs or ICOOPs or ICOBIs.
            percentage (float): The percentage to scale extremum.

        Returns:
            float: Min value of ICOHPs, or max value of ICOOPs/ICOBIs.
        """

        which_extr = min if not self.are_coops and not self.are_cobis else max
        return which_extr(list_icohps) * percentage

    def _get_limit_from_extremum(
        self,
        icohpcollection: IcohpCollection,
        percentage: float = 0.15,
        adapt_extremum_to_add_cond: bool = False,
        additional_condition: Literal[0, 1, 2, 3, 4, 5, 6] = 0,
    ) -> tuple[float, float] | None:
        """Get range for the ICOHP values from an IcohpCollection.

        Currently only work for ICOHPs.

        Args:
            icohpcollection (IcohpCollection): IcohpCollection object.
            percentage (float): Determine which ICOHPs/ICOOP/ICOBI will be considered.
            adapt_extremum_to_add_cond (bool): Whether the extrumum be adapted to
                the additional condition.
            additional_condition (int): Additional condition to determine which bonds to include.

        Returns:
            tuple[float, float]: [-inf, min(strongest_icohp*0.15, -noise_cutoff)]
                or [max(strongest_icohp*0.15, noise_cutoff), inf].
        """
        extremum_based = None

        if self.valences is None:
            raise ValueError(f"{self.valences=}")

        if not adapt_extremum_to_add_cond or additional_condition == 0:
            extremum_based = icohpcollection.extremum_icohpvalue(summed_spin_channels=True) * percentage

        elif additional_condition == 1:
            # ONLY_ANION_CATION_BONDS
            list_icohps = []
            for value in icohpcollection._icohplist.values():
                atomnr1 = type(self)._get_atomnumber(value._atom1)
                atomnr2 = type(self)._get_atomnumber(value._atom2)

                val1 = self.valences[atomnr1]
                val2 = self.valences[atomnr2]
                if (val1 < 0.0 < val2) or (val2 < 0.0 < val1):
                    list_icohps.append(value.summed_icohp)

            extremum_based = self._adapt_extremum_to_add_cond(list_icohps, percentage)

        elif additional_condition == 2:
            # NO_ELEMENT_TO_SAME_ELEMENT_BONDS
            list_icohps = []
            for value in icohpcollection._icohplist.values():
                if value._atom1.rstrip("0123456789") != value._atom2.rstrip("0123456789"):
                    list_icohps.append(value.summed_icohp)

            extremum_based = self._adapt_extremum_to_add_cond(list_icohps, percentage)

        elif additional_condition == 3:
            # ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS
            list_icohps = []
            for value in icohpcollection._icohplist.values():
                atomnr1 = type(self)._get_atomnumber(value._atom1)
                atomnr2 = type(self)._get_atomnumber(value._atom2)
                val1 = self.valences[atomnr1]
                val2 = self.valences[atomnr2]

                if ((val1 < 0.0 < val2) or (val2 < 0.0 < val1)) and value._atom1.rstrip(
                    "0123456789"
                ) != value._atom2.rstrip("0123456789"):
                    list_icohps.append(value.summed_icohp)

            extremum_based = self._adapt_extremum_to_add_cond(list_icohps, percentage)

        elif additional_condition == 4:
            # ONLY_ELEMENT_TO_OXYGEN_BONDS
            list_icohps = []
            for value in icohpcollection._icohplist.values():
                if value._atom1.rstrip("0123456789") == "O" or value._atom2.rstrip("0123456789") == "O":
                    list_icohps.append(value.summed_icohp)

            extremum_based = self._adapt_extremum_to_add_cond(list_icohps, percentage)

        elif additional_condition == 5:
            # DO_NOT_CONSIDER_ANION_CATION_BONDS
            list_icohps = []
            for value in icohpcollection._icohplist.values():
                atomnr1 = type(self)._get_atomnumber(value._atom1)
                atomnr2 = type(self)._get_atomnumber(value._atom2)
                val1 = self.valences[atomnr1]
                val2 = self.valences[atomnr2]

                if (val1 > 0.0 and val2 > 0.0) or (val1 < 0.0 and val2 < 0.0):
                    list_icohps.append(value.summed_icohp)

            extremum_based = self._adapt_extremum_to_add_cond(list_icohps, percentage)

        elif additional_condition == 6:
            # ONLY_CATION_CATION_BONDS
            list_icohps = []
            for value in icohpcollection._icohplist.values():
                atomnr1 = type(self)._get_atomnumber(value._atom1)
                atomnr2 = type(self)._get_atomnumber(value._atom2)
                val1 = self.valences[atomnr1]
                val2 = self.valences[atomnr2]

                if val1 > 0.0 and val2 > 0.0:
                    list_icohps.append(value.summed_icohp)

            extremum_based = self._adapt_extremum_to_add_cond(list_icohps, percentage)

        if not self.are_coops and not self.are_cobis:
            max_here = min(extremum_based, -self.noise_cutoff) if self.noise_cutoff is not None else extremum_based
            return -float("inf"), max_here

        if self.are_coops or self.are_cobis:
            min_here = max(extremum_based, self.noise_cutoff) if self.noise_cutoff is not None else extremum_based
            return min_here, float("inf")

        return None


class LobsterLightStructureEnvironments(LightStructureEnvironments):
    """Store LightStructureEnvironments based on LOBSTER outputs."""

    @classmethod
    def from_Lobster(
        cls,
        list_ce_symbol: list[str],
        list_csm: list[float],
        list_permutation: list,
        list_neighsite: list[PeriodicSite],
        list_neighisite: list[list[int]],
        structure: Structure,
        valences: list[float] | None = None,
    ) -> Self:
        """Set up a LightStructureEnvironments from LOBSTER.

        Args:
            list_ce_symbol (list[str]): Coordination environments symbols.
            list_csm (list[float]): Continuous symmetry measures.
            list_permutation (list): Permutations.
            list_neighsite (list[PeriodicSite]): Neighboring sites.
            list_neighisite (list[list[int]]): Neighboring sites indexes.
            structure (Structure): Structure object.
            valences (list[float]): Valences.

        Returns:
            LobsterLightStructureEnvironments
        """
        strategy = None
        valences_origin = "user-defined"
        coordination_environments = []
        all_nbs_sites = []
        all_nbs_sites_indices = []
        neighbors_sets = []
        counter = 0

        for site_idx in range(len(structure)):
            # Coordination environment
            if list_ce_symbol is not None:
                ce_dict = {
                    "ce_symbol": list_ce_symbol[site_idx],
                    "ce_fraction": 1.0,
                    "csm": list_csm[site_idx],
                    "permutation": list_permutation[site_idx],
                }
            else:
                ce_dict = None

            if list_neighisite[site_idx] is not None:
                all_nbs_sites_indices_here = []
                for neigh_site_idx, neigh_site in enumerate(list_neighsite[site_idx]):
                    diff = neigh_site.frac_coords - structure[list_neighisite[site_idx][neigh_site_idx]].frac_coords
                    round_diff = np.round(diff)
                    if not np.allclose(diff, round_diff):
                        raise ValueError(
                            "Weird, differences between one site in a periodic image cell is not integer ..."
                        )
                    nb_image_cell = np.array(round_diff, int)

                    all_nbs_sites_indices_here.append(counter)

                    neighbor = {
                        "site": neigh_site,
                        "index": list_neighisite[site_idx][neigh_site_idx],
                        "image_cell": nb_image_cell,
                    }
                    all_nbs_sites.append(neighbor)
                    counter += 1

                all_nbs_sites_indices.append(all_nbs_sites_indices_here)

            else:
                all_nbs_sites.append({"site": None, "index": None, "image_cell": None})
                all_nbs_sites_indices.append([])

            if list_neighisite[site_idx] is not None:
                nb_set = cls.NeighborsSet(
                    structure=structure,
                    isite=site_idx,
                    all_nbs_sites=all_nbs_sites,
                    all_nbs_sites_indices=all_nbs_sites_indices[site_idx],
                )

            else:
                nb_set = cls.NeighborsSet(
                    structure=structure,
                    isite=site_idx,
                    all_nbs_sites=[],
                    all_nbs_sites_indices=[],
                )

            coordination_environments.append([ce_dict])
            neighbors_sets.append([nb_set])

        return cls(
            strategy=strategy,
            coordination_environments=coordination_environments,
            all_nbs_sites=all_nbs_sites,
            neighbors_sets=neighbors_sets,
            structure=structure,
            valences=valences,
            valences_origin=valences_origin,
        )

    @property
    def uniquely_determines_coordination_environments(self) -> Literal[True]:
        """Whether the coordination environments are uniquely determined."""
        return True

    def as_dict(self) -> dict[str, Any]:
        """Bson-serializable dict representation of the object.

        Returns:
            Bson-serializable dict representation.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "strategy": self.strategy,
            "structure": self.structure.as_dict(),
            "coordination_environments": self.coordination_environments,
            "all_nbs_sites": [
                {
                    "site": nb_site["site"].as_dict(),
                    "index": nb_site["index"],
                    "image_cell": [int(ii) for ii in nb_site["image_cell"]],
                }
                for nb_site in self._all_nbs_sites
            ],
            "neighbors_sets": [
                [nb_set.as_dict() for nb_set in site_nb_sets] or None for site_nb_sets in self.neighbors_sets
            ],
            "valences": self.valences,
        }


class ICOHPNeighborsInfo(NamedTuple):
    """Tuple to record information on relevant bonds.

    Args:
        total_icohp (float): Sum of ICOHP values of neighbors to the selected
            sites (given by the index in structure).
        list_icohps (list): Summed ICOHP values for all identified interactions with neighbors.
        n_bonds (int): Number of identified bonds to the selected sites.
        labels (list[str]): Labels (from ICOHPLIST) for all identified bonds.
        atoms (list[list[str]]): Lists describing the species present (from ICOHPLIST)
            in the identified interactions , e.g. ["Ag3", "O5"].
        central_isites (list[int]): The central site indexes for each identified interaction.
    """

    total_icohp: float
    list_icohps: list[float]
    n_bonds: int
    labels: list[str]
    atoms: list[list[str]]
    central_isites: list[int] | None
