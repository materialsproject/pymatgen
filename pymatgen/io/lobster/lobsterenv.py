# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to perform analyses of
the local environments (e.g., finding near neighbors)
of single sites in molecules and structures based on
bonding analysis with Lobster.
"""

import collections
import copy
import math
import os

import numpy as np
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.analysis.local_env import NearNeighbors
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import CohpPlotter
from pymatgen.io.lobster import Charge, Icohplist

__author__ = "Janine George"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "1.0"
__maintainer__ = "J. George"
__email__ = "janinegeorge.ulfen@gmail.com"
__status__ = "Production"
__date__ = "February 2, 2021"


class LobsterNeighbors(NearNeighbors):
    """
    This class combines capabilities from LocalEnv and ChemEnv to determine coordination environments based on
    bonding analysis
    """

    def __init__(
        self,
        are_coops=False,
        filename_ICOHP=None,
        valences=None,
        limits=None,
        structure=None,
        additional_condition=0,
        only_bonds_to=None,
        perc_strength_ICOHP=0.15,
        valences_from_charges=False,
        filename_CHARGE=None,
    ):
        """

        Args:
            are_coops: (Bool) if True, the file is a ICOOPLIST.lobster and not a ICOHPLIST.lobster; only tested for
            ICOHPLIST.lobster so far
            filename_ICOHP: (str) Path to ICOOPLIST.lobster
            valences: (list of integers/floats) gives valence/charge for each element
            limits: limit to decide which ICOHPs should be considered
            structure: (Structure Object) typically constructed by: Structure.from_file("POSCAR") (Structure object
            from pymatgen.core.structure)
            additional_condition:   Additional condition that decides which kind of bonds will be considered
                                    NO_ADDITIONAL_CONDITION = 0
                                    ONLY_ANION_CATION_BONDS = 1
                                    NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
                                    ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
                                    ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
                                    DO_NOT_CONSIDER_ANION_CATION_BONDS=5
                                    ONLY_CATION_CATION_BONDS=6
            only_bonds_to: (list of str) will only consider bonds to certain elements (e.g. ["O"] for oxygen)
            perc_strength_ICOHP: if no limits are given, this will decide which icohps will still be considered (
            relative to
            the strongest ICOHP)
            valences_from_charges: if True and path to CHARGE.lobster is provided, will use Lobster charges (
            Mulliken) instead of valences
            filename_CHARGE: (str) Path to Charge.lobster
        """

        self.ICOHP = Icohplist(are_coops=are_coops, filename=filename_ICOHP)
        self.Icohpcollection = self.ICOHP.icohpcollection
        self.structure = structure
        self.limits = limits
        self.only_bonds_to = only_bonds_to
        self.are_coops = are_coops
        if are_coops:
            raise ValueError("Algorithm only works correctly for ICOHPLIST.lobster")

        # will check if the additional condition is correctly delivered
        if additional_condition in range(0, 7):
            self.additional_condition = additional_condition
        else:
            raise ValueError("No correct additional condition")

        # will read in valences, will prefer manual setting of valences
        if valences is None:
            if valences_from_charges and filename_CHARGE is not None:
                chg = Charge(filename=filename_CHARGE)
                self.valences = chg.Mulliken
            else:
                bv_analyzer = BVAnalyzer()
                try:
                    self.valences = bv_analyzer.get_valences(structure=self.structure)
                except ValueError:
                    self.valences = None
                    if additional_condition in [1, 3, 5, 6]:
                        print("Valences cannot be assigned, additional_conditions 1 and 3 and 5 and 6 will not work")
        else:
            self.valences = valences

        if limits is None:
            self.lowerlimit = None
            self.upperlimit = None

        else:
            self.lowerlimit = limits[0]
            self.upperlimit = limits[1]

        # will evaluate coordination environments
        self._evaluate_ce(
            lowerlimit=self.lowerlimit,
            upperlimit=self.upperlimit,
            only_bonds_to=only_bonds_to,
            additional_condition=self.additional_condition,
            perc_strength_ICOHP=perc_strength_ICOHP,
        )

    @property
    def structures_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Structure
        objects?
        """
        return True

    @property
    def molecules_allowed(self):
        """
        Boolean property: can this NearNeighbors class be used with Molecule
        objects?
        """
        return False

    def get_nn_info(self, structure, n, use_weights=False):
        """
        Get coordination number, CN, of site with index n in structure.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine CN.
            use_weights (boolean): flag indicating whether (True)
                to use weights for computing the coordination number
                or not (False, default: each coordinated site has equal
                weight).
                True is not implemented for LobsterNeighbors
        Returns:
            cn (integer or float): coordination number.
        """
        if use_weights:
            raise ValueError("LobsterEnv cannot use weights")
        if len(structure) != len(self.structure):
            raise ValueError("The wrong structure was provided")
        return self.sg_list[n]

    def get_light_structure_environment(self, only_cation_environments=False):
        """
        will return a LobsterLightStructureEnvironments object
        if the structure only contains coordination environments smaller 13
        Args:
            only_cation_environments: only data for cations will be returned

        Returns: LobsterLightStructureEnvironments Object

        """

        lgf = LocalGeometryFinder()
        lgf.setup_structure(structure=self.structure)
        list_ce_symbols = []
        list_csm = []
        list_permut = []
        for ival, _neigh_coords in enumerate(self.list_coords):

            if (len(_neigh_coords)) > 13:
                raise ValueError("Environment cannot be determined. Number of neighbors is larger than 13.")
            lgf.setup_local_geometry(isite=ival, coords=_neigh_coords, optimization=2)
            cncgsm = lgf.get_coordination_symmetry_measures(optimization=2)

            list_ce_symbols.append(min(cncgsm.items(), key=lambda t: t[1]["csm_wcs_ctwcc"])[0])
            list_csm.append(min(cncgsm.items(), key=lambda t: t[1]["csm_wcs_ctwcc"])[1]["csm_wcs_ctwcc"])
            list_permut.append(min(cncgsm.items(), key=lambda t: t[1]["csm_wcs_ctwcc"])[1]["indices"])

        if not only_cation_environments:
            lse = LobsterLightStructureEnvironments.from_Lobster(
                list_ce_symbol=list_ce_symbols,
                list_csm=list_csm,
                list_permutation=list_permut,
                list_neighsite=self.list_neighsite,
                list_neighisite=self.list_neighisite,
                structure=self.structure,
                valences=self.valences,
            )
        else:
            new_list_ce_symbols = []
            new_list_csm = []
            new_list_permut = []
            new_list_neighsite = []
            new_list_neighisite = []

            for ival, val in enumerate(self.valences):

                if val >= 0.0:

                    new_list_ce_symbols.append(list_ce_symbols[ival])
                    new_list_csm.append(list_csm[ival])
                    new_list_permut.append(list_permut[ival])
                    new_list_neighisite.append(self.list_neighisite[ival])
                    new_list_neighsite.append(self.list_neighsite[ival])
                else:
                    new_list_ce_symbols.append(None)
                    new_list_csm.append(None)
                    new_list_permut.append([])
                    new_list_neighisite.append([])
                    new_list_neighsite.append([])

            lse = LobsterLightStructureEnvironments.from_Lobster(
                list_ce_symbol=new_list_ce_symbols,
                list_csm=new_list_csm,
                list_permutation=new_list_permut,
                list_neighsite=new_list_neighsite,
                list_neighisite=new_list_neighisite,
                structure=self.structure,
                valences=self.valences,
            )

        return lse

    def get_info_icohps_to_neighbors(self, isites=[], onlycation_isites=True):
        """
        this method will return information of cohps of neighbors
        Args:
            isites: list of site ids, if isite==[], all isites will be used to add the icohps of the neighbors
            onlycation_isites: will only use cations, if isite==[]


        Returns:
            sum of icohps of neighbors to certain sites [given by the id in structure], number of bonds to this site,
            labels (from ICOHPLIST) for
            these bonds
            [the latter is useful for plotting summed COHP plots]
        """

        if self.valences is None and onlycation_isites:
            raise ValueError("No valences are provided")
        if isites == []:
            if onlycation_isites:
                isites = [i for i in range(len(self.structure)) if self.valences[i] >= 0.0]
            else:
                isites = list(range(len(self.structure)))

        summed_icohps = 0.0
        list_icohps = []
        number_bonds = 0
        labels = []
        atoms = []
        for ival, site in enumerate(self.structure):
            if ival in isites:
                for keys, icohpsum in zip(self.list_keys[ival], self.list_icohps[ival]):
                    summed_icohps += icohpsum
                    list_icohps.append(icohpsum)
                    labels.append(keys)
                    atoms.append(
                        [
                            self.Icohpcollection._list_atom1[int(keys) - 1],
                            self.Icohpcollection._list_atom2[int(keys) - 1],
                        ]
                    )
                    number_bonds += 1

        return summed_icohps, list_icohps, number_bonds, labels, atoms

    def plot_cohps_of_neighbors(
        self,
        path_to_COHPCAR="COHPCAR.lobster",
        isites=[],
        onlycation_isites=True,
        only_bonds_to=None,
        per_bond=False,
        summed_spin_channels=False,
        xlim=None,
        ylim=[-10, 6],
        integrated=False,
    ):

        """
        will plot summed cohps (please be careful in the spin polarized case (plots might overlap (exactly!))
        Args:
            isites: list of site ids, if isite==[], all isites will be used to add the icohps of the neighbors
            onlycation_isites: bool, will only use cations, if isite==[]
            only_bonds_to: list of str, only anions in this list will be considered
            per_bond: bool, will lead to a normalization of the plotted COHP per number of bond if True,
            otherwise the sum
            will be plotted
            xlim: list of float, limits of x values
            ylim: list of float, limits of y values
            integrated: bool, if true will show integrated cohp instead of cohp

        Returns:
            plt of the cohps

        """

        # include COHPPlotter and plot a sum of these COHPs
        # might include option to add Spin channels
        # implement only_bonds_to
        cp = CohpPlotter()

        plotlabel, summed_cohp = self.get_info_cohps_to_neighbors(
            path_to_COHPCAR,
            isites,
            only_bonds_to,
            onlycation_isites,
            per_bond,
            summed_spin_channels=summed_spin_channels,
        )

        cp.add_cohp(plotlabel, summed_cohp)
        plot = cp.get_plot(integrated=integrated)
        if xlim is not None:
            plot.xlim(xlim)

        if ylim is not None:
            plot.ylim(ylim)

        return plot

    def get_info_cohps_to_neighbors(
        self,
        path_to_COHPCAR="COHPCAR.lobster",
        isites=[],
        only_bonds_to=None,
        onlycation_isites=True,
        per_bond=True,
        summed_spin_channels=False,
    ):
        """
        will return info about the cohps from all sites mentioned in isites with neighbors
        Args:
            path_to_COHPCAR: str, path to COHPCAR
            isites: list of int that indicate the number of the site
            only_bonds_to: list of str, e.g. ["O"] to only show cohps of anything to oxygen
            onlycation_isites: if isites=[], only cation sites will be returned
            per_bond: will normalize per bond
            summed_spin_channels: will sum all spin channels

        Returns: label for cohp (str), CompleteCohp object which describes all cohps of the sites as given by isites
        and the other parameters

        """
        # TODO: add options for orbital-resolved cohps
        summed_icohps, list_icohps, number_bonds, labels, atoms = self.get_info_icohps_to_neighbors(
            isites=isites, onlycation_isites=onlycation_isites
        )
        import tempfile

        with tempfile.TemporaryDirectory() as t:
            path = os.path.join(t, "POSCAR.vasp")

            self.structure.to(filename=path, fmt="POSCAR")

            completecohp = CompleteCohp.from_file(fmt="LOBSTER", filename=path_to_COHPCAR, structure_file=path)

        # will check that the number of bonds in ICOHPLIST and COHPCAR are identical
        # further checks could be implemented
        if len(self.Icohpcollection._list_atom1) != len(completecohp.bonds.keys()):
            raise ValueError("COHPCAR and ICOHPLIST do not fit together")
        is_spin_completecohp = Spin.down in completecohp.get_cohp_by_label("1").cohp
        if self.Icohpcollection.is_spin_polarized != is_spin_completecohp:
            raise ValueError("COHPCAR and ICOHPLIST do not fit together")

        if only_bonds_to is None:
            # sort by anion type
            if per_bond:
                divisor = len(labels)
            else:
                divisor = 1

            plotlabel = self._get_plot_label(atoms, per_bond)
            summed_cohp = completecohp.get_summed_cohp_by_label_list(
                label_list=labels, divisor=divisor, summed_spin_channels=summed_spin_channels
            )

        else:
            # labels of the COHPs that will be summed!
            # iterate through labels and atoms and check which bonds can be included
            new_labels = []
            new_atoms = []
            for label, atompair in zip(labels, atoms):
                # durchlaufe only_bonds_to=[] und sage ja, falls eines der Labels in atompair ist, dann speichere
                # new_label
                present = False
                for atomtype in only_bonds_to:
                    if atomtype in (self._split_string(atompair[0])[0], self._split_string(atompair[1])[0]):
                        present = True
                if present:
                    new_labels.append(label)
                    new_atoms.append(atompair)

            if len(new_labels) > 0:
                if per_bond:
                    divisor = len(new_labels)
                else:
                    divisor = 1

                plotlabel = self._get_plot_label(new_atoms, per_bond)
                summed_cohp = completecohp.get_summed_cohp_by_label_list(
                    label_list=new_labels, divisor=divisor, summed_spin_channels=summed_spin_channels
                )
            else:
                plotlabel = None
                summed_cohp = None

        return plotlabel, summed_cohp

    def _get_plot_label(self, atoms, per_bond):
        # count the types of bonds and append a label:
        all_labels = []
        for atomsnames in atoms:
            new = [self._split_string(atomsnames[0])[0], self._split_string(atomsnames[1])[0]]
            new.sort()
            # print(new2)
            string_here = new[0] + "-" + new[1]
            all_labels.append(string_here)
        count = collections.Counter(all_labels)
        plotlabels = []
        for key, item in count.items():
            plotlabels.append(str(item) + " x " + str(key))
        plotlabel = ", ".join(plotlabels)
        if per_bond:
            plotlabel = plotlabel + " (per bond)"
        return plotlabel

    def get_info_icohps_between_neighbors(self, isites=[], onlycation_isites=True):

        """
        will return infos about interactions between neighbors of a certain atom
        Args:
            isites: list of site ids, if isite==[], all isites will be used
            onlycation_isites: will only use cations, if isite==[]

        Returns:

        """

        lowerlimit = self.lowerlimit
        upperlimit = self.upperlimit

        if self.valences is None and onlycation_isites:
            raise ValueError("No valences are provided")
        if isites == []:
            if onlycation_isites:
                isites = [i for i in range(len(self.structure)) if self.valences[i] >= 0.0]
            else:
                isites = list(range(len(self.structure)))

        summed_icohps = 0.0
        list_icohps = []
        number_bonds = 0
        label_list = []
        atoms_list = []
        for iisite, isite in enumerate(isites):
            for in_site, n_site in enumerate(self.list_neighsite[isite]):
                for in_site2, n_site2 in enumerate(self.list_neighsite[isite]):
                    if in_site < in_site2:
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
                            isite=index_n_site,
                            lowerlimit=lowerlimit,
                            upperlimit=upperlimit,
                            only_bonds_to=self.only_bonds_to,
                        )

                        done = False
                        for key, icohp in icohps.items():

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
                                        label_list.append(label)
                                        atoms_list.append(
                                            [
                                                self.Icohpcollection._list_atom1[int(label) - 1],
                                                self.Icohpcollection._list_atom2[int(label) - 1],
                                            ]
                                        )

                                else:
                                    if not done:
                                        if (np.all(np.asarray(translation) == np.asarray(icohp._translation))) or (
                                            np.all(
                                                np.asarray(translation)
                                                == np.asarray(
                                                    [
                                                        -icohp._translation[0],
                                                        -icohp._translation[1],
                                                        -icohp._translation[2],
                                                    ]
                                                )
                                            )
                                        ):
                                            summed_icohps += icohp.summed_icohp
                                            list_icohps.append(icohp.summed_icohp)
                                            number_bonds += 1
                                            label_list.append(label)
                                            atoms_list.append(
                                                [
                                                    self.Icohpcollection._list_atom1[int(label) - 1],
                                                    self.Icohpcollection._list_atom2[int(label) - 1],
                                                ]
                                            )
                                            done = True

        return summed_icohps, list_icohps, number_bonds, label_list, atoms_list

    def _evaluate_ce(
        self, lowerlimit, upperlimit, only_bonds_to=None, additional_condition=0, perc_strength_ICOHP=0.15
    ):
        """

        Args:
            lowerlimit: lower limit which determines the ICOHPs that are considered for the determination of the
            neighbors
            upperlimit: upper limit which determines the ICOHPs that are considered for the determination of the
            neighbors
            only_bonds_to: restricts the types of bonds that will be considered
            additional_condition: Additional condition for the evaluation
            perc_strength_ICOHP: will be used to determine how strong the ICOHPs (percentage*strongest ICOHP) will be
            that are still considered for the evalulation

        Returns:

        """
        # get extremum
        if lowerlimit is None and upperlimit is None:
            lowerlimit, upperlimit = self._get_limit_from_extremum(self.Icohpcollection, percentage=perc_strength_ICOHP)
        elif lowerlimit is None and upperlimit is not None:
            raise ValueError("Please give two limits or leave them both at None")
        elif upperlimit is None and lowerlimit is not None:
            raise ValueError("Please give two limits or leave them both at None")

        # find environments based on ICOHP values
        list_icohps, list_keys, list_lengths, list_neighisite, list_neighsite, list_coords = self._find_environments(
            additional_condition, lowerlimit, upperlimit, only_bonds_to
        )

        self.list_icohps = list_icohps
        self.list_lengths = list_lengths
        self.list_keys = list_keys
        self.list_neighsite = list_neighsite
        self.list_neighisite = list_neighisite
        self.list_coords = list_coords

        # make a structure graph
        # make sure everything is relative to the given Structure and not just the atoms in the unit cell
        self.sg_list = [
            [
                {
                    "site": neighbor,
                    "image": tuple(
                        int(round(i))
                        for i in (
                            neighbor.frac_coords
                            - self.structure[
                                [
                                    isite
                                    for isite, site in enumerate(self.structure)
                                    if neighbor.is_periodic_image(site)
                                ][0]
                            ].frac_coords
                        )
                    ),
                    "weight": 1,
                    "site_index": [
                        isite for isite, site in enumerate(self.structure) if neighbor.is_periodic_image(site)
                    ][0],
                }
                for neighbor in neighbors
            ]
            for neighbors in self.list_neighsite
        ]

    def _find_environments(self, additional_condition, lowerlimit, upperlimit, only_bonds_to):
        """
        will find all relevant neighbors based on certain restrictions
        Args:
            additional_condition (int): additional condition (see above)
            lowerlimit (float): lower limit that tells you which ICOHPs are considered
            upperlimit (float): upper limit that tells you which ICOHPs are considerd
            only_bonds_to (list): list of str, e.g. ["O"] that will ensure that only bonds to "O" will be considered

        Returns:

        """
        # run over structure
        list_neighsite = []
        list_neighisite = []
        list_coords = []
        list_icohps = []
        list_lengths = []
        list_keys = []
        for isite, site in enumerate(self.structure):

            icohps = self._get_icohps(
                icohpcollection=self.Icohpcollection,
                isite=isite,
                lowerlimit=lowerlimit,
                upperlimit=upperlimit,
                only_bonds_to=only_bonds_to,
            )

            (
                keys_from_ICOHPs,
                lengths_from_ICOHPs,
                neighbors_from_ICOHPs,
                selected_ICOHPs,
            ) = self._find_relevant_atoms_additional_condition(isite, icohps, additional_condition)

            if len(neighbors_from_ICOHPs) > 0:
                centralsite = self.structure.sites[isite]

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
                    newcoords = [
                        site_here.frac_coords[0] + np.float(cell_here[0]),
                        site_here.frac_coords[1] + np.float(cell_here[1]),
                        site_here.frac_coords[2] + np.float(cell_here[2]),
                    ]
                    coords.append(site_here.lattice.get_cartesian_coords(newcoords))

                    # new_site = PeriodicSite(species=site_here.species_string,
                    #                         coords=site_here.lattice.get_cartesian_coords(newcoords),
                    #                         lattice=site_here.lattice, to_unit_cell=False, coords_are_cartesian=True)
                    neighbors_by_distance.append(neigh_new[0])
                    list_distances.append(neigh_new[1])
                _list_neighsite = []
                _list_neighisite = []
                copied_neighbors_from_ICOHPs = copy.copy(neighbors_from_ICOHPs)
                copied_distances_from_ICOHPs = copy.copy(lengths_from_ICOHPs)
                _neigh_coords = []
                _neigh_frac_coords = []

                for ineigh, neigh in enumerate(neighbors_by_distance):
                    index_here2 = index_here_list[ineigh]

                    for idist, dist in enumerate(copied_distances_from_ICOHPs):
                        if (
                            np.isclose(dist, list_distances[ineigh], rtol=1e-4)
                            and copied_neighbors_from_ICOHPs[idist] == index_here2
                        ):
                            _list_neighsite.append(neigh)
                            _list_neighisite.append(index_here2)
                            _neigh_coords.append(coords[ineigh])
                            _neigh_frac_coords.append(neigh.frac_coords)
                            del copied_distances_from_ICOHPs[idist]
                            del copied_neighbors_from_ICOHPs[idist]
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
        return list_icohps, list_keys, list_lengths, list_neighisite, list_neighsite, list_coords

    def _find_relevant_atoms_additional_condition(self, isite, icohps, additional_condition):
        """
        will find all relevant atoms that fulfill the additional_conditions
        Args:
            isite: number of site in structure (starts with 0)
            icohps: icohps
            additional_condition (int): additonal condition

        Returns:

        """
        neighbors_from_ICOHPs = []
        lengths_from_ICOHPs = []
        icohps_from_ICOHPs = []
        keys_from_ICOHPs = []

        for key, icohp in icohps.items():
            atomnr1 = self._get_atomnumber(icohp._atom1)
            atomnr2 = self._get_atomnumber(icohp._atom2)

            # test additional conditions
            if additional_condition in (1, 3, 5, 6):
                val1 = self.valences[atomnr1]
                val2 = self.valences[atomnr2]

            if additional_condition == 0:
                # NO_ADDITIONAL_CONDITION
                if atomnr1 == isite:
                    neighbors_from_ICOHPs.append(atomnr2)
                    lengths_from_ICOHPs.append(icohp._length)
                    icohps_from_ICOHPs.append(icohp.summed_icohp)
                    keys_from_ICOHPs.append(key)
                elif atomnr2 == isite:
                    neighbors_from_ICOHPs.append(atomnr1)
                    lengths_from_ICOHPs.append(icohp._length)
                    icohps_from_ICOHPs.append(icohp.summed_icohp)
                    keys_from_ICOHPs.append(key)

            elif additional_condition == 1:
                # ONLY_ANION_CATION_BONDS
                if (val1 < 0.0 < val2) or (val2 < 0.0 < val1):
                    if atomnr1 == isite:
                        neighbors_from_ICOHPs.append(atomnr2)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

                    elif atomnr2 == isite:
                        neighbors_from_ICOHPs.append(atomnr1)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

            elif additional_condition == 2:
                # NO_ELEMENT_TO_SAME_ELEMENT_BONDS
                if icohp._atom1.rstrip("0123456789") != icohp._atom2.rstrip("0123456789"):
                    if atomnr1 == isite:
                        neighbors_from_ICOHPs.append(atomnr2)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

                    elif atomnr2 == isite:
                        neighbors_from_ICOHPs.append(atomnr1)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

            elif additional_condition == 3:
                # ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
                if (val1 < 0.0 < val2) or (val2 < 0.0 < val1):
                    if icohp._atom1.rstrip("0123456789") != icohp._atom2.rstrip("0123456789"):
                        if atomnr1 == isite:
                            neighbors_from_ICOHPs.append(atomnr2)
                            lengths_from_ICOHPs.append(icohp._length)
                            icohps_from_ICOHPs.append(icohp.summed_icohp)
                            keys_from_ICOHPs.append(key)

                        elif atomnr2 == isite:
                            neighbors_from_ICOHPs.append(atomnr1)
                            lengths_from_ICOHPs.append(icohp._length)
                            icohps_from_ICOHPs.append(icohp.summed_icohp)
                            keys_from_ICOHPs.append(key)

            elif additional_condition == 4:
                # ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
                if icohp._atom1.rstrip("0123456789") == "O" or icohp._atom2.rstrip("0123456789") == "O":

                    if atomnr1 == isite:
                        neighbors_from_ICOHPs.append(atomnr2)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

                    elif atomnr2 == isite:
                        neighbors_from_ICOHPs.append(atomnr1)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

            elif additional_condition == 5:
                # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
                if (val1 > 0.0 and val2 > 0.0) or (val1 < 0.0 and val2 < 0.0):
                    if atomnr1 == isite:
                        neighbors_from_ICOHPs.append(atomnr2)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

                    elif atomnr2 == isite:
                        neighbors_from_ICOHPs.append(atomnr1)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

            elif additional_condition == 6:
                # ONLY_CATION_CATION_BONDS=6
                if val1 > 0.0 and val2 > 0.0:
                    if atomnr1 == isite:
                        neighbors_from_ICOHPs.append(atomnr2)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

                    elif atomnr2 == isite:
                        neighbors_from_ICOHPs.append(atomnr1)
                        lengths_from_ICOHPs.append(icohp._length)
                        icohps_from_ICOHPs.append(icohp.summed_icohp)
                        keys_from_ICOHPs.append(key)

        return keys_from_ICOHPs, lengths_from_ICOHPs, neighbors_from_ICOHPs, icohps_from_ICOHPs

    @staticmethod
    def _get_icohps(icohpcollection, isite, lowerlimit, upperlimit, only_bonds_to):
        """
        will return icohp dict for certain site
        Args:
            icohpcollection: Icohpcollection object
            isite (int): number of a site
            lowerlimit (float): lower limit that tells you which ICOHPs are considered
            upperlimit (float): upper limit that tells you which ICOHPs are considerd
            only_bonds_to (list): list of str, e.g. ["O"] that will ensure that only bonds to "O" will be considered

        Returns:

        """
        icohps = icohpcollection.get_icohp_dict_of_site(
            site=isite,
            maxbondlength=6.0,
            minsummedicohp=lowerlimit,
            maxsummedicohp=upperlimit,
            only_bonds_to=only_bonds_to,
        )
        return icohps

    def _get_atomnumber(self, atomstring):
        """
        will return the number of the atom within the initial POSCAR (e.g., will return 0 for "Na1")
        Args:
            atomstring: string such as "Na1"

        Returns: integer indicating the position in the POSCAR

        """
        return int(self._split_string(atomstring)[1]) - 1

    @staticmethod
    def _split_string(s):
        """
        will split strings such as "Na1" in "Na" and "1" and return "1"
        Args:
            s (str): string

        Returns:

        """
        head = s.rstrip("0123456789")
        tail = s[len(head) :]
        return head, tail

    @staticmethod
    def _determine_unit_cell(site):
        """
        based on the site it will determine the unit cell, in which this site is based
        Args:
            site: site object

        Returns:

        """
        unitcell = []
        for coord in site.frac_coords:
            value = math.floor(round(coord, 4))
            unitcell.append(value)

        return unitcell

    @staticmethod
    def _get_limit_from_extremum(icohpcollection, percentage=0.15):
        """
        will return limits for the evaluation of the icohp values from an icohpcollection
        will return -100000, min(max_icohp*0.15,-0.1)
        Args:
            icohpcollection: icohpcollection object
            percentage: will determine which ICOHPs will be considered (only 0.15 from the maximum value)

        Returns: [-100000, min(max_icohp*0.15,-0.1)]

        """
        # TODO: make it work for COOPs
        extremum_based = icohpcollection.extremum_icohpvalue(summed_spin_channels=True) * percentage
        # if not self.are_coops:
        max_here = min(extremum_based, -0.1)
        return -100000, max_here
        # else:
        #    return extremum_based, 100000


class LobsterLightStructureEnvironments(LightStructureEnvironments):
    """
    Class to store LightStructureEnvironments based on Lobster outputs
    """

    @classmethod
    def from_Lobster(
        cls, list_ce_symbol, list_csm, list_permutation, list_neighsite, list_neighisite, structure, valences=None
    ):
        """
        will set up a LightStructureEnvironments from Lobster
        Args:
            structure: Structure object
            list_ce_symbol: list of symbols for coordination environments
            list_csm: list of continous symmetry measures
            list_permutation: list of permutations
            list_neighsite: list of neighboring sites
            list_neighisite: list of neighboring isites (number of a site)
            valences: list of valences

        Returns: LobsterLightStructureEnvironments

        """
        strategy = None
        valences = valences
        valences_origin = "user-defined"
        structure = structure

        coordination_environments = []

        all_nbs_sites = []
        all_nbs_sites_indices = []
        neighbors_sets = []
        counter = 0
        for isite, site in enumerate(structure):

            # all_nbs_sites_here=[]
            all_nbs_sites_indices_here = []
            # Coordination environment
            if list_ce_symbol is not None:
                ce_dict = {
                    "ce_symbol": list_ce_symbol[isite],
                    "ce_fraction": 1.0,
                    "csm": list_csm[isite],
                    "permutation": list_permutation[isite],
                }
            else:
                ce_dict = None

            if list_neighisite[isite] is not None:
                for ineighsite, neighsite in enumerate(list_neighsite[isite]):
                    diff = neighsite.frac_coords - structure[list_neighisite[isite][ineighsite]].frac_coords
                    rounddiff = np.round(diff)
                    if not np.allclose(diff, rounddiff):
                        raise ValueError(
                            "Weird, differences between one site in a periodic image cell is not " "integer ..."
                        )
                    nb_image_cell = np.array(rounddiff, np.int)

                    all_nbs_sites_indices_here.append(counter)

                    all_nbs_sites.append(
                        {"site": neighsite, "index": list_neighisite[isite][ineighsite], "image_cell": nb_image_cell}
                    )
                    counter = counter + 1

                all_nbs_sites_indices.append(all_nbs_sites_indices_here)
            else:
                all_nbs_sites.append({"site": None, "index": None, "image_cell": None})  # all_nbs_sites_here)
                all_nbs_sites_indices.append([])  # all_nbs_sites_indices_here)

            if list_neighisite[isite] is not None:
                nb_set = cls.NeighborsSet(
                    structure=structure,
                    isite=isite,
                    all_nbs_sites=all_nbs_sites,
                    all_nbs_sites_indices=all_nbs_sites_indices[isite],
                )

            else:
                nb_set = cls.NeighborsSet(structure=structure, isite=isite, all_nbs_sites=[], all_nbs_sites_indices=[])

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
    def uniquely_determines_coordination_environments(self):
        """
        True if the coordination environments are uniquely determined.
        """
        return True

    def as_dict(self):
        """
        Bson-serializable dict representation of the LightStructureEnvironments object.
        :return: Bson-serializable dict representation of the LightStructureEnvironments object.
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
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
                [nb_set.as_dict() for nb_set in site_nb_sets] if site_nb_sets is not None else None
                for site_nb_sets in self.neighbors_sets
            ],
            "valences": self.valences,
        }
