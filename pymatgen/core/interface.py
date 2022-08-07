# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides classes to store, generate, and manipulate material interfaces.
"""

from __future__ import annotations

from itertools import chain, combinations, product

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Site, Structure
from pymatgen.core.surface import Slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class Interface(Structure):
    """
    This class stores data for defining an interface between two structures.
    It is a subclass of pymatgen.core.structure.Structure.
    """

    def __init__(
        self,
        lattice,
        species,
        coords,
        site_properties,
        validate_proximity=False,
        to_unit_cell=False,
        coords_are_cartesian=False,
        in_plane_offset: tuple[float, float] = (0, 0),
        gap: float = 0,
        vacuum_over_film: float = 0.0,
        interface_properties: dict | None = None,
    ):
        """
        Makes an interface structure, a structure object with additional information
        and methods pertaining to interfaces.

        Args:
            lattice (Lattice/3x3 array): The lattice, either as a
                :class:`pymatgen.core.lattice.Lattice` or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species ([Species]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / species specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Species objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of
                each species.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in Cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
            in_plane_offset: fractional shift in plane for the film with respect
                to the substrate
            gap: gap between substrate and film in Angstroms; zero corresponds to
                the original distance between substrate and film sites
            vacuum_over_film: vacuum space above the film in Angstroms
        """

        assert (
            "interface_label" in site_properties
        ), "Must provide labeling of substrate and film sites in site properties"

        self._in_plane_offset = np.array(in_plane_offset, dtype="float")
        self._gap = gap
        self._vacuum_over_film = vacuum_over_film
        self.interface_properties = interface_properties or {}

        super().__init__(
            lattice,
            species,
            coords,
            validate_proximity=validate_proximity,
            to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties,
        )

        self.sort()

    @property
    def in_plane_offset(self) -> np.ndarray:
        """
        The shift between the film and substrate in fractional
        coordinates
        """
        return self._in_plane_offset

    @in_plane_offset.setter
    def in_plane_offset(self, new_shift: np.ndarray) -> None:
        if len(new_shift) != 2:
            raise ValueError("In-plane shifts require two floats for a and b vectors")
        new_shift = np.mod(new_shift, 1)
        delta = new_shift - np.array(self.in_plane_offset)
        self._in_plane_offset = new_shift
        self.translate_sites(self.film_indices, [delta[0], delta[1], 0], to_unit_cell=True)

    @property
    def gap(self) -> float:
        """
        The gap in Cartesian units between the film and the substrate
        """
        return self._gap

    @gap.setter
    def gap(self, new_gap: float) -> None:
        if new_gap < 0:
            raise ValueError("Can't reduce interface gap below 0")

        delta = new_gap - self.gap
        self._gap = new_gap

        self.__update_c(self.lattice.c + delta)
        self.translate_sites(self.film_indices, [0, 0, delta], frac_coords=False, to_unit_cell=True)

    @property
    def vacuum_over_film(self) -> float:
        """
        The vacuum space over the film in Cartesian units
        """
        return self._vacuum_over_film

    @vacuum_over_film.setter
    def vacuum_over_film(self, new_vacuum: float) -> None:
        if new_vacuum < 0:
            raise ValueError("The vacuum over the film can not be less then 0")

        delta = new_vacuum - self.vacuum_over_film
        self._vacuum_over_film = new_vacuum

        self.__update_c(self.lattice.c + delta)

    @property
    def substrate_indices(self) -> list[int]:
        """
        Site indices for the substrate atoms
        """
        sub_indices = [i for i, tag in enumerate(self.site_properties["interface_label"]) if "substrate" in tag]
        return sub_indices

    @property
    def substrate_sites(self) -> list[Site]:
        """
        The site objects in the substrate
        """
        sub_sites = [site for site, tag in zip(self, self.site_properties["interface_label"]) if "substrate" in tag]
        return sub_sites

    @property
    def substrate(self) -> Structure:
        """
        A pymatgen Structure for just the substrate
        """
        return Structure.from_sites(self.substrate_sites)

    @property
    def film_indices(self) -> list[int]:
        """
        Site indices of the film sites
        """
        f_indices = [i for i, tag in enumerate(self.site_properties["interface_label"]) if "film" in tag]
        return f_indices

    @property
    def film_sites(self) -> list[Site]:
        """
        Return the film sites of the interface.
        """
        film_sites = [site for site, tag in zip(self, self.site_properties["interface_label"]) if "film" in tag]
        return film_sites

    @property
    def film(self) -> Structure:
        """
        A pymatgen Structure for just the film
        """
        return Structure.from_sites(self.film_sites)

    def copy(self):
        """
        Returns:
            Interface: A copy of the Interface.
        """

        return Interface.from_dict(self.as_dict())

    def get_sorted_structure(self, key=None, reverse=False) -> Structure:
        """
        Get a sorted structure for the interface. The parameters have the same
        meaning as in list.sort. By default, sites are sorted by the
        electronegativity of the species.

        Args:
            key: Specifies a function of one argument that is used to extract
                a comparison key from each list element: key=str.lower. The
                default value is None (compare the elements directly).
            reverse (bool): If set to True, then the list elements are sorted
                as if each comparison were reversed.
        """
        struct_copy = Structure.from_sites(self)
        struct_copy.sort(key=key, reverse=reverse)
        return struct_copy

    def get_shifts_based_on_adsorbate_sites(self, tolerance: float = 0.1) -> list[tuple[float, float]]:
        """
        Computes possible in-plane shifts based on an adsorbate site  algorithm

        Args:
            tolerance: tolerance for "uniqueness" for shifts in Cartesian unit
                This is usually Angstroms.
        """
        substrate = self.substrate
        film = self.film

        substrate_surface_sites = np.dot(
            list(chain.from_iterable(AdsorbateSiteFinder(substrate).find_adsorption_sites().values())),
            substrate.lattice.inv_matrix,
        )

        # Film gets forced into substrate lattice anyways, so shifts can be computed in fractional coords
        film_surface_sites = np.dot(
            list(chain.from_iterable(AdsorbateSiteFinder(film).find_adsorption_sites().values())),
            film.lattice.inv_matrix,
        )
        pos_shift = np.array(
            [
                np.add(np.multiply(-1, film_shift), sub_shift)
                for film_shift, sub_shift in product(film_surface_sites, substrate_surface_sites)
            ]
        )

        def _base_round(x, base=0.05):
            return base * (np.array(x) / base).round()

        # Round shifts to tolerance
        pos_shift[:, 0] = _base_round(pos_shift[:, 0], base=tolerance / substrate.lattice.a)
        pos_shift[:, 1] = _base_round(pos_shift[:, 1], base=tolerance / substrate.lattice.b)
        # C-axis is not useful
        pos_shift = pos_shift[:, 0:2]

        return list(np.unique(pos_shift, axis=0))

    @property
    def film_termination(self) -> str:
        """Label for the film termination chemistry"""
        return label_termination(self.film)

    @property
    def substrate_termination(self) -> str:
        """Label for the substrate termination chemistry"""
        return label_termination(self.substrate)

    @property
    def film_layers(self) -> int:
        """Number of layers of the minimum element in the film composition"""
        sorted_element_list = sorted(
            self.film.composition.element_composition.items(), key=lambda x: x[1], reverse=True
        )
        return count_layers(self.film, sorted_element_list[0][0])

    @property
    def substrate_layers(self) -> int:
        """Number of layers of the minimum element in the substrate composition"""
        sorted_element_list = sorted(
            self.substrate.composition.element_composition.items(), key=lambda x: x[1], reverse=True
        )
        return count_layers(self.substrate, sorted_element_list[0][0])

    def __update_c(self, new_c: float) -> None:
        """
        Modifies the c-direction of the lattice without changing the site Cartesian coordinates
        Be careful you can mess up the interface by setting a c-length that can't accommodate all the sites
        """
        if new_c <= 0:
            raise ValueError("New c-length must be greater than 0")

        new_latt_matrix = self.lattice.matrix[:2].tolist() + [[0, 0, new_c]]
        new_latice = Lattice(new_latt_matrix)
        self._lattice = new_latice

        for site, c_coords in zip(self, self.cart_coords):
            site._lattice = new_latice  # Update the lattice
            site.coords = c_coords  # Put back into original Cartesian space

    def as_dict(self):
        """
        :return: MSONAble dict
        """
        d = super().as_dict()
        d["in_plane_offset"] = self.in_plane_offset.tolist()
        d["gap"] = self.gap
        d["vacuum_over_film"] = self.vacuum_over_film
        d["interface_properties"] = self.interface_properties
        return d

    @classmethod
    def from_dict(cls, d):
        """
        :param d: dict
        :return: Creates slab from dict.
        """
        lattice = Lattice.from_dict(d["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in d["sites"]]
        s = Structure.from_sites(sites)

        optional = dict(
            in_plane_offset=d.get("in_plane_offset"),
            gap=d.get("gap"),
            vacuum_over_film=d.get("vacuum_over_film"),
            interface_properties=d.get("interface_properties"),
        )
        return Interface(
            lattice=lattice,
            species=s.species_and_occu,
            coords=s.frac_coords,
            site_properties=s.site_properties,
            **{k: v for k, v in optional.items() if v is not None},
        )

    @classmethod
    def from_slabs(
        cls,
        substrate_slab: Slab,
        film_slab: Slab,
        in_plane_offset: tuple[float, float] = (0, 0),
        gap: float = 1.6,
        vacuum_over_film: float = 0.0,
        interface_properties: dict | None = None,
        center_slab: bool = True,
    ) -> Interface:
        """
        Makes an interface structure by merging a substrate and film slabs
        The film a- and b-vectors will be forced to be the substrate slab's
        a- and b-vectors.

        For now, it's suggested to use a factory method that will ensure the
        appropriate interface structure is already met.

        Args:
            sub_slab: slab for the substrate
            film_slab: slab for the film
            in_plane_offset: fractional shift in plane
                for the film with respect to the substrate
            gap: gap between substrate and film in Angstroms
            vacuum_over_film: vacuum space above the film in Angstroms
            structure_properties: dictionary of misc properties for this structure
            center_slab: center the slab
        """
        interface_properties = interface_properties or {}

        # Ensure c-axis is orthogonal to a/b plane
        if isinstance(substrate_slab, Slab):
            substrate_slab = substrate_slab.get_orthogonal_c_slab()
        if isinstance(film_slab, Slab):
            film_slab = film_slab.get_orthogonal_c_slab()
        assert np.allclose(film_slab.lattice.alpha, 90, 0.1)
        assert np.allclose(film_slab.lattice.beta, 90, 0.1)
        assert np.allclose(substrate_slab.lattice.alpha, 90, 0.1)
        assert np.allclose(substrate_slab.lattice.beta, 90, 0.1)

        # Ensure sub is right-handed
        # IE sub has surface facing "up"
        sub_vecs = substrate_slab.lattice.matrix.copy()
        if np.dot(np.cross(*sub_vecs[:2]), sub_vecs[2]) < 0:
            sub_vecs[2] *= -1.0
            substrate_slab.lattice = Lattice(sub_vecs)

        # Find the limits of C-coords
        sub_coords = substrate_slab.frac_coords
        film_coords = film_slab.frac_coords
        sub_min_c = np.min(sub_coords[:, 2]) * substrate_slab.lattice.c
        sub_max_c = np.max(sub_coords[:, 2]) * substrate_slab.lattice.c
        film_min_c = np.min(film_coords[:, 2]) * film_slab.lattice.c
        film_max_c = np.max(film_coords[:, 2]) * film_slab.lattice.c
        min_height = np.abs(film_max_c - film_min_c) + np.abs(sub_max_c - sub_min_c)

        # construct new lattice
        abc = substrate_slab.lattice.abc[:2] + (min_height + gap + vacuum_over_film,)
        angles = substrate_slab.lattice.angles
        lattice = Lattice.from_parameters(*abc, *angles)

        # Get the species
        species = substrate_slab.species + film_slab.species

        # Get the coords
        # Shift substrate to bottom in new lattice
        sub_coords = np.subtract(sub_coords, [0, 0, np.min(sub_coords[:, 2])])
        sub_coords[:, 2] *= substrate_slab.lattice.c / lattice.c

        # Flip the film over
        film_coords[:, 2] *= -1.0
        film_coords[:, 2] *= film_slab.lattice.c / lattice.c

        # Shift the film coords to right over the substrate + gap
        film_coords = np.subtract(film_coords, [0, 0, np.min(film_coords[:, 2])])
        film_coords = np.add(film_coords, [0, 0, gap / lattice.c + np.max(sub_coords[:, 2])])

        # Build coords
        coords = np.concatenate([sub_coords, film_coords])

        # Shift coords to center
        if center_slab:
            coords = np.add(coords, [0, 0, 0.5 - np.average(coords[:, 2])])

        # Only merge site properties in both slabs
        site_properties = {}
        site_props_in_both = set(substrate_slab.site_properties) & set(film_slab.site_properties)

        for key in site_props_in_both:
            site_properties[key] = [
                *substrate_slab.site_properties[key],
                *film_slab.site_properties[key],
            ]

        site_properties["interface_label"] = ["substrate"] * len(substrate_slab) + ["film"] * len(film_slab)

        iface = cls(
            lattice=lattice,
            species=species,
            coords=coords,
            to_unit_cell=False,
            coords_are_cartesian=False,
            site_properties=site_properties,
            validate_proximity=False,
            in_plane_offset=in_plane_offset,
            gap=gap,
            vacuum_over_film=vacuum_over_film,
            interface_properties=interface_properties,
        )

        iface.sort()
        return iface


def label_termination(slab: Structure) -> str:
    """Labels the slab surface termination"""
    frac_coords = slab.frac_coords
    n = len(frac_coords)

    if n == 1:
        # Clustering does not work when there is only one data point.
        form = slab.composition.reduced_formula
        sp_symbol = SpacegroupAnalyzer(slab, symprec=0.1).get_space_group_symbol()
        return f"{form}_{sp_symbol}_{len(slab)}"

    dist_matrix = np.zeros((n, n))
    h = slab.lattice.c
    # Projection of c lattice vector in
    # direction of surface normal.
    for i, j in combinations(list(range(n)), 2):
        if i != j:
            cdist = frac_coords[i][2] - frac_coords[j][2]
            cdist = abs(cdist - round(cdist)) * h
            dist_matrix[i, j] = cdist
            dist_matrix[j, i] = cdist

    condensed_m = squareform(dist_matrix)
    z = linkage(condensed_m)
    clusters = fcluster(z, 0.25, criterion="distance")

    clustered_sites: dict[int, list[Site]] = {c: [] for c in clusters}
    for i, c in enumerate(clusters):
        clustered_sites[c].append(slab[i])

    plane_heights = {
        np.average(np.mod([s.frac_coords[2] for s in sites], 1)): c for c, sites in clustered_sites.items()
    }
    top_plane_cluster = sorted(plane_heights.items(), key=lambda x: x[0])[-1][1]
    top_plane_sites = clustered_sites[top_plane_cluster]
    top_plane = Structure.from_sites(top_plane_sites)

    sp_symbol = SpacegroupAnalyzer(top_plane, symprec=0.1).get_space_group_symbol()
    form = top_plane.composition.reduced_formula
    return f"{form}_{sp_symbol}_{len(top_plane)}"


def count_layers(struc: Structure, el=None) -> int:
    """
    Counts the number of 'layers' along the c-axis
    """
    el = el if el else struc.composition.elements[0]
    frac_coords = [site.frac_coords for site in struc if site.species_string == str(el)]
    n = len(frac_coords)

    if n == 1:
        return 1

    dist_matrix = np.zeros((n, n))
    h = struc.lattice.c
    # Projection of c lattice vector in
    # direction of surface normal.
    for i, j in combinations(list(range(n)), 2):
        if i != j:
            cdist = frac_coords[i][2] - frac_coords[j][2]
            cdist = abs(cdist - round(cdist)) * h
            dist_matrix[i, j] = cdist
            dist_matrix[j, i] = cdist

    condensed_m = squareform(dist_matrix)
    z = linkage(condensed_m)
    clusters = fcluster(z, 0.25, criterion="distance")

    clustered_sites: dict[int, list[Site]] = {c: [] for c in clusters}
    for i, c in enumerate(clusters):
        clustered_sites[c].append(struc[i])

    plane_heights = {
        np.average(np.mod([s.frac_coords[2] for s in sites], 1)): c for c, sites in clustered_sites.items()
    }

    return len(plane_heights)
