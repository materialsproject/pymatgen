"""This module implements representation of Slab, SlabGenerator
for generating Slabs, ReconstructionGenerator to generate
reconstructed Slabs, and some related utility functions.

If you use this module, please consider citing the following work:

    R. Tran, Z. Xu, B. Radhakrishnan, D. Winston, W. Sun, K. A. Persson,
    S. P. Ong, "Surface Energies of Elemental Crystals", Scientific Data,
    2016, 3:160080, doi: 10.1038/sdata.2016.80.

    Sun, W.; Ceder, G. Efficient creation and convergence of surface slabs,
    Surface Science, 2013, 617, 53-59, doi:10.1016/j.susc.2013.05.016.
"""

from __future__ import annotations

import copy
import itertools
import json
import logging
import math
import os
import warnings
from functools import reduce
from math import gcd, isclose
from typing import TYPE_CHECKING, cast

import numpy as np
from monty.fractions import lcm
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice, PeriodicSite, Structure, get_el_sp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import in_coord_list
from pymatgen.util.due import Doi, due
from pymatgen.util.typing import Tuple3Ints

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from numpy.typing import ArrayLike, NDArray
    from typing_extensions import Self

    from pymatgen.core.composition import Element, Species
    from pymatgen.symmetry.groups import CrystalSystem
    from pymatgen.util.typing import MillerIndex

__author__ = "Richard Tran, Wenhao Sun, Zihan Xu, Shyue Ping Ong"

due.dcite(
    Doi("10.1038/sdata.2016.80"),
    description="Surface Energies of Elemental Crystals",
)
due.dcite(
    Doi("10.1016/j.susc.2013.05.016"),
    description="Efficient creation and convergence of surface slabs",
)

logger = logging.getLogger(__name__)


class Slab(Structure):
    """Hold information for a Slab, with additional
    attributes pertaining to slabs, but the init method does not
    actually create a slab. Also has additional methods that returns other information
    about a Slab such as the surface area, normal, and atom adsorption.

    Note that all Slabs have the surface normal oriented perpendicular to the
    a and b lattice vectors. This means the lattice vectors a and b are in the
    surface plane and the c vector is out of the surface plane (though not
    necessarily perpendicular to the surface).
    """

    def __init__(
        self,
        lattice: Lattice | np.ndarray,
        species: Sequence[Any],
        coords: np.ndarray,
        miller_index: MillerIndex,
        oriented_unit_cell: Structure,
        shift: float,
        scale_factor: np.ndarray,
        reorient_lattice: bool = True,
        validate_proximity: bool = False,
        to_unit_cell: bool = False,
        reconstruction: str | None = None,
        coords_are_cartesian: bool = False,
        site_properties: dict | None = None,
        energy: float | None = None,
    ) -> None:
        """A Structure object with additional information
        and methods pertaining to Slabs.

        Args:
            lattice (Lattice/3x3 array): The lattice, either as a
                pymatgen.core.Lattice or simply as any 2D array.
                Each row should correspond to a lattice
                vector. e.g. [[10,0,0], [20,10,0], [0,0,30]].
            species ([Species]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / species specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g. (3, 56, ...) or actual Element or Species objects.

                ii. List of dict of elements/species and occupancies, e.g.
                    [{"Fe": 0.5, "Mn": 0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of each species.
            miller_index (MillerIndex): Miller index of plane parallel to
                surface. Note that this is referenced to the input structure. If
                you need this to be based on the conventional cell,
                you should supply the conventional structure.
            oriented_unit_cell (Structure): The oriented_unit_cell from which
                this Slab is created (by scaling in the c-direction).
            shift (float): The NEGATIVE of shift in the c-direction applied
                to get the termination.
            scale_factor (np.ndarray): scale_factor Final computed scale factor
                that brings the parent cell to the surface cell.
            reorient_lattice (bool): reorients the lattice parameters such that
                the c direction is along the z axis.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            reconstruction (str): Type of reconstruction. Defaults to None if
                the slab is not reconstructed.
            to_unit_cell (bool): Translates fractional coordinates into the
                unit cell. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in Cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g. {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
            energy (float): A value for the energy.
        """
        self.oriented_unit_cell = oriented_unit_cell
        self.miller_index = miller_index
        self.shift = shift
        self.reconstruction = reconstruction
        self.scale_factor = scale_factor
        self.energy = energy
        self.reorient_lattice = reorient_lattice

        if self.reorient_lattice:
            if coords_are_cartesian:
                coords = lattice.get_fractional_coords(coords)
                coords_are_cartesian = False
            lattice = Lattice.from_parameters(
                lattice.a,
                lattice.b,
                lattice.c,
                lattice.alpha,
                lattice.beta,
                lattice.gamma,
            )

        super().__init__(
            lattice,
            species,
            coords,
            validate_proximity=validate_proximity,
            to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties,
        )

    def __str__(self) -> str:
        outs = [
            f"Slab Summary ({self.composition.formula})",
            f"Reduced Formula: {self.composition.reduced_formula}",
            f"Miller index: {self.miller_index}",
            f"Shift: {self.shift:.4f}, Scale Factor: {self.scale_factor}",
            f"abc   : {' '.join(f'{i:0.6f}'.rjust(10) for i in self.lattice.abc)}",
            f"angles: {' '.join(f'{i:0.6f}'.rjust(10) for i in self.lattice.angles)}",
            f"Sites ({len(self)})",
        ]

        for idx, site in enumerate(self):
            outs.append(f"{idx + 1} {site.species_string} {' '.join(f'{j:0.6f}'.rjust(12) for j in site.frac_coords)}")

        return "\n".join(outs)

    @property
    def center_of_mass(self) -> np.ndarray:
        """The center of mass of the Slab in fractional coordinates."""
        weights = [site.species.weight for site in self]
        return np.average(self.frac_coords, weights=weights, axis=0)

    @property
    def dipole(self) -> np.ndarray:
        """The dipole moment of the Slab in the direction of the surface normal.

        Note that the Slab must be oxidation state decorated for this to work properly.
        Otherwise, the Slab will always have a dipole moment of 0.
        """
        centroid = np.sum(self.cart_coords, axis=0) / len(self)

        dipole = np.zeros(3)
        for site in self:
            charge = sum(getattr(sp, "oxi_state", 0) * amt for sp, amt in site.species.items())
            dipole += charge * np.dot(site.coords - centroid, self.normal) * self.normal
        return dipole

    @property
    def normal(self) -> np.ndarray:
        """The surface normal vector of the Slab, normalized to unit length."""
        normal = np.cross(self.lattice.matrix[0], self.lattice.matrix[1])
        normal /= np.linalg.norm(normal)
        return normal

    @property
    def surface_area(self) -> float:
        """The surface area of the Slab."""
        matrix = self.lattice.matrix
        return np.linalg.norm(np.cross(matrix[0], matrix[1]))

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct: dict.

        Returns:
            Slab: Created from dict.
        """
        lattice = Lattice.from_dict(dct["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in dct["sites"]]
        struct = Structure.from_sites(sites)

        return cls(
            lattice=lattice,
            species=struct.species_and_occu,
            coords=struct.frac_coords,
            miller_index=dct["miller_index"],
            oriented_unit_cell=Structure.from_dict(dct["oriented_unit_cell"]),
            shift=dct["shift"],
            scale_factor=np.array(dct["scale_factor"]),
            site_properties=struct.site_properties,
            energy=dct["energy"],
        )

    def as_dict(self, **kwargs) -> dict:
        """MSONable dict."""
        dct = super().as_dict(**kwargs)
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["oriented_unit_cell"] = self.oriented_unit_cell.as_dict()
        dct["miller_index"] = self.miller_index
        dct["shift"] = self.shift
        dct["scale_factor"] = self.scale_factor.tolist()  # np.ndarray is not JSON serializable
        dct["reconstruction"] = self.reconstruction
        dct["energy"] = self.energy
        return dct

    def copy(self, site_properties: dict[str, Any] | None = None) -> Self:
        """Get a copy of the Slab, with options to update site properties.

        Args:
            site_properties (dict): Properties to update. The
                properties are specified in the same way as the constructor,
                i.e., as a dict of the form {property: [values]}.

        Returns:
            A copy of the Structure, with optionally new site_properties
        """
        props = self.site_properties
        if site_properties:
            props.update(site_properties)

        return type(self)(
            self.lattice,
            self.species_and_occu,
            self.frac_coords,
            self.miller_index,
            self.oriented_unit_cell,
            self.shift,
            self.scale_factor,
            site_properties=props,
            reorient_lattice=self.reorient_lattice,
        )

    def is_symmetric(self, symprec: float = 0.1) -> bool:
        """Check if Slab is symmetric, i.e., contains inversion, mirror on (hkl) plane,
            or screw axis (rotation and translation) about [hkl].

        Args:
            symprec (float): Symmetry precision used for SpaceGroup analyzer.

        Returns:
            bool: True if surfaces are symmetric.
        """
        spg_analyzer = SpacegroupAnalyzer(self, symprec=symprec)
        symm_ops = spg_analyzer.get_point_group_operations()

        # Check for inversion symmetry. Or if sites from surface (a) can be translated
        # to surface (b) along the [hkl]-axis, surfaces are symmetric. Or because the
        # two surfaces of our slabs are always parallel to the (hkl) plane,
        # any operation where there's an (hkl) mirror plane has surface symmetry
        return bool(
            spg_analyzer.is_laue()
            or any(op.translation_vector[2] != 0 for op in symm_ops)
            or any(np.all(op.rotation_matrix[2] == np.array([0, 0, -1])) for op in symm_ops)
        )

    def is_polar(self, tol_dipole_per_unit_area: float = 1e-3) -> bool:
        """Check if the Slab is polar by computing the normalized dipole per unit area.
        Normalized dipole per unit area is used as it is more reliable than
        using the absolute value, which varies with surface area.

        Note that the Slab must be oxidation state decorated for this to work properly.
        Otherwise, the Slab will always have a dipole moment of 0.

        Args:
            tol_dipole_per_unit_area (float): A tolerance above which the Slab is
                considered polar.
        """
        dip_per_unit_area = self.dipole / self.surface_area
        return bool(np.linalg.norm(dip_per_unit_area) > tol_dipole_per_unit_area)

    def get_surface_sites(self, tag: bool = False) -> dict[str, list]:
        """Get the surface sites and their indices in a dictionary.
        Useful for analysis involving broken bonds and for finding adsorption sites.

        The oriented unit cell of the slab will determine the
        coordination number of a typical site.
        We use VoronoiNN to determine the coordination number of sites.
        Due to the pathological error resulting from some surface sites in the
        VoronoiNN, we assume any site that has this error is a surface
        site as well. This will only work for single-element systems for now.

        Args:
            tag (bool): Add attribute "is_surf_site" (bool)
                to all sites of the Slab. Defaults to False.

        Returns:
            A dictionary grouping sites on top and bottom of the slab together.
                {"top": [sites with indices], "bottom": [sites with indices]}

        Todo:
            Is there a way to determine site equivalence between sites in a slab
            and bulk system? This would allow us get the coordination number of
            a specific site for multi-elemental systems or systems with more
            than one inequivalent site. This will allow us to use this for
            compound systems.
        """
        from pymatgen.analysis.local_env import VoronoiNN

        # Get a dictionary of coordination numbers for each distinct site in the structure
        spg_analyzer = SpacegroupAnalyzer(self.oriented_unit_cell)
        u_cell = spg_analyzer.get_symmetrized_structure()
        cn_dict: dict = {}
        voronoi_nn = VoronoiNN()
        unique_indices = [equ[0] for equ in u_cell.equivalent_indices]

        for idx in unique_indices:
            el = u_cell[idx].species_string
            if el not in cn_dict:
                cn_dict[el] = []
            # Since this will get the CN as a result of the weighted polyhedra, the
            # slightest difference in CN will indicate a different environment for a
            # species, eg. bond distance of each neighbor or neighbor species. The
            # decimal place to get some CN to be equal.
            cn = voronoi_nn.get_cn(u_cell, idx, use_weights=True)
            cn = float(f"{round(cn, 5):.5f}")
            if cn not in cn_dict[el]:
                cn_dict[el].append(cn)

        voronoi_nn = VoronoiNN()

        surf_sites_dict: dict = {"top": [], "bottom": []}
        properties: list = []
        for idx, site in enumerate(self):
            # Determine if site is closer to the top or bottom of the slab
            is_top: bool = site.frac_coords[2] > self.center_of_mass[2]

            try:
                # A site is a surface site, if its environment does
                # not fit the environment of other sites
                cn = float(f"{round(voronoi_nn.get_cn(self, idx, use_weights=True), 5):.5f}")
                if cn < min(cn_dict[site.species_string]):
                    properties.append(True)
                    key = "top" if is_top else "bottom"
                    surf_sites_dict[key].append([site, idx])
                else:
                    properties.append(False)
            except RuntimeError:
                # or if pathological error is returned, indicating a surface site
                properties.append(True)
                key = "top" if is_top else "bottom"
                surf_sites_dict[key].append([site, idx])

        if tag:
            self.add_site_property("is_surf_site", properties)
        return surf_sites_dict

    def get_symmetric_site(
        self,
        point: ArrayLike,
        cartesian: bool = False,
    ) -> ArrayLike:
        """Use symmetry operations to find an equivalent site on the other side of
        the slab. Works mainly for slabs with Laue symmetry.

        This is useful for retaining the non-polar and
        symmetric properties of a slab when creating adsorbed
        structures or symmetric reconstructions.

        Args:
            point (ArrayLike): Fractional coordinate of the original site.
            cartesian (bool): Use Cartesian coordinates.

        Returns:
            ArrayLike: Fractional coordinate. A site equivalent to the
                original site, but on the other side of the slab
        """
        spg_analyzer = SpacegroupAnalyzer(self)
        ops = spg_analyzer.get_symmetry_operations(cartesian=cartesian)

        # Each operation on a site will return an equivalent site.
        # We want to find the site on the other side of the slab.
        site_other = None
        for op in ops:
            slab = self.copy()
            site_other = op.operate(point)
            if isclose(site_other[2], point[2], abs_tol=1e-6):
                continue

            # Add dummy sites to check if the overall structure is symmetric
            slab.append("O", point, coords_are_cartesian=cartesian)
            slab.append("O", site_other, coords_are_cartesian=cartesian)
            if SpacegroupAnalyzer(slab).is_laue():
                break

            # If not symmetric, remove the two added
            # sites and try another symmetry operator
            slab.remove_sites([len(slab) - 1])
            slab.remove_sites([len(slab) - 1])

        if site_other is None:
            raise RuntimeError("Failed to get symmetric site.")

        return site_other

    def get_orthogonal_c_slab(self) -> Self:
        """Generate a Slab where the normal (c lattice vector) is
        forced to be orthogonal to the surface a and b lattice vectors.

        **Note that this breaks inherent symmetries in the slab.**

        It should be pointed out that orthogonality is not required to get good
        surface energies, but it can be useful in cases where the slabs are
        subsequently used for postprocessing of some kind, e.g. generating
        grain boundaries or interfaces.
        """
        a, b, c = self.lattice.matrix
        _new_c = np.cross(a, b)
        _new_c /= np.linalg.norm(_new_c)
        new_c = np.dot(c, _new_c) * _new_c
        new_latt = Lattice([a, b, new_c])

        return type(self)(
            lattice=new_latt,
            species=self.species_and_occu,
            coords=self.cart_coords,
            miller_index=self.miller_index,
            oriented_unit_cell=self.oriented_unit_cell,
            shift=self.shift,
            scale_factor=self.scale_factor,
            coords_are_cartesian=True,
            energy=self.energy,
            reorient_lattice=self.reorient_lattice,
            site_properties=self.site_properties,
        )

    def get_tasker2_slabs(
        self,
        tol: float = 0.01,
        same_species_only: bool = True,
    ) -> list[Self]:
        """Get a list of slabs that have been Tasker 2 corrected.

        Args:
            tol (float): Fractional tolerance to determine if atoms are within same plane.
            same_species_only (bool): If True, only those are of the exact same
                species as the atom at the outermost surface are considered for moving.
                Otherwise, all atoms regardless of species within tol are considered for moving.
                Default is True (usually the desired behavior).

        Returns:
            list[Slab]: Tasker 2 corrected slabs.
        """

        def get_equi_index(site: PeriodicSite) -> int:
            """Get the index of the equivalent site for a given site."""
            for idx, equi_sites in enumerate(symm_structure.equivalent_sites):
                if site in equi_sites:
                    return idx
            raise ValueError("Cannot determine equi index!")

        sites = list(self.sites)
        slabs = []

        sorted_csites = sorted(sites, key=lambda site: site.c)

        # Determine what fraction the slab is of the total cell size in the
        # c direction. Round to nearest rational number.
        n_layers_total = int(round(self.lattice.c / self.oriented_unit_cell.lattice.c))
        n_layers_slab = int(round((sorted_csites[-1].c - sorted_csites[0].c) * n_layers_total))
        slab_ratio = n_layers_slab / n_layers_total

        spg_analyzer = SpacegroupAnalyzer(self)
        symm_structure = spg_analyzer.get_symmetrized_structure()

        for surface_site, shift in [
            (sorted_csites[0], slab_ratio),
            (sorted_csites[-1], -slab_ratio),
        ]:
            to_move = []
            fixed = []
            for site in sites:
                if abs(site.c - surface_site.c) < tol and (
                    (not same_species_only) or site.species == surface_site.species
                ):
                    to_move.append(site)
                else:
                    fixed.append(site)

            # Sort and group the sites by the species and symmetry equivalence
            to_move = sorted(to_move, key=get_equi_index)

            grouped = [list(sites) for k, sites in itertools.groupby(to_move, key=get_equi_index)]

            if len(to_move) == 0 or any(len(g) % 2 != 0 for g in grouped):
                warnings.warn(
                    "Odd number of sites to divide! Try changing "
                    "the tolerance to ensure even division of "
                    "sites or create supercells in a or b directions "
                    "to allow for atoms to be moved!"
                )
                continue
            combinations = []
            for g in grouped:
                combinations.append(list(itertools.combinations(g, int(len(g) / 2))))

            for selection in itertools.product(*combinations):
                species = [site.species for site in fixed]
                frac_coords = [site.frac_coords for site in fixed]

                for struct_matcher in to_move:
                    species.append(struct_matcher.species)
                    for group in selection:
                        if struct_matcher in group:
                            frac_coords.append(struct_matcher.frac_coords)
                            break
                    else:
                        # Move unselected atom to the opposite surface.
                        frac_coords.append(struct_matcher.frac_coords + np.array([0, 0, shift]))

                # sort by species to put all similar species together.
                sp_fcoord = sorted(zip(species, frac_coords, strict=True), key=lambda x: x[0])
                species = [x[0] for x in sp_fcoord]
                frac_coords = [x[1] for x in sp_fcoord]
                slab = type(self)(
                    self.lattice,
                    species,
                    frac_coords,
                    self.miller_index,
                    self.oriented_unit_cell,
                    self.shift,
                    self.scale_factor,
                    energy=self.energy,
                    reorient_lattice=self.reorient_lattice,
                )
                slabs.append(slab)
        struct_matcher = StructureMatcher()
        return [ss[0] for ss in struct_matcher.group_structures(slabs)]

    def get_sorted_structure(self, key=None, reverse: bool = False) -> Self:
        """Get a sorted copy of the structure. The parameters have the same
        meaning as in list.sort. By default, sites are sorted by the
        electronegativity of the species. Note that Slab has to override this
        because of the different __init__ args.

        Args:
            key: Specifies a function of one argument that is used to extract
                a comparison key from each list element: key=str.lower. The
                default value is None (compare the elements directly).
            reverse (bool): If set to True, then the list elements are sorted
                as if each comparison were reversed.
        """
        sites = sorted(self, key=key, reverse=reverse)
        struct = Structure.from_sites(sites)
        return type(self)(
            struct.lattice,
            struct.species_and_occu,
            struct.frac_coords,
            self.miller_index,
            self.oriented_unit_cell,
            self.shift,
            self.scale_factor,
            site_properties=struct.site_properties,
            reorient_lattice=self.reorient_lattice,
        )

    def add_adsorbate_atom(
        self,
        indices: list[int],
        species: str | Element | Species,
        distance: float,
        specie: Species | Element | str | None = None,
    ) -> Self:
        """Add adsorbate onto the Slab, along the c lattice vector.

        Args:
            indices (list[int]): Indices of sites on which to put the adsorbate.
                Adsorbate will be placed relative to the center of these sites.
            species (str | Element | Species): The species to add.
            distance (float): between centers of the adsorbed atom and the
                given site in Angstroms, along the c lattice vector.
            specie: Deprecated argument in #3691. Use 'species' instead.

        Returns:
            Slab: self with adsorbed atom.
        """
        # Check if deprecated argument is used
        if specie is not None:
            warnings.warn(
                "The argument 'specie' is deprecated. Use 'species' instead.",
                DeprecationWarning,
            )
            species = specie

        # Calculate target site as the center of sites
        center = np.sum([self[idx].coords for idx in indices], axis=0) / len(indices)

        coords = center + self.normal * distance

        self.append(species, coords, coords_are_cartesian=True)

        return self

    def symmetrically_add_atom(
        self,
        species: str | Element | Species,
        point: ArrayLike,
        specie: str | Element | Species | None = None,
        coords_are_cartesian: bool = False,
    ) -> None:
        """Add a species at a selected site in a Slab. Will also add an
        equivalent site on the other side to maintain symmetry.

        Args:
            species (str | Element | Species): The species to add.
            point (ArrayLike): The coordinate of the target site.
            specie: Deprecated argument name in #3691. Use 'species' instead.
            coords_are_cartesian (bool): If the site is in Cartesian coordinates.
        """
        # Check if deprecated argument is used
        if specie is not None:
            warnings.warn(
                "The argument 'specie' is deprecated. Use 'species' instead.",
                DeprecationWarning,
            )
            species = specie

        # Get the index of the equivalent site on the other side
        equi_site = self.get_symmetric_site(point, cartesian=coords_are_cartesian)

        self.append(species, point, coords_are_cartesian=coords_are_cartesian)
        self.append(species, equi_site, coords_are_cartesian=coords_are_cartesian)

    def symmetrically_remove_atoms(self, indices: list[int]) -> None:
        """Remove sites from a list of indices. Will also remove the
        equivalent site on the other side of the slab to maintain symmetry.

        Args:
            indices (list[int]): The indices of the sites to remove.

        TODO(@DanielYang59):
        1. Reuse public method get_symmetric_site to get equi sites?
        2. If not 1, get_equi_sites has multiple nested loops
        """

        def get_equi_sites(slab: Slab, sites: list[int]) -> list[int]:
            """
            Get the indices of the equivalent sites of given sites.

            Parameters:
                slab (Slab): The slab structure.
                sites (list[int]): Original indices of sites.

            Returns:
                list[int]: Indices of the equivalent sites.
            """
            equi_sites = []
            eq_indices = []
            eq_sites: list = []

            for pt in sites:
                # Get the index of the original site
                cart_point = slab.lattice.get_cartesian_coords(pt)
                dist = [site.distance_from_point(cart_point) for site in slab]
                site1 = dist.index(min(dist))

                # Get the index of the equivalent site on the other side
                for i, eq_sites in enumerate(slab.equivalent_sites):
                    if slab[site1] in eq_sites:
                        eq_indices = slab.equivalent_indices[i]
                        break
                i1 = eq_indices[eq_sites.index(slab[site1])]

                for i2 in eq_indices:
                    if i2 == i1:
                        continue
                    if slab[i2].frac_coords[2] == slab[i1].frac_coords[2]:
                        continue
                    # Test site remove to see if it results in symmetric slab
                    slab = self.copy()
                    slab.remove_sites([i1, i2])
                    if slab.is_symmetric():
                        equi_sites.append(i2)
                        break

            return equi_sites

        # Generate the equivalent sites of the original sites
        slab_copy = SpacegroupAnalyzer(self.copy()).get_symmetrized_structure()
        sites = [slab_copy[i].frac_coords for i in indices]

        equi_sites = get_equi_sites(slab_copy, sites)

        # Check if found any equivalent sites
        if len(equi_sites) == len(indices):
            self.remove_sites(indices)
            self.remove_sites(equi_sites)

        else:
            warnings.warn("Equivalent sites could not be found for some indices. Surface unchanged.")


def center_slab(slab: Structure) -> Structure:
    """Relocate the slab to the center such that its center
    (the slab region) is close to z=0.5.

    This makes it easier to find surface sites and apply
    operations like doping.

    There are two possible cases:
        1. When the slab region is completely positioned between
        two vacuum layers in the cell but is not centered, we simply
        shift the slab to the center along z-axis.
        2. If the slab completely resides outside the cell either
        from the bottom or the top, we iterate through all sites that
        spill over and shift all sites such that it is now
        on the other side. An edge case being, either the top
        of the slab is at z = 0 or the bottom is at z = 1.

    Args:
        slab (Structure): The slab to center.

    Returns:
        Structure: The centered slab.
    """
    # Get all site indices
    all_indices = list(range(len(slab)))

    # Get a reasonable cutoff radius to sample neighbors
    bond_dists = sorted(nn[1] for nn in slab.get_neighbors(slab[0], 10) if nn[1] > 0)
    # TODO (@DanielYang59): magic number for cutoff radius (would 3 be too large?)
    cutoff_radius = bond_dists[0] * 3

    # TODO (@DanielYang59): do we need the following complex method?
    # Why don't we just calculate the center of the Slab and move it to z=0.5?
    # Before moving we need to ensure there is only one Slab layer though

    # If structure is case 2, shift all the sites
    # to the other side until it is case 1
    for site in slab:  # DEBUG (@DanielYang59): Slab position changes during loop?
        # DEBUG (@DanielYang59): sites below z=0 is not considered (only check coord > c)
        if any(nn[1] >= slab.lattice.c for nn in slab.get_neighbors(site, cutoff_radius)):
            # TODO (@DanielYang59): the magic offset "0.05" seems unnecessary,
            # as the Slab would be centered later anyway
            shift = 1 - site.frac_coords[2] + 0.05
            slab.translate_sites(all_indices, [0, 0, shift])

    # Now the slab is case 1, move it to the center
    weights = [site.species.weight for site in slab]
    center_of_mass = np.average(slab.frac_coords, weights=weights, axis=0)
    shift = 0.5 - center_of_mass[2]

    slab.translate_sites(all_indices, [0, 0, shift])

    return slab


def get_slab_regions(
    slab: Slab,
    blength: float = 3.5,
) -> list[tuple[float, float]]:
    """Find the z-ranges for the slab region.

    Useful for discerning where the slab ends and vacuum begins
    if the slab is not fully within the cell.

    Args:
        slab (Slab): The Slab to analyse.
        blength (float): The bond length between atoms in Angstrom.
            You generally want this value to be larger than the actual
            bond length in order to find atoms that are part of the slab.

    TODO (@DanielYang59): this should be a method for `Slab`?
    TODO (@DanielYang59): maybe project all z coordinates to 1D?
    """
    frac_coords: list = []  # TODO (@DanielYang59): zip site and coords?
    indices: list = []

    all_indices: list = []

    for site in slab:
        neighbors = slab.get_neighbors(site, blength)
        for nn in neighbors:
            # TODO (@DanielYang59): use z coordinate (z<0) to check
            # if a Slab is contiguous is suspicious (Slab could locate
            # entirely below z=0)

            # Find sites with z < 0 (sites noncontiguous within cell)
            if nn[0].frac_coords[2] < 0:
                frac_coords.append(nn[0].frac_coords[2])
                indices.append(nn[-2])

                if nn[-2] not in all_indices:
                    all_indices.append(nn[-2])

    # If slab is noncontiguous
    if frac_coords:
        # Locate the lowest site within the upper Slab
        last_frac_coords = []
        last_indices = []
        while frac_coords:
            last_frac_coords = copy.copy(frac_coords)
            last_indices = copy.copy(indices)

            site = slab[indices[frac_coords.index(min(frac_coords))]]
            neighbors = slab.get_neighbors(site, blength, include_index=True, include_image=True)
            frac_coords, indices = [], []
            for nn in neighbors:
                if 1 > nn[0].frac_coords[2] > 0 and nn[0].frac_coords[2] < site.frac_coords[2]:
                    # Sites are noncontiguous within cell
                    frac_coords.append(nn[0].frac_coords[2])
                    indices.append(nn[-2])
                    if nn[-2] not in all_indices:
                        all_indices.append(nn[-2])

        # Locate the highest site within the lower Slab
        upper_fcoords: list = []
        for site in slab:
            if all(nn.index not in all_indices for nn in slab.get_neighbors(site, blength)):
                upper_fcoords.append(site.frac_coords[2])
        coords: list = copy.copy(frac_coords) if frac_coords else copy.copy(last_frac_coords)
        min_top = slab[last_indices[coords.index(min(coords))]].frac_coords[2]
        return [(0, max(upper_fcoords)), (min_top, 1)]

    # If the entire slab region is within the cell, just
    # set the range as the highest and lowest site in the Slab
    sorted_sites = sorted(slab, key=lambda site: site.frac_coords[2])
    return [(sorted_sites[0].frac_coords[2], sorted_sites[-1].frac_coords[2])]


class SlabGenerator:
    """Generate different slabs using shift values determined by where
    a unique termination can be found, along with other criteria such as where a
    termination doesn't break a polyhedral bond. The shift value then indicates
    where the slab layer will begin and terminate in the slab-vacuum system.

    Attributes:
        oriented_unit_cell (Structure): An oriented unit cell of the parent structure.
        parent (Structure): Parent structure from which Slab was derived.
        lll_reduce (bool): Whether the slabs will be orthogonalized.
        center_slab (bool): Whether the slabs will be centered in the slab-vacuum system.
        slab_scale_factor (float): Scale factor that brings
            the parent cell to the surface cell.
        miller_index (tuple): Miller index of plane parallel to surface.
        min_slab_size (float): Minimum size of layers containing atoms, in angstroms.
        min_vac_size (float): Minimum vacuum layer size, in angstroms.
    """

    def __init__(
        self,
        initial_structure: Structure,
        miller_index: MillerIndex,
        min_slab_size: float,
        min_vacuum_size: float,
        lll_reduce: bool = False,
        center_slab: bool = False,
        in_unit_planes: bool = False,
        primitive: bool = True,
        max_normal_search: int | None = None,
        reorient_lattice: bool = True,
    ) -> None:
        """Calculate the slab scale factor and uses it to generate an
        oriented unit cell (OUC) of the initial structure.
        Also stores the initial information needed later on to generate a slab.

        Args:
            initial_structure (Structure): Initial input structure. Note that to
                ensure that the Miller indices correspond to usual
                crystallographic definitions, you should supply a conventional
                unit cell structure.
            miller_index ([h, k, l]): Miller index of the plane parallel to
                the surface. Note that this is referenced to the input structure.
                If you need this to be based on the conventional cell,
                you should supply the conventional structure.
            min_slab_size (float): In Angstroms or number of hkl planes
            min_vacuum_size (float): In Angstroms or number of hkl planes
            lll_reduce (bool): Whether to perform an LLL reduction on the
                final structure.
            center_slab (bool): Whether to center the slab in the cell with
                equal vacuum spacing from the top and bottom.
            in_unit_planes (bool): Whether to set min_slab_size and min_vac_size
                in number of hkl planes or Angstrom (default).
                Setting in units of planes is useful to ensure some slabs
                to have a certain number of layers, e.g. for Cs(100), 10 Ang
                will result in a slab with only 2 layers, whereas
                Fe(100) will have more layers. The slab thickness
                will be in min_slab_size/math.ceil(self._proj_height/dhkl)
                multiples of oriented unit cells.
            primitive (bool): Whether to reduce generated slabs to
                primitive cell. Note this does NOT generate a slab
                from a primitive cell, it means that after slab
                generation, we attempt to reduce the generated slab to
                primitive cell.
            max_normal_search (int): If set to a positive integer, the code
                will search for a normal lattice vector that is as
                perpendicular to the surface as possible, by considering
                multiple linear combinations of lattice vectors up to
                this value. This has no bearing on surface energies,
                but may be useful as a preliminary step to generate slabs
                for absorption or other sizes. It may not be the smallest possible
                cell for simulation. Normality is not guaranteed, but the oriented
                cell will have the c vector as normal as possible to the surface.
                The max absolute Miller index is usually sufficient.
            reorient_lattice (bool): reorient the lattice such that
                the c direction is parallel to the third lattice vector
        """

        def reduce_vector(vector: MillerIndex) -> MillerIndex:
            """Helper function to reduce vectors."""
            divisor = abs(reduce(gcd, vector))  # type: ignore[arg-type]
            return cast(Tuple3Ints, tuple(int(idx / divisor) for idx in vector))

        def add_site_types() -> None:
            """Add Wyckoff symbols and equivalent sites to the initial structure."""
            if (
                "bulk_wyckoff" not in initial_structure.site_properties
                or "bulk_equivalent" not in initial_structure.site_properties
            ):
                spg_analyzer = SpacegroupAnalyzer(initial_structure)
                initial_structure.add_site_property("bulk_wyckoff", spg_analyzer.get_symmetry_dataset().wyckoffs)
                initial_structure.add_site_property(
                    "bulk_equivalent",
                    spg_analyzer.get_symmetry_dataset().equivalent_atoms.tolist(),
                )

        def calculate_surface_normal() -> np.ndarray:
            """Calculate the unit surface normal vector using the reciprocal
            lattice vector.
            """
            recip_lattice = lattice.reciprocal_lattice_crystallographic

            normal = recip_lattice.get_cartesian_coords(miller_index)
            normal /= np.linalg.norm(normal)
            return normal

        def calculate_scaling_factor() -> np.ndarray:
            """Calculate scaling factor.

            # TODO (@DanielYang59): revise docstring to add more details.
            """
            slab_scale_factor = []
            non_orth_ind = []
            eye = np.eye(3, dtype=np.int64)
            for idx, miller_idx in enumerate(miller_index):
                if miller_idx == 0:
                    # If lattice vector is perpendicular to surface normal, i.e.,
                    # in plane of surface. We will simply choose this lattice
                    # vector as the basis vector
                    slab_scale_factor.append(eye[idx])

                else:
                    # Calculate projection of lattice vector onto surface normal.
                    d = abs(np.dot(normal, lattice.matrix[idx])) / lattice.abc[idx]
                    non_orth_ind.append((idx, d))

            # We want the vector that has maximum magnitude in the
            # direction of the surface normal as the c-direction.
            # Results in a more "orthogonal" unit cell.
            c_index, _dist = max(non_orth_ind, key=lambda t: t[1])

            if len(non_orth_ind) > 1:
                lcm_miller = lcm(*(miller_index[i] for i, _d in non_orth_ind))
                for (ii, _di), (jj, _dj) in itertools.combinations(non_orth_ind, 2):
                    scale_factor = [0, 0, 0]
                    scale_factor[ii] = -int(round(lcm_miller / miller_index[ii]))
                    scale_factor[jj] = int(round(lcm_miller / miller_index[jj]))
                    slab_scale_factor.append(scale_factor)
                    if len(slab_scale_factor) == 2:
                        break

            if max_normal_search is None:
                slab_scale_factor.append(eye[c_index])
            else:
                index_range = sorted(
                    range(-max_normal_search, max_normal_search + 1),
                    key=lambda x: -abs(x),
                )
                candidates = []
                for uvw in itertools.product(index_range, index_range, index_range):
                    if (not any(uvw)) or abs(np.linalg.det([*slab_scale_factor, uvw])) < 1e-8:
                        continue
                    vec = lattice.get_cartesian_coords(uvw)
                    osdm = np.linalg.norm(vec)
                    cosine = abs(np.dot(vec, normal) / osdm)
                    candidates.append((uvw, cosine, osdm))
                    # Stop searching if cosine equals 1 or -1
                    if isclose(abs(cosine), 1, abs_tol=1e-8):
                        break
                # We want the indices with the maximum absolute cosine,
                # but smallest possible length.
                uvw, cosine, osdm = max(candidates, key=lambda x: (x[1], -x[2]))
                slab_scale_factor.append(uvw)

            slab_scale_factor = np.array(slab_scale_factor)

            # Let's make sure we have a left-handed crystallographic system
            if np.linalg.det(slab_scale_factor) < 0:
                slab_scale_factor *= -1

            # Make sure the slab_scale_factor is reduced to avoid
            # unnecessarily large slabs
            reduced_scale_factor = [reduce_vector(v) for v in slab_scale_factor]
            return np.array(reduced_scale_factor)

        # Add Wyckoff symbols and equivalent sites to the initial structure,
        # to help identify types of sites in the generated slab
        add_site_types()

        # Calculate the surface normal
        lattice = initial_structure.lattice
        miller_index = reduce_vector(miller_index)
        normal = calculate_surface_normal()

        # Calculate scale factor
        slab_scale_factor = calculate_scaling_factor()

        single = initial_structure.copy()
        single.make_supercell(slab_scale_factor)

        # Calculate the most reduced structure as OUC to minimize calculations
        self.oriented_unit_cell = Structure.from_sites(single, to_unit_cell=True)

        self.max_normal_search = max_normal_search
        self.parent = initial_structure
        self.lll_reduce = lll_reduce
        self.center_slab = center_slab
        self.slab_scale_factor = slab_scale_factor
        self.miller_index = miller_index
        self.min_vac_size = min_vacuum_size
        self.min_slab_size = min_slab_size
        self.in_unit_planes = in_unit_planes
        self.primitive = primitive
        self._normal = normal  # TODO (@DanielYang59): used only in unit test
        self.reorient_lattice = reorient_lattice

        _a, _b, c = self.oriented_unit_cell.lattice.matrix
        self._proj_height = abs(np.dot(normal, c))

    def get_slab(
        self,
        shift: float = 0,
        tol: float = 0.1,
        energy: float | None = None,
    ) -> Slab:
        """[Private method] Generate a slab based on a given termination
            coordinate along the lattice c direction.

        You should RARELY use this method directly.

        Args:
            shift (float): The termination coordinate along the lattice c
                direction in fractional coordinates.
            tol (float): Tolerance to determine primitive cell.
            energy (float): The energy to assign to the slab.

        Returns:
            Slab: from a shifted oriented unit cell.
        """
        # Calculate total number of layers
        height = self._proj_height
        height_per_layer = round(height / self.parent.lattice.d_hkl(self.miller_index), 8)

        if self.in_unit_planes:
            n_layers_slab = math.ceil(self.min_slab_size / height_per_layer)
            n_layers_vac = math.ceil(self.min_vac_size / height_per_layer)
        else:
            n_layers_slab = math.ceil(self.min_slab_size / height)
            n_layers_vac = math.ceil(self.min_vac_size / height)

        n_layers = n_layers_slab + n_layers_vac

        # Prepare for Slab generation: lattice, species, coords and site_properties
        a, b, c = self.oriented_unit_cell.lattice.matrix
        new_lattice = [a, b, n_layers * c]

        species = self.oriented_unit_cell.species_and_occu

        # Shift all atoms to the termination
        frac_coords = self.oriented_unit_cell.frac_coords
        frac_coords = np.array(frac_coords) + np.array([0, 0, -shift])[None, :]
        frac_coords -= np.floor(frac_coords)  # wrap to the [0, 1) range

        # Scale down z-coordinate by the number of layers
        frac_coords[:, 2] /= n_layers

        # Duplicate atom layers by stacking along the z-axis
        all_coords = []
        for idx in range(n_layers_slab):
            _frac_coords = frac_coords.copy()
            _frac_coords[:, 2] += idx / n_layers
            all_coords.extend(_frac_coords)

        # Scale properties by number of atom layers (excluding vacuum)
        props = self.oriented_unit_cell.site_properties
        props = {k: v * n_layers_slab for k, v in props.items()}

        # Generate Slab
        struct: Structure = Structure(new_lattice, species * n_layers_slab, all_coords, site_properties=props)

        # (Optionally) Post-process the Slab
        # Orthogonalize the structure (through LLL lattice basis reduction)
        scale_factor = self.slab_scale_factor
        if self.lll_reduce:
            # Sanitize Slab (LLL reduction + site sorting + map frac_coords)
            lll_slab = struct.copy(sanitize=True)

            # Apply reduction on the scaling factor
            mapping = lll_slab.lattice.find_mapping(struct.lattice)
            struct = lll_slab
            if mapping is None:
                raise RuntimeError("LLL reduction has failed")
            scale_factor = np.dot(mapping[2], scale_factor)

        # Center the slab layer around the vacuum
        if self.center_slab:
            struct = center_slab(struct)

        # Reduce to primitive cell
        if self.primitive:
            prim_slab = struct.get_primitive_structure(tolerance=tol)
            struct = prim_slab

            if energy is not None:
                energy *= prim_slab.volume / struct.volume

        # Reorient the lattice to get the correctly reduced cell
        ouc = self.oriented_unit_cell.copy()
        if self.primitive:
            # Find a reduced OUC
            slab_l = struct.lattice
            ouc = ouc.get_primitive_structure(
                constrain_latt={
                    "a": slab_l.a,
                    "b": slab_l.b,
                    "alpha": slab_l.alpha,
                    "beta": slab_l.beta,
                    "gamma": slab_l.gamma,
                }
            )

            # Ensure lattice a and b are consistent between the OUC and the Slab
            ouc = ouc if (slab_l.a == ouc.lattice.a and slab_l.b == ouc.lattice.b) else self.oriented_unit_cell

        return Slab(
            struct.lattice,
            struct.species_and_occu,
            struct.frac_coords,
            self.miller_index,
            ouc,
            shift,
            scale_factor,
            reorient_lattice=self.reorient_lattice,
            site_properties=struct.site_properties,
            energy=energy,
        )

    def get_slabs(
        self,
        bonds: dict[tuple[Species | Element, Species | Element], float] | None = None,
        ftol: float = 0.1,
        tol: float = 0.1,
        max_broken_bonds: int = 0,
        symmetrize: bool = False,
        repair: bool = False,
        ztol: float = 0,
        filter_out_sym_slabs: bool = True,
    ) -> list[Slab]:
        """Generate slabs with shift values calculated from the internal
        gen_possible_terminations func. If the user decide to avoid breaking
        any polyhedral bond (by setting `bonds`), any shift value that do so
        would be filtered out.

        Args:
            bonds (dict): A {(species1, species2): max_bond_dist} dict.
                For example, PO4 groups may be defined as {("P", "O"): 3}.
            tol (float): Fractional tolerance for getting primitive cells
                and matching structures.
            ftol (float): Threshold for fcluster to check if two atoms are
                on the same plane. Default to 0.1 Angstrom in the direction of
                the surface normal.
            max_broken_bonds (int): Maximum number of allowable broken bonds
                for the slab. Use this to limit number of slabs. Defaults to 0,
                which means no bonds could be broken.
            symmetrize (bool): Whether to enforce the equivalency of slab surfaces.
            repair (bool): Whether to repair terminations with broken bonds (True)
                or just omit them (False). Default to False as repairing terminations
                can lead to many more possible slabs.
            ztol (float): Fractional tolerance for determine overlapping z-ranges,
                smaller ztol might result in more possible Slabs.
            filter_out_sym_slabs (bool): If True filter out identical slabs with different terminations.

        Returns:
            list[Slab]: All possible Slabs of a particular surface,
                sorted by the number of bonds broken.
        """

        def gen_possible_terminations(ftol: float) -> list[float]:
            """Generate possible terminations by clustering z coordinates.

            Args:
                ftol (float): Threshold for fcluster to check if
                    two atoms are on the same plane.
            """
            frac_coords = self.oriented_unit_cell.frac_coords
            n_atoms: int = len(frac_coords)

            # Skip clustering when there is only one atom
            if n_atoms == 1:
                # Put the atom to the center
                termination = frac_coords[0][2] + 0.5
                return [termination - math.floor(termination)]

            # Compute a Cartesian z-coordinate distance matrix
            # TODO (@DanielYang59): account for periodic boundary condition
            dist_matrix: NDArray = np.zeros((n_atoms, n_atoms))
            for i, j in itertools.combinations(list(range(n_atoms)), 2):
                if i != j:
                    z_dist = frac_coords[i][2] - frac_coords[j][2]
                    z_dist = abs(z_dist - round(z_dist)) * self._proj_height
                    dist_matrix[i, j] = z_dist
                    dist_matrix[j, i] = z_dist

            # Cluster the sites by z coordinates
            z_matrix = linkage(squareform(dist_matrix))
            clusters = fcluster(z_matrix, ftol, criterion="distance")

            # Generate cluster to z-coordinate mapping
            clst_loc: dict[Any, float] = {clst: frac_coords[idx][2] for idx, clst in enumerate(clusters)}

            # Wrap all clusters into the unit cell ([0, 1) range)
            possible_clst: list[float] = [coord - math.floor(coord) for coord in sorted(clst_loc.values())]

            # Calculate terminations
            n_terms: int = len(possible_clst)
            terminations: list[float] = []
            for idx in range(n_terms):
                # Handle the special case for the first-last pair of
                # z coordinates (because of periodic boundary condition)
                if idx == n_terms - 1:
                    termination = (possible_clst[0] + 1 + possible_clst[idx]) * 0.5
                else:
                    termination = (possible_clst[idx] + possible_clst[idx + 1]) * 0.5

                # Wrap termination to [0, 1) range
                terminations.append(termination - math.floor(termination))

            return sorted(terminations)

        def get_z_ranges(
            bonds: dict[tuple[Species | Element, Species | Element], float],
            ztol: float,
        ) -> list[tuple[float, float]]:
            """Collect occupied z ranges where each range is a (lower_z, upper_z) tuple.

            This method examines all sites in the oriented unit cell (OUC)
            and considers all neighboring sites within the specified bond distance
            for each site. If a site and its neighbor meet bonding and species
            requirements, their respective z-ranges will be collected.

            Args:
                bonds (dict): A {(species1, species2): max_bond_dist} dict.
                ztol (float): Fractional tolerance for determine overlapping z-ranges.
            """
            # Sanitize species in dict keys
            bonds = {(get_el_sp(s1), get_el_sp(s2)): dist for (s1, s2), dist in bonds.items()}

            z_ranges = []
            for (sp1, sp2), bond_dist in bonds.items():
                for site in self.oriented_unit_cell:
                    if sp1 in site.species:
                        for nn in self.oriented_unit_cell.get_neighbors(site, bond_dist):
                            if sp2 in nn.species:
                                z_range = tuple(sorted([site.frac_coords[2], nn.frac_coords[2]]))

                                # Handle cases when z coordinate of site goes
                                # beyond the upper boundary
                                if z_range[1] > 1:
                                    z_ranges.extend([(z_range[0], 1), (0, z_range[1] - 1)])

                                # When z coordinate is below the lower boundary
                                elif z_range[0] < 0:
                                    z_ranges.extend([(0, z_range[1]), (z_range[0] + 1, 1)])

                                # Neglect overlapping positions
                                elif not isclose(z_range[0], z_range[1], abs_tol=ztol):
                                    z_ranges.append(z_range)

            return z_ranges

        # Get occupied z_ranges
        z_ranges = [] if bonds is None else get_z_ranges(bonds, ztol)

        slabs = []
        for termination in gen_possible_terminations(ftol=ftol):
            # Calculate total number of bonds broken (how often the
            # termination fall within the z_range occupied by a bond)
            bonds_broken = 0
            for z_range in z_ranges:
                if z_range[0] <= termination <= z_range[1]:
                    bonds_broken += 1

            # DEBUG(@DanielYang59): number of bonds broken passed to energy
            # As per the docstring this is to sort final Slabs by number
            # of bonds broken, but this may very likely lead to errors
            # if the "energy" is used literally (Maybe reset energy to None?)
            slab = self.get_slab(shift=termination, tol=tol, energy=bonds_broken)

            if bonds_broken <= max_broken_bonds:
                slabs.append(slab)

            # If the number of broken bonds is exceeded, repair the broken bonds
            elif repair and bonds is not None:
                slabs.append(self.repair_broken_bonds(slab=slab, bonds=bonds))

        # Filter out surfaces that might be the same
        if filter_out_sym_slabs:
            matcher = StructureMatcher(ltol=tol, stol=tol, primitive_cell=False, scale=False)

            final_slabs: list[Slab] = []
            for group in matcher.group_structures(slabs):
                # For each unique slab, symmetrize the
                # surfaces by removing sites from the bottom
                if symmetrize:
                    sym_slabs = self.nonstoichiometric_symmetrized_slab(group[0])
                    final_slabs.extend(sym_slabs)
                else:
                    final_slabs.append(group[0])

            # Filter out similar surfaces generated by symmetrization
            if symmetrize:
                matcher_sym = StructureMatcher(ltol=tol, stol=tol, primitive_cell=False, scale=False)
                final_slabs = [group[0] for group in matcher_sym.group_structures(final_slabs)]
        else:
            final_slabs = slabs

        return cast(list[Slab], sorted(final_slabs, key=lambda slab: slab.energy))

    def repair_broken_bonds(
        self,
        slab: Slab,
        bonds: dict[tuple[Species | Element, Species | Element], float],
    ) -> Slab:
        """Repair broken bonds (specified by the bonds parameter) due to
        slab cleaving, and repair them by moving undercoordinated atoms
        to the other surface.

        How it works:
            For example a P-O4 bond may have P and O(4-x) on one side
            of the surface, and Ox on the other side, this method would
            first move P (the reference atom) to the other side,
            find its missing nearest neighbours (Ox), and move P
            and Ox back together.

        Args:
            slab (Slab): The Slab to repair.
            bonds (dict): A {(species1, species2): max_bond_dist} dict.
                For example, PO4 groups may be defined as {("P", "O"): 3}.

        Returns:
            Slab: The repaired Slab.
        """
        for species_pair, bond_dist in bonds.items():
            # Determine which element should be the reference (center)
            # element for determining broken bonds, e.g. P for PO4 bond.
            cn_dict = {}
            for idx, ele in enumerate(species_pair):
                cn_list = []
                for site in self.oriented_unit_cell:
                    # Find integer coordination numbers for element pairs
                    ref_cn = 0
                    if site.species_string == ele:
                        for nn in self.oriented_unit_cell.get_neighbors(site, bond_dist):
                            if nn[0].species_string == species_pair[idx - 1]:
                                ref_cn += 1

                    cn_list.append(ref_cn)
                cn_dict[ele] = cn_list

            # Make the element with higher coordination the reference
            if max(cn_dict[species_pair[0]]) > max(cn_dict[species_pair[1]]):
                ele_ref, ele_other = species_pair
            else:
                ele_other, ele_ref = species_pair

            for idx, site in enumerate(slab):
                # Determine the coordination of the reference
                if site.species_string == ele_ref:
                    ref_cn = sum(
                        1 if neighbor.species_string == ele_other else 0
                        for neighbor in slab.get_neighbors(site, bond_dist)
                    )

                    # Suppose we find an undercoordinated reference atom
                    # TODO (@DanielYang59): maybe use the following to
                    # check if the reference atom is "undercoordinated"
                    # if ref_cn < min(cn_dict[ele_ref]):
                    if ref_cn not in cn_dict[ele_ref]:
                        # Move this reference atom to the other side
                        slab = self.move_to_other_side(slab, [idx])

                        # Find its NNs (with right species) it should bond to
                        neighbors = slab.get_neighbors(slab[idx], r=bond_dist)
                        to_move = [nn[2] for nn in neighbors if nn[0].species_string == ele_other]
                        to_move.append(idx)

                        # Move those NNs along with the reference
                        # atom back to the other side of the slab
                        slab = self.move_to_other_side(slab, to_move)

        return slab

    def move_to_other_side(
        self,
        init_slab: Slab,
        index_of_sites: list[int],
    ) -> Slab:
        """Move surface sites to the opposite surface of the Slab.

        If a selected site resides on the top half of the Slab,
        it would be moved to the bottom side, and vice versa.
        The distance moved is equal to the thickness of the Slab.

        Note:
            You should only use this method on sites close to the
            surface, otherwise it would end up deep inside the
            vacuum layer.

        Args:
            init_slab (Slab): The Slab whose sites would be moved.
            index_of_sites (list[int]): Indices representing
                the sites to move.

        Returns:
            Slab: The Slab with selected sites moved.
        """
        # Calculate Slab height
        height: float = self._proj_height
        # Scale height if using number of hkl planes
        if self.in_unit_planes:
            height /= self.parent.lattice.d_hkl(self.miller_index)

        # Calculate the moving distance as the fractional height
        # of the Slab inside the cell
        # DEBUG(@DanielYang59): use actual sizes for slab/vac
        # instead of the input arg (min_slab/vac_size)
        n_layers_slab: int = math.ceil(self.min_slab_size / height)
        n_layers_vac: int = math.ceil(self.min_vac_size / height)
        n_layers: int = n_layers_slab + n_layers_vac

        frac_dist: float = n_layers_slab / n_layers

        # Separate selected sites into top and bottom
        top_site_index: list[int] = []
        bottom_site_index: list[int] = []
        for idx in index_of_sites:
            if init_slab[idx].frac_coords[2] >= init_slab.center_of_mass[2]:
                top_site_index.append(idx)
            else:
                bottom_site_index.append(idx)

        # Move sites to the opposite surface
        slab = init_slab.copy()
        slab.translate_sites(top_site_index, vector=[0, 0, -frac_dist], frac_coords=True)
        slab.translate_sites(bottom_site_index, vector=[0, 0, frac_dist], frac_coords=True)

        return Slab(
            init_slab.lattice,
            slab.species,
            slab.frac_coords,
            init_slab.miller_index,
            init_slab.oriented_unit_cell,
            init_slab.shift,
            init_slab.scale_factor,
            energy=init_slab.energy,
        )

    def nonstoichiometric_symmetrized_slab(self, init_slab: Slab) -> list[Slab]:
        """Symmetrize the two surfaces of a Slab, but may break the stoichiometry.

        How it works:
            1. Check whether two surfaces of the slab are equivalent.
            If the point group of the slab has an inversion symmetry (
            ie. belong to one of the Laue groups), then it's assumed that the
            surfaces are equivalent.

            2.If not symmetrical, sites at the bottom of the slab will be removed
            until the slab is symmetric, which may break the stoichiometry.

        Args:
            init_slab (Slab): The initial Slab.

        Returns:
            list[Slabs]: The symmetrized Slabs.
        """
        if init_slab.is_symmetric():
            return [init_slab]

        non_stoich_slabs = []
        # Build a symmetrical surface slab for each of the different surfaces
        for surface in ("top", "bottom"):
            is_sym: bool = False
            slab = init_slab.copy()
            slab.energy = init_slab.energy

            while not is_sym:
                # Keep removing sites from the bottom until surfaces are
                # symmetric or the number of sites removed has
                # exceeded 10 percent of the original slab
                # TODO: (@DanielYang59) comment differs from implementation:
                # no "exceeded 10 percent" check
                z_coords: list[float] = [site[2] for site in slab.frac_coords]

                if surface == "top":
                    slab.remove_sites([z_coords.index(max(z_coords))])
                else:
                    slab.remove_sites([z_coords.index(min(z_coords))])

                if len(slab) <= len(self.parent):
                    warnings.warn("Too many sites removed, please use a larger slab.")
                    break

                # Check if the new Slab is symmetric
                # TODO: (@DanielYang59): should have some feedback (warning)
                # if cannot symmetrize the Slab
                if slab.is_symmetric():
                    is_sym = True
                    non_stoich_slabs.append(slab)

        return non_stoich_slabs


def generate_all_slabs(
    structure: Structure,
    max_index: int,
    min_slab_size: float,
    min_vacuum_size: float,
    bonds: dict | None = None,
    tol: float = 0.1,
    ftol: float = 0.1,
    max_broken_bonds: int = 0,
    lll_reduce: bool = False,
    center_slab: bool = False,
    primitive: bool = True,
    max_normal_search: int | None = None,
    symmetrize: bool = False,
    repair: bool = False,
    include_reconstructions: bool = False,
    in_unit_planes: bool = False,
) -> list[Slab]:
    """Find all unique Slabs up to a given Miller index.

    Slabs oriented along certain Miller indices may be equivalent to
    other Miller indices under symmetry operations. To avoid
    duplication, such equivalent slabs would be filtered out.
    For instance, CsCl has equivalent slabs in the (0,0,1),
    (0,1,0), and (1,0,0) directions under symmetry operations.

    Args:
        structure (Structure): Initial input structure. To
            ensure that the Miller indices correspond to usual
            crystallographic definitions, you should supply a
            conventional unit cell.
        max_index (int): The maximum Miller index to go up to.
        min_slab_size (float): The minimum slab size in Angstrom.
        min_vacuum_size (float): The minimum vacuum layer thickness in Angstrom.
        bonds (dict): A {(species1, species2): max_bond_dist} dict.
                For example, PO4 groups may be defined as {("P", "O"): 3}.
        tol (float): Tolerance for getting primitive cells and
            matching structures.
        ftol (float): Tolerance in Angstrom for fcluster to check
            if two atoms are on the same plane. Default to 0.1 Angstrom
            in the direction of the surface normal.
        max_broken_bonds (int): Maximum number of allowable broken bonds
            for the slab. Use this to limit the number of slabs.
            Defaults to zero, which means no bond can be broken.
        lll_reduce (bool): Whether to perform an LLL reduction on the
            final Slab.
        center_slab (bool): Whether to center the slab in the cell with
            equal vacuum spacing from the top and bottom.
        primitive (bool): Whether to reduce generated slabs to
            primitive cell. Note this does NOT generate a slab
            from a primitive cell, it means that after slab
            generation, we attempt to reduce the generated slab to
            primitive cell.
        max_normal_search (int): If set to a positive integer, the code
            will search for a normal lattice vector that is as
            perpendicular to the surface as possible, by considering
            multiple linear combinations of lattice vectors up to
            this value. This has no bearing on surface energies,
            but may be useful as a preliminary step to generate slabs
            for absorption or other sizes. It may not be the smallest possible
            cell for simulation. Normality is not guaranteed, but the oriented
            cell will have the c vector as normal as possible to the surface.
            The max absolute Miller index is usually sufficient.
        symmetrize (bool): Whether to ensure the surfaces of the
            slabs are equivalent.
        repair (bool): Whether to repair terminations with broken bonds
            or just omit them.
        include_reconstructions (bool): Whether to include reconstructed
            slabs available in the reconstructions_archive.json file. Defaults to False.
        in_unit_planes (bool): Whether to set min_slab_size and min_vac_size
            in number of hkl planes or Angstrom (default).
            Setting in units of planes is useful to ensure some slabs
            to have a certain number of layers, e.g. for Cs(100), 10 Ang
            will result in a slab with only 2 layers, whereas
            Fe(100) will have more layers. The slab thickness
            will be in min_slab_size/math.ceil(self._proj_height/dhkl)
            multiples of oriented unit cells.
    """
    all_slabs: list[Slab] = []

    for miller in get_symmetrically_distinct_miller_indices(structure, max_index):
        gen = SlabGenerator(
            structure,
            miller,
            min_slab_size,
            min_vacuum_size,
            lll_reduce=lll_reduce,
            center_slab=center_slab,
            primitive=primitive,
            max_normal_search=max_normal_search,
            in_unit_planes=in_unit_planes,
        )
        slabs = gen.get_slabs(
            bonds=bonds,
            tol=tol,
            ftol=ftol,
            symmetrize=symmetrize,
            max_broken_bonds=max_broken_bonds,
            repair=repair,
        )

        if len(slabs) > 0:
            logger.debug(f"{miller} has {len(slabs)} slabs... ")
            all_slabs.extend(slabs)

    if include_reconstructions:
        symbol = SpacegroupAnalyzer(structure).get_space_group_symbol()
        # Enumerate through all reconstructions in the
        # archive available for this particular spacegroup
        for name, instructions in RECONSTRUCTIONS_ARCHIVE.items():
            if "base_reconstruction" in instructions:
                instructions = RECONSTRUCTIONS_ARCHIVE[instructions["base_reconstruction"]]

            if instructions["spacegroup"]["symbol"] == symbol:
                # Make sure this reconstruction has a max index
                # equal or less than the given max index
                if max(instructions["miller_index"]) > max_index:
                    continue
                recon = ReconstructionGenerator(structure, min_slab_size, min_vacuum_size, name)
                all_slabs.extend(recon.build_slabs())

    return all_slabs


# Load the reconstructions_archive JSON file
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
with open(f"{MODULE_DIR}/reconstructions_archive.json", encoding="utf-8") as data_file:
    RECONSTRUCTIONS_ARCHIVE = json.load(data_file)


def get_d(slab: Slab) -> float:
    """Determine the z-spacing between the bottom two layers for a Slab."""
    # Sort all sites by z-coordinates
    sorted_sites = sorted(slab, key=lambda site: site.frac_coords[2])

    distance = None
    for site, next_site in itertools.pairwise(sorted_sites):
        if not isclose(site.frac_coords[2], next_site.frac_coords[2], abs_tol=1e-6):
            distance = next_site.frac_coords[2] - site.frac_coords[2]
            break

    if distance is None:
        raise RuntimeError("Cannot identify any layer.")

    return slab.lattice.get_cartesian_coords([0, 0, distance])[2]


class ReconstructionGenerator:
    """Build a reconstructed Slab from a given initial Structure.

    This class needs a pre-defined dictionary specifying the parameters
    needed such as the SlabGenerator parameters, transformation matrix,
    sites to remove/add and slab/vacuum sizes.

    Attributes:
        slabgen_params (dict): Parameters for the SlabGenerator.
        trans_matrix (np.ndarray): A 3x3 transformation matrix to generate
            the reconstructed slab. Only the a and b lattice vectors are
            actually changed while the c vector remains the same.
            This matrix is what the Wood's notation is based on.
        reconstruction_json (dict): The full JSON or dictionary containing
            the instructions for building the slab.

    Todo:
        - Right now there is no way to specify what atom is being added.
            Use basis sets in the future?
    """

    def __init__(
        self,
        initial_structure: Structure,
        min_slab_size: float,
        min_vacuum_size: float,
        reconstruction_name: str,
    ) -> None:
        """Generate reconstructed slabs from a set of instructions.

        Args:
            initial_structure (Structure): Initial input structure. Note
                that to ensure that the Miller indices correspond to usual
                crystallographic definitions, you should supply a conventional
                unit cell structure.
            min_slab_size (float): Minimum Slab size in Angstrom.
            min_vacuum_size (float): Minimum vacuum layer size in Angstrom.
            reconstruction_name (str): Name of the dict containing the build
                instructions. The dictionary can contain any item, however
                any instructions archived in pymatgen for public use need
                to contain the following keys and items to ensure
                compatibility with the ReconstructionGenerator:

                    "name" (str): A descriptive name for the reconstruction,
                        typically including the type of structure,
                        the Miller index, the Wood's notation and additional
                        descriptors for the reconstruction.
                        Example: "fcc_110_missing_row_1x2"
                    "description" (str): A detailed description of the
                        reconstruction, intended to assist future contributors
                        in avoiding duplicate entries. Please read the description
                        carefully before adding to prevent duplications.
                    "reference" (str): Optional reference to the source of
                        the reconstruction.
                    "spacegroup" (dict): A dictionary indicating the space group
                        of the reconstruction. e.g. {"symbol": "Fm-3m", "number": 225}.
                    "miller_index" ([h, k, l]): Miller index of the reconstruction
                    "Woods_notation" (str): For a reconstruction, the a and b
                        lattice may change to accommodate the symmetry.
                        This notation indicates the change in
                        the vectors relative to the primitive (p) or
                        conventional (c) slab cell. E.g. p(2x1).

                        Reference: Wood, E. A. (1964). Vocabulary of surface
                        crystallography. Journal of Applied Physics, 35(4),
                        1306-1312.
                    "transformation_matrix" (numpy array): A 3x3 matrix to
                        transform the slab. Only the a and b lattice vectors
                        should change while the c vector remains the same.
                    "SlabGenerator_parameters" (dict): A dictionary containing
                        the parameters for the SlabGenerator, excluding the
                        miller_index, min_slab_size and min_vac_size. As the
                        Miller index is already specified and the min_slab_size
                        and min_vac_size can be changed regardless of the
                        reconstruction type. Having a consistent set of
                        SlabGenerator parameters allows for the instructions to
                        be reused.
                    "points_to_remove" (list[site]): A list of sites to
                        remove where the first two indices are fractional (in a
                        and b) and the third index is in units of 1/d (in c),
                        see the below "Notes" for details.
                    "points_to_add" (list[site]): A list of sites to add
                        where the first two indices are fractional (in a an b) and
                        the third index is in units of 1/d (in c), see the below
                        "Notes" for details.
                    "base_reconstruction" (dict, Optional): A dictionary specifying
                        an existing reconstruction model upon which the current
                        reconstruction is built to avoid repetition. E.g. the
                        alpha reconstruction of halites is based on the octopolar
                        reconstruction but with the topmost atom removed. The dictionary
                        for the alpha reconstruction would therefore contain the item
                        "reconstruction_base": "halite_111_octopolar_2x2", and
                        additional sites can be added by "points_to_add".

        Notes:
            1. For "points_to_remove" and "points_to_add", the third index
                for the c vector is specified in units of 1/d, where d represents
                the spacing between atoms along the hkl (the c vector), relative
                to the topmost site in the unreconstructed slab. For instance,
                a point of [0.5, 0.25, 1] corresponds to the 0.5 fractional
                coordinate of a, 0.25 fractional coordinate of b, and a
                distance of 1 atomic layer above the topmost site. Similarly,
                [0.5, 0.25, -0.5] corresponds to a point half an atomic layer
                below the topmost site, and [0.5, 0.25, 0] corresponds to a
                point at the same position along c as the topmost site.
                This approach is employed because while the primitive units
                of a and b remain constant, the user can vary the length
                of the c direction by adjusting the slab layer or the vacuum layer.

            2. The dictionary should only provide "points_to_remove" and
                "points_to_add" for the top surface. The ReconstructionGenerator
                will modify the bottom surface accordingly to return a symmetric Slab.
        """

        def build_recon_json() -> dict:
            """Build reconstruction instructions, optionally upon a base instruction set."""
            # Check if reconstruction instruction exists
            # TODO (@DanielYang59): can we avoid asking user to modify the source file?
            if reconstruction_name not in RECONSTRUCTIONS_ARCHIVE:
                raise KeyError(
                    f"{reconstruction_name=} does not exist in the archive. "
                    "Please select from one of the following: "
                    f"{list(RECONSTRUCTIONS_ARCHIVE)} or add it to the "
                    "archive file 'reconstructions_archive.json'."
                )

            # Get the reconstruction instructions from the archive file
            recon_json: dict = copy.deepcopy(RECONSTRUCTIONS_ARCHIVE[reconstruction_name])

            # Build new instructions from a base reconstruction
            if "base_reconstruction" in recon_json:
                new_points_to_add: list = []
                new_points_to_remove: list = []

                if "points_to_add" in recon_json:
                    new_points_to_add = recon_json["points_to_add"]
                if "points_to_remove" in recon_json:
                    new_points_to_remove = recon_json["points_to_remove"]

                # DEBUG (@DanielYang59): the following overwrites previously
                # loaded "recon_json", use condition to avoid this
                recon_json = copy.deepcopy(RECONSTRUCTIONS_ARCHIVE[recon_json["base_reconstruction"]])

                # TODO (@DanielYang59): use "site" over "point" for consistency?
                if "points_to_add" in recon_json:
                    del recon_json["points_to_add"]
                if new_points_to_add:
                    recon_json["points_to_add"] = new_points_to_add

                if "points_to_remove" in recon_json:
                    del recon_json["points_to_remove"]
                if new_points_to_remove:
                    recon_json["points_to_remove"] = new_points_to_remove

            return recon_json

        def build_slabgen_params() -> dict:
            """Build SlabGenerator parameters."""
            slabgen_params: dict = copy.deepcopy(recon_json["SlabGenerator_parameters"])
            slabgen_params["initial_structure"] = initial_structure.copy()
            slabgen_params["miller_index"] = recon_json["miller_index"]
            slabgen_params["min_slab_size"] = min_slab_size
            slabgen_params["min_vacuum_size"] = min_vacuum_size

            return slabgen_params

        # Build reconstruction instructions
        recon_json = build_recon_json()

        # Build SlabGenerator parameters
        slabgen_params = build_slabgen_params()

        self.name = reconstruction_name
        self.slabgen_params = slabgen_params
        self.reconstruction_json = recon_json
        self.trans_matrix = recon_json["transformation_matrix"]

    def build_slabs(self) -> list[Slab]:
        """Build reconstructed Slabs by:
            (1) Obtaining the unreconstructed Slab using the specified
                parameters for the SlabGenerator.
            (2) Applying the appropriate lattice transformation to the
                a and b lattice vectors.
            (3) Remove and then add specified sites from both surfaces.

        Returns:
            list[Slab]: The reconstructed slabs.
        """
        slabs = self.get_unreconstructed_slabs()

        recon_slabs = []

        for slab in slabs:
            z_spacing = get_d(slab)
            top_site = max(slab, key=lambda site: site.frac_coords[2]).coords

            # Remove specified sites
            if "points_to_remove" in self.reconstruction_json:
                sites_to_rm: list = copy.deepcopy(self.reconstruction_json["points_to_remove"])
                for site in sites_to_rm:
                    site[2] = slab.lattice.get_fractional_coords(
                        [top_site[0], top_site[1], top_site[2] + site[2] * z_spacing]
                    )[2]

                    # Find and remove nearest site
                    cart_point = slab.lattice.get_cartesian_coords(site)
                    distances: list[float] = [site.distance_from_point(cart_point) for site in slab]
                    nearest_site = distances.index(min(distances))
                    slab.symmetrically_remove_atoms(indices=[nearest_site])

            # Add specified sites
            if "points_to_add" in self.reconstruction_json:
                sites_to_add: list = copy.deepcopy(self.reconstruction_json["points_to_add"])
                for site in sites_to_add:
                    site[2] = slab.lattice.get_fractional_coords(
                        [top_site[0], top_site[1], top_site[2] + site[2] * z_spacing]
                    )[2]
                    # TODO: see ReconstructionGenerator docstring:
                    # cannot specify species to add
                    slab.symmetrically_add_atom(species=slab[0].specie, point=site)

            slab.reconstruction = self.name
            slab.recon_trans_matrix = self.trans_matrix

            # Get the oriented unit cell with the same a*b area
            ouc = slab.oriented_unit_cell.copy()
            ouc.make_supercell(self.trans_matrix)
            slab.oriented_unit_cell = ouc
            recon_slabs.append(slab)

        return recon_slabs

    def get_unreconstructed_slabs(self) -> list[Slab]:
        """Generate the unreconstructed (super) Slabs.

        TODO (@DanielYang59): this should be a private method.
        """
        return [slab.make_supercell(self.trans_matrix) for slab in SlabGenerator(**self.slabgen_params).get_slabs()]


def get_symmetrically_equivalent_miller_indices(
    structure: Structure,
    miller_index: tuple[int, ...],
    return_hkil: bool = True,
    system: CrystalSystem | None = None,
) -> list:
    """Get indices for all equivalent sites within a given structure.
    Analysis is based on the symmetry of its reciprocal lattice.

    Args:
        structure (Structure): Structure to analyze.
        miller_index (tuple): Designates the family of Miller indices
            to find. Can be hkl or hkil for hexagonal systems.
        return_hkil (bool): Whether to return hkil (True) form of Miller
            index for hexagonal systems, or hkl (False).
        system: The crystal system of the structure.
    """
    # Convert to hkl if hkil, because in_coord_list only handles tuples of 3
    if len(miller_index) >= 3:
        _miller_index: MillerIndex = (
            miller_index[0],
            miller_index[1],
            miller_index[-1],
        )
    else:
        _miller_index = (miller_index[0], miller_index[1], miller_index[2])

    max_idx = max(np.abs(miller_index))
    idx_range = list(range(-max_idx, max_idx + 1))
    idx_range.reverse()

    # Skip crystal system analysis if already given
    if system:
        spg_analyzer = None
    else:
        spg_analyzer = SpacegroupAnalyzer(structure)
        system = spg_analyzer.get_crystal_system()

    # Get distinct hkl planes from the rhombohedral setting if trigonal
    if system == "trigonal":
        if not spg_analyzer:
            spg_analyzer = SpacegroupAnalyzer(structure)
        prim_structure = spg_analyzer.get_primitive_standard_structure()
        symm_ops = prim_structure.lattice.get_recp_symmetry_operation()

    else:
        symm_ops = structure.lattice.get_recp_symmetry_operation()

    equivalent_millers: list[Tuple3Ints] = [_miller_index]
    for miller in itertools.product(idx_range, idx_range, idx_range):
        if miller == _miller_index:
            continue

        if any(idx != 0 for idx in miller):
            if _is_in_miller_family(miller, equivalent_millers, symm_ops):
                equivalent_millers += [miller]

            # Include larger Miller indices in the family of planes
            if (
                all(max_idx > i for i in np.abs(miller))
                and not in_coord_list(equivalent_millers, miller)
                and _is_in_miller_family(max_idx * np.array(miller), equivalent_millers, symm_ops)
            ):
                equivalent_millers += [miller]

    # Convert hkl to hkil if necessary
    if return_hkil and system in {"trigonal", "hexagonal"}:
        return [(hkl[0], hkl[1], -1 * hkl[0] - hkl[1], hkl[2]) for hkl in equivalent_millers]

    return equivalent_millers


def get_symmetrically_distinct_miller_indices(
    structure: Structure,
    max_index: int,
    return_hkil: bool = False,
) -> list:
    """Find all symmetrically distinct indices below a certain max-index
    for a given structure. Analysis is based on the symmetry of the
    reciprocal lattice of the structure.

    Args:
        structure (Structure): The input structure.
        max_index (int): The maximum index. For example, 1 means that
            (100), (110), and (111) are returned for the cubic structure.
            All other indices are equivalent to one of these.
        return_hkil (bool): Whether to return hkil (True) form of Miller
            index for hexagonal systems, or hkl (False).
    """
    # Get a list of all hkls for conventional (including equivalent)
    rng = list(range(-max_index, max_index + 1))[::-1]
    conv_hkl_list = [miller for miller in itertools.product(rng, rng, rng) if any(i != 0 for i in miller)]

    # Sort by the maximum absolute values of Miller indices so that
    # low-index planes come first. This is important for trigonal systems.
    conv_hkl_list = sorted(conv_hkl_list, key=lambda x: max(np.abs(x)))

    # Get distinct hkl planes from the rhombohedral setting if trigonal
    spg_analyzer = SpacegroupAnalyzer(structure)
    if spg_analyzer.get_crystal_system() == "trigonal":
        transf = spg_analyzer.get_conventional_to_primitive_transformation_matrix()
        miller_list: list[Tuple3Ints] = [hkl_transformation(transf, hkl) for hkl in conv_hkl_list]
        prim_structure = SpacegroupAnalyzer(structure).get_primitive_standard_structure()
        symm_ops = prim_structure.lattice.get_recp_symmetry_operation()

    else:
        miller_list = conv_hkl_list
        symm_ops = structure.lattice.get_recp_symmetry_operation()

    unique_millers: list = []
    unique_millers_conv: list = []

    for idx, miller in enumerate(miller_list):
        denom = abs(reduce(gcd, miller))  # type: ignore[arg-type]
        miller = cast(Tuple3Ints, tuple(int(idx / denom) for idx in miller))
        if not _is_in_miller_family(miller, unique_millers, symm_ops):
            if spg_analyzer.get_crystal_system() == "trigonal":
                # Now we find the distinct primitive hkls using
                # the primitive symmetry operations and their
                # corresponding hkls in the conventional setting
                unique_millers.append(miller)
                denom = abs(reduce(gcd, conv_hkl_list[idx]))  # type: ignore[arg-type]
                cmiller = tuple(int(idx / denom) for idx in conv_hkl_list[idx])
                unique_millers_conv.append(cmiller)
            else:
                unique_millers.append(miller)
                unique_millers_conv.append(miller)

    if return_hkil and spg_analyzer.get_crystal_system() in {"trigonal", "hexagonal"}:
        return [(hkl[0], hkl[1], -1 * hkl[0] - hkl[1], hkl[2]) for hkl in unique_millers_conv]

    return unique_millers_conv


def _is_in_miller_family(
    miller_index: MillerIndex,
    miller_list: list[MillerIndex],
    symm_ops: list,
) -> bool:
    """Helper function to check if the given Miller index belongs
    to the same family of any index in the provided list.

    Args:
        miller_index (MillerIndex): The Miller index to analyze.
        miller_list (list): List of Miller indices.
        symm_ops (list): Symmetry operations for a lattice,
            used to define the indices family.
    """
    return any(in_coord_list(miller_list, op.operate(miller_index)) for op in symm_ops)


def hkl_transformation(
    transf: np.ndarray,
    miller_index: MillerIndex,
) -> Tuple3Ints:
    """Transform the Miller index from setting A to B with a transformation matrix.

    Args:
        transf (3x3 array): The matrix that transforms a lattice from A to B.
        miller_index (MillerIndex): The Miller index [h, k, l] to transform.
    """

    def math_lcm(a: int, b: int) -> int:
        """Calculate the least common multiple."""
        return a * b // math.gcd(a, b)

    # Convert the elements of the transformation matrix to integers
    reduced_transf = reduce(math_lcm, [int(1 / i) for i in itertools.chain(*transf) if i != 0]) * transf
    reduced_transf = reduced_transf.astype(int)

    # Perform the transformation
    transf_hkl = np.dot(reduced_transf, miller_index)
    divisor = abs(reduce(gcd, transf_hkl))  # type: ignore[arg-type]
    transf_hkl = np.array([idx // divisor for idx in transf_hkl])

    # Get positive Miller index
    if sum(idx < 0 for idx in transf_hkl) > 1:
        transf_hkl *= -1

    return tuple(transf_hkl)


def miller_index_from_sites(
    lattice: Lattice | ArrayLike,
    coords: ArrayLike,
    coords_are_cartesian: bool = True,
    round_dp: int = 4,
    verbose: bool = True,
) -> Tuple3Ints:
    """Get the Miller index of a plane, determined by a given set of coordinates.

    A minimum of 3 sets of coordinates are required. If more than 3
    coordinates are given, the plane that minimises the distance to all
    sites will be calculated.

    Args:
        lattice (matrix or Lattice): A 3x3 lattice matrix or `Lattice` object.
        coords (ArrayLike): A list or numpy array of coordinates. Can be
            Cartesian or fractional coordinates.
        coords_are_cartesian (bool, optional): Whether the coordinates are
            in Cartesian coordinates, or fractional (False).
        round_dp (int, optional): The number of decimal places to round the
            Miller index to.
        verbose (bool, optional): Whether to print warnings.

    Returns:
        tuple[int]: The Miller index.
    """
    if not isinstance(lattice, Lattice):
        lattice = Lattice(lattice)

    return lattice.get_miller_index_from_coords(
        coords,
        coords_are_cartesian=coords_are_cartesian,
        round_dp=round_dp,
        verbose=verbose,
    )
