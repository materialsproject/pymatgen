# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements representations of slabs and surfaces, as well as
algorithms for generating them. If you use this module, please consider
citing the following work::

    R. Tran, Z. Xu, B. Radhakrishnan, D. Winston, W. Sun, K. A. Persson,
    S. P. Ong, "Surface Energies of Elemental Crystals", Scientific Data,
    2016, 3:160080, doi: 10.1038/sdata.2016.80.

as well as::

    Sun, W.; Ceder, G. Efficient creation and convergence of surface slabs,
    Surface Science, 2013, 617, 53–59, doi:10.1016/j.susc.2013.05.016.
"""

import copy
import itertools
import json
import logging
import math
import os
import warnings
from functools import reduce
from math import gcd

import numpy as np
from monty.fractions import lcm
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import in_coord_list

__author__ = "Richard Tran, Wenhao Sun, Zihan Xu, Shyue Ping Ong"


logger = logging.getLogger(__name__)


class Slab(Structure):
    """
    Subclass of Structure representing a Slab. Implements additional
    attributes pertaining to slabs, but the init method does not
    actually implement any algorithm that creates a slab. This is a
    DUMMY class who's init method only holds information about the
    slab. Also has additional methods that returns other information
    about a slab such as the surface area, normal, and atom adsorption.

    Note that all Slabs have the surface normal oriented perpendicular to the a
    and b lattice vectors. This means the lattice vectors a and b are in the
    surface plane and the c vector is out of the surface plane (though not
    necessarily perpendicular to the surface).

    .. attribute:: miller_index

        Miller index of plane parallel to surface.

    .. attribute:: scale_factor

        Final computed scale factor that brings the parent cell to the
        surface cell.

    .. attribute:: shift

        The shift value in Angstrom that indicates how much this
        slab has been shifted.
    """

    def __init__(
        self,
        lattice,
        species,
        coords,
        miller_index,
        oriented_unit_cell,
        shift,
        scale_factor,
        reorient_lattice=True,
        validate_proximity=False,
        to_unit_cell=False,
        reconstruction=None,
        coords_are_cartesian=False,
        site_properties=None,
        energy=None,
    ):
        """
        Makes a Slab structure, a structure object with additional information
        and methods pertaining to slabs.

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
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface. Note that this is referenced to the input structure. If
                you need this to be based on the conventional cell,
                you should supply the conventional structure.
            oriented_unit_cell (Structure): The oriented_unit_cell from which
                this Slab is created (by scaling in the c-direction).
            shift (float): The shift in the c-direction applied to get the
                termination.
            scale_factor (np.ndarray): scale_factor Final computed scale factor
                that brings the parent cell to the surface cell.
            reorient_lattice (bool): reorients the lattice parameters such that
                the c direction is along the z axis.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            reconstruction (str): Type of reconstruction. Defaults to None if
                the slab is not reconstructed.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
            energy (float): A value for the energy.
        """
        self.oriented_unit_cell = oriented_unit_cell
        self.miller_index = tuple(miller_index)
        self.shift = shift
        self.reconstruction = reconstruction
        self.scale_factor = np.array(scale_factor)
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

    def get_orthogonal_c_slab(self):
        """
        This method returns a Slab where the normal (c lattice vector) is
        "forced" to be exactly orthogonal to the surface a and b lattice
        vectors. **Note that this breaks inherent symmetries in the slab.**
        It should be pointed out that orthogonality is not required to get good
        surface energies, but it can be useful in cases where the slabs are
        subsequently used for postprocessing of some kind, e.g. generating
        GBs or interfaces.
        """
        a, b, c = self.lattice.matrix
        new_c = np.cross(a, b)
        new_c /= np.linalg.norm(new_c)
        new_c = np.dot(c, new_c) * new_c
        new_latt = Lattice([a, b, new_c])
        return Slab(
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

    def get_tasker2_slabs(self, tol=0.01, same_species_only=True):
        """
        Get a list of slabs that have been Tasker 2 corrected.

        Args:
            tol (float): Tolerance to determine if atoms are within same plane.
                This is a fractional tolerance, not an absolute one.
            same_species_only (bool): If True, only that are of the exact same
                species as the atom at the outermost surface are considered for
                moving. Otherwise, all atoms regardless of species that is
                within tol are considered for moving. Default is True (usually
                the desired behavior).

        Returns:
            ([Slab]) List of tasker 2 corrected slabs.
        """
        sites = list(self.sites)
        slabs = []

        sortedcsites = sorted(sites, key=lambda site: site.c)

        # Determine what fraction the slab is of the total cell size in the
        # c direction. Round to nearest rational number.
        nlayers_total = int(round(self.lattice.c / self.oriented_unit_cell.lattice.c))
        nlayers_slab = int(round((sortedcsites[-1].c - sortedcsites[0].c) * nlayers_total))
        slab_ratio = nlayers_slab / nlayers_total

        a = SpacegroupAnalyzer(self)
        symm_structure = a.get_symmetrized_structure()

        def equi_index(site):
            for i, equi_sites in enumerate(symm_structure.equivalent_sites):
                if site in equi_sites:
                    return i
            raise ValueError("Cannot determine equi index!")

        for surface_site, shift in [
            (sortedcsites[0], slab_ratio),
            (sortedcsites[-1], -slab_ratio),
        ]:
            tomove = []
            fixed = []
            for site in sites:
                if abs(site.c - surface_site.c) < tol and (
                    (not same_species_only) or site.species == surface_site.species
                ):
                    tomove.append(site)
                else:
                    fixed.append(site)

            # Sort and group the sites by the species and symmetry equivalence
            tomove = sorted(tomove, key=lambda s: equi_index(s))

            grouped = [list(sites) for k, sites in itertools.groupby(tomove, key=lambda s: equi_index(s))]

            if len(tomove) == 0 or any(len(g) % 2 != 0 for g in grouped):
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
                fcoords = [site.frac_coords for site in fixed]

                for s in tomove:
                    species.append(s.species)
                    for group in selection:
                        if s in group:
                            fcoords.append(s.frac_coords)
                            break
                    else:
                        # Move unselected atom to the opposite surface.
                        fcoords.append(s.frac_coords + [0, 0, shift])

                # sort by species to put all similar species together.
                sp_fcoord = sorted(zip(species, fcoords), key=lambda x: x[0])
                species = [x[0] for x in sp_fcoord]
                fcoords = [x[1] for x in sp_fcoord]
                slab = Slab(
                    self.lattice,
                    species,
                    fcoords,
                    self.miller_index,
                    self.oriented_unit_cell,
                    self.shift,
                    self.scale_factor,
                    energy=self.energy,
                    reorient_lattice=self.reorient_lattice,
                )
                slabs.append(slab)
        s = StructureMatcher()
        unique = [ss[0] for ss in s.group_structures(slabs)]
        return unique

    def is_symmetric(self, symprec=0.1):
        """
        Checks if slab is symmetric, i.e., contains inversion symmetry.

        Args:
            symprec (float): Symmetry precision used for SpaceGroup analyzer.

        Returns:
            (bool) Whether slab contains inversion symmetry.
        """

        sg = SpacegroupAnalyzer(self, symprec=symprec)
        return sg.is_laue()

    def get_sorted_structure(self, key=None, reverse=False):
        """
        Get a sorted copy of the structure. The parameters have the same
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
        s = Structure.from_sites(sites)
        return Slab(
            s.lattice,
            s.species_and_occu,
            s.frac_coords,
            self.miller_index,
            self.oriented_unit_cell,
            self.shift,
            self.scale_factor,
            site_properties=s.site_properties,
            reorient_lattice=self.reorient_lattice,
        )

    def copy(self, site_properties=None, sanitize=False):
        """
        Convenience method to get a copy of the structure, with options to add
        site properties.

        Args:
            site_properties (dict): Properties to add or override. The
                properties are specified in the same way as the constructor,
                i.e., as a dict of the form {property: [values]}. The
                properties should be in the order of the *original* structure
                if you are performing sanitization.
            sanitize (bool): If True, this method will return a sanitized
                structure. Sanitization performs a few things: (i) The sites are
                sorted by electronegativity, (ii) a LLL lattice reduction is
                carried out to obtain a relatively orthogonalized cell,
                (iii) all fractional coords for sites are mapped into the
                unit cell.

        Returns:
            A copy of the Structure, with optionally new site_properties and
            optionally sanitized.
        """
        props = self.site_properties
        if site_properties:
            props.update(site_properties)
        return Slab(
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

    @property
    def dipole(self):
        """
        Calculates the dipole of the Slab in the direction of the surface
        normal. Note that the Slab must be oxidation state-decorated for this
        to work properly. Otherwise, the Slab will always have a dipole of 0.
        """
        dipole = np.zeros(3)
        mid_pt = np.sum(self.cart_coords, axis=0) / len(self)
        normal = self.normal
        for site in self:
            charge = sum([getattr(sp, "oxi_state", 0) * amt for sp, amt in site.species.items()])
            dipole += charge * np.dot(site.coords - mid_pt, normal) * normal
        return dipole

    def is_polar(self, tol_dipole_per_unit_area=1e-3):
        """
        Checks whether the surface is polar by computing the dipole per unit
        area. Note that the Slab must be oxidation state-decorated for this
        to work properly. Otherwise, the Slab will always be non-polar.

        Args:
            tol_dipole_per_unit_area (float): A tolerance. If the dipole
                magnitude per unit area is less than this value, the Slab is
                considered non-polar. Defaults to 1e-3, which is usually
                pretty good. Normalized dipole per unit area is used as it is
                more reliable than using the total, which tends to be larger for
                slabs with larger surface areas.
        """
        dip_per_unit_area = self.dipole / self.surface_area
        return np.linalg.norm(dip_per_unit_area) > tol_dipole_per_unit_area

    @property
    def normal(self):
        """
        Calculates the surface normal vector of the slab
        """
        normal = np.cross(self.lattice.matrix[0], self.lattice.matrix[1])
        normal /= np.linalg.norm(normal)
        return normal

    @property
    def surface_area(self):
        """
        Calculates the surface area of the slab
        """
        m = self.lattice.matrix
        return np.linalg.norm(np.cross(m[0], m[1]))

    @property
    def center_of_mass(self):
        """
        Calculates the center of mass of the slab
        """
        weights = [s.species.weight for s in self]
        center_of_mass = np.average(self.frac_coords, weights=weights, axis=0)
        return center_of_mass

    def add_adsorbate_atom(self, indices, specie, distance):
        """
        Gets the structure of single atom adsorption.
        slab structure from the Slab class(in [0, 0, 1])

        Args:
            indices ([int]): Indices of sites on which to put the absorbate.
                Absorbed atom will be displaced relative to the center of
                these sites.
            specie (Species/Element/str): adsorbed atom species
            distance (float): between centers of the adsorbed atom and the
                given site in Angstroms.
        """
        # Let's do the work in cartesian coords
        center = np.sum([self[i].coords for i in indices], axis=0) / len(indices)

        coords = center + self.normal * distance / np.linalg.norm(self.normal)

        self.append(specie, coords, coords_are_cartesian=True)

    def __str__(self):
        comp = self.composition
        outs = [
            "Slab Summary (%s)" % comp.formula,
            "Reduced Formula: %s" % comp.reduced_formula,
            "Miller index: %s" % (self.miller_index,),
            "Shift: %.4f, Scale Factor: %s" % (self.shift, self.scale_factor.__str__()),
        ]

        def to_s(x):
            return "%0.6f" % x

        outs.append("abc   : " + " ".join([to_s(i).rjust(10) for i in self.lattice.abc]))
        outs.append("angles: " + " ".join([to_s(i).rjust(10) for i in self.lattice.angles]))
        outs.append("Sites ({i})".format(i=len(self)))
        for i, site in enumerate(self):
            outs.append(
                " ".join(
                    [
                        str(i + 1),
                        site.species_string,
                        " ".join([to_s(j).rjust(12) for j in site.frac_coords]),
                    ]
                )
            )
        return "\n".join(outs)

    def as_dict(self):
        """
        :return: MSONAble dict
        """
        d = super().as_dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["oriented_unit_cell"] = self.oriented_unit_cell.as_dict()
        d["miller_index"] = self.miller_index
        d["shift"] = self.shift
        d["scale_factor"] = self.scale_factor.tolist()
        d["reconstruction"] = self.reconstruction
        d["energy"] = self.energy
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

        return Slab(
            lattice=lattice,
            species=s.species_and_occu,
            coords=s.frac_coords,
            miller_index=d["miller_index"],
            oriented_unit_cell=Structure.from_dict(d["oriented_unit_cell"]),
            shift=d["shift"],
            scale_factor=d["scale_factor"],
            site_properties=s.site_properties,
            energy=d["energy"],
        )

    def get_surface_sites(self, tag=False):
        """
        Returns the surface sites and their indices in a dictionary. The
        oriented unit cell of the slab will determine the coordination number
        of a typical site. We use VoronoiNN to determine the
        coordination number of bulk sites and slab sites. Due to the
        pathological error resulting from some surface sites in the
        VoronoiNN, we assume any site that has this error is a surface
        site as well. This will work for elemental systems only for now. Useful
        for analysis involving broken bonds and for finding adsorption sites.

            Args:
                tag (bool): Option to adds site attribute "is_surfsite" (bool)
                    to all sites of slab. Defaults to False

            Returns:
                A dictionary grouping sites on top and bottom of the slab
                together.
                {"top": [sites with indices], "bottom": [sites with indices}

        TODO:
            Is there a way to determine site equivalence between sites in a slab
            and bulk system? This would allow us get the coordination number of
            a specific site for multi-elemental systems or systems with more
            than one unequivalent site. This will allow us to use this for
            compound systems.
        """

        from pymatgen.analysis.local_env import VoronoiNN

        # Get a dictionary of coordination numbers
        # for each distinct site in the structure
        a = SpacegroupAnalyzer(self.oriented_unit_cell)
        ucell = a.get_symmetrized_structure()
        cn_dict = {}
        v = VoronoiNN()
        unique_indices = [equ[0] for equ in ucell.equivalent_indices]

        for i in unique_indices:
            el = ucell[i].species_string
            if el not in cn_dict.keys():
                cn_dict[el] = []
            # Since this will get the cn as a result of the weighted polyhedra, the
            # slightest difference in cn will indicate a different environment for a
            # species, eg. bond distance of each neighbor or neighbor species. The
            # decimal place to get some cn to be equal.
            cn = v.get_cn(ucell, i, use_weights=True)
            cn = float("%.5f" % (round(cn, 5)))
            if cn not in cn_dict[el]:
                cn_dict[el].append(cn)

        v = VoronoiNN()

        surf_sites_dict, properties = {"top": [], "bottom": []}, []
        for i, site in enumerate(self):
            # Determine if site is closer to the top or bottom of the slab
            top = site.frac_coords[2] > self.center_of_mass[2]

            try:
                # A site is a surface site, if its environment does
                # not fit the environment of other sites
                cn = float("%.5f" % (round(v.get_cn(self, i, use_weights=True), 5)))
                if cn < min(cn_dict[site.species_string]):
                    properties.append(True)
                    key = "top" if top else "bottom"
                    surf_sites_dict[key].append([site, i])
                else:
                    properties.append(False)
            except RuntimeError:
                # or if pathological error is returned, indicating a surface site
                properties.append(True)
                key = "top" if top else "bottom"
                surf_sites_dict[key].append([site, i])

        if tag:
            self.add_site_property("is_surf_site", properties)
        return surf_sites_dict

    def have_equivalent_surfaces(self):
        """
        Check if we have same number of equivalent sites on both surfaces.
        This is an alternative to checking Laue symmetry (is_symmetric())
        if we want to ensure both surfaces in the slab are the same
        """

        # tag the sites as either surface sites or not
        self.get_surface_sites(tag=True)

        a = SpacegroupAnalyzer(self)
        symm_structure = a.get_symmetrized_structure()

        # ensure each site on one surface has a
        # corresponding equivalent site on the other
        equal_surf_sites = []
        for equ in symm_structure.equivalent_sites:
            # Top and bottom are arbitrary, we will just determine
            # if one site is on one side of the slab or the other
            top, bottom = 0, 0
            for s in equ:
                if s.is_surf_site:
                    if s.frac_coords[2] > self.center_of_mass[2]:
                        top += 1
                    else:
                        bottom += 1
            # Check to see if the number of equivalent sites
            # on one side of the slab are equal to the other
            equal_surf_sites.append(top == bottom)

        return all(equal_surf_sites)

    def get_symmetric_site(self, point, cartesian=False):
        """
        This method uses symmetry operations to find equivalent sites on
            both sides of the slab. Works mainly for slabs with Laue
            symmetry. This is useful for retaining the non-polar and
            symmetric properties of a slab when creating adsorbed
            structures or symmetric reconstructions.

        Arg:
            point: Fractional coordinate.

        Returns:
            point: Fractional coordinate. A point equivalent to the
                parameter point, but on the other side of the slab
        """

        sg = SpacegroupAnalyzer(self)
        ops = sg.get_symmetry_operations(cartesian=cartesian)

        # Each operation on a point will return an equivalent point.
        # We want to find the point on the other side of the slab.
        for op in ops:
            slab = self.copy()
            site2 = op.operate(point)
            if "%.6f" % (site2[2]) == "%.6f" % (point[2]):
                continue

            # Add dummy site to check the overall structure is symmetric
            slab.append("O", point, coords_are_cartesian=cartesian)
            slab.append("O", site2, coords_are_cartesian=cartesian)
            sg = SpacegroupAnalyzer(slab)
            if sg.is_laue():
                break

            # If not symmetric, remove the two added
            # sites and try another symmetry operator
            slab.remove_sites([len(slab) - 1])
            slab.remove_sites([len(slab) - 1])

        return site2

    def symmetrically_add_atom(self, specie, point, coords_are_cartesian=False):
        """
        Class method for adding a site at a specified point in a slab.
            Will add the corresponding site on the other side of the
            slab to maintain equivalent surfaces.

        Arg:
            specie (str): The specie to add
            point (coords): The coordinate of the site in the slab to add.
            coords_are_cartesian (bool): Is the point in cartesian coordinates

        Returns:
            (Slab): The modified slab
        """

        # For now just use the species of the
        # surface atom as the element to add

        # Get the index of the corresponding site at the bottom
        point2 = self.get_symmetric_site(point, cartesian=coords_are_cartesian)

        self.append(specie, point, coords_are_cartesian=coords_are_cartesian)
        self.append(specie, point2, coords_are_cartesian=coords_are_cartesian)

    def symmetrically_remove_atoms(self, indices):
        """
        Class method for removing sites corresponding to a list of indices.
            Will remove the corresponding site on the other side of the
            slab to maintain equivalent surfaces.

        Arg:
            indices ([indices]): The indices of the sites
                in the slab to remove.
        """

        slabcopy = SpacegroupAnalyzer(self.copy()).get_symmetrized_structure()
        points = [slabcopy[i].frac_coords for i in indices]
        removal_list = []

        for pt in points:
            # Get the index of the original site on top
            cart_point = slabcopy.lattice.get_cartesian_coords(pt)
            dist = [site.distance_from_point(cart_point) for site in slabcopy]
            site1 = dist.index(min(dist))

            # Get the index of the corresponding site at the bottom
            for i, eq_sites in enumerate(slabcopy.equivalent_sites):
                if slabcopy[site1] in eq_sites:
                    eq_indices = slabcopy.equivalent_indices[i]
                    break
            i1 = eq_indices[eq_sites.index(slabcopy[site1])]

            for i2 in eq_indices:
                if i2 == i1:
                    continue
                if slabcopy[i2].frac_coords[2] == slabcopy[i1].frac_coords[2]:
                    continue
                # Test site remove to see if it results in symmetric slab
                s = self.copy()
                s.remove_sites([i1, i2])
                if s.is_symmetric():
                    removal_list.extend([i1, i2])
                    break

        # If expected, 2 atoms are removed per index
        if len(removal_list) == 2 * len(indices):
            self.remove_sites(removal_list)
        else:
            warnings.warn("Equivalent sites could not be found for removal for all indices. Surface unchanged.")


class SlabGenerator:
    """
    This class generates different slabs using shift values determined by where
    a unique termination can be found along with other criterias such as where a
    termination doesn't break a polyhedral bond. The shift value then indicates
    where the slab layer will begin and terminate in the slab-vacuum system.

    .. attribute:: oriented_unit_cell

        A unit cell of the parent structure with the miller
        index of plane parallel to surface

    .. attribute:: parent

        Parent structure from which Slab was derived.

    .. attribute:: lll_reduce

        Whether or not the slabs will be orthogonalized

    .. attribute:: center_slab

        Whether or not the slabs will be centered between
        the vacuum layer

    .. attribute:: slab_scale_factor

        Final computed scale factor that brings the parent cell to the
        surface cell.

    .. attribute:: miller_index

        Miller index of plane parallel to surface.

    .. attribute:: min_slab_size

        Minimum size in angstroms of layers containing atoms

    .. attribute:: min_vac_size

        Minimize size in angstroms of layers containing vacuum

    """

    def __init__(
        self,
        initial_structure,
        miller_index,
        min_slab_size,
        min_vacuum_size,
        lll_reduce=False,
        center_slab=False,
        in_unit_planes=False,
        primitive=True,
        max_normal_search=None,
        reorient_lattice=True,
    ):
        """
        Calculates the slab scale factor and uses it to generate a unit cell
        of the initial structure that has been oriented by its miller index.
        Also stores the initial information needed later on to generate a slab.

        Args:
            initial_structure (Structure): Initial input structure. Note that to
                ensure that the miller indices correspond to usual
                crystallographic definitions, you should supply a conventional
                unit cell structure.
            miller_index ([h, k, l]): Miller index of plane parallel to
                surface. Note that this is referenced to the input structure. If
                you need this to be based on the conventional cell,
                you should supply the conventional structure.
            min_slab_size (float): In Angstroms or number of hkl planes
            min_vacuum_size (float): In Angstroms or number of hkl planes
            lll_reduce (bool): Whether to perform an LLL reduction on the
                eventual structure.
            center_slab (bool): Whether to center the slab in the cell with
                equal vacuum spacing from the top and bottom.
            in_unit_planes (bool): Whether to set min_slab_size and min_vac_size
                in units of hkl planes (True) or Angstrom (False/default).
                Setting in units of planes is useful for ensuring some slabs
                have a certain nlayer of atoms. e.g. for Cs (100), a 10 Ang
                slab will result in a slab with only 2 layer of atoms, whereas
                Fe (100) will have more layer of atoms. By using units of hkl
                planes instead, we ensure both slabs
                have the same number of atoms. The slab thickness will be in
                min_slab_size/math.ceil(self._proj_height/dhkl)
                multiples of oriented unit cells.
            primitive (bool): Whether to reduce any generated slabs to a
                primitive cell (this does **not** mean the slab is generated
                from a primitive cell, it simply means that after slab
                generation, we attempt to find shorter lattice vectors,
                which lead to less surface area and smaller cells).
            max_normal_search (int): If set to a positive integer, the code will
                conduct a search for a normal lattice vector that is as
                perpendicular to the surface as possible by considering
                multiples linear combinations of lattice vectors up to
                max_normal_search. This has no bearing on surface energies,
                but may be useful as a preliminary step to generating slabs
                for absorption and other sizes. It is typical that this will
                not be the smallest possible cell for simulation. Normality
                is not guaranteed, but the oriented cell will have the c
                vector as normal as possible (within the search range) to the
                surface. A value of up to the max absolute Miller index is
                usually sufficient.
            reorient_lattice (bool): reorients the lattice parameters such that
                the c direction is the third vector of the lattice matrix

        """
        # pylint: disable=E1130
        # Add Wyckoff symbols of the bulk, will help with
        # identfying types of sites in the slab system
        sg = SpacegroupAnalyzer(initial_structure)
        initial_structure.add_site_property("bulk_wyckoff", sg.get_symmetry_dataset()["wyckoffs"])
        initial_structure.add_site_property("bulk_equivalent", sg.get_symmetry_dataset()["equivalent_atoms"].tolist())
        latt = initial_structure.lattice
        miller_index = _reduce_vector(miller_index)
        # Calculate the surface normal using the reciprocal lattice vector.
        recp = latt.reciprocal_lattice_crystallographic
        normal = recp.get_cartesian_coords(miller_index)
        normal /= np.linalg.norm(normal)

        slab_scale_factor = []
        non_orth_ind = []
        eye = np.eye(3, dtype=np.int_)
        for i, j in enumerate(miller_index):
            if j == 0:
                # Lattice vector is perpendicular to surface normal, i.e.,
                # in plane of surface. We will simply choose this lattice
                # vector as one of the basis vectors.
                slab_scale_factor.append(eye[i])
            else:
                # Calculate projection of lattice vector onto surface normal.
                d = abs(np.dot(normal, latt.matrix[i])) / latt.abc[i]
                non_orth_ind.append((i, d))

        # We want the vector that has maximum magnitude in the
        # direction of the surface normal as the c-direction.
        # Results in a more "orthogonal" unit cell.
        c_index, dist = max(non_orth_ind, key=lambda t: t[1])

        if len(non_orth_ind) > 1:
            lcm_miller = lcm(*[miller_index[i] for i, d in non_orth_ind])
            for (i, di), (j, dj) in itertools.combinations(non_orth_ind, 2):
                l = [0, 0, 0]
                l[i] = -int(round(lcm_miller / miller_index[i]))
                l[j] = int(round(lcm_miller / miller_index[j]))
                slab_scale_factor.append(l)
                if len(slab_scale_factor) == 2:
                    break

        if max_normal_search is None:
            slab_scale_factor.append(eye[c_index])
        else:

            index_range = sorted(
                reversed(range(-max_normal_search, max_normal_search + 1)),
                key=lambda x: abs(x),
            )
            candidates = []
            for uvw in itertools.product(index_range, index_range, index_range):
                if (not any(uvw)) or abs(np.linalg.det(slab_scale_factor + [uvw])) < 1e-8:
                    continue
                vec = latt.get_cartesian_coords(uvw)
                l = np.linalg.norm(vec)
                cosine = abs(np.dot(vec, normal) / l)
                candidates.append((uvw, cosine, l))
                if abs(abs(cosine) - 1) < 1e-8:
                    # If cosine of 1 is found, no need to search further.
                    break
            # We want the indices with the maximum absolute cosine,
            # but smallest possible length.
            uvw, cosine, l = max(candidates, key=lambda x: (x[1], -x[2]))
            slab_scale_factor.append(uvw)

        slab_scale_factor = np.array(slab_scale_factor)

        # Let's make sure we have a left-handed crystallographic system
        if np.linalg.det(slab_scale_factor) < 0:
            slab_scale_factor *= -1

        # Make sure the slab_scale_factor is reduced to avoid
        # unnecessarily large slabs

        reduced_scale_factor = [_reduce_vector(v) for v in slab_scale_factor]
        slab_scale_factor = np.array(reduced_scale_factor)

        single = initial_structure.copy()
        single.make_supercell(slab_scale_factor)

        # When getting the OUC, lets return the most reduced
        # structure as possible to reduce calculations
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
        self._normal = normal
        a, b, c = self.oriented_unit_cell.lattice.matrix
        self._proj_height = abs(np.dot(normal, c))
        self.reorient_lattice = reorient_lattice

    def get_slab(self, shift=0, tol=0.1, energy=None):
        """
        This method takes in shift value for the c lattice direction and
        generates a slab based on the given shift. You should rarely use this
        method. Instead, it is used by other generation algorithms to obtain
        all slabs.

        Arg:
            shift (float): A shift value in Angstrom that determines how much a
                slab should be shifted.
            tol (float): Tolerance to determine primitive cell.
            energy (float): An energy to assign to the slab.

        Returns:
            (Slab) A Slab object with a particular shifted oriented unit cell.
        """

        h = self._proj_height
        p = round(h / self.parent.lattice.d_hkl(self.miller_index), 8)
        if self.in_unit_planes:
            nlayers_slab = int(math.ceil(self.min_slab_size / p))
            nlayers_vac = int(math.ceil(self.min_vac_size / p))
        else:
            nlayers_slab = int(math.ceil(self.min_slab_size / h))
            nlayers_vac = int(math.ceil(self.min_vac_size / h))
        nlayers = nlayers_slab + nlayers_vac

        species = self.oriented_unit_cell.species_and_occu
        props = self.oriented_unit_cell.site_properties
        props = {k: v * nlayers_slab for k, v in props.items()}
        frac_coords = self.oriented_unit_cell.frac_coords
        frac_coords = np.array(frac_coords) + np.array([0, 0, -shift])[None, :]
        frac_coords -= np.floor(frac_coords)
        a, b, c = self.oriented_unit_cell.lattice.matrix
        new_lattice = [a, b, nlayers * c]
        frac_coords[:, 2] = frac_coords[:, 2] / nlayers
        all_coords = []
        for i in range(nlayers_slab):
            fcoords = frac_coords.copy()
            fcoords[:, 2] += i / nlayers
            all_coords.extend(fcoords)

        slab = Structure(new_lattice, species * nlayers_slab, all_coords, site_properties=props)

        scale_factor = self.slab_scale_factor
        # Whether or not to orthogonalize the structure
        if self.lll_reduce:
            lll_slab = slab.copy(sanitize=True)
            mapping = lll_slab.lattice.find_mapping(slab.lattice)
            scale_factor = np.dot(mapping[2], scale_factor)
            slab = lll_slab

        # Whether or not to center the slab layer around the vacuum
        if self.center_slab:
            avg_c = np.average([c[2] for c in slab.frac_coords])
            slab.translate_sites(list(range(len(slab))), [0, 0, 0.5 - avg_c])

        if self.primitive:
            prim = slab.get_primitive_structure(tolerance=tol)
            if energy is not None:
                energy = prim.volume / slab.volume * energy
            slab = prim

        # Reorient the lattice to get the correct reduced cell
        ouc = self.oriented_unit_cell.copy()
        if self.primitive:
            # find a reduced ouc
            slab_l = slab.lattice
            ouc = ouc.get_primitive_structure(
                constrain_latt={
                    "a": slab_l.a,
                    "b": slab_l.b,
                    "alpha": slab_l.alpha,
                    "beta": slab_l.beta,
                    "gamma": slab_l.gamma,
                }
            )
            # Check this is the correct oriented unit cell
            ouc = self.oriented_unit_cell if slab_l.a != ouc.lattice.a or slab_l.b != ouc.lattice.b else ouc

        return Slab(
            slab.lattice,
            slab.species_and_occu,
            slab.frac_coords,
            self.miller_index,
            ouc,
            shift,
            scale_factor,
            energy=energy,
            site_properties=slab.site_properties,
            reorient_lattice=self.reorient_lattice,
        )

    def _calculate_possible_shifts(self, tol=0.1):
        frac_coords = self.oriented_unit_cell.frac_coords
        n = len(frac_coords)

        if n == 1:
            # Clustering does not work when there is only one data point.
            shift = frac_coords[0][2] + 0.5
            return [shift - math.floor(shift)]

        # We cluster the sites according to the c coordinates. But we need to
        # take into account PBC. Let's compute a fractional c-coordinate
        # distance matrix that accounts for PBC.
        dist_matrix = np.zeros((n, n))
        h = self._proj_height
        # Projection of c lattice vector in
        # direction of surface normal.
        for i, j in itertools.combinations(list(range(n)), 2):
            if i != j:
                cdist = frac_coords[i][2] - frac_coords[j][2]
                cdist = abs(cdist - round(cdist)) * h
                dist_matrix[i, j] = cdist
                dist_matrix[j, i] = cdist

        condensed_m = squareform(dist_matrix)
        z = linkage(condensed_m)
        clusters = fcluster(z, tol, criterion="distance")

        # Generate dict of cluster# to c val - doesn't matter what the c is.
        c_loc = {c: frac_coords[i][2] for i, c in enumerate(clusters)}

        # Put all c into the unit cell.
        possible_c = [c - math.floor(c) for c in sorted(c_loc.values())]

        # Calculate the shifts
        nshifts = len(possible_c)
        shifts = []
        for i in range(nshifts):
            if i == nshifts - 1:
                # There is an additional shift between the first and last c
                # coordinate. But this needs special handling because of PBC.
                shift = (possible_c[0] + 1 + possible_c[i]) * 0.5
                if shift > 1:
                    shift -= 1
            else:
                shift = (possible_c[i] + possible_c[i + 1]) * 0.5
            shifts.append(shift - math.floor(shift))
        shifts = sorted(shifts)
        return shifts

    def _get_c_ranges(self, bonds):
        c_ranges = []
        bonds = {(get_el_sp(s1), get_el_sp(s2)): dist for (s1, s2), dist in bonds.items()}
        for (sp1, sp2), bond_dist in bonds.items():
            for site in self.oriented_unit_cell:
                if sp1 in site.species:
                    for nn in self.oriented_unit_cell.get_neighbors(site, bond_dist):
                        if sp2 in nn.species:
                            c_range = tuple(sorted([site.frac_coords[2], nn.frac_coords[2]]))
                            if c_range[1] > 1:
                                # Takes care of PBC when c coordinate of site
                                # goes beyond the upper boundary of the cell
                                c_ranges.append((c_range[0], 1))
                                c_ranges.append((0, c_range[1] - 1))
                            elif c_range[0] < 0:
                                # Takes care of PBC when c coordinate of site
                                # is below the lower boundary of the unit cell
                                c_ranges.append((0, c_range[1]))
                                c_ranges.append((c_range[0] + 1, 1))
                            elif c_range[0] != c_range[1]:
                                c_ranges.append((c_range[0], c_range[1]))
        return c_ranges

    def get_slabs(
        self,
        bonds=None,
        ftol=0.1,
        tol=0.1,
        max_broken_bonds=0,
        symmetrize=False,
        repair=False,
    ):
        """
        This method returns a list of slabs that are generated using the list of
        shift values from the method, _calculate_possible_shifts(). Before the
        shifts are used to create the slabs however, if the user decides to take
        into account whether or not a termination will break any polyhedral
        structure (bonds is not None), this method will filter out any shift
        values that do so.

        Args:
            bonds ({(specie1, specie2): max_bond_dist}: bonds are
                specified as a dict of tuples: float of specie1, specie2
                and the max bonding distance. For example, PO4 groups may be
                defined as {("P", "O"): 3}.
            tol (float): General tolerance paramter for getting primitive
                cells and matching structures
            ftol (float): Threshold parameter in fcluster in order to check
                if two atoms are lying on the same plane. Default thresh set
                to 0.1 Angstrom in the direction of the surface normal.
            max_broken_bonds (int): Maximum number of allowable broken bonds
                for the slab. Use this to limit # of slabs (some structures
                may have a lot of slabs). Defaults to zero, which means no
                defined bonds must be broken.
            symmetrize (bool): Whether or not to ensure the surfaces of the
                slabs are equivalent.
            repair (bool): Whether to repair terminations with broken bonds
                or just omit them. Set to False as repairing terminations can
                lead to many possible slabs as oppose to just omitting them.

        Returns:
            ([Slab]) List of all possible terminations of a particular surface.
            Slabs are sorted by the # of bonds broken.
        """
        c_ranges = [] if bonds is None else self._get_c_ranges(bonds)

        slabs = []
        for shift in self._calculate_possible_shifts(tol=ftol):
            bonds_broken = 0
            for r in c_ranges:
                if r[0] <= shift <= r[1]:
                    bonds_broken += 1
            slab = self.get_slab(shift, tol=tol, energy=bonds_broken)
            if bonds_broken <= max_broken_bonds:
                slabs.append(slab)
            elif repair:
                # If the number of broken bonds is exceeded,
                # we repair the broken bonds on the slab
                slabs.append(self.repair_broken_bonds(slab, bonds))

        # Further filters out any surfaces made that might be the same
        m = StructureMatcher(ltol=tol, stol=tol, primitive_cell=False, scale=False)

        new_slabs = []
        for g in m.group_structures(slabs):
            # For each unique termination, symmetrize the
            # surfaces by removing sites from the bottom.
            if symmetrize:
                slabs = self.nonstoichiometric_symmetrized_slab(g[0])
                new_slabs.extend(slabs)
            else:
                new_slabs.append(g[0])

        match = StructureMatcher(ltol=tol, stol=tol, primitive_cell=False, scale=False)
        new_slabs = [g[0] for g in match.group_structures(new_slabs)]

        return sorted(new_slabs, key=lambda s: s.energy)

    def repair_broken_bonds(self, slab, bonds):
        """
        This method will find undercoordinated atoms due to slab
        cleaving specified by the bonds parameter and move them
        to the other surface to make sure the bond is kept intact.
        In a future release of surface.py, the ghost_sites will be
        used to tell us how the repair bonds should look like.

        Arg:
            slab (structure): A structure object representing a slab.
            bonds ({(specie1, specie2): max_bond_dist}: bonds are
                specified as a dict of tuples: float of specie1, specie2
                and the max bonding distance. For example, PO4 groups may be
                defined as {("P", "O"): 3}.

        Returns:
            (Slab) A Slab object with a particular shifted oriented unit cell.
        """

        for pair in bonds.keys():
            blength = bonds[pair]

            # First lets determine which element should be the
            # reference (center element) to determine broken bonds.
            # e.g. P for a PO4 bond. Find integer coordination
            # numbers of the pair of elements wrt to each other
            cn_dict = {}
            for i, el in enumerate(pair):
                cnlist = []
                for site in self.oriented_unit_cell:
                    poly_coord = 0
                    if site.species_string == el:

                        for nn in self.oriented_unit_cell.get_neighbors(site, blength):
                            if nn[0].species_string == pair[i - 1]:
                                poly_coord += 1
                    cnlist.append(poly_coord)
                cn_dict[el] = cnlist

            # We make the element with the higher coordination our reference
            if max(cn_dict[pair[0]]) > max(cn_dict[pair[1]]):
                element1, element2 = pair
            else:
                element2, element1 = pair

            for i, site in enumerate(slab):
                # Determine the coordination of our reference
                if site.species_string == element1:
                    poly_coord = 0
                    for neighbor in slab.get_neighbors(site, blength):
                        poly_coord += 1 if neighbor.species_string == element2 else 0

                    # suppose we find an undercoordinated reference atom
                    if poly_coord not in cn_dict[element1]:
                        # We get the reference atom of the broken bonds
                        # (undercoordinated), move it to the other surface
                        slab = self.move_to_other_side(slab, [i])

                        # find its NNs with the corresponding
                        # species it should be coordinated with
                        neighbors = slab.get_neighbors(slab[i], blength, include_index=True)
                        tomove = [nn[2] for nn in neighbors if nn[0].species_string == element2]
                        tomove.append(i)
                        # and then move those NNs along with the central
                        # atom back to the other side of the slab again
                        slab = self.move_to_other_side(slab, tomove)

        return slab

    def move_to_other_side(self, init_slab, index_of_sites):
        """
        This method will Move a set of sites to the
        other side of the slab (opposite surface).

        Arg:
            init_slab (structure): A structure object representing a slab.
            index_of_sites (list of ints): The list of indices representing
                the sites we want to move to the other side.

        Returns:
            (Slab) A Slab object with a particular shifted oriented unit cell.
        """

        slab = init_slab.copy()

        # Determine what fraction the slab is of the total cell size
        # in the c direction. Round to nearest rational number.
        h = self._proj_height
        p = h / self.parent.lattice.d_hkl(self.miller_index)
        if self.in_unit_planes:
            nlayers_slab = int(math.ceil(self.min_slab_size / p))
            nlayers_vac = int(math.ceil(self.min_vac_size / p))
        else:
            nlayers_slab = int(math.ceil(self.min_slab_size / h))
            nlayers_vac = int(math.ceil(self.min_vac_size / h))
        nlayers = nlayers_slab + nlayers_vac
        slab_ratio = nlayers_slab / nlayers

        # Sort the index of sites based on which side they are on
        top_site_index = [i for i in index_of_sites if slab[i].frac_coords[2] > slab.center_of_mass[2]]
        bottom_site_index = [i for i in index_of_sites if slab[i].frac_coords[2] < slab.center_of_mass[2]]

        # Translate sites to the opposite surfaces
        slab.translate_sites(top_site_index, [0, 0, slab_ratio])
        slab.translate_sites(bottom_site_index, [0, 0, -slab_ratio])

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

    def nonstoichiometric_symmetrized_slab(self, init_slab, tol=1e-3):
        """
        This method checks whether or not the two surfaces of the slab are
        equivalent. If the point group of the slab has an inversion symmetry (
        ie. belong to one of the Laue groups), then it is assumed that the
        surfaces should be equivalent. Otherwise, sites at the bottom of the
        slab will be removed until the slab is symmetric. Note the removal of sites
        can destroy the stoichiometry of the slab. For non-elemental
        structures, the chemical potential will be needed to calculate surface energy.

        Arg:
            init_slab (Structure): A single slab structure
            tol (float): Tolerance for SpaceGroupanalyzer.

        Returns:
            Slab (structure): A symmetrized Slab object.
        """

        sg = SpacegroupAnalyzer(init_slab, symprec=tol)

        if sg.is_laue():
            return [init_slab]

        nonstoich_slabs = []
        # Build an equivalent surface slab for each of the different surfaces
        for top in [True, False]:
            asym = True
            slab = init_slab.copy()
            slab.energy = init_slab.energy

            while asym:
                # Keep removing sites from the bottom one by one until both
                # surfaces are symmetric or the number of sites removed has
                # exceeded 10 percent of the original slab

                c_dir = [site[2] for i, site in enumerate(slab.frac_coords)]

                if top:
                    slab.remove_sites([c_dir.index(max(c_dir))])
                else:
                    slab.remove_sites([c_dir.index(min(c_dir))])
                if len(slab) <= len(self.parent):
                    break

                # Check if the altered surface is symmetric
                sg = SpacegroupAnalyzer(slab, symprec=tol)
                if sg.is_laue():
                    asym = False
                    nonstoich_slabs.append(slab)

        if len(slab) <= len(self.parent):
            warnings.warn("Too many sites removed, please use a larger slab " "size.")

        return nonstoich_slabs


module_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(module_dir, "reconstructions_archive.json")) as data_file:
    reconstructions_archive = json.load(data_file)


class ReconstructionGenerator:
    """
    This class takes in a pre-defined dictionary specifying the parameters
    need to build a reconstructed slab such as the SlabGenerator parameters,
    transformation matrix, sites to remove/add and slab/vacuum size. It will
    then use the formatted instructions provided by the dictionary to build
    the desired reconstructed slab from the initial structure.

    .. attribute:: slabgen_params

        Parameters for the SlabGenerator

    .. trans_matrix::

        A 3x3 transformation matrix to generate the reconstructed
            slab. Only the a and b lattice vectors are actually
            changed while the c vector remains the same. This
            matrix is what the Wood's notation is based on.

    .. reconstruction_json::

        The full json or dictionary containing the instructions
            for building the reconstructed slab

    .. termination::

        The index of the termination of the slab

    TODO:
    - Right now there is no way to specify what atom is being
        added. In the future, use basis sets?
    """

    def __init__(self, initial_structure, min_slab_size, min_vacuum_size, reconstruction_name):
        """
        Generates reconstructed slabs from a set of instructions
            specified by a dictionary or json file.

        Args:
            initial_structure (Structure): Initial input structure. Note
                that to ensure that the miller indices correspond to usual
                crystallographic definitions, you should supply a conventional
                unit cell structure.
            min_slab_size (float): In Angstroms
            min_vacuum_size (float): In Angstroms

            reconstruction (str): Name of the dict containing the instructions
                for building a reconstructed slab. The dictionary can contain
                any item the creator deems relevant, however any instructions
                archived in pymatgen for public use needs to contain the
                following keys and items to ensure compatibility with the
                ReconstructionGenerator:

                    "name" (str): A descriptive name for the type of
                        reconstruction. Typically the name will have the type
                        of structure the reconstruction is for, the Miller
                        index, and Wood's notation along with anything to
                        describe the reconstruction: e.g.:
                        "fcc_110_missing_row_1x2"
                    "description" (str): A longer description of your
                        reconstruction. This is to help future contributors who
                        want to add other types of reconstructions to the
                        archive on pymatgen to check if the reconstruction
                        already exists. Please read the descriptions carefully
                        before adding a new type of reconstruction to ensure it
                        is not in the archive yet.
                    "reference" (str): Optional reference to where the
                        reconstruction was taken from or first observed.
                    "spacegroup" (dict): e.g. {"symbol": "Fm-3m", "number": 225}
                        Indicates what kind of structure is this reconstruction.
                    "miller_index" ([h,k,l]): Miller index of your reconstruction
                    "Woods_notation" (str): For a reconstruction, the a and b
                        lattice may change to accomodate the symmetry of the
                        reconstruction. This notation indicates the change in
                        the vectors relative to the primitive (p) or
                        conventional (c) slab cell. E.g. p(2x1):

                        Wood, E. A. (1964). Vocabulary of surface
                        crystallography. Journal of Applied Physics, 35(4),
                        1306–1312.

                    "transformation_matrix" (numpy array): A 3x3 matrix to
                        transform the slab. Only the a and b lattice vectors
                        should change while the c vector remains the same.
                    "SlabGenerator_parameters" (dict): A dictionary containing
                        the parameters for the SlabGenerator class excluding the
                        miller_index, min_slab_size and min_vac_size as the
                        Miller index is already specified and the min_slab_size
                        and min_vac_size can be changed regardless of what type
                        of reconstruction is used. Having a consistent set of
                        SlabGenerator parameters allows for the instructions to
                        be reused to consistently build a reconstructed slab.
                    "points_to_remove" (list of coords): A list of sites to
                        remove where the first two indices are fraction (in a
                        and b) and the third index is in units of 1/d (in c).
                    "points_to_add" (list of frac_coords): A list of sites to add
                        where the first two indices are fraction (in a an b) and
                        the third index is in units of 1/d (in c).

                    "base_reconstruction" (dict): Option to base a reconstruction on
                        an existing reconstruction model also exists to easily build
                        the instructions without repeating previous work. E.g. the
                        alpha reconstruction of halites is based on the octopolar
                        reconstruction but with the topmost atom removed. The dictionary
                        for the alpha reconstruction would therefore contain the item
                        "reconstruction_base": "halite_111_octopolar_2x2", and
                        additional sites for "points_to_remove" and "points_to_add"
                        can be added to modify this reconstruction.

                    For "points_to_remove" and "points_to_add", the third index for
                        the c vector is in units of 1/d where d is the spacing
                        between atoms along hkl (the c vector) and is relative to
                        the topmost site in the unreconstructed slab. e.g. a point
                        of [0.5, 0.25, 1] corresponds to the 0.5 frac_coord of a,
                        0.25 frac_coord of b and a distance of 1 atomic layer above
                        the topmost site. [0.5, 0.25, -0.5] where the third index
                        corresponds to a point half a atomic layer below the topmost
                        site. [0.5, 0.25, 0] corresponds to a point in the same
                        position along c as the topmost site. This is done because
                        while the primitive units of a and b will remain constant,
                        the user can vary the length of the c direction by changing
                        the slab layer or the vacuum layer.

            NOTE: THE DICTIONARY SHOULD ONLY CONTAIN "points_to_remove" AND
            "points_to_add" FOR THE TOP SURFACE. THE ReconstructionGenerator
            WILL MODIFY THE BOTTOM SURFACE ACCORDINGLY TO RETURN A SLAB WITH
            EQUIVALENT SURFACES.
        """

        if reconstruction_name not in reconstructions_archive.keys():
            raise KeyError(
                "The reconstruction_name entered (%s) does not exist in the "
                "archive. Please select from one of the following reconstructions: %s "
                "or add the appropriate dictionary to the archive file "
                "reconstructions_archive.json." % (reconstruction_name, list(reconstructions_archive.keys()))
            )

        # Get the instructions to build the reconstruction
        # from the reconstruction_archive
        recon_json = copy.deepcopy(reconstructions_archive[reconstruction_name])
        new_points_to_add, new_points_to_remove = [], []
        if "base_reconstruction" in recon_json.keys():
            if "points_to_add" in recon_json.keys():
                new_points_to_add = recon_json["points_to_add"]
            if "points_to_remove" in recon_json.keys():
                new_points_to_remove = recon_json["points_to_remove"]

            # Build new instructions from a base reconstruction
            recon_json = copy.deepcopy(reconstructions_archive[recon_json["base_reconstruction"]])
            if "points_to_add" in recon_json.keys():
                del recon_json["points_to_add"]
            if "points_to_remove" in recon_json.keys():
                del recon_json["points_to_remove"]
            if new_points_to_add:
                recon_json["points_to_add"] = new_points_to_add
            if new_points_to_remove:
                recon_json["points_to_remove"] = new_points_to_remove

        slabgen_params = copy.deepcopy(recon_json["SlabGenerator_parameters"])
        slabgen_params["initial_structure"] = initial_structure.copy()
        slabgen_params["miller_index"] = recon_json["miller_index"]
        slabgen_params["min_slab_size"] = min_slab_size
        slabgen_params["min_vacuum_size"] = min_vacuum_size

        self.slabgen_params = slabgen_params
        self.trans_matrix = recon_json["transformation_matrix"]
        self.reconstruction_json = recon_json
        self.name = reconstruction_name

    def build_slabs(self):
        """
        Builds the reconstructed slab by:
            (1) Obtaining the unreconstructed slab using the specified
                parameters for the SlabGenerator.
            (2) Applying the appropriate lattice transformation in the
                a and b lattice vectors.
            (3) Remove any specified sites from both surfaces.
            (4) Add any specified sites to both surfaces.

        Returns:
            (Slab): The reconstructed slab.
        """

        slabs = self.get_unreconstructed_slabs()
        recon_slabs = []

        for slab in slabs:
            d = get_d(slab)
            top_site = sorted(slab, key=lambda site: site.frac_coords[2])[-1].coords

            # Remove any specified sites
            if "points_to_remove" in self.reconstruction_json.keys():
                pts_to_rm = copy.deepcopy(self.reconstruction_json["points_to_remove"])
                for p in pts_to_rm:
                    p[2] = slab.lattice.get_fractional_coords([top_site[0], top_site[1], top_site[2] + p[2] * d])[2]
                    cart_point = slab.lattice.get_cartesian_coords(p)
                    dist = [site.distance_from_point(cart_point) for site in slab]
                    site1 = dist.index(min(dist))
                    slab.symmetrically_remove_atoms([site1])

            # Add any specified sites
            if "points_to_add" in self.reconstruction_json.keys():
                pts_to_add = copy.deepcopy(self.reconstruction_json["points_to_add"])
                for p in pts_to_add:
                    p[2] = slab.lattice.get_fractional_coords([top_site[0], top_site[1], top_site[2] + p[2] * d])[2]
                    slab.symmetrically_add_atom(slab[0].specie, p)

            slab.reconstruction = self.name
            setattr(slab, "recon_trans_matrix", self.trans_matrix)

            # Get the oriented_unit_cell with the same axb area.
            ouc = slab.oriented_unit_cell.copy()
            ouc.make_supercell(self.trans_matrix)
            slab.oriented_unit_cell = ouc
            recon_slabs.append(slab)

        return recon_slabs

    def get_unreconstructed_slabs(self):
        """
        Generates the unreconstructed or pristine super slab.
        """
        slabs = []
        for slab in SlabGenerator(**self.slabgen_params).get_slabs():
            slab.make_supercell(self.trans_matrix)
            slabs.append(slab)
        return slabs


def get_d(slab):
    """
    Determine the distance of space between
    each layer of atoms along c
    """
    sorted_sites = sorted(slab, key=lambda site: site.frac_coords[2])
    for i, site in enumerate(sorted_sites):
        if not "%.6f" % (site.frac_coords[2]) == "%.6f" % (sorted_sites[i + 1].frac_coords[2]):
            d = abs(site.frac_coords[2] - sorted_sites[i + 1].frac_coords[2])
            break
    return slab.lattice.get_cartesian_coords([0, 0, d])[2]


def is_already_analyzed(miller_index: tuple, miller_list: list, symm_ops: list) -> bool:
    """
    Helper function to check if a given Miller index is
    part of the family of indices of any index in a list

    Args:
        miller_index (tuple): The Miller index to analyze
        miller_list (list): List of Miller indices. If the given
            Miller index belongs in the same family as any of the
            indices in this list, return True, else return False
        symm_ops (list): Symmetry operations of a
            lattice, used to define family of indices
    """
    for op in symm_ops:
        if in_coord_list(miller_list, op.operate(miller_index)):
            return True
    return False


def get_symmetrically_equivalent_miller_indices(structure, miller_index, return_hkil=True):
    """
    Returns all symmetrically equivalent indices for a given structure. Analysis
    is based on the symmetry of the reciprocal lattice of the structure.

    Args:
        miller_index (tuple): Designates the family of Miller indices
            to find. Can be hkl or hkil for hexagonal systems
        return_hkil (bool): If true, return hkil form of Miller
            index for hexagonal systems, otherwise return hkl
    """

    # Change to hkl if hkil because in_coord_list only handles tuples of 3
    miller_index = (miller_index[0], miller_index[1], miller_index[3]) if len(miller_index) == 4 else miller_index
    mmi = max(np.abs(miller_index))
    r = list(range(-mmi, mmi + 1))
    r.reverse()

    sg = SpacegroupAnalyzer(structure)
    # Get distinct hkl planes from the rhombohedral setting if trigonal
    if sg.get_crystal_system() == "trigonal":
        prim_structure = SpacegroupAnalyzer(structure).get_primitive_standard_structure()
        symm_ops = prim_structure.lattice.get_recp_symmetry_operation()
    else:
        symm_ops = structure.lattice.get_recp_symmetry_operation()

    equivalent_millers = [miller_index]
    for miller in itertools.product(r, r, r):
        if miller == miller_index:
            continue
        if any(i != 0 for i in miller):
            if is_already_analyzed(miller, equivalent_millers, symm_ops):
                equivalent_millers.append(miller)

            # include larger Miller indices in the family of planes
            if all(mmi > i for i in np.abs(miller)) and not in_coord_list(equivalent_millers, miller):
                if is_already_analyzed(mmi * np.array(miller), equivalent_millers, symm_ops):
                    equivalent_millers.append(miller)

    if return_hkil and sg.get_crystal_system() in ["trigonal", "hexagonal"]:
        return [(hkl[0], hkl[1], -1 * hkl[0] - hkl[1], hkl[2]) for hkl in equivalent_millers]
    return equivalent_millers


def get_symmetrically_distinct_miller_indices(structure, max_index, return_hkil=False):
    """
    Returns all symmetrically distinct indices below a certain max-index for
    a given structure. Analysis is based on the symmetry of the reciprocal
    lattice of the structure.
    Args:
        structure (Structure): input structure.
        max_index (int): The maximum index. For example, a max_index of 1
            means that (100), (110), and (111) are returned for the cubic
            structure. All other indices are equivalent to one of these.
        return_hkil (bool): If true, return hkil form of Miller
            index for hexagonal systems, otherwise return hkl
    """

    r = list(range(-max_index, max_index + 1))
    r.reverse()

    # First we get a list of all hkls for conventional (including equivalent)
    conv_hkl_list = [miller for miller in itertools.product(r, r, r) if any(i != 0 for i in miller)]

    sg = SpacegroupAnalyzer(structure)
    # Get distinct hkl planes from the rhombohedral setting if trigonal
    if sg.get_crystal_system() == "trigonal":
        transf = sg.get_conventional_to_primitive_transformation_matrix()
        miller_list = [hkl_transformation(transf, hkl) for hkl in conv_hkl_list]
        prim_structure = SpacegroupAnalyzer(structure).get_primitive_standard_structure()
        symm_ops = prim_structure.lattice.get_recp_symmetry_operation()
    else:
        miller_list = conv_hkl_list
        symm_ops = structure.lattice.get_recp_symmetry_operation()

    unique_millers, unique_millers_conv = [], []

    for i, miller in enumerate(miller_list):
        d = abs(reduce(gcd, miller))
        miller = tuple(int(i / d) for i in miller)
        if not is_already_analyzed(miller, unique_millers, symm_ops):
            if sg.get_crystal_system() == "trigonal":
                # Now we find the distinct primitive hkls using
                # the primitive symmetry operations and their
                # corresponding hkls in the conventional setting
                unique_millers.append(miller)
                d = abs(reduce(gcd, conv_hkl_list[i]))
                cmiller = tuple(int(i / d) for i in conv_hkl_list[i])
                unique_millers_conv.append(cmiller)
            else:
                unique_millers.append(miller)
                unique_millers_conv.append(miller)

    if return_hkil and sg.get_crystal_system() in ["trigonal", "hexagonal"]:
        return [(hkl[0], hkl[1], -1 * hkl[0] - hkl[1], hkl[2]) for hkl in unique_millers_conv]
    return unique_millers_conv


def hkl_transformation(transf, miller_index):
    """
    Returns the Miller index from setting
    A to B using a transformation matrix
    Args:
        transf (3x3 array): The transformation matrix
            that transforms a lattice of A to B
        miller_index ([h, k, l]): Miller index to transform to setting B
    """
    # Get a matrix of whole numbers (ints)

    def lcm(a, b):
        return a * b // math.gcd(a, b)

    reduced_transf = reduce(lcm, [int(1 / i) for i in itertools.chain(*transf) if i != 0]) * transf
    reduced_transf = reduced_transf.astype(int)

    # perform the transformation
    t_hkl = np.dot(reduced_transf, miller_index)
    d = abs(reduce(gcd, t_hkl))
    t_hkl = np.array([int(i / d) for i in t_hkl])

    # get mostly positive oriented Miller index
    if len([i for i in t_hkl if i < 0]) > 1:
        t_hkl *= -1

    return tuple(t_hkl)


def generate_all_slabs(
    structure,
    max_index,
    min_slab_size,
    min_vacuum_size,
    bonds=None,
    tol=0.1,
    ftol=0.1,
    max_broken_bonds=0,
    lll_reduce=False,
    center_slab=False,
    primitive=True,
    max_normal_search=None,
    symmetrize=False,
    repair=False,
    include_reconstructions=False,
    in_unit_planes=False,
):
    """
    A function that finds all different slabs up to a certain miller index.
    Slabs oriented under certain Miller indices that are equivalent to other
    slabs in other Miller indices are filtered out using symmetry operations
    to get rid of any repetitive slabs. For example, under symmetry operations,
    CsCl has equivalent slabs in the (0,0,1), (0,1,0), and (1,0,0) direction.

    Args:
        structure (Structure): Initial input structure. Note that to
                ensure that the miller indices correspond to usual
                crystallographic definitions, you should supply a conventional
                unit cell structure.
        max_index (int): The maximum Miller index to go up to.
        min_slab_size (float): In Angstroms
        min_vacuum_size (float): In Angstroms
        bonds ({(specie1, specie2): max_bond_dist}: bonds are
            specified as a dict of tuples: float of specie1, specie2
            and the max bonding distance. For example, PO4 groups may be
            defined as {("P", "O"): 3}.
        tol (float): Threshold parameter in fcluster in order to check
            if two atoms are lying on the same plane. Default thresh set
            to 0.1 Angstrom in the direction of the surface normal.
        max_broken_bonds (int): Maximum number of allowable broken bonds
            for the slab. Use this to limit # of slabs (some structures
            may have a lot of slabs). Defaults to zero, which means no
            defined bonds must be broken.
        lll_reduce (bool): Whether to perform an LLL reduction on the
            eventual structure.
        center_slab (bool): Whether to center the slab in the cell with
            equal vacuum spacing from the top and bottom.
        primitive (bool): Whether to reduce any generated slabs to a
            primitive cell (this does **not** mean the slab is generated
            from a primitive cell, it simply means that after slab
            generation, we attempt to find shorter lattice vectors,
            which lead to less surface area and smaller cells).
        max_normal_search (int): If set to a positive integer, the code will
            conduct a search for a normal lattice vector that is as
            perpendicular to the surface as possible by considering
            multiples linear combinations of lattice vectors up to
            max_normal_search. This has no bearing on surface energies,
            but may be useful as a preliminary step to generating slabs
            for absorption and other sizes. It is typical that this will
            not be the smallest possible cell for simulation. Normality
            is not guaranteed, but the oriented cell will have the c
            vector as normal as possible (within the search range) to the
            surface. A value of up to the max absolute Miller index is
            usually sufficient.
        symmetrize (bool): Whether or not to ensure the surfaces of the
            slabs are equivalent.
        repair (bool): Whether to repair terminations with broken bonds
            or just omit them
        include_reconstructions (bool): Whether to include reconstructed
            slabs available in the reconstructions_archive.json file.
    """
    all_slabs = []

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
            logger.debug("%s has %d slabs... " % (miller, len(slabs)))
            all_slabs.extend(slabs)

    if include_reconstructions:
        sg = SpacegroupAnalyzer(structure)
        symbol = sg.get_space_group_symbol()
        # enumerate through all posisble reconstructions in the
        # archive available for this particular structure (spacegroup)
        for name, instructions in reconstructions_archive.items():
            if "base_reconstruction" in instructions.keys():
                instructions = reconstructions_archive[instructions["base_reconstruction"]]
            if instructions["spacegroup"]["symbol"] == symbol:
                # check if this reconstruction has a max index
                # equal or less than the given max index
                if max(instructions["miller_index"]) > max_index:
                    continue
                recon = ReconstructionGenerator(structure, min_slab_size, min_vacuum_size, name)
                all_slabs.extend(recon.build_slabs())

    return all_slabs


def get_slab_regions(slab, blength=3.5):
    """
    Function to get the ranges of the slab regions. Useful for discerning where
    the slab ends and vacuum begins if the slab is not fully within the cell
    Args:
        slab (Structure): Structure object modelling the surface
        blength (float, Ang): The bondlength between atoms. You generally
            want this value to be larger than the actual bondlengths in
            order to find atoms that are part of the slab
    """

    fcoords, indices, all_indices = [], [], []
    for site in slab:
        # find sites with c < 0 (noncontiguous)
        neighbors = slab.get_neighbors(site, blength, include_index=True, include_image=True)
        for nn in neighbors:
            if nn[0].frac_coords[2] < 0:
                # sites are noncontiguous within cell
                fcoords.append(nn[0].frac_coords[2])
                indices.append(nn[-2])
                if nn[-2] not in all_indices:
                    all_indices.append(nn[-2])

    if fcoords:
        # If slab is noncontiguous, locate the lowest
        # site within the upper region of the slab
        while fcoords:
            last_fcoords = copy.copy(fcoords)
            last_indices = copy.copy(indices)
            site = slab[indices[fcoords.index(min(fcoords))]]
            neighbors = slab.get_neighbors(site, blength, include_index=True, include_image=True)
            fcoords, indices = [], []
            for nn in neighbors:
                if 1 > nn[0].frac_coords[2] > 0 and nn[0].frac_coords[2] < site.frac_coords[2]:
                    # sites are noncontiguous within cell
                    fcoords.append(nn[0].frac_coords[2])
                    indices.append(nn[-2])
                    if nn[-2] not in all_indices:
                        all_indices.append(nn[-2])

        # Now locate the highest site within the lower region of the slab
        upper_fcoords = []
        for site in slab:
            if all(nn.index not in all_indices for nn in slab.get_neighbors(site, blength)):
                upper_fcoords.append(site.frac_coords[2])
        coords = copy.copy(last_fcoords) if not fcoords else copy.copy(fcoords)
        min_top = slab[last_indices[coords.index(min(coords))]].frac_coords[2]
        ranges = [[0, max(upper_fcoords)], [min_top, 1]]
    else:
        # If the entire slab region is within the slab cell, just
        # set the range as the highest and lowest site in the slab
        sorted_sites = sorted(slab, key=lambda site: site.frac_coords[2])
        ranges = [[sorted_sites[0].frac_coords[2], sorted_sites[-1].frac_coords[2]]]

    return ranges


def miller_index_from_sites(lattice, coords, coords_are_cartesian=True, round_dp=4, verbose=True):
    """
    Get the Miller index of a plane from a list of site coordinates.

    A minimum of 3 sets of coordinates are required. If more than 3 sets of
    coordinates are given, the best plane that minimises the distance to all
    points will be calculated.

    Args:
        lattice (list or Lattice): A 3x3 lattice matrix or `Lattice` object (for
            example obtained from Structure.lattice).
        coords (iterable): A list or numpy array of coordinates. Can be
            cartesian or fractional coordinates. If more than three sets of
            coordinates are provided, the best plane that minimises the
            distance to all sites will be calculated.
        coords_are_cartesian (bool, optional): Whether the coordinates are
            in cartesian space. If using fractional coordinates set to False.
        round_dp (int, optional): The number of decimal places to round the
            miller index to.
        verbose (bool, optional): Whether to print warnings.

    Returns:
        (tuple): The Miller index.
    """
    if not isinstance(lattice, Lattice):
        lattice = Lattice(lattice)

    return lattice.get_miller_index_from_coords(
        coords,
        coords_are_cartesian=coords_are_cartesian,
        round_dp=round_dp,
        verbose=verbose,
    )


def center_slab(slab):
    """
    The goal here is to ensure the center of the slab region
        is centered close to c=0.5. This makes it easier to
        find the surface sites and apply operations like doping.

    There are three cases where the slab in not centered:

    1. The slab region is completely between two vacuums in the
    box but not necessarily centered. We simply shift the
    slab by the difference in its center of mass and 0.5
    along the c direction.

    2. The slab completely spills outside the box from the bottom
    and into the top. This makes it incredibly difficult to
    locate surface sites. We iterate through all sites that
    spill over (z>c) and shift all sites such that this specific
    site is now on the other side. Repeat for all sites with z>c.

    3. This is a simpler case of scenario 2. Either the top or bottom
    slab sites are at c=0 or c=1. Treat as scenario 2.

    Args:
        slab (Slab): Slab structure to center
    Returns:
        Returns a centered slab structure
    """

    # get a reasonable r cutoff to sample neighbors
    bdists = sorted([nn[1] for nn in slab.get_neighbors(slab[0], 10) if nn[1] > 0])
    r = bdists[0] * 3

    all_indices = [i for i, site in enumerate(slab)]

    # check if structure is case 2 or 3, shift all the
    # sites up to the other side until it is case 1
    for site in slab:
        if any(nn[1] > slab.lattice.c for nn in slab.get_neighbors(site, r)):
            shift = 1 - site.frac_coords[2] + 0.05
            slab.translate_sites(all_indices, [0, 0, shift])

    # now the slab is case 1, shift the center of mass of the slab to 0.5
    weights = [s.species.weight for s in slab]
    center_of_mass = np.average(slab.frac_coords, weights=weights, axis=0)
    shift = 0.5 - center_of_mass[2]
    slab.translate_sites(all_indices, [0, 0, shift])

    return slab


def _reduce_vector(vector):
    # small function to reduce vectors

    d = abs(reduce(gcd, vector))
    vector = tuple(int(i / d) for i in vector)

    return vector
