# coding: utf-8

# should go into pymatgen.core.structure module

from __future__ import division, print_function, unicode_literals

from math import fabs, pi, cos, sin
import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import Element
from pymatgen.transformations.standard_transformations import \
        SupercellTransformation

def get_borders_of_enveloping_cuboid(lattice):
    """
    Determines the borders of the cuboid (rectangular parallelepiped)
    that envelops the unitcell of the lattice provided.
    Could go into Lattice class itself?!

    Args:
        lattice (Lattice): lattice object for which the enveloping
            cuboid is to be determined.

    Returns:
        borders (2 x 3 array): borders of the cuboid that envelops the
            unitcell of the lattice provided; for example, borders[0][0]
            stores the left (smallest value) border in x direction,
            borders[1][0] the right (largest value) border in x,
            and borders[1][2] the right border in z direction.
    """
    corners = []
    corners.append(lattice.get_cartesian_coords([0.0, 0.0, 0.0]))
    corners.append(lattice.get_cartesian_coords([0.0, 1.0, 0.0]))
    corners.append(lattice.get_cartesian_coords([1.0, 0.0, 0.0]))
    corners.append(lattice.get_cartesian_coords([1.0, 1.0, 0.0]))
    corners.append(lattice.get_cartesian_coords([0.0, 0.0, 1.0]))
    corners.append(lattice.get_cartesian_coords([0.0, 1.0, 1.0]))
    corners.append(lattice.get_cartesian_coords([1.0, 0.0, 1.0]))
    corners.append(lattice.get_cartesian_coords([1.0, 1.0, 1.0]))
    borders = np.zeros((2, 3), dtype=float)
    for d in range(3):
        for i, corner in enumerate(corners):
            if i == 0 or corner[d] < borders[0][d]:
                borders[0][d] = corner[d]
            if i == 0 or corner[d] > borders[1][d]:
                borders[1][d] = corner[d]
    return borders

class PrototypeStructure(Structure):
    """
    Enable facile generation of Structure objects of common prototypes
    [e.g., primitive cubic (pc), face-centered cubic (fcc), hexagonal
    closed packed (hcp)].  Useful together with typical site
    perturbations and structure rotations to test local order parameters
    in both ways: conceptual performance and concrete implementation
    in pymatgen.
    """

    __types__ = {"c":   ["cubic", "primitive cubic", "pc", "c"],\
                 "bcc": ["body-centered cubic", "body centered cubic", "bcc"],\
                 "fcc": ["face-centered cubic", "face centered cubic", "fcc"],\
                 "hcp": ["hexagonal closed packed", "hcp"],\
                 "d":   ["diamond", "d"],\
                }

    def __init__(self, prototype, spec="H", lengths=[], angles=[], \
                 supercell_scaling=[1, 1, 1], enveloping_box_factor=0.0):
        """
        Create a PrototypeStructure object.

        Args:
            prototype (string): target prototype;  following
                self-explanatory values are recognized:
                    "cubic" or "primitive cubic" or "pc" or "c";
                    "body-centered cubic" or "body centered cubic" or
                        "bcc";
                    "face-centered cubic" or "face centered cubic" or
                        "fcc";
                    "hexagonal closed packed" or "hcp";
                    "diamond" or "d".
            spec (string): specie symbol for all sites in prototype
                structure.
            lengths ([3 x float]): lengths of unit-cell vectors a, b,
                and c.  Default is 1.0 for each axis.  Note that all
                prototype structures are setup with fractional coordinates.
                This permits the setup of a "body-centered tetragonal"
                structure, for example, by specifying
                lengths=[1.0, 1.0, 2.0].
            angles ([3 x float]): lattice angles.  If specified, the
                typical lattice angles of the chosen prototype structure
                will be overwritten.
            supercell_scaling ([3 x int]): number of
                times the unit cell is to be repeated along each of the
                three crystallographic axes.
            enveloping_box_factor (float): determines the size of the
                enveloping cube around the supercell. 
                Note that the factor is only considered if it is greater
                than unity.  In such case, the largest extension of the
                supercell along any one of the three Cartesian coordinates
                is multiplied by enveloping_box_factor.  The resulting value
                is used to define the final cubic lattice of the prototype
                structure.  The sites from the supercell structure are then
                placed in the center of the cube.  If the factor is smaller or
                equal unity, no enveloping box is placed around the
                supercell structure, thus, preserving the natural lattice
                size and shape from supercell scaling.  A reasonable value
                seems to be 2.0 because that allows any rotation of the
                structure around the cube center without sites crossing
                the boundaries of the enveloping cube.
        """
        self.prototype = ""
        for type_short, type_name_list in self.__types__.items():
            if prototype in type_name_list:
                self.prototype = type_short
                break
        if not self.prototype:
            raise ValueError("could not find input prototype \"" \
                             +prototype+"\"")

        # Generate unit-cell structure.
        coords = []
        specs  = []
        local_lengths = [1.0, 1.0, 1.0]
        local_angles = [90.0, 90.0, 90.0]
        if self.prototype == "c":
            coords = [[0.0, 0.0, 0.0]]
            specs  = [spec]
        elif self.prototype == "bcc":
            coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
            specs  = [spec, spec]
        elif self.prototype == "fcc":
            coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
                      [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
            specs  = [spec, spec, spec, spec]
        elif self.prototype == "hcp":
            coords = [[0.3333, 0.6667, 0.25], [0.6667, 0.3333, 0.75]]
            specs  = [spec, spec]
            local_lengths = [1.0, 1.0, 1.633]
            local_angles = [90.0, 90.0, 120.0]
        elif self.prototype == "d":
            coords = [[0.0, 0.0, 0.5], [0.75, 0.75, 0.75],
                      [0.0, 0.5, 0.0], [0.75, 0.25, 0.25],
                      [0.5, 0.0, 0.0], [0.25, 0.75, 0.25],
                      [0.5, 0.5, 0.5], [0.25, 0.25, 0.75]]
            specs  = [spec, spec, spec, spec, spec, spec, spec, spec]
        else:
            raise RuntimeError("unrecognized prototype!")
        if not coords or not specs:
            raise RuntimeError("could not setup unitcell for unknown reason!")
        if lengths:
            local_lengths = lengths
        if angles:
            local_angles = angles
        self.unitcell = Structure(
                Lattice.from_lengths_and_angles(
                        local_lengths, local_angles),
                specs, coords, validate_proximity=False, to_unit_cell=False,
                coords_are_cartesian=False, site_properties=None)

        # Generate supercell structure.
        if supercell_scaling[0] < 1 or supercell_scaling[1] < 1 or \
                supercell_scaling[2] < 1:
            raise ValueError("supercell scaling entries must all be"
                             " larger than zero!")
        self.supercell = SupercellTransformation(
                scaling_matrix=((supercell_scaling[0], 0, 0),
                        (0, supercell_scaling[1], 0),
                        (0, 0, supercell_scaling[2]))).apply_transformation(
                                self.unitcell)

        # If desired, place a bounding cube around the supercell structure.
        if enveloping_box_factor > 1.0:
            # Determine smallest orthogonal box that envelops the unit cell.
            borders = get_borders_of_enveloping_cuboid(self.supercell.lattice)
            enveloping_box_lengths = borders[1] - borders[0]
            max_enveloping_box_length = max(enveloping_box_lengths)
            shift = - borders[0] - 0.5*enveloping_box_lengths + 0.5*np.array(
                    [enveloping_box_factor*max_enveloping_box_length,
                     enveloping_box_factor*max_enveloping_box_length,
                     enveloping_box_factor*max_enveloping_box_length])

            # Use largest length of enveloping box,
            # max_enveloping_box_length,
            # in order to define final enveloping cubic box;
            # recommended: 2 * max_enveloping_box_length.
            tmp_coords = self.supercell.sites[0].coords + shift
            structure = Structure(
                    Lattice.cubic(
                            enveloping_box_factor*max_enveloping_box_length),
                    [self.supercell.sites[0].specie], [tmp_coords],
                    validate_proximity=False, to_unit_cell=False,
                    coords_are_cartesian=True, site_properties=None)
            for i in range(1, len(self.supercell.sites)):
                tmp_coords = self.supercell.sites[i].coords + shift
                structure.append(
                        self.supercell.sites[i].specie, tmp_coords,
                        coords_are_cartesian=True)

        # Because bounding cube is not wanted, we simply copy the supercell
        # structure.
        else:
            structure = self.supercell.copy()

        specs = [site.specie for site in structure.sites]
        coords = [site.coords for site in structure.sites]
        super(PrototypeStructure, self).__init__(
                structure.lattice, specs, coords,
                validate_proximity=False, to_unit_cell=False,
                coords_are_cartesian=True, site_properties=None)


    def get_index_of_site_closest_to_center_of_enveloping_cuboid(self):
        """
        Determine (one of) the site(s) nearest to the center of the
        enveloping cuboid.  Note that periodic boundary conditions
        are not applied.

        Return:
            index (int): index of site that is nearest to the center of
                the enveloping cuboid.  In the case that there are
                multiple sites equally far away from the center, the first
                site obtained from the standard iterator is used.
        """
        borders = get_borders_of_enveloping_cuboid(self._lattice)
        center = borders[0] + 0.5 * (borders[1] - borders[0])

        closest_site = None
        index = -1
        nearest2 = 0.0
        for i, site in enumerate(self._sites):
            distance_vector = site.coords - center
            distance2 = np.dot(distance_vector, distance_vector)
            if i == 0 or distance2 < nearest2:
                index = i
                nearest2 = distance2
                closest_site = site
        if not closest_site:
            raise RuntimeError("could not find closest site to enveloping" \
                               " cuboid center!")

        return index

    def rotate_sites_around_center_of_enveloping_cuboid(self, axis, angle):
        """
        Rotate the sites of the structure by angle degrees around the
        vector axis with the origin of rotation located at the center
        of the unit-cell enveloping cuboid.  The rotation matrix is taken
        from the "Rotation matrix from axis and angle" subsection on
        https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations.

        Args:
            axis ([3x1 float array]): axis about which the structure is to
                be rotated.
            angle (float): angle in degrees by which the structure is to
                be rotated around the axis vector.

        Return:
            rot_struct (Structure): rotated structure.
        """
        # Copy this structure for returning, determine center of enveloping
        # cuboid around which rotation is to be performed, convert angle,
        # and normalize rotation axis.
        borders = get_borders_of_enveloping_cuboid(self._lattice)
        center = borders[0] + 0.5 * (borders[1] - borders[0])
        theta = angle * pi / 180.0
        cos_theta = cos(theta)
        sin_theta = sin(theta)
        if len(axis) != 3:
            raise TypeError("expected a 3D vector for the rotation axis!")
        axis_length = np.linalg.norm(np.array(axis))
        if fabs(axis_length) <= 1.0e-8:
            raise ValueError("zero-length rotation axis encountered!")
        u   = np.array(axis) / axis_length
        ux  = u[0]
        ux2 = ux * ux
        uy  = u[1]
        uy2 = uy * uy
        uz  = u[2]
        uz2 = uz * uz
        R = np.array([[cos_theta + ux2 * (1.0 - cos_theta),
                       ux * uy * (1.0 - cos_theta) - uz * sin_theta,
                       ux * uz * (1.0 - cos_theta) + uy * sin_theta],
                      [uy * ux * (1.0 - cos_theta) + uz * sin_theta,
                       cos_theta + uy2 * (1.0 - cos_theta),
                       uy * uz * (1.0 - cos_theta) - ux * sin_theta],
                      [uz * ux * (1.0 - cos_theta) - uy * sin_theta,
                       uz * uy * (1.0 - cos_theta) + ux * sin_theta,
                       cos_theta + uz2 * (1.0 - cos_theta)]])

        # Shift each site so that center of enveloping cuboid is located
        # at the origin, rotate site around axis, and, finally,
        # shift it back and store it.
        rot_struct_sites = []
        x = np.array((3), float)
        for site in self._sites:
            x = site.coords - center
            x = np.dot(R, x)
            x = x + center
            rot_struct_sites.append(
                    PeriodicSite(site.specie, x, self._lattice,
                            to_unit_cell=False, coords_are_cartesian=True,
                            properties=None))
        rot_struct = Structure.from_sites(rot_struct_sites)
        return rot_struct

    def perturb_einstein_crystal_style(self, sqrt_kBT_over_kspring=1.0):
        """
        Perturb the position of each site so that the distribution around
        the equilibrium position yields a normal distribution for each
        Cartesian component.  The perturbation complies thus with the
        expectation for an Einstein crystal, in which the potential is
        given by V(dr) = 1/2 * kspring * (dr)^2.  kspring denotes
        the spring constant with which the sites are tethered to their
        equilibrium position, and dr is the distance of the site under
        consideration from its equilibrium position.  The displacements
        are obtained with numpy's random.randn() method, the values of
        which are scaled by sqrt_kBT_over_kspring.
        TODO: move to Structure class!

        Args:
            sqrt_kBT_over_kspring (float): width of the underlying normal
                distribution, which is sigma = (kB*T/kspring)^0.5.
        """
        for i in range(len(self._sites)):
            displ = np.random.randn(3)
            displ = sqrt_kBT_over_kspring * displ
            self.translate_sites([i], displ, frac_coords=False)
