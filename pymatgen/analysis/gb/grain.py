# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import numpy as np
from fractions import Fraction
from math import gcd, floor, cos
from functools import reduce
from pymatgen import Structure, Lattice
from pymatgen.core.sites import PeriodicSite
from monty.fractions import lcm
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import itertools

import logging
import warnings

# This module implements representations of grain boundaries, as well as
# algorithms for generating them.

__author__ = "Xiang-Guo Li"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Xiang-Guo Li"
__email__ = "xil110@ucsd.edu"
__date__ = "7/30/18"

logger = logging.getLogger(__name__)


class GrainBoundary(Structure):
    """
    Subclass of Structure representing a GrainBoundary (gb) object.
    Implements additional attributes pertaining to gbs, but the
    init method does not actually implement any algorithm that
    creates a gb. This is a DUMMY class who's init method only holds
    information about the gb. Also has additional methods that returns
    other information about a gb such as sigma value.

    Note that all gbs have the gb surface normal oriented in the c-direction.
    This means the lattice vectors a and b are in the gb surface plane (at
     least for one grain) and the c vector is out of the surface plane
     (though not necessary perpendicular to the surface.)
    """

    def __init__(self, lattice, species, coords, rotation_axis, rotation_angle,
                 gb_plane, join_plane, init_cell, vacuum_thickness, ab_shift,
                 site_properties, oriented_unit_cell, validate_proximity=False,
                 coords_are_cartesian=False):
        """
        Makes a gb structure, a structure object with additional information
        and methods pertaining to gbs.

        Args:
            lattice (Lattice/3x3 array): The lattice, either as a
                :class:`pymatgen.core.lattice.Lattice` or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species ([Specie]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of
                each species.
            rotation_axis (list): Rotation axis of GB in the form of a list of
                integers, e.g. [1, 1, 0].
            rotation_angle (float, in unit of degree): rotation angle of GB.
            gb_plane (list): Grain boundary plane in the form of a list of integers
                e.g.: [1, 2, 3].
            join_plane (list): Joining plane of the second grain in the form of a list of
                integers. e.g.: [1, 2, 3].
            init_cell (Structure): initial bulk structure to form the GB.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, The sequences have to be the same length as
                the atomic species and fractional_coords. For gb, you should
                have the 'grain_label' properties to classify the sites as 'top',
                'bottom', 'top_incident', or 'bottom_incident'.
            vacuum_thickness (float in angstrom): The thickness of vacuum inserted
                between two grains of the GB.
            ab_shift (list of float, in unit of crystal vector a, b): The relative
                shift along a, b vectors.
            oriented_unit_cell (Structure): oriented unit cell of the bulk init_cell.
                Help to accurate calculate the bulk properties that are consistent
                with gb calculations.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in cartesian coordinates. Defaults to False.
        """
        self.oriented_unit_cell = oriented_unit_cell
        self.rotation_axis = rotation_axis
        self.rotation_angle = rotation_angle
        self.gb_plane = gb_plane
        self.join_plane = join_plane
        self.init_cell = init_cell
        self.vacuum_thickness = vacuum_thickness
        self.ab_shift = ab_shift
        super().__init__(
            lattice, species, coords, validate_proximity=validate_proximity,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties)

    def copy(self):
        """
        Convenience method to get a copy of the structure, with options to add
        site properties.

        Returns:
            A copy of the Structure, with optionally new site_properties and
            optionally sanitized.
        """
        return GrainBoundary(self.lattice, self.species_and_occu, self.frac_coords,
                             self.rotation_axis, self.rotation_angle, self.gb_plane,
                             self.join_plane, self.init_cell, self.vacuum_thickness,
                             self.ab_shift, self.site_properties, self.oriented_unit_cell)

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
        return GrainBoundary(s.lattice, s.species_and_occu, s.frac_coords,
                             self.rotation_axis, self.rotation_angle, self.gb_plane,
                             self.join_plane, self.init_cell, self.vacuum_thickness,
                             self.ab_shift, self.site_properties, self.oriented_unit_cell)

    @property
    def sigma(self):
        """
        This method returns the sigma value of the gb.
        If using 'quick_gen' to generate GB, this value is not valid.
        """
        return int(round(self.oriented_unit_cell.volume / self.init_cell.volume))

    @property
    def sigma_from_site_prop(self):
        """
        This method returns the sigma value of the gb from site properties.
        If the GB structure merge some atoms due to the atoms too closer with
        each other, this property will not work.
        """
        num_coi = 0
        if None in self.site_properties['grain_label']:
            raise RuntimeError('Site were merged, this property do not work')
        for tag in self.site_properties['grain_label']:
            if 'incident' in tag:
                num_coi += 1
        return int(round(self.num_sites / num_coi))

    @property
    def top_grain(self):
        """
        return the top grain (Structure) of the GB.
        """
        top_sites = []
        for i, tag in enumerate(self.site_properties['grain_label']):
            if 'top' in tag:
                top_sites.append(self.sites[i])
        return Structure.from_sites(top_sites)

    @property
    def bottom_grain(self):
        """
        return the bottom grain (Structure) of the GB.
        """
        bottom_sites = []
        for i, tag in enumerate(self.site_properties['grain_label']):
            if 'bottom' in tag:
                bottom_sites.append(self.sites[i])
        return Structure.from_sites(bottom_sites)

    @property
    def coincidents(self):
        """
        return the a list of coincident sites.
        """
        coincident_sites = []
        for i, tag in enumerate(self.site_properties['grain_label']):
            if 'incident' in tag:
                coincident_sites.append(self.sites[i])
        return coincident_sites

    def __str__(self):
        comp = self.composition
        outs = [
            "Gb Summary (%s)" % comp.formula,
            "Reduced Formula: %s" % comp.reduced_formula,
            "Rotation axis: %s" % (self.rotation_axis,),
            "Rotation angle: %s" % (self.rotation_angle,),
            "GB plane: %s" % (self.gb_plane,),
            "Join plane: %s" % (self.join_plane,),
            "vacuum thickness: %s" % (self.vacuum_thickness,),
            "ab_shift: %s" % (self.ab_shift,), ]
        to_s = lambda x: "%0.6f" % x
        outs.append("abc   : " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.abc]))
        outs.append("angles: " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.angles]))
        outs.append("Sites ({i})".format(i=len(self)))
        for i, site in enumerate(self):
            outs.append(" ".join([str(i + 1), site.species_string,
                                  " ".join([to_s(j).rjust(12)
                                            for j in site.frac_coords])]))
        return "\n".join(outs)

    def as_dict(self):
        d = super().as_dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["init_cell"] = self.init_cell.as_dict()
        d["rotation_axis"] = self.rotation_axis
        d["rotation_angle"] = self.rotation_angle
        d["gb_plane"] = self.gb_plane
        d["join_plane"] = self.join_plane
        d["vacuum_thickness"] = self.vacuum_thickness
        d["ab_shift"] = self.ab_shift
        d["oriented_unit_cell"] = self.oriented_unit_cell.as_dict()
        return d

    @classmethod
    def from_dict(cls, d):
        lattice = Lattice.from_dict(d["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in d["sites"]]
        s = Structure.from_sites(sites)

        return GrainBoundary(
            lattice=lattice,
            species=s.species_and_occu, coords=s.frac_coords,
            rotation_axis=d["rotation_axis"],
            rotation_angle=d["rotation_angle"],
            gb_plane=d["gb_plane"],
            join_plane=d["join_plane"],
            init_cell=Structure.from_dict(d["init_cell"]),
            vacuum_thickness=d["vacuum_thickness"],
            ab_shift=d["ab_shift"],
            oriented_unit_cell=Structure.from_dict(d["oriented_unit_cell"]),
            site_properties=s.site_properties)


class GrainBoundaryGenerator:
    """
    This class is to generate grain boundaries (GBs) from bulk
    conventional cell (fcc, bcc can from the primitive cell), and works for Cubic,
    Tetragonal, Orthorhombic, Rhombohedral, and Hexagonal systems.
    It generate GBs from given parameters, which includes
    GB plane, rotation axis, rotation angle.

    This class works for any general GB, including twist, tilt and mixed GBs.
    The three parameters, rotation axis, GB plane and rotation angle, are
    sufficient to identify one unique GB. While sometimes, users may not be able
    to tell what exactly rotation angle is but prefer to use sigma as an parameter,
    this class also provides the function that is able to return all possible
    rotation angles for a specific sigma value.
    The same sigma value (with rotation axis fixed) can correspond to
    multiple rotation angles.
    Users can use structure matcher in pymatgen to get rid of the redundant structures.
    """

    def __init__(self, initial_structure, symprec=0.1, angle_tolerance=1):

        """
        initial_structure (Structure): Initial input structure. It can
               be conventional or primitive cell (primitive cell works for bcc and fcc).
               For fcc and bcc, using conventional cell can lead to a non-primitive
               grain boundary structure.
               This code supplies Cubic, Tetragonal, Orthorhombic, Rhombohedral, and
               Hexagonal systems.
        symprec (float): Tolerance for symmetry finding. Defaults to 0.1 (the value used
                in Materials Project), which is for structures with slight deviations
                from their proper atomic positions (e.g., structures relaxed with
                electronic structure codes).
                A smaller value of 0.01 is often used for properly refined
                structures with atoms in the proper symmetry coordinates.
                User should make sure the symmetry is what you want.
        angle_tolerance (float): Angle tolerance for symmetry finding.
        """
        analyzer = SpacegroupAnalyzer(initial_structure, symprec, angle_tolerance)
        self.lat_type = analyzer.get_lattice_type()[0]
        if (self.lat_type == 't'):
            # need to use the conventional cell for tetragonal
            initial_structure = analyzer.get_conventional_standard_structure()
            a, b, c = initial_structure.lattice.abc
            # c axis of tetragonal structure not in the third direction
            if abs(a - b) > symprec:
                # a == c, rotate b to the third direction
                if abs(a - c) < symprec:
                    initial_structure.make_supercell([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
                # b == c, rotate a to the third direction
                else:
                    initial_structure.make_supercell([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        elif (self.lat_type == 'h'):
            alpha, beta, gamma = initial_structure.lattice.angles
            # c axis is not in the third direction
            if (abs(gamma - 90) < angle_tolerance):
                # alpha = 120 or 60, rotate b, c to a, b vectors
                if (abs(alpha - 90) > angle_tolerance):
                    initial_structure.make_supercell([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
                # beta = 120 or 60, rotate c, a to a, b vectors
                elif (abs(beta - 90) > angle_tolerance):
                    initial_structure.make_supercell([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
        elif (self.lat_type == 'r'):
            # need to use primitive cell for rhombohedra
            initial_structure = analyzer.get_primitive_standard_structure()
        elif (self.lat_type == 'o'):
            # need to use the conventional cell for orthorombic
            initial_structure = analyzer.get_conventional_standard_structure()
        self.initial_structure = initial_structure

    def gb_from_parameters(self, rotation_axis, rotation_angle, expand_times=4, vacuum_thickness=0.0,
                           ab_shift=[0, 0], normal=False, ratio=None, plane=None, max_search=20,
                           tol_coi=1.e-8, rm_ratio=0.7, quick_gen=False):

        """
        Args:
           rotation_axis (list): Rotation axis of GB in the form of a list of integer
                e.g.: [1, 1, 0]
           rotation_angle (float, in unit of degree): rotation angle used to generate GB.
                Make sure the angle is accurate enough. You can use the enum* functions
                in this class to extract the accurate angle.
                e.g.: The rotation angle of sigma 3 twist GB with the rotation axis
                [1, 1, 1] and GB plane (1, 1, 1) can be 60.000000000 degree.
                If you do not know the rotation angle, but know the sigma value, we have
                provide the function get_rotation_angle_from_sigma which is able to return
                all the rotation angles of sigma value you provided.
           expand_times (int): The multiple times used to expand one unit grain to larger grain.
                This is used to tune the grain length of GB to warrant that the two GBs in one
                cell do not interact with each other. Default set to 4.
           vacuum_thickness (float, in angstrom): The thickness of vacuum that you want to insert
                between two grains of the GB. Default to 0.
            ab_shift (list of float, in unit of a, b vectors of Gb): in plane shift of two grains
            normal (logic):
                determine if need to require the c axis of top grain (first transformation matrix)
                perperdicular to the surface or not.
                default to false.
            ratio (list of integers):
                    lattice axial ratio.
                    For cubic system, ratio is not needed.
                    For tetragonal system, ratio = [mu, mv], list of two integers,
                    that is, mu/mv = c2/a2. If it is irrational, set it to none.
                    For orthorhombic system, ratio = [mu, lam, mv], list of three integers,
                    that is, mu:lam:mv = c2:b2:a2. If irrational for one axis, set it to None.
                    e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
                    For rhombohedral system, ratio = [mu, mv], list of two integers,
                    that is, mu/mv is the ratio of (1+2*cos(alpha))/cos(alpha).
                    If irrational, set it to None.
                    For hexagonal system, ratio = [mu, mv], list of two integers,
                    that is, mu/mv = c2/a2. If it is irrational, set it to none.
                    This code also supplies a class method to generate the ratio from the
                    structure (get_ratio). User can also make their own approximation and
                    input the ratio directly.
            plane (list): Grain boundary plane in the form of a list of integers
                e.g.: [1, 2, 3]. If none, we set it as twist GB. The plane will be perpendicular
                to the rotation axis.
            max_search (int): max search for the GB lattice vectors that give the smallest GB
                lattice. If normal is true, also max search the GB c vector that perpendicular
                to the plane. For complex GB, if you want to speed up, you can reduce this value.
                But too small of this value may lead to error.
            tol_coi (float): tolerance to find the coincidence sites. When making approximations to
                the ratio needed to generate the GB, you probably need to increase this tolerance to
                obtain the correct number of coincidence sites. To check the number of coincidence
                sites are correct or not, you can compare the generated Gb object's sigma_from_site_prop
                with enum* sigma values (what user expected by input).
            rm_ratio (float): the criteria to remove the atoms which are too close with each other.
                rm_ratio*bond_length of bulk system is the criteria of bond length, below which the atom
                will be removed. Default to 0.7.
            quick_gen (bool): whether to quickly generate a supercell, if set to true, no need to
                find the smallest cell.

        Returns:
           Grain boundary structure (gb object).
               """
        lat_type = self.lat_type
        # if the initial structure is primitive cell in cubic system,
        # calculate the transformation matrix from its conventional cell
        # to primitive cell, basically for bcc and fcc systems.
        trans_cry = np.eye(3)
        if lat_type == 'c':
            analyzer = SpacegroupAnalyzer(self.initial_structure)
            convention_cell = analyzer.get_conventional_standard_structure()
            vol_ratio = self.initial_structure.volume / convention_cell.volume
            # bcc primitive cell, belong to cubic system
            if abs(vol_ratio - 0.5) < 1.e-3:
                trans_cry = np.array([[0.5, 0.5, -0.5], [-0.5, 0.5, 0.5], [0.5, -0.5, 0.5]])
                logger.info("Make sure this is for cubic with bcc primitive cell")
            # fcc primitive cell, belong to cubic system
            elif abs(vol_ratio - 0.25) < 1.e-3:
                trans_cry = np.array([[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]])
                logger.info("Make sure this is for cubic with fcc primitive cell")
            else:
                logger.info("Make sure this is for cubic with conventional cell")
        elif lat_type == 't':
            logger.info("Make sure this is for tetragonal system")
            if ratio is None:
                logger.info('Make sure this is for irrational c2/a2')
            elif len(ratio) != 2:
                raise RuntimeError('Tetragonal system needs correct c2/a2 ratio')
        elif lat_type == 'o':
            logger.info('Make sure this is for orthorhombic system')
            if ratio is None:
                raise RuntimeError('CSL does not exist if all axial ratios are irrational '
                                   'for an orthorhombic system')
            elif len(ratio) != 3:
                raise RuntimeError('Orthorhombic system needs correct c2:b2:a2 ratio')
        elif lat_type == 'h':
            logger.info('Make sure this is for hexagonal system')
            if ratio is None:
                logger.info('Make sure this is for irrational c2/a2')
            elif len(ratio) != 2:
                raise RuntimeError('Hexagonal system needs correct c2/a2 ratio')
        elif lat_type == 'r':
            logger.info('Make sure this is for rhombohedral system')
            if ratio is None:
                logger.info('Make sure this is for irrational (1+2*cos(alpha)/cos(alpha) ratio')
            elif len(ratio) != 2:
                raise RuntimeError('Rhombohedral system needs correct '
                                   '(1+2*cos(alpha)/cos(alpha) ratio')
        else:
            raise RuntimeError('Lattice type not implemented. This code works for cubic, '
                               'tetragonal, orthorhombic, rhombehedral, hexagonal systems')

        # transform four index notation to three index notation for hexagonal and rhombohedral
        if len(rotation_axis) == 4:
            u1 = rotation_axis[0]
            v1 = rotation_axis[1]
            w1 = rotation_axis[3]
            if lat_type.lower() == 'h':
                u = 2 * u1 + v1
                v = 2 * v1 + u1
                w = w1
                rotation_axis = [u, v, w]
            elif lat_type.lower() == 'r':
                u = 2 * u1 + v1 + w1
                v = v1 + w1 - u1
                w = w1 - 2 * v1 - u1
                rotation_axis = [u, v, w]

        # make sure gcd(rotation_axis)==1
        if reduce(gcd, rotation_axis) != 1:
            rotation_axis = [int(round(x / reduce(gcd, rotation_axis))) for x in rotation_axis]
        # transform four index notation to three index notation for plane
        if plane is not None:
            if len(plane) == 4:
                u1 = plane[0]
                v1 = plane[1]
                w1 = plane[3]
                plane = [u1, v1, w1]
        # set the plane for grain boundary when plane is None.
        if plane is None:
            if lat_type.lower() == 'c':
                plane = rotation_axis
            else:
                if lat_type.lower() == 'h':
                    if ratio is None:
                        c2_a2_ratio = 1
                    else:
                        c2_a2_ratio = ratio[0] / ratio[1]
                    metric = np.array([[1, -0.5, 0], [-0.5, 1, 0], [0, 0, c2_a2_ratio]])
                elif lat_type.lower() == 'r':
                    if ratio is None:
                        cos_alpha = 0.5
                    else:
                        cos_alpha = 1.0 / (ratio[0] / ratio[1] - 2)
                    metric = np.array([[1, cos_alpha, cos_alpha], [cos_alpha, 1, cos_alpha],
                                       [cos_alpha, cos_alpha, 1]])
                elif lat_type.lower() == 't':
                    if ratio is None:
                        c2_a2_ratio = 1
                    else:
                        c2_a2_ratio = ratio[0] / ratio[1]
                    metric = np.array([[1, 0, 0], [0, 1, 0], [0, 0, c2_a2_ratio]])
                elif lat_type.lower() == 'o':
                    for i in range(3):
                        if ratio[i] is None:
                            ratio[i] = 1
                    metric = np.array([[1, 0, 0], [0, ratio[1] / ratio[2], 0],
                                       [0, 0, ratio[0] / ratio[2]]])
                else:
                    raise RuntimeError('Lattice type has not implemented.')

                plane = np.matmul(rotation_axis, metric)
                fractions = [Fraction(x).limit_denominator() for x in plane]
                least_mul = reduce(lcm, [f.denominator for f in fractions])
                plane = [int(round(x * least_mul)) for x in plane]

        if reduce(gcd, plane) != 1:
            index = reduce(gcd, plane)
            plane = [int(round(x / index)) for x in plane]

        t1, t2 = self.get_trans_mat(r_axis=rotation_axis, angle=rotation_angle, normal=normal,
                                    trans_cry=trans_cry, lat_type=lat_type, ratio=ratio,
                                    surface=plane, max_search=max_search, quick_gen=quick_gen)

        # find the join_plane
        if lat_type.lower() != 'c':
            if lat_type.lower() == 'h':
                if ratio is None:
                    mu, mv = [1, 1]
                else:
                    mu, mv = ratio
                trans_cry1 = np.array([[1, 0, 0], [-0.5, np.sqrt(3.0) / 2.0, 0],
                                       [0, 0, np.sqrt(mu / mv)]])
            elif lat_type.lower() == 'r':
                if ratio is None:
                    c2_a2_ratio = 1
                else:
                    mu, mv = ratio
                    c2_a2_ratio = 3.0 / (2 - 6 * mv / mu)
                trans_cry1 = np.array([[0.5, np.sqrt(3.0) / 6.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)],
                                       [-0.5, np.sqrt(3.0) / 6.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)],
                                       [0, -1 * np.sqrt(3.0) / 3.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)]])
            else:
                if lat_type.lower() == 't':
                    if ratio is None:
                        mu, mv = [1, 1]
                    else:
                        mu, mv = ratio
                    lam = mv
                elif lat_type.lower() == 'o':
                    new_ratio = [1 if v is None else v for v in ratio]
                    mu, lam, mv = new_ratio
                trans_cry1 = np.array([[1, 0, 0], [0, np.sqrt(lam / mv), 0], [0, 0, np.sqrt(mu / mv)]])
        else:
            trans_cry1 = trans_cry
        grain_matrix = np.dot(t2, trans_cry1)
        plane_init = np.cross(grain_matrix[0], grain_matrix[1])
        if lat_type.lower() != 'c':
            plane_init = np.dot(plane_init, trans_cry1.T)
        join_plane = self.vec_to_surface(plane_init)

        parent_structure = self.initial_structure.copy()
        # calculate the bond_length in bulk system.
        if len(parent_structure) == 1:
            temp_str = parent_structure.copy()
            temp_str.make_supercell([1, 1, 2])
            distance = temp_str.distance_matrix
        else:
            distance = parent_structure.distance_matrix
        bond_length = np.min(distance[np.nonzero(distance)])

        # top grain
        top_grain = fix_pbc(parent_structure * t1)

        # obtain the smallest oriended cell
        if normal and not quick_gen:
            t_temp = self.get_trans_mat(r_axis=rotation_axis, angle=rotation_angle, normal=False,
                                        trans_cry=trans_cry, lat_type=lat_type, ratio=ratio,
                                        surface=plane, max_search=max_search)
            oriended_unit_cell = fix_pbc(parent_structure * t_temp[0])
            t_matrix = oriended_unit_cell.lattice.matrix
            normal_v_plane = np.cross(t_matrix[0], t_matrix[1])
            unit_normal_v = normal_v_plane / np.linalg.norm(normal_v_plane)
            unit_ab_adjust = (t_matrix[2] - np.dot(unit_normal_v, t_matrix[2]) * unit_normal_v) \
                             / np.dot(unit_normal_v, t_matrix[2])
        else:
            oriended_unit_cell = top_grain.copy()
            unit_ab_adjust = 0.0

        # bottom grain, using top grain's lattice matrix
        bottom_grain = fix_pbc(parent_structure * t2, top_grain.lattice.matrix)

        # label both grains with 'top','bottom','top_incident','bottom_incident'
        n_sites = top_grain.num_sites
        t_and_b = Structure(top_grain.lattice, top_grain.species + bottom_grain.species,
                            list(top_grain.frac_coords) + list(bottom_grain.frac_coords))
        t_and_b_dis = t_and_b.lattice.get_all_distances(t_and_b.frac_coords[0:n_sites],
                                                        t_and_b.frac_coords[n_sites:n_sites * 2])
        index_incident = np.nonzero(t_and_b_dis < np.min(t_and_b_dis) + tol_coi)

        top_labels = []
        for i in range(n_sites):
            if i in index_incident[0]:
                top_labels.append('top_incident')
            else:
                top_labels.append('top')
        bottom_labels = []
        for i in range(n_sites):
            if i in index_incident[1]:
                bottom_labels.append('bottom_incident')
            else:
                bottom_labels.append('bottom')
        top_grain = Structure(Lattice(top_grain.lattice.matrix), top_grain.species,
                              top_grain.frac_coords, site_properties={'grain_label': top_labels})
        bottom_grain = Structure(Lattice(bottom_grain.lattice.matrix), bottom_grain.species,
                                 bottom_grain.frac_coords, site_properties={'grain_label': bottom_labels})

        # expand both grains
        top_grain.make_supercell([1, 1, expand_times])
        bottom_grain.make_supercell([1, 1, expand_times])
        top_grain = fix_pbc(top_grain)
        bottom_grain = fix_pbc(bottom_grain)

        # determine the top-grain location.
        edge_b = 1.0 - max(bottom_grain.frac_coords[:, 2])
        edge_t = 1.0 - max(top_grain.frac_coords[:, 2])
        c_adjust = (edge_t - edge_b) / 2.0

        # construct all species
        all_species = []
        all_species.extend([site.specie for site in bottom_grain])
        all_species.extend([site.specie for site in top_grain])

        half_lattice = top_grain.lattice
        # calculate translation vector, perpendicular to the plane
        normal_v_plane = np.cross(half_lattice.matrix[0], half_lattice.matrix[1])
        unit_normal_v = normal_v_plane / np.linalg.norm(normal_v_plane)
        translation_v = unit_normal_v * vacuum_thickness

        # construct the final lattice
        whole_matrix_no_vac = np.array(half_lattice.matrix)
        whole_matrix_no_vac[2] = half_lattice.matrix[2] * 2
        whole_matrix_with_vac = whole_matrix_no_vac.copy()
        whole_matrix_with_vac[2] = whole_matrix_no_vac[2] + translation_v * 2
        whole_lat = Lattice(whole_matrix_with_vac)

        # construct the coords, move top grain with translation_v
        all_coords = []
        grain_labels = bottom_grain.site_properties['grain_label'] \
                       + top_grain.site_properties['grain_label']
        for site in bottom_grain:
            all_coords.append(site.coords)
        for site in top_grain:
            all_coords.append(site.coords + half_lattice.matrix[2] * (1 + c_adjust) +
                              unit_ab_adjust * np.linalg.norm(half_lattice.matrix[2] * (1 + c_adjust)) +
                              translation_v + ab_shift[0] * whole_matrix_with_vac[0] +
                              ab_shift[1] * whole_matrix_with_vac[1])

        gb_with_vac = Structure(whole_lat, all_species, all_coords,
                                coords_are_cartesian=True,
                                site_properties={'grain_label': grain_labels})
        # merge closer atoms. extract near gb atoms.
        cos_c_norm_plane = np.dot(unit_normal_v, whole_matrix_with_vac[2]) / whole_lat.c
        range_c_len = abs(bond_length / cos_c_norm_plane / whole_lat.c)
        sites_near_gb = []
        sites_away_gb = []
        for site in gb_with_vac.sites:
            if site.frac_coords[2] < range_c_len or site.frac_coords[2] > 1 - range_c_len \
                    or (site.frac_coords[2] > 0.5 - range_c_len and site.frac_coords[2] < 0.5 + range_c_len):
                sites_near_gb.append(site)
            else:
                sites_away_gb.append(site)
        if len(sites_near_gb) >= 1:
            s_near_gb = Structure.from_sites(sites_near_gb)
            s_near_gb.merge_sites(tol=bond_length * rm_ratio, mode='d')
            all_sites = sites_away_gb + s_near_gb.sites
            gb_with_vac = Structure.from_sites(all_sites)

        return GrainBoundary(whole_lat, gb_with_vac.species, gb_with_vac.cart_coords, rotation_axis,
                             rotation_angle, plane, join_plane, self.initial_structure,
                             vacuum_thickness, ab_shift, site_properties=gb_with_vac.site_properties,
                             oriented_unit_cell=oriended_unit_cell,
                             coords_are_cartesian=True)

    def get_ratio(self, max_denominator=5, index_none=None):
        """
        find the axial ratio needed for GB generator input.
        Args:
            max_denominator (int): the maximum denominator for
                the computed ratio, default to be 5.
            index_none (int): specify the irrational axis.
                0-a, 1-b, 2-c. Only may be needed for orthorombic system.
        Returns:
               axial ratio needed for GB generator (list of integers).

        """
        structure = self.initial_structure
        lat_type = self.lat_type
        if lat_type == 't' or lat_type == 'h':
            # For tetragonal and hexagonal system, ratio = c2 / a2.
            a, c = (structure.lattice.a, structure.lattice.c)
            if c > a:
                frac = Fraction(c ** 2 / a ** 2).limit_denominator(max_denominator)
                ratio = [frac.numerator, frac.denominator]
            else:
                frac = Fraction(a ** 2 / c ** 2).limit_denominator(max_denominator)
                ratio = [frac.denominator, frac.numerator]
        elif lat_type == 'r':
            # For rhombohedral system, ratio = (1 + 2 * cos(alpha)) / cos(alpha).
            cos_alpha = cos(structure.lattice.alpha / 180 * np.pi)
            frac = Fraction((1 + 2 * cos_alpha) / cos_alpha).limit_denominator(max_denominator)
            ratio = [frac.numerator, frac.denominator]
        elif lat_type == 'o':
            # For orthorhombic system, ratio = c2:b2:a2.If irrational for one axis, set it to None.
            ratio = [None] * 3
            lat = (structure.lattice.c, structure.lattice.b, structure.lattice.a)
            index = [0, 1, 2]
            if index_none is None:
                min_index = np.argmin(lat)
                index.pop(min_index)
                frac1 = Fraction(lat[index[0]] ** 2 / lat[min_index] ** 2).limit_denominator(max_denominator)
                frac2 = Fraction(lat[index[1]] ** 2 / lat[min_index] ** 2).limit_denominator(max_denominator)
                com_lcm = lcm(frac1.denominator, frac2.denominator)
                ratio[min_index] = com_lcm
                ratio[index[0]] = frac1.numerator * int(round((com_lcm / frac1.denominator)))
                ratio[index[1]] = frac2.numerator * int(round((com_lcm / frac2.denominator)))
            else:
                index.pop(index_none)
                if (lat[index[0]] > lat[index[1]]):
                    frac = Fraction(lat[index[0]] ** 2 / lat[index[1]] ** 2).limit_denominator(max_denominator)
                    ratio[index[0]] = frac.numerator
                    ratio[index[1]] = frac.denominator
                else:
                    frac = Fraction(lat[index[1]] ** 2 / lat[index[0]] ** 2).limit_denominator(max_denominator)
                    ratio[index[1]] = frac.numerator
                    ratio[index[0]] = frac.denominator
        elif lat_type == 'c':
            raise RuntimeError('Cubic system does not need axial ratio.')
        else:
            raise RuntimeError('Lattice type not implemented.')
        return ratio

    @staticmethod
    def get_trans_mat(r_axis, angle, normal=False, trans_cry=np.eye(3), lat_type='c',
                      ratio=None, surface=None, max_search=20, quick_gen=False):
        """
        Find the two transformation matrix for each grain from given rotation axis,
        GB plane, rotation angle and corresponding ratio (see explanation for ratio
        below).
        The structure of each grain can be obtained by applying the corresponding
        transformation matrix to the conventional cell.
        The algorithm for this code is from reference, Acta Cryst, A32,783(1976).

        Args:
            r_axis (list of three integers, e.g. u, v, w
                    or four integers, e.g. u, v, t, w for hex/rho system only):
                    the rotation axis of the grain boundary.
            angle (float, in unit of degree) :
                    the rotation angle of the grain boundary
            normal (logic):
                    determine if need to require the c axis of one grain associated with
                    the first transformation matrix perperdicular to the surface or not.
                    default to false.
            trans_cry (3 by 3 array):
                    if the structure given are primitive cell in cubic system, e.g.
                    bcc or fcc system, trans_cry is the transformation matrix from its
                    conventional cell to the primitive cell.
            lat_type ( one character):
                    'c' or 'C': cubic system
                     't' or 'T': tetragonal system
                     'o' or 'O': orthorhombic system
                     'h' or 'H': hexagonal system
                     'r' or 'R': rhombohedral system
                     default to cubic system
            ratio (list of integers):
                    lattice axial ratio.
                    For cubic system, ratio is not needed.
                    For tetragonal system, ratio = [mu, mv], list of two integers,
                    that is, mu/mv = c2/a2. If it is irrational, set it to none.
                    For orthorhombic system, ratio = [mu, lam, mv], list of three integers,
                    that is, mu:lam:mv = c2:b2:a2. If irrational for one axis, set it to None.
                    e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
                    For rhombohedral system, ratio = [mu, mv], list of two integers,
                    that is, mu/mv is the ratio of (1+2*cos(alpha)/cos(alpha).
                    If irrational, set it to None.
                    For hexagonal system, ratio = [mu, mv], list of two integers,
                    that is, mu/mv = c2/a2. If it is irrational, set it to none.
            surface (list of three integers, e.g. h, k, l
                     or four integers, e.g. h, k, i, l for hex/rho system only):
                    the miller index of grain boundary plane, with the format of [h,k,l]
                    if surface is not given, the default is perpendicular to r_axis, which is
                    a twist grain boundary.
            max_search (int): max search for the GB lattice vectors that give the smallest GB
                lattice. If normal is true, also max search the GB c vector that perpendicular
                to the plane.
            quick_gen (bool): whether to quickly generate a supercell, if set to true, no need to
                find the smallest cell.
        Returns:
            t1 (3 by 3 integer array):
                    The transformation array for one grain.
            t2 (3 by 3 integer array):
                    The transformation array for the other grain
        """
        # transform four index notation to three index notation
        if len(r_axis) == 4:
            u1 = r_axis[0]
            v1 = r_axis[1]
            w1 = r_axis[3]
            if lat_type.lower() == 'h':
                u = 2 * u1 + v1
                v = 2 * v1 + u1
                w = w1
                r_axis = [u, v, w]
            elif lat_type.lower() == 'r':
                u = 2 * u1 + v1 + w1
                v = v1 + w1 - u1
                w = w1 - 2 * v1 - u1
                r_axis = [u, v, w]

        # make sure gcd(r_axis)==1
        if reduce(gcd, r_axis) != 1:
            r_axis = [int(round(x / reduce(gcd, r_axis))) for x in r_axis]

        if surface is not None:
            if len(surface) == 4:
                u1 = surface[0]
                v1 = surface[1]
                w1 = surface[3]
                surface = [u1, v1, w1]
        # set the surface for grain boundary.
        if surface is None:
            if lat_type.lower() == 'c':
                surface = r_axis
            else:
                if lat_type.lower() == 'h':
                    if ratio is None:
                        c2_a2_ratio = 1
                    else:
                        c2_a2_ratio = ratio[0] / ratio[1]
                    metric = np.array([[1, -0.5, 0], [-0.5, 1, 0], [0, 0, c2_a2_ratio]])
                elif lat_type.lower() == 'r':
                    if ratio is None:
                        cos_alpha = 0.5
                    else:
                        cos_alpha = 1.0 / (ratio[0] / ratio[1] - 2)
                    metric = np.array([[1, cos_alpha, cos_alpha], [cos_alpha, 1, cos_alpha],
                                       [cos_alpha, cos_alpha, 1]])
                elif lat_type.lower() == 't':
                    if ratio is None:
                        c2_a2_ratio = 1
                    else:
                        c2_a2_ratio = ratio[0] / ratio[1]
                    metric = np.array([[1, 0, 0], [0, 1, 0], [0, 0, c2_a2_ratio]])
                elif lat_type.lower() == 'o':
                    for i in range(3):
                        if ratio[i] is None:
                            ratio[i] = 1
                    metric = np.array([[1, 0, 0], [0, ratio[1] / ratio[2], 0],
                                       [0, 0, ratio[0] / ratio[2]]])
                else:
                    raise RuntimeError('Lattice type has not implemented.')

                surface = np.matmul(r_axis, metric)
                fractions = [Fraction(x).limit_denominator() for x in surface]
                least_mul = reduce(lcm, [f.denominator for f in fractions])
                surface = [int(round(x * least_mul)) for x in surface]

        if reduce(gcd, surface) != 1:
            index = reduce(gcd, surface)
            surface = [int(round(x / index)) for x in surface]

        if lat_type.lower() == 'h':
            # set the value for u,v,w,mu,mv,m,n,d,x
            # check the reference for the meaning of these parameters
            u, v, w = r_axis
            # make sure mu, mv are coprime integers.
            if ratio is None:
                mu, mv = [1, 1]
                if w != 0:
                    if u != 0 or (v != 0):
                        raise RuntimeError('For irrational c2/a2, CSL only exist for [0,0,1] '
                                           'or [u,v,0] and m = 0')
            else:
                mu, mv = ratio
            if gcd(mu, mv) != 1:
                temp = gcd(mu, mv)
                mu = int(round(mu / temp))
                mv = int(round(mv / temp))
            d = (u ** 2 + v ** 2 - u * v) * mv + w ** 2 * mu
            if abs(angle - 180.0) < 1.e0:
                m = 0
                n = 1
            else:
                fraction = Fraction(np.tan(angle / 2 / 180.0 * np.pi) /
                                    np.sqrt(float(d) / 3.0 / mu)).limit_denominator()
                m = fraction.denominator
                n = fraction.numerator

            # construct the rotation matrix, check reference for details
            r_list = [(u ** 2 * mv - v ** 2 * mv - w ** 2 * mu) * n ** 2 +
                      2 * w * mu * m * n + 3 * mu * m ** 2,
                      (2 * v - u) * u * mv * n ** 2 - 4 * w * mu * m * n,
                      2 * u * w * mu * n ** 2 + 2 * (2 * v - u) * mu * m * n,
                      (2 * u - v) * v * mv * n ** 2 + 4 * w * mu * m * n,
                      (v ** 2 * mv - u ** 2 * mv - w ** 2 * mu) * n ** 2 -
                      2 * w * mu * m * n + 3 * mu * m ** 2,
                      2 * v * w * mu * n ** 2 - 2 * (2 * u - v) * mu * m * n,
                      (2 * u - v) * w * mv * n ** 2 - 3 * v * mv * m * n,
                      (2 * v - u) * w * mv * n ** 2 + 3 * u * mv * m * n,
                      (w ** 2 * mu - u ** 2 * mv - v ** 2 * mv + u * v * mv) *
                      n ** 2 + 3 * mu * m ** 2]
            m = -1 * m
            r_list_inv = [(u ** 2 * mv - v ** 2 * mv - w ** 2 * mu) * n ** 2 +
                          2 * w * mu * m * n + 3 * mu * m ** 2,
                          (2 * v - u) * u * mv * n ** 2 - 4 * w * mu * m * n,
                          2 * u * w * mu * n ** 2 + 2 * (2 * v - u) * mu * m * n,
                          (2 * u - v) * v * mv * n ** 2 + 4 * w * mu * m * n,
                          (v ** 2 * mv - u ** 2 * mv - w ** 2 * mu) * n ** 2 -
                          2 * w * mu * m * n + 3 * mu * m ** 2,
                          2 * v * w * mu * n ** 2 - 2 * (2 * u - v) * mu * m * n,
                          (2 * u - v) * w * mv * n ** 2 - 3 * v * mv * m * n,
                          (2 * v - u) * w * mv * n ** 2 + 3 * u * mv * m * n,
                          (w ** 2 * mu - u ** 2 * mv - v ** 2 * mv + u * v * mv) *
                          n ** 2 + 3 * mu * m ** 2]
            m = -1 * m
            F = 3 * mu * m ** 2 + d * n ** 2
            all_list = r_list + r_list_inv + [F]
            com_fac = reduce(gcd, all_list)
            sigma = F / com_fac
            r_matrix = (np.array(r_list) / com_fac / sigma).reshape(3, 3)
        elif lat_type.lower() == 'r':
            # set the value for u,v,w,mu,mv,m,n,d
            # check the reference for the meaning of these parameters
            u, v, w = r_axis
            # make sure mu, mv are coprime integers.
            if ratio is None:
                mu, mv = [1, 1]
                if u + v + w != 0:
                    if u != v or u != w:
                        raise RuntimeError('For irrational ratio_alpha, CSL only exist for [1,1,1]'
                                           'or [u, v, -(u+v)] and m =0')
            else:
                mu, mv = ratio
            if gcd(mu, mv) != 1:
                temp = gcd(mu, mv)
                mu = int(round(mu / temp))
                mv = int(round(mv / temp))
            d = (u ** 2 + v ** 2 + w ** 2) * (mu - 2 * mv) + \
                2 * mv * (v * w + w * u + u * v)
            if abs(angle - 180.0) < 1.e0:
                m = 0
                n = 1
            else:
                fraction = Fraction(np.tan(angle / 2 / 180.0 * np.pi) /
                                    np.sqrt(float(d) / mu)).limit_denominator()
                m = fraction.denominator
                n = fraction.numerator

            # construct the rotation matrix, check reference for details
            r_list = [(mu - 2 * mv) * (u ** 2 - v ** 2 - w ** 2) * n ** 2 +
                      2 * mv * (v - w) * m * n - 2 * mv * v * w * n ** 2 +
                      mu * m ** 2,
                      2 * (mv * u * n * (w * n + u * n - m) - (mu - mv) *
                           m * w * n + (mu - 2 * mv) * u * v * n ** 2),
                      2 * (mv * u * n * (v * n + u * n + m) + (mu - mv) *
                           m * v * n + (mu - 2 * mv) * w * u * n ** 2),
                      2 * (mv * v * n * (w * n + v * n + m) + (mu - mv) *
                           m * w * n + (mu - 2 * mv) * u * v * n ** 2),
                      (mu - 2 * mv) * (v ** 2 - w ** 2 - u ** 2) * n ** 2 +
                      2 * mv * (w - u) * m * n - 2 * mv * u * w * n ** 2 +
                      mu * m ** 2,
                      2 * (mv * v * n * (v * n + u * n - m) - (mu - mv) *
                           m * u * n + (mu - 2 * mv) * w * v * n ** 2),
                      2 * (mv * w * n * (w * n + v * n - m) - (mu - mv) *
                           m * v * n + (mu - 2 * mv) * w * u * n ** 2),
                      2 * (mv * w * n * (w * n + u * n + m) + (mu - mv) *
                           m * u * n + (mu - 2 * mv) * w * v * n ** 2),
                      (mu - 2 * mv) * (w ** 2 - u ** 2 - v ** 2) * n ** 2 +
                      2 * mv * (u - v) * m * n - 2 * mv * u * v * n ** 2 +
                      mu * m ** 2]
            m = -1 * m
            r_list_inv = [(mu - 2 * mv) * (u ** 2 - v ** 2 - w ** 2) * n ** 2 +
                          2 * mv * (v - w) * m * n - 2 * mv * v * w * n ** 2 +
                          mu * m ** 2,
                          2 * (mv * u * n * (w * n + u * n - m) - (mu - mv) *
                               m * w * n + (mu - 2 * mv) * u * v * n ** 2),
                          2 * (mv * u * n * (v * n + u * n + m) + (mu - mv) *
                               m * v * n + (mu - 2 * mv) * w * u * n ** 2),
                          2 * (mv * v * n * (w * n + v * n + m) + (mu - mv) *
                               m * w * n + (mu - 2 * mv) * u * v * n ** 2),
                          (mu - 2 * mv) * (v ** 2 - w ** 2 - u ** 2) * n ** 2 +
                          2 * mv * (w - u) * m * n - 2 * mv * u * w * n ** 2 +
                          mu * m ** 2,
                          2 * (mv * v * n * (v * n + u * n - m) - (mu - mv) *
                               m * u * n + (mu - 2 * mv) * w * v * n ** 2),
                          2 * (mv * w * n * (w * n + v * n - m) - (mu - mv) *
                               m * v * n + (mu - 2 * mv) * w * u * n ** 2),
                          2 * (mv * w * n * (w * n + u * n + m) + (mu - mv) *
                               m * u * n + (mu - 2 * mv) * w * v * n ** 2),
                          (mu - 2 * mv) * (w ** 2 - u ** 2 - v ** 2) * n ** 2 +
                          2 * mv * (u - v) * m * n - 2 * mv * u * v * n ** 2 +
                          mu * m ** 2]
            m = -1 * m
            F = mu * m ** 2 + d * n ** 2
            all_list = r_list_inv + r_list + [F]
            com_fac = reduce(gcd, all_list)
            sigma = F / com_fac
            r_matrix = (np.array(r_list) / com_fac / sigma).reshape(3, 3)
        else:
            u, v, w = r_axis
            if lat_type.lower() == 'c':
                mu = 1
                lam = 1
                mv = 1
            elif lat_type.lower() == 't':
                if ratio is None:
                    mu, mv = [1, 1]
                    if w != 0:
                        if u != 0 or (v != 0):
                            raise RuntimeError('For irrational c2/a2, CSL only exist for [0,0,1] '
                                               'or [u,v,0] and m = 0')
                else:
                    mu, mv = ratio
                lam = mv
            elif lat_type.lower() == 'o':
                if None in ratio:
                    mu, lam, mv = ratio
                    non_none = [i for i in ratio if i is not None]
                    if len(non_none) < 2:
                        raise RuntimeError('No CSL exist for two irrational numbers')
                    non1, non2 = non_none
                    if mu is None:
                        lam = non1
                        mv = non2
                        mu = 1
                        if w != 0:
                            if u != 0 or (v != 0):
                                raise RuntimeError('For irrational c2, CSL only exist for [0,0,1] '
                                                   'or [u,v,0] and m = 0')
                    elif lam is None:
                        mu = non1
                        mv = non2
                        lam = 1
                        if v != 0:
                            if u != 0 or (w != 0):
                                raise RuntimeError('For irrational b2, CSL only exist for [0,1,0] '
                                                   'or [u,0,w] and m = 0')
                    elif mv is None:
                        mu = non1
                        lam = non2
                        mv = 1
                        if u != 0:
                            if w != 0 or (v != 0):
                                raise RuntimeError('For irrational a2, CSL only exist for [1,0,0] '
                                                   'or [0,v,w] and m = 0')
                else:
                    mu, lam, mv = ratio
                    if u == 0 and v == 0:
                        mu = 1
                    if u == 0 and w == 0:
                        lam = 1
                    if v == 0 and w == 0:
                        mv = 1

            # make sure mu, lambda, mv are coprime integers.
            if reduce(gcd, [mu, lam, mv]) != 1:
                temp = reduce(gcd, [mu, lam, mv])
                mu = int(round(mu / temp))
                mv = int(round(mv / temp))
                lam = int(round(lam / temp))
            d = (mv * u ** 2 + lam * v ** 2) * mv + w ** 2 * mu * mv
            if abs(angle - 180.0) < 1.e0:
                m = 0
                n = 1
            else:
                fraction = Fraction(np.tan(angle / 2 / 180.0 * np.pi) /
                                    np.sqrt(d / mu / lam)).limit_denominator()
                m = fraction.denominator
                n = fraction.numerator
            r_list = [(u ** 2 * mv * mv - lam * v ** 2 * mv -
                       w ** 2 * mu * mv) * n ** 2 + lam * mu * m ** 2,
                      2 * lam * (v * u * mv * n ** 2 - w * mu * m * n),
                      2 * mu * (u * w * mv * n ** 2 + v * lam * m * n),
                      2 * mv * (u * v * mv * n ** 2 + w * mu * m * n),
                      (v ** 2 * mv * lam - u ** 2 * mv * mv -
                       w ** 2 * mu * mv) * n ** 2 + lam * mu * m ** 2,
                      2 * mv * mu * (v * w * n ** 2 - u * m * n),
                      2 * mv * (u * w * mv * n ** 2 - v * lam * m * n),
                      2 * lam * mv * (v * w * n ** 2 + u * m * n),
                      (w ** 2 * mu * mv - u ** 2 * mv * mv -
                       v ** 2 * mv * lam) * n ** 2 + lam * mu * m ** 2]
            m = -1 * m
            r_list_inv = [(u ** 2 * mv * mv - lam * v ** 2 * mv -
                           w ** 2 * mu * mv) * n ** 2 + lam * mu * m ** 2,
                          2 * lam * (v * u * mv * n ** 2 - w * mu * m * n),
                          2 * mu * (u * w * mv * n ** 2 + v * lam * m * n),
                          2 * mv * (u * v * mv * n ** 2 + w * mu * m * n),
                          (v ** 2 * mv * lam - u ** 2 * mv * mv -
                           w ** 2 * mu * mv) * n ** 2 + lam * mu * m ** 2,
                          2 * mv * mu * (v * w * n ** 2 - u * m * n),
                          2 * mv * (u * w * mv * n ** 2 - v * lam * m * n),
                          2 * lam * mv * (v * w * n ** 2 + u * m * n),
                          (w ** 2 * mu * mv - u ** 2 * mv * mv -
                           v ** 2 * mv * lam) * n ** 2 + lam * mu * m ** 2]
            m = -1 * m
            F = mu * lam * m ** 2 + d * n ** 2
            all_list = r_list + r_list_inv + [F]
            com_fac = reduce(gcd, all_list)
            sigma = F / com_fac
            r_matrix = (np.array(r_list) / com_fac / sigma).reshape(3, 3)

        if (sigma > 1000):
            raise RuntimeError('Sigma >1000 too large. Are you sure what you are doing, '
                               'Please check the GB if exist')
        # transform surface, r_axis, r_matrix in terms of primitive lattice
        surface = np.matmul(surface, np.transpose(trans_cry))
        fractions = [Fraction(x).limit_denominator() for x in surface]
        least_mul = reduce(lcm, [f.denominator for f in fractions])
        surface = [int(round(x * least_mul)) for x in surface]
        if reduce(gcd, surface) != 1:
            index = reduce(gcd, surface)
            surface = [int(round(x / index)) for x in surface]
        r_axis = np.rint(np.matmul(r_axis, np.linalg.inv(trans_cry))).astype(int)
        if reduce(gcd, r_axis) != 1:
            r_axis = [int(round(x / reduce(gcd, r_axis))) for x in r_axis]
        r_matrix = np.dot(np.dot(np.linalg.inv(trans_cry.T), r_matrix), trans_cry.T)
        # set one vector of the basis to the rotation axis direction, and
        # obtain the corresponding transform matrix
        eye = np.eye(3, dtype=np.int)
        for h in range(3):
            if abs(r_axis[h]) != 0:
                eye[h] = np.array(r_axis)
                k = h + 1 if h + 1 < 3 else abs(2 - h)
                l = h + 2 if h + 2 < 3 else abs(1 - h)
                break
        trans = eye.T
        new_rot = np.array(r_matrix)

        # with the rotation matrix to construct the CSL lattice, check reference for details
        fractions = [Fraction(x).limit_denominator() for x in new_rot[:, k]]
        least_mul = reduce(lcm, [f.denominator for f in fractions])
        scale = np.zeros((3, 3))
        scale[h, h] = 1
        scale[k, k] = least_mul
        scale[l, l] = sigma / least_mul
        for i in range(least_mul):
            check_int = i * new_rot[:, k] + (sigma / least_mul) * new_rot[:, l]
            if all([np.round(x, 5).is_integer() for x in list(check_int)]):
                n_final = i
                break
        try:
            n_final
        except NameError:
            raise RuntimeError('Something is wrong. Check if this GB exists or not')
        scale[k, l] = n_final
        # each row of mat_csl is the CSL lattice vector
        csl_init = np.rint(np.dot(np.dot(r_matrix, trans), scale)).astype(int).T
        if abs(r_axis[h]) > 1:
            csl_init = GrainBoundaryGenerator.reduce_mat(np.array(csl_init), r_axis[h], r_matrix)
        csl = np.rint(Lattice(csl_init).get_niggli_reduced_lattice().matrix).astype(int)

        # find the best slab supercell in terms of the conventional cell from the csl lattice,
        # which is the transformation matrix

        # now trans_cry is the transformation matrix from crystal to cartesian coordinates.
        # for cubic, do not need to change.
        if lat_type.lower() != 'c':
            if lat_type.lower() == 'h':
                trans_cry = np.array([[1, 0, 0], [-0.5, np.sqrt(3.0) / 2.0, 0],
                                      [0, 0, np.sqrt(mu / mv)]])
            elif lat_type.lower() == 'r':
                if ratio is None:
                    c2_a2_ratio = 1
                else:
                    c2_a2_ratio = 3.0 / (2 - 6 * mv / mu)
                trans_cry = np.array([[0.5, np.sqrt(3.0) / 6.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)],
                                      [-0.5, np.sqrt(3.0) / 6.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)],
                                      [0, -1 * np.sqrt(3.0) / 3.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)]])
            else:
                trans_cry = np.array([[1, 0, 0], [0, np.sqrt(lam / mv), 0], [0, 0, np.sqrt(mu / mv)]])
        t1_final = GrainBoundaryGenerator.slab_from_csl(csl, surface, normal, trans_cry, max_search=max_search,
                                                        quick_gen=quick_gen)
        t2_final = np.array(np.rint(np.dot(t1_final, np.linalg.inv(r_matrix.T)))).astype(int)
        return t1_final, t2_final

    @staticmethod
    def enum_sigma_cubic(cutoff, r_axis):
        """
        Find all possible sigma values and corresponding rotation angles
        within a sigma value cutoff with known rotation axis in cubic system.
        The algorithm for this code is from reference, Acta Cryst, A40,108(1984)
        Args:
            cutoff (integer): the cutoff of sigma values.
            r_axis (list of three integers, e.g. u, v, w):
                    the rotation axis of the grain boundary, with the format of [u,v,w].
        Returns:
            sigmas (dict):
                    dictionary with keys as the possible integer sigma values
                    and values as list of the possible rotation angles to the
                    corresponding sigma values.
                    e.g. the format as
                    {sigma1: [angle11,angle12,...], sigma2: [angle21, angle22,...],...}
                    Note: the angles are the rotation angles of one grain respect to
                    the other grain.
                    When generate the microstructures of the grain boundary using these angles,
                    you need to analyze the symmetry of the structure. Different angles may
                    result in equivalent microstructures.

        """
        sigmas = {}
        # make sure gcd(r_axis)==1
        if reduce(gcd, r_axis) != 1:
            r_axis = [int(round(x / reduce(gcd, r_axis))) for x in r_axis]

        # count the number of odds in r_axis
        odd_r = len(list(filter(lambda x: x % 2 == 1, r_axis)))
        # Compute the max n we need to enumerate.
        if odd_r == 3:
            a_max = 4
        elif odd_r == 0:
            a_max = 1
        else:
            a_max = 2
        n_max = int(np.sqrt(cutoff * a_max / sum(np.array(r_axis) ** 2)))
        # enumerate all possible n, m to give possible sigmas within the cutoff.
        for n_loop in range(1, n_max + 1):
            n = n_loop
            m_max = int(np.sqrt(cutoff * a_max - n ** 2 * sum(np.array(r_axis) ** 2)))
            for m in range(0, m_max + 1):
                if gcd(m, n) == 1 or m == 0:
                    if m == 0:
                        n = 1
                    else:
                        n = n_loop
                    # construct the quadruple [m, U,V,W], count the number of odds in
                    # quadruple to determine the parameter a, refer to the reference
                    quadruple = [m] + [x * n for x in r_axis]
                    odd_qua = len(list(filter(lambda x: x % 2 == 1, quadruple)))
                    if odd_qua == 4:
                        a = 4
                    elif odd_qua == 2:
                        a = 2
                    else:
                        a = 1
                    sigma = int(round((m ** 2 + n ** 2 * sum(np.array(r_axis) ** 2)) / a))
                    if (sigma <= cutoff) and (sigma > 1):
                        if sigma not in list(sigmas.keys()):
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n * np.sqrt(sum(np.array(r_axis) ** 2)) / m) \
                                        / np.pi * 180
                            sigmas[sigma] = [angle]
                        else:
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n * np.sqrt(sum(np.array(r_axis) ** 2)) / m) \
                                        / np.pi * 180
                            if angle not in sigmas[sigma]:
                                sigmas[sigma].append(angle)
        return sigmas

    @staticmethod
    def enum_sigma_hex(cutoff, r_axis, c2_a2_ratio):
        """
        Find all possible sigma values and corresponding rotation angles
        within a sigma value cutoff with known rotation axis in hexagonal system.
        The algorithm for this code is from reference, Acta Cryst, A38,550(1982)

        Args:
            cutoff (integer): the cutoff of sigma values.
            r_axis (list of three integers, e.g. u, v, w
                    or four integers, e.g. u, v, t, w):
                    the rotation axis of the grain boundary.
            c2_a2_ratio (list of two integers, e.g. mu, mv):
                    mu/mv is the square of the hexagonal axial ratio, which is rational
                    number. If irrational, set c2_a2_ratio = None
        Returns:
            sigmas (dict):
                    dictionary with keys as the possible integer sigma values
                    and values as list of the possible rotation angles to the
                    corresponding sigma values.
                    e.g. the format as
                    {sigma1: [angle11,angle12,...], sigma2: [angle21, angle22,...],...}
                    Note: the angles are the rotation angle of one grain respect to the
                    other grain.
                    When generate the microstructure of the grain boundary using these
                    angles, you need to analyze the symmetry of the structure. Different
                    angles may result in equivalent microstructures.

        """
        sigmas = {}
        # make sure gcd(r_axis)==1
        if reduce(gcd, r_axis) != 1:
            r_axis = [int(round(x / reduce(gcd, r_axis))) for x in r_axis]
        # transform four index notation to three index notation
        if len(r_axis) == 4:
            u1 = r_axis[0]
            v1 = r_axis[1]
            w1 = r_axis[3]
            u = 2 * u1 + v1
            v = 2 * v1 + u1
            w = w1
        else:
            u, v, w = r_axis

        # make sure mu, mv are coprime integers.
        if c2_a2_ratio is None:
            mu, mv = [1, 1]
            if w != 0:
                if u != 0 or (v != 0):
                    raise RuntimeError('For irrational c2/a2, CSL only exist for [0,0,1] '
                                       'or [u,v,0] and m = 0')
        else:
            mu, mv = c2_a2_ratio
            if gcd(mu, mv) != 1:
                temp = gcd(mu, mv)
                mu = int(round(mu / temp))
                mv = int(round(mv / temp))

        # refer to the meaning of d in reference
        d = (u ** 2 + v ** 2 - u * v) * mv + w ** 2 * mu

        # Compute the max n we need to enumerate.
        n_max = int(np.sqrt((cutoff * 12 * mu * mv) / abs(d)))

        # Enumerate all possible n, m to give possible sigmas within the cutoff.
        for n in range(1, n_max + 1):
            if (c2_a2_ratio is None) and w == 0:
                m_max = 0
            else:
                m_max = int(np.sqrt((cutoff * 12 * mu * mv - n ** 2 * d) / (3 * mu)))
            for m in range(0, m_max + 1):
                if gcd(m, n) == 1 or m == 0:
                    # construct the rotation matrix, refer to the reference
                    R_list = [(u ** 2 * mv - v ** 2 * mv - w ** 2 * mu) * n ** 2 +
                              2 * w * mu * m * n + 3 * mu * m ** 2,
                              (2 * v - u) * u * mv * n ** 2 - 4 * w * mu * m * n,
                              2 * u * w * mu * n ** 2 + 2 * (2 * v - u) * mu * m * n,
                              (2 * u - v) * v * mv * n ** 2 + 4 * w * mu * m * n,
                              (v ** 2 * mv - u ** 2 * mv - w ** 2 * mu) * n ** 2 -
                              2 * w * mu * m * n + 3 * mu * m ** 2,
                              2 * v * w * mu * n ** 2 - 2 * (2 * u - v) * mu * m * n,
                              (2 * u - v) * w * mv * n ** 2 - 3 * v * mv * m * n,
                              (2 * v - u) * w * mv * n ** 2 + 3 * u * mv * m * n,
                              (w ** 2 * mu - u ** 2 * mv - v ** 2 * mv + u * v * mv) *
                              n ** 2 + 3 * mu * m ** 2]
                    m = -1 * m
                    # inverse of the rotation matrix
                    R_list_inv = [(u ** 2 * mv - v ** 2 * mv - w ** 2 * mu) * n ** 2 +
                                  2 * w * mu * m * n + 3 * mu * m ** 2,
                                  (2 * v - u) * u * mv * n ** 2 - 4 * w * mu * m * n,
                                  2 * u * w * mu * n ** 2 + 2 * (2 * v - u) * mu * m * n,
                                  (2 * u - v) * v * mv * n ** 2 + 4 * w * mu * m * n,
                                  (v ** 2 * mv - u ** 2 * mv - w ** 2 * mu) * n ** 2 -
                                  2 * w * mu * m * n + 3 * mu * m ** 2,
                                  2 * v * w * mu * n ** 2 - 2 * (2 * u - v) * mu * m * n,
                                  (2 * u - v) * w * mv * n ** 2 - 3 * v * mv * m * n,
                                  (2 * v - u) * w * mv * n ** 2 + 3 * u * mv * m * n,
                                  (w ** 2 * mu - u ** 2 * mv - v ** 2 * mv + u * v * mv) *
                                  n ** 2 + 3 * mu * m ** 2]
                    m = -1 * m
                    F = 3 * mu * m ** 2 + d * n ** 2
                    all_list = R_list_inv + R_list + [F]
                    # Compute the max common factors for the elements of the rotation matrix
                    # and its inverse.
                    com_fac = reduce(gcd, all_list)
                    sigma = int(round((3 * mu * m ** 2 + d * n ** 2) / com_fac))
                    if (sigma <= cutoff) and (sigma > 1):
                        if sigma not in list(sigmas.keys()):
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n / m * np.sqrt(d / 3.0 / mu)) \
                                        / np.pi * 180
                            sigmas[sigma] = [angle]
                        else:
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n / m * np.sqrt(d / 3.0 / mu)) \
                                        / np.pi * 180
                            if angle not in sigmas[sigma]:
                                sigmas[sigma].append(angle)
            if m_max == 0:
                break
        return sigmas

    @staticmethod
    def enum_sigma_rho(cutoff, r_axis, ratio_alpha):
        """
        Find all possible sigma values and corresponding rotation angles
        within a sigma value cutoff with known rotation axis in rhombohedral system.
        The algorithm for this code is from reference, Acta Cryst, A45,505(1989).

        Args:
            cutoff (integer): the cutoff of sigma values.
            r_axis (list of three integers, e.g. u, v, w
                    or four integers, e.g. u, v, t, w):
                    the rotation axis of the grain boundary, with the format of [u,v,w]
                    or Weber indices [u, v, t, w].
            ratio_alpha (list of two integers, e.g. mu, mv):
                    mu/mv is the ratio of (1+2*cos(alpha))/cos(alpha) with rational number.
                    If irrational, set ratio_alpha = None.
        Returns:
            sigmas (dict):
                    dictionary with keys as the possible integer sigma values
                    and values as list of the possible rotation angles to the
                    corresponding sigma values.
                    e.g. the format as
                    {sigma1: [angle11,angle12,...], sigma2: [angle21, angle22,...],...}
                    Note: the angles are the rotation angle of one grain respect to the
                    other grain.
                    When generate the microstructure of the grain boundary using these
                    angles, you need to analyze the symmetry of the structure. Different
                    angles may result in equivalent microstructures.

        """
        sigmas = {}
        # transform four index notation to three index notation
        if len(r_axis) == 4:
            u1 = r_axis[0]
            v1 = r_axis[1]
            w1 = r_axis[3]
            u = 2 * u1 + v1 + w1
            v = v1 + w1 - u1
            w = w1 - 2 * v1 - u1
            r_axis = [u, v, w]
        # make sure gcd(r_axis)==1
        if reduce(gcd, r_axis) != 1:
            r_axis = [int(round(x / reduce(gcd, r_axis))) for x in r_axis]
        u, v, w = r_axis
        # make sure mu, mv are coprime integers.
        if ratio_alpha is None:
            mu, mv = [1, 1]
            if u + v + w != 0:
                if u != v or u != w:
                    raise RuntimeError('For irrational ratio_alpha, CSL only exist for [1,1,1]'
                                       'or [u, v, -(u+v)] and m =0')
        else:
            mu, mv = ratio_alpha
            if gcd(mu, mv) != 1:
                temp = gcd(mu, mv)
                mu = int(round(mu / temp))
                mv = int(round(mv / temp))

        # refer to the meaning of d in reference
        d = (u ** 2 + v ** 2 + w ** 2) * (mu - 2 * mv) + \
            2 * mv * (v * w + w * u + u * v)
        # Compute the max n we need to enumerate.
        n_max = int(np.sqrt((cutoff * abs(4 * mu * (mu - 3 * mv))) / abs(d)))

        # Enumerate all possible n, m to give possible sigmas within the cutoff.
        for n in range(1, n_max + 1):
            if ratio_alpha is None and u + v + w == 0:
                m_max = 0
            else:
                m_max = int(np.sqrt((cutoff * abs(4 * mu * (mu - 3 * mv)) - n ** 2 * d) / (mu)))
            for m in range(0, m_max + 1):
                if gcd(m, n) == 1 or m == 0:
                    # construct the rotation matrix, refer to the reference
                    R_list = [(mu - 2 * mv) * (u ** 2 - v ** 2 - w ** 2) * n ** 2 +
                              2 * mv * (v - w) * m * n - 2 * mv * v * w * n ** 2 +
                              mu * m ** 2,
                              2 * (mv * u * n * (w * n + u * n - m) - (mu - mv) *
                                   m * w * n + (mu - 2 * mv) * u * v * n ** 2),
                              2 * (mv * u * n * (v * n + u * n + m) + (mu - mv) *
                                   m * v * n + (mu - 2 * mv) * w * u * n ** 2),
                              2 * (mv * v * n * (w * n + v * n + m) + (mu - mv) *
                                   m * w * n + (mu - 2 * mv) * u * v * n ** 2),
                              (mu - 2 * mv) * (v ** 2 - w ** 2 - u ** 2) * n ** 2 +
                              2 * mv * (w - u) * m * n - 2 * mv * u * w * n ** 2 +
                              mu * m ** 2,
                              2 * (mv * v * n * (v * n + u * n - m) - (mu - mv) *
                                   m * u * n + (mu - 2 * mv) * w * v * n ** 2),
                              2 * (mv * w * n * (w * n + v * n - m) - (mu - mv) *
                                   m * v * n + (mu - 2 * mv) * w * u * n ** 2),
                              2 * (mv * w * n * (w * n + u * n + m) + (mu - mv) *
                                   m * u * n + (mu - 2 * mv) * w * v * n ** 2),
                              (mu - 2 * mv) * (w ** 2 - u ** 2 - v ** 2) * n ** 2 +
                              2 * mv * (u - v) * m * n - 2 * mv * u * v * n ** 2 +
                              mu * m ** 2]
                    m = -1 * m
                    # inverse of the rotation matrix
                    R_list_inv = [(mu - 2 * mv) * (u ** 2 - v ** 2 - w ** 2) * n ** 2 +
                                  2 * mv * (v - w) * m * n - 2 * mv * v * w * n ** 2 +
                                  mu * m ** 2,
                                  2 * (mv * u * n * (w * n + u * n - m) - (mu - mv) *
                                       m * w * n + (mu - 2 * mv) * u * v * n ** 2),
                                  2 * (mv * u * n * (v * n + u * n + m) + (mu - mv) *
                                       m * v * n + (mu - 2 * mv) * w * u * n ** 2),
                                  2 * (mv * v * n * (w * n + v * n + m) + (mu - mv) *
                                       m * w * n + (mu - 2 * mv) * u * v * n ** 2),
                                  (mu - 2 * mv) * (v ** 2 - w ** 2 - u ** 2) * n ** 2 +
                                  2 * mv * (w - u) * m * n - 2 * mv * u * w * n ** 2 +
                                  mu * m ** 2,
                                  2 * (mv * v * n * (v * n + u * n - m) - (mu - mv) *
                                       m * u * n + (mu - 2 * mv) * w * v * n ** 2),
                                  2 * (mv * w * n * (w * n + v * n - m) - (mu - mv) *
                                       m * v * n + (mu - 2 * mv) * w * u * n ** 2),
                                  2 * (mv * w * n * (w * n + u * n + m) + (mu - mv) *
                                       m * u * n + (mu - 2 * mv) * w * v * n ** 2),
                                  (mu - 2 * mv) * (w ** 2 - u ** 2 - v ** 2) * n ** 2 +
                                  2 * mv * (u - v) * m * n - 2 * mv * u * v * n ** 2 +
                                  mu * m ** 2]
                    m = -1 * m
                    F = mu * m ** 2 + d * n ** 2
                    all_list = R_list_inv + R_list + [F]
                    # Compute the max common factors for the elements of the rotation matrix
                    #  and its inverse.
                    com_fac = reduce(gcd, all_list)
                    sigma = int(round(abs(F / com_fac)))
                    if (sigma <= cutoff) and (sigma > 1):
                        if sigma not in list(sigmas.keys()):
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n / m * np.sqrt(d / mu)) \
                                        / np.pi * 180
                            sigmas[sigma] = [angle]
                        else:
                            if m == 0:
                                angle = 180
                            else:
                                angle = 2 * np.arctan(n / m * np.sqrt(d / mu)) \
                                        / np.pi * 180.0
                            if angle not in sigmas[sigma]:
                                sigmas[sigma].append(angle)
            if m_max == 0:
                break
        return sigmas

    @staticmethod
    def enum_sigma_tet(cutoff, r_axis, c2_a2_ratio):
        """
        Find all possible sigma values and corresponding rotation angles
        within a sigma value cutoff with known rotation axis in tetragonal system.
        The algorithm for this code is from reference, Acta Cryst, B46,117(1990)

        Args:
            cutoff (integer): the cutoff of sigma values.
            r_axis (list of three integers, e.g. u, v, w):
                    the rotation axis of the grain boundary, with the format of [u,v,w].
            c2_a2_ratio (list of two integers, e.g. mu, mv):
                    mu/mv is the square of the tetragonal axial ratio with rational number.
                    if irrational, set c2_a2_ratio = None
        Returns:
            sigmas (dict):
                    dictionary with keys as the possible integer sigma values
                    and values as list of the possible rotation angles to the
                    corresponding sigma values.
                    e.g. the format as
                    {sigma1: [angle11,angle12,...], sigma2: [angle21, angle22,...],...}
                    Note: the angles are the rotation angle of one grain respect to the
                    other grain.
                    When generate the microstructure of the grain boundary using these
                    angles, you need to analyze the symmetry of the structure. Different
                    angles may result in equivalent microstructures.

        """
        sigmas = {}
        # make sure gcd(r_axis)==1
        if reduce(gcd, r_axis) != 1:
            r_axis = [int(round(x / reduce(gcd, r_axis))) for x in r_axis]

        u, v, w = r_axis

        # make sure mu, mv are coprime integers.
        if c2_a2_ratio is None:
            mu, mv = [1, 1]
            if w != 0:
                if u != 0 or (v != 0):
                    raise RuntimeError('For irrational c2/a2, CSL only exist for [0,0,1] '
                                       'or [u,v,0] and m = 0')
        else:
            mu, mv = c2_a2_ratio
            if gcd(mu, mv) != 1:
                temp = gcd(mu, mv)
                mu = int(round(mu / temp))
                mv = int(round(mv / temp))

        # refer to the meaning of d in reference
        d = (u ** 2 + v ** 2) * mv + w ** 2 * mu

        # Compute the max n we need to enumerate.
        n_max = int(np.sqrt((cutoff * 4 * mu * mv) / d))

        # Enumerate all possible n, m to give possible sigmas within the cutoff.
        for n in range(1, n_max + 1):
            if c2_a2_ratio is None and w == 0:
                m_max = 0
            else:
                m_max = int(np.sqrt((cutoff * 4 * mu * mv - n ** 2 * d) / mu))
            for m in range(0, m_max + 1):
                if gcd(m, n) == 1 or m == 0:
                    # construct the rotation matrix, refer to the reference
                    R_list = [(u ** 2 * mv - v ** 2 * mv - w ** 2 * mu) * n ** 2 +
                              mu * m ** 2,
                              2 * v * u * mv * n ** 2 - 2 * w * mu * m * n,
                              2 * u * w * mu * n ** 2 + 2 * v * mu * m * n,
                              2 * u * v * mv * n ** 2 + 2 * w * mu * m * n,
                              (v ** 2 * mv - u ** 2 * mv - w ** 2 * mu) * n ** 2 +
                              mu * m ** 2,
                              2 * v * w * mu * n ** 2 - 2 * u * mu * m * n,
                              2 * u * w * mv * n ** 2 - 2 * v * mv * m * n,
                              2 * v * w * mv * n ** 2 + 2 * u * mv * m * n,
                              (w ** 2 * mu - u ** 2 * mv - v ** 2 * mv) * n ** 2 +
                              mu * m ** 2]
                    m = -1 * m
                    # inverse of rotation matrix
                    R_list_inv = [(u ** 2 * mv - v ** 2 * mv - w ** 2 * mu) * n ** 2 +
                                  mu * m ** 2,
                                  2 * v * u * mv * n ** 2 - 2 * w * mu * m * n,
                                  2 * u * w * mu * n ** 2 + 2 * v * mu * m * n,
                                  2 * u * v * mv * n ** 2 + 2 * w * mu * m * n,
                                  (v ** 2 * mv - u ** 2 * mv - w ** 2 * mu) * n ** 2 +
                                  mu * m ** 2,
                                  2 * v * w * mu * n ** 2 - 2 * u * mu * m * n,
                                  2 * u * w * mv * n ** 2 - 2 * v * mv * m * n,
                                  2 * v * w * mv * n ** 2 + 2 * u * mv * m * n,
                                  (w ** 2 * mu - u ** 2 * mv - v ** 2 * mv) * n ** 2 +
                                  mu * m ** 2]
                    m = -1 * m
                    F = mu * m ** 2 + d * n ** 2
                    all_list = R_list + R_list_inv + [F]
                    # Compute the max common factors for the elements of the rotation matrix
                    #  and its inverse.
                    com_fac = reduce(gcd, all_list)
                    sigma = int(round((mu * m ** 2 + d * n ** 2) / com_fac))
                    if (sigma <= cutoff) and (sigma > 1):
                        if sigma not in list(sigmas.keys()):
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n / m * np.sqrt(d / mu)) \
                                        / np.pi * 180
                            sigmas[sigma] = [angle]
                        else:
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n / m * np.sqrt(d / mu)) \
                                        / np.pi * 180
                            if angle not in sigmas[sigma]:
                                sigmas[sigma].append(angle)
            if m_max == 0:
                break

        return sigmas

    @staticmethod
    def enum_sigma_ort(cutoff, r_axis, c2_b2_a2_ratio):
        """
        Find all possible sigma values and corresponding rotation angles
        within a sigma value cutoff with known rotation axis in orthorhombic system.
        The algorithm for this code is from reference, Scipta Metallurgica 27, 291(1992)

        Args:
            cutoff (integer): the cutoff of sigma values.
            r_axis (list of three integers, e.g. u, v, w):
                    the rotation axis of the grain boundary, with the format of [u,v,w].
            c2_b2_a2_ratio (list of three integers, e.g. mu,lamda, mv):
                    mu:lam:mv is the square of the orthorhombic axial ratio with rational
                    numbers. If irrational for one axis, set it to None.
                    e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
        Returns:
            sigmas (dict):
                    dictionary with keys as the possible integer sigma values
                    and values as list of the possible rotation angles to the
                    corresponding sigma values.
                    e.g. the format as
                    {sigma1: [angle11,angle12,...], sigma2: [angle21, angle22,...],...}
                    Note: the angles are the rotation angle of one grain respect to the
                    other grain.
                    When generate the microstructure of the grain boundary using these
                    angles, you need to analyze the symmetry of the structure. Different
                    angles may result in equivalent microstructures.

        """
        sigmas = {}
        # make sure gcd(r_axis)==1
        if reduce(gcd, r_axis) != 1:
            r_axis = [int(round(x / reduce(gcd, r_axis))) for x in r_axis]

        u, v, w = r_axis
        # make sure mu, lambda, mv are coprime integers.
        if None in c2_b2_a2_ratio:
            mu, lam, mv = c2_b2_a2_ratio
            non_none = [i for i in c2_b2_a2_ratio if i is not None]
            if len(non_none) < 2:
                raise RuntimeError('No CSL exist for two irrational numbers')
            non1, non2 = non_none
            if reduce(gcd, non_none) != 1:
                temp = reduce(gcd, non_none)
                non1 = int(round(non1 / temp))
                non2 = int(round(non2 / temp))
            if mu is None:
                lam = non1
                mv = non2
                mu = 1
                if w != 0:
                    if u != 0 or (v != 0):
                        raise RuntimeError('For irrational c2, CSL only exist for [0,0,1] '
                                           'or [u,v,0] and m = 0')
            elif lam is None:
                mu = non1
                mv = non2
                lam = 1
                if v != 0:
                    if u != 0 or (w != 0):
                        raise RuntimeError('For irrational b2, CSL only exist for [0,1,0] '
                                           'or [u,0,w] and m = 0')
            elif mv is None:
                mu = non1
                lam = non2
                mv = 1
                if u != 0:
                    if w != 0 or (v != 0):
                        raise RuntimeError('For irrational a2, CSL only exist for [1,0,0] '
                                           'or [0,v,w] and m = 0')
        else:
            mu, lam, mv = c2_b2_a2_ratio
            if reduce(gcd, c2_b2_a2_ratio) != 1:
                temp = reduce(gcd, c2_b2_a2_ratio)
                mu = int(round(mu / temp))
                mv = int(round(mv / temp))
                lam = int(round(lam / temp))
            if u == 0 and v == 0:
                mu = 1
            if u == 0 and w == 0:
                lam = 1
            if v == 0 and w == 0:
                mv = 1
        # refer to the meaning of d in reference
        d = (mv * u ** 2 + lam * v ** 2) * mv + w ** 2 * mu * mv

        # Compute the max n we need to enumerate.
        n_max = int(np.sqrt((cutoff * 4 * mu * mv * mv * lam) / d))
        # Enumerate all possible n, m to give possible sigmas within the cutoff.
        for n in range(1, n_max + 1):
            mu_temp, lam_temp, mv_temp = c2_b2_a2_ratio
            if (mu_temp is None and w == 0) or (lam_temp is None and v == 0) \
                    or (mv_temp is None and u == 0):
                m_max = 0
            else:
                m_max = int(np.sqrt((cutoff * 4 * mu * mv * lam * mv -
                                     n ** 2 * d) / mu / lam))
            for m in range(0, m_max + 1):

                if gcd(m, n) == 1 or m == 0:
                    # construct the rotation matrix, refer to the reference
                    R_list = [(u ** 2 * mv * mv - lam * v ** 2 * mv -
                               w ** 2 * mu * mv) * n ** 2 + lam * mu * m ** 2,
                              2 * lam * (v * u * mv * n ** 2 - w * mu * m * n),
                              2 * mu * (u * w * mv * n ** 2 + v * lam * m * n),
                              2 * mv * (u * v * mv * n ** 2 + w * mu * m * n),
                              (v ** 2 * mv * lam - u ** 2 * mv * mv -
                               w ** 2 * mu * mv) * n ** 2 + lam * mu * m ** 2,
                              2 * mv * mu * (v * w * n ** 2 - u * m * n),
                              2 * mv * (u * w * mv * n ** 2 - v * lam * m * n),
                              2 * lam * mv * (v * w * n ** 2 + u * m * n),
                              (w ** 2 * mu * mv - u ** 2 * mv * mv -
                               v ** 2 * mv * lam) * n ** 2 + lam * mu * m ** 2]
                    m = -1 * m
                    # inverse of rotation matrix
                    R_list_inv = [(u ** 2 * mv * mv - lam * v ** 2 * mv -
                                   w ** 2 * mu * mv) * n ** 2 + lam * mu * m ** 2,
                                  2 * lam * (v * u * mv * n ** 2 - w * mu * m * n),
                                  2 * mu * (u * w * mv * n ** 2 + v * lam * m * n),
                                  2 * mv * (u * v * mv * n ** 2 + w * mu * m * n),
                                  (v ** 2 * mv * lam - u ** 2 * mv * mv -
                                   w ** 2 * mu * mv) * n ** 2 + lam * mu * m ** 2,
                                  2 * mv * mu * (v * w * n ** 2 - u * m * n),
                                  2 * mv * (u * w * mv * n ** 2 - v * lam * m * n),
                                  2 * lam * mv * (v * w * n ** 2 + u * m * n),
                                  (w ** 2 * mu * mv - u ** 2 * mv * mv -
                                   v ** 2 * mv * lam) * n ** 2 + lam * mu * m ** 2]
                    m = -1 * m
                    F = mu * lam * m ** 2 + d * n ** 2
                    all_list = R_list + R_list_inv + [F]
                    # Compute the max common factors for the elements of the rotation matrix
                    #  and its inverse.
                    com_fac = reduce(gcd, all_list)
                    sigma = int(round((mu * lam * m ** 2 + d * n ** 2) / com_fac))
                    if (sigma <= cutoff) and (sigma > 1):
                        if sigma not in list(sigmas.keys()):
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n / m * np.sqrt(d / mu / lam)) \
                                        / np.pi * 180
                            sigmas[sigma] = [angle]
                        else:
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n / m * np.sqrt(d / mu / lam)) \
                                        / np.pi * 180
                            if angle not in sigmas[sigma]:
                                sigmas[sigma].append(angle)
            if m_max == 0:
                break

        return sigmas

    @staticmethod
    def enum_possible_plane_cubic(plane_cutoff, r_axis, r_angle):
        """
        Find all possible plane combinations for GBs given a rotaion axis and angle for
        cubic system, and classify them to different categories, including 'Twist',
        'Symmetric tilt', 'Normal tilt', 'Mixed' GBs.

        Args:
            plane_cutoff (integer): the cutoff of plane miller index.
            r_axis (list of three integers, e.g. u, v, w):
                    the rotation axis of the grain boundary, with the format of [u,v,w].
            r_angle (float): rotation angle of the GBs.
        Returns:
            all_combinations (dict):
                    dictionary with keys as GB type, e.g. 'Twist','Symmetric tilt',etc.
                    and values as the combination of the two plane miller index
                     (GB plane and joining plane).
        """
        all_combinations = {}
        all_combinations['Symmetric tilt'] = []
        all_combinations['Twist'] = []
        all_combinations['Normal tilt'] = []
        all_combinations['Mixed'] = []
        sym_plane = symm_group_cubic([[1, 0, 0], [1, 1, 0]])
        j = np.arange(0, plane_cutoff + 1)
        combination = []
        for i in itertools.product(j, repeat=3):
            if sum(abs(np.array(i))) != 0:
                combination.append(list(i))
            if len(np.nonzero(i)[0]) == 3:
                for i1 in range(3):
                    new_i = list(i).copy()
                    new_i[i1] = -1 * new_i[i1]
                    combination.append(new_i)
            elif len(np.nonzero(i)[0]) == 2:
                new_i = list(i).copy()
                new_i[np.nonzero(i)[0][0]] = -1 * new_i[np.nonzero(i)[0][0]]
                combination.append(new_i)
        miller = np.array(combination)
        miller = miller[np.argsort(np.linalg.norm(miller, axis=1))]
        for i, val in enumerate(miller):
            if reduce(gcd, val) == 1:
                matrix = GrainBoundaryGenerator.get_trans_mat(r_axis, r_angle, surface=val, quick_gen=True)
                vec = np.cross(matrix[1][0], matrix[1][1])
                miller2 = GrainBoundaryGenerator.vec_to_surface(vec)
                if np.all(np.abs(np.array(miller2)) <= plane_cutoff):
                    cos_1 = abs(np.dot(val, r_axis) / np.linalg.norm(val) / np.linalg.norm(r_axis))
                    if 1 - cos_1 < 1.e-5:
                        all_combinations['Twist'].append([list(val), miller2])
                    elif cos_1 < 1.e-8:
                        sym_tilt = False
                        if np.sum(np.abs(val)) == np.sum(np.abs(miller2)):
                            ave = (np.array(val) + np.array(miller2)) / 2
                            ave1 = (np.array(val) - np.array(miller2)) / 2
                            for plane in sym_plane:
                                cos_2 = abs(np.dot(ave, plane) / np.linalg.norm(ave) / np.linalg.norm(plane))
                                cos_3 = abs(np.dot(ave1, plane) / np.linalg.norm(ave1) / np.linalg.norm(plane))
                                if 1 - cos_2 < 1.e-5 or 1 - cos_3 < 1.e-5:
                                    all_combinations['Symmetric tilt'].append([list(val), miller2])
                                    sym_tilt = True
                                    break
                        if not sym_tilt:
                            all_combinations['Normal tilt'].append([list(val), miller2])
                    else:
                        all_combinations['Mixed'].append([list(val), miller2])
        return all_combinations

    @staticmethod
    def get_rotation_angle_from_sigma(sigma, r_axis, lat_type='C', ratio=None):
        """
        Find all possible rotation angle for the given sigma value.

        Args:
            sigma (integer):
                    sigma value provided
            r_axis (list of three integers, e.g. u, v, w
                    or four integers, e.g. u, v, t, w for hex/rho system only):
                    the rotation axis of the grain boundary.
            lat_type ( one character):
                    'c' or 'C': cubic system
                     't' or 'T': tetragonal system
                     'o' or 'O': orthorhombic system
                     'h' or 'H': hexagonal system
                     'r' or 'R': rhombohedral system
                     default to cubic system
            ratio (list of integers):
                    lattice axial ratio.
                    For cubic system, ratio is not needed.
                    For tetragonal system, ratio = [mu, mv], list of two integers,
                    that is, mu/mv = c2/a2. If it is irrational, set it to none.
                    For orthorhombic system, ratio = [mu, lam, mv], list of three integers,
                    that is, mu:lam:mv = c2:b2:a2. If irrational for one axis, set it to None.
                    e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
                    For rhombohedral system, ratio = [mu, mv], list of two integers,
                    that is, mu/mv is the ratio of (1+2*cos(alpha)/cos(alpha).
                    If irrational, set it to None.
                    For hexagonal system, ratio = [mu, mv], list of two integers,
                    that is, mu/mv = c2/a2. If it is irrational, set it to none.

        Returns:
            rotation_angles corresponding to the provided sigma value.
            If the sigma value is not correct, return the rotation angle corresponding
            to the correct possible sigma value right smaller than the wrong sigma value provided.
        """
        if lat_type.lower() == 'c':
            logger.info('Make sure this is for cubic system')
            sigma_dict = GrainBoundaryGenerator.enum_sigma_cubic(cutoff=sigma, r_axis=r_axis)
        elif lat_type.lower() == 't':
            logger.info('Make sure this is for tetragonal system')
            if ratio is None:
                logger.info('Make sure this is for irrational c2/a2 ratio')
            elif len(ratio) != 2:
                raise RuntimeError('Tetragonal system needs correct c2/a2 ratio')
            sigma_dict = GrainBoundaryGenerator.enum_sigma_tet(cutoff=sigma, r_axis=r_axis, c2_a2_ratio=ratio)
        elif lat_type.lower() == 'o':
            logger.info('Make sure this is for orthorhombic system')
            if len(ratio) != 3:
                raise RuntimeError('Orthorhombic system needs correct c2:b2:a2 ratio')
            sigma_dict = GrainBoundaryGenerator.enum_sigma_ort(cutoff=sigma, r_axis=r_axis, c2_b2_a2_ratio=ratio)
        elif lat_type.lower() == 'h':
            logger.info('Make sure this is for hexagonal system')
            if ratio is None:
                logger.info('Make sure this is for irrational c2/a2 ratio')
            elif len(ratio) != 2:
                raise RuntimeError('Hexagonal system needs correct c2/a2 ratio')
            sigma_dict = GrainBoundaryGenerator.enum_sigma_hex(cutoff=sigma, r_axis=r_axis, c2_a2_ratio=ratio)
        elif lat_type.lower() == 'r':
            logger.info('Make sure this is for rhombohedral system')
            if ratio is None:
                logger.info('Make sure this is for irrational (1+2*cos(alpha)/cos(alpha) ratio')
            elif len(ratio) != 2:
                raise RuntimeError('Rhombohedral system needs correct '
                                   '(1+2*cos(alpha)/cos(alpha) ratio')
            sigma_dict = GrainBoundaryGenerator.enum_sigma_rho(cutoff=sigma, r_axis=r_axis, ratio_alpha=ratio)
        else:
            raise RuntimeError('Lattice type not implemented')

        sigmas = list(sigma_dict.keys())
        if not sigmas:
            raise RuntimeError('This is a wriong sigma value, and no sigma exists smaller than this value.')
        if sigma in sigmas:
            rotation_angles = sigma_dict[sigma]
        else:
            sigmas.sort()
            warnings.warn("This is not the possible sigma value according to the rotation axis!"
                          "The nearest neighbor sigma and its corresponding angle are returned")
            rotation_angles = sigma_dict[sigmas[-1]]
        rotation_angles.sort()
        return rotation_angles

    @staticmethod
    def slab_from_csl(csl, surface, normal, trans_cry, max_search=20, quick_gen=False):
        """
        By linear operation of csl lattice vectors to get the best corresponding
        slab lattice. That is the area of a,b vectors (within the surface plane)
        is the smallest, the c vector first, has shortest length perpendicular
        to surface [h,k,l], second, has shortest length itself.

        Args:
            csl (3 by 3 integer array):
                    input csl lattice.
            surface (list of three integers, e.g. h, k, l):
                    the miller index of the surface, with the format of [h,k,l]
            normal (logic):
                    determine if the c vector needs to perpendicular to surface
            trans_cry (3 by 3 array):
                    transform matrix from crystal system to orthogonal system
            max_search (int): max search for the GB lattice vectors that give the smallest GB
                lattice. If normal is true, also max search the GB c vector that perpendicular
                to the plane.
            quick_gen (bool): whether to quickly generate a supercell, no need to find the smallest
                cell if set to true.

        Returns:
            t_matrix: a slab lattice ( 3 by 3 integer array):
        """

        # set the transform matrix in real space
        trans = trans_cry
        # transform matrix in reciprocal space
        ctrans = np.linalg.inv(trans.T)

        t_matrix = csl.copy()
        # vectors constructed from csl that perpendicular to surface
        ab_vector = []
        # obtain the miller index of surface in terms of csl.
        miller = np.matmul(surface, csl.T)
        if reduce(gcd, miller) != 1:
            miller = [int(round(x / reduce(gcd, miller))) for x in miller]
        miller_nonzero = []
        # quickly generate a supercell, normal is not work in this way
        if quick_gen:
            scale_factor = []
            eye = np.eye(3, dtype=np.int)
            for i, j in enumerate(miller):
                if j == 0:
                    scale_factor.append(eye[i])
                else:
                    miller_nonzero.append(i)
            if len(scale_factor) < 2:
                index_len = len(miller_nonzero)
                for i in range(index_len):
                    for j in range(i + 1, index_len):
                        lcm_miller = lcm(miller[miller_nonzero[i]], miller[miller_nonzero[j]])
                        l = [0, 0, 0]
                        l[miller_nonzero[i]] = -int(round(lcm_miller / miller[miller_nonzero[i]]))
                        l[miller_nonzero[j]] = int(round(lcm_miller / miller[miller_nonzero[j]]))
                        scale_factor.append(l)
                        if len(scale_factor) == 2:
                            break
            t_matrix[0] = np.array(np.dot(scale_factor[0], csl))
            t_matrix[1] = np.array(np.dot(scale_factor[1], csl))
            t_matrix[2] = csl[miller_nonzero[0]]
            if abs(np.linalg.det(t_matrix)) > 1000:
                warnings.warn('Too large matrix. Suggest to use quick_gen=False')
            return t_matrix

        for i, j in enumerate(miller):
            if j == 0:
                ab_vector.append(csl[i])
            else:
                c_index = i
                miller_nonzero.append(j)

        if len(miller_nonzero) > 1:
            t_matrix[2] = csl[c_index]
            index_len = len(miller_nonzero)
            lcm_miller = []
            for i in range(index_len):
                for j in range(i + 1, index_len):
                    com_gcd = gcd(miller_nonzero[i], miller_nonzero[j])
                    mil1 = int(round(miller_nonzero[i] / com_gcd))
                    mil2 = int(round(miller_nonzero[j] / com_gcd))
                    lcm_miller.append(max(abs(mil1), abs(mil2)))
            lcm_sorted = sorted(lcm_miller)
            if index_len == 2:
                max_j = lcm_sorted[0]
            else:
                max_j = lcm_sorted[1]
        else:
            if not normal:
                t_matrix[0] = ab_vector[0]
                t_matrix[1] = ab_vector[1]
                t_matrix[2] = csl[c_index]
                return t_matrix
            else:
                max_j = abs(miller_nonzero[0])
        if max_j > max_search:
            max_j = max_search
        # area of a, b vectors
        area = None
        # length of c vector
        c_norm = np.linalg.norm(np.matmul(t_matrix[2], trans))
        # c vector length along the direction perpendicular to surface
        c_length = np.abs(np.dot(t_matrix[2], surface))
        # check if the init c vector perpendicular to the surface
        if normal:
            c_cross = np.cross(np.matmul(t_matrix[2], trans), np.matmul(surface, ctrans))
            if np.linalg.norm(c_cross) < 1.e-8:
                normal_init = True
            else:
                normal_init = False

        j = np.arange(0, max_j + 1)
        combination = []
        for i in itertools.product(j, repeat=3):
            if sum(abs(np.array(i))) != 0:
                combination.append(list(i))
            if len(np.nonzero(i)[0]) == 3:
                for i1 in range(3):
                    new_i = list(i).copy()
                    new_i[i1] = -1 * new_i[i1]
                    combination.append(new_i)
            elif len(np.nonzero(i)[0]) == 2:
                new_i = list(i).copy()
                new_i[np.nonzero(i)[0][0]] = -1 * new_i[np.nonzero(i)[0][0]]
                combination.append(new_i)
        for i in combination:
            if reduce(gcd, i) == 1:
                temp = np.dot(np.array(i), csl)
                if abs(np.dot(temp, surface) - 0) < 1.e-8:
                    ab_vector.append(temp)
                else:
                    # c vector length along the direction perpendicular to surface
                    c_len_temp = np.abs(np.dot(temp, surface))
                    # c vector length itself
                    c_norm_temp = np.linalg.norm(np.matmul(temp, trans))
                    if normal:
                        c_cross = np.cross(np.matmul(temp, trans), np.matmul(surface, ctrans))
                        if np.linalg.norm(c_cross) < 1.e-8:
                            if normal_init:
                                if c_norm_temp < c_norm:
                                    t_matrix[2] = temp
                                    c_norm = c_norm_temp
                            else:
                                c_norm = c_norm_temp
                                normal_init = True
                                t_matrix[2] = temp
                    else:
                        if c_len_temp < c_length or \
                                (abs(c_len_temp - c_length) < 1.e-8 and c_norm_temp < c_norm):
                            t_matrix[2] = temp
                            c_norm = c_norm_temp
                            c_length = c_len_temp

        if normal and (not normal_init):
            logger.info('Did not find the perpendicular c vector, increase max_j')
            while (not normal_init):
                if max_j == max_search:
                    warnings.warn('Cannot find the perpendicular c vector, please increase max_search')
                    break
                max_j = 3 * max_j
                if max_j > max_search:
                    max_j = max_search
                j = np.arange(0, max_j + 1)
                combination = []
                for i in itertools.product(j, repeat=3):
                    if sum(abs(np.array(i))) != 0:
                        combination.append(list(i))
                    if len(np.nonzero(i)[0]) == 3:
                        for i1 in range(3):
                            new_i = list(i).copy()
                            new_i[i1] = -1 * new_i[i1]
                            combination.append(new_i)
                    elif len(np.nonzero(i)[0]) == 2:
                        new_i = list(i).copy()
                        new_i[np.nonzero(i)[0][0]] = -1 * new_i[np.nonzero(i)[0][0]]
                        combination.append(new_i)
                for i in combination:
                    if reduce(gcd, i) == 1:
                        temp = np.dot(np.array(i), csl)
                        if abs(np.dot(temp, surface) - 0) > 1.e-8:
                            c_cross = np.cross(np.matmul(temp, trans), np.matmul(surface, ctrans))
                            if np.linalg.norm(c_cross) < 1.e-8:
                                # c vetor length itself
                                c_norm_temp = np.linalg.norm(np.matmul(temp, trans))
                                if normal_init:
                                    if c_norm_temp < c_norm:
                                        t_matrix[2] = temp
                                        c_norm = c_norm_temp
                                else:
                                    c_norm = c_norm_temp
                                    normal_init = True
                                    t_matrix[2] = temp
                if normal_init:
                    logger.info('Found perpendicular c vector')

        # find the best a, b vectors with their formed area smallest and average norm of a,b smallest.
        for i in itertools.combinations(ab_vector, 2):
            area_temp = np.linalg.norm(np.cross(np.matmul(i[0], trans),
                                                np.matmul(i[1], trans)))
            if abs(area_temp - 0) > 1.e-8:
                ab_norm_temp = np.linalg.norm(np.matmul(i[0], trans)) + \
                               np.linalg.norm(np.matmul(i[1], trans))
                if area is None:
                    area = area_temp
                    ab_norm = ab_norm_temp
                    t_matrix[0] = i[0]
                    t_matrix[1] = i[1]
                elif area_temp < area:
                    t_matrix[0] = i[0]
                    t_matrix[1] = i[1]
                    area = area_temp
                    ab_norm = ab_norm_temp
                elif abs(area - area_temp) < 1.e-8 and ab_norm_temp < ab_norm:
                    t_matrix[0] = i[0]
                    t_matrix[1] = i[1]
                    area = area_temp
                    ab_norm = ab_norm_temp

        # make sure we have a left-handed crystallographic system
        if np.linalg.det(np.matmul(t_matrix, trans)) < 0:
            t_matrix *= -1

        if normal and abs(np.linalg.det(t_matrix)) > 1000:
            warnings.warn('Too large matrix. Suggest to use Normal=False')
        return t_matrix

    @staticmethod
    def reduce_mat(mat, mag, r_matrix):
        """
        Reduce integer array mat's determinant mag times by linear combination
        of its row vectors, so that the new array after rotation (r_matrix) is
        still an integer array

        Args:
            mat (3 by 3 array): input matrix
            mag (integer): reduce times for the determinant
            r_matrix (3 by 3 array): rotation matrix
        Return:
            the reduced integer array
        """
        max_j = abs(int(round(np.linalg.det(mat) / mag)))
        reduced = False
        for h in range(3):
            k = h + 1 if h + 1 < 3 else abs(2 - h)
            l = h + 2 if h + 2 < 3 else abs(1 - h)
            j = np.arange(-max_j, max_j + 1)
            for j1, j2 in itertools.product(j, repeat=2):
                temp = mat[h] + j1 * mat[k] + j2 * mat[l]
                if all([np.round(x, 5).is_integer() for x in list(temp / mag)]):
                    mat_copy = mat.copy()
                    mat_copy[h] = np.array([int(round(ele / mag)) for ele in temp])
                    new_mat = np.dot(mat_copy, np.linalg.inv(r_matrix.T))
                    if all([np.round(x, 5).is_integer() for x in list(np.ravel(new_mat))]):
                        reduced = True
                        mat[h] = np.array([int(round(ele / mag)) for ele in temp])
                        break
            if reduced:
                break

        if not reduced:
            warnings.warn("Matrix reduction not performed, may lead to non-primitive gb cell.")
        return mat

    @staticmethod
    def vec_to_surface(vec):
        """
        Transform a float vector to a surface miller index with integers.

        Args:
            vec (1 by 3 array float vector): input float vector
        Return:
            the surface miller index of the input vector.
        """
        miller = [None] * 3
        index = []
        for i, value in enumerate(vec):
            if abs(value) < 1.e-8:
                miller[i] = 0
            else:
                index.append(i)
        if len(index) == 1:
            miller[index[0]] = 1
        else:
            min_index = np.argmin([i for i in vec if i != 0])
            true_index = index[min_index]
            index.pop(min_index)
            frac = []
            for i, value in enumerate(index):
                frac.append(Fraction(vec[value] / vec[true_index]).limit_denominator(100))
            if len(index) == 1:
                miller[true_index] = frac[0].denominator
                miller[index[0]] = frac[0].numerator
            else:
                com_lcm = lcm(frac[0].denominator, frac[1].denominator)
                miller[true_index] = com_lcm
                miller[index[0]] = frac[0].numerator * int(round((com_lcm / frac[0].denominator)))
                miller[index[1]] = frac[1].numerator * int(round((com_lcm / frac[1].denominator)))
        return miller


def factors(n):
    """
    Compute the factors of a integer.
    Args:
        n: the input integer

    Returns:
        a set of integers that are the factors of the input integer.
    """
    return set(reduce(list.__add__,
                      ([i, n // i] for i in range(1, int(np.sqrt(n)) + 1) if n % i == 0)))


def fix_pbc(structure, matrix=None):
    """
    Set all frac_coords of the input structure within [0,1].

    Args:
        structure (pymatgen structure object):
            input structure
        matrix (lattice matrix, 3 by 3 array/matrix)
            new structure's lattice matrix, if none, use
            input structure's matrix

    Return:
        new structure with fixed frac_coords and lattice matrix
    """

    spec = []
    coords = []
    if matrix is None:
        latte = Lattice(structure.lattice.matrix)
    else:
        latte = Lattice(matrix)

    for site in structure:
        spec.append(site.specie)
        coord = np.array(site.frac_coords)
        for i in range(3):
            coord[i] -= floor(coord[i])
            if np.allclose(coord[i], 1):
                coord[i] = 0
            elif np.allclose(coord[i], 0):
                coord[i] = 0
            else:
                coord[i] = round(coord[i], 7)
        coords.append(coord)

    return Structure(latte, spec, coords, site_properties=structure.site_properties)


def symm_group_cubic(mat):
    """
     obtain cubic symmetric eqivalents of the list of vectors.

    Args:
        matrix (lattice matrix, n by 3 array/matrix)

    Return:
        cubic symmetric eqivalents of the list of vectors.
    """
    sym_group = np.zeros([24, 3, 3])
    sym_group[0, :] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    sym_group[1, :] = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    sym_group[2, :] = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
    sym_group[3, :] = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
    sym_group[4, :] = [[0, -1, 0], [-1, 0, 0], [0, 0, -1]]
    sym_group[5, :] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    sym_group[6, :] = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
    sym_group[7, :] = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
    sym_group[8, :] = [[-1, 0, 0], [0, 0, -1], [0, -1, 0]]
    sym_group[9, :] = [[-1, 0, 0], [0, 0, 1], [0, 1, 0]]
    sym_group[10, :] = [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
    sym_group[11, :] = [[1, 0, 0], [0, 0, 1], [0, -1, 0]]
    sym_group[12, :] = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    sym_group[13, :] = [[0, 1, 0], [0, 0, -1], [-1, 0, 0]]
    sym_group[14, :] = [[0, -1, 0], [0, 0, 1], [-1, 0, 0]]
    sym_group[15, :] = [[0, -1, 0], [0, 0, -1], [1, 0, 0]]
    sym_group[16, :] = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
    sym_group[17, :] = [[0, 0, 1], [-1, 0, 0], [0, -1, 0]]
    sym_group[18, :] = [[0, 0, -1], [1, 0, 0], [0, -1, 0]]
    sym_group[19, :] = [[0, 0, -1], [-1, 0, 0], [0, 1, 0]]
    sym_group[20, :] = [[0, 0, -1], [0, -1, 0], [-1, 0, 0]]
    sym_group[21, :] = [[0, 0, -1], [0, 1, 0], [1, 0, 0]]
    sym_group[22, :] = [[0, 0, 1], [0, -1, 0], [1, 0, 0]]
    sym_group[23, :] = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]

    mat = np.atleast_2d(mat)
    all_vectors = []
    for sym in sym_group:
        for vec in mat:
            all_vectors.append(np.dot(sym, vec))
    return np.unique(np.array(all_vectors), axis=0)
