"""This module provides classes to store, generate,
and manipulate interfaces, including grain boundaries.
"""

from __future__ import annotations

import logging
import math
import warnings
from fractions import Fraction
from functools import reduce
from itertools import chain, combinations, product
from typing import TYPE_CHECKING, Literal, cast

import numpy as np
from monty.fractions import lcm
from numpy.testing import assert_allclose
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite, Site
from pymatgen.core.structure import Structure
from pymatgen.core.surface import Slab
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.typing import Tuple3Ints

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from typing import Any

    from numpy.typing import ArrayLike, NDArray
    from typing_extensions import Self

    from pymatgen.core import Element
    from pymatgen.util.typing import CompositionLike, Matrix3D, MillerIndex, Tuple3Floats, Vector3D

Tuple4Ints = tuple[int, int, int, int]
logger = logging.getLogger(__name__)


__author__ = "Xiang-Guo Li"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Xiang-Guo Li"
__email__ = "xil110@ucsd.edu"
__date__ = "7/30/18"


class GrainBoundary(Structure):
    """
    Representation of grain boundary (GB). Implements additional
    attributes pertaining to GBs, but the init method does not actually implement any
    algorithm that creates a GB. This is a DUMMY class who's init method only holds
    information about the GB. Also has additional methods that returns other information
    about a GB such as sigma value.

    Note that all GBs have their surface normal oriented in the c-direction. This means
    the lattice vectors a and b are in the GB surface plane (at least for one grain) and
    the c vector is out of the surface plane (though not necessarily perpendicular to the
    surface).
    """

    def __init__(
        self,
        lattice: np.ndarray | Lattice,
        species: Sequence[CompositionLike],
        coords: Sequence[ArrayLike],
        rotation_axis: Tuple3Ints | Tuple4Ints,
        rotation_angle: float,
        gb_plane: Tuple3Ints,
        join_plane: Tuple3Ints,
        init_cell: Structure,
        vacuum_thickness: float,
        ab_shift: tuple[float, float],
        site_properties: dict[str, Any],
        oriented_unit_cell: Structure,
        validate_proximity: bool = False,
        coords_are_cartesian: bool = False,
        properties: dict | None = None,
    ) -> None:
        """A Structure with additional information and methods pertaining to GBs.

        Args:
            lattice (Lattice | np.ndarray): The lattice, either as an instance or
                a 3x3 array. Each row should correspond to a lattice vector.
            species ([Species]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / species specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g. (3, 56, ...) or actual Element or Species objects.

                ii. List of dict of elements/species and occupancies, e.g.
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates for each species.
            rotation_axis (list[int]): Rotation axis of GB in the form of a list of integers, e.g. [1, 1, 0].
            rotation_angle (float, in unit of degree): rotation angle of GB.
            gb_plane (list): Grain boundary plane in the form of a list of integers
                e.g.: [1, 2, 3].
            join_plane (list): Joining plane of the second grain in the form of a list of
                integers. e.g.: [1, 2, 3].
            init_cell (Structure): initial bulk structure to form the GB.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, The sequences have to be the same length as
                the atomic species and fractional_coords. For GB, you should
                have the 'grain_label' properties to classify the sites as 'top',
                'bottom', 'top_incident', or 'bottom_incident'.
            vacuum_thickness (float in angstrom): The thickness of vacuum inserted
                between two grains of the GB.
            ab_shift (list of float, in unit of crystal vector a, b): The relative
                shift along a, b vectors.
            oriented_unit_cell (Structure): oriented unit cell of the bulk init_cell.
                Helps to accurately calculate the bulk properties that are consistent
                with GB calculations.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in Cartesian coordinates. Defaults to False.
            properties (dict): dictionary containing properties associated
                with the whole GrainBoundary.
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
            lattice,
            species,
            coords,
            validate_proximity=validate_proximity,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties,
            properties=properties,
        )

    def __str__(self) -> str:
        comp = self.composition
        outs = [
            f"Gb Summary ({comp.formula})",
            f"Reduced Formula: {comp.reduced_formula}",
            f"Rotation axis: {self.rotation_axis}",
            f"Rotation angle: {self.rotation_angle}",
            f"GB plane: {self.gb_plane}",
            f"Join plane: {self.join_plane}",
            f"vacuum thickness: {self.vacuum_thickness}",
            f"ab_shift: {self.ab_shift}",
        ]

        def to_str(number: float, rjust: int = 10) -> str:
            """Convert a float to string and right justify."""
            return (f"{number:0.6f}").rjust(rjust)

        outs += (
            f"abc   : {' '.join(to_str(i) for i in self.lattice.abc)}",
            f"angles: {' '.join(to_str(i) for i in self.lattice.angles)}",
            f"Sites ({len(self)})",
        )
        for idx, site in enumerate(self, start=1):
            outs.append(f"{idx} {site.species_string} {' '.join(to_str(coord, 12) for coord in site.frac_coords)}")
        return "\n".join(outs)

    def copy(self) -> Self:
        """Make a copy of the GrainBoundary."""
        return type(self)(
            self.lattice,
            self.species_and_occu,
            self.frac_coords,
            self.rotation_axis,
            self.rotation_angle,
            self.gb_plane,
            self.join_plane,
            self.init_cell,
            self.vacuum_thickness,
            self.ab_shift,
            self.site_properties,
            self.oriented_unit_cell,
        )

    def get_sorted_structure(
        self,
        key: Callable | None = None,
        reverse: bool = False,
    ) -> Self:
        """Get a sorted copy of the Structure. The parameters have the same
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
            self.rotation_axis,
            self.rotation_angle,
            self.gb_plane,
            self.join_plane,
            self.init_cell,
            self.vacuum_thickness,
            self.ab_shift,
            self.site_properties,
            self.oriented_unit_cell,
        )

    @property
    def sigma(self) -> int:
        """The sigma value of the GB. If using 'quick_gen' to generate GB, this value is not valid."""
        return int(round(self.oriented_unit_cell.volume / self.init_cell.volume))

    @property
    def sigma_from_site_prop(self) -> int:
        """The sigma value of the GB from site properties.
        If the GB structure merge some atoms due to the atoms too close with
        each other, this property will not work.
        """
        if None in self.site_properties["grain_label"]:
            raise ValueError("Sites were merged, this property does not work")

        n_coi = sum("incident" in tag for tag in self.site_properties["grain_label"])
        return round(len(self) / n_coi)

    @property
    def top_grain(self) -> Structure:
        """The top grain (Structure) of the GB."""
        top_sites = []
        for i, tag in enumerate(self.site_properties["grain_label"]):
            if "top" in tag:
                top_sites.append(self.sites[i])
        return Structure.from_sites(top_sites)

    @property
    def bottom_grain(self) -> Structure:
        """The bottom grain (Structure) of the GB."""
        bottom_sites = []
        for i, tag in enumerate(self.site_properties["grain_label"]):
            if "bottom" in tag:
                bottom_sites.append(self.sites[i])
        return Structure.from_sites(bottom_sites)

    @property
    def coincidents(self) -> list[Site]:
        """A list of coincident sites."""
        coincident_sites = []
        for idx, tag in enumerate(self.site_properties["grain_label"]):
            if "incident" in tag:
                coincident_sites.append(self.sites[idx])
        return coincident_sites

    def as_dict(self) -> dict:
        """
        Returns:
            Dictionary representation of GrainBoundary object.
        """
        return {
            **super().as_dict(),
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "init_cell": self.init_cell.as_dict(),
            "rotation_axis": self.rotation_axis,
            "rotation_angle": self.rotation_angle,
            "gb_plane": self.gb_plane,
            "join_plane": self.join_plane,
            "vacuum_thickness": self.vacuum_thickness,
            "ab_shift": self.ab_shift,
            "oriented_unit_cell": self.oriented_unit_cell.as_dict(),
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Generate GrainBoundary from a dict created by as_dict().

        Args:
            dct: dict

        Returns:
            GrainBoundary object
        """
        lattice = Lattice.from_dict(dct["lattice"])
        sites = [PeriodicSite.from_dict(site_dict, lattice) for site_dict in dct["sites"]]
        struct = Structure.from_sites(sites)

        return cls(
            lattice=lattice,
            species=struct.species_and_occu,
            coords=struct.frac_coords,
            rotation_axis=dct["rotation_axis"],
            rotation_angle=dct["rotation_angle"],
            gb_plane=dct["gb_plane"],
            join_plane=dct["join_plane"],
            init_cell=Structure.from_dict(dct["init_cell"]),
            vacuum_thickness=dct["vacuum_thickness"],
            ab_shift=dct["ab_shift"],
            oriented_unit_cell=Structure.from_dict(dct["oriented_unit_cell"]),
            site_properties=struct.site_properties,
        )


class GrainBoundaryGenerator:
    """
    Generate grain boundaries (GBs) from bulk conventional cell (FCC, BCC can
    from the primitive cell), and works for Cubic, Tetragonal, Orthorhombic,
    Rhombohedral, and Hexagonal systems. It generate GBs from given parameters,
    which includes GB plane, rotation axis, rotation angle.

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

    def __init__(
        self,
        initial_structure: Structure,
        symprec: float = 0.1,
        angle_tolerance: float = 1.0,
    ) -> None:
        """
        Args:
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

        # Use the conventional cell for tetragonal
        if self.lat_type == "t":
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

        elif self.lat_type == "h":
            alpha, beta, gamma = initial_structure.lattice.angles
            # c axis is not in the third direction
            if abs(gamma - 90) < angle_tolerance:
                # alpha = 120 or 60, rotate b, c to a, b vectors
                if abs(alpha - 90) > angle_tolerance:
                    initial_structure.make_supercell([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
                # beta = 120 or 60, rotate c, a to a, b vectors
                elif abs(beta - 90) > angle_tolerance:
                    initial_structure.make_supercell([[0, 0, 1], [1, 0, 0], [0, 1, 0]])

        # Use primitive cell for rhombohedra
        elif self.lat_type == "r":
            initial_structure = analyzer.get_primitive_standard_structure()

        # Use the conventional cell for orthorhombic
        elif self.lat_type == "o":
            initial_structure = analyzer.get_conventional_standard_structure()

        self.initial_structure = initial_structure

    def gb_from_parameters(
        self,
        rotation_axis: Tuple3Ints,
        rotation_angle: float,
        expand_times: int = 4,
        vacuum_thickness: float = 0.0,
        ab_shift: tuple[float, float] = (0, 0),
        normal: bool = False,
        ratio: list[int] | None = None,
        plane: Tuple3Ints | None = None,
        max_search: int = 20,
        tol_coi: float = 1.0e-8,
        rm_ratio: float = 0.7,
        quick_gen: bool = False,
    ) -> GrainBoundary:
        """
        Args:
            rotation_axis (tuple of 3): Rotation axis of GB e.g.: (1, 1, 0).
            rotation_angle (float, in unit of degree): rotation angle used to generate GB.
                Make sure the angle is accurate enough. You can use the enum* functions
                in this class to extract the accurate angle.
                e.g.: The rotation angle of sigma 3 twist GB with the rotation axis
                (1, 1, 1) and GB plane (1, 1, 1) can be 60 degree.
                If you do not know the rotation angle, but know the sigma value, we have
                provide the function get_rotation_angle_from_sigma which is able to return
                all the rotation angles of sigma value you provided.
            expand_times (int): The multiple times used to expand one unit grain to larger grain.
                This is used to tune the grain length of GB to warrant that the two GBs in one
                cell do not interact with each other. Default set to 4.
            vacuum_thickness (float, in angstrom): The thickness of vacuum that you want to insert
                between two grains of the GB. Default to 0.
            ab_shift (list of float, in unit of a, b vectors of Gb): in plane shift of two grains
            normal (bool):
                determine if need to require the c axis of top grain (first transformation matrix)
                perpendicular to the surface or not.
                default to false.
            ratio (list[int]): lattice axial ratio.
                For cubic system, ratio is not needed.
                For tetragonal system, ratio = [mu, mv], list of two integers,
                that is, mu/mv = c2/a2. If it is irrational, set it to none.
                For orthorhombic system, ratio = [mu, lam, mv], list of 3 integers,
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
            plane (tuple of 3): Grain boundary plane. If none, we set it as twist GB.
                The plane will be perpendicular to the rotation axis.
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
            GrainBoundary object
        """
        lat_type = self.lat_type.lower()
        # If the initial structure is primitive cell in cubic system,
        # calculate the transformation matrix from its conventional cell
        # to primitive cell, basically for BCC and FCC systems.
        trans_cry = np.eye(3)
        if lat_type == "c":
            analyzer = SpacegroupAnalyzer(self.initial_structure)
            convention_cell = analyzer.get_conventional_standard_structure()
            vol_ratio = self.initial_structure.volume / convention_cell.volume
            # BCC primitive cell, belong to cubic system
            if abs(vol_ratio - 0.5) < 1.0e-3:
                trans_cry = np.array([[0.5, 0.5, -0.5], [-0.5, 0.5, 0.5], [0.5, -0.5, 0.5]])
                logger.info("Make sure this is for cubic with bcc primitive cell")
            # FCC primitive cell, belong to cubic system
            elif abs(vol_ratio - 0.25) < 1.0e-3:
                trans_cry = np.array([[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]])
                logger.info("Make sure this is for cubic with fcc primitive cell")
            else:
                logger.info("Make sure this is for cubic with conventional cell")

        elif lat_type == "t":
            logger.info("Make sure this is for tetragonal system")
            if ratio is None:
                logger.info("Make sure this is for irrational c2/a2")
            elif len(ratio) != 2:
                raise RuntimeError("Tetragonal system needs correct c2/a2 ratio")

        elif lat_type == "o":
            logger.info("Make sure this is for orthorhombic system")
            if ratio is None:
                raise RuntimeError("CSL does not exist if all axial ratios are irrational for an orthorhombic system")
            if len(ratio) != 3:
                raise RuntimeError("Orthorhombic system needs correct c2:b2:a2 ratio")

        elif lat_type == "h":
            logger.info("Make sure this is for hexagonal system")
            if ratio is None:
                logger.info("Make sure this is for irrational c2/a2")
            elif len(ratio) != 2:
                raise RuntimeError("Hexagonal system needs correct c2/a2 ratio")

        elif lat_type == "r":
            logger.info("Make sure this is for rhombohedral system")
            if ratio is None:
                logger.info("Make sure this is for irrational (1+2*cos(alpha)/cos(alpha) ratio")
            elif len(ratio) != 2:
                raise RuntimeError("Rhombohedral system needs correct (1+2*cos(alpha)/cos(alpha) ratio")

        else:
            raise RuntimeError(
                "Lattice type not implemented. This code works for cubic, "
                "tetragonal, orthorhombic, rhombohedral, hexagonal systems"
            )

        # Transform four index notation to three index notation for hexagonal and rhombohedral
        if len(rotation_axis) == 4:
            u1 = rotation_axis[0]
            v1 = rotation_axis[1]
            w1 = rotation_axis[3]
            if lat_type == "h":
                u = 2 * u1 + v1
                v = 2 * v1 + u1
                w = w1
                _rotation_axis: Tuple3Ints | Tuple4Ints = (u, v, w)
            elif lat_type == "r":
                u = 2 * u1 + v1 + w1
                v = v1 + w1 - u1
                w = w1 - 2 * v1 - u1
                _rotation_axis = (u, v, w)
            else:
                _rotation_axis = cast(Tuple4Ints, tuple(rotation_axis))

        elif len(rotation_axis) == 3:
            _rotation_axis = cast(Tuple3Ints, tuple(rotation_axis))

        else:
            raise ValueError("Invalid length of rotation axis.")

        # Make sure math.gcd(rotation_axis) == 1
        if reduce(math.gcd, _rotation_axis) != 1:
            _rotation_axis = tuple(round(x / reduce(math.gcd, _rotation_axis)) for x in _rotation_axis)  # type: ignore[assignment]

        # Transform four index notation to three index notation for plane
        if plane is not None:
            if len(plane) == 4:
                u1, v1, w1 = plane[0], plane[1], plane[3]
                plane = (u1, v1, w1)
            elif len(plane) == 3:
                plane = cast(Tuple3Ints, tuple(plane))

        # Set the plane for grain boundary when plane is None
        if plane is None:
            if lat_type == "c" and len(_rotation_axis) == 3:
                _plane = _rotation_axis
            else:
                if lat_type == "h":
                    c2_a2_ratio = 1.0 if ratio is None else ratio[0] / ratio[1]
                    metric = np.array([[1, -0.5, 0], [-0.5, 1, 0], [0, 0, c2_a2_ratio]])
                elif lat_type == "r":
                    cos_alpha = 0.5 if ratio is None else 1.0 / (ratio[0] / ratio[1] - 2)
                    metric = np.array(
                        [
                            [1, cos_alpha, cos_alpha],
                            [cos_alpha, 1, cos_alpha],
                            [cos_alpha, cos_alpha, 1],
                        ]
                    )
                elif lat_type == "t":
                    c2_a2_ratio = 1.0 if ratio is None else ratio[0] / ratio[1]
                    metric = np.array([[1, 0, 0], [0, 1, 0], [0, 0, c2_a2_ratio]])
                elif lat_type == "o" and ratio is not None:
                    for idx in range(3):
                        if ratio[idx] is None:
                            ratio[idx] = 1
                    metric = np.array(
                        [
                            [1, 0, 0],
                            [0, ratio[1] / ratio[2], 0],
                            [0, 0, ratio[0] / ratio[2]],
                        ]
                    )
                else:
                    raise RuntimeError("Lattice type is not implemented.")

                _plane = np.matmul(_rotation_axis, metric)
                fractions = [Fraction(x).limit_denominator() for x in _plane]
                least_mul = reduce(lcm, [fraction.denominator for fraction in fractions])
                _plane = cast(Tuple3Ints, tuple(round(x * least_mul) for x in _plane))

        else:
            _plane = plane

        if reduce(math.gcd, _plane) != 1:
            index = reduce(math.gcd, _plane)
            _plane = cast(Tuple3Ints, tuple(round(x / index) for x in _plane))

        t1, t2 = self.get_trans_mat(
            r_axis=_rotation_axis,
            angle=rotation_angle,
            normal=normal,
            trans_cry=trans_cry,
            lat_type=lat_type,
            ratio=ratio,
            surface=_plane,
            max_search=max_search,
            quick_gen=quick_gen,
        )

        # Find the join_plane
        if lat_type != "c":
            if lat_type == "h":
                if ratio is None:
                    mu, mv = [1, 1]
                else:
                    mu, mv = ratio
                trans_cry1 = np.array([[1, 0, 0], [-0.5, np.sqrt(3.0) / 2.0, 0], [0, 0, np.sqrt(mu / mv)]])

            elif lat_type == "r":
                if ratio is None:
                    c2_a2_ratio = 1.0
                else:
                    mu, mv = ratio
                    c2_a2_ratio = 3 / (2 - 6 * mv / mu)
                trans_cry1 = np.array(
                    [
                        [0.5, np.sqrt(3.0) / 6.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)],
                        [-0.5, np.sqrt(3.0) / 6.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)],
                        [0, -1 * np.sqrt(3.0) / 3.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)],
                    ]
                )

            else:
                if lat_type == "t":
                    if ratio is None:
                        mu, mv = [1, 1]
                    else:
                        mu, mv = ratio
                    lam = mv

                elif lat_type == "o" and ratio is not None:
                    new_ratio = [1 if v is None else v for v in ratio]
                    mu, lam, mv = new_ratio

                else:
                    raise RuntimeError("Invalid lattice type.")

                trans_cry1 = np.array([[1, 0, 0], [0, np.sqrt(lam / mv), 0], [0, 0, np.sqrt(mu / mv)]])

        else:
            trans_cry1 = trans_cry

        grain_matrix = np.dot(t2, trans_cry1)
        plane_init = np.cross(grain_matrix[0], grain_matrix[1])
        if lat_type != "c":
            plane_init = np.dot(plane_init, trans_cry1.T)
        join_plane = self.vec_to_surface(plane_init)

        parent_structure = self.initial_structure.copy()
        # Calculate the bond_length in bulk system
        if len(parent_structure) == 1:
            temp_str = parent_structure.copy()
            temp_str.make_supercell([1, 1, 2])
            distance = temp_str.distance_matrix
        else:
            distance = parent_structure.distance_matrix
        bond_length = np.min(distance[np.nonzero(distance)])

        # Top grain
        top_grain = fix_pbc(parent_structure * t1)

        # Obtain the smallest oriented cell
        if normal and not quick_gen:
            t_temp = self.get_trans_mat(
                r_axis=_rotation_axis,
                angle=rotation_angle,
                normal=False,
                trans_cry=trans_cry,
                lat_type=lat_type,
                ratio=ratio,
                surface=_plane,
                max_search=max_search,
            )
            oriented_unit_cell = fix_pbc(parent_structure * t_temp[0])
            t_matrix = oriented_unit_cell.lattice.matrix
            normal_v_plane = np.cross(t_matrix[0], t_matrix[1])
            unit_normal_v = normal_v_plane / np.linalg.norm(normal_v_plane)
            unit_ab_adjust = (t_matrix[2] - np.dot(unit_normal_v, t_matrix[2]) * unit_normal_v) / np.dot(
                unit_normal_v, t_matrix[2]
            )
        else:
            oriented_unit_cell = top_grain.copy()
            unit_ab_adjust = 0.0

        # Bottom grain, using top grain's lattice matrix
        bottom_grain = fix_pbc(parent_structure * t2, top_grain.lattice.matrix)

        # Label both grains with 'top', 'bottom', 'top_incident', 'bottom_incident'
        n_sites = len(top_grain)
        t_and_b = Structure(
            top_grain.lattice,
            top_grain.species + bottom_grain.species,
            list(top_grain.frac_coords) + list(bottom_grain.frac_coords),
        )
        t_and_b_dis = t_and_b.lattice.get_all_distances(
            t_and_b.frac_coords[:n_sites], t_and_b.frac_coords[n_sites : n_sites * 2]
        )
        index_incident = np.nonzero(t_and_b_dis < np.min(t_and_b_dis) + tol_coi)

        top_labels = []
        for idx in range(n_sites):
            if idx in index_incident[0]:
                top_labels.append("top_incident")
            else:
                top_labels.append("top")
        bottom_labels = []
        for idx in range(n_sites):
            if idx in index_incident[1]:
                bottom_labels.append("bottom_incident")
            else:
                bottom_labels.append("bottom")
        top_grain = Structure(
            Lattice(top_grain.lattice.matrix),
            top_grain.species,
            top_grain.frac_coords,
            site_properties={"grain_label": top_labels},
        )
        bottom_grain = Structure(
            Lattice(bottom_grain.lattice.matrix),
            bottom_grain.species,
            bottom_grain.frac_coords,
            site_properties={"grain_label": bottom_labels},
        )

        # Expand both grains
        top_grain.make_supercell([1, 1, expand_times])
        bottom_grain.make_supercell([1, 1, expand_times])
        top_grain = fix_pbc(top_grain)
        bottom_grain = fix_pbc(bottom_grain)

        # Determine the top-grain location
        edge_b = 1.0 - max(bottom_grain.frac_coords[:, 2])
        edge_t = 1.0 - max(top_grain.frac_coords[:, 2])
        c_adjust = (edge_t - edge_b) / 2.0

        # Construct all species
        all_species = []
        all_species.extend([site.specie for site in bottom_grain])
        all_species.extend([site.specie for site in top_grain])

        half_lattice = top_grain.lattice
        # Calculate translation vector, perpendicular to the plane
        normal_v_plane = np.cross(half_lattice.matrix[0], half_lattice.matrix[1])
        unit_normal_v = normal_v_plane / np.linalg.norm(normal_v_plane)
        translation_v = unit_normal_v * vacuum_thickness

        # Construct the final lattice
        whole_matrix_no_vac = np.array(half_lattice.matrix)
        whole_matrix_no_vac[2] = half_lattice.matrix[2] * 2
        whole_matrix_with_vac = whole_matrix_no_vac.copy()
        whole_matrix_with_vac[2] = whole_matrix_no_vac[2] + translation_v * 2
        whole_lat = Lattice(whole_matrix_with_vac)

        # Construct the coords, move top grain with translation_v
        all_coords = []
        grain_labels = bottom_grain.site_properties["grain_label"] + top_grain.site_properties["grain_label"]  # type: ignore[operator]
        for site in bottom_grain:
            all_coords.append(site.coords)
        for site in top_grain:
            all_coords.append(
                site.coords
                + half_lattice.matrix[2] * (1 + c_adjust)
                + unit_ab_adjust * np.linalg.norm(half_lattice.matrix[2] * (1 + c_adjust))
                + translation_v
                + ab_shift[0] * whole_matrix_with_vac[0]
                + ab_shift[1] * whole_matrix_with_vac[1]
            )

        gb_with_vac = Structure(
            whole_lat,
            all_species,
            all_coords,
            coords_are_cartesian=True,
            site_properties={"grain_label": grain_labels},
        )
        # Merge closer atoms. extract near GB atoms.
        cos_c_norm_plane = np.dot(unit_normal_v, whole_matrix_with_vac[2]) / whole_lat.c
        range_c_len = abs(bond_length / cos_c_norm_plane / whole_lat.c)
        sites_near_gb = []
        sites_away_gb: list[PeriodicSite] = []
        for site in gb_with_vac:
            if (
                site.frac_coords[2] < range_c_len
                or site.frac_coords[2] > 1 - range_c_len
                or (site.frac_coords[2] > 0.5 - range_c_len and site.frac_coords[2] < 0.5 + range_c_len)
            ):
                sites_near_gb.append(site)
            else:
                sites_away_gb.append(site)
        if len(sites_near_gb) >= 1:
            s_near_gb = Structure.from_sites(sites_near_gb)
            s_near_gb.merge_sites(tol=bond_length * rm_ratio, mode="delete")
            all_sites = sites_away_gb + s_near_gb.sites  # type: ignore[operator]
            gb_with_vac = Structure.from_sites(all_sites)

        # Move coordinates into the periodic cell
        gb_with_vac = fix_pbc(gb_with_vac, whole_lat.matrix)
        return GrainBoundary(
            whole_lat,
            gb_with_vac.species,
            gb_with_vac.cart_coords,  # type: ignore[arg-type]
            _rotation_axis,
            rotation_angle,
            _plane,
            join_plane,
            self.initial_structure,
            vacuum_thickness,
            ab_shift,
            site_properties=gb_with_vac.site_properties,
            oriented_unit_cell=oriented_unit_cell,
            coords_are_cartesian=True,
        )

    def get_ratio(
        self,
        max_denominator: int = 5,
        index_none: int | None = None,
    ) -> list[int] | None:
        """Find the axial ratio needed for GB generator input.

        Args:
            max_denominator (int): the maximum denominator for
                the computed ratio, default to be 5.
            index_none (int): specify the irrational axis.
                0-a, 1-b, 2-c. Only may be needed for orthorhombic system.

        Returns:
            axial ratio needed for GB generator (list of integers).
        """
        structure = self.initial_structure
        lat_type = self.lat_type

        # For tetragonal and hexagonal systems, ratio = c2 / a2
        if lat_type in {"t", "h"}:
            a, _, c = structure.lattice.lengths
            if c > a:
                frac = Fraction(c**2 / a**2).limit_denominator(max_denominator)
                ratio: list[int | None] = [frac.numerator, frac.denominator]
            else:
                frac = Fraction(a**2 / c**2).limit_denominator(max_denominator)
                ratio = [frac.denominator, frac.numerator]

        # For rhombohedral system, ratio = (1 + 2 * cos(alpha)) / cos(alpha)
        elif lat_type == "r":
            cos_alpha = math.cos(structure.lattice.alpha / 180 * np.pi)
            frac = Fraction((1 + 2 * cos_alpha) / cos_alpha).limit_denominator(max_denominator)
            ratio = [frac.numerator, frac.denominator]

        # For orthorhombic system, ratio = c2:b2:a2. If irrational for one axis, set it to None
        elif lat_type == "o":
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
                ratio[index[0]] = frac1.numerator * round(com_lcm / frac1.denominator)
                ratio[index[1]] = frac2.numerator * round(com_lcm / frac2.denominator)

            else:
                index.pop(index_none)
                if lat[index[0]] > lat[index[1]]:
                    frac = Fraction(lat[index[0]] ** 2 / lat[index[1]] ** 2).limit_denominator(max_denominator)
                    ratio[index[0]] = frac.numerator
                    ratio[index[1]] = frac.denominator
                else:
                    frac = Fraction(lat[index[1]] ** 2 / lat[index[0]] ** 2).limit_denominator(max_denominator)
                    ratio[index[1]] = frac.numerator
                    ratio[index[0]] = frac.denominator

        # Cubic system does not need axial ratio
        elif lat_type == "c":
            return None

        else:
            raise RuntimeError("Lattice type not implemented.")
        return cast(list[int], ratio)

    @staticmethod
    def get_trans_mat(
        r_axis: Tuple3Ints | Tuple4Ints,
        angle: float,
        normal: bool = False,
        trans_cry: NDArray | None = None,
        lat_type: str = "c",
        ratio: list[int] | None = None,
        surface: Tuple3Ints | Tuple4Ints | None = None,
        max_search: int = 20,
        quick_gen: bool = False,
    ):
        """Find the two transformation matrix for each grain from given rotation axis,
        GB plane, rotation angle and corresponding ratio (see explanation for ratio
        below).
        The structure of each grain can be obtained by applying the corresponding
        transformation matrix to the conventional cell.
        The algorithm for this code is from reference, Acta Cryst, A32,783(1976).

        Args:
            r_axis ((u, v, w) or (u, v, t, w) for hex/rho systems):
                the rotation axis of the grain boundary.
            angle (float, in unit of degree): the rotation angle of the grain boundary
            normal (bool): determine if need to require the c axis of one grain associated with
                the first transformation matrix perpendicular to the surface or not.
                default to false.
            trans_cry (np.array): shape 3x3. If the structure given are primitive cell in cubic system, e.g.
                bcc or fcc system, trans_cry is the transformation matrix from its
                conventional cell to the primitive cell.
            lat_type (str): one character to specify the lattice type. Defaults to 'c' for cubic.
                'c' or 'C': cubic system
                't' or 'T': tetragonal system
                'o' or 'O': orthorhombic system
                'h' or 'H': hexagonal system
                'r' or 'R': rhombohedral system
            ratio (list[int]): lattice axial ratio.
                For cubic system, ratio is not needed.
                For tetragonal system, ratio = [mu, mv], list of two integers, that is, mu/mv = c2/a2. If it is
                irrational, set it to none.
                For orthorhombic system, ratio = [mu, lam, mv], list of 3 integers, that is, mu:lam:mv = c2:b2:a2.
                If irrational for one axis, set it to None. e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
                For rhombohedral system, ratio = [mu, mv], list of two integers,
                that is, mu/mv is the ratio of (1+2*cos(alpha)/cos(alpha).
                If irrational, set it to None.
                For hexagonal system, ratio = [mu, mv], list of two integers,
                that is, mu/mv = c2/a2. If it is irrational, set it to none.
            surface ((h, k, l) or (h, k, i, l) for hex/rho systems): The
                miller index of grain boundary plane, with the format of (h, k, l) if surface
                is not given, the default is perpendicular to r_axis, which is a twist grain boundary.
            max_search (int): max search for the GB lattice vectors that give the smallest GB
                lattice. If normal is true, also max search the GB c vector that perpendicular
                to the plane.
            quick_gen (bool): whether to quickly generate a supercell, if set to true, no need to
                find the smallest cell.

        Returns:
            t1 (3 by 3 integer array): The transformation array for one grain.
            t2 (3 by 3 integer array): The transformation array for the other grain
        """
        trans_cry = np.eye(3) if trans_cry is None else trans_cry
        lat_type = lat_type.lower()

        # Transform four index notation to three index notation
        if len(r_axis) == 4:
            u1 = r_axis[0]
            v1 = r_axis[1]
            w1 = r_axis[3]
            if lat_type == "h":
                u = 2 * u1 + v1
                v = 2 * v1 + u1
                w = w1
                r_axis = (u, v, w)
            elif lat_type == "r":
                u = 2 * u1 + v1 + w1
                v = v1 + w1 - u1
                w = w1 - 2 * v1 - u1
                r_axis = (u, v, w)

        # Make sure gcd(r_axis) == 1
        if reduce(math.gcd, r_axis) != 1:
            r_axis = cast(
                Tuple3Ints | Tuple4Ints,
                tuple(round(x / reduce(math.gcd, r_axis)) for x in r_axis),
            )

        if surface is not None and len(surface) == 4:
            u1 = surface[0]
            v1 = surface[1]
            w1 = surface[3]
            surface = (u1, v1, w1)

        # Set the surface for grain boundary.
        if surface is None:
            if lat_type == "c":
                surface = r_axis
            else:
                if lat_type == "h":
                    c2_a2_ratio = 1.0 if ratio is None else ratio[0] / ratio[1]
                    metric = np.array([[1, -0.5, 0], [-0.5, 1, 0], [0, 0, c2_a2_ratio]])
                elif lat_type == "r":
                    cos_alpha = 0.5 if ratio is None else 1.0 / (ratio[0] / ratio[1] - 2)
                    metric = np.array(
                        [
                            [1, cos_alpha, cos_alpha],
                            [cos_alpha, 1, cos_alpha],
                            [cos_alpha, cos_alpha, 1],
                        ]
                    )
                elif lat_type == "t":
                    c2_a2_ratio = 1.0 if ratio is None else ratio[0] / ratio[1]
                    metric = np.array([[1, 0, 0], [0, 1, 0], [0, 0, c2_a2_ratio]])
                elif lat_type == "o":
                    if ratio is None:
                        raise ValueError(f"Invalid {ratio=} for orthorhombic system")
                    for idx in range(3):
                        if ratio is not None and ratio[idx] is None:
                            ratio[idx] = 1
                    metric = np.array(
                        [
                            [1, 0, 0],
                            [0, ratio[1] / ratio[2], 0],
                            [0, 0, ratio[0] / ratio[2]],
                        ]
                    )
                else:
                    raise RuntimeError("Lattice type has not implemented.")

                surface = np.matmul(r_axis, metric)
                fractions = [Fraction(x).limit_denominator() for x in surface]
                least_mul = reduce(lcm, [fraction.denominator for fraction in fractions])
                surface = cast(
                    Tuple3Ints | Tuple4Ints,
                    tuple(round(x * least_mul) for x in surface),
                )

        if reduce(math.gcd, surface) != 1:
            index = reduce(math.gcd, surface)
            surface = cast(Tuple3Ints | Tuple4Ints, tuple(round(x / index) for x in surface))

        lam = None
        if lat_type == "h":
            # Set the values for u, v, w, mu, mv, m, n, d, x
            # Check the reference for the meaning of these parameters
            u, v, w = cast(Tuple3Ints, r_axis)
            # Make sure mu, mv are coprime integers
            if ratio is None:
                mu, mv = (1, 1)
                if w != 0 and (u != 0 or (v != 0)):
                    raise RuntimeError("For irrational c2/a2, CSL only exist for [0,0,1] or [u,v,0] and m = 0")
            else:
                mu, mv = ratio
            if math.gcd(mu, mv) != 1:
                temp = math.gcd(mu, mv)
                mu = round(mu / temp)
                mv = round(mv / temp)
            d = (u**2 + v**2 - u * v) * mv + w**2 * mu
            if abs(angle - 180.0) < 1.0e0:
                m = 0
                n = 1
            else:
                fraction = Fraction(
                    np.tan(angle / 2 / 180.0 * np.pi) / np.sqrt(float(d) / 3.0 / mu)
                ).limit_denominator()
                m = fraction.denominator
                n = fraction.numerator

            # Construct the rotation matrix, check reference for details
            r_list = [
                (u**2 * mv - v**2 * mv - w**2 * mu) * n**2 + 2 * w * mu * m * n + 3 * mu * m**2,
                (2 * v - u) * u * mv * n**2 - 4 * w * mu * m * n,
                2 * u * w * mu * n**2 + 2 * (2 * v - u) * mu * m * n,
                (2 * u - v) * v * mv * n**2 + 4 * w * mu * m * n,
                (v**2 * mv - u**2 * mv - w**2 * mu) * n**2 - 2 * w * mu * m * n + 3 * mu * m**2,
                2 * v * w * mu * n**2 - 2 * (2 * u - v) * mu * m * n,
                (2 * u - v) * w * mv * n**2 - 3 * v * mv * m * n,
                (2 * v - u) * w * mv * n**2 + 3 * u * mv * m * n,
                (w**2 * mu - u**2 * mv - v**2 * mv + u * v * mv) * n**2 + 3 * mu * m**2,
            ]
            m = -1 * m
            r_list_inv = [
                (u**2 * mv - v**2 * mv - w**2 * mu) * n**2 + 2 * w * mu * m * n + 3 * mu * m**2,
                (2 * v - u) * u * mv * n**2 - 4 * w * mu * m * n,
                2 * u * w * mu * n**2 + 2 * (2 * v - u) * mu * m * n,
                (2 * u - v) * v * mv * n**2 + 4 * w * mu * m * n,
                (v**2 * mv - u**2 * mv - w**2 * mu) * n**2 - 2 * w * mu * m * n + 3 * mu * m**2,
                2 * v * w * mu * n**2 - 2 * (2 * u - v) * mu * m * n,
                (2 * u - v) * w * mv * n**2 - 3 * v * mv * m * n,
                (2 * v - u) * w * mv * n**2 + 3 * u * mv * m * n,
                (w**2 * mu - u**2 * mv - v**2 * mv + u * v * mv) * n**2 + 3 * mu * m**2,
            ]
            m = -1 * m
            F = 3 * mu * m**2 + d * n**2
            all_list = r_list + r_list_inv + [F]
            com_fac = reduce(math.gcd, all_list)
            sigma = F / com_fac
            r_matrix = (np.array(r_list) / com_fac / sigma).reshape(3, 3)

        elif lat_type == "r":
            # Set the values for u, v, w, mu, mv ,m, n, d
            # Check the reference for the meaning of these parameters
            u, v, w = cast(Tuple3Ints, r_axis)
            # make sure mu, mv are coprime integers
            if ratio is None:
                mu, mv = (1, 1)
                if u + v + w != 0 and (u != v or u != w):
                    raise RuntimeError(
                        "For irrational ratio_alpha, CSL only exist for [1,1,1] or [u, v, -(u+v)] and m =0"
                    )
            else:
                mu, mv = ratio
            if math.gcd(mu, mv) != 1:
                temp = math.gcd(mu, mv)
                mu = round(mu / temp)
                mv = round(mv / temp)
            d = (u**2 + v**2 + w**2) * (mu - 2 * mv) + 2 * mv * (v * w + w * u + u * v)
            if abs(angle - 180.0) < 1.0e0:
                m = 0
                n = 1
            else:
                fraction = Fraction(np.tan(angle / 2 / 180.0 * np.pi) / np.sqrt(float(d) / mu)).limit_denominator()
                m = fraction.denominator
                n = fraction.numerator

            # Construct the rotation matrix, check reference for details
            r_list = [
                (mu - 2 * mv) * (u**2 - v**2 - w**2) * n**2
                + 2 * mv * (v - w) * m * n
                - 2 * mv * v * w * n**2
                + mu * m**2,
                2 * (mv * u * n * (w * n + u * n - m) - (mu - mv) * m * w * n + (mu - 2 * mv) * u * v * n**2),
                2 * (mv * u * n * (v * n + u * n + m) + (mu - mv) * m * v * n + (mu - 2 * mv) * w * u * n**2),
                2 * (mv * v * n * (w * n + v * n + m) + (mu - mv) * m * w * n + (mu - 2 * mv) * u * v * n**2),
                (mu - 2 * mv) * (v**2 - w**2 - u**2) * n**2
                + 2 * mv * (w - u) * m * n
                - 2 * mv * u * w * n**2
                + mu * m**2,
                2 * (mv * v * n * (v * n + u * n - m) - (mu - mv) * m * u * n + (mu - 2 * mv) * w * v * n**2),
                2 * (mv * w * n * (w * n + v * n - m) - (mu - mv) * m * v * n + (mu - 2 * mv) * w * u * n**2),
                2 * (mv * w * n * (w * n + u * n + m) + (mu - mv) * m * u * n + (mu - 2 * mv) * w * v * n**2),
                (mu - 2 * mv) * (w**2 - u**2 - v**2) * n**2
                + 2 * mv * (u - v) * m * n
                - 2 * mv * u * v * n**2
                + mu * m**2,
            ]
            m = -1 * m
            r_list_inv = [
                (mu - 2 * mv) * (u**2 - v**2 - w**2) * n**2
                + 2 * mv * (v - w) * m * n
                - 2 * mv * v * w * n**2
                + mu * m**2,
                2 * (mv * u * n * (w * n + u * n - m) - (mu - mv) * m * w * n + (mu - 2 * mv) * u * v * n**2),
                2 * (mv * u * n * (v * n + u * n + m) + (mu - mv) * m * v * n + (mu - 2 * mv) * w * u * n**2),
                2 * (mv * v * n * (w * n + v * n + m) + (mu - mv) * m * w * n + (mu - 2 * mv) * u * v * n**2),
                (mu - 2 * mv) * (v**2 - w**2 - u**2) * n**2
                + 2 * mv * (w - u) * m * n
                - 2 * mv * u * w * n**2
                + mu * m**2,
                2 * (mv * v * n * (v * n + u * n - m) - (mu - mv) * m * u * n + (mu - 2 * mv) * w * v * n**2),
                2 * (mv * w * n * (w * n + v * n - m) - (mu - mv) * m * v * n + (mu - 2 * mv) * w * u * n**2),
                2 * (mv * w * n * (w * n + u * n + m) + (mu - mv) * m * u * n + (mu - 2 * mv) * w * v * n**2),
                (mu - 2 * mv) * (w**2 - u**2 - v**2) * n**2
                + 2 * mv * (u - v) * m * n
                - 2 * mv * u * v * n**2
                + mu * m**2,
            ]
            m = -1 * m
            F = mu * m**2 + d * n**2
            all_list = r_list_inv + r_list + [F]
            com_fac = reduce(math.gcd, all_list)
            sigma = F / com_fac
            r_matrix = (np.array(r_list) / com_fac / sigma).reshape(3, 3)

        else:
            u, v, w = cast(Tuple3Ints, r_axis)
            mu = mv = None  # type: ignore[assignment]
            if lat_type == "c":
                mu = lam = mv = 1
            elif lat_type == "t":
                if ratio is None:
                    mu, mv = (1, 1)
                    if w != 0 and (u != 0 or (v != 0)):
                        raise RuntimeError("For irrational c2/a2, CSL only exist for [0,0,1] or [u,v,0] and m = 0")
                else:
                    mu, mv = ratio
                lam = mv
            elif lat_type == "o":
                if ratio is not None and None in ratio:
                    mu, lam, mv = ratio
                    non_none = [i for i in ratio if i is not None]
                    if len(non_none) < 2:
                        raise RuntimeError("No CSL exist for two irrational numbers")
                    non1, non2 = non_none
                    if mu is None:
                        lam = non1
                        mv = non2
                        mu = 1
                        if w != 0 and (u != 0 or (v != 0)):
                            raise RuntimeError("For irrational c2, CSL only exist for [0,0,1] or [u,v,0] and m = 0")
                    elif lam is None:
                        mu = non1
                        mv = non2
                        lam = 1
                        if v != 0 and (u != 0 or (w != 0)):
                            raise RuntimeError("For irrational b2, CSL only exist for [0,1,0] or [u,0,w] and m = 0")
                    elif mv is None:
                        mu = non1
                        lam = non2
                        mv = 1
                        if u != 0 and (w != 0 or (v != 0)):
                            raise RuntimeError("For irrational a2, CSL only exist for [1,0,0] or [0,v,w] and m = 0")
                else:
                    mu, lam, mv = cast(list[int], ratio)
                    if u == 0 and v == 0:
                        mu = 1
                    if u == 0 and w == 0:
                        lam = 1
                    if v == 0 and w == 0:
                        mv = 1

            # Make sure mu, lambda, mv are coprime integers
            if mu is None:
                raise ValueError("mu is None.")
            if lam is None:
                raise ValueError("lambda is None.")
            if mv is None:
                raise ValueError("mv is None.")

            if reduce(math.gcd, [mu, lam, mv]) != 1:
                temp = cast(int, reduce(math.gcd, [mu, lam, mv]))
                mu = round(mu / temp)
                mv = round(mv / temp)
                lam = round(lam / temp)
            d = (mv * u**2 + lam * v**2) * mv + w**2 * mu * mv
            if abs(angle - 180.0) < 1.0e0:
                m = 0
                n = 1
            else:
                fraction = Fraction(np.tan(angle / 2 / 180.0 * np.pi) / np.sqrt(d / mu / lam)).limit_denominator()
                m = fraction.denominator
                n = fraction.numerator
            r_list = [
                (u**2 * mv * mv - lam * v**2 * mv - w**2 * mu * mv) * n**2 + lam * mu * m**2,
                2 * lam * (v * u * mv * n**2 - w * mu * m * n),
                2 * mu * (u * w * mv * n**2 + v * lam * m * n),
                2 * mv * (u * v * mv * n**2 + w * mu * m * n),
                (v**2 * mv * lam - u**2 * mv * mv - w**2 * mu * mv) * n**2 + lam * mu * m**2,
                2 * mv * mu * (v * w * n**2 - u * m * n),
                2 * mv * (u * w * mv * n**2 - v * lam * m * n),
                2 * lam * mv * (v * w * n**2 + u * m * n),
                (w**2 * mu * mv - u**2 * mv * mv - v**2 * mv * lam) * n**2 + lam * mu * m**2,
            ]
            m = -1 * m
            r_list_inv = [
                (u**2 * mv * mv - lam * v**2 * mv - w**2 * mu * mv) * n**2 + lam * mu * m**2,
                2 * lam * (v * u * mv * n**2 - w * mu * m * n),
                2 * mu * (u * w * mv * n**2 + v * lam * m * n),
                2 * mv * (u * v * mv * n**2 + w * mu * m * n),
                (v**2 * mv * lam - u**2 * mv * mv - w**2 * mu * mv) * n**2 + lam * mu * m**2,
                2 * mv * mu * (v * w * n**2 - u * m * n),
                2 * mv * (u * w * mv * n**2 - v * lam * m * n),
                2 * lam * mv * (v * w * n**2 + u * m * n),
                (w**2 * mu * mv - u**2 * mv * mv - v**2 * mv * lam) * n**2 + lam * mu * m**2,
            ]
            m = -1 * m
            F = mu * lam * m**2 + d * n**2
            all_list = r_list + r_list_inv + [F]
            com_fac = reduce(math.gcd, all_list)
            sigma = F / com_fac
            r_matrix = (np.array(r_list) / com_fac / sigma).reshape(3, 3)

        if sigma > 1000:
            raise RuntimeError("Sigma >1000 too large. Are you sure what you are doing, Please check the GB if exist")
        # Transform surface, r_axis, r_matrix in terms of primitive lattice
        surface = np.matmul(surface, np.transpose(trans_cry))
        if surface is None:
            raise ValueError("surface is None.")
        fractions = [Fraction(x).limit_denominator() for x in surface]
        least_mul = reduce(lcm, [fraction.denominator for fraction in fractions])
        surface = cast(Tuple3Ints, tuple(round(x * least_mul) for x in surface))
        if reduce(math.gcd, surface) != 1:
            index = reduce(math.gcd, surface)
            surface = cast(Tuple3Ints, tuple(round(x / index) for x in surface))
        r_axis = np.rint(np.matmul(r_axis, np.linalg.inv(trans_cry))).astype(int)
        if reduce(math.gcd, r_axis) != 1:
            r_axis = cast(
                Tuple3Ints,
                tuple(round(x / reduce(math.gcd, r_axis)) for x in r_axis),
            )
        r_matrix = np.dot(np.dot(np.linalg.inv(trans_cry.T), r_matrix), trans_cry.T)
        # Set one vector of the basis to the rotation axis direction, and
        # obtain the corresponding transform matrix
        eye = np.eye(3, dtype=np.int64)
        hh = kk = ll = None
        for hh in range(3):
            if abs(r_axis[hh]) != 0:
                eye[hh] = np.array(r_axis)
                kk = hh + 1 if hh + 1 < 3 else abs(2 - hh)
                ll = hh + 2 if hh + 2 < 3 else abs(1 - hh)
                break
        trans = eye.T
        new_rot = np.array(r_matrix)

        # With the rotation matrix to construct the CSL lattice, check reference for details
        fractions = [Fraction(x).limit_denominator() for x in new_rot[:, kk]]
        least_mul = reduce(lcm, [fraction.denominator for fraction in fractions])
        scale = np.zeros((3, 3))
        scale[hh, hh] = 1
        scale[kk, kk] = least_mul
        scale[ll, ll] = sigma / least_mul
        n_final = None
        for idx in range(least_mul):
            check_int = idx * new_rot[:, kk] + (sigma / least_mul) * new_rot[:, ll]
            if all(np.round(x, 5).is_integer() for x in list(check_int)):
                n_final = idx
                break

        if n_final is None:
            raise RuntimeError("Something is wrong. Check if this GB exists or not")
        scale[kk, ll] = n_final
        # Each row of mat_csl is the CSL lattice vector
        csl_init = np.rint(np.dot(np.dot(r_matrix, trans), scale)).astype(int).T
        if abs(r_axis[hh]) > 1:
            csl_init = GrainBoundaryGenerator.reduce_mat(np.array(csl_init), r_axis[hh], r_matrix)
        csl = np.rint(Lattice(csl_init).get_niggli_reduced_lattice().matrix).astype(int)

        # Find the best slab supercell in terms of the conventional cell from the csl lattice,
        # which is the transformation matrix

        # Now trans_cry is the transformation matrix from crystal to Cartesian coordinates.
        # for cubic, do not need to change.
        if lat_type != "c":
            if lat_type == "h":
                trans_cry = np.array([[1, 0, 0], [-0.5, np.sqrt(3.0) / 2.0, 0], [0, 0, np.sqrt(mu / mv)]])  # type: ignore[operator]
            elif lat_type == "r":
                c2_a2_ratio = 1.0 if ratio is None else 3.0 / (2 - 6 * mv / mu)  # type: ignore[operator]
                trans_cry = np.array(
                    [
                        [0.5, np.sqrt(3.0) / 6.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)],
                        [-0.5, np.sqrt(3.0) / 6.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)],
                        [0, -1 * np.sqrt(3.0) / 3.0, 1.0 / 3 * np.sqrt(c2_a2_ratio)],
                    ]
                )
            else:
                trans_cry = np.array([[1, 0, 0], [0, np.sqrt(lam / mv), 0], [0, 0, np.sqrt(mu / mv)]])  # type: ignore[operator]
        t1_final = GrainBoundaryGenerator.slab_from_csl(
            csl, surface, normal, trans_cry, max_search=max_search, quick_gen=quick_gen
        )
        t2_final = np.array(np.rint(np.dot(t1_final, np.linalg.inv(r_matrix.T)))).astype(int)
        return t1_final, t2_final

    @staticmethod
    def enum_sigma_cubic(
        cutoff: int,
        r_axis: Tuple3Ints,
    ) -> dict[int, list[float]]:
        """Find all possible sigma values and corresponding rotation angles
        within a sigma value cutoff with known rotation axis in cubic system.
        The algorithm for this code is from reference, Acta Cryst, A40,108(1984).

        Args:
            cutoff (int): the cutoff of sigma values.
            r_axis ((u, v, w)): the rotation axis of the grain boundary.

        Returns:
            dict: sigmas dictionary with keys as the possible integer sigma values
                and values as list of the possible rotation angles to the
                corresponding sigma values. e.g. the format as
                {sigma1: [angle11,angle12,...], sigma2: [angle21, angle22,...],...}
                Note: the angles are the rotation angles of one grain respect to
                the other grain.
                When generating the microstructures of the grain boundary using these angles,
                you need to analyze the symmetry of the structure. Different angles may
                result in equivalent microstructures.
        """
        # Make sure math.gcd(r_axis) == 1
        if reduce(math.gcd, r_axis) != 1:
            r_axis = cast(Tuple3Ints, tuple(round(x / reduce(math.gcd, r_axis)) for x in r_axis))

        # Count the number of odds in r_axis
        odd_r = len(list(filter(lambda x: x % 2 == 1, r_axis)))
        # Compute the max n we need to enumerate
        if odd_r == 3:
            a_max = 4
        elif odd_r == 0:
            a_max = 1
        else:
            a_max = 2
        n_max = int(np.sqrt(cutoff * a_max / sum(np.array(r_axis) ** 2)))
        # Enumerate all possible n, m to give possible sigmas within the cutoff
        sigmas: dict[int, list[float]] = {}
        for n_loop in range(1, n_max + 1):
            n = n_loop
            m_max = int(np.sqrt(cutoff * a_max - n**2 * sum(np.array(r_axis) ** 2)))
            for m in range(m_max + 1):
                if math.gcd(m, n) == 1 or m == 0:
                    n = 1 if m == 0 else n_loop
                    # Construct the quadruple [m, U,V,W], count the number of odds in
                    # quadruple to determine the parameter a, refer to the reference
                    quadruple = [m] + [x * n for x in r_axis]
                    odd_qua = len(list(filter(lambda x: x % 2 == 1, quadruple)))
                    if odd_qua == 4:
                        a = 4
                    elif odd_qua == 2:
                        a = 2
                    else:
                        a = 1
                    sigma = round((m**2 + n**2 * sum(np.array(r_axis) ** 2)) / a)
                    if 1 < sigma <= cutoff:
                        if sigma not in list(sigmas):
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n * np.sqrt(sum(np.array(r_axis) ** 2)) / m) / np.pi * 180
                            sigmas[sigma] = [angle]
                        else:
                            if m == 0:
                                angle = 180.0
                            else:
                                angle = 2 * np.arctan(n * np.sqrt(sum(np.array(r_axis) ** 2)) / m) / np.pi * 180
                            if angle not in sigmas[sigma]:
                                sigmas[sigma].append(angle)
        return sigmas

    @staticmethod
    def enum_sigma_hex(
        cutoff: int,
        r_axis: Tuple3Ints | Tuple4Ints,
        c2_a2_ratio: tuple[int, int],
    ) -> dict[int, list[float]]:
        """Find all possible sigma values and corresponding rotation angles
        within a sigma value cutoff with known rotation axis in hexagonal system.
        The algorithm for this code is from reference, Acta Cryst, A38,550(1982).

        Args:
            cutoff (int): the cutoff of sigma values.
            r_axis ((u, v, w) or (u, v, t, w)): the rotation axis of the grain boundary.
            c2_a2_ratio ((mu, mv)): mu/mv is the square of the hexagonal axial ratio,
                which is rational number. If irrational, set c2_a2_ratio = None

        Returns:
            dict: sigmas dictionary with keys as the possible integer sigma values
                and values as list of the possible rotation angles to the
                corresponding sigma values. e.g. the format as
                {sigma1: [angle11, angle12, ...], sigma2: [angle21, angle22, ...], ...}
                Note: the angles are the rotation angles of one grain respect to
                the other grain.
                When generating the microstructures of the grain boundary using these angles,
                you need to analyze the symmetry of the structure. Different angles may
                result in equivalent microstructures.
        """
        # Make sure math.gcd(r_axis) == 1
        if reduce(math.gcd, r_axis) != 1:
            r_axis = cast(
                Tuple3Ints | Tuple4Ints,
                tuple(round(x / reduce(math.gcd, r_axis)) for x in r_axis),
            )

        # Transform four index notation to three index notation
        if len(r_axis) == 4:
            u1 = r_axis[0]
            v1 = r_axis[1]
            w1 = r_axis[3]
            u = 2 * u1 + v1
            v = 2 * v1 + u1
            w = w1
        else:
            u, v, w = r_axis  # type: ignore[misc]

        # Make sure mu, mv are coprime integers
        if c2_a2_ratio is None:
            mu, mv = [1, 1]
            if w != 0 and (u != 0 or (v != 0)):
                raise RuntimeError("For irrational c2/a2, CSL only exist for [0,0,1] or [u,v,0] and m = 0")
        else:
            mu, mv = c2_a2_ratio
            if math.gcd(mu, mv) != 1:
                temp = math.gcd(mu, mv)
                mu = round(mu / temp)
                mv = round(mv / temp)

        # Refer to the meaning of d in reference
        d = (u**2 + v**2 - u * v) * mv + w**2 * mu

        # Compute the max n we need to enumerate
        n_max = int(np.sqrt((cutoff * 12 * mu * mv) / abs(d)))

        # Enumerate all possible n, m to give possible sigmas within the cutoff
        sigmas: dict[int, list[float]] = {}
        for n in range(1, n_max + 1):
            if (c2_a2_ratio is None) and w == 0:
                m_max = 0
            else:
                m_max = int(np.sqrt((cutoff * 12 * mu * mv - n**2 * d) / (3 * mu)))
            for m in range(m_max + 1):
                if math.gcd(m, n) == 1 or m == 0:
                    # Construct the rotation matrix, refer to the reference
                    R_list = [
                        (u**2 * mv - v**2 * mv - w**2 * mu) * n**2 + 2 * w * mu * m * n + 3 * mu * m**2,
                        (2 * v - u) * u * mv * n**2 - 4 * w * mu * m * n,
                        2 * u * w * mu * n**2 + 2 * (2 * v - u) * mu * m * n,
                        (2 * u - v) * v * mv * n**2 + 4 * w * mu * m * n,
                        (v**2 * mv - u**2 * mv - w**2 * mu) * n**2 - 2 * w * mu * m * n + 3 * mu * m**2,
                        2 * v * w * mu * n**2 - 2 * (2 * u - v) * mu * m * n,
                        (2 * u - v) * w * mv * n**2 - 3 * v * mv * m * n,
                        (2 * v - u) * w * mv * n**2 + 3 * u * mv * m * n,
                        (w**2 * mu - u**2 * mv - v**2 * mv + u * v * mv) * n**2 + 3 * mu * m**2,
                    ]
                    m = -1 * m
                    # Inverse of the rotation matrix
                    R_list_inv = [
                        (u**2 * mv - v**2 * mv - w**2 * mu) * n**2 + 2 * w * mu * m * n + 3 * mu * m**2,
                        (2 * v - u) * u * mv * n**2 - 4 * w * mu * m * n,
                        2 * u * w * mu * n**2 + 2 * (2 * v - u) * mu * m * n,
                        (2 * u - v) * v * mv * n**2 + 4 * w * mu * m * n,
                        (v**2 * mv - u**2 * mv - w**2 * mu) * n**2 - 2 * w * mu * m * n + 3 * mu * m**2,
                        2 * v * w * mu * n**2 - 2 * (2 * u - v) * mu * m * n,
                        (2 * u - v) * w * mv * n**2 - 3 * v * mv * m * n,
                        (2 * v - u) * w * mv * n**2 + 3 * u * mv * m * n,
                        (w**2 * mu - u**2 * mv - v**2 * mv + u * v * mv) * n**2 + 3 * mu * m**2,
                    ]
                    m = -1 * m
                    F = 3 * mu * m**2 + d * n**2
                    all_list = R_list_inv + R_list + [F]
                    # Compute the max common factors for the elements of the rotation matrix
                    # and its inverse.
                    com_fac = reduce(math.gcd, all_list)
                    sigma = round((3 * mu * m**2 + d * n**2) / com_fac)
                    if 1 < sigma <= cutoff:
                        if sigma not in list(sigmas):
                            angle = 180.0 if m == 0 else 2 * np.arctan(n / m * np.sqrt(d / 3.0 / mu)) / np.pi * 180
                            sigmas[sigma] = [angle]
                        else:
                            angle = 180.0 if m == 0 else 2 * np.arctan(n / m * np.sqrt(d / 3.0 / mu)) / np.pi * 180
                            if angle not in sigmas[sigma]:
                                sigmas[sigma].append(angle)
            if m_max == 0:
                break
        return sigmas

    @staticmethod
    def enum_sigma_rho(
        cutoff: int,
        r_axis: Tuple3Ints | Tuple4Ints,
        ratio_alpha: tuple[int, int],
    ) -> dict[int, list[float]]:
        """Find all possible sigma values and corresponding rotation angles
        within a sigma value cutoff with known rotation axis in rhombohedral system.
        The algorithm for this code is from reference, Acta Cryst, A45,505(1989).

        Args:
            cutoff (int): the cutoff of sigma values.
            r_axis ((u, v, w) or (u, v, t, w)): the rotation axis of the grain boundary.
            ratio_alpha (tuple of two integers, e.g. mu, mv):
                    mu/mv is the ratio of (1+2*cos(alpha))/cos(alpha) with rational number.
                    If irrational, set ratio_alpha = None.

        Returns:
            dict[int, list[float]]: keys are possible integer sigma values
                and values are lists of possible rotation angles to the
                {sigma1: [angle11, angle12,...], sigma2: [angle21, angle22,...],...}
                Note: the angles are the rotation angle of one grain respect to the
                other grain.
                When generating the microstructure of the grain boundary using these
                angles, you need to analyze the symmetry of the structure. Different
                angles may result in equivalent microstructures.
        """
        # Transform four index notation to three index notation
        if len(r_axis) == 4:
            u1 = r_axis[0]
            v1 = r_axis[1]
            w1 = r_axis[3]
            u = 2 * u1 + v1 + w1
            v = v1 + w1 - u1
            w = w1 - 2 * v1 - u1
            r_axis = (u, v, w)

        # Make sure math.(r_axis) == 1
        if reduce(math.gcd, r_axis) != 1:
            r_axis = cast(Tuple3Ints, tuple(round(x / reduce(math.gcd, r_axis)) for x in r_axis))
        u, v, w = r_axis  # type: ignore[misc]

        # Make sure mu, mv are coprime integers
        if ratio_alpha is None:
            mu, mv = [1, 1]
            if u + v + w != 0 and (u != v or u != w):
                raise RuntimeError("For irrational ratio_alpha, CSL only exist for [1,1,1] or [u, v, -(u+v)] and m =0")
        else:
            mu, mv = ratio_alpha
            if math.gcd(mu, mv) != 1:
                temp = math.gcd(mu, mv)
                mu = round(mu / temp)
                mv = round(mv / temp)

        # Refer to the meaning of d in reference
        d = (u**2 + v**2 + w**2) * (mu - 2 * mv) + 2 * mv * (v * w + w * u + u * v)
        # Compute the max n we need to enumerate
        n_max = int(np.sqrt((cutoff * abs(4 * mu * (mu - 3 * mv))) / abs(d)))

        # Enumerate all possible n, m to give possible sigmas within the cutoff
        sigmas: dict[int, list[float]] = {}
        for n in range(1, n_max + 1):
            if ratio_alpha is None and u + v + w == 0:
                m_max = 0
            else:
                m_max = int(np.sqrt((cutoff * abs(4 * mu * (mu - 3 * mv)) - n**2 * d) / (mu)))
            for m in range(m_max + 1):
                if math.gcd(m, n) == 1 or m == 0:
                    # Construct the rotation matrix, refer to the reference
                    R_list = [
                        (mu - 2 * mv) * (u**2 - v**2 - w**2) * n**2
                        + 2 * mv * (v - w) * m * n
                        - 2 * mv * v * w * n**2
                        + mu * m**2,
                        2 * (mv * u * n * (w * n + u * n - m) - (mu - mv) * m * w * n + (mu - 2 * mv) * u * v * n**2),
                        2 * (mv * u * n * (v * n + u * n + m) + (mu - mv) * m * v * n + (mu - 2 * mv) * w * u * n**2),
                        2 * (mv * v * n * (w * n + v * n + m) + (mu - mv) * m * w * n + (mu - 2 * mv) * u * v * n**2),
                        (mu - 2 * mv) * (v**2 - w**2 - u**2) * n**2
                        + 2 * mv * (w - u) * m * n
                        - 2 * mv * u * w * n**2
                        + mu * m**2,
                        2 * (mv * v * n * (v * n + u * n - m) - (mu - mv) * m * u * n + (mu - 2 * mv) * w * v * n**2),
                        2 * (mv * w * n * (w * n + v * n - m) - (mu - mv) * m * v * n + (mu - 2 * mv) * w * u * n**2),
                        2 * (mv * w * n * (w * n + u * n + m) + (mu - mv) * m * u * n + (mu - 2 * mv) * w * v * n**2),
                        (mu - 2 * mv) * (w**2 - u**2 - v**2) * n**2
                        + 2 * mv * (u - v) * m * n
                        - 2 * mv * u * v * n**2
                        + mu * m**2,
                    ]
                    m = -1 * m
                    # Inverse of the rotation matrix
                    R_list_inv = [
                        (mu - 2 * mv) * (u**2 - v**2 - w**2) * n**2
                        + 2 * mv * (v - w) * m * n
                        - 2 * mv * v * w * n**2
                        + mu * m**2,
                        2 * (mv * u * n * (w * n + u * n - m) - (mu - mv) * m * w * n + (mu - 2 * mv) * u * v * n**2),
                        2 * (mv * u * n * (v * n + u * n + m) + (mu - mv) * m * v * n + (mu - 2 * mv) * w * u * n**2),
                        2 * (mv * v * n * (w * n + v * n + m) + (mu - mv) * m * w * n + (mu - 2 * mv) * u * v * n**2),
                        (mu - 2 * mv) * (v**2 - w**2 - u**2) * n**2
                        + 2 * mv * (w - u) * m * n
                        - 2 * mv * u * w * n**2
                        + mu * m**2,
                        2 * (mv * v * n * (v * n + u * n - m) - (mu - mv) * m * u * n + (mu - 2 * mv) * w * v * n**2),
                        2 * (mv * w * n * (w * n + v * n - m) - (mu - mv) * m * v * n + (mu - 2 * mv) * w * u * n**2),
                        2 * (mv * w * n * (w * n + u * n + m) + (mu - mv) * m * u * n + (mu - 2 * mv) * w * v * n**2),
                        (mu - 2 * mv) * (w**2 - u**2 - v**2) * n**2
                        + 2 * mv * (u - v) * m * n
                        - 2 * mv * u * v * n**2
                        + mu * m**2,
                    ]
                    m = -1 * m
                    F = mu * m**2 + d * n**2
                    all_list = R_list_inv + R_list + [F]
                    # Compute the max common factors for the elements of the rotation matrix and its inverse.
                    com_fac = reduce(math.gcd, all_list)
                    sigma = round(abs(F / com_fac))
                    if 1 < sigma <= cutoff:
                        if sigma not in list(sigmas):
                            angle = 180.0 if m == 0 else 2 * np.arctan(n / m * np.sqrt(d / mu)) / np.pi * 180
                            sigmas[sigma] = [angle]
                        else:
                            angle = 180 if m == 0 else 2 * np.arctan(n / m * np.sqrt(d / mu)) / np.pi * 180.0
                            if angle not in sigmas[sigma]:
                                sigmas[sigma].append(angle)
            if m_max == 0:
                break
        return sigmas

    @staticmethod
    def enum_sigma_tet(
        cutoff: int,
        r_axis: Tuple3Ints,
        c2_a2_ratio: tuple[int, int],
    ) -> dict[int, list[float]]:
        """Find all possible sigma values and corresponding rotation angles
        within a sigma value cutoff with known rotation axis in tetragonal system.
        The algorithm for this code is from reference, Acta Cryst, B46,117(1990).

        Args:
            cutoff (int): the cutoff of sigma values.
            r_axis ((u, v, w)): the rotation axis of the grain boundary.
            c2_a2_ratio ((mu, mv)): mu/mv is the square of the tetragonal axial
                ratio with rational number. If irrational, set c2_a2_ratio = None.

        Returns:
            dict: sigmas dictionary with keys as the possible integer sigma values
                and values as list of the possible rotation angles to the
                corresponding sigma values. e.g. the format as
                {sigma1: [angle11, angle12, ...], sigma2: [angle21, angle22, ...], ...}
                Note: the angles are the rotation angle of one grain respect to the
                other grain.
                When generating the microstructure of the grain boundary using these
                angles, you need to analyze the symmetry of the structure. Different
                angles may result in equivalent microstructures.
        """
        # Make sure math.gcd(r_axis) == 1
        if reduce(math.gcd, r_axis) != 1:
            r_axis = cast(Tuple3Ints, tuple(round(x / reduce(math.gcd, r_axis)) for x in r_axis))

        u, v, w = r_axis

        # Make sure mu, mv are coprime integers
        if c2_a2_ratio is None:
            mu, mv = [1, 1]
            if w != 0 and (u != 0 or (v != 0)):
                raise RuntimeError("For irrational c2/a2, CSL only exist for [0,0,1] or [u,v,0] and m = 0")
        else:
            mu, mv = c2_a2_ratio
            if math.gcd(mu, mv) != 1:
                temp = math.gcd(mu, mv)
                mu = round(mu / temp)
                mv = round(mv / temp)

        # Refer to the meaning of d in reference
        d = (u**2 + v**2) * mv + w**2 * mu

        # Compute the max n we need to enumerate
        n_max = int(np.sqrt((cutoff * 4 * mu * mv) / d))

        # Enumerate all possible n, m to give possible sigmas within the cutoff
        sigmas: dict[int, list[float]] = {}
        for n in range(1, n_max + 1):
            m_max = 0 if c2_a2_ratio is None and w == 0 else int(np.sqrt((cutoff * 4 * mu * mv - n**2 * d) / mu))
            for m in range(m_max + 1):
                if math.gcd(m, n) == 1 or m == 0:
                    # Construct the rotation matrix, refer to the reference
                    R_list = [
                        (u**2 * mv - v**2 * mv - w**2 * mu) * n**2 + mu * m**2,
                        2 * v * u * mv * n**2 - 2 * w * mu * m * n,
                        2 * u * w * mu * n**2 + 2 * v * mu * m * n,
                        2 * u * v * mv * n**2 + 2 * w * mu * m * n,
                        (v**2 * mv - u**2 * mv - w**2 * mu) * n**2 + mu * m**2,
                        2 * v * w * mu * n**2 - 2 * u * mu * m * n,
                        2 * u * w * mv * n**2 - 2 * v * mv * m * n,
                        2 * v * w * mv * n**2 + 2 * u * mv * m * n,
                        (w**2 * mu - u**2 * mv - v**2 * mv) * n**2 + mu * m**2,
                    ]
                    m = -1 * m
                    # Inverse of rotation matrix
                    R_list_inv = [
                        (u**2 * mv - v**2 * mv - w**2 * mu) * n**2 + mu * m**2,
                        2 * v * u * mv * n**2 - 2 * w * mu * m * n,
                        2 * u * w * mu * n**2 + 2 * v * mu * m * n,
                        2 * u * v * mv * n**2 + 2 * w * mu * m * n,
                        (v**2 * mv - u**2 * mv - w**2 * mu) * n**2 + mu * m**2,
                        2 * v * w * mu * n**2 - 2 * u * mu * m * n,
                        2 * u * w * mv * n**2 - 2 * v * mv * m * n,
                        2 * v * w * mv * n**2 + 2 * u * mv * m * n,
                        (w**2 * mu - u**2 * mv - v**2 * mv) * n**2 + mu * m**2,
                    ]
                    m = -1 * m
                    F = mu * m**2 + d * n**2
                    all_list = R_list + R_list_inv + [F]
                    # Compute the max common factors for the elements of the rotation matrix
                    # and its inverse
                    com_fac = reduce(math.gcd, all_list)
                    sigma = round((mu * m**2 + d * n**2) / com_fac)
                    if 1 < sigma <= cutoff:
                        if sigma not in list(sigmas):
                            angle = 180.0 if m == 0 else 2 * np.arctan(n / m * np.sqrt(d / mu)) / np.pi * 180
                            sigmas[sigma] = [angle]
                        else:
                            angle = 180.0 if m == 0 else 2 * np.arctan(n / m * np.sqrt(d / mu)) / np.pi * 180
                            if angle not in sigmas[sigma]:
                                sigmas[sigma].append(angle)
            if m_max == 0:
                break

        return sigmas

    @staticmethod
    def enum_sigma_ort(
        cutoff: int,
        r_axis: Tuple3Ints,
        c2_b2_a2_ratio: Tuple3Floats,
    ) -> dict[int, list[float]]:
        """Find all possible sigma values and corresponding rotation angles
        within a sigma value cutoff with known rotation axis in orthorhombic system.

        Reference: Scipta Metallurgica 27, 291(1992).

        Args:
            cutoff (int): the cutoff of sigma values.
            r_axis ((u, v, w)): the rotation axis of the grain boundary.
            c2_b2_a2_ratio ((mu, lambda, mv)):
                mu:lam:mv is the square of the orthorhombic axial ratio with rational
                numbers. If irrational for one axis, set it to None.
                e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.

        Returns:
            dict: sigmas dictionary with keys as the possible integer sigma values
                and values as list of the possible rotation angles to the
                corresponding sigma values. e.g. the format as
                {sigma1: [angle11,angle12,...], sigma2: [angle21, angle22,...],...}
                Note: the angles are the rotation angle of one grain respect to the
                other grain.
                When generating the microstructure of the grain boundary using these
                angles, you need to analyze the symmetry of the structure. Different
                angles may result in equivalent microstructures.
        """
        # Make sure math.gcd(r_axis) == 1
        if reduce(math.gcd, r_axis) != 1:
            r_axis = cast(Tuple3Ints, tuple(round(x / reduce(math.gcd, r_axis)) for x in r_axis))

        u, v, w = r_axis

        # Make sure mu, lambda, mv are coprime integers
        if None in c2_b2_a2_ratio:
            mu, lam, mv = c2_b2_a2_ratio
            non_none = [i for i in c2_b2_a2_ratio if i is not None]
            if len(non_none) < 2:
                raise RuntimeError("No CSL exist for two irrational numbers")
            non1, non2 = non_none
            if reduce(math.gcd, non_none) != 1:  # type: ignore[arg-type]
                temp = reduce(math.gcd, non_none)  # type: ignore[arg-type]
                non1 = round(non1 / temp)
                non2 = round(non2 / temp)
            if mu is None:
                lam = non1
                mv = non2
                mu = 1
                if w != 0 and (u != 0 or (v != 0)):
                    raise RuntimeError("For irrational c2, CSL only exist for [0,0,1] or [u,v,0] and m = 0")
            elif lam is None:
                mu = non1
                mv = non2
                lam = 1
                if v != 0 and (u != 0 or (w != 0)):
                    raise RuntimeError("For irrational b2, CSL only exist for [0,1,0] or [u,0,w] and m = 0")
            elif mv is None:
                mu = non1
                lam = non2
                mv = 1
                if u != 0 and (w != 0 or (v != 0)):
                    raise RuntimeError("For irrational a2, CSL only exist for [1,0,0] or [0,v,w] and m = 0")
        else:
            mu, lam, mv = c2_b2_a2_ratio
            if reduce(math.gcd, c2_b2_a2_ratio) != 1:  # type: ignore[arg-type]
                temp = reduce(math.gcd, c2_b2_a2_ratio)  # type: ignore[arg-type]
                mu = round(mu / temp)
                mv = round(mv / temp)
                lam = round(lam / temp)
            if u == 0 and v == 0:
                mu = 1
            if u == 0 and w == 0:
                lam = 1
            if v == 0 and w == 0:
                mv = 1
        # Refer to the reference for the meaning of d
        d = (mv * u**2 + lam * v**2) * mv + w**2 * mu * mv

        # Compute the max n we need to enumerate
        n_max = int(np.sqrt((cutoff * 4 * mu * mv * mv * lam) / d))
        # Enumerate all possible n, m to give possible sigmas within the cutoff
        sigmas: dict[int, list[float]] = {}
        for n in range(1, n_max + 1):
            mu_temp, lam_temp, mv_temp = c2_b2_a2_ratio
            if (mu_temp is None and w == 0) or (lam_temp is None and v == 0) or (mv_temp is None and u == 0):
                m_max = 0
            else:
                m_max = int(np.sqrt((cutoff * 4 * mu * mv * lam * mv - n**2 * d) / mu / lam))
            for m in range(m_max + 1):
                if math.gcd(m, n) == 1 or m == 0:
                    # Construct the rotation matrix, refer to the reference
                    R_list = [
                        (u**2 * mv * mv - lam * v**2 * mv - w**2 * mu * mv) * n**2 + lam * mu * m**2,
                        2 * lam * (v * u * mv * n**2 - w * mu * m * n),
                        2 * mu * (u * w * mv * n**2 + v * lam * m * n),
                        2 * mv * (u * v * mv * n**2 + w * mu * m * n),
                        (v**2 * mv * lam - u**2 * mv * mv - w**2 * mu * mv) * n**2 + lam * mu * m**2,
                        2 * mv * mu * (v * w * n**2 - u * m * n),
                        2 * mv * (u * w * mv * n**2 - v * lam * m * n),
                        2 * lam * mv * (v * w * n**2 + u * m * n),
                        (w**2 * mu * mv - u**2 * mv * mv - v**2 * mv * lam) * n**2 + lam * mu * m**2,
                    ]
                    m = -1 * m
                    # Inverse of rotation matrix
                    R_list_inv = [
                        (u**2 * mv * mv - lam * v**2 * mv - w**2 * mu * mv) * n**2 + lam * mu * m**2,
                        2 * lam * (v * u * mv * n**2 - w * mu * m * n),
                        2 * mu * (u * w * mv * n**2 + v * lam * m * n),
                        2 * mv * (u * v * mv * n**2 + w * mu * m * n),
                        (v**2 * mv * lam - u**2 * mv * mv - w**2 * mu * mv) * n**2 + lam * mu * m**2,
                        2 * mv * mu * (v * w * n**2 - u * m * n),
                        2 * mv * (u * w * mv * n**2 - v * lam * m * n),
                        2 * lam * mv * (v * w * n**2 + u * m * n),
                        (w**2 * mu * mv - u**2 * mv * mv - v**2 * mv * lam) * n**2 + lam * mu * m**2,
                    ]
                    m = -1 * m
                    F = mu * lam * m**2 + d * n**2
                    all_list = R_list + R_list_inv + [F]
                    # Compute the max common factors for the elements of the rotation matrix
                    # and its inverse.
                    com_fac = reduce(math.gcd, all_list)  # type: ignore[arg-type]
                    sigma = round((mu * lam * m**2 + d * n**2) / com_fac)
                    if 1 < sigma <= cutoff:
                        if sigma not in list(sigmas):
                            angle = 180.0 if m == 0 else 2 * np.arctan(n / m * np.sqrt(d / mu / lam)) / np.pi * 180
                            sigmas[sigma] = [angle]
                        else:
                            angle = 180.0 if m == 0 else 2 * np.arctan(n / m * np.sqrt(d / mu / lam)) / np.pi * 180
                            if angle not in sigmas[sigma]:
                                sigmas[sigma].append(angle)
            if m_max == 0:
                break

        return sigmas

    @staticmethod
    def enum_possible_plane_cubic(
        plane_cutoff: int,
        r_axis: Tuple3Ints,
        r_angle: float,
    ) -> dict[Literal["Twist", "Symmetric tilt", "Normal tilt", "Mixed"], list[list]]:
        """Find all possible plane combinations for GBs given a rotation axis and angle for
        cubic system, and classify them to different categories, including "Twist",
        "Symmetric tilt", "Normal tilt", "Mixed" GBs.

        Args:
            plane_cutoff (int): the cutoff of plane miller index.
            r_axis ((u, v, w])): the rotation axis of the grain boundary.
            r_angle (float): rotation angle of the GBs.

        Returns:
            dict: all combinations with keys as GB type, e.g. "Twist", "Symmetric tilt", etc.
                and values as the combination of the two plane miller index (GB plane and joining plane).
        """
        all_combinations: dict[Literal["Twist", "Symmetric tilt", "Normal tilt", "Mixed"], list[list]] = {
            "Symmetric tilt": [],
            "Twist": [],
            "Normal tilt": [],
            "Mixed": [],
        }
        sym_plane = symm_group_cubic([[1, 0, 0], [1, 1, 0]])
        j = np.arange(0, plane_cutoff + 1)
        combination = []
        for idx in product(j, repeat=3):
            if sum(abs(np.array(idx))) != 0:
                combination.append(list(idx))
            if len(np.nonzero(idx)[0]) == 3:
                for i1 in range(3):
                    new_i = list(idx).copy()
                    new_i[i1] = -1 * new_i[i1]
                    combination.append(new_i)
            elif len(np.nonzero(idx)[0]) == 2:
                new_i = list(idx).copy()
                new_i[np.nonzero(idx)[0][0]] = -1 * new_i[np.nonzero(idx)[0][0]]
                combination.append(new_i)
        miller = np.array(combination)
        miller = miller[np.argsort(np.linalg.norm(miller, axis=1))]
        for val in miller:
            if reduce(math.gcd, val) == 1:
                matrix = GrainBoundaryGenerator.get_trans_mat(r_axis, r_angle, surface=val, quick_gen=True)
                vec = np.cross(matrix[1][0], matrix[1][1])
                miller2 = GrainBoundaryGenerator.vec_to_surface(vec)
                if np.all(np.abs(np.array(miller2)) <= plane_cutoff):
                    cos_1 = abs(np.dot(val, r_axis) / np.linalg.norm(val) / np.linalg.norm(r_axis))
                    if 1 - cos_1 < 1.0e-5:
                        all_combinations["Twist"].append([list(val), miller2])
                    elif cos_1 < 1.0e-8:
                        sym_tilt = False
                        if np.sum(np.abs(val)) == np.sum(np.abs(miller2)):
                            ave = (np.array(val) + np.array(miller2)) / 2
                            ave1 = (np.array(val) - np.array(miller2)) / 2
                            for plane in sym_plane:
                                cos_2 = abs(np.dot(ave, plane) / np.linalg.norm(ave) / np.linalg.norm(plane))
                                cos_3 = abs(np.dot(ave1, plane) / np.linalg.norm(ave1) / np.linalg.norm(plane))
                                if 1 - cos_2 < 1.0e-5 or 1 - cos_3 < 1.0e-5:
                                    all_combinations["Symmetric tilt"].append([list(val), miller2])
                                    sym_tilt = True
                                    break
                        if not sym_tilt:
                            all_combinations["Normal tilt"].append([list(val), miller2])
                    else:
                        all_combinations["Mixed"].append([list(val), miller2])
        return all_combinations

    @staticmethod
    def get_rotation_angle_from_sigma(
        sigma: int,
        r_axis: Tuple3Ints | Tuple4Ints,
        lat_type: str = "c",
        ratio: tuple[int, int] | Tuple3Ints | None = None,
    ) -> list[float]:
        """Find all possible rotation angles for the given sigma value.

        Args:
            sigma (int): sigma value provided
            r_axis ((u, v, w) or (u, v, t, w) for hex/rho systems): the
                rotation axis of the grain boundary.
            lat_type (str): one character to specify the lattice type. Defaults to 'c' for cubic.
                'c' or 'C': cubic system
                't' or 'T': tetragonal system
                'o' or 'O': orthorhombic system
                'h' or 'H': hexagonal system
                'r' or 'R': rhombohedral system
            ratio (tuple[int, ...]): lattice axial ratio.
                For cubic system, ratio is not needed.
                For tetragonal system, ratio = (mu, mv), tuple of two integers,
                that is, mu/mv = c2/a2. If it is irrational, set it to none.
                For orthorhombic system, ratio = (mu, lam, mv), tuple of 3 integers,
                that is, mu:lam:mv = c2:b2:a2. If irrational for one axis, set it to None.
                e.g. mu:lam:mv = c2,None,a2, means b2 is irrational.
                For rhombohedral system, ratio = (mu, mv), tuple of two integers,
                that is, mu/mv is the ratio of (1+2*cos(alpha)/cos(alpha).
                If irrational, set it to None.
                For hexagonal system, ratio = (mu, mv), tuple of two integers,
                that is, mu/mv = c2/a2. If it is irrational, set it to none.

        Returns:
            rotation_angles corresponding to the provided sigma value.
            If the sigma value is not correct, return the rotation angle corresponding
            to the correct possible sigma value right smaller than the wrong sigma value provided.
        """
        lat_type = lat_type.lower()

        # Check r_axis length
        if lat_type in {"c", "t"} and len(r_axis) != 3:
            raise ValueError(f"expect r_axis length 3 for selected lattice system, got {len(r_axis)}")

        # Check lattice axial ratio length
        if lat_type == "o" and (ratio is None or len(ratio) != 3):
            raise RuntimeError("Orthorhombic system needs correct c2:b2:a2 ratio")

        if lat_type in {"t", "h"} and (ratio is None or len(ratio) != 2):
            raise RuntimeError("System needs correct c2/a2 ratio.")

        if lat_type == "r" and (ratio is None or len(ratio) != 2):
            raise RuntimeError("Rhombohedral system needs correct (1+2*cos(alpha)/cos(alpha) ratio")

        if lat_type == "c":
            logger.info("Make sure this is for cubic system")
            sigma_dict = GrainBoundaryGenerator.enum_sigma_cubic(cutoff=sigma, r_axis=cast(Tuple3Ints, r_axis))

        elif lat_type == "t":
            logger.info("Make sure this is for tetragonal system")
            if ratio is None:
                logger.info("Make sure this is for irrational c2/a2 ratio")
            sigma_dict = GrainBoundaryGenerator.enum_sigma_tet(
                cutoff=sigma,
                r_axis=cast(Tuple3Ints, r_axis),
                c2_a2_ratio=cast(tuple[int, int], ratio),
            )

        elif lat_type == "o":
            logger.info("Make sure this is for orthorhombic system")
            sigma_dict = GrainBoundaryGenerator.enum_sigma_ort(
                cutoff=sigma,
                r_axis=cast(Tuple3Ints, r_axis),
                c2_b2_a2_ratio=cast(Tuple3Ints, ratio),
            )

        elif lat_type == "h":
            logger.info("Make sure this is for hexagonal system")
            if ratio is None:
                logger.info("Make sure this is for irrational c2/a2 ratio")
            sigma_dict = GrainBoundaryGenerator.enum_sigma_hex(
                cutoff=sigma, r_axis=r_axis, c2_a2_ratio=cast(tuple[int, int], ratio)
            )

        elif lat_type == "r":
            logger.info("Make sure this is for rhombohedral system")
            if ratio is None:
                logger.info("Make sure this is for irrational (1+2*cos(alpha)/cos(alpha) ratio")
            sigma_dict = GrainBoundaryGenerator.enum_sigma_rho(
                cutoff=sigma,
                r_axis=cast(Tuple3Ints, r_axis),
                ratio_alpha=cast(tuple[int, int], ratio),
            )

        else:
            raise RuntimeError("Lattice type not implemented")

        sigmas = list(sigma_dict)
        if not sigmas:
            raise RuntimeError("This is a wrong sigma value, and no sigma exists smaller than this value.")
        if sigma in sigmas:
            rotation_angles = sigma_dict[sigma]
        else:
            sigmas.sort()
            warnings.warn(
                "This is not the possible sigma value according to the rotation axis!"
                "The nearest neighbor sigma and its corresponding angle are returned"
            )
            rotation_angles = sigma_dict[sigmas[-1]]
        rotation_angles.sort()
        return rotation_angles

    @staticmethod
    def slab_from_csl(
        csl: NDArray,
        surface: Tuple3Ints,
        normal: bool,
        trans_cry: NDArray,
        max_search: int = 20,
        quick_gen: bool = False,
    ) -> Matrix3D:
        """By linear operation of csl lattice vectors to get the best corresponding
        slab lattice. That is the area of a,b vectors (within the surface plane)
        is the smallest, the c vector first, has shortest length perpendicular
        to surface [h,k,l], second, has shortest length itself.

        Args:
            csl (3 by 3 integer array): input csl lattice.
            surface ((h, k, l)): the miller index of the surface.
                    determine if the c vector needs to perpendicular to surface.
            normal (bool): search the GB c vector that is perpendicular to the plane.
            trans_cry (3 by 3 array): transform matrix from crystal system to orthogonal system.
            max_search (int): max search for the GB lattice vectors that give the smallest GB
                lattice. If normal is True, also max search the GB c vector that is perpendicular
                to the plane.
            quick_gen (bool): whether to quickly generate a supercell, no need to find the smallest
                cell if set to true.

        Returns:
            t_matrix: a slab lattice (3 by 3 integer array)
        """
        # Set the transform matrix in real space
        trans = trans_cry
        # Transform matrix in reciprocal space
        ctrans = np.linalg.inv(trans.T)

        t_matrix = csl.copy()
        # Vectors constructed from csl that perpendicular to surface
        ab_vector = []
        # Obtain the miller index of surface in terms of csl
        miller = np.matmul(surface, csl.T)
        if reduce(math.gcd, miller) != 1:
            miller = [round(x / reduce(math.gcd, miller)) for x in miller]
        miller_nonzero = []
        # Quickly generate a supercell, normal is not working in this way
        if quick_gen:
            scale_factor = []
            eye = np.eye(3, dtype=np.int64)
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
                        scl_factor = [0, 0, 0]
                        scl_factor[miller_nonzero[i]] = -round(lcm_miller / miller[miller_nonzero[i]])
                        scl_factor[miller_nonzero[j]] = round(lcm_miller / miller[miller_nonzero[j]])
                        scale_factor.append(scl_factor)
                        if len(scale_factor) == 2:
                            break
            t_matrix[0] = np.array(np.dot(scale_factor[0], csl))
            t_matrix[1] = np.array(np.dot(scale_factor[1], csl))
            t_matrix[2] = csl[miller_nonzero[0]]
            if abs(np.linalg.det(t_matrix)) > 1000:
                warnings.warn("Too large matrix. Suggest to use quick_gen=False")
            return t_matrix

        c_index = 0
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
                    com_gcd = math.gcd(miller_nonzero[i], miller_nonzero[j])
                    mil1 = round(miller_nonzero[i] / com_gcd)
                    mil2 = round(miller_nonzero[j] / com_gcd)
                    lcm_miller.append(max(abs(mil1), abs(mil2)))
            lcm_sorted = sorted(lcm_miller)
            max_j = lcm_sorted[0] if index_len == 2 else lcm_sorted[1]
        else:
            if not normal:
                t_matrix[0] = ab_vector[0]
                t_matrix[1] = ab_vector[1]
                t_matrix[2] = csl[c_index]
                return t_matrix
            max_j = abs(miller_nonzero[0])
        max_j = min(max_j, max_search)
        # Area of a, b vectors
        area = None
        # Length of c vector
        c_norm = np.linalg.norm(np.matmul(t_matrix[2], trans))
        # c vector length along the direction perpendicular to surface
        c_length = np.abs(np.dot(t_matrix[2], surface))
        # Check if the init c vector perpendicular to the surface
        if normal:
            c_cross = np.cross(np.matmul(t_matrix[2], trans), np.matmul(surface, ctrans))
            normal_init = np.linalg.norm(c_cross) < 1e-8
        else:
            normal_init = False

        jj = np.arange(0, max_j + 1)
        combination = []
        for ii in product(jj, repeat=3):
            if sum(abs(np.array(ii))) != 0:
                combination.append(list(ii))
            if len(np.nonzero(ii)[0]) == 3:
                for i1 in range(3):
                    new_i = list(ii).copy()
                    new_i[i1] = -1 * new_i[i1]
                    combination.append(new_i)
            elif len(np.nonzero(ii)[0]) == 2:
                new_i = list(ii).copy()
                new_i[np.nonzero(ii)[0][0]] = -1 * new_i[np.nonzero(ii)[0][0]]
                combination.append(new_i)

        for ii in combination:  # type: ignore[assignment]
            if reduce(math.gcd, ii) == 1:
                temp = np.dot(np.array(ii), csl)
                if abs(np.dot(temp, surface) - 0) < 1.0e-8:
                    ab_vector.append(temp)
                else:
                    # c vector length along the direction perpendicular to surface
                    c_len_temp = np.abs(np.dot(temp, surface))
                    # c vector length itself
                    c_norm_temp = np.linalg.norm(np.matmul(temp, trans))
                    if normal:
                        c_cross = np.cross(np.matmul(temp, trans), np.matmul(surface, ctrans))
                        if np.linalg.norm(c_cross) < 1.0e-8:
                            if normal_init:
                                if c_norm_temp < c_norm:
                                    t_matrix[2] = temp
                                    c_norm = c_norm_temp
                            else:
                                c_norm = c_norm_temp
                                normal_init = True
                                t_matrix[2] = temp
                    elif c_len_temp < c_length or (abs(c_len_temp - c_length) < 1.0e-8 and c_norm_temp < c_norm):
                        t_matrix[2] = temp
                        c_norm = c_norm_temp
                        c_length = c_len_temp

        if normal and (not normal_init):
            logger.info("Did not find the perpendicular c vector, increase max_j")
            while not normal_init:
                if max_j == max_search:
                    warnings.warn("Cannot find the perpendicular c vector, please increase max_search")
                    break
                max_j *= 3
                max_j = min(max_j, max_search)
                jj = np.arange(0, max_j + 1)
                combination = []
                for ii in product(jj, repeat=3):
                    if sum(abs(np.array(ii))) != 0:
                        combination.append(list(ii))
                    if len(np.nonzero(ii)[0]) == 3:
                        for i1 in range(3):
                            new_i = list(ii).copy()
                            new_i[i1] = -1 * new_i[i1]
                            combination.append(new_i)
                    elif len(np.nonzero(ii)[0]) == 2:
                        new_i = list(ii).copy()
                        new_i[np.nonzero(ii)[0][0]] = -1 * new_i[np.nonzero(ii)[0][0]]
                        combination.append(new_i)

                for ii in combination:  # type: ignore[assignment]
                    if reduce(math.gcd, ii) == 1:
                        temp = np.dot(np.array(ii), csl)
                        if abs(np.dot(temp, surface) - 0) > 1.0e-8:
                            c_cross = np.cross(np.matmul(temp, trans), np.matmul(surface, ctrans))
                            if np.linalg.norm(c_cross) < 1.0e-8:
                                # c vector length itself
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
                    logger.info("Found perpendicular c vector")

        # Find the best a, b vectors with their formed area smallest and average norm of a,b smallest
        ab_norm = None
        for ii in combinations(ab_vector, 2):
            area_temp = np.linalg.norm(np.cross(np.matmul(ii[0], trans), np.matmul(ii[1], trans)))
            if abs(area_temp - 0) > 1.0e-8:
                ab_norm_temp = np.linalg.norm(np.matmul(ii[0], trans)) + np.linalg.norm(np.matmul(ii[1], trans))
                if area is None:
                    area = area_temp
                    ab_norm = ab_norm_temp
                    t_matrix[0] = ii[0]
                    t_matrix[1] = ii[1]

                elif area_temp < area or (abs(area - area_temp) < 1.0e-8 and ab_norm_temp < ab_norm):
                    t_matrix[0] = ii[0]
                    t_matrix[1] = ii[1]
                    area = area_temp
                    ab_norm = ab_norm_temp

        # Make sure we have a left-handed crystallographic system
        if np.linalg.det(np.matmul(t_matrix, trans)) < 0:
            t_matrix *= -1

        if normal and abs(np.linalg.det(t_matrix)) > 1000:
            warnings.warn("Too large matrix. Suggest to use Normal=False")
        return t_matrix

    @staticmethod
    def reduce_mat(mat: NDArray, mag: int, r_matrix: NDArray) -> NDArray:
        """Reduce integer array mat's determinant mag times by linear combination
        of its row vectors, so that the new array after rotation (r_matrix) is
        still an integer array.

        Args:
            mat (3 by 3 array): input matrix
            mag (int): reduce times for the determinant
            r_matrix (3 by 3 array): rotation matrix

        Returns:
            the reduced integer array
        """
        max_j = abs(round(np.linalg.det(mat) / mag))
        reduced = False
        for h in range(3):
            kk = h + 1 if h < 2 else abs(2 - h)
            ll = h + 2 if h < 1 else abs(1 - h)
            jj = np.arange(-max_j, max_j + 1)
            for j1, j2 in product(jj, repeat=2):
                temp = mat[h] + j1 * mat[kk] + j2 * mat[ll]
                if all(np.round(x, 5).is_integer() for x in list(temp / mag)):
                    mat_copy = mat.copy()
                    mat_copy[h] = np.array([round(ele / mag) for ele in temp])
                    new_mat = np.dot(mat_copy, np.linalg.inv(r_matrix.T))
                    if all(np.round(x, 5).is_integer() for x in list(np.ravel(new_mat))):
                        reduced = True
                        mat[h] = np.array([round(ele / mag) for ele in temp])
                        break
            if reduced:
                break

        if not reduced:
            warnings.warn("Matrix reduction not performed, may lead to non-primitive GB cell.")
        return mat

    @staticmethod
    def vec_to_surface(vec: Vector3D) -> MillerIndex:
        """Transform a float vector to a surface miller index with integers.

        Args:
            vec (Vector3D): input float vector

        Returns:
            the surface miller index of the input vector.
        """
        miller: list[None | int] = [None] * 3
        index = []
        for idx, value in enumerate(vec):
            if abs(value) < 1.0e-8:
                miller[idx] = 0
            else:
                index.append(idx)

        if len(index) == 1:
            miller[index[0]] = 1
        else:
            min_index = np.argmin([i for i in vec if i != 0])
            true_index = index[min_index]
            index.pop(min_index)
            frac = []
            for value in index:
                frac.append(Fraction(vec[value] / vec[true_index]).limit_denominator(100))
            if len(index) == 1:
                miller[true_index] = frac[0].denominator
                miller[index[0]] = frac[0].numerator
            else:
                com_lcm = lcm(frac[0].denominator, frac[1].denominator)
                miller[true_index] = com_lcm
                miller[index[0]] = frac[0].numerator * round(com_lcm / frac[0].denominator)
                miller[index[1]] = frac[1].numerator * round(com_lcm / frac[1].denominator)
        return cast(Tuple3Ints, miller)


def fix_pbc(structure: Structure, matrix: NDArray = None) -> Structure:
    """Wrap all frac_coords of the input structure within [0, 1].

    Args:
        structure (pymatgen structure object): input structure
        matrix (lattice matrix, 3 by 3 array/matrix): new structure's lattice matrix,
            If None, use input structure's matrix.


    Returns:
        new structure with fixed frac_coords and lattice matrix
    """
    spec = []
    coords = []
    latte = Lattice(structure.lattice.matrix) if matrix is None else Lattice(matrix)

    for site in structure:
        spec.append(site.specie)
        coord = np.array(site.frac_coords)
        for i in range(3):
            coord[i] -= math.floor(coord[i])
            if np.allclose(coord[i], 1) or np.allclose(coord[i], 0):
                coord[i] = 0
            else:
                coord[i] = round(coord[i], 7)
        coords.append(coord)

    return Structure(latte, spec, coords, site_properties=structure.site_properties)


def symm_group_cubic(mat: NDArray) -> list:
    """Obtain cubic symmetric equivalents of the list of vectors.

    Args:
        mat (np.ndarray): n x 3 lattice matrix


    Returns:
        cubic symmetric equivalents of the list of vectors.
    """
    sym_group = np.zeros([24, 3, 3])
    sym_group[0, :] = np.eye(3)
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


class Interface(Structure):
    """Store data for defining an interface between two Structures."""

    def __init__(
        self,
        lattice: Lattice | NDArray,
        species: list[Any],
        coords: NDArray,
        site_properties: dict[str, Any],
        validate_proximity: bool = False,
        to_unit_cell: bool = False,
        coords_are_cartesian: bool = False,
        in_plane_offset: tuple[float, float] = (0, 0),
        gap: float = 0,
        vacuum_over_film: float = 0,
        interface_properties: dict | None = None,
    ) -> None:
        """Make an Interface, a Structure with additional information
        and methods pertaining to interfaces.

        Args:
            lattice (Lattice | np.ndarray): The lattice, either as a pymatgen.core.Lattice
                or a 3x3 array. Each row should correspond to a lattice
                vector. e.g. [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species ([Species]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / species specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g. (3, 56, ...) or actual Element or Species objects.

                ii. List of dict of elements/species and occupancies, e.g.
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            coords (Nx3 array): list of fractional/cartesian coordinates of
                each species.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            to_unit_cell (bool): Whether to translate sites into the unit cell.
                Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in Cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g. {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
            in_plane_offset: fractional shift in plane for the film with respect
                to the substrate
            gap: gap between substrate and film in Angstroms; zero corresponds to
                the original distance between substrate and film sites
            vacuum_over_film: vacuum space above the film in Angstroms. Defaults to 0.
            interface_properties: properties associated with the Interface. Defaults to None.
        """
        if "interface_label" not in site_properties:
            raise RuntimeError("Must provide labeling of substrate and film sites in site properties")

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
        """The shift between the film and substrate in fractional coordinates."""
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
        """The gap in Cartesian units between the film and the substrate."""
        return self._gap

    @gap.setter
    def gap(self, new_gap: float) -> None:
        if new_gap < 0:
            raise ValueError("Can't reduce interface gap below 0")

        delta = new_gap - self.gap
        self._gap = new_gap

        self._update_c(self.lattice.c + delta)
        self.translate_sites(self.film_indices, [0, 0, delta], frac_coords=False, to_unit_cell=True)

    @property
    def vacuum_over_film(self) -> float:
        """The vacuum space over the film in Cartesian units."""
        return self._vacuum_over_film

    @vacuum_over_film.setter
    def vacuum_over_film(self, new_vacuum: float) -> None:
        if new_vacuum < 0:
            raise ValueError("The vacuum over the film can not be less then 0")

        delta = new_vacuum - self.vacuum_over_film
        self._vacuum_over_film = new_vacuum

        self._update_c(self.lattice.c + delta)

    @property
    def substrate_indices(self) -> list[int]:
        """Site indices for the substrate atoms."""
        return [i for i, tag in enumerate(self.site_properties["interface_label"]) if "substrate" in tag]

    @property
    def substrate_sites(self) -> list[Site]:
        """The site objects in the substrate."""
        return [
            site for site, tag in zip(self, self.site_properties["interface_label"], strict=True) if "substrate" in tag
        ]

    @property
    def substrate(self) -> Structure:
        """A Structure for just the substrate."""
        return Structure.from_sites(self.substrate_sites)

    @property
    def film_indices(self) -> list[int]:
        """Site indices of the film sites."""
        return [i for i, tag in enumerate(self.site_properties["interface_label"]) if "film" in tag]

    @property
    def film_sites(self) -> list[Site]:
        """The film sites of the interface."""
        return [site for site, tag in zip(self, self.site_properties["interface_label"], strict=True) if "film" in tag]

    @property
    def film(self) -> Structure:
        """A Structure for just the film."""
        return Structure.from_sites(self.film_sites)

    def copy(self) -> Self:
        """Make a copy of the Interface."""
        return type(self).from_dict(self.as_dict())

    def get_sorted_structure(
        self,
        key: Callable | None = None,
        reverse: bool = False,
    ) -> Structure:
        """Get a sorted structure for the Interface. The parameters have the same
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

    def get_shifts_based_on_adsorbate_sites(
        self,
        tolerance: float = 0.1,
    ) -> list[tuple[float, float]]:
        """Compute possible in-plane shifts based on an adsorbate site algorithm.

        Args:
            tolerance: tolerance for "uniqueness" for shifts in Cartesian unit
                This is usually Angstroms.
        """
        substrate, film = self.substrate, self.film

        substrate_surface_sites = np.dot(
            list(chain.from_iterable(AdsorbateSiteFinder(substrate).find_adsorption_sites().values())),
            substrate.lattice.inv_matrix,
        )

        # Film gets forced into substrate lattice anyways, so shifts can be
        # computed in fractional coords
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

        def _base_round(x, base: float = 0.05):
            return base * (np.array(x) / base).round()

        # Round shifts to tolerance
        pos_shift[:, 0] = _base_round(pos_shift[:, 0], base=tolerance / substrate.lattice.a)
        pos_shift[:, 1] = _base_round(pos_shift[:, 1], base=tolerance / substrate.lattice.b)
        # C-axis is not useful
        pos_shift = pos_shift[:, 0:2]

        return list(np.unique(pos_shift, axis=0))

    @property
    def film_termination(self) -> str:
        """Label for the film termination chemistry."""
        return label_termination(self.film)

    @property
    def substrate_termination(self) -> str:
        """Label for the substrate termination chemistry."""
        return label_termination(self.substrate)

    @property
    def film_layers(self) -> int:
        """Number of layers of the minimum element in the film composition."""
        sorted_element_list = sorted(
            self.film.composition.element_composition.items(),
            key=lambda x: x[1],
            reverse=True,
        )
        return count_layers(self.film, sorted_element_list[0][0])

    @property
    def substrate_layers(self) -> int:
        """Number of layers of the minimum element in the substrate composition."""
        sorted_element_list = sorted(
            self.substrate.composition.element_composition.items(),
            key=lambda x: x[1],
            reverse=True,
        )
        return count_layers(self.substrate, sorted_element_list[0][0])

    def _update_c(self, new_c: float) -> None:
        """Modify the c-direction of the lattice without changing the site
        Cartesian coordinates.

        WARNING: you can mess up the Interface by setting a c-length that
        can't accommodate all the sites.
        """
        if new_c <= 0:
            raise ValueError("New c-length must be greater than 0")

        new_latt_matrix = [*self.lattice.matrix[:2].tolist(), [0, 0, new_c]]
        new_lattice = Lattice(new_latt_matrix)
        self._lattice = new_lattice

        for site, c_coords in zip(self, self.cart_coords, strict=True):
            site._lattice = new_lattice  # Update the lattice
            site.coords = c_coords  # Put back into original Cartesian space

    def as_dict(self) -> dict:
        """MSONable dict."""
        return {
            **super().as_dict(),
            "in_plane_offset": self.in_plane_offset.tolist(),
            "gap": self.gap,
            "vacuum_over_film": self.vacuum_over_film,
            "interface_properties": self.interface_properties,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct: dict.

        Returns:
            Creates slab from dict.
        """
        lattice = Lattice.from_dict(dct["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in dct["sites"]]
        struct = Structure.from_sites(sites)

        optional = {
            "in_plane_offset": dct.get("in_plane_offset"),
            "gap": dct.get("gap"),
            "vacuum_over_film": dct.get("vacuum_over_film"),
            "interface_properties": dct.get("interface_properties"),
        }
        return cls(
            lattice=lattice,
            species=struct.species_and_occu,
            coords=struct.frac_coords,
            site_properties=struct.site_properties,
            **{k: v for k, v in optional.items() if v is not None},
        )

    @classmethod
    def from_slabs(
        cls,
        substrate_slab: Slab,
        film_slab: Slab,
        in_plane_offset: tuple[float, float] = (0, 0),
        gap: float = 1.6,
        vacuum_over_film: float = 0,
        interface_properties: dict | None = None,
        center_slab: bool = True,
    ) -> Self:
        """Make an Interface by merging a substrate and film slabs
        The film a- and b-vectors will be forced to be the substrate slab's
        a- and b-vectors.

        For now, it's suggested to use a factory method that will ensure the
        appropriate Interface is already met.

        Args:
            substrate_slab (Slab): slab for the substrate.
            film_slab (Slab): slab for the film.
            in_plane_offset (tuple): fractional shift in plane for the film with
                respect to the substrate. For example, (0.5, 0.5) will shift the
                film by half the substrate's a- and b-vectors. Defaults to (0, 0).
            gap (float): gap between substrate and film in Angstroms. Defaults to 1.6.
            vacuum_over_film (float): vacuum space above the film in Angstroms.
                Defaults to 0.
            interface_properties (dict): misc properties to assign to the Interface.
                Defaults to None.
            center_slab (bool): center the slab. Defaults to True.
        """
        # Ensure c-axis is orthogonal to a/b plane
        if isinstance(substrate_slab, Slab):
            substrate_slab = substrate_slab.get_orthogonal_c_slab()
        if isinstance(film_slab, Slab):
            film_slab = film_slab.get_orthogonal_c_slab()
        assert_allclose(film_slab.lattice.alpha, 90, 0.1)
        assert_allclose(film_slab.lattice.beta, 90, 0.1)
        assert_allclose(substrate_slab.lattice.alpha, 90, 0.1)
        assert_allclose(substrate_slab.lattice.beta, 90, 0.1)

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

        # Construct new lattice
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
            coords = np.add(coords, [0, 0, 0.5 - np.mean(coords[:, 2])])

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
            interface_properties=interface_properties or {},
        )

        iface.sort()
        return iface


def label_termination(slab: Structure, ftol: float = 0.25, t_idx: int | None = None) -> str:
    """Label the slab surface termination.

    Args:
        slab (Slab): film or substrate slab to label termination for
        ftol (float): tolerance for terminating position hierarchical clustering
        t_idx (None | int): if not None, adding an extra index to the termination label output
    """
    frac_coords = slab.frac_coords
    n = len(frac_coords)

    if n == 1:
        # Clustering does not work when there is only one data point.
        form = slab.reduced_formula
        sp_symbol = SpacegroupAnalyzer(slab, symprec=0.1).get_space_group_symbol()
        return f"{form}_{sp_symbol}_{len(slab)}"

    dist_matrix = np.zeros((n, n))
    h = slab.lattice.c
    # Projection of c lattice vector in
    # direction of surface normal.
    for ii, jj in combinations(list(range(n)), 2):
        if ii != jj:
            cdist = frac_coords[ii][2] - frac_coords[jj][2]
            cdist = abs(cdist - round(cdist)) * h
            dist_matrix[ii, jj] = cdist
            dist_matrix[jj, ii] = cdist

    condensed_m = squareform(dist_matrix)
    z = linkage(condensed_m)
    clusters = fcluster(z, ftol, criterion="distance")

    clustered_sites: dict[int, list[Site]] = {c: [] for c in clusters}
    for idx, cluster in enumerate(clusters):
        clustered_sites[cluster].append(slab[idx])

    plane_heights = {np.mean(np.mod([s.frac_coords[2] for s in sites], 1)): c for c, sites in clustered_sites.items()}
    top_plane_cluster = max(plane_heights.items(), key=lambda x: x[0])[1]
    top_plane_sites = clustered_sites[top_plane_cluster]
    top_plane = Structure.from_sites(top_plane_sites)

    sp_symbol = SpacegroupAnalyzer(top_plane, symprec=0.1).get_space_group_symbol()
    form = top_plane.reduced_formula

    if t_idx is None:
        return f"{form}_{sp_symbol}_{len(top_plane)}"

    return f"{t_idx}_{form}_{sp_symbol}_{len(top_plane)}"


def count_layers(struct: Structure, el: Element | None = None) -> int:
    """Count the number of layers along the c-axis."""
    el = el or struct.elements[0]
    frac_coords = [site.frac_coords for site in struct if site.species_string == str(el)]
    n_el_sites = len(frac_coords)

    if n_el_sites == 1:
        return 1

    dist_matrix = np.zeros((n_el_sites, n_el_sites))
    h = struct.lattice.c
    # Projection of c lattice vector in
    # direction of surface normal.
    for ii, jj in combinations(list(range(n_el_sites)), 2):
        if ii != jj:
            cdist = frac_coords[ii][2] - frac_coords[jj][2]
            cdist = abs(cdist - round(cdist)) * h
            dist_matrix[ii, jj] = cdist
            dist_matrix[jj, ii] = cdist

    condensed_m = squareform(dist_matrix)
    z = linkage(condensed_m)
    clusters = fcluster(z, 0.25, criterion="distance")

    clustered_sites: dict[int, list[Site]] = {c: [] for c in clusters}
    for idx, cluster in enumerate(clusters):
        clustered_sites[cluster].append(struct[idx])

    plane_heights = {np.mean(np.mod([s.frac_coords[2] for s in sites], 1)): c for c, sites in clustered_sites.items()}

    return len(plane_heights)
