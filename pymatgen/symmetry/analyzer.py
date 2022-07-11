# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
An interface to the excellent spglib library by Atsushi Togo
(http://spglib.sourceforge.net/) for pymatgen.

v1.0 - Now works with both ordered and disordered structure.
v2.0 - Updated for spglib 1.6.
v3.0 - pymatgen no longer ships with spglib. Instead, spglib (the python
       version) is now a dependency and the SpacegroupAnalyzer merely serves
       as an interface to spglib for pymatgen Structures.
"""

import copy
import itertools
import logging
import math
import warnings
from collections import defaultdict
from fractions import Fraction
from math import cos, sin
from typing import Literal

import os

import numpy as np
import spglib

from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Molecule, PeriodicSite, Structure
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.coord import find_in_coord_list, pbc_diff
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.core.operations import MagSymmOp

logger = logging.getLogger(__name__)


class SpacegroupAnalyzer:
    """
    Takes a pymatgen.core.structure.Structure object and a symprec.
    Uses spglib to perform various symmetry finding operations.
    """

    def __init__(self, structure, symprec=0.01, angle_tolerance=5.0):
        """
        Args:
            structure (Structure/IStructure): Structure to find symmetry
            symprec (float): Tolerance for symmetry finding. Defaults to 0.01,
                which is fairly strict and works well for properly refined
                structures with atoms in the proper symmetry coordinates. For
                structures with slight deviations from their proper atomic
                positions (e.g., structures relaxed with electronic structure
                codes), a looser tolerance of 0.1 (the value used in Materials
                Project) is often needed.
            angle_tolerance (float): Angle tolerance for symmetry finding.
        """
        self._symprec = symprec
        self._angle_tol = angle_tolerance
        self._structure = structure
        self._siteprops = structure.site_properties
        latt = structure.lattice.matrix
        positions = structure.frac_coords
        unique_species = []
        zs = []
        magmoms = []

        for species, g in itertools.groupby(structure, key=lambda s: s.species):
            if species in unique_species:
                ind = unique_species.index(species)
                zs.extend([ind + 1] * len(tuple(g)))
            else:
                unique_species.append(species)
                zs.extend([len(unique_species)] * len(tuple(g)))

        for site in structure:
            if hasattr(site, "magmom"):
                magmoms.append(site.magmom)
            elif site.is_ordered and hasattr(site.specie, "spin"):
                magmoms.append(site.specie.spin)
            else:
                magmoms.append(0)  # needed for spglib

        self._unique_species = unique_species
        self._numbers = zs
        self._cell = latt, positions, zs, magmoms

        self._space_group_data = spglib.get_symmetry_dataset(
            self._cell, symprec=self._symprec, angle_tolerance=angle_tolerance
        )

    def get_space_group_symbol(self):
        """
        Get the spacegroup symbol (e.g., Pnma) for structure.

        Returns:
            (str): Spacegroup symbol for structure.
        """
        return self._space_group_data["international"]

    def get_space_group_number(self):
        """
        Get the international spacegroup number (e.g., 62) for structure.

        Returns:
            (int): International spacegroup number for structure.
        """
        return int(self._space_group_data["number"])

    def get_space_group_operations(self):
        """
        Get the SpacegroupOperations for the Structure.

        Returns:
            SpacegroupOperations object.
        """
        return SpacegroupOperations(
            self.get_space_group_symbol(),
            self.get_space_group_number(),
            self.get_symmetry_operations(),
        )

    def get_hall(self):
        """
        Returns Hall symbol for structure.

        Returns:
            (str): Hall symbol
        """
        return self._space_group_data["hall"]

    def get_point_group_symbol(self):
        """
        Get the point group associated with the structure.

        Returns:
            (Pointgroup): Point group for structure.
        """
        rotations = self._space_group_data["rotations"]
        # passing a 0-length rotations list to spglib can segfault
        if len(rotations) == 0:
            return "1"
        return spglib.get_pointgroup(rotations)[0].strip()

    def get_crystal_system(
        self,
    ) -> Literal["triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"]:
        """
        Get the crystal system for the structure, e.g., (triclinic,
        orthorhombic, cubic, etc.).

        Raises:
            ValueError: on invalid space group numbers < 1 or > 230.

        Returns:
            (str): Crystal system for structure
        """
        n = self._space_group_data["number"]

        # not using isinstance(n, int) to allow 0-decimal floats
        if not (n == int(n) and 0 < n < 231):
            raise ValueError(f"Received invalid space group {n}")

        if 0 < n < 3:
            return "triclinic"
        if n < 16:
            return "monoclinic"
        if n < 75:
            return "orthorhombic"
        if n < 143:
            return "tetragonal"
        if n < 168:
            return "trigonal"
        if n < 195:
            return "hexagonal"
        return "cubic"

    def get_lattice_type(
        self,
    ) -> Literal["triclinic", "monoclinic", "orthorhombic", "tetragonal", "rhombohedral", "hexagonal", "cubic"]:
        """
        Get the lattice for the structure, e.g., (triclinic, orthorhombic, cubic, etc.).This is
        the same as the crystal system with the exception of the hexagonal/rhombohedral lattice.

        Raises:
            ValueError: on invalid space group numbers < 1 or > 230.

        Returns:
            (str): Lattice type for structure
        """
        n = self._space_group_data["number"]
        system = self.get_crystal_system()
        if n in [146, 148, 155, 160, 161, 166, 167]:
            return "rhombohedral"
        if system == "trigonal":
            return "hexagonal"
        return system

    def get_symmetry_dataset(self):
        """
        Returns the symmetry dataset as a dict.

        Returns:
            (dict): With the following properties:
                number: International space group number
                international: International symbol
                hall: Hall symbol
                transformation_matrix: Transformation matrix from lattice of
                input cell to Bravais lattice L^bravais = L^original * Tmat
                origin shift: Origin shift in the setting of "Bravais lattice"
                rotations, translations: Rotation matrices and translation
                vectors. Space group operations are obtained by
                [(r,t) for r, t in zip(rotations, translations)]
                wyckoffs: Wyckoff letters
        """
        return self._space_group_data

    def _get_symmetry(self):
        """
        Get the symmetry operations associated with the structure.

        Returns:
            Symmetry operations as a tuple of two equal length sequences.
            (rotations, translations). "rotations" is the numpy integer array
            of the rotation matrices for scaled positions
            "translations" gives the numpy float64 array of the translation
            vectors in scaled positions.
        """
        d = spglib.get_symmetry(self._cell, symprec=self._symprec, angle_tolerance=self._angle_tol)
        # Sometimes spglib returns small translation vectors, e.g.
        # [1e-4, 2e-4, 1e-4]
        # (these are in fractional coordinates, so should be small denominator
        # fractions)
        trans = []
        for t in d["translations"]:
            trans.append([float(Fraction.from_float(c).limit_denominator(1000)) for c in t])
        trans = np.array(trans)

        # fractional translations of 1 are more simply 0
        trans[np.abs(trans) == 1] = 0
        return d["rotations"], trans

    def get_symmetry_operations(self, cartesian=False):
        """
        Return symmetry operations as a list of SymmOp objects.
        By default returns fractional coord symmops.
        But Cartesian can be returned too.

        Returns:
            ([SymmOp]): List of symmetry operations.
        """
        rotation, translation = self._get_symmetry()
        symmops = []
        mat = self._structure.lattice.matrix.T
        invmat = np.linalg.inv(mat)
        for rot, trans in zip(rotation, translation):
            if cartesian:
                rot = np.dot(mat, np.dot(rot, invmat))
                trans = np.dot(trans, self._structure.lattice.matrix)
            op = SymmOp.from_rotation_and_translation(rot, trans)
            symmops.append(op)
        return symmops

    def get_point_group_operations(self, cartesian=False):
        """
        Return symmetry operations as a list of SymmOp objects.
        By default returns fractional coord symmops.
        But Cartesian can be returned too.

        Args:
            cartesian (bool): Whether to return SymmOps as Cartesian or
                direct coordinate operations.

        Returns:
            ([SymmOp]): List of point group symmetry operations.
        """
        rotation, translation = self._get_symmetry()
        symmops = []
        mat = self._structure.lattice.matrix.T
        invmat = np.linalg.inv(mat)
        for rot in rotation:
            if cartesian:
                rot = np.dot(mat, np.dot(rot, invmat))
            op = SymmOp.from_rotation_and_translation(rot, np.array([0, 0, 0]))
            symmops.append(op)
        return symmops

    def get_symmetrized_structure(self):
        """
        Get a symmetrized structure. A symmetrized structure is one where the
        sites have been grouped into symmetrically equivalent groups.

        Returns:
            :class:`pymatgen.symmetry.structure.SymmetrizedStructure` object.
        """
        ds = self.get_symmetry_dataset()
        sg = SpacegroupOperations(
            self.get_space_group_symbol(),
            self.get_space_group_number(),
            self.get_symmetry_operations(),
        )
        return SymmetrizedStructure(self._structure, sg, ds["equivalent_atoms"], ds["wyckoffs"])

    def get_refined_structure(self, keep_site_properties=False):
        """
        Get the refined structure based on detected symmetry. The refined
        structure is a *conventional* cell setting with atoms moved to the
        expected symmetry positions.

        Args:
            keep_site_properties (bool): Whether to keep the input site properties (including
                magnetic moments) on the sites that are still present after the refinement. Note:
                This is disabled by default because the magnetic moments are not always directly
                transferable between unit cell definitions. For instance, long-range magnetic
                ordering or antiferromagnetic character may no longer be present (or exist in
                the same way) in the returned structure. If keep_site_properties is True,
                each site retains the same site property as in the original structure without
                further adjustment.

        Returns:
            Refined structure.
        """
        # Atomic positions have to be specified by scaled positions for spglib.
        lattice, scaled_positions, numbers = spglib.refine_cell(self._cell, self._symprec, self._angle_tol)
        species = [self._unique_species[i - 1] for i in numbers]
        if keep_site_properties:
            site_properties = {}
            for k, v in self._siteprops.items():
                site_properties[k] = [v[i - 1] for i in numbers]
        else:
            site_properties = None
        s = Structure(lattice, species, scaled_positions, site_properties=site_properties)
        return s.get_sorted_structure()

    def find_primitive(self, keep_site_properties=False):
        """
        Find a primitive version of the unit cell.

        Args:
            keep_site_properties (bool): Whether to keep the input site properties (including
                magnetic moments) on the sites that are still present after the refinement. Note:
                This is disabled by default because the magnetic moments are not always directly
                transferable between unit cell definitions. For instance, long-range magnetic
                ordering or antiferromagnetic character may no longer be present (or exist in
                the same way) in the returned structure. If keep_site_properties is True,
                each site retains the same site property as in the original structure without
                further adjustment.

        Returns:
            A primitive cell in the input cell is searched and returned
            as an Structure object. If no primitive cell is found, None is
            returned.
        """
        lattice, scaled_positions, numbers = spglib.find_primitive(self._cell, symprec=self._symprec)
        species = [self._unique_species[i - 1] for i in numbers]
        if keep_site_properties:
            site_properties = {}
            for k, v in self._siteprops.items():
                site_properties[k] = [v[i - 1] for i in numbers]
        else:
            site_properties = None

        return Structure(
            lattice, species, scaled_positions, to_unit_cell=True, site_properties=site_properties
        ).get_reduced_structure()

    def get_ir_reciprocal_mesh(self, mesh=(10, 10, 10), is_shift=(0, 0, 0)):
        """
        k-point mesh of the Brillouin zone generated taken into account
        symmetry.The method returns the irreducible kpoints of the mesh
        and their weights

        Args:
            mesh (3x1 array): The number of kpoint for the mesh needed in
                each direction
            is_shift (3x1 array): Whether to shift the kpoint grid. (1, 1,
            1) means all points are shifted by 0.5, 0.5, 0.5.

        Returns:
            A list of irreducible kpoints and their weights as a list of
            tuples [(ir_kpoint, weight)], with ir_kpoint given
            in fractional coordinates
        """
        shift = np.array([1 if i else 0 for i in is_shift])
        mapping, grid = spglib.get_ir_reciprocal_mesh(np.array(mesh), self._cell, is_shift=shift, symprec=self._symprec)

        results = []
        for i, count in zip(*np.unique(mapping, return_counts=True)):
            results.append(((grid[i] + shift * (0.5, 0.5, 0.5)) / mesh, count))
        return results

    def get_conventional_to_primitive_transformation_matrix(self, international_monoclinic=True):
        """
        Gives the transformation matrix to transform a conventional
        unit cell to a primitive cell according to certain standards
        the standards are defined in Setyawan, W., & Curtarolo, S. (2010).
        High-throughput electronic band structure calculations:
        Challenges and tools. Computational Materials Science,
        49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010

        Args:
            international_monoclinic (bool): Whether to convert to proper international convention
                such that beta is the non-right angle.

        Returns:
            Transformation matrix to go from conventional to primitive cell
        """
        conv = self.get_conventional_standard_structure(international_monoclinic=international_monoclinic)
        lattice = self.get_lattice_type()

        if "P" in self.get_space_group_symbol() or lattice == "hexagonal":
            return np.eye(3)

        if lattice == "rhombohedral":
            # check if the conventional representation is hexagonal or
            # rhombohedral
            lengths = conv.lattice.lengths
            if abs(lengths[0] - lengths[2]) < 0.0001:
                transf = np.eye
            else:
                transf = np.array([[-1, 1, 1], [2, 1, 1], [-1, -2, 1]], dtype=np.float_) / 3

        elif "I" in self.get_space_group_symbol():
            transf = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype=np.float_) / 2
        elif "F" in self.get_space_group_symbol():
            transf = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=np.float_) / 2
        elif "C" in self.get_space_group_symbol() or "A" in self.get_space_group_symbol():
            if self.get_crystal_system() == "monoclinic":
                transf = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 2]], dtype=np.float_) / 2
            else:
                transf = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 2]], dtype=np.float_) / 2
        else:
            transf = np.eye(3)

        return transf

    def get_primitive_standard_structure(self, international_monoclinic=True, keep_site_properties=False):
        """
        Gives a structure with a primitive cell according to certain standards
        the standards are defined in Setyawan, W., & Curtarolo, S. (2010).
        High-throughput electronic band structure calculations:
        Challenges and tools. Computational Materials Science,
        49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010

        Args:
            international_monoclinic (bool): Whether to convert to proper international convention
                such that beta is the non-right angle.
            keep_site_properties (bool): Whether to keep the input site properties (including
                magnetic moments) on the sites that are still present after the refinement. Note:
                This is disabled by default because the magnetic moments are not always directly
                transferable between unit cell definitions. For instance, long-range magnetic
                ordering or antiferromagnetic character may no longer be present (or exist in
                the same way) in the returned structure. If keep_site_properties is True,
                each site retains the same site property as in the original structure without
                further adjustment.

        Returns:
            The structure in a primitive standardized cell
        """
        conv = self.get_conventional_standard_structure(
            international_monoclinic=international_monoclinic, keep_site_properties=keep_site_properties
        )
        lattice = self.get_lattice_type()

        if "P" in self.get_space_group_symbol() or lattice == "hexagonal":
            return conv

        transf = self.get_conventional_to_primitive_transformation_matrix(
            international_monoclinic=international_monoclinic
        )

        new_sites = []
        latt = Lattice(np.dot(transf, conv.lattice.matrix))
        for s in conv:
            new_s = PeriodicSite(
                s.specie,
                s.coords,
                latt,
                to_unit_cell=True,
                coords_are_cartesian=True,
                properties=s.properties,
            )
            if not any(map(new_s.is_periodic_image, new_sites)):
                new_sites.append(new_s)

        if lattice == "rhombohedral":
            prim = Structure.from_sites(new_sites)
            lengths = prim.lattice.lengths
            angles = prim.lattice.angles
            a = lengths[0]
            alpha = math.pi * angles[0] / 180
            new_matrix = [
                [a * cos(alpha / 2), -a * sin(alpha / 2), 0],
                [a * cos(alpha / 2), a * sin(alpha / 2), 0],
                [
                    a * cos(alpha) / cos(alpha / 2),
                    0,
                    a * math.sqrt(1 - (cos(alpha) ** 2 / (cos(alpha / 2) ** 2))),
                ],
            ]
            new_sites = []
            latt = Lattice(new_matrix)
            for s in prim:
                new_s = PeriodicSite(
                    s.specie,
                    s.frac_coords,
                    latt,
                    to_unit_cell=True,
                    properties=s.properties,
                )
                if not any(map(new_s.is_periodic_image, new_sites)):
                    new_sites.append(new_s)
            return Structure.from_sites(new_sites)

        return Structure.from_sites(new_sites)

    def get_conventional_standard_structure(self, international_monoclinic=True, keep_site_properties=False):
        """
        Gives a structure with a conventional cell according to certain
        standards. The standards are defined in Setyawan, W., & Curtarolo,
        S. (2010). High-throughput electronic band structure calculations:
        Challenges and tools. Computational Materials Science,
        49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010
        They basically enforce as much as possible
        norm(a1)<norm(a2)<norm(a3). NB This is not necessarily the same as the
        standard settings within the International Tables of Crystallography,
        for which get_refined_structure should be used instead.

        Args:
            international_monoclinic (bool): Whether to convert to proper international convention
                such that beta is the non-right angle.
            keep_site_properties (bool): Whether to keep the input site properties (including
                magnetic moments) on the sites that are still present after the refinement. Note:
                This is disabled by default because the magnetic moments are not always directly
                transferable between unit cell definitions. For instance, long-range magnetic
                ordering or antiferromagnetic character may no longer be present (or exist in
                the same way) in the returned structure. If keep_site_properties is True,
                each site retains the same site property as in the original structure without
                further adjustment.

        Returns:
            The structure in a conventional standardized cell
        """
        tol = 1e-5
        struct = self.get_refined_structure(keep_site_properties=keep_site_properties)
        latt = struct.lattice
        latt_type = self.get_lattice_type()
        sorted_lengths = sorted(latt.abc)
        sorted_dic = sorted(
            ({"vec": latt.matrix[i], "length": latt.abc[i], "orig_index": i} for i in [0, 1, 2]),
            key=lambda k: k["length"],
        )

        if latt_type in ("orthorhombic", "cubic"):
            # you want to keep the c axis where it is
            # to keep the C- settings
            transf = np.zeros(shape=(3, 3))
            if self.get_space_group_symbol().startswith("C"):
                transf[2] = [0, 0, 1]
                a, b = sorted(latt.abc[:2])
                sorted_dic = sorted(
                    ({"vec": latt.matrix[i], "length": latt.abc[i], "orig_index": i} for i in [0, 1]),
                    key=lambda k: k["length"],
                )
                for i in range(2):
                    transf[i][sorted_dic[i]["orig_index"]] = 1
                c = latt.abc[2]
            elif self.get_space_group_symbol().startswith(
                "A"
            ):  # change to C-centering to match Setyawan/Curtarolo convention
                transf[2] = [1, 0, 0]
                a, b = sorted(latt.abc[1:])
                sorted_dic = sorted(
                    ({"vec": latt.matrix[i], "length": latt.abc[i], "orig_index": i} for i in [1, 2]),
                    key=lambda k: k["length"],
                )
                for i in range(2):
                    transf[i][sorted_dic[i]["orig_index"]] = 1
                c = latt.abc[0]
            else:
                for i, d in enumerate(sorted_dic):
                    transf[i][d["orig_index"]] = 1
                a, b, c = sorted_lengths
            latt = Lattice.orthorhombic(a, b, c)

        elif latt_type == "tetragonal":
            # find the "a" vectors
            # it is basically the vector repeated two times
            transf = np.zeros(shape=(3, 3))
            a, b, c = sorted_lengths
            for i, d in enumerate(sorted_dic):
                transf[i][d["orig_index"]] = 1

            if abs(b - c) < tol < abs(a - c):
                a, c = c, a
                transf = np.dot([[0, 0, 1], [0, 1, 0], [1, 0, 0]], transf)
            latt = Lattice.tetragonal(a, c)
        elif latt_type in ("hexagonal", "rhombohedral"):
            # for the conventional cell representation,
            # we always show the rhombohedral lattices as hexagonal

            # check first if we have the refined structure shows a rhombohedral
            # cell
            # if so, make a supercell
            a, b, c = latt.abc
            if np.all(np.abs([a - b, c - b, a - c]) < 0.001):
                struct.make_supercell(((1, -1, 0), (0, 1, -1), (1, 1, 1)))
                a, b, c = sorted(struct.lattice.abc)

            if abs(b - c) < 0.001:
                a, c = c, a
            new_matrix = [
                [a / 2, -a * math.sqrt(3) / 2, 0],
                [a / 2, a * math.sqrt(3) / 2, 0],
                [0, 0, c],
            ]
            latt = Lattice(new_matrix)
            transf = np.eye(3, 3)

        elif latt_type == "monoclinic":
            # You want to keep the c axis where it is to keep the C- settings

            if self.get_space_group_operations().int_symbol.startswith("C"):
                transf = np.zeros(shape=(3, 3))
                transf[2] = [0, 0, 1]
                sorted_dic = sorted(
                    ({"vec": latt.matrix[i], "length": latt.abc[i], "orig_index": i} for i in [0, 1]),
                    key=lambda k: k["length"],
                )
                a = sorted_dic[0]["length"]
                b = sorted_dic[1]["length"]
                c = latt.abc[2]
                new_matrix = None
                for t in itertools.permutations(list(range(2)), 2):
                    m = latt.matrix
                    latt2 = Lattice([m[t[0]], m[t[1]], m[2]])
                    lengths = latt2.lengths
                    angles = latt2.angles
                    if angles[0] > 90:
                        # if the angle is > 90 we invert a and b to get
                        # an angle < 90
                        a, b, c, alpha, beta, gamma = Lattice([-m[t[0]], -m[t[1]], m[2]]).parameters
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = -1
                        transf[1][t[1]] = -1
                        transf[2][2] = 1
                        alpha = math.pi * alpha / 180
                        new_matrix = [
                            [a, 0, 0],
                            [0, b, 0],
                            [0, c * cos(alpha), c * sin(alpha)],
                        ]
                        continue

                    if angles[0] < 90:
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = 1
                        transf[1][t[1]] = 1
                        transf[2][2] = 1
                        a, b, c = lengths
                        alpha = math.pi * angles[0] / 180
                        new_matrix = [
                            [a, 0, 0],
                            [0, b, 0],
                            [0, c * cos(alpha), c * sin(alpha)],
                        ]

                if new_matrix is None:
                    # this if is to treat the case
                    # where alpha==90 (but we still have a monoclinic sg
                    new_matrix = [[a, 0, 0], [0, b, 0], [0, 0, c]]
                    transf = np.zeros(shape=(3, 3))
                    transf[2] = [0, 0, 1]  # see issue #1929
                    for i, d in enumerate(sorted_dic):
                        transf[i][d["orig_index"]] = 1
            # if not C-setting
            else:
                # try all permutations of the axis
                # keep the ones with the non-90 angle=alpha
                # and b<c
                new_matrix = None
                for t in itertools.permutations(list(range(3)), 3):
                    m = latt.matrix
                    a, b, c, alpha, beta, gamma = Lattice([m[t[0]], m[t[1]], m[t[2]]]).parameters
                    if alpha > 90 and b < c:
                        a, b, c, alpha, beta, gamma = Lattice([-m[t[0]], -m[t[1]], m[t[2]]]).parameters
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = -1
                        transf[1][t[1]] = -1
                        transf[2][t[2]] = 1
                        alpha = math.pi * alpha / 180
                        new_matrix = [
                            [a, 0, 0],
                            [0, b, 0],
                            [0, c * cos(alpha), c * sin(alpha)],
                        ]
                        continue
                    if alpha < 90 and b < c:
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = 1
                        transf[1][t[1]] = 1
                        transf[2][t[2]] = 1
                        alpha = math.pi * alpha / 180
                        new_matrix = [
                            [a, 0, 0],
                            [0, b, 0],
                            [0, c * cos(alpha), c * sin(alpha)],
                        ]
                if new_matrix is None:
                    # this if is to treat the case
                    # where alpha==90 (but we still have a monoclinic sg
                    new_matrix = [
                        [sorted_lengths[0], 0, 0],
                        [0, sorted_lengths[1], 0],
                        [0, 0, sorted_lengths[2]],
                    ]
                    transf = np.zeros(shape=(3, 3))
                    for i, d in enumerate(sorted_dic):
                        transf[i][d["orig_index"]] = 1

            if international_monoclinic:
                # The above code makes alpha the non-right angle.
                # The following will convert to proper international convention
                # that beta is the non-right angle.
                op = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
                transf = np.dot(op, transf)
                new_matrix = np.dot(op, new_matrix)
                beta = Lattice(new_matrix).beta
                if beta < 90:
                    op = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
                    transf = np.dot(op, transf)
                    new_matrix = np.dot(op, new_matrix)

            latt = Lattice(new_matrix)

        elif latt_type == "triclinic":
            # we use a LLL Minkowski-like reduction for the triclinic cells
            struct = struct.get_reduced_structure("LLL")
            latt = struct.lattice

            a, b, c = latt.lengths
            alpha, beta, gamma = (math.pi * i / 180 for i in latt.angles)
            new_matrix = None
            test_matrix = [
                [a, 0, 0],
                [b * cos(gamma), b * sin(gamma), 0.0],
                [
                    c * cos(beta),
                    c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma),
                    c
                    * math.sqrt(
                        sin(gamma) ** 2 - cos(alpha) ** 2 - cos(beta) ** 2 + 2 * cos(alpha) * cos(beta) * cos(gamma)
                    )
                    / sin(gamma),
                ],
            ]

            def is_all_acute_or_obtuse(m):
                recp_angles = np.array(Lattice(m).reciprocal_lattice.angles)
                return np.all(recp_angles <= 90) or np.all(recp_angles > 90)

            if is_all_acute_or_obtuse(test_matrix):
                transf = np.eye(3)
                new_matrix = test_matrix

            test_matrix = [
                [-a, 0, 0],
                [b * cos(gamma), b * sin(gamma), 0.0],
                [
                    -c * cos(beta),
                    -c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma),
                    -c
                    * math.sqrt(
                        sin(gamma) ** 2 - cos(alpha) ** 2 - cos(beta) ** 2 + 2 * cos(alpha) * cos(beta) * cos(gamma)
                    )
                    / sin(gamma),
                ],
            ]

            if is_all_acute_or_obtuse(test_matrix):
                transf = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
                new_matrix = test_matrix

            test_matrix = [
                [-a, 0, 0],
                [-b * cos(gamma), -b * sin(gamma), 0.0],
                [
                    c * cos(beta),
                    c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma),
                    c
                    * math.sqrt(
                        sin(gamma) ** 2 - cos(alpha) ** 2 - cos(beta) ** 2 + 2 * cos(alpha) * cos(beta) * cos(gamma)
                    )
                    / sin(gamma),
                ],
            ]

            if is_all_acute_or_obtuse(test_matrix):
                transf = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
                new_matrix = test_matrix

            test_matrix = [
                [a, 0, 0],
                [-b * cos(gamma), -b * sin(gamma), 0.0],
                [
                    -c * cos(beta),
                    -c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma),
                    -c
                    * math.sqrt(
                        sin(gamma) ** 2 - cos(alpha) ** 2 - cos(beta) ** 2 + 2 * cos(alpha) * cos(beta) * cos(gamma)
                    )
                    / sin(gamma),
                ],
            ]
            if is_all_acute_or_obtuse(test_matrix):
                transf = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
                new_matrix = test_matrix

            latt = Lattice(new_matrix)

        new_coords = np.dot(transf, np.transpose(struct.frac_coords)).T
        new_struct = Structure(
            latt,
            struct.species_and_occu,
            new_coords,
            site_properties=struct.site_properties,
            to_unit_cell=True,
        )
        return new_struct.get_sorted_structure()

    def get_kpoint_weights(self, kpoints, atol=1e-5):
        """
        Calculate the weights for a list of kpoints.

        Args:
            kpoints (Sequence): Sequence of kpoints. np.arrays is fine. Note
                that the code does not check that the list of kpoints
                provided does not contain duplicates.
            atol (float): Tolerance for fractional coordinates comparisons.

        Returns:
            List of weights, in the SAME order as kpoints.
        """
        kpts = np.array(kpoints)
        shift = []
        mesh = []
        for i in range(3):
            nonzero = [i for i in kpts[:, i] if abs(i) > 1e-5]
            if len(nonzero) != len(kpts):
                # gamma centered
                if not nonzero:
                    mesh.append(1)
                else:
                    m = np.abs(np.round(1 / np.array(nonzero)))
                    mesh.append(int(max(m)))
                shift.append(0)
            else:
                # Monk
                m = np.abs(np.round(0.5 / np.array(nonzero)))
                mesh.append(int(max(m)))
                shift.append(1)

        mapping, grid = spglib.get_ir_reciprocal_mesh(np.array(mesh), self._cell, is_shift=shift, symprec=self._symprec)
        mapping = list(mapping)
        grid = (np.array(grid) + np.array(shift) * (0.5, 0.5, 0.5)) / mesh
        weights = []
        mapped = defaultdict(int)
        for k in kpoints:
            for i, g in enumerate(grid):
                if np.allclose(pbc_diff(k, g), (0, 0, 0), atol=atol):
                    mapped[tuple(g)] += 1
                    weights.append(mapping.count(mapping[i]))
                    break
        if (len(mapped) != len(set(mapping))) or (not all(v == 1 for v in mapped.values())):
            raise ValueError("Unable to find 1:1 corresponding between input kpoints and irreducible grid!")
        return [w / sum(weights) for w in weights]

    def is_laue(self):
        """
        Check if the point group of the structure
            has Laue symmetry (centrosymmetry)
        """

        laue = [
            "-1",
            "2/m",
            "mmm",
            "4/m",
            "4/mmm",
            "-3",
            "-3m",
            "6/m",
            "6/mmm",
            "m-3",
            "m-3m",
        ]

        return str(self.get_point_group_symbol()) in laue


class PointGroupAnalyzer:
    """
    A class to analyze the point group of a molecule. The general outline of
    the algorithm is as follows:

    1. Center the molecule around its center of mass.
    2. Compute the inertia tensor and the eigenvalues and eigenvectors.
    3. Handle the symmetry detection based on eigenvalues.

        a. Linear molecules have one zero eigenvalue. Possible symmetry
           operations are C*v or D*v
        b. Asymmetric top molecules have all different eigenvalues. The
           maximum rotational symmetry in such molecules is 2
        c. Symmetric top molecules have 1 unique eigenvalue, which gives a
           unique rotation axis. All axial point groups are possible
           except the cubic groups (T & O) and I.
        d. Spherical top molecules have all three eigenvalues equal. They
           have the rare T, O or I point groups.

    .. attribute:: sch_symbol

        Schoenflies symbol of the detected point group.
    """

    inversion_op = SymmOp.inversion()

    def __init__(self, mol, tolerance=0.3, eigen_tolerance=0.01, matrix_tolerance=0.1):
        """
        The default settings are usually sufficient.

        Args:
            mol (Molecule): Molecule to determine point group for.
            tolerance (float): Distance tolerance to consider sites as
                symmetrically equivalent. Defaults to 0.3 Angstrom.
            eigen_tolerance (float): Tolerance to compare eigen values of
                the inertia tensor. Defaults to 0.01.
            matrix_tolerance (float): Tolerance used to generate the full set of
                symmetry operations of the point group.
        """
        self.mol = mol
        self.centered_mol = mol.get_centered_molecule()
        self.tol = tolerance
        self.eig_tol = eigen_tolerance
        self.mat_tol = matrix_tolerance
        self._analyze()
        if self.sch_symbol in ["C1v", "C1h"]:
            self.sch_symbol = "Cs"

    def _analyze(self):
        if len(self.centered_mol) == 1:
            self.sch_symbol = "Kh"
        else:
            inertia_tensor = np.zeros((3, 3))
            total_inertia = 0
            for site in self.centered_mol:
                c = site.coords
                wt = site.species.weight
                for i in range(3):
                    inertia_tensor[i, i] += wt * (c[(i + 1) % 3] ** 2 + c[(i + 2) % 3] ** 2)
                for i, j in [(0, 1), (1, 2), (0, 2)]:
                    inertia_tensor[i, j] += -wt * c[i] * c[j]
                    inertia_tensor[j, i] += -wt * c[j] * c[i]
                total_inertia += wt * np.dot(c, c)

            # Normalize the inertia tensor so that it does not scale with size
            # of the system. This mitigates the problem of choosing a proper
            # comparison tolerance for the eigenvalues.
            inertia_tensor /= total_inertia
            eigvals, eigvecs = np.linalg.eig(inertia_tensor)
            self.principal_axes = eigvecs.T
            self.eigvals = eigvals
            v1, v2, v3 = eigvals
            eig_zero = abs(v1 * v2 * v3) < self.eig_tol
            eig_all_same = abs(v1 - v2) < self.eig_tol and abs(v1 - v3) < self.eig_tol
            eig_all_diff = abs(v1 - v2) > self.eig_tol and abs(v1 - v3) > self.eig_tol and abs(v2 - v3) > self.eig_tol

            self.rot_sym = []
            self.symmops = [SymmOp(np.eye(4))]
            if eig_zero:
                logger.debug("Linear molecule detected")
                self._proc_linear()
            elif eig_all_same:
                logger.debug("Spherical top molecule detected")
                self._proc_sph_top()
            elif eig_all_diff:
                logger.debug("Asymmetric top molecule detected")
                self._proc_asym_top()
            else:
                logger.debug("Symmetric top molecule detected")
                self._proc_sym_top()

    def _proc_linear(self):
        if self.is_valid_op(PointGroupAnalyzer.inversion_op):
            self.sch_symbol = "D*h"
            self.symmops.append(PointGroupAnalyzer.inversion_op)
        else:
            self.sch_symbol = "C*v"

    def _proc_asym_top(self):
        """
        Handles asymmetric top molecules, which cannot contain rotational
        symmetry larger than 2.
        """
        self._check_R2_axes_asym()
        if len(self.rot_sym) == 0:
            logger.debug("No rotation symmetries detected.")
            self._proc_no_rot_sym()
        elif len(self.rot_sym) == 3:
            logger.debug("Dihedral group detected.")
            self._proc_dihedral()
        else:
            logger.debug("Cyclic group detected.")
            self._proc_cyclic()

    def _proc_sym_top(self):
        """
        Handles symmetric top molecules which has one unique eigenvalue whose
        corresponding principal axis is a unique rotational axis. More complex
        handling required to look for R2 axes perpendicular to this unique
        axis.
        """
        if abs(self.eigvals[0] - self.eigvals[1]) < self.eig_tol:
            ind = 2
        elif abs(self.eigvals[1] - self.eigvals[2]) < self.eig_tol:
            ind = 0
        else:
            ind = 1
        logger.debug(f"Eigenvalues = {self.eigvals}.")
        unique_axis = self.principal_axes[ind]
        self._check_rot_sym(unique_axis)
        logger.debug(f"Rotation symmetries = {self.rot_sym}")
        if len(self.rot_sym) > 0:
            self._check_perpendicular_r2_axis(unique_axis)

        if len(self.rot_sym) >= 2:
            self._proc_dihedral()
        elif len(self.rot_sym) == 1:
            self._proc_cyclic()
        else:
            self._proc_no_rot_sym()

    def _proc_no_rot_sym(self):
        """
        Handles molecules with no rotational symmetry. Only possible point
        groups are C1, Cs and Ci.
        """
        self.sch_symbol = "C1"
        if self.is_valid_op(PointGroupAnalyzer.inversion_op):
            self.sch_symbol = "Ci"
            self.symmops.append(PointGroupAnalyzer.inversion_op)
        else:
            for v in self.principal_axes:
                mirror_type = self._find_mirror(v)
                if not mirror_type == "":
                    self.sch_symbol = "Cs"
                    break

    def _proc_cyclic(self):
        """
        Handles cyclic group molecules.
        """
        main_axis, rot = max(self.rot_sym, key=lambda v: v[1])
        self.sch_symbol = f"C{rot}"
        mirror_type = self._find_mirror(main_axis)
        if mirror_type == "h":
            self.sch_symbol += "h"
        elif mirror_type == "v":
            self.sch_symbol += "v"
        elif mirror_type == "":
            if self.is_valid_op(SymmOp.rotoreflection(main_axis, angle=180 / rot)):
                self.sch_symbol = f"S{2 * rot}"

    def _proc_dihedral(self):
        """
        Handles dihedral group molecules, i.e those with intersecting R2 axes
        and a main axis.
        """
        main_axis, rot = max(self.rot_sym, key=lambda v: v[1])
        self.sch_symbol = f"D{rot}"
        mirror_type = self._find_mirror(main_axis)
        if mirror_type == "h":
            self.sch_symbol += "h"
        elif not mirror_type == "":
            self.sch_symbol += "d"

    def _check_R2_axes_asym(self):
        """
        Test for 2-fold rotation along the principal axes. Used to handle
        asymmetric top molecules.
        """
        for v in self.principal_axes:
            op = SymmOp.from_axis_angle_and_translation(v, 180)
            if self.is_valid_op(op):
                self.symmops.append(op)
                self.rot_sym.append((v, 2))

    def _find_mirror(self, axis):
        """
        Looks for mirror symmetry of specified type about axis. Possible
        types are "h" or "vd". Horizontal (h) mirrors are perpendicular to
        the axis while vertical (v) or diagonal (d) mirrors are parallel. v
        mirrors has atoms lying on the mirror plane while d mirrors do
        not.
        """
        mirror_type = ""

        # First test whether the axis itself is the normal to a mirror plane.
        if self.is_valid_op(SymmOp.reflection(axis)):
            self.symmops.append(SymmOp.reflection(axis))
            mirror_type = "h"
        else:
            # Iterate through all pairs of atoms to find mirror
            for s1, s2 in itertools.combinations(self.centered_mol, 2):
                if s1.species == s2.species:
                    normal = s1.coords - s2.coords
                    if np.dot(normal, axis) < self.tol:
                        op = SymmOp.reflection(normal)
                        if self.is_valid_op(op):
                            self.symmops.append(op)
                            if len(self.rot_sym) > 1:
                                mirror_type = "d"
                                for v, r in self.rot_sym:
                                    if np.linalg.norm(v - axis) >= self.tol:
                                        if np.dot(v, normal) < self.tol:
                                            mirror_type = "v"
                                            break
                            else:
                                mirror_type = "v"
                            break

        return mirror_type

    def _get_smallest_set_not_on_axis(self, axis):
        """
        Returns the smallest list of atoms with the same species and
        distance from origin AND does not lie on the specified axis. This
        maximal set limits the possible rotational symmetry operations,
        since atoms lying on a test axis is irrelevant in testing rotational
        symmetryOperations.
        """

        def not_on_axis(site):
            v = np.cross(site.coords, axis)
            return np.linalg.norm(v) > self.tol

        valid_sets = []
        origin_site, dist_el_sites = cluster_sites(self.centered_mol, self.tol)
        for test_set in dist_el_sites.values():
            valid_set = list(filter(not_on_axis, test_set))
            if len(valid_set) > 0:
                valid_sets.append(valid_set)

        return min(valid_sets, key=lambda s: len(s))

    def _check_rot_sym(self, axis):
        """
        Determines the rotational symmetry about supplied axis. Used only for
        symmetric top molecules which has possible rotational symmetry
        operations > 2.
        """
        min_set = self._get_smallest_set_not_on_axis(axis)
        max_sym = len(min_set)
        for i in range(max_sym, 0, -1):
            if max_sym % i != 0:
                continue
            op = SymmOp.from_axis_angle_and_translation(axis, 360 / i)
            rotvalid = self.is_valid_op(op)
            if rotvalid:
                self.symmops.append(op)
                self.rot_sym.append((axis, i))
                return i
        return 1

    def _check_perpendicular_r2_axis(self, axis):
        """
        Checks for R2 axes perpendicular to unique axis. For handling
        symmetric top molecules.
        """
        min_set = self._get_smallest_set_not_on_axis(axis)
        for s1, s2 in itertools.combinations(min_set, 2):
            test_axis = np.cross(s1.coords - s2.coords, axis)
            if np.linalg.norm(test_axis) > self.tol:
                op = SymmOp.from_axis_angle_and_translation(test_axis, 180)
                r2present = self.is_valid_op(op)
                if r2present:
                    self.symmops.append(op)
                    self.rot_sym.append((test_axis, 2))
                    return True
        return None

    def _proc_sph_top(self):
        """
        Handles Sperhical Top Molecules, which belongs to the T, O or I point
        groups.
        """
        self._find_spherical_axes()
        if len(self.rot_sym) == 0:
            logger.debug("Accidental speherical top!")
            self._proc_sym_top()
        main_axis, rot = max(self.rot_sym, key=lambda v: v[1])
        if rot < 3:
            logger.debug("Accidental speherical top!")
            self._proc_sym_top()
        elif rot == 3:
            mirror_type = self._find_mirror(main_axis)
            if mirror_type != "":
                if self.is_valid_op(PointGroupAnalyzer.inversion_op):
                    self.symmops.append(PointGroupAnalyzer.inversion_op)
                    self.sch_symbol = "Th"
                else:
                    self.sch_symbol = "Td"
            else:
                self.sch_symbol = "T"
        elif rot == 4:
            if self.is_valid_op(PointGroupAnalyzer.inversion_op):
                self.symmops.append(PointGroupAnalyzer.inversion_op)
                self.sch_symbol = "Oh"
            else:
                self.sch_symbol = "O"
        elif rot == 5:
            if self.is_valid_op(PointGroupAnalyzer.inversion_op):
                self.symmops.append(PointGroupAnalyzer.inversion_op)
                self.sch_symbol = "Ih"
            else:
                self.sch_symbol = "I"

    def _find_spherical_axes(self):
        """
        Looks for R5, R4, R3 and R2 axes in spherical top molecules. Point
        group T molecules have only one unique 3-fold and one unique 2-fold
        axis. O molecules have one unique 4, 3 and 2-fold axes. I molecules
        have a unique 5-fold axis.
        """
        rot_present = defaultdict(bool)
        origin_site, dist_el_sites = cluster_sites(self.centered_mol, self.tol)
        test_set = min(dist_el_sites.values(), key=lambda s: len(s))
        coords = [s.coords for s in test_set]
        for c1, c2, c3 in itertools.combinations(coords, 3):
            for cc1, cc2 in itertools.combinations([c1, c2, c3], 2):
                if not rot_present[2]:
                    test_axis = cc1 + cc2
                    if np.linalg.norm(test_axis) > self.tol:
                        op = SymmOp.from_axis_angle_and_translation(test_axis, 180)
                        rot_present[2] = self.is_valid_op(op)
                        if rot_present[2]:
                            self.symmops.append(op)
                            self.rot_sym.append((test_axis, 2))

            test_axis = np.cross(c2 - c1, c3 - c1)
            if np.linalg.norm(test_axis) > self.tol:
                for r in (3, 4, 5):
                    if not rot_present[r]:
                        op = SymmOp.from_axis_angle_and_translation(test_axis, 360 / r)
                        rot_present[r] = self.is_valid_op(op)
                        if rot_present[r]:
                            self.symmops.append(op)
                            self.rot_sym.append((test_axis, r))
                            break
            if rot_present[2] and rot_present[3] and (rot_present[4] or rot_present[5]):
                break

    def get_pointgroup(self):
        """
        Returns a PointGroup object for the molecule.
        """
        return PointGroupOperations(self.sch_symbol, self.symmops, self.mat_tol)

    def get_symmetry_operations(self):
        """
        Return symmetry operations as a list of SymmOp objects.
        Returns Cartesian coord symmops.

        Returns:
            ([SymmOp]): List of symmetry operations.
        """
        return generate_full_symmops(self.symmops, self.tol)

    def get_rotational_symmetry_number(self):
        """
        Return the rotational symmetry number.
        """
        symm_ops = self.get_symmetry_operations()
        symm_number = 0
        for symm in symm_ops:
            rot = symm.rotation_matrix
            if np.abs(np.linalg.det(rot) - 1) < 1e-4:
                symm_number += 1
        return symm_number

    def is_valid_op(self, symmop):
        """
        Check if a particular symmetry operation is a valid symmetry operation
        for a molecule, i.e., the operation maps all atoms to another
        equivalent atom.

        Args:
            symmop (SymmOp): Symmetry operation to test.

        Returns:
            (bool): Whether SymmOp is valid for Molecule.
        """
        coords = self.centered_mol.cart_coords
        for site in self.centered_mol:
            coord = symmop.operate(site.coords)
            ind = find_in_coord_list(coords, coord, self.tol)
            if not (len(ind) == 1 and self.centered_mol[ind[0]].species == site.species):
                return False
        return True

    def _get_eq_sets(self):
        """
        Calculates the dictionary for mapping equivalent atoms onto each other.

        Args:
            None

        Returns:
            dict: The returned dictionary has two possible keys:

            ``eq_sets``:
            A dictionary of indices mapping to sets of indices,
            each key maps to indices of all equivalent atoms.
            The keys are guaranteed to be not equivalent.

            ``sym_ops``:
            Twofold nested dictionary.
            ``operations[i][j]`` gives the symmetry operation
            that maps atom ``i`` unto ``j``.
        """
        UNIT = np.eye(3)
        eq_sets, operations = defaultdict(set), defaultdict(dict)
        symm_ops = [op.rotation_matrix for op in generate_full_symmops(self.symmops, self.tol)]

        def get_clustered_indices():
            indices = cluster_sites(self.centered_mol, self.tol, give_only_index=True)
            out = list(indices[1].values())
            if indices[0] is not None:
                out.append([indices[0]])
            return out

        for index in get_clustered_indices():
            sites = self.centered_mol.cart_coords[index]
            for i, reference in zip(index, sites):
                for op in symm_ops:
                    rotated = np.dot(op, sites.T).T
                    matched_indices = find_in_coord_list(rotated, reference, self.tol)
                    matched_indices = {dict(enumerate(index))[i] for i in matched_indices}
                    eq_sets[i] |= matched_indices

                    if i not in operations:
                        operations[i] = {j: op.T if j != i else UNIT for j in matched_indices}
                    else:
                        for j in matched_indices:
                            if j not in operations[i]:
                                operations[i][j] = op.T if j != i else UNIT
                    for j in matched_indices:
                        if j not in operations:
                            operations[j] = {i: op if j != i else UNIT}
                        elif i not in operations[j]:
                            operations[j][i] = op if j != i else UNIT

        return {"eq_sets": eq_sets, "sym_ops": operations}

    @staticmethod
    def _combine_eq_sets(eq_sets, operations):
        """Combines the dicts of _get_equivalent_atom_dicts into one

        Args:
            eq_sets (dict)
            operations (dict)

        Returns:
            dict: The returned dictionary has two possible keys:

            ``eq_sets``:
            A dictionary of indices mapping to sets of indices,
            each key maps to indices of all equivalent atoms.
            The keys are guaranteed to be not equivalent.

            ``sym_ops``:
            Twofold nested dictionary.
            ``operations[i][j]`` gives the symmetry operation
            that maps atom ``i`` unto ``j``.
        """
        UNIT = np.eye(3)

        def all_equivalent_atoms_of_i(i, eq_sets, ops):
            """WORKS INPLACE on operations"""
            visited = {i}
            tmp_eq_sets = {j: (eq_sets[j] - visited) for j in eq_sets[i]}

            while tmp_eq_sets:
                new_tmp_eq_sets = {}
                for j in tmp_eq_sets:
                    if j in visited:
                        continue
                    visited.add(j)
                    for k in tmp_eq_sets[j]:
                        new_tmp_eq_sets[k] = eq_sets[k] - visited
                        if i not in ops[k]:
                            ops[k][i] = np.dot(ops[j][i], ops[k][j]) if k != i else UNIT
                        ops[i][k] = ops[k][i].T
                tmp_eq_sets = new_tmp_eq_sets
            return visited, ops

        eq_sets = copy.deepcopy(eq_sets)
        ops = copy.deepcopy(operations)
        to_be_deleted = set()
        for i in eq_sets:
            if i in to_be_deleted:
                continue
            visited, ops = all_equivalent_atoms_of_i(i, eq_sets, ops)
            to_be_deleted |= visited - {i}

        for k in to_be_deleted:
            eq_sets.pop(k, None)
        return {"eq_sets": eq_sets, "sym_ops": ops}

    def get_equivalent_atoms(self):
        """Returns sets of equivalent atoms with symmetry operations

        Args:
            None

        Returns:
            dict: The returned dictionary has two possible keys:

            ``eq_sets``:
            A dictionary of indices mapping to sets of indices,
            each key maps to indices of all equivalent atoms.
            The keys are guaranteed to be not equivalent.

            ``sym_ops``:
            Twofold nested dictionary.
            ``operations[i][j]`` gives the symmetry operation
            that maps atom ``i`` unto ``j``.
        """
        eq = self._get_eq_sets()
        return self._combine_eq_sets(eq["eq_sets"], eq["sym_ops"])

    def symmetrize_molecule(self):
        """Returns a symmetrized molecule

        The equivalent atoms obtained via
        :meth:`~pymatgen.symmetry.analyzer.PointGroupAnalyzer.get_equivalent_atoms`
        are rotated, mirrored... unto one position.
        Then the average position is calculated.
        The average position is rotated, mirrored... back with the inverse
        of the previous symmetry operations, which gives the
        symmetrized molecule

        Args:
            None

        Returns:
            dict: The returned dictionary has three possible keys:

            ``sym_mol``:
            A symmetrized molecule instance.

            ``eq_sets``:
            A dictionary of indices mapping to sets of indices,
            each key maps to indices of all equivalent atoms.
            The keys are guaranteed to be not equivalent.

            ``sym_ops``:
            Twofold nested dictionary.
            ``operations[i][j]`` gives the symmetry operation
            that maps atom ``i`` unto ``j``.
        """
        eq = self.get_equivalent_atoms()
        eq_sets, ops = eq["eq_sets"], eq["sym_ops"]
        coords = self.centered_mol.cart_coords.copy()
        for i, eq_indices in eq_sets.items():
            for j in eq_indices:
                coords[j] = np.dot(ops[j][i], coords[j])
            coords[i] = np.mean(coords[list(eq_indices)], axis=0)
            for j in eq_indices:
                if j == i:
                    continue
                coords[j] = np.dot(ops[i][j], coords[i])
                coords[j] = np.dot(ops[i][j], coords[i])
        molecule = Molecule(species=self.centered_mol.species_and_occu, coords=coords)
        return {"sym_mol": molecule, "eq_sets": eq_sets, "sym_ops": ops}


def iterative_symmetrize(mol, max_n=10, tolerance=0.3, epsilon=1e-2):
    """Returns a symmetrized molecule

    The equivalent atoms obtained via
    :meth:`~pymatgen.symmetry.analyzer.PointGroupAnalyzer.get_equivalent_atoms`
    are rotated, mirrored... unto one position.
    Then the average position is calculated.
    The average position is rotated, mirrored... back with the inverse
    of the previous symmetry operations, which gives the
    symmetrized molecule

    Args:
        mol (Molecule): A pymatgen Molecule instance.
        max_n (int): Maximum number of iterations.
        tolerance (float): Tolerance for detecting symmetry.
            Gets passed as Argument into
            :class:`~pymatgen.analyzer.symmetry.PointGroupAnalyzer`.
        epsilon (float): If the elementwise absolute difference of two
            subsequently symmetrized structures is smaller epsilon,
            the iteration stops before ``max_n`` is reached.


    Returns:
        dict: The returned dictionary has three possible keys:

        ``sym_mol``:
        A symmetrized molecule instance.

        ``eq_sets``:
        A dictionary of indices mapping to sets of indices,
        each key maps to indices of all equivalent atoms.
        The keys are guaranteed to be not equivalent.

        ``sym_ops``:
        Twofold nested dictionary.
        ``operations[i][j]`` gives the symmetry operation
        that maps atom ``i`` unto ``j``.
    """
    new = mol
    n = 0
    finished = False
    while not finished and n <= max_n:
        previous = new
        PA = PointGroupAnalyzer(previous, tolerance=tolerance)
        eq = PA.symmetrize_molecule()
        new = eq["sym_mol"]
        finished = np.allclose(new.cart_coords, previous.cart_coords, atol=epsilon)
        n += 1
    return eq


def cluster_sites(mol, tol, give_only_index=False):
    """
    Cluster sites based on distance and species type.

    Args:
        mol (Molecule): Molecule **with origin at center of mass**.
        tol (float): Tolerance to use.

    Returns:
        (origin_site, clustered_sites): origin_site is a site at the center
        of mass (None if there are no origin atoms). clustered_sites is a
        dict of {(avg_dist, species_and_occu): [list of sites]}
    """
    # Cluster works for dim > 2 data. We just add a dummy 0 for second
    # coordinate.
    dists = [[np.linalg.norm(site.coords), 0] for site in mol]
    import scipy.cluster as spcluster

    f = spcluster.hierarchy.fclusterdata(dists, tol, criterion="distance")
    clustered_dists = defaultdict(list)
    for i, site in enumerate(mol):
        clustered_dists[f[i]].append(dists[i])
    avg_dist = {label: np.mean(val) for label, val in clustered_dists.items()}
    clustered_sites = defaultdict(list)
    origin_site = None
    for i, site in enumerate(mol):
        if avg_dist[f[i]] < tol:
            if give_only_index:
                origin_site = i
            else:
                origin_site = site
        else:
            if give_only_index:
                clustered_sites[(avg_dist[f[i]], site.species)].append(i)
            else:
                clustered_sites[(avg_dist[f[i]], site.species)].append(site)
    return origin_site, clustered_sites


def generate_full_symmops(symmops, tol):
    """
    Recursive algorithm to permute through all possible combinations of the
    initially supplied symmetry operations to arrive at a complete set of
    operations mapping a single atom to all other equivalent atoms in the
    point group. This assumes that the initial number already uniquely
    identifies all operations.

    Args:
        symmops ([SymmOp]): Initial set of symmetry operations.

    Returns:
        Full set of symmetry operations.
    """
    # Uses an algorithm described in:
    # Gregory Butler. Fundamental Algorithms for Permutation Groups.
    # Lecture Notes in Computer Science (Book 559). Springer, 1991. page 15
    UNIT = np.eye(4)
    generators = [op.affine_matrix for op in symmops if not np.allclose(op.affine_matrix, UNIT)]
    if not generators:
        # C1 symmetry breaks assumptions in the algorithm afterwards
        return symmops

    full = list(generators)

    for g in full:
        for s in generators:
            op = np.dot(g, s)
            d = np.abs(full - op) < tol
            if not np.any(np.all(np.all(d, axis=2), axis=1)):
                full.append(op)
            if len(full) > 1000:
                warnings.warn(
                    f"{len(full)} matrices have been generated. The tol may be too small. Please terminate"
                    f" and rerun with a different tolerance."
                )

    d = np.abs(full - UNIT) < tol
    if not np.any(np.all(np.all(d, axis=2), axis=1)):
        full.append(UNIT)
    return [SymmOp(op) for op in full]


class SpacegroupOperations(list):
    """
    Represents a space group, which is a collection of symmetry operations.
    """

    def __init__(self, int_symbol, int_number, symmops):
        """
        Args:
            int_symbol (str): International symbol of the spacegroup.
            int_number (int): International number of the spacegroup.
            symmops ([SymmOp]): Symmetry operations associated with the
                spacegroup.
        """
        self.int_symbol = int_symbol
        self.int_number = int_number
        super().__init__(symmops)

    def are_symmetrically_equivalent(self, sites1, sites2, symm_prec=1e-3):
        """
        Given two sets of PeriodicSites, test if they are actually
        symmetrically equivalent under this space group. Useful, for example,
        if you want to test if selecting atoms 1 and 2 out of a set of 4 atoms
        are symmetrically the same as selecting atoms 3 and 4, etc.

        One use is in PartialRemoveSpecie transformation to return only
        symmetrically distinct arrangements of atoms.

        Args:
            sites1 ([PeriodicSite]): 1st set of sites
            sites2 ([PeriodicSite]): 2nd set of sites
            symm_prec (float): Tolerance in atomic distance to test if atoms
                are symmetrically similar.

        Returns:
            (bool): Whether the two sets of sites are symmetrically
            equivalent.
        """

        def in_sites(site):
            for test_site in sites1:
                if test_site.is_periodic_image(site, symm_prec, False):
                    return True
            return False

        for op in self:
            newsites2 = [PeriodicSite(site.species, op.operate(site.frac_coords), site.lattice) for site in sites2]
            for site in newsites2:
                if not in_sites(site):
                    break
            else:
                return True
        return False

    def __str__(self):
        return f"{self.int_symbol} ({self.int_number}) spacegroup"


class PointGroupOperations(list):
    """
    Defines a point group, which is essentially a sequence of symmetry
    operations.

    .. attribute:: sch_symbol

        Schoenflies symbol of the point group.
    """

    def __init__(self, sch_symbol, operations, tol=0.1):
        """
        Args:
            sch_symbol (str): Schoenflies symbol of the point group.
            operations ([SymmOp]): Initial set of symmetry operations. It is
                sufficient to provide only just enough operations to generate
                the full set of symmetries.
            tol (float): Tolerance to generate the full set of symmetry
                operations.
        """
        self.sch_symbol = sch_symbol
        super().__init__(generate_full_symmops(operations, tol))

    def __str__(self):
        return self.sch_symbol

    def __repr__(self):
        return self.__str__()


class MagneticSpacegroupAnalyzer:
    """
    Takes a pymatgen.core.structure.Structure object, symprec and angle_tolarance.

    Determines magnetic spacegroup operations and magnetic spacegroup. 
    Depends on pymatgen.symmetry.analyzer.SpacegroupAnalyzer, which interfaces
    with spglib.

    Args:
        structure (Structure/IStructure): Structure to find symmetry
        symprec (float): Tolerance for symmetry finding of SpacegroupAnalyzer. 
            Defaults to 0.01, which is fairly strict and works well for properly 
            refined structures with atoms in the proper symmetry coordinates. For
            structures with slight deviations from their proper atomic
            positions (e.g., structures relaxed with electronic structure
            codes), a looser tolerance of 0.1 (the value used in Materials
            Project) is often needed.
        angle_tolerance (float): Angle tolerance for symmetry finding.
        pos_tolerance (float): Tolerance for fractional coordinates in magnetic 
            structures. Defaults to 0.0015.
        magmom_tolerance (float): Tolerance for the magnetic moment. 
            Defaults to 0.3, which is about the uncertainty one would get in an
            ab initio calculation.
    """

    def __init__(self, structure, symprec=0.01, angle_tolerance=5, pos_tolerance=0.0015, magmom_tolerance=0.3):
        # TODO: structure code such that the circular dependencies btw pymatgen.transformations 
        #       and pymatgen.symmetry are resolved
      
        from pymatgen.transformations.standard_transformations import ApplyMagSymmOpTransformation
        from pymatgen.alchemy.materials import TransformedStructure
      
        self._structure = structure
        self._symprec = symprec
        self._angle_tol = angle_tolerance
        self._pos_tol = pos_tolerance
        self._magmom_tol = magmom_tolerance
        self._domains = None
      
        self._spacegroupAnalyzer = SpacegroupAnalyzer(self._structure, 
                                              symprec = self._symprec, 
                                              angle_tolerance = self._angle_tol
                                              )
      
        if np.any(np.abs(self._structure._lattice.matrix - self._spacegroupAnalyzer._cell[0]) > 1e-12): 
            exit("Error in the definitions of the unit cell.")

        self._paramag_symmops = self._paramag_symmetry_operations()

        # apply all paramagnetic spacegroup operations and check if they are a symmetry
        mspg_list = [mop for mop in self._paramag_symmops if self._check_mop(mop)]

        self.symmetry_ops = mspg_list
        self._MagneticSpacegroup = self._search4MagneticSpacegroup()

    def _check_mop(self, mop):
        from pymatgen.transformations.standard_transformations import ApplyMagSymmOpTransformation
        from pymatgen.alchemy.materials import TransformedStructure

        trans = [ApplyMagSymmOpTransformation(mop)]
        s_new = TransformedStructure(self._structure, trans)
        same = self._same_MagneticStructure(s_new.final_structure,  self._structure)[0]
        return same

    def _get_domains(self):
        from pymatgen.transformations.standard_transformations import ApplyMagSymmOpTransformation
        from pymatgen.alchemy.materials import TransformedStructure

        domains = {MagSymmOp.from_xyzt_string('x, y, z, +1') : self._structure}
        relevant_mops = [mop for mop in self._paramag_symmops if not np.any([mop == mop2 for mop2 in self.symmetry_ops])]
      
        for mop in relevant_mops:
            trans = [ApplyMagSymmOpTransformation(mop)]
            s_new = TransformedStructure(self._structure, trans)
            new_domain = True
            for do in domains.values():
                if self._same_MagneticStructure(s_new.final_structure, do)[0]:
                    new_domain = False
                    break
            if new_domain: 
                domains[mop] = s_new.final_structure
        return domains

    def get_parent_spacegroup(self):
        """
        Returns analyzer of the parent spacegroup, i.e. the spacegroup of 
        the nonmagnetic structure.

        Returns:
            pymatgen.symmetry.analyzer.SpacegroupAnalyzer
        """
        return self._spacegroupAnalyzer

    def get_paramagnetic_symmetry_operations(self):
        """
        Returns all magnetic symmetry operations of the nonmagnetic structure, which comprises all
        crystallographic spacegroup operations with time reversal symmetry even and odd.
      
        Returns:
            list(pymatgen.core.operations.MagSymmOp): List of symmetry operations.
        """
        return self._paramag_symmetry_ops

    def _paramag_symmetry_operations(self, cartesian=False):
        """
        Computes all magnetic symmetry operations of the nonmagnetic structure, which comprises all
        crystallographic spacegroup operations with time reversal symmetry even and odd.
        By default returns symmetry operations in fractional coordinates for a given structure,
        but can also return symmetry operations in Cartesian coordinates.

        Returns:
            list(pymatgen.core.operations.MagSymmOp): List of symmetry operations.
        """
        if np.any(np.abs(self._structure._lattice.matrix - self._spacegroupAnalyzer._cell[0]) > 1e-12):
            exit("Error in the definitions of the unit cell.")

        rotation = self._spacegroupAnalyzer._space_group_data['rotations']
        translation = self._spacegroupAnalyzer._space_group_data['translations']
        
        # if cell is conventional prim2cart will be conv2cart
        prim2cart = np.transpose(self._spacegroupAnalyzer._cell[0])
        prim2cart_inv = np.linalg.inv(prim2cart)
        # if cell is already conventional prim2conv will be unit matrix
        prim2conv = self._spacegroupAnalyzer._space_group_data["transformation_matrix"]
        prim2conv_inv = np.linalg.inv(prim2conv)

        symmops = []
        for rot, trans in zip(rotation, translation):
            if cartesian:
                rot = prim2cart.dot(rot.dot(prim2cart_inv))
                trans = prim2cart.dot(trans)
            else:
                rot = prim2conv.dot(rot.dot(prim2conv_inv)).round(8)
                trans = prim2conv.dot(trans).round(8)
          
            op = SymmOp.from_rotation_and_translation(rot, trans)
            symmops.append(op)
        
        paramag_symmops = [MagSymmOp.from_symmop(op, time_reversal) for op in symmops for time_reversal in [+1, -1]]
        
        return paramag_symmops

    def get_magnetic_pointgroup(self):
        """
        Get magnetic pointgroup information as tuple (label, number).

        Returns:
            tuple(str, list(int)): Magnetic point group for structure.
        """
        MAGSPACEGROUP2MAGPOINTGROUP = os.path.join(os.path.dirname(__file__), "mspg2mpg_map.txt")
      
        glob_og_no = ".".join([str(x) for x in self._MagneticSpacegroup._data["og_number"]])

        with open(MAGSPACEGROUP2MAGPOINTGROUP, 'r') as f:
            for line in f:
                og_no, og_label, pointgroup_no, pointgroup_label = line.split()
                if og_no == glob_og_no:
                    return (pointgroup_label, pointgroup_no)

    def get_magnetic_spacegroup(self):
        """
        Get the magnetic spacegroup.

        Returns:
            pymatgen.symmetry.maggroups.MagneticSpaceGroup 
        """
        return self._MagneticSpacegroup

    def get_domains(self):
        """
        Get all possible magnetic domains as a dictionary.

        A magnetic domain is obtained when a symmetry operation of the crystallographic 
        spacegroup is applied to a magnetic structure with either time reversal odd or 
        even and the resulting magnetic structure is not covering. 

        The keys() are the applied magnetic symmetry operations. 
        The values() are the resulting magnetic structure.

        Returns:
            dict(pymatgen.core.operations.MagSymmOp : pymatgen.core.structure.Structure) 
        """
      
        if not self._domains:
            self._domains = self._get_domains()
      
        return self._domains

    def get_symmetry_ops(self):
        """
        Get magnetic spacegroup operations for the magnetic structure.

        Returns:
            list(pymatgen.core.operations.MagSymmOp)
        """ 
        return self._MagneticSpacegroup.symmetry_ops

    def get_bns(self):
        """
        Get the magnetic spacegroup in the BNS notation [1], e.g. (P1, 1.1). 

        [1] Belov et al. (1957). Sov. Phys. Crystallogr. 2, 311-322 
      
        Returns:
            (str, list(int)): Magnetic space group in a tuple of (symbol, number) for structure.
        """
        bns_label = self._MagneticSpacegroup._data["bns_label"]
        bns_number = ".".join([str(x) for x in self._MagneticSpacegroup._data["bns_number"]])
        return (bns_label, bns_number)  

    def get_og(self):
        """
        Get the magnetic spacegroup in the OG notation [1], e.g. (P1, 1.1.1). 

        [1] Opechowski and Guccione (1965). Magnetism, edited by G. T. Rado and H. Suhl, Vol. II, 
            Part A, pp. 105-165. New York: Academic Press
      
        Returns:
            (str, list(int)): Magnetic space group in a tuple of (symbol, number) for structure.
        """
        og_label = self._MagneticSpacegroup._data["og_label"]
        og_number = ".".join([str(x) for x in self._MagneticSpacegroup._data["og_number"]])
        return (og_label, og_number)  

    def _same_MagneticStructure(self, struc1, struc2, pos_tolerance=0.0015, magmom_tolerance=0.3):
        """
        Compare two pymatgen.core.structure.Structure objects, struc1 and struc2,
        with respect to their sitesi, species and magnetic moment.

        Takes:
          struc1           pymatgen.core.structure.Structure 
          struc1           pymatgen.core.structure.Structure 
          pos_tolerance    float (default: 0.0015)
          magmom_tolerance float (default: 0.3)

        Returns:
          boolean          True if struc1 covers struc2
          set(tuple(tuple(str), tuple(float, float, float), tuple(float, float, float))) 
                         Residue: Set of species, fractional coordinates 
                         and magnetic moments that are distinct. 
                         Empty set if boolean is True.
        """
        same = False
        nsite = len(struc1.sites)
        if nsite != len(struc2.sites):
            res = "number of sites differ"
            pass
        else:
            sites1_prim = [np.mod(i.frac_coords + 10 + 1e-9, 1) - 1e-9 for i in struc1._sites]
            sites1_prim = [tuple(k) for k in sites1_prim]
            sites2_prim = [np.mod(i.frac_coords + 10 + 1e-9, 1) - 1e-9 for i in struc2._sites]
            sites2_prim = [tuple(k) for k in sites2_prim]
        
            magms1_cart = [i.properties["magmom"].moment for i in struc1._sites]
            magms1_cart = [tuple(k.round(3)) for k in magms1_cart]
            magms2_cart = [i.properties["magmom"].moment for i in struc2._sites]
            magms2_cart = [tuple(k.round(3)) for k in magms2_cart]
            species1 = [i._species for i in struc1._sites]
            species2 = [i._species for i in struc2._sites]
      
            set1 = set(tuple([tuple([species1[i],sites1_prim[i],magms1_cart[i]]) for i in range(nsite)]))
            set2 = set(tuple([tuple([species2[i],sites2_prim[i],magms2_cart[i]]) for i in range(nsite)]))
            res = self._difference(set1, set2, pos_tolerance=pos_tolerance, magmom_tolerance=magmom_tolerance)
            if res == set():
                same = True
      
        return same, res

    def _difference(self, sitesmagms_set, sitesmagms_set_i, pos_tolerance=0.0015, magmom_tolerance=0.3):
        """
        Returns the difference between two sets of sites representing magnetic structures.

        Takes:
            sitesmagms_set1  set(tuple(tuple(str), tuple(float, float, float), tuple(float, float, float)))
            sitesmagms_set2  set(tuple(tuple(str), tuple(float, float, float), tuple(float, float, float)))
            pos_tolerance    float (default: 0.0015)
            magmom_tolerance float (default: 0.3)

        Returns:
            set(tuple(tuple(str), tuple(float, float, float), tuple(float, float, float)))
        """
        res_sitesmagms_i = sitesmagms_set_i - sitesmagms_set
        res_sitesmagms_new = sitesmagms_set - sitesmagms_set_i
        res = []
        for t1 in res_sitesmagms_new:
            t1_equal2_t2 = False
            for t2 in res_sitesmagms_i:
                if t1[0] == t2[0]:
                    if np.all(np.mod(np.abs(np.array(t1[1]) - np.array(t2[1])), 1) <= pos_tolerance):
                        if np.all(np.abs(np.array(t1[2]) - np.array(t2[2])) <= magmom_tolerance):
                            t1_equal2_t2 = True
            if not t1_equal2_t2: 
                res.append(t1)
        return set(res) 

    def _search4MagneticSpacegroup(self):
        """
        Returns the magnetic spacegroup of the magnetic structure.
        Compares the magnetic spacegroup operations of the magnetic structure
        to the magnetic spacegroup operations of all magnetic spacegroups that
        have the same crystal system and number of operations.

        Returns:
            pymatgen.symmetry.MagneticSpaceGroup
        """
      
        mop_list = self.symmetry_ops

        # get all potential magnetic spacegroups 
        # according to the crystal system
        from pymatgen.symmetry.maggroups import MagneticSpaceGroup
        cs = self._spacegroupAnalyzer.get_crystal_system()

        cs_range = {
            "triclinic": (1, 8),
            "monoclinic": (8, 99),
            "orthorhombic": (99, 661),
            "tetragonal": (661, 1231),
            "trigonal": (1231, 1339),
            "hexagonal": (1339, 1503),
            "cubic": (1503, 1652),
        }
      
        # start1 = time.time()    
        # get all potential magnetic spacegroups 
        # according to the order of the magnetic spacegroup
        ORDER2MSPG = {
                1 : [1],
                2 : [2, 3, 4, 6, 8, 10, 15, 17, 25, 27, 32, 34],
                3 : [1231, 1234, 1237],
                4 : [5, 7, 9, 11, 12, 14, 16, 18, 19, 21, 23, 24, 26, 28, 29, 
                        31, 33, 35, 36, 38, 40, 42, 44, 45, 47, 48, 49, 51, 52, 
                        53, 59, 61, 62, 63, 77, 79, 80, 81, 86, 88, 89, 90, 99, 
                        101, 106, 108, 109, 113, 115, 116, 119, 121, 155, 157, 
                        158, 168, 170, 171, 172, 178, 180, 181, 185, 187, 188, 
                        189, 198, 200, 201, 202, 205, 207, 208, 209, 212, 214, 
                        215, 216, 219, 221, 222, 226, 228, 229, 230, 231, 233, 
                        234, 661, 663, 668, 670, 672, 674, 679, 681, 693, 695
                        ],
                6 : [1232, 1233, 1235, 1236, 1238, 1239, 1243, 1245, 1251, 1253, 
                        1255, 1257, 1259, 1261, 1263, 1265, 1267, 1269, 1271, 1273, 
                        1279, 1281, 1284, 1286, 1289, 1291, 1292, 1294, 1339, 1341, 
                        1344, 1346, 1347, 1349, 1350, 1352, 1355, 1357, 1360, 1362, 
                        1363, 1365
                        ],
                8 : [13, 20, 22, 30, 37, 39, 41, 43, 46, 50, 54, 55, 57, 58, 60, 64, 
                        65, 66, 68, 69, 70, 72, 74, 75, 76, 78, 82, 83, 85, 87, 91, 
                        92, 94, 95, 96, 97, 98, 100, 102, 105, 107, 110, 112, 114, 117, 
                        118, 120, 122, 124, 125, 126, 127, 128, 129, 131, 132, 134, 
                        137, 138, 145, 147, 148, 149, 150, 152, 153, 154, 156, 159, 
                        160, 164, 165, 166, 169, 173, 174, 176, 177, 179, 182, 184, 
                        186, 190, 191, 193, 194, 195, 196, 199, 203, 204, 206, 210, 
                        211, 213, 217, 218, 220, 223, 224, 225, 227, 232, 236, 238, 
                        239, 241, 245, 246, 249, 251, 252, 253, 254, 255, 256, 257, 
                        258, 260, 261, 262, 263, 264, 265, 267, 268, 269, 271, 274, 
                        275, 276, 278, 280, 281, 282, 284, 287, 288, 289, 291, 293, 
                        294, 295, 296, 297, 298, 299, 300, 302, 303, 304, 305, 306, 
                        307, 308, 324, 326, 327, 328, 329, 330, 331, 333, 334, 335, 
                        336, 337, 338, 340, 341, 342, 343, 344, 345, 346, 347, 349, 
                        350, 351, 358, 360, 361, 362, 364, 366, 367, 368, 369, 370, 
                        377, 379, 380, 381, 382, 383, 387, 389, 390, 391, 392, 393, 
                        394, 395, 406, 408, 409, 410, 411, 412, 413, 414, 415, 417, 
                        418, 419, 420, 421, 422, 423, 428, 430, 431, 432, 433, 434, 
                        435, 436, 441, 443, 444, 445, 446, 447, 451, 453, 454, 455, 
                        456, 457, 458, 460, 461, 462, 463, 464, 465, 466, 471, 473, 
                        474, 475, 476, 477, 478, 480, 481, 482, 483, 484, 488, 490, 
                        491, 492, 493, 494, 495, 496, 497, 499, 500, 501, 502, 504, 
                        505, 506, 507, 508, 509, 510, 662, 664, 665, 667, 669, 671, 
                        673, 675, 676, 678, 680, 682, 683, 685, 686, 687, 688, 690, 
                        691, 692, 694, 696, 697, 699, 701, 702, 703, 705, 706, 707, 
                        713, 715, 716, 717, 720, 722, 723, 724, 727, 729, 730, 731, 
                        747, 749, 750, 751, 757, 759, 760, 761, 764, 766, 767, 768, 
                        771, 773, 774, 775, 776, 778, 779, 780, 786, 788, 789, 790, 
                        793, 795, 796, 797, 800, 802, 803, 804, 823, 825, 826, 827, 
                        836, 838, 839, 840, 845, 847, 848, 849, 852, 854, 855, 856, 
                        859, 861, 862, 863, 866, 868, 869, 870, 871, 873, 874, 875, 
                        878, 880, 881, 882, 911, 913, 914, 915, 922, 924, 925, 926, 
                        929, 931, 932, 933, 936, 938, 939, 940, 941, 943, 944, 945, 
                        951, 953, 954, 955, 958, 960, 961, 962, 965, 967, 968, 969
                        ],
                9 : [1240],
                12 : [1244, 1246, 1252, 1254, 1256, 1258, 1260, 1262, 1264, 1266, 1268, 
                        1270, 1272, 1274, 1280, 1282, 1283, 1285, 1287, 1288, 1290, 1293, 
                        1303, 1305, 1306, 1307, 1310, 1312, 1313, 1314, 1315, 1317, 1318, 
                        1319, 1322, 1324, 1325, 1326, 1340, 1342, 1343, 1345, 1348, 1351, 
                        1353, 1354, 1356, 1358, 1359, 1361, 1364, 1366, 1367, 1369, 1370, 
                        1371, 1374, 1376, 1377, 1378, 1379, 1381, 1382, 1383, 1386, 1388, 
                        1389, 1390, 1391, 1393, 1394, 1395, 1396, 1398, 1399, 1400, 1403, 
                        1405, 1406, 1407, 1410, 1412, 1413, 1414, 1415, 1417, 1418, 1419, 
                        1424, 1426, 1427, 1428, 1429, 1431, 1432, 1433, 1434, 1436, 1437, 
                        1438, 1439, 1441, 1442, 1443, 1446, 1448, 1449, 1450, 1451, 1453, 
                        1454, 1455, 1458, 1460, 1461, 1462, 1503, 1511
                        ],
                16 : [56, 67, 71, 73, 84, 93, 103, 111, 123, 130, 133, 135, 136, 139, 140, 
                        142, 143, 144, 146, 151, 161, 162, 167, 175, 183, 192, 197, 237, 
                        240, 242, 243, 244, 247, 248, 250, 259, 266, 270, 272, 273, 277, 
                        279, 283, 285, 286, 290, 292, 301, 309, 311, 312, 313, 314, 315, 
                        316, 317, 318, 319, 320, 322, 323, 325, 332, 339, 348, 352, 355, 
                        356, 359, 365, 371, 373, 374, 375, 378, 384, 385, 386, 388, 396, 
                        397, 399, 400, 401, 402, 403, 404, 407, 416, 424, 425, 426, 427, 
                        429, 437, 438, 439, 440, 442, 448, 449, 450, 452, 459, 467, 468, 
                        469, 470, 472, 479, 485, 486, 487, 489, 498, 503, 511, 513, 514, 
                        515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 
                        528, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 
                        542, 543, 544, 545, 547, 548, 549, 550, 551, 553, 557, 558, 559, 
                        560, 561, 564, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 
                        576, 577, 579, 580, 581, 582, 583, 585, 589, 590, 591, 594, 596, 
                        597, 598, 599, 600, 601, 602, 603, 604, 621, 623, 624, 625, 626, 
                        627, 628, 629, 630, 632, 633, 634, 635, 636, 637, 638, 639, 640, 
                        641, 642, 643, 645, 646, 647, 648, 649, 650, 652, 653, 654, 655, 
                        656, 657, 658, 659, 660, 666, 677, 684, 689, 698, 700, 704, 708, 
                        709, 711, 712, 714, 718, 719, 721, 725, 726, 728, 733, 735, 736, 
                        737, 738, 739, 740, 741, 742, 744, 745, 746, 748, 752, 753, 755,
                        756, 758, 762, 763, 765, 769, 770, 772, 777, 781, 782, 784, 785, 
                        787, 791, 792, 794, 798, 799, 801, 805, 807, 808, 809, 810, 811, 
                        812, 813, 814, 816, 817, 818, 819, 820, 821, 822, 824, 828, 829, 
                        831, 832, 833, 834, 837, 841, 842, 843, 844, 846, 850, 851, 853, 
                        860, 864, 865, 867, 872, 876, 877, 879, 883, 885, 886, 887, 888, 
                        889, 890, 891, 892, 894, 895, 896, 897, 898, 899, 900, 901, 903, 
                        904, 905, 906, 908, 909, 910, 912, 916, 917, 919, 920, 923, 927, 
                        928, 930, 934, 935, 937, 942, 946, 947, 949, 950, 952, 956, 957, 
                        959, 963, 964, 966, 971, 973, 974, 975, 976, 977, 978, 980, 981, 
                        982, 983, 984, 985, 987, 988, 989, 990, 991, 992, 993, 994, 996, 
                        997, 998, 999, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1018, 
                        1020, 1021, 1022, 1023, 1024, 1025, 1026, 1031, 1033, 1034, 1035, 
                        1036, 1037, 1038, 1039, 1044, 1046, 1047, 1048, 1049, 1050, 1051, 
                        1052, 1053, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1066, 1068, 
                        1069, 1070, 1071, 1072, 1073, 1074, 1075, 1077, 1078, 1079, 1080, 
                        1081, 1082, 1083, 1088, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 
                        1097, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1110, 1112, 1113, 
                        1114, 1115, 1116, 1117, 1118, 1123, 1125, 1126, 1127, 1128, 1129, 
                        1130, 1131, 1132, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1143, 
                        1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1154, 1155, 1156, 
                        1157, 1158, 1159, 1160, 1161, 1163, 1164, 1165, 1166, 1167, 1168, 
                        1169, 1170, 1172, 1173, 1174, 1175, 1176, 1177, 1178
                        ],
                18 : [1241, 1242, 1247, 1249, 1275, 1277, 1295, 1297, 1300, 1302] ,
                24 : [1304, 1308, 1309, 1311, 1316, 1320, 1321, 1323, 1368, 1372, 1373, 1375, 
                        1380, 1384, 1385, 1387, 1392, 1397, 1401, 1402, 1404, 1408, 1409, 
                        1411, 1416, 1420, 1421, 1422, 1423, 1425, 1430, 1435, 1440, 1444, 
                        1445, 1447, 1452, 1456, 1457, 1459, 1463, 1465, 1466, 1467, 1468, 
                        1469, 1470, 1471, 1476, 1478, 1479, 1480, 1481, 1482, 1483, 1484, 
                        1485, 1487, 1488, 1489, 1490, 1491, 1492, 1493, 1494, 1496, 1497, 
                        1498, 1499, 1500, 1501, 1502, 1504, 1508, 1510, 1512, 1513, 1515, 
                        1516, 1518, 1520, 1522, 1535, 1537, 1542, 1544, 1546, 1548, 1561, 
                        1563, 1564, 1566, 1572, 1574, 1585, 1587
                        ],
                32 : [104, 141, 163, 235, 310, 321, 353, 357, 372, 376, 398, 405, 512, 529, 
                        546, 552, 554, 555, 556, 562, 563, 565, 578, 584, 586, 587, 588, 592, 
                        593, 595, 605, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 618, 
                        619, 620, 622, 631, 644, 651, 710, 732, 734, 743, 754, 783, 806, 815, 
                        830, 835, 857, 858, 884, 893, 902, 907, 918, 921, 948, 970, 972, 979, 
                        986, 995, 1000, 1008, 1009, 1011, 1012, 1013, 1014, 1015, 1016, 1019, 
                        1027, 1028, 1029, 1030, 1032, 1040, 1041, 1042, 1043, 1045, 1054, 1062, 
                        1063, 1064, 1065, 1067, 1076, 1084, 1085, 1086, 1087, 1089, 1098, 1106, 
                        1107, 1108, 1109, 1111, 1119, 1120, 1121, 1122, 1124, 1133, 1144, 1153, 
                        1162, 1171, 1179, 1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 
                        1190, 1191, 1192, 1193, 1194, 1195, 1196, 1198, 1199, 1200, 1201, 1202, 
                        1203, 1204, 1205, 1206, 1207, 1208, 1209, 1210, 1211, 1212, 1213, 1215, 
                        1216, 1217, 1218, 1219, 1220, 1221, 1222, 1224, 1225, 1226, 1227, 1228, 
                        1229, 1230
                        ],
                36 : [1248, 1250, 1276, 1278, 1296, 1298, 1299, 1301, 1327, 1329, 1330, 1331, 
                        1334, 1336, 1337, 1338
                        ],
                48 : [1464, 1472, 1473, 1474, 1475, 1477, 1486, 1495, 1506, 1509, 1514, 1517, 
                        1521, 1530, 1532, 1533, 1534, 1536, 1538, 1540, 1541, 1543, 1547, 1556, 
                        1558, 1559, 1560, 1562, 1565, 1567, 1569, 1570, 1571, 1573, 1580, 1582, 
                        1583, 1584, 1586, 1591, 1593, 1594, 1596, 1597, 1598, 1601, 1603, 1604, 
                        1605, 1606, 1608, 1609, 1610, 1611, 1613, 1614, 1615
                        ],
                64 : [354, 363, 606, 617, 1010, 1017, 1141, 1142, 1180, 1197, 1214, 1223],
                72 : [1328, 1332, 1333, 1335],
                96 : [1505, 1507, 1524, 1526, 1527, 1529, 1531, 1539, 1550, 1552, 1553, 1555, 
                        1557, 1568, 1577, 1579, 1581, 1588, 1590, 1592, 1595, 1602, 1607, 1612, 
                        1638, 1640, 1641, 1642, 1643, 1644, 1645, 1646, 1647, 1649, 1650, 1651
                        ],
                192 : [1519, 1523, 1525, 1528, 1545, 1549, 1551, 1554, 1575, 1576, 1578, 1589, 
                        1618, 1620, 1621, 1622, 1623, 1625, 1626, 1627, 1628, 1630, 1631, 1632, 
                        1633, 1635, 1636, 1637, 1639, 1648
                        ],
                384 : [1599, 1600, 1616, 1617, 1619, 1624, 1629, 1634],
                }
        
        all_mspg = [MagneticSpaceGroup(i) for i in range(*cs_range[cs]) if i in ORDER2MSPG[len(mop_list)]]
        # print("init mspg with dict: ", time.time() - start1)

        # slower alternative that doesn't need above dictionary
        # start1 = time.time()
        # all_mspg = [MagneticSpaceGroup(i) for i in range(*cs_range[cs])]
        # all_mspg = [mspg for mspg in all_mspg if len(mspg.symmetry_ops)==len(mop_list)]
        # print("init mspg: ", time.time() - start1)

        # compare list of symm operations to each magnetic spacegroup
        for tmp_mspg in all_mspg:  
            tmp_mspg_minus_input = []
            input_minus_tmp_mspg = []
        
            for mop1 in tmp_mspg.symmetry_ops:
                if not np.any([mop1 == mop2 for mop2 in mop_list]):
                    tmp_mspg_minus_input.append(mop1)

            if tmp_mspg_minus_input == []:
                for mop2 in mop_list:
                    if not np.any([mop1 == mop2 for mop1 in tmp_mspg.symmetry_ops]):
                        input_minus_tmp_mspg.append(mop2)

            if tmp_mspg_minus_input == [] and input_minus_tmp_mspg == []:
                return tmp_mspg
        
        return exit("Error in the list of symmetry operations.")

    def origin_centering(self):
      
        international_symbol = self._spacegroupAnalyzer._space_group_data['international']
        origins = []

        if international_symbol[0] == "F":
            origings = ["x, y+1/2, z+1/2, +1", "x+1/2, y, z+1/2, +1", "x+1/2, y+1/2, z, +1"]
        elif international_symbol[0] == "I":
            origings = ["x+1/2, y+1/2, z+1/2, +1"]
        elif international_symbol[0] == "R":
            origings = ["x+1/3, y+2/3, z+2/3, +1", "x+2/3, y+1/3, z+1/3, +1"]
        elif international_symbol[0] == "C":
            origings = ["x+1/2, y+1/2, z, +1"]
        elif international_symbol[0] == "A":
            origings = ["x, y+1/2, z+1/2, +1"]
        elif international_symbol[0] == "P":
            pass
        else:
            exit("Error: International symbol wrong. Argument given ", international_symbol)

        return [MagSymmOp.from_xyzt_string(xyzt_string) for xyzt_string in origings]

    # TODO: add a function get_refined_structure(self) analogous to SpacegroupAnalyzer
    # TODO: add function to return Magnetic Laue group
    # TODO: Use Magnetic Laue Group + Seemann et al (2015) to determine symmetry-imposed 
    #       shape of linear response tensor.
