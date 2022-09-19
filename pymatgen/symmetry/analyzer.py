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

from __future__ import annotations

import copy
import itertools
import logging
import math
import warnings
from collections import defaultdict
from fractions import Fraction
from math import cos, sin
from typing import Literal

import numpy as np
import spglib

from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.core.periodic_table import Element, Species
from pymatgen.core.structure import Molecule, PeriodicSite, Structure
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.coord import find_in_coord_list, pbc_diff

logger = logging.getLogger(__name__)


class SpacegroupAnalyzer:
    """
    Takes a pymatgen.core.structure.Structure object and a symprec.
    Uses spglib to perform various symmetry finding operations.
    """

    def __init__(self, structure: Structure, symprec: float = 0.01, angle_tolerance=5.0):
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
        unique_species: list[Element | Species] = []
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

    def get_ir_reciprocal_mesh_map(self, mesh=(10, 10, 10), is_shift=(0, 0, 0)):
        """
        Same as 'get_ir_reciprocal_mesh' but the full grid together with
        the mapping that maps a reducible to an irreducible kpoint is
        returned.

        Args:
            mesh (3x1 array): The number of kpoint for the mesh needed in
                each direction
            is_shift (3x1 array): Whether to shift the kpoint grid. (1, 1,
            1) means all points are shifted by 0.5, 0.5, 0.5.

        Returns:
            A tuple containing two numpy.ndarray. The first is the mesh in
            fractional coordinates and the second is an array of integers
            that maps all the reducible kpoints from to irreducible ones.
        """
        shift = np.array([1 if i else 0 for i in is_shift])
        mapping, grid = spglib.get_ir_reciprocal_mesh(np.array(mesh), self._cell, is_shift=shift, symprec=self._symprec)

        grid_fractional_coords = (grid + shift * (0.5, 0.5, 0.5)) / mesh

        return grid_fractional_coords, mapping

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
            ({"vec": latt.matrix[i], "length": latt.abc[i], "orig_index": i} for i in range(3)),
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
                                for v, _ in self.rot_sym:
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
    for idx in range(len(mol)):
        clustered_dists[f[idx]].append(dists[idx])
    avg_dist = {label: np.mean(val) for label, val in clustered_dists.items()}
    clustered_sites = defaultdict(list)
    origin_site = None
    for idx, site in enumerate(mol):
        if avg_dist[f[idx]] < tol:
            if give_only_index:
                origin_site = idx
            else:
                origin_site = site
        else:
            if give_only_index:
                clustered_sites[(avg_dist[f[idx]], site.species)].append(idx)
            else:
                clustered_sites[(avg_dist[f[idx]], site.species)].append(site)
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

    def __init__(self, sch_symbol, operations, tol: float = 0.1):
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
        return str(self)
