"""An interface to the excellent spglib library by Atsushi Togo
(https://github.com/spglib/spglib) for pymatgen.

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
from collections.abc import Sequence
from fractions import Fraction
from functools import lru_cache
from math import cos, sin
from typing import TYPE_CHECKING

import numpy as np
import scipy.cluster
import spglib

from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Molecule, PeriodicSite, Structure
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.coord import find_in_coord_list, pbc_diff
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing import Any, Literal

    from numpy.typing import NDArray
    from spglib import SpglibDataset

    from pymatgen.core import Element, Species
    from pymatgen.core.sites import Site
    from pymatgen.symmetry.groups import CrystalSystem
    from pymatgen.util.typing import Kpoint

    LatticeType = Literal[
        "cubic",
        "hexagonal",
        "monoclinic",
        "orthorhombic",
        "rhombohedral",
        "tetragonal",
        "triclinic",
    ]

logger = logging.getLogger(__name__)


cite_conventional_cell_algo = due.dcite(
    Doi("10.1016/j.commatsci.2010.05.010"),
    description="High-throughput electronic band structure calculations: Challenges and tools",
)


class SymmetryUndeterminedError(ValueError):
    """
    An Exception for when symmetry cannot be determined. This might happen
    when, for example, atoms are very close together.
    """


@lru_cache(maxsize=32)
def _get_symmetry_dataset(cell, symprec, angle_tolerance):
    """Simple wrapper to cache results of spglib.get_symmetry_dataset since this call is
    expensive.
    """
    dataset = spglib.get_symmetry_dataset(cell, symprec=symprec, angle_tolerance=angle_tolerance)
    if dataset is None:
        raise SymmetryUndeterminedError
    return dataset


class SpacegroupAnalyzer:
    """Takes a pymatgen Structure object and a symprec.

    Uses spglib to perform various symmetry finding operations.
    """

    def __init__(
        self,
        structure: Structure,
        symprec: float | None = 0.01,
        angle_tolerance: float = 5,
    ) -> None:
        """
        Args:
            structure (Structure | IStructure): Structure to find symmetry
            symprec (float): Tolerance for symmetry finding. Defaults to 0.01,
                which is fairly strict and works well for properly refined
                structures with atoms in the proper symmetry coordinates. For
                structures with slight deviations from their proper atomic
                positions (e.g., structures relaxed with electronic structure
                codes), a looser tolerance of 0.1 (the value used in Materials
                Project) is often needed.
            angle_tolerance (float): Angle tolerance for symmetry finding. Defaults to 5 degrees.
        """
        self._symprec = symprec
        self._angle_tol = angle_tolerance
        self._structure = structure
        self._site_props = structure.site_properties
        unique_species: list[Element | Species] = []
        zs = []
        for species, group in itertools.groupby(structure, key=lambda s: s.species):
            if species in unique_species:
                ind = unique_species.index(species)
                zs.extend([ind + 1] * len(tuple(group)))
            else:
                unique_species.append(species)
                zs.extend([len(unique_species)] * len(tuple(group)))

        has_explicit_magmoms = "magmom" in structure.site_properties or any(
            getattr(specie, "spin", None) is not None for specie in structure.types_of_species
        )

        magmoms = []
        for site in structure:
            if hasattr(site, "magmom"):
                magmoms.append(site.magmom)
            elif site.is_ordered and getattr(site.specie, "spin", None) is not None:
                magmoms.append(site.specie.spin)
            elif has_explicit_magmoms:  # if any site has a magmom, all sites must have magmoms
                magmoms.append(0)

        self._unique_species = unique_species
        self._numbers = zs

        if len(magmoms) > 0:
            self._cell: tuple[Any, ...] = (
                tuple(map(tuple, structure.lattice.matrix.tolist())),
                tuple(map(tuple, structure.frac_coords.tolist())),
                tuple(zs),
                tuple(map(tuple, magmoms) if isinstance(magmoms[0], Sequence) else magmoms),
            )
        else:  # if no magmoms given do not add to cell
            self._cell = (
                tuple(map(tuple, structure.lattice.matrix.tolist())),
                tuple(map(tuple, structure.frac_coords.tolist())),
                tuple(zs),
            )

        self._space_group_data = _get_symmetry_dataset(self._cell, symprec, angle_tolerance)

    def get_space_group_symbol(self) -> str:
        """Get the spacegroup symbol (e.g., Pnma) for structure.

        Returns:
            str: Spacegroup symbol for structure.
        """
        return self._space_group_data.international

    def get_space_group_number(self) -> int:
        """Get the international spacegroup number (e.g., 62) for structure.

        Returns:
            int: International spacegroup number for structure.
        """
        return int(self._space_group_data.number)

    def get_space_group_operations(self) -> SpacegroupOperations:
        """Get the SpacegroupOperations for the Structure.

        Returns:
            SpacegroupOperations object.
        """
        return SpacegroupOperations(
            self.get_space_group_symbol(),
            self.get_space_group_number(),
            self.get_symmetry_operations(),
        )

    def get_hall(self) -> str:
        """Get Hall symbol for structure.

        Returns:
            str: Hall symbol
        """
        return self._space_group_data.hall

    def get_point_group_symbol(self) -> str:
        """Get the point group associated with the structure.

        Returns:
            Pointgroup: Point group for structure.
        """
        rotations = self._space_group_data.rotations
        # passing a 0-length rotations list to spglib can segfault
        if len(rotations) == 0:
            return "1"
        return spglib.get_pointgroup(rotations)[0].strip()

    def get_crystal_system(self) -> CrystalSystem:
        """Get the crystal system for the structure, e.g. (triclinic, orthorhombic,
        cubic, etc.).

        Raises:
            ValueError: on invalid space group numbers < 1 or > 230.

        Returns:
            str: Crystal system for structure
        """
        n = self._space_group_data.number

        # Not using isinstance(n, int) to allow 0-decimal floats
        if n != int(n) or not 0 < n < 231:
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

    def get_lattice_type(self) -> LatticeType:
        """Get the lattice for the structure, e.g. (triclinic, orthorhombic, cubic,
        etc.).This is the same as the crystal system with the exception of the
        hexagonal/rhombohedral lattice.

        Raises:
            ValueError: on invalid space group numbers < 1 or > 230.

        Returns:
            str: Lattice type for structure
        """
        spg_num = self._space_group_data.number
        system = self.get_crystal_system()
        if spg_num in {146, 148, 155, 160, 161, 166, 167}:
            return "rhombohedral"
        return "hexagonal" if system == "trigonal" else system

    def get_symmetry_dataset(self) -> SpglibDataset:
        """Get the symmetry dataset as a SpglibDataset.

        Returns:
            frozen dict: With the following properties:
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

    def _get_symmetry(self) -> tuple[NDArray, NDArray]:
        """Get the symmetry operations associated with the structure.

        Returns:
            Symmetry operations as a tuple of two equal length sequences.
            (rotations, translations). "rotations" is the numpy integer array
            of the rotation matrices for scaled positions
            "translations" gives the numpy float64 array of the translation
            vectors in scaled positions.
        """
        dct = spglib.get_symmetry(self._cell, symprec=self._symprec, angle_tolerance=self._angle_tol)
        if dct is None:
            symprec = self._symprec
            raise ValueError(
                f"Symmetry detection failed for structure with formula {self._structure.formula}. "
                f"Try setting {symprec=} to a different value."
            )
        # Sometimes spglib returns small translation vectors, e.g.
        # [1e-4, 2e-4, 1e-4]
        # (these are in fractional coordinates, so should be small denominator
        # fractions)
        translations: NDArray = np.array(
            [[float(Fraction(c).limit_denominator(1000)) for c in trans] for trans in dct["translations"]]
        )

        # Fractional translations of 1 are more simply 0
        translations[np.abs(translations) == 1] = 0
        return dct["rotations"], translations

    def get_symmetry_operations(self, cartesian: bool = False) -> list[SymmOp]:
        """Return symmetry operations as a list of SymmOp objects. By default returns
        fractional coord sym_ops. But Cartesian can be returned too.

        Returns:
            list[SymmOp]: symmetry operations.
        """
        rotation, translation = self._get_symmetry()
        sym_ops = []
        mat = self._structure.lattice.matrix.T
        inv_mat = np.linalg.inv(mat)
        for rot, trans in zip(rotation, translation, strict=True):
            if cartesian:
                rot = np.dot(mat, np.dot(rot, inv_mat))
                trans = np.dot(trans, self._structure.lattice.matrix)
            op = SymmOp.from_rotation_and_translation(rot, trans)
            sym_ops.append(op)
        return sym_ops

    def get_point_group_operations(self, cartesian: bool = False) -> list[SymmOp]:
        """Return symmetry operations as a list of SymmOp objects. By default returns
        fractional coord symm ops. But Cartesian can be returned too.

        Args:
            cartesian (bool): Whether to return SymmOps as Cartesian or
                direct coordinate operations.

        Returns:
            list[SymmOp]: Point group symmetry operations.
        """
        rotation, _translation = self._get_symmetry()
        symm_ops = []
        seen = set()
        mat = self._structure.lattice.matrix.T
        inv_mat = self._structure.lattice.inv_matrix.T
        for rot in rotation:
            rot_hash = rot.tobytes()
            if rot_hash in seen:
                continue
            seen.add(rot_hash)
            if cartesian:
                rot = np.dot(mat, np.dot(rot, inv_mat))
            op = SymmOp.from_rotation_and_translation(rot, np.array([0, 0, 0]))
            symm_ops.append(op)
        return symm_ops

    def get_symmetrized_structure(self) -> SymmetrizedStructure:
        """Get a symmetrized structure. A symmetrized structure is one where the sites
        have been grouped into symmetrically equivalent groups.

        Returns:
            pymatgen.symmetry.structure.SymmetrizedStructure object.
        """
        sym_dataset = self.get_symmetry_dataset()
        spg_ops = SpacegroupOperations(
            self.get_space_group_symbol(),
            self.get_space_group_number(),
            self.get_symmetry_operations(),
        )
        return SymmetrizedStructure(self._structure, spg_ops, sym_dataset.equivalent_atoms, sym_dataset.wyckoffs)

    def get_refined_structure(self, keep_site_properties: bool = False) -> Structure:
        """Get the refined structure based on detected symmetry. The refined structure is
        a *conventional* cell setting with atoms moved to the expected symmetry positions.

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
            for k, v in self._site_props.items():
                site_properties[k] = [v[self._numbers.index(i)] for i in numbers]
        else:
            site_properties = None
        struct = Structure(lattice, species, scaled_positions, site_properties=site_properties)
        return struct.get_sorted_structure()

    def find_primitive(self, keep_site_properties: bool = False) -> Structure:
        """Find a primitive version of the unit cell.

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
            as a Structure object. If no primitive cell is found, None is
            returned.
        """
        lattice, scaled_positions, numbers = spglib.find_primitive(self._cell, symprec=self._symprec)
        species = [self._unique_species[i - 1] for i in numbers]
        if keep_site_properties:
            site_properties = {}
            for key, val in self._site_props.items():
                site_properties[key] = [val[self._numbers.index(i)] for i in numbers]
        else:
            site_properties = None

        return Structure(
            lattice,
            species,
            scaled_positions,
            to_unit_cell=True,
            site_properties=site_properties,
        ).get_reduced_structure()

    def get_ir_reciprocal_mesh(
        self,
        mesh: tuple[int, int, int] = (10, 10, 10),
        is_shift: tuple[float, float, float] = (0, 0, 0),
    ) -> list[tuple[Kpoint, float]]:
        """k-point mesh of the Brillouin zone generated taken into account symmetry.
        The method returns the irreducible kpoints of the mesh and their weights.

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
        for idx, count in zip(*np.unique(mapping, return_counts=True), strict=True):
            results.append(((grid[idx] + shift * (0.5, 0.5, 0.5)) / mesh, count))
        return results

    def get_ir_reciprocal_mesh_map(
        self,
        mesh: tuple[int, int, int] = (10, 10, 10),
        is_shift: tuple[float, float, float] = (0, 0, 0),
    ) -> tuple[NDArray, NDArray]:
        """Same as 'get_ir_reciprocal_mesh' but the full grid together with the mapping
        that maps a reducible to an irreducible kpoint is returned.

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

    @cite_conventional_cell_algo
    def get_conventional_to_primitive_transformation_matrix(
        self,
        international_monoclinic: bool = True,
    ) -> NDArray:
        """Get the transformation matrix to transform a conventional unit cell to a
        primitive cell according to certain standards the standards are defined in
        Setyawan, W., & Curtarolo, S. (2010). High-throughput electronic band structure
        calculations: Challenges and tools. Computational Materials Science, 49(2),
        299-312. doi:10.1016/j.commatsci.2010.05.010.

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
            # Check if the conventional representation is hexagonal or
            # rhombohedral
            lengths = conv.lattice.lengths
            if abs(lengths[0] - lengths[2]) < 0.0001:
                return np.eye
            return np.array([[-1, 1, 1], [2, 1, 1], [-1, -2, 1]], dtype=np.float64) / 3

        if "I" in self.get_space_group_symbol():
            return np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype=np.float64) / 2

        if "F" in self.get_space_group_symbol():
            return np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=np.float64) / 2

        if "C" in self.get_space_group_symbol() or "A" in self.get_space_group_symbol():
            if self.get_crystal_system() == "monoclinic":
                return np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 2]], dtype=np.float64) / 2
            return np.array([[1, -1, 0], [1, 1, 0], [0, 0, 2]], dtype=np.float64) / 2

        return np.eye(3)

    @cite_conventional_cell_algo
    def get_primitive_standard_structure(
        self,
        international_monoclinic: bool = True,
        keep_site_properties: bool = False,
    ) -> Structure:
        """Get a structure with a primitive cell according to certain standards. The
        standards are defined in Setyawan, W., & Curtarolo, S. (2010). High-throughput
        electronic band structure calculations: Challenges and tools. Computational
        Materials Science, 49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010.

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
            international_monoclinic=international_monoclinic,
            keep_site_properties=keep_site_properties,
        )
        lattice = self.get_lattice_type()

        if "P" in self.get_space_group_symbol() or lattice == "hexagonal":
            return conv

        transf = self.get_conventional_to_primitive_transformation_matrix(
            international_monoclinic=international_monoclinic
        )

        new_sites: list[PeriodicSite] = []
        lattice = Lattice(np.dot(transf, conv.lattice.matrix))
        for site in conv:
            new_s = PeriodicSite(
                site.specie,
                site.coords,
                lattice,
                to_unit_cell=True,
                coords_are_cartesian=True,
                properties=site.properties,
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
            lattice = Lattice(new_matrix)
            for site in prim:
                new_s = PeriodicSite(
                    site.specie,
                    site.frac_coords,
                    lattice,
                    to_unit_cell=True,
                    properties=site.properties,
                )
                if not any(map(new_s.is_periodic_image, new_sites)):
                    new_sites.append(new_s)
            return Structure.from_sites(new_sites)

        return Structure.from_sites(new_sites)

    @cite_conventional_cell_algo
    def get_conventional_standard_structure(
        self,
        international_monoclinic: bool = True,
        keep_site_properties: bool = False,
    ) -> Structure:
        """Get a structure with a conventional cell according to certain standards. The
        standards are defined in Setyawan, W., & Curtarolo, S. (2010). High-throughput
        electronic band structure calculations: Challenges and tools. Computational
        Materials Science, 49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010 They
        basically enforce as much as possible norm(a1)<norm(a2)<norm(a3). NB This is not
        necessarily the same as the standard settings within the International Tables of
        Crystallography, for which get_refined_structure should be used instead.

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
        transf = None
        struct = self.get_refined_structure(keep_site_properties=keep_site_properties)
        lattice = struct.lattice
        latt_type = self.get_lattice_type()
        sorted_lengths = sorted(lattice.abc)
        sorted_dic = sorted(
            ({"vec": lattice.matrix[i], "length": lattice.abc[i], "orig_index": i} for i in range(3)),
            key=lambda k: k["length"],
        )

        if latt_type in {"orthorhombic", "cubic"}:
            # you want to keep the c axis where it is
            # to keep the C- settings
            transf = np.zeros(shape=(3, 3))
            if self.get_space_group_symbol().startswith("C"):
                transf[2] = [0, 0, 1]
                a, b = sorted(lattice.abc[:2])
                sorted_dic = sorted(
                    (
                        {
                            "vec": lattice.matrix[i],
                            "length": lattice.abc[i],
                            "orig_index": i,
                        }
                        for i in (0, 1)
                    ),
                    key=lambda k: k["length"],
                )
                for idx in range(2):
                    transf[idx][sorted_dic[idx]["orig_index"]] = 1
                c = lattice.abc[2]
            elif self.get_space_group_symbol().startswith(
                "A"
            ):  # change to C-centering to match Setyawan/Curtarolo convention
                transf[2] = [1, 0, 0]
                a, b = sorted(lattice.abc[1:])
                sorted_dic = sorted(
                    (
                        {
                            "vec": lattice.matrix[i],
                            "length": lattice.abc[i],
                            "orig_index": i,
                        }
                        for i in (1, 2)
                    ),
                    key=lambda k: k["length"],
                )
                for idx in range(2):
                    transf[idx][sorted_dic[idx]["orig_index"]] = 1
                c = lattice.abc[0]
            else:
                for idx, dct in enumerate(sorted_dic):
                    transf[idx][dct["orig_index"]] = 1
                a, b, c = sorted_lengths
            lattice = Lattice.orthorhombic(a, b, c)

        elif latt_type == "tetragonal":
            # find the "a" vectors
            # it is basically the vector repeated two times
            transf = np.zeros(shape=(3, 3))
            a, b, c = sorted_lengths
            for idx, dct in enumerate(sorted_dic):
                transf[idx][dct["orig_index"]] = 1

            if abs(b - c) < tol < abs(a - c):
                a, c = c, a
                transf = np.dot([[0, 0, 1], [0, 1, 0], [1, 0, 0]], transf)
            lattice = Lattice.tetragonal(a, c)

        elif latt_type in {"hexagonal", "rhombohedral"}:
            # for the conventional cell representation,
            # we always show the rhombohedral lattices as hexagonal

            # check first if we have the refined structure shows a rhombohedral
            # cell
            # if so, make a supercell
            a, b, c = lattice.abc
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
            lattice = Lattice(new_matrix)
            transf = np.eye(3, 3)

        elif latt_type == "monoclinic":
            # You want to keep the c axis where it is to keep the C- settings

            if self.get_space_group_operations().int_symbol.startswith("C"):
                transf = np.zeros(shape=(3, 3))
                transf[2] = [0, 0, 1]
                sorted_dic = sorted(
                    (
                        {
                            "vec": lattice.matrix[i],
                            "length": lattice.abc[i],
                            "orig_index": i,
                        }
                        for i in (0, 1)
                    ),
                    key=lambda k: k["length"],
                )
                a = sorted_dic[0]["length"]
                b = sorted_dic[1]["length"]
                c = lattice.abc[2]
                new_matrix = None
                for tp2 in itertools.permutations(list(range(2)), 2):
                    m = lattice.matrix
                    latt2 = Lattice([m[tp2[0]], m[tp2[1]], m[2]])
                    lengths = latt2.lengths
                    angles = latt2.angles
                    if angles[0] > 90:
                        # if the angle is > 90 we invert a and b to get
                        # an angle < 90
                        a, b, c, alpha, beta, gamma = Lattice([-m[tp2[0]], -m[tp2[1]], m[2]]).parameters
                        transf = np.zeros(shape=(3, 3))
                        transf[0][tp2[0]] = -1
                        transf[1][tp2[1]] = -1
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
                        transf[0][tp2[0]] = 1
                        transf[1][tp2[1]] = 1
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
                    for idx, dct in enumerate(sorted_dic):
                        transf[idx][dct["orig_index"]] = 1
            # if not C-setting
            else:
                # try all permutations of the axis
                # keep the ones with the non-90 angle=alpha
                # and b<c
                new_matrix = None

                for tp3 in itertools.permutations(list(range(3)), 3):
                    m = lattice.matrix
                    a, b, c, alpha, beta, gamma = Lattice([m[tp3[0]], m[tp3[1]], m[tp3[2]]]).parameters
                    if alpha > 90 and b < c:
                        a, b, c, alpha, beta, gamma = Lattice([-m[tp3[0]], -m[tp3[1]], m[tp3[2]]]).parameters
                        transf = np.zeros(shape=(3, 3))
                        transf[0][tp3[0]] = -1
                        transf[1][tp3[1]] = -1
                        transf[2][tp3[2]] = 1
                        alpha = math.pi * alpha / 180
                        new_matrix = [
                            [a, 0, 0],
                            [0, b, 0],
                            [0, c * cos(alpha), c * sin(alpha)],
                        ]
                        continue

                    if alpha < 90 and b < c:
                        transf = np.zeros(shape=(3, 3))
                        transf[0][tp3[0]] = 1
                        transf[1][tp3[1]] = 1
                        transf[2][tp3[2]] = 1
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
                    for idx, dct in enumerate(sorted_dic):
                        transf[idx][dct["orig_index"]] = 1

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

            lattice = Lattice(new_matrix)

        elif latt_type == "triclinic":
            # we use a LLL Minkowski-like reduction for the triclinic cells
            struct = struct.get_reduced_structure("LLL")
            lattice = struct.lattice

            a, b, c = lattice.lengths
            alpha, beta, gamma = (math.pi * i / 180 for i in lattice.angles)
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

            def is_all_acute_or_obtuse(matrix) -> bool:
                recp_angles = np.array(Lattice(matrix).reciprocal_lattice.angles)
                return all(recp_angles <= 90) or all(recp_angles > 90)

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

            lattice = Lattice(new_matrix)

        new_coords = np.dot(transf, np.transpose(struct.frac_coords)).T
        new_struct = Structure(
            lattice,
            struct.species_and_occu,
            new_coords,
            site_properties=struct.site_properties,
            to_unit_cell=True,
        )
        return new_struct.get_sorted_structure()

    def get_kpoint_weights(self, kpoints: Sequence[Kpoint], atol: float = 1e-5) -> list[float]:
        """Calculate the weights for a list of kpoints.

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
        for idx in range(3):
            nonzero = [i for i in kpts[:, idx] if abs(i) > 1e-5]
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
        mapped: dict[tuple, int] = defaultdict(int)
        for kpt in kpoints:
            for idx, g in enumerate(grid):
                if np.allclose(pbc_diff(kpt, g), (0, 0, 0), atol=atol):
                    mapped[tuple(g)] += 1
                    weights.append(mapping.count(mapping[idx]))
                    break
        if (len(mapped) != len(set(mapping))) or any(v != 1 for v in mapped.values()):
            raise ValueError("Unable to find 1:1 corresponding between input kpoints and irreducible grid!")
        return [w / sum(weights) for w in weights]

    def is_laue(self) -> bool:
        """Check if the point group of the structure has Laue symmetry (centrosymmetry)."""
        laue = {
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
        }

        return str(self.get_point_group_symbol()) in laue


class PointGroupAnalyzer:
    """A class to analyze the point group of a molecule.

    The general outline of the algorithm is as follows:

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

    Attribute:
        sch_symbol (str): Schoenflies symbol of the detected point group.
    """

    inversion_op = SymmOp.inversion()

    def __init__(
        self,
        mol: Molecule,
        tolerance: float = 0.3,
        eigen_tolerance: float = 0.01,
        matrix_tolerance: float = 0.1,
    ) -> None:
        """The default settings are usually sufficient.

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
        if self.sch_symbol in {"C1v", "C1h"}:
            self.sch_symbol: str = "Cs"

    def _analyze(self) -> None:
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
                for i, j in ((0, 1), (1, 2), (0, 2)):
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

            self.rot_sym: list = []
            self.symmops: list[SymmOp] = [SymmOp(np.eye(4))]
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

    def _proc_linear(self) -> None:
        if self.is_valid_op(PointGroupAnalyzer.inversion_op):
            self.sch_symbol = "D*h"
            self.symmops.append(PointGroupAnalyzer.inversion_op)
        else:
            self.sch_symbol = "C*v"

    def _proc_asym_top(self) -> None:
        """Handles asymmetric top molecules, which cannot contain rotational symmetry
        larger than 2.
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

    def _proc_sym_top(self) -> None:
        """Handles symmetric top molecules which has one unique eigenvalue whose
        corresponding principal axis is a unique rotational axis.

        More complex handling required to look for R2 axes perpendicular to this unique
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

    def _proc_no_rot_sym(self) -> None:
        """Handles molecules with no rotational symmetry.

        Only possible point groups are C1, Cs and Ci.
        """
        self.sch_symbol = "C1"
        if self.is_valid_op(PointGroupAnalyzer.inversion_op):
            self.sch_symbol = "Ci"
            self.symmops.append(PointGroupAnalyzer.inversion_op)
        else:
            for v in self.principal_axes:
                mirror_type = self._find_mirror(v)
                if mirror_type != "":
                    self.sch_symbol = "Cs"
                    break

    def _proc_cyclic(self) -> None:
        """Handles cyclic group molecules."""
        main_axis, rot = max(self.rot_sym, key=lambda v: v[1])
        self.sch_symbol = f"C{rot}"
        mirror_type = self._find_mirror(main_axis)
        if mirror_type == "h":
            self.sch_symbol += "h"
        elif mirror_type == "v":
            self.sch_symbol += "v"
        elif mirror_type == "" and self.is_valid_op(SymmOp.rotoreflection(main_axis, angle=180 / rot)):
            self.sch_symbol = f"S{2 * rot}"

    def _proc_dihedral(self) -> None:
        """Handles dihedral group molecules, i.e those with intersecting R2 axes and a
        main axis.
        """
        main_axis, rot = max(self.rot_sym, key=lambda v: v[1])
        self.sch_symbol = f"D{rot}"
        mirror_type = self._find_mirror(main_axis)
        if mirror_type == "h":
            self.sch_symbol += "h"
        elif mirror_type != "":
            self.sch_symbol += "d"

    def _check_R2_axes_asym(self) -> None:
        """Test for 2-fold rotation along the principal axes.

        Used to handle asymmetric top molecules.
        """
        for v in self.principal_axes:
            op = SymmOp.from_axis_angle_and_translation(v, 180)
            if self.is_valid_op(op):
                self.symmops.append(op)
                self.rot_sym.append((v, 2))

    def _find_mirror(self, axis: NDArray) -> Literal["h", "d", "v", ""]:
        """Looks for mirror symmetry of specified type about axis.

        Possible types are "h" or "vd". Horizontal (h) mirrors are perpendicular to the
        axis while vertical (v) or diagonal (d) mirrors are parallel. v mirrors has atoms
        lying on the mirror plane while d mirrors do not.
        """
        mirror_type: Literal["h", "d", "v", ""] = ""

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
                                    if np.linalg.norm(v - axis) >= self.tol and np.dot(v, normal) < self.tol:
                                        mirror_type = "v"
                                        break
                            else:
                                mirror_type = "v"
                            break

        return mirror_type

    def _get_smallest_set_not_on_axis(self, axis: NDArray) -> list:
        """Get the smallest list of atoms with the same species and distance from
        origin AND does not lie on the specified axis.

        This maximal set limits the possible rotational symmetry operations, since atoms
        lying on a test axis is irrelevant in testing rotational symmetryOperations.
        """

        def not_on_axis(site):
            return np.linalg.norm(np.cross(site.coords, axis)) > self.tol

        valid_sets = []
        _origin_site, dist_el_sites = cluster_sites(self.centered_mol, self.tol)
        for test_set in dist_el_sites.values():
            valid_set = list(filter(not_on_axis, test_set))
            if len(valid_set) > 0:
                valid_sets.append(valid_set)

        return min(valid_sets, key=len)

    def _check_rot_sym(self, axis: NDArray) -> int:
        """Determine the rotational symmetry about supplied axis.

        Used only for symmetric top molecules which has possible rotational symmetry
        operations > 2.
        """
        min_set = self._get_smallest_set_not_on_axis(axis)
        max_sym = len(min_set)
        for idx in range(max_sym, 0, -1):
            if max_sym % idx != 0:
                continue
            op = SymmOp.from_axis_angle_and_translation(axis, 360 / idx)
            if self.is_valid_op(op):
                self.symmops.append(op)
                self.rot_sym.append((axis, idx))
                return idx
        return 1

    def _check_perpendicular_r2_axis(self, axis: NDArray) -> None | Literal[True]:
        """Check for R2 axes perpendicular to unique axis.

        For handling symmetric top molecules.
        """
        min_set = self._get_smallest_set_not_on_axis(axis)
        for s1, s2 in itertools.combinations(min_set, 2):
            test_axis = np.cross(s1.coords - s2.coords, axis)
            if np.linalg.norm(test_axis) > self.tol:
                op = SymmOp.from_axis_angle_and_translation(test_axis, 180)
                if self.is_valid_op(op):
                    self.symmops.append(op)
                    self.rot_sym.append((test_axis, 2))
                    return True
        return None

    def _proc_sph_top(self) -> None:
        """Handles Spherical Top Molecules, which belongs to the T, O or I point
        groups.
        """
        self._find_spherical_axes()
        if len(self.rot_sym) == 0:
            logger.debug("Accidental spherical top!")
            self._proc_sym_top()
        main_axis, rot = max(self.rot_sym, key=lambda v: v[1])
        if rot < 3:
            logger.debug("Accidental spherical top!")
            self._proc_sym_top()

        elif rot == 3:
            mirror_type = self._find_mirror(main_axis)
            if mirror_type == "":
                self.sch_symbol = "T"
            elif self.is_valid_op(PointGroupAnalyzer.inversion_op):
                self.symmops.append(PointGroupAnalyzer.inversion_op)
                self.sch_symbol = "Th"
            else:
                self.sch_symbol = "Td"

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

    def _find_spherical_axes(self) -> None:
        """Looks for R5, R4, R3 and R2 axes in spherical top molecules.

        Point group T molecules have only one unique 3-fold and one unique 2-fold axis. O
        molecules have one unique 4, 3 and 2-fold axes. I molecules have a unique 5-fold
        axis.
        """
        rot_present: dict[int, bool] = defaultdict(bool)
        _origin_site, dist_el_sites = cluster_sites(self.centered_mol, self.tol)
        test_set = min(dist_el_sites.values(), key=len)
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

    def get_pointgroup(self) -> PointGroupOperations:
        """Get a PointGroup object for the molecule."""
        return PointGroupOperations(self.sch_symbol, self.symmops, self.mat_tol)

    def get_symmetry_operations(self) -> Sequence[SymmOp]:
        """Get symmetry operations.

        Returns:
            list[SymmOp]: symmetry operations in Cartesian coord.
        """
        return generate_full_symmops(self.symmops, self.tol)

    def get_rotational_symmetry_number(self) -> int:
        """Get rotational symmetry number.

        Returns:
            int: Rotational symmetry number.
        """
        if self.sch_symbol == "D*h":
            # Special case. H2 for example has rotational symmetry number 2
            return 2

        """Get the rotational symmetry number."""
        symm_ops = self.get_symmetry_operations()
        symm_number = 0
        for symm in symm_ops:
            rot = symm.rotation_matrix
            if np.abs(np.linalg.det(rot) - 1) < 1e-4:
                symm_number += 1
        return symm_number

    def is_valid_op(self, symm_op: SymmOp) -> bool:
        """Check if a particular symmetry operation is a valid symmetry operation for a
        molecule, i.e., the operation maps all atoms to another equivalent atom.

        Args:
            symm_op (SymmOp): Symmetry operation to test.

        Returns:
            bool: True if SymmOp is valid for Molecule.
        """
        coords = self.centered_mol.cart_coords
        for site in self.centered_mol:
            coord = symm_op.operate(site.coords)
            ind = find_in_coord_list(coords, coord, self.tol)
            if len(ind) != 1 or self.centered_mol[ind[0]].species != site.species:
                return False
        return True

    def _get_eq_sets(self) -> dict[Literal["eq_sets", "sym_ops"], Any]:
        """Calculate the dictionary for mapping equivalent atoms onto each other.

        Returns:
            dict: with two possible keys:
                eq_sets: A dictionary of indices mapping to sets of indices, each key maps to
                    indices of all equivalent atoms. The keys are guaranteed to be not equivalent.
                sym_ops: Twofold nested dictionary. operations[i][j] gives the symmetry
                    operation that maps atom i unto j.
        """
        unit_matrix = np.eye(3)
        eq_sets: dict[int, set] = defaultdict(set)
        operations: dict[int, dict] = defaultdict(dict)
        symm_ops = [op.rotation_matrix for op in generate_full_symmops(self.symmops, self.tol)]

        def get_clustered_indices():
            indices = cluster_sites(self.centered_mol, self.tol, give_only_index=True)
            out = list(indices[1].values())
            if indices[0] is not None:
                out.append([indices[0]])
            return out

        for index in get_clustered_indices():
            sites = self.centered_mol.cart_coords[index]
            for idx, reference in zip(index, sites, strict=True):
                for op in symm_ops:
                    rotated = np.dot(op, sites.T).T
                    matched_indices = find_in_coord_list(rotated, reference, self.tol)
                    matched_indices = {dict(enumerate(index))[i] for i in matched_indices}
                    eq_sets[idx] |= matched_indices

                    if idx not in operations:
                        operations[idx] = {j: op.T if j != idx else unit_matrix for j in matched_indices}
                    else:
                        for j in matched_indices:
                            if j not in operations[idx]:
                                operations[idx][j] = op.T if j != idx else unit_matrix
                    for j in matched_indices:
                        if j not in operations:
                            operations[j] = {idx: op if j != idx else unit_matrix}
                        elif idx not in operations[j]:
                            operations[j][idx] = op if j != idx else unit_matrix

        return {"eq_sets": eq_sets, "sym_ops": operations}

    @staticmethod
    def _combine_eq_sets(equiv_sets: dict, sym_ops: dict) -> dict:
        """Combines the dicts of _get_equivalent_atom_dicts into one.

        Args:
            equiv_sets (dict): Map of equivalent atoms onto each other (i.e. indices to indices).
            sym_ops (dict): Map of symmetry operations that map atoms onto each other.

        Returns:
            dict: with two possible keys:
                eq_sets: A dictionary of indices mapping to sets of indices, each key maps to
                    indices of all equivalent atoms. The keys are guaranteed to be not equivalent.
                sym_ops: Twofold nested dictionary. operations[i][j] gives the symmetry
                    operation that maps atom i unto j.
        """
        unit_matrix = np.eye(3)

        def all_equivalent_atoms_of_i(idx, eq_sets, ops):
            """WORKS INPLACE on operations."""
            visited = {idx}
            tmp_eq_sets = {j: (eq_sets[j] - visited) for j in eq_sets[idx]}

            while tmp_eq_sets:
                new_tmp_eq_sets = {}
                for j, eq_set in tmp_eq_sets.items():
                    if j in visited:
                        continue
                    visited.add(j)
                    for k in eq_set:
                        new_tmp_eq_sets[k] = eq_sets[k] - visited
                        if idx not in ops[k]:
                            ops[k][idx] = np.dot(ops[j][idx], ops[k][j]) if k != idx else unit_matrix
                        ops[idx][k] = ops[k][idx].T
                tmp_eq_sets = new_tmp_eq_sets
            return visited, ops

        equiv_sets = copy.deepcopy(equiv_sets)
        ops = copy.deepcopy(sym_ops)
        to_be_deleted = set()
        for idx in equiv_sets:
            if idx in to_be_deleted:
                continue
            visited, ops = all_equivalent_atoms_of_i(idx, equiv_sets, ops)
            to_be_deleted |= visited - {idx}

        for key in to_be_deleted:
            equiv_sets.pop(key, None)
        return {"eq_sets": equiv_sets, "sym_ops": ops}

    def get_equivalent_atoms(self):
        """Get sets of equivalent atoms with symmetry operations.

        Returns:
            dict: with two possible keys:
                eq_sets: A dictionary of indices mapping to sets of indices, each key maps to
                    indices of all equivalent atoms. The keys are guaranteed to be not equivalent.
                sym_ops: Twofold nested dictionary. operations[i][j] gives the symmetry
                    operation that maps atom i unto j.
        """
        eq = self._get_eq_sets()
        return self._combine_eq_sets(eq["eq_sets"], eq["sym_ops"])

    def symmetrize_molecule(self) -> dict:
        """Get a symmetrized molecule.

        The equivalent atoms obtained via
        :meth:`~pymatgen.symmetry.analyzer.PointGroupAnalyzer.get_equivalent_atoms`
        are rotated, mirrored... unto one position.
        Then the average position is calculated.
        The average position is rotated, mirrored... back with the inverse
        of the previous symmetry operations, which gives the
        symmetrized molecule

        Returns:
            dict: with three possible keys:
                sym_mol: A symmetrized molecule instance.
                eq_sets: A dictionary of indices mapping to sets of indices, each key maps to indices
                    of all equivalent atoms. The keys are guaranteed to be not equivalent.
                sym_ops: Twofold nested dictionary. operations[i][j] gives the symmetry operation
                    that maps atom i unto j.
        """
        eq = self.get_equivalent_atoms()
        eq_sets, ops = eq["eq_sets"], eq["sym_ops"]
        coords = self.centered_mol.cart_coords.copy()
        for idx, eq_indices in eq_sets.items():
            for j in eq_indices:
                coords[j] = np.dot(ops[j][idx], coords[j])
            coords[idx] = np.mean(coords[list(eq_indices)], axis=0)
            for j in eq_indices:
                if j == idx:
                    continue
                coords[j] = np.dot(ops[idx][j], coords[idx])
                coords[j] = np.dot(ops[idx][j], coords[idx])
        molecule = Molecule(species=self.centered_mol.species_and_occu, coords=coords)
        return {"sym_mol": molecule, "eq_sets": eq_sets, "sym_ops": ops}


def iterative_symmetrize(
    mol: Molecule,
    max_n: int = 10,
    tolerance: float = 0.3,
    epsilon: float = 1e-2,
) -> dict[Literal["sym_mol", "eq_sets", "sym_ops"], Molecule | dict]:
    """Get a symmetrized molecule.

    The equivalent atoms obtained via `PointGroupAnalyzer.get_equivalent_atoms`
    are rotated, mirrored... unto one position.
    Then the average position is calculated, which is rotated, mirrored...
    back with the inverse of the previous symmetry operations, giving the
    symmetrized molecule.

    Args:
        mol (Molecule): A pymatgen Molecule instance.
        max_n (int): Maximum number of iterations.
        tolerance (float): Tolerance for detecting symmetry with PointGroupAnalyzer.
        epsilon (float): If the element-wise absolute difference of two
            subsequently symmetrized structures is smaller epsilon,
            the iteration stops before max_n is reached.

    Returns:
        dict with three keys:
            sym_mol: A symmetrized Molecule instance.
            eq_sets: A dictionary of indices mapping to sets of indices, each key maps to indices
                of all equivalent atoms. The keys are guaranteed to be not equivalent.
            sym_ops: Two-fold nested dictionary. operations[i][j] gives the symmetry operation
                that maps atom i unto j.
    """
    new_mol: Molecule = mol
    sym_mol: dict = {"sym_mol": new_mol, "eq_sets": {}, "sym_ops": {}}
    for _ in range(max_n):
        prev_mol: Molecule = new_mol
        sym_mol = PointGroupAnalyzer(prev_mol, tolerance=tolerance).symmetrize_molecule()
        new_mol = sym_mol["sym_mol"]

        if np.allclose(new_mol.cart_coords, prev_mol.cart_coords, atol=epsilon):
            break
    return sym_mol


def cluster_sites(
    mol: Molecule,
    tol: float,
    give_only_index: bool = False,
) -> tuple[Site | None, dict]:
    """Cluster sites based on distance and species type.

    Args:
        mol (Molecule): Molecule **with origin at center of mass**.
        tol (float): Tolerance to use.
        give_only_index (bool): Whether to return only the index of the
            origin site, instead of the site itself. Defaults to False.

    Returns:
        tuple[Site | None, dict]: origin_site is a site at the center
            of mass (None if there are no origin atoms). clustered_sites is a
            dict of {(avg_dist, species_and_occu): [list of sites]}
    """
    # Cluster works for dim > 2 data. We just add a dummy 0 for second
    # coordinate.
    dists: list[list[float]] = [[float(np.linalg.norm(site.coords)), 0] for site in mol]

    f_cluster = scipy.cluster.hierarchy.fclusterdata(dists, tol, criterion="distance")
    clustered_dists: dict[str, list[list[float]]] = defaultdict(list)
    for idx in range(len(mol)):
        clustered_dists[f_cluster[idx]].append(dists[idx])
    avg_dist = {key: np.mean(val) for key, val in clustered_dists.items()}
    clustered_sites = defaultdict(list)
    origin_site = None
    for idx, site in enumerate(mol):
        if avg_dist[f_cluster[idx]] < tol:
            origin_site = idx if give_only_index else site
        elif give_only_index:
            clustered_sites[avg_dist[f_cluster[idx]], site.species].append(idx)
        else:
            clustered_sites[avg_dist[f_cluster[idx]], site.species].append(site)
    return origin_site, clustered_sites


def generate_full_symmops(
    symmops: Sequence[SymmOp],
    tol: float,
) -> Sequence[SymmOp]:
    """Recursive algorithm to permute through all possible combinations of the initially
    supplied symmetry operations to arrive at a complete set of operations mapping a
    single atom to all other equivalent atoms in the point group. This assumes that the
    initial number already uniquely identifies all operations.

    Args:
        symmops (list[SymmOp]): Initial set of symmetry operations.
        tol (float): Tolerance for detecting symmetry.

    Returns:
        list[SymmOp]: Full set of symmetry operations.
    """
    # Uses an algorithm described in:
    # Gregory Butler. Fundamental Algorithms for Permutation Groups.
    # Lecture Notes in Computer Science (Book 559). Springer, 1991. page 15
    identity = np.eye(4)
    generators = [op.affine_matrix for op in symmops if not np.allclose(op.affine_matrix, identity)]
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
                    " and rerun with a different tolerance."
                )

    d = np.abs(full - identity) < tol
    if not np.any(np.all(np.all(d, axis=2), axis=1)):
        full.append(identity)
    return [SymmOp(op) for op in full]


class SpacegroupOperations(list):
    """Represents a space group, which is a collection of symmetry operations."""

    def __init__(
        self,
        int_symbol: str,
        int_number: int,
        symmops: Sequence[SymmOp],
    ) -> None:
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

    def __str__(self) -> str:
        return f"{self.int_symbol} ({self.int_number}) spacegroup"

    def are_symmetrically_equivalent(
        self,
        sites1: set[PeriodicSite],
        sites2: set[PeriodicSite],
        symm_prec: float = 1e-3,
    ) -> bool:
        """Given two sets of PeriodicSites, test if they are actually symmetrically
        equivalent under this space group. Useful, for example, if you want to test if
        selecting atoms 1 and 2 out of a set of 4 atoms are symmetrically the same as
        selecting atoms 3 and 4, etc.

        One use is in PartialRemoveSpecie transformation to return only
        symmetrically distinct arrangements of atoms.

        Args:
            sites1 ([PeriodicSite]): 1st set of sites
            sites2 ([PeriodicSite]): 2nd set of sites
            symm_prec (float): Tolerance in atomic distance to test if atoms
                are symmetrically similar.

        Returns:
            bool: True if the two sets of sites are symmetrically equivalent.
        """

        def in_sites(site):
            return any(test_site.is_periodic_image(site, symm_prec, check_lattice=False) for test_site in sites1)

        for op in self:
            new_sites2 = [PeriodicSite(site.species, op.operate(site.frac_coords), site.lattice) for site in sites2]
            for site in new_sites2:
                if not in_sites(site):
                    break
            else:
                return True
        return False


class PointGroupOperations(list):
    """Represents a point group, which is a sequence of symmetry operations.

    Attributes:
        sch_symbol (str): Schoenflies symbol of the point group.
    """

    def __init__(
        self,
        sch_symbol: str,
        operations: Sequence[SymmOp],
        tol: float = 0.1,
    ) -> None:
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

    def __repr__(self) -> str:
        return self.sch_symbol
