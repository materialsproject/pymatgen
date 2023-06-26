"""This module provides classes to store, generate, and manipulate material interfaces."""

from __future__ import annotations

from itertools import product
from typing import TYPE_CHECKING, Iterator, Sequence

import numpy as np
from scipy.linalg import polar

from pymatgen.analysis.elasticity.strain import Deformation
from pymatgen.analysis.interfaces.zsl import ZSLGenerator, fast_norm
from pymatgen.core.interface import Interface, label_termination
from pymatgen.core.surface import SlabGenerator

if TYPE_CHECKING:
    from pymatgen.core import Structure


class CoherentInterfaceBuilder:
    """
    This class constructs the coherent interfaces between two crystalline slabs
    Coherency is defined by matching lattices not sub-planes.
    """

    def __init__(
        self,
        substrate_structure: Structure,
        film_structure: Structure,
        film_miller: tuple[int, int, int],
        substrate_miller: tuple[int, int, int],
        zslgen: ZSLGenerator | None = None,
    ):
        """
        Args:
            substrate_structure: structure of substrate
            film_structure: structure of film
            film_miller: miller index of the film layer
            substrate_miller: miller index for the substrate layer
            zslgen: BiDirectionalZSL if you want custom lattice matching tolerances for coherency.
        """
        # Bulk structures
        self.substrate_structure = substrate_structure
        self.film_structure = film_structure
        self.film_miller = film_miller
        self.substrate_miller = substrate_miller
        self.zslgen = zslgen or ZSLGenerator(bidirectional=True)

        self._find_matches()
        self._find_terminations()

    def _find_matches(self) -> None:
        """Finds and stores the ZSL matches."""
        self.zsl_matches = []

        film_sg = SlabGenerator(
            self.film_structure,
            self.film_miller,
            min_slab_size=1,
            min_vacuum_size=3,
            in_unit_planes=True,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        sub_sg = SlabGenerator(
            self.substrate_structure,
            self.substrate_miller,
            min_slab_size=1,
            min_vacuum_size=3,
            in_unit_planes=True,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        film_slab = film_sg.get_slab(shift=0)
        sub_slab = sub_sg.get_slab(shift=0)

        film_vectors = film_slab.lattice.matrix
        substrate_vectors = sub_slab.lattice.matrix

        # Generate all possible interface matches
        self.zsl_matches = list(self.zslgen(film_vectors[:2], substrate_vectors[:2], lowest=False))

        for match in self.zsl_matches:
            xform = get_2d_transform(film_vectors, match.film_vectors)
            strain, rot = polar(xform)
            assert np.allclose(
                strain, np.round(strain)
            ), "Film lattice vectors changed during ZSL match, check your ZSL Generator parameters"

            xform = get_2d_transform(substrate_vectors, match.substrate_vectors)
            strain, rot = polar(xform)
            assert np.allclose(
                strain, strain.astype(int)
            ), "Substrate lattice vectors changed during ZSL match, check your ZSL Generator parameters"

    def _find_terminations(self):
        """Finds all terminations."""
        film_sg = SlabGenerator(
            self.film_structure,
            self.film_miller,
            min_slab_size=1,
            min_vacuum_size=3,
            in_unit_planes=True,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        sub_sg = SlabGenerator(
            self.substrate_structure,
            self.substrate_miller,
            min_slab_size=1,
            min_vacuum_size=3,
            in_unit_planes=True,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        film_slabs = film_sg.get_slabs()
        sub_slabs = sub_sg.get_slabs()

        film_shits = [s.shift for s in film_slabs]
        film_terminations = [label_termination(s) for s in film_slabs]

        sub_shifts = [s.shift for s in sub_slabs]
        sub_terminations = [label_termination(s) for s in sub_slabs]

        self._terminations = {
            (film_label, sub_label): (film_shift, sub_shift)
            for (film_label, film_shift), (sub_label, sub_shift) in product(
                zip(film_terminations, film_shits), zip(sub_terminations, sub_shifts)
            )
        }
        self.terminations = list(self._terminations)

    def get_interfaces(
        self,
        termination: tuple[str, str],
        gap: float = 2.0,
        vacuum_over_film: float = 20.0,
        film_thickness: float | int = 1,
        substrate_thickness: float | int = 1,
        in_layers: bool = True,
    ) -> Iterator[Interface]:
        """Generates interface structures given the film and substrate structure
        as well as the desired terminations.

        Args:
            termination (tuple[str, str]): termination from self.termination list
            gap (float, optional): gap between film and substrate. Defaults to 2.0.
            vacuum_over_film (float, optional): vacuum over the top of the film. Defaults to 20.0.
            film_thickness (float | int, optional): the film thickness. Defaults to 1.
            substrate_thickness (float | int, optional): substrate thickness. Defaults to 1.
            in_layers (bool, optional): set the thickness in layer units. Defaults to True.

        Yields:
            Iterator[Interface]: interfaces from slabs
        """
        film_sg = SlabGenerator(
            self.film_structure,
            self.film_miller,
            min_slab_size=film_thickness,
            min_vacuum_size=3,
            in_unit_planes=in_layers,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        sub_sg = SlabGenerator(
            self.substrate_structure,
            self.substrate_miller,
            min_slab_size=substrate_thickness,
            min_vacuum_size=3,
            in_unit_planes=in_layers,
            center_slab=True,
            primitive=True,
            reorient_lattice=False,  # This is necessary to not screw up the lattice
        )

        film_shift, sub_shift = self._terminations[termination]

        film_slab = film_sg.get_slab(shift=film_shift)
        sub_slab = sub_sg.get_slab(shift=sub_shift)

        for match in self.zsl_matches:
            # Build film superlattice
            super_film_transform = np.round(
                from_2d_to_3d(get_2d_transform(film_slab.lattice.matrix[:2], match.film_sl_vectors))
            ).astype(int)
            film_sl_slab = film_slab.copy()
            film_sl_slab.make_supercell(super_film_transform)
            assert np.allclose(
                film_sl_slab.lattice.matrix[2], film_slab.lattice.matrix[2]
            ), "2D transformation affected C-axis for Film transformation"
            assert np.allclose(
                film_sl_slab.lattice.matrix[:2], match.film_sl_vectors
            ), "Transformation didn't make proper supercell for film"

            # Build substrate superlattice
            super_sub_transform = np.round(
                from_2d_to_3d(get_2d_transform(sub_slab.lattice.matrix[:2], match.substrate_sl_vectors))
            ).astype(int)
            sub_sl_slab = sub_slab.copy()
            sub_sl_slab.make_supercell(super_sub_transform)
            assert np.allclose(
                sub_sl_slab.lattice.matrix[2], sub_slab.lattice.matrix[2]
            ), "2D transformation affected C-axis for Film transformation"
            assert np.allclose(
                sub_sl_slab.lattice.matrix[:2], match.substrate_sl_vectors
            ), "Transformation didn't make proper supercell for substrate"

            # Add extra info
            match_dict = match.as_dict()
            interface_properties = {k: match_dict[k] for k in match_dict if not k.startswith("@")}

            dfm = Deformation(match.match_transformation)

            strain = dfm.green_lagrange_strain
            interface_properties["strain"] = strain
            interface_properties["von_mises_strain"] = strain.von_mises_strain
            interface_properties["termination"] = termination
            interface_properties["film_thickness"] = film_thickness
            interface_properties["substrate_thickness"] = substrate_thickness

            yield (
                Interface.from_slabs(
                    substrate_slab=sub_sl_slab,
                    film_slab=film_sl_slab,
                    gap=gap,
                    vacuum_over_film=vacuum_over_film,
                    interface_properties=interface_properties,
                )
            )


def get_rot_3d_for_2d(film_matrix, sub_matrix) -> np.ndarray:
    """Find transformation matrix that will rotate and strain the film to the substrate while preserving the c-axis."""
    film_matrix = np.array(film_matrix)
    film_matrix = film_matrix.tolist()[:2]
    film_matrix.append(np.cross(film_matrix[0], film_matrix[1]))

    # Generate 3D lattice vectors for substrate super lattice
    # Out of plane substrate super lattice has to be same length as
    # Film out of plane vector to ensure no extra deformation in that
    # direction
    sub_matrix = np.array(sub_matrix)
    sub_matrix = sub_matrix.tolist()[:2]
    temp_sub = np.cross(sub_matrix[0], sub_matrix[1]).astype(float)  # conversion to float necessary if using numba
    temp_sub = temp_sub * fast_norm(
        np.array(film_matrix[2], dtype=float)
    )  # conversion to float necessary if using numba
    sub_matrix.append(temp_sub)

    transform_matrix = np.transpose(np.linalg.solve(film_matrix, sub_matrix))

    rot, _ = polar(transform_matrix)

    return rot


def get_2d_transform(start: Sequence, end: Sequence) -> np.ndarray:
    """
    Gets a 2d transformation matrix
    that converts start to end.
    """
    return np.dot(end, np.linalg.pinv(start))


def from_2d_to_3d(mat: np.ndarray) -> np.ndarray:
    """Converts a 2D matrix to a 3D matrix."""
    new_mat = np.diag([1.0, 1.0, 1.0])
    new_mat[:2, :2] = mat
    return new_mat
