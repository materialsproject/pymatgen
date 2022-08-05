"""
Provides analysis of site symmetries.
"""

from __future__ import annotations

import numpy as np

from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga


def get_site_symmetries(struct: Structure, precision: float = 0.1) -> list[list[SymmOp]]:
    """
    Get all the point group operations centered on each atomic site
    in the form [[point operations of site index 1]...[[point operations of site index N]]]

    Args:
        struct: Pymatgen structure
        precision (float): tolerance to find symmetry operations

    Return:
        list of lists of point operations for each atomic site
    """

    point_ops: list[list[SymmOp]] = []

    # Point symmetries of each atom
    for site1 in range(len(struct.sites)):
        temp_struct = struct.copy()

        # Place the origin of the cell at each atomic site
        point_ops.append([])

        for site2 in range(len(struct.sites)):  # pylint: disable=C0200
            temp_struct.replace(
                site2,
                temp_struct.sites[site2].specie,
                temp_struct.frac_coords[site2] - struct.frac_coords[site1],
            )

        sgastruc = sga(temp_struct, symprec=precision)
        ops = sgastruc.get_symmetry_operations(cartesian=True)
        for op in ops:
            if all(op.translation_vector == [0, 0, 0]):
                point_ops[site1].append(op)
    return point_ops


def get_shared_symmetry_operations(struct: Structure, pointops: list[list[SymmOp]], tol=0.1):
    """
    Get all the point group operations shared by a pair of atomic sites
    in the form [[point operations of site index 1],[],...,[]]

    Args:
        struct: Pymatgen structure
        pointops: list of point group operations from get_site_symmetries method

    Return:
        list of lists of shared point operations for each pair of atomic sites
    """
    num_sites = len(struct)
    shared_ops = np.zeros((num_sites, num_sites), dtype=object)
    for site1 in range(num_sites):
        for site2 in range(num_sites):
            shared_ops[site1][site2] = []
            for pop1 in pointops[site1]:
                for pop2 in pointops[site2]:
                    if np.allclose(pop1.rotation_matrix, pop2.rotation_matrix):
                        shared_ops[site1][site2].append(pop1)

    for site1, sops in enumerate(shared_ops):
        for site2, sop in enumerate(sops):
            unique_ops = []
            for ops in sop:
                op = SymmOp.from_rotation_and_translation(
                    rotation_matrix=ops.rotation_matrix,
                    translation_vec=(0, 0, 0),
                    tol=tol,
                )
                if op not in unique_ops:
                    unique_ops.append(op)

            shared_ops[site1][site2] = unique_ops

    return shared_ops
