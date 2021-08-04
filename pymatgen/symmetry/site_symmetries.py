"""
Provides analysis of site symmetries.
"""

import numpy as np

from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga


def get_site_symmetries(struc, precision=0.1):
    """
    Get all the point group operations centered on each atomic site
    in the form [[point operations of site index 1]...[[point operations of site index N]]]

    Args:
        struc: Pymatgen structure
        precision (float): tolerance to find symmetry operaitons

    Return:
        list of lists of point operations for each atomic site
    """

    pointops = []

    # Point symmetries of each atom
    for site1 in range(len(struc.sites)):
        tempstruc = struc.copy()

        # Place the origin of the cell at each atomic site
        pointops.append([])

        for site2 in range(len(struc.sites)):  # pylint: disable=C0200
            tempstruc.replace(
                site2,
                tempstruc.sites[site2].specie,
                tempstruc.frac_coords[site2] - struc.frac_coords[site1],
            )

        sgastruc = sga(tempstruc, symprec=precision)
        ops = sgastruc.get_symmetry_operations(cartesian=True)
        for site2, op in enumerate(ops):
            if all(op.translation_vector == [0, 0, 0]):
                pointops[site1].append(op)
    return pointops


def get_shared_symmetry_operations(struc, pointops, tol=0.1):
    """
    Get all the point group operations shared by a pair of atomic sites
    in the form [[point operations of site index 1],[],...,[]]

    Args:
        struc: Pymatgen structure
        pointops: list of point group operations from get_site_symmetries method

    Return:
        list of lists of shared point operations for each pair of atomic sites
    """
    numsites = len(struc)
    sharedops = [[0 for x in range(numsites)] for y in range(numsites)]
    for site1 in range(numsites):
        for site2 in range(numsites):
            sharedops[site1][site2] = []
            for op1, pop1 in enumerate(pointops[site1]):
                for op2, pop2 in enumerate(pointops[site2]):
                    if np.allclose(pop1.rotation_matrix, pop2.rotation_matrix):
                        sharedops[site1][site2].append(pop1)

    for site1, sops in enumerate(sharedops):
        for site2, sop in enumerate(sops):
            uniqueops = []
            for ops in sop:
                op = SymmOp.from_rotation_and_translation(
                    rotation_matrix=ops.rotation_matrix,
                    translation_vec=(0, 0, 0),
                    tol=tol,
                )
                if op not in uniqueops:
                    uniqueops.append(op)

            sharedops[site1][site2] = uniqueops

    return sharedops
