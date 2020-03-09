"""
Provides analysis of site symmetries.
"""

import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
from pymatgen.core.operations import SymmOp


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

        for site2 in range(len(struc.sites)):
            tempstruc.replace(
                site2,
                tempstruc.sites[site2].specie,
                tempstruc.frac_coords[site2] - struc.frac_coords[site1],
            )

        sgastruc = sga(tempstruc, symprec=precision)
        ops = sgastruc.get_symmetry_operations(cartesian=True)
        for site2 in range(len(ops)):
            if all(ops[site2].translation_vector == [0, 0, 0]):
                pointops[site1].append(ops[site2])
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
            for op1 in range(len(pointops[site1])):
                for op2 in range(len(pointops[site2])):
                    if np.allclose(
                        pointops[site1][op1].rotation_matrix,
                        pointops[site2][op2].rotation_matrix,
                    ):
                        sharedops[site1][site2].append(pointops[site1][op1])

    for site1 in range(len(sharedops)):
        for site2 in range(len(sharedops[site1])):
            uniqueops = []
            for ops in range(len(sharedops[site1][site2])):
                op = SymmOp.from_rotation_and_translation(
                    rotation_matrix=sharedops[site1][site2][ops].rotation_matrix,
                    translation_vec=(0, 0, 0),
                    tol=tol,
                )
                if op in uniqueops:
                    continue
                else:
                    uniqueops.append(op)
            sharedops[site1][site2] = uniqueops

    return sharedops
