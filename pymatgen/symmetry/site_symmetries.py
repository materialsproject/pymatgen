import pymatgen
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
from pymatgen.core.operations import SymmOp
from pymatgen import Element
from pymatgen.analysis.elasticity.tensors import Tensor

# function that returns the point group symmetries centered on each site
def get_site_symmetries(struc, precision = 0.1):
    numsites = len(struc.sites)
    sgastruc = sga(struc)
 
    pointops = []

    # Point symmetries of each atom
    for site1 in range(len(struc.sites)):
        tempstruc = struc.copy()

        # Place the origin of the cell at each atomic site
        pointops.append([])
        for site2 in range(len(struc.sites)):
            tempstruc.replace(site2, tempstruc.sites[site2].specie, tempstruc.frac_coords[site2]-struc.frac_coords[site1])


        sgastruc = sga(tempstruc, symprec = precision)
        ops = sgastruc.get_symmetry_operations(cartesian = True)
        for site2 in range(len(ops)):
            if all(ops[site2].translation_vector == [0,0,0]):
                pointops[site1].append(ops[site2])
    return pointops


# Function that returns the symmetry operations shared by a pair of atoms 
def get_shared_symmetry_operations(struc, pointops, tol = 0.1):
    numsites = len(struc)
    sharedops = [[0 for x in range(numsites)] for y in range(numsites)] 
    for site1 in range(numsites):
        for site2 in range(numsites):
            sharedops[site1][site2] = ([])
            for op1 in range(len(pointops[site1])):
                for op2 in range(len(pointops[site2])):
                    if np.allclose(pointops[site1][op1].rotation_matrix, pointops[site2][op2].rotation_matrix):
                        sharedops[site1][site2].append(pointops[site1][op1])

    for site1 in range(len(sharedops)):
        for site2 in range(len(sharedops[site1])):
            uniqueops = []
            for ops in range(len(sharedops[site1][site2])):
                op = SymmOp.from_rotation_and_translation(
                        rotation_matrix=sharedops[site1][site2][ops].rotation_matrix,
                        translation_vec=(0, 0, 0), tol=tol)
                if op in uniqueops:
                    continue
                else:
                    uniqueops.append(op)
            sharedops[site1][site2] = uniqueops

    return sharedops



                    

