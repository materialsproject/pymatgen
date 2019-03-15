import pymatgen
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
from pymatgen.core.operations import SymmOp
from pymatgen import Element
import matplotlib.pyplot as plt
from pymatgen.analysis.elasticity.tensors import Tensor

# function that returns the point group symmetries centered on each site
def get_site_symmetries(struc, precision = 0.1):
    numsites = len(struc.sites)
    tempstruc1 = struc.copy()

    sgastruc = sga(struc)
    maptoprim = sgastruc.get_symmetry_dataset()['mapping_to_primitive']
 
    symops = []
    symopsorig = sgastruc.get_symmetry_operations(cartesian = True)

    # point symmetries of each atom
    for j in range(len(struc.sites)):
        tempstruc = struc.copy()

        symops.append([])
        for i in range(len(strucs.sites)):
            tempstruc.replace(i, tempstruc.sites[i].specie, tempstruc.frac_coords[i]-struc.frac_coords[j])


        sgastruc = sga(tempstruc, symprec = precision)
        ops = sgastruc.get_symmetry_operations(cartesian = True)
        for k in range(len(ops)):
            if all(ops[k].translation_vector == [0,0,0]):
                symops[j].append(ops[k])
    return symops

# Function that get's symmetry operations that map one atomic site onto another
def get_equivalent_atom_symmetries(struc, tol = 1e-03):
    #get equivalent atoms
    symstruc = sga(struc).get_symmetrized_structure()
    eq = []
    for i in range(len(symstruc.sites)):
        eqsites = symstruc.find_equivalent_sites(symstruc.sites[i])
        if any(symstruc.sites[i] in sublist for sublist in eq):
            eq.append([])
            continue
        else:
            eq.append(eqsites)

    order = []
    for i in range(len(eq)):
        order.append([])
        for j in range(len(symstruc.sites)):
            if symstruc.sites[j] in eq[i]:
                order[i].append(j)

    while [] in order:
        order.remove([])


    # equivalent atom symmetry operations
    symopsfrac = sga(struc).get_symmetry_operations()
    eqop = []
    eqopfrac = []
    for i in range(numsites):
        eqop.append([])
        eqopfrac.append([])
    for i in range(len(order)): 
        for j in range(len(order[i])):

            site = symstruc.frac_coords[order[i][j]]
            for l in range(len(order[i])):

                eqsite = symstruc.frac_coords[order[i][l]]

                for k in range(len(symopsfrac)):
                    count = 0
                    siteupdate = eqsite
                    dist = np.array([1,1,1])

                    siteupdate = symopsfrac[k].operate(siteupdate)
                    while not all(dist < tol):
                        count = count + 1

                        rounded2 = np.round(np.mod(np.abs(siteupdate - site), 1))
                        dist2 = np.abs(rounded2 - np.mod(np.abs(siteupdate - site), 1))

                        if all(dist2 < tol):
                            eqop[order[i][j]].append([order[i][l], symopsorig[k], count])
                            eqopfrac[order[i][j]].append([order[i][l], symopsfrac[k], count])

                        siteupdate = symopsfrac[k].operate(siteupdate)
                        rounded = np.round(np.mod(np.abs(siteupdate - eqsite), 1))
                        dist = np.abs(rounded - np.mod(np.abs(siteupdate - eqsite), 1))
    return eqop


# Function that returns the symmetry operations shared by a pair of atoms 
def get_shared_symmetry_operations(strucs, symops, tol = 0.1):
    sharedops = [[0 for x in range(numsites)] for y in range(numsites)] 
    for i in range(numsites):
        for j in range(numsites):
            sharedops[i][j] = ([])
            for k in range(len(symops[i])):
                for l in range(len(symops[j])):
                    if np.allclose(symops[i][k].rotation_matrix, symops[j][l].rotation_matrix):
                        sharedops[i][j].append(symops[i][k])

    for i in range(len(sharedops)):
        for k in range(len(sharedops[i])):
            uniqueops = []
            for j in range(len(sharedops[i][k])):
                op = SymmOp.from_rotation_and_translation(
                        rotation_matrix=sharedops[i][k][j].rotation_matrix,
                        translation_vec=(0, 0, 0), tol=tol)
                if op in uniqueops:
                    continue
                else:
                    uniqueops.append(op)
            sharedops[i][k] = uniqueops

    return sharedops



                    

