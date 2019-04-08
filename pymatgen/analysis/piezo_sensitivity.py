from pymatgen.core.tensors import Tensor
from pymatgen.analysis.piezo import BornEffectiveChargeTensor, InternalStrainTensor
from pymatgen.analysis.symmetry import sitesymmeties as ss
import numpy as np
import warnings

try:
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
    from phonopy.file_IO import write_disp_yaml
except ImportError:
    Phonopy = None

__author__ = "Handong Ling"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Handong Ling"
__email__ = "hling@lbl.gov"
__status__ = "Development"
__date__ = "Feb, 2019"

class BornEffectiveChargeTensor(Tensor):
    """
    This class describes the Nx3x3 born effective charge tensor
    """

    def __new__(cls, input_array, tol=1e-3):
        """
        Create an BornEffectiveChargeTensor object.  The constructor throws 
        an error if the input_matrix argument does not satisfy the acoustic 
        sum rule. Note that the constructor uses __new__ rather than __init__ 
        according to the standard method of subclassing numpy ndarrays.

        Args:
            input_matrix (Nx3x3 array-like): the Nx3x63 array-like
                representing the born effective charge tensor
        """
        obj = super(PiezoTensor, cls).__new__(cls, input_array, check_rank=4)
        if not (np.sum(obj)) < tol).all():
            warnings.warn("Input born effective charge tensor does "
                          "not satisfy acoustic sum rule")
        return obj.view(cls)

class InternalStrainTensor(Tensor):
    """
    This class describes the Nx3x3x3 born effective charge tensor
    """

    def __new__(cls, input_array, tol=1e-3):
        """
        Create an InternalStrainTensor object.  The constructor throws an error if
        the shape of the input_matrix argument is not Nx3x3x3, i. e. in true
        tensor notation. Note that the constructor uses __new__ rather than
        __init__ according to the standard method of subclassing numpy
        ndarrays.

        Args:
            input_matrix (Nx3x3x3 array-like): the Nx3x3x3 array-like
                representing the internal strain tensor
        """
        obj = super(PiezoTensor, cls).__new__(cls, input_array, check_rank=4)
        if not (obj - np.transpose(obj, (0, 1, 3, 2)) < tol).all():
            warnings.warn("Input internal strain tensor does "
                          "not satisfy standard symmetries")
        return obj.view(cls)

class ForceConstantMatrix(Tensor):
    """
    This class describes the NxNx3x3 born effective charge tensor
    """

    def __new__(cls, input_array, tol=1e-3):
        """
        Create an ForceConstantMatrix object.  The constructor throws an error if
        the shape of the input_matrix argument is not NxNx3x3, i. e. in true
        tensor notation. Note that the constructor uses __new__ rather than
        __init__ according to the standard method of subclassing numpy
        ndarrays.

        Args:
            input_matrix (NxNx3x3 array-like): the NxNx3x3 array-like
                representing the force constant matrix
        """
        obj = super(PiezoTensor, cls).__new__(cls, input_array, check_rank=4)

        return obj.view(cls)


def get_rand_BEC(struc, eqop, symops, max_charge=1):
    """
    Generate a random born effective charge tensor which obeys a structure's 
    symmetry and the acoustic sum rule

    Args:
        struc (pymatgen structure): 
        eqop: operations which map atoms to equivalent atoms 
            in the structure
        symops: point symmetry operations of each atomic site
        max_charge (float): maximum born effective charge value

    Return:
        BornEffectiveChargeTensor object
    """

    symstruc = sga(struc)
    maptoprim = sgastruc.get_symmetry_dataset()['mapping_to_primitive']
    eqatoms = sgastruc.get_symmetry_dataset()['equivalent_atoms']
    
    l = len(struc)
    BEC = np.zeros((l,3,3))
    primsites = []
    for i in range(l):
        good = 0    
        for j in reversed(range(i)):
            if maptoprim[j] == maptoprim[i]:
                BEC[i] = BEC[j]
                good = 1
                break
            if eqatoms[j] == eqatoms[i] and maptoprim[j] != maptoprim[i]:
                for k in range(len(eqop[i])):
                    if eqop[i][k][0] == j:
                        tempbec = BEC[j]
                        for m in range(eqop[i][k][2]):
                            tempbec = eqop[i][k][1].transform_tensor(tempbec)
                            good = 1

                        
                        BEC[i] = tempbec
                    if good == 1:
                        break
            if good == 1:
                break
        if good == 1:
            continue
                
                
        symsite = 0 
        temp_tensor = Tensor(np.random.rand(3,3)-0.5)         
        temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in symops[i]]) \
            /len(symops[i])
            
        BEC[i] = np.array(temp_tensor)

    count = 0

    
#     Enforce Acoustic Sum
    disp_charge = np.einsum("ijk->jk",BEC)/l
    
    add = np.zeros([l,3,3])
            
    for i in range(l):
        add[i] = disp_charge
        
        eq = 0
        for j in range(i):
            if maptoprim[i] == maptoprim[j] and i != j:
                add[i] = add[j]
                eq = 1
                break

        if eq == 1:
            continue

        good = 0    
        for j in reversed(range(i)):

            if eqatoms[j] == eqatoms[i]:
                for k in range(len(eqop[i])):
                    if eqop[i][k][0] == j:
                        tempbec = add[j]
                        for m in range(eqop[i][k][2]):
                            tempbec = eqop[i][k][1].transform_tensor(tempbec)
                            good = 1
                            break
                        add[i] = tempbec
                    if good == 1:
                        break
            if good == 1:
                break
        if good == 1:
            continue

        add[i] = sum([Tensor(add[i]).transform(symm_op) for symm_op in symops[i]]) / len(symops[i])
    
    BEC = BEC - add
    
    disp_charge = np.einsum("ijk->jk",BEC)/l
 
    return BornEffectiveChargeTensor(BEC*max_charge)

def get_rand_IST(struc, eqop, symops, max_force=1):
    """
    Generate a random internal strain tensor which obeys a structure's 
    symmetry and the acoustic sum rule

    Args:
        struc (pymatgen structure): 
        eqop: operations which map atoms to equivalent atoms 
            in the structure
        symops: point symmetry operations of each atomic site
        max_charge (float): maximum born effective charge value

    Return:
        InternalStrainTensor object
    """
    symstruc = sga(struc)
    maptoprim = sgastruc.get_symmetry_dataset()['mapping_to_primitive']
    eqatoms = sgastruc.get_symmetry_dataset()['equivalent_atoms']
    
    l = len(struc)
    IST = np.zeros((l,3,3,3))
    for i in range(l):
        good = 0    
        for j in reversed(range(i)):
            if maptoprim[j] == maptoprim[i]:
                IST[i] = IST[j]
                good = 1
                break
            if eqatoms[j] == eqatoms[i]:
                
                for k in range(len(eqop[i])):
                    if eqop[i][k][0] == j:
                        tempist = IST[j]
                        for m in range(eqop[i][k][2]):
                            tempist = eqop[i][k][1].transform_tensor(tempist)
                            good = 1
                        IST[i] = tempist
                    if good == 1:
                        break
            if good == 1:
                break
        if good == 1:
            continue
                
        temp_tensor = Tensor(np.random.rand(3,3,3)-0.5)    
        for m in range(3):
            for n in range(3):
                for o in range(n,3):
                    temp_tensor[m][o][n] = temp_tensor[m][n][o]
        temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in symops[i]]) \
            /len(symops[i])
 
        IST[i] = np.array(temp_tensor)
    count = 0

    return InternalStrainTensor(IST*max_force)


def get_FCM_symmetry(fcm, operations, sharedops, max_force = 1):
    """
    Generate a symmeterized force constant matrix from an unsymmeterized matrix

    Args:
        fcm (numpy array): unsymmeterized force constant matrix 
        operations: list of operation mappings for indexed sites in the force
            constant matrix
        sharedops: symmetry operations share by a pair of atomic sites
        max_charge (float): maximum born effective charge value

    Return:
        numpy array representing the force constant matrix
    """

    # set max force in reciprocal space
    numsites = int(len(fcm/3))
    D = fcm
    passed = []
    #symmetry operations
    for i in range(len(operations)):
        if type(operations[i][4]) == pymatgen.core.operations.SymmOp:

            tempfcm = D[3*operations[i][2]:3*operations[i][2]+3, 3*operations[i][3]:3*operations[i][3]+3]
            tempfcm = operations[i][4].transform_tensor(tempfcm)

            D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] = tempfcm
            D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = tempfcm.T

            continue

        temp_tensor = Tensor(D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3])
        temp_tensor_sum = sum([temp_tensor.transform(symm_op) for symm_op in \
                               sharedops[operations[i][0]][operations[i][1]]])
        temp_tensor_sum = temp_tensor_sum/(len(sharedops[operations[i][0]][operations[i][1]]))
        
        for pair in range(len(operations[i][5])):
            
            temp_tensor2 = temp_tensor_sum.T
            temp_tensor2 = operations[i][5][pair].transform_tensor(temp_tensor2)
            temp_tensor_sum = (temp_tensor_sum + temp_tensor2)/2

        D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] = temp_tensor_sum
        D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = temp_tensor_sum.T

    return(D)


def get_stable_fcm(fcm, struc, operations, sharedops, asum = 10):
    
    """
    Generate a symmeterized force constant matrix 
    that has no unstable modes and also obeys the acoustic sum rule through an
    iterative procedure

    Args:
        fcm (numpy array): unsymmeterized force constant matrix 
        operations: list of operation mappings for indexed sites in the force
            constant matrix
        sharedops: 
        max_charge (float): maximum born effective charge value
        asum (int): number of iterations to attempt to obey the acoustic sum
            rule

    Return:
        numpy array representing the force constant matrix
    """

    numsites = int(len(fcm)/3)
    
    
    check = 0
    count = 0
    while check == 0:
        if count > 20:
            check = 1

        eigs, vecs = np.linalg.eig(fcm)

        maxeig = np.max(np.abs(eigs))
        mineig = np.min(eigs)
        eigsort = np.argsort(np.abs(eigs))
        for i in range(3,len(eigs)):
            if eigs[eigsort[i]] > 1e-06:
                eigs[eigsort[i]] = -1*maxeig*np.random.rand()
        diag = np.real(np.eye(len(fcm))*eigs)

        fcm = np.matmul(np.matmul(vecs, diag), vecs.T)
        fcm = get_FCM_symmetry(fcm, operations, sharedops)
        fcm = get_asum_fcm(fcm, operations, sharedops, asum)
        #symmetry operations
        eigs, vecs = np.linalg.eig(fcm)
        a = 0
        eigsort = np.argsort(np.abs(eigs))
        for i in range(3,len(eigs)):
            if eigs[eigsort[i]] > 1e-06:
                a = 1
        if a == 1:
            count = count + 1 
            continue
        else: 
            check = 1
        
        
    return fcm


    #acoustic sum

def get_asum_fcm(fcm, operations, sharedops, numiter):
    """
    Generate a symmeterized force constant matrix that also obeys the acoustic 
    sum rule through an iterative procedure

    Args:
        fcm (numpy array): unsymmeterized force constant matrix 
        operations: list of operation mappings for indexed sites in the force
            constant matrix
        sharedops: 
        max_charge (float): maximum born effective charge value
        numiter (int): number of iterations to attempt to obey the acoustic sum
            rule

    Return:
        numpy array representing the force constant matrix
    """

    # set max force in reciprocal space 
    numsites = int(len(fcm)/3)
    
    D = np.ones([numsites*3, numsites*3])
    for num in range(numiter):
        X = np.real(fcm)
        passed = []
        
        #symmetry operations
        pastrow = 0 
        total = np.zeros([3,3])
        for col in range(numsites):
            total = total + X[0:3, col*3:col*3+3]

        total = total/(numsites)
        for i in range(len(operations)):
            prim = 0
            if type(operations[i][4]) == pymatgen.core.operations.SymmOp:
                
                tempfcm = D[3*operations[i][2]:3*operations[i][2]+3, 3*operations[i][3]:3*operations[i][3]+3]
                
                tempfcm = operations[i][4].transform_tensor(tempfcm)

                D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] = tempfcm
                D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = tempfcm.T
                
                continue
                
            currrow = operations[i][0]
            if currrow != pastrow:
                total = np.zeros([3,3])
                for col in range(numsites):
                    total = total + X[currrow*3:currrow*3+3, col*3:col*3+3]
                for col in range(currrow):
                    total = total - D[currrow*3:currrow*3+3, col*3:col*3+3]
                total = total/(numsites-currrow)
            pastrow = currrow 
            

            temp_tensor = Tensor(total)
            temp_tensor_sum = sum([temp_tensor.transform(symm_op) for symm_op in \
                                   sharedops[operations[i][0]][operations[i][1]]])
        
            if len(sharedops[operations[i][0]][operations[i][1]]) != 0:
                temp_tensor_sum = temp_tensor_sum/(len(sharedops[operations[i][0]][operations[i][1]]))
                
            numops = 1
            count = 0
            for pair in range(len(operations[i][5])):

                temp_tensor2 = temp_tensor_sum.T
                temp_tensor2 = operations[i][5][pair].transform_tensor(temp_tensor2)
                temp_tensor_sum = (temp_tensor_sum + temp_tensor2)/2

            D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] = temp_tensor_sum
            D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = temp_tensor_sum.T
        fcm = fcm - D
            
    return(fcm)


def get_fullfcm(struc, fcm, force = 1, asum = 10):    
    """
    Generate a symmeterized force constant matrix from an unsymmeterized matrix
    that has no unstable modes and also obeys the acoustic sum rule through an 
    iterative procedure

    Args:
        fcm (numpy array): unsymmeterized force constant matrix 
        operations: list of operation mappings for indexed sites in the force
            constant matrix
        sharedops: 
        force (float): maximum force constant
        asum (int): number of iterations to attempt to obey the acoustic sum
            rule

    Return:
        force constant matrix representing the force constant matrix
    """
    numsites = len(struc.sites)
    structure = phonopy.get_phonopy_structure(struc)
    operations = get_fcm_operations(struc, fcm)
    pnstruc = pn.Phonopy(structure, np.eye(3), np.eye(3))
    symops = ss.get_site_symmetries(struc)
    sharedops = ss.get_shared_symmetry_operations(strucs, symops)

    fcm  = get_FCM2(struc, sharedops, operations, force)
    fcm = get_asum_fcm(fcm, operations, sharedops, numiter = asum)
    fcm = get_stable_fcm(fcm, struc, operations, sharedops)
    
    
    fcm = np.reshape(fcm, (numsites,3,numsites,3)).swapaxes(1,2)
    dynmass= np.zeros([len(struc),len(struc),3,3])
    masses = []
    for j in range(numsites):
        masses.append(struc.sites[j].specie.atomic_mass)
    fcmsmass = np.zeros([numsites,numsites,3,3])
    for m in range(numsites):
        for n in range(numsites):
            dynmass[m][n] = fcm[m][n]*np.sqrt(masses[m])*np.sqrt(masses[n])

    fcm = dynmass
    fcm = np.reshape(np.swapaxes(fcm,1,2),(numsites*3,numsites*3))
    supercell = pnstruc.get_supercell()
    primitive = pnstruc.get_primitive()
    
    converter = dyntofc.DynmatToForceConstants(primitive, supercell)
    
    converter.set_dynamical_matrices(dynmat = [fcm])

    try:
        converter.run()
        fc = converter.get_force_constants()
        fc = np.reshape(np.swapaxes(fc,1,2),(len(fc)*3,len(fc)*3))
        fc = get_FCM_symmetry(fc, operations, sharedops)
    except:
        fc = converter.get_force_constants()
        fc = np.reshape(np.swapaxes(fc,1,2),(len(fc)*3,len(fc)*3))
        fc = get_FCM_symmetry(fc, operations, sharedops)
    return ForceConstantMatrix(fc)

def get_fcm_operations(struc, fcm):
    """
    Generate a mapping of symmetry operations which the force constant matrix must obey
    Args:
        struc (pymatgen structure): structure whose symmetry operations the piezo
            tensor must obey
        fcm (np.array): example force constant matrix whose symmetry operations 
            the structure must obey
    Return:
        PiezoTensor object
    """

#     fcm = np.reshape(fcm, (len(struc),3,len(struc),3)).swapaxes(1,2)
    ops = sga(struc).get_symmetry_operations(cartesian = True)
    structure = phonopy.get_phonopy_structure(struc)

    pnstruc = pn.Phonopy(structure, np.eye(3), np.eye(3))
    supercell = pnstruc.get_supercell()
    primitive = pnstruc.get_primitive()

    
    converter = dyntofc.DynmatToForceConstants(primitive, supercell)
    pnstruc.set_force_constants(fcm)
    dyn = pnstruc.get_dynamical_matrix_at_q([0,0,0])
    dyn = np.reshape(dyn, (len(struc),3,len(struc),3)).swapaxes(1,2)
    passed = []
    relations = []
    for i in range(len(dyn)):
        for j in range(i, len(dyn)):
            eig1, vecs1 = np.linalg.eig(dyn[i][j])
    #         eig2, vecs2 = np.linalg.eig(dyn[j][i])
#             print(i,j, eig)
            index = np.argsort(eig1)
            neweig = np.real([eig1[index[0]],eig1[index[1]],eig1[index[2]]])
            passed.append([i,j,np.real(neweig)])
    #         passed.append([j,i,eig2])
            for k in range(len(passed)):
                if np.allclose(neweig,passed[k][2], atol = 1e-05):
                    relations.append([i,j, passed[k][0], passed[k][1]])
                    break
    operations = []
    for i in range(len(relations)):
        good = 0
        operations.append(relations[i])
        if relations[i][2] == relations[i][0] and relations[i][1] == relations[i][3]:
            operations[i].append([])
            transpose = []
            
            for j in range(len(ops)):
                new = ops[j].transform_tensor(dyn[relations[i][2]][relations[i][3]])
                if np.allclose(new, dyn[relations[i][1]][relations[i][0]], atol = 1e-05):
                    transpose.append(ops[j])
            operations[i].append(transpose)        
            transpose = []
            continue
        
        for j in range(len(ops)):
            new = ops[j].transform_tensor(dyn[relations[i][2]][relations[i][3]])
            if np.allclose(new, dyn[relations[i][0]][relations[i][1]], atol = 1e-05):
                operations[i].append(ops[j])
                good = 1

            elif np.allclose(new, dyn[relations[i][0]][relations[i][1]].T, atol = 1e-05):
                operations[i].append(ops[j])
                temp1 = operations[i][2]
                temp2 = operations[i][3]
                operations[i][2] = temp2
                operations[i][3] = temp1
                good = 1
            if good ==1:
                transpose = []
            
                for j in range(len(ops)):
                    new = ops[j].transform_tensor(dyn[relations[i][2]][relations[i][3]])
                    if np.allclose(new, dyn[relations[i][1]][relations[i][0]], atol = 1e-05):
                        transpose.append(ops[j])
                operations[i].append(transpose)
                break
        
    return operations

                    
def rand_piezo(struc, fcm, anumiter = 10):
    """
    Generate a random piezoelectric tensor based on a structure and corresponding
    symmetry

    Args:
        struc (pymatgen structure): structure whose symmetry operations the piezo
            tensor must obey
        fcm (np.array): example force constant matrix whose symmetry operations 
            the structure must obey
    Return:
        PiezoTensor object
    """
    symops = ss.get_site_symmetries(struc)
    eqop = ss.get_equivalent_atom_symmetrxies(struc)
    sharedops = ss.get_shared_symmetry_operations(strucs, symops)
    numsites = len(struc.sites)
    BEC = get_rand_BEC(struc, eqop, symops, max_charge=20)
    IST = get_rand_IST(struc, eqop, symops, max_force=10)
    FCM = get_fullfcm(struc, fcm, sharedops, anumiter)
    FCM = np.reshape(FCM, (numsites,3,numsites,3)).swapaxes(1,2)
    P = PiezoTensor(get_piezo(BEC,FCM,IST)*16.0216559424/struc.volume)
    return (BEC,IST,FCM,P)
