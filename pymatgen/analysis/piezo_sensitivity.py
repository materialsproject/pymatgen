from pymatgen.core.tensors import Tensor
from pymatgen.symmetry import site_symmetries as ss
from pymatgen.analysis.piezo import PiezoTensor
import numpy as np
import warnings
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
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
        obj = super(Tensor, cls).__new__(cls, input_array, check_rank=4)
        # if not (np.sum(obj)) < tol).all():
        #     warnings.warn("Input born effective charge tensor does "
        #                   "not satisfy charge neutrality")
        return obj.view(cls)

    def get_BEC_IST_operations(struc, bec, pointops,eigtol = 1e-05, opstol = 1e-03):
        ops = sga(struc).get_symmetry_operations(cartesian = True)
        uniquepointops = [] 
        for i in range(len(ops)):
            uniquepointops.append(ops[i])
        
        for i in range(len(pointops)):
            for j in range(len(pointops[i])):
                if pointops[i][j] not in uniquepointops:
                    uniquepointops.append(pointops[i][j])
                
        structure = phonopy.get_phonopy_structure(struc)

        pnstruc = pn.Phonopy(structure, np.eye(3), np.eye(3))
        supercell = pnstruc.get_supercell()
        primitive = pnstruc.get_primitive()

        
        converter = dyntofc.DynmatToForceConstants(primitive, supercell)
        pnstruc.set_force_constants(fcm)
         
        passed = []
        relations = []
        for i in range(len(bec)):
            unique = 1
            eig1, vecs1 = np.linalg.eig(bec[i])
            index = np.argsort(eig1)
            neweig = np.real([eig1[index[0]],eig1[index[1]],eig1[index[2]]])
            for k in range(len(passed)):

                if np.allclose(neweig, passed[k][1], atol = eigtol):
                    relations.append([i,k])
                    unique = 0 
                    passed.append([i,passed[k][0], neweig])
                    break
            if unique == 1:
                relations.append([i, i])
                passed.append([i, neweig])
        BEC_IST_operations = []
        for i in range(len(relations)):
            good = 0
            BEC_IST_operations.append(relations[i])
            BEC_IST_operations[i].append([])

            good = 0
            for j in range(len(uniquepointops)):
                new = uniquepointops[j].transform_tensor(bec[relations[i][1]])
                
                ## Check the matrix it references
                if np.allclose(new, bec[relations[i][0]], atol = opstol):
                    BEC_IST_operations[i][2].append(uniquepointops[j])

        
        return BEC_IST_operations

                        

    def get_rand_BEC(struc, symops, BEC_operations, max_charge=1):
        """
        Generate a random born effective charge tensor which obeys a structure's 
        symmetry and the acoustic sum rule

        Args:
            struc (pymatgen structure): 
            symops: point symmetry operations of each atomic site
            BEC_operations: list of operations which map atomic sites onto each other
            max_charge (float): maximum born effective charge value

        Return:
            BornEffectiveChargeTensor object
        """
        symstruc = sga(struc)
        symstruc = symstruc.get_symmetrized_structure()
        
        
        l = len(struc)
        BEC = np.zeros((l,3,3))
        primsites = []
        for i in range(len(BEC_operations)):
            if BEC_operations[i][0] == BEC_operations[i][1]:
                temp_tensor = Tensor(np.random.rand(3,3)-0.5)         
                temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in symops[i]]) \
                    /len(symops[i])
                BEC[i] = temp_tensor
            else: 
                tempfcm = np.zeros([3,3])
                for j in range(len(BEC_operations[i][2])):

                    tempfcm += BEC_operations[i][2][j].transform_tensor(BEC[BEC_operations[i][1]])
                BEC[BEC_operations[i][0]] = tempfcm
                if len(BEC_operations[i][2]) != 0:
                    BEC[BEC_operations[i][0]] = BEC[BEC_operations[i][0]]/len(BEC_operations[i][2])
        
    #     Enforce Acoustic Sum
        disp_charge = np.einsum("ijk->jk",BEC)/l
        add = np.zeros([l,3,3])
        
        for i in range(len(BEC_operations)):  
            
            if BEC_operations[i][0] == BEC_operations[i][1]:
                temp_tensor = Tensor(disp_charge)
                temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in symops[i]]) \
                    /len(symops[i])
                add[BEC_operations[i][0]] = temp_tensor
            else: 
                temp_tensor = np.zeros([3,3])
                for j in range(len(BEC_operations[i][2])):

                    temp_tensor += BEC_operations[i][2][j].transform_tensor(add[BEC_operations[i][1]])
                
                add[BEC_operations[i][0]] = temp_tensor
                
                if len(BEC_operations[i][2]) != 0:
                    add[BEC_operations[i][0]] = add[BEC_operations[i][0]]/len(BEC_operations[i][2])
            
            BEC = BEC - add
         

        return BornEffectiveChargeTensor(BEC*max_charge)




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
        obj = super(Tensor, cls).__new__(cls, input_array, check_rank=4)
        if not (obj - np.transpose(obj, (0, 1, 3, 2)) < tol).all():
            warnings.warn("Input internal strain tensor does "
                          "not satisfy standard symmetries")
        return obj.view(cls)

    def get_BEC_IST_operations(struc, bec, pointops,eigtol = 1e-05, opstol = 1e-03):
        ops = sga(struc).get_symmetry_operations(cartesian = True)
        uniquepointops = [] 
        for i in range(len(ops)):
            uniquepointops.append(ops[i])
        
        for i in range(len(pointops)):
            for j in range(len(pointops[i])):
                if pointops[i][j] not in uniquepointops:
                    uniquepointops.append(pointops[i][j])
                
        structure = phonopy.get_phonopy_structure(struc)

        pnstruc = pn.Phonopy(structure, np.eye(3), np.eye(3))
        supercell = pnstruc.get_supercell()
        primitive = pnstruc.get_primitive()

        
        converter = dyntofc.DynmatToForceConstants(primitive, supercell)
        pnstruc.set_force_constants(fcm)
         
        passed = []
        relations = []
        for i in range(len(bec)):
            unique = 1
            eig1, vecs1 = np.linalg.eig(bec[i])
            index = np.argsort(eig1)
            neweig = np.real([eig1[index[0]],eig1[index[1]],eig1[index[2]]])
            for k in range(len(passed)):

                if np.allclose(neweig, passed[k][1], atol = eigtol):
                    relations.append([i,k])
                    unique = 0 
                    passed.append([i,passed[k][0], neweig])
                    break
            if unique == 1:
                relations.append([i, i])
                passed.append([i, neweig])
        BEC_IST_operations = []
        for i in range(len(relations)):
            good = 0
            BEC_IST_operations.append(relations[i])
            BEC_IST_operations[i].append([])

            good = 0
            for j in range(len(uniquepointops)):
                new = uniquepointops[j].transform_tensor(bec[relations[i][1]])
                
                ## Check the matrix it references
                if np.allclose(new, bec[relations[i][0]], atol = opstol):
                    BEC_IST_operations[i][2].append(uniquepointops[j])

        
        return BEC_IST_operations    

    def get_rand_IST(struc, symops, IST_operations, max_force=1):
        """
        Generate a random internal strain tensor which obeys a structure's 
        symmetry and the acoustic sum rule

        Args:
            struc (pymatgen structure): 
            symops: point symmetry operations of each atomic site
            IST_operations: list of operations which map atomic sites onto each other
            max_charge (float): maximum born effective charge value

        Return:
            InternalStrainTensor object
        """
        symstruc = sga(struc)
        symstruc = symstruc.get_symmetrized_structure()
        
        
        l = len(struc)
        IST = np.zeros((l,3,3,3))
        primsites = []
        for i in range(len(IST_operations)):
            if IST_operations[i][0] == IST_operations[i][1]:
                temp_tensor = Tensor(np.random.rand(3,3,3)-0.5) 
                for dim in range(3):
                    temp_tensor[dim] = (temp_tensor[dim]+temp_tensor[dim].T)/2
                temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in symops[i]]) \
                    /len(symops[i])
                IST[IST_operations[i][0]] = temp_tensor
            else: 
                temp_tensor = np.zeros([3,3,3])
                for j in range(len(IST_operations[i][2])):

                    temp_tensor += IST_operations[i][2][j].transform_tensor(IST[IST_operations[i][1]])
                
                IST[IST_operations[i][0]] = temp_tensor
                if len(IST_operations[i][2]) != 0:
                    IST[IST_operations[i][0]] = IST[IST_operations[i][0]]/len(IST_operations[i][2])
        
        return IST*max_force



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
        obj = super(Tensor, cls).__new__(cls, input_array, check_rank=3)

        return obj.view(cls)

    def get_FCM_operations(struc, fcm, pointops,eigtol = 1e-05, opstol = 1e-05):
        ops = sga(struc).get_symmetry_operations(cartesian = True)
        uniquepointops = [] 
        for i in range(len(ops)):
            uniquepointops.append(ops[i])
        
        for i in range(len(pointops)):
            for j in range(len(pointops[i])):
                if pointops[i][j] not in uniquepointops:
                    uniquepointops.append(pointops[i][j])
                
        dyn = fcm
        passed = []
        relations = []
        for i in range(len(dyn)):
            for j in range(i, len(dyn)):
                unique = 1
                eig1, vecs1 = np.linalg.eig(dyn[i][j])
                index = np.argsort(eig1)
                neweig = np.real([eig1[index[0]],eig1[index[1]],eig1[index[2]]])
                
                for k in range(len(passed)):
                    if np.allclose(neweig,passed[k][2], atol = eigtol):
                        relations.append([i,j, passed[k][0], passed[k][1]])
                        unique = 0 
                        break
                if unique == 1:
                    relations.append([i,j, j, i])
                    passed.append([i,j,np.real(neweig)])
        FCM_operations = []
        for i in range(len(relations)):
            good = 0
            FCM_operations.append(relations[i])
            FCM_operations[i].append([])

            good = 0
            for j in range(len(uniquepointops)):
                new = uniquepointops[j].transform_tensor(dyn[relations[i][2]][relations[i][3]])
                
                if np.allclose(new, dyn[relations[i][0]][relations[i][1]], atol = opstol):
                    FCM_operations[i][4].append(uniquepointops[j])
                    good = 1
            if relations[i][0] == relations[i][3] and relations[i][1] == relations[i][2]:
                good = 1
            if relations[i][0] == relations[i][2] and relations[i][1] == relations[i][3]:
                good = 1
            if good == 0:
                FCM_operations[i] = [relations[i][0], relations[i][1], relations[i][3], relations[i][2]]
                FCM_operations[i].append([])
                for j in range(len(uniquepointops)):
                    new = uniquepointops[j].transform_tensor(dyn[relations[i][2]][relations[i][3]])
                    if np.allclose(new.T, dyn[relations[i][0]][relations[i][1]], atol = opstol):
                        FCM_operations[i][4].append(uniquepointops[j])
            
        return FCM_operations

                    

    def get_new_FCM(struc, sharedops, operations, max_force = 1):
        # set max force in reciprocal space
        numsites = len(struc.sites)
        D = (1/max_force)*2*(np.ones([numsites*3, numsites*3]))
        for i in range(len(operations)):
            same = 0
            transpose = 0
            if operations[i][0] == operations[i][1] and operations[i][0] == operations[i][2] and operations[i][0] == operations[i][3]:
                same = 1
            if operations[i][0] == operations[i][3] and operations[i][1] == operations[i][2]:
                transpose = 1
            if transpose ==0 and same == 0:
                D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] = np.zeros([3,3])
                D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = np.zeros([3,3])

                for j in range(len(operations[i][4])):

                    tempfcm = D[3*operations[i][2]:3*operations[i][2]+3, 3*operations[i][3]:3*operations[i][3]+3]
                    tempfcm = operations[i][4][j].transform_tensor(tempfcm)
                    D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] += tempfcm

                if len(operations[i][4]) != 0:
                    D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] =\
                        D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3]/len(operations[i][4])
                
                D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = \
                D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3].T
                continue
            else:
                temp_tensor = Tensor(np.random.rand(3,3)-0.5)*max_force

                temp_tensor_sum = sum([temp_tensor.transform(symm_op) for symm_op in \
                                       sharedops[operations[i][0]][operations[i][1]]])
                temp_tensor_sum = temp_tensor_sum/(len(sharedops[operations[i][0]][operations[i][1]]))
                if operations[i][0] != operations[i][1]:
                    for pair in range(len(operations[i][4])):

                        temp_tensor2 = temp_tensor_sum.T
                        temp_tensor2 = operations[i][4][pair].transform_tensor(temp_tensor2)
                        temp_tensor_sum = (temp_tensor_sum + temp_tensor2)/2
                        
                else:
                    temp_tensor_sum = (temp_tensor_sum + temp_tensor_sum.T)/2

                D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] = temp_tensor_sum
                D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = temp_tensor_sum.T
                

        return D
    def get_symmetrized_FCM(fcm, operations, sharedops, max_force = 1):
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


    def get_stable_FCM(fcm, struc, operations, sharedops, asum = 10):
        
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
            if count > 100:
                check = 1
            
            eigs, vecs = np.linalg.eig(fcm)

            maxeig = np.max(-1*eigs)
            eigsort = np.argsort(np.abs(eigs))
            for i in range(3,len(eigs)):
                if eigs[eigsort[i]] > 1e-06:
                    eigs[eigsort[i]] = -1*maxeig*np.random.rand()
            diag = np.real(np.eye(len(fcm))*eigs)

            fcm = np.matmul(np.matmul(vecs, diag), vecs.T)
            fcm = get_FCM_symmetry(fcm, operations, sharedops)
            fcm = get_asum3(fcm, operations, sharedops, asum)
            #symmetry operations
            eigs, vecs = np.linalg.eig(fcm)
            unstable_modes = 0
            eigsort = np.argsort(np.abs(eigs))
            for i in range(3,len(eigs)):
                if eigs[eigsort[i]] > 1e-06:
                    unstable_modes = 1
            if unstable_modes == 1:
                count = count + 1 
                continue
            else: 
                check = 1
            
            
        return fcm


    #acoustic sum

    def get_asum_FCM(fcm, operations, sharedops, numiter):
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
                same = 0
                transpose = 0
                if operations[i][0] == operations[i][1] and operations[i][0] == operations[2] and operations[i][0] == operations[i][3]:
                    same = 1
                if operations[i][0] == operations[i][3] and operations[i][1] == operations[i][2]:
                    transpose = 1
                if transpose ==0 and same == 0:
                    D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] = np.zeros([3,3])
                
                    for j in range(len(operations[i][4])):

                        tempfcm = D[3*operations[i][2]:3*operations[i][2]+3, 3*operations[i][3]:3*operations[i][3]+3]
                        tempfcm = operations[i][4][j].transform_tensor(tempfcm)

                        D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] += tempfcm
                    
                    if len(operations[i][4]) != 0:
                        D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] =\
                            D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3]/len(operations[i][4])
                    D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = \
                        D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3].T
                    continue
                else:
                    ## Get the difference in the sum up to this point
                    currrow = operations[i][0]
                    if currrow != pastrow:
                        total = np.zeros([3,3])
                        for col in range(numsites):
                            total = total + X[currrow*3:currrow*3+3, col*3:col*3+3]
                        for col in range(currrow):
                            total = total - D[currrow*3:currrow*3+3, col*3:col*3+3]
                        total = total/(numsites-currrow)
                    pastrow = currrow 

                    ## Apply the point symmetry operations of the site
                    temp_tensor = Tensor(total)
                    temp_tensor_sum = sum([temp_tensor.transform(symm_op) for symm_op in \
                                           sharedops[operations[i][0]][operations[i][1]]])

                    if len(sharedops[operations[i][0]][operations[i][1]]) != 0:
                        temp_tensor_sum = temp_tensor_sum/(len(sharedops[operations[i][0]][operations[i][1]]))



                    ## Apply the proper transformation if there is an equivalent already
                    if operations[i][0] != operations[i][1]:

                        for pair in range(len(operations[i][4])):

                            temp_tensor2 = temp_tensor_sum.T
                            temp_tensor2 = operations[i][4][pair].transform_tensor(temp_tensor2)
                            temp_tensor_sum = (temp_tensor_sum + temp_tensor2)/2
                        
                    else:
                        temp_tensor_sum = (temp_tensor_sum + temp_tensor_sum.T)/2

                    D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3] = temp_tensor_sum
                    D[3*operations[i][1]:3*operations[i][1]+3, 3*operations[i][0]:3*operations[i][0]+3] = temp_tensor_sum.T
            fcm = fcm - D

        
        return(fcm)

    def get_FCM(struc, fcm, force = 1, asum = 10):    
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
        operations = FCM_operations
        pnstruc = pn.Phonopy(structure, np.eye(3), np.eye(3))
        
        dyn = get_new_FCM(struc, sharedops, operations, force)
        dyn = get_asum_FCM(dyn, operations, sharedops, numiter = asum)
        dyn = get_stable_FCM(dyn, struc, operations, sharedops)
        
        dyn = np.reshape(dyn, (numsites,3,numsites,3)).swapaxes(1,2)
        

        dynmass= np.zeros([len(struc),len(struc),3,3])
        masses = []
        for j in range(numsites):
            masses.append(struc.sites[j].specie.atomic_mass)
        dynmass = np.zeros([numsites,numsites,3,3])
        for m in range(numsites):
            for n in range(numsites):
                dynmass[m][n] = dyn[m][n]*np.sqrt(masses[m])*np.sqrt(masses[n])
    #     pnstruc.set_force_constants(fcm2)

        fcm = dynmass
        fcm = np.reshape(np.swapaxes(fcm,1,2),(numsites*3,numsites*3))

        fcm = dynmass
        fcm = np.reshape(np.swapaxes(fcm,1,2),(numsites*3,numsites*3))
        supercell = pnstruc.get_supercell()
        primitive = pnstruc.get_primitive()
        
        converter = dyntofc.DynmatToForceConstants(primitive, supercell)
        
        converter.set_dynamical_matrices(dynmat = [fcm])

        converter.run()
        fc = converter.get_force_constants()
        fc = np.reshape(np.swapaxes(fc,1,2),(len(fc)*3,len(fc)*3))
        fc = get_symmetrized_FCM(fc, operations, sharedops)

        return ForceConstantMatrix(fc)
 
                    
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
