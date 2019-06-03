from pymatgen.core.tensors import Tensor
from pymatgen.symmetry import site_symmetries as ss
from pymatgen.analysis.piezo import PiezoTensor
import pymatgen.io.phonopy

import numpy as np
import warnings
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
from phonopy.harmonic import dynmat_to_fc as dyntofc

try:
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms
except ImportError:
    Phonopy = None

__author__ = "Handong Ling"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Handong Ling"
__email__ = "hling@lbl.gov"
__status__ = "Development"
__date__ = "Feb, 2019"

class BornEffectiveCharge:
    """
    This class describes the Nx3x3 born effective charge tensor
    """

    def __init__(self, structure, bec, pointops, tol=1e-3):
        """
        Create an BornEffectiveChargeTensor object.  The constructor throws 
        an error if the input_matrix argument does not satisfy the charge 
        sum rule. Note that the constructor uses __new__ rather than __init__ 
        according to the standard method of subclassing numpy ndarrays.

        Args:
            input_matrix (Nx3x3 array-like): the Nx3x3 array-like
                representing the born effective charge tensor
        """
        self.structure = structure
        self.bec = bec
        self.pointops = pointops
        self.BEC_operations = None
        if not np.sum(self.bec) < tol:
            warnings.warn("Input born effective charge tensor does "
                          "not satisfy charge neutrality")

    def get_BEC_IST_operations(self, eigtol = 1e-05, opstol = 1e-03):
        bec = self.bec
        struc = self.structure
        ops = sga(struc).get_symmetry_operations(cartesian = True)
        uniquepointops = [] 
        for i in range(len(ops)):
            uniquepointops.append(ops[i])
        
        for i in range(len(self.pointops)):
            for j in range(len(self.pointops[i])):
                if self.pointops[i][j] not in uniquepointops:
                    uniquepointops.append(self.pointops[i][j])
         
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
                new = uniquepointops[j].transform_tensor(self.bec[relations[i][1]])
                
                ## Check the matrix it references
                if np.allclose(new, self.bec[relations[i][0]], atol = opstol):
                    BEC_IST_operations[i][2].append(uniquepointops[j])

        self.BEC_operations = BEC_IST_operations

                        

    def get_rand_BEC(self, max_charge=1):
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
        struc = self.structure
        symstruc = sga(struc)
        symstruc = symstruc.get_symmetrized_structure()
        
        
        l = len(struc)
        BEC = np.zeros((l,3,3))
        primsites = []
        for i in range(len(self.BEC_operations)):
            if self.BEC_operations[i][0] == self.BEC_operations[i][1]:
                temp_tensor = Tensor(np.random.rand(3,3)-0.5)         
                temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in self.pointops[i]]) \
                    /len(self.pointops[i])
                BEC[i] = temp_tensor
            else: 
                tempfcm = np.zeros([3,3])
                for j in range(len(self.BEC_operations[i][2])):

                    tempfcm += self.BEC_operations[i][2][j].transform_tensor(BEC[self.BEC_operations[i][1]])
                BEC[self.BEC_operations[i][0]] = tempfcm
                if len(self.BEC_operations[i][2]) != 0:
                    BEC[self.BEC_operations[i][0]] = BEC[self.BEC_operations[i][0]]/len(self.BEC_operations[i][2])
        
    #     Enforce Acoustic Sum
        disp_charge = np.einsum("ijk->jk",BEC)/l
        add = np.zeros([l,3,3])
        
        for i in range(len(self.BEC_operations)):  
            
            if self.BEC_operations[i][0] == self.BEC_operations[i][1]:
                temp_tensor = Tensor(disp_charge)
                temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in self.pointops[i]]) \
                    /len(self.pointops[i])
                add[self.BEC_operations[i][0]] = temp_tensor
            else: 
                temp_tensor = np.zeros([3,3])
                for j in range(len(self.BEC_operations[i][2])):

                    temp_tensor += self.BEC_operations[i][2][j].transform_tensor(add[self.BEC_operations[i][1]])
                
                add[self.BEC_operations[i][0]] = temp_tensor
                
                if len(self.BEC_operations[i][2]) != 0:
                    add[self.BEC_operations[i][0]] = add[self.BEC_operations[i][0]]/len(self.BEC_operations[i][2])
            
        BEC = BEC - add
         

        return BEC*max_charge




class InternalStrainTensor:
    """
    This class describes the Nx3x3x3 born effective charge tensor
    """

    def __init__(self, structure, ist, pointops,  tol=1e-3):
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

        self.structure = structure
        self.ist = ist
        self.pointops = pointops
        self.IST_operations = None 

        obj = self.ist
        if not (obj - np.transpose(obj, (0, 1, 3, 2)) < tol).all():
            warnings.warn("Input internal strain tensor does "
                          "not satisfy standard symmetries")

    def get_BEC_IST_operations(self, bec, eigtol = 1e-05, opstol = 1e-03):
        struc = self.structure
        ops = sga(struc).get_symmetry_operations(cartesian = True)
        uniquepointops = [] 
        for i in range(len(ops)):
            uniquepointops.append(ops[i])
        
        for i in range(len(self.pointops)):
            for j in range(len(self.pointops[i])):
                if self.pointops[i][j] not in uniquepointops:
                    uniquepointops.append(self.pointops[i][j])
                
         
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

        self.IST_operations = BEC_IST_operations    

    def get_rand_IST(self, max_force=1):
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
        struc = self.structure
        symstruc = sga(struc)
        symstruc = symstruc.get_symmetrized_structure()
        
        
        l = len(struc)
        IST = np.zeros((l,3,3,3))
        primsites = []
        for i in range(len(self.IST_operations)):
            if self.IST_operations[i][0] == self.IST_operations[i][1]:
                temp_tensor = Tensor(np.random.rand(3,3,3)-0.5) 
                for dim in range(3):
                    temp_tensor[dim] = (temp_tensor[dim]+temp_tensor[dim].T)/2
                temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in self.pointops[i]]) \
                    /len(self.pointops[i])
                IST[self.IST_operations[i][0]] = temp_tensor
            else: 
                temp_tensor = np.zeros([3,3,3])
                for j in range(len(self.IST_operations[i][2])):

                    temp_tensor += self.IST_operations[i][2][j].transform_tensor(IST[self.IST_operations[i][1]])
                
                IST[self.IST_operations[i][0]] = temp_tensor
                if len(self.IST_operations[i][2]) != 0:
                    IST[self.IST_operations[i][0]] = IST[self.IST_operations[i][0]]/len(self.IST_operations[i][2])
        
        return IST*max_force



class ForceConstantMatrix:
    """
    This class describes the NxNx3x3 born effective charge tensor
    """

    def __init__(self, structure, fcm, pointops, sharedops, tol=1e-3):
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

        self.structure = structure
        self.fcm = fcm 
        self.pointops = pointops
        self.sharedops = sharedops
        self.FCM_operations = None

    def get_FCM_operations(self, eigtol = 1e-05, opstol = 1e-05):
        struc = self.structure
        ops = sga(struc).get_symmetry_operations(cartesian = True)
        uniquepointops = [] 
        for i in range(len(ops)):
            uniquepointops.append(ops[i])
        
        for i in range(len(self.pointops)):
            for j in range(len(self.pointops[i])):
                if self.pointops[i][j] not in uniquepointops:
                    uniquepointops.append(self.pointops[i][j])
                
        passed = []
        relations = []
        for i in range(len(self.fcm)):
            for j in range(i, len(self.fcm)):
                unique = 1
                eig1, vecs1 = np.linalg.eig(self.fcm[i][j])
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
                new = uniquepointops[j].transform_tensor(self.fcm[relations[i][2]][relations[i][3]])
                
                if np.allclose(new, self.fcm[relations[i][0]][relations[i][1]], atol = opstol):
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
                    new = uniquepointops[j].transform_tensor(self.fcm[relations[i][2]][relations[i][3]])
                    if np.allclose(new.T, self.fcm[relations[i][0]][relations[i][1]], atol = opstol):
                        FCM_operations[i][4].append(uniquepointops[j])
         
        self.FCM_operations = FCM_operations    
        return FCM_operations

                    

    def get_unstable_FCM(self, max_force = 1):

        struc = self.structure
        operations = self.FCM_operations
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
                                       self.sharedops[operations[i][0]][operations[i][1]]])
                temp_tensor_sum = temp_tensor_sum/(len(self.sharedops[operations[i][0]][operations[i][1]]))
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
    def get_symmetrized_FCM(self, unsymmetrized_fcm, max_force = 1):
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

        operations = self.FCM_operations
        numsites = len(self.structure)
        D = unsymmetrized_fcm
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
    #         print(5, operations[i][0], operations[i][1])

                temp_tensor = Tensor(D[3*operations[i][0]:3*operations[i][0]+3, 3*operations[i][1]:3*operations[i][1]+3])
                temp_tensor_sum = sum([temp_tensor.transform(symm_op) for symm_op in \
                                       self.sharedops[operations[i][0]][operations[i][1]]])
                if len(self.sharedops[operations[i][0]][operations[i][1]]) != 0:
                    temp_tensor_sum = temp_tensor_sum/(len(self.sharedops[operations[i][0]][operations[i][1]]))



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

        return(D)


    def get_stable_FCM(self, fcm, fcmasum = 10):
        
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

        numsites = len(self.structure)
        check = 0
        count = 0
        # fcm = self.get_asum_FCM(fcm)
        while check == 0:
            if count > 100:
                check = 1
                break
            
            eigs, vecs = np.linalg.eig(fcm)

            maxeig = np.max(-1*eigs)
            eigsort = np.argsort(np.abs(eigs))
            for i in range(3,len(eigs)):
                if eigs[eigsort[i]] > 1e-06:
                    eigs[eigsort[i]] = -1*maxeig*np.random.rand()
            diag = np.real(np.eye(len(fcm))*eigs)

            fcm = np.real(np.matmul(np.matmul(vecs, diag), vecs.T))
            fcm = self.get_symmetrized_FCM(fcm)
            fcm = self.get_asum_FCM(fcm)
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

    def get_asum_FCM(self, fcm, numiter = 15):
    # set max force in reciprocal space 
        operations = self.FCM_operations
        numsites = len(self.structure)
        
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
                if operations[i][0] == operations[i][1] and operations[i][0] == operations[i][2] and operations[i][0] == operations[i][3]:
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
                                           self.sharedops[operations[i][0]][operations[i][1]]])

                    if len(self.sharedops[operations[i][0]][operations[i][1]]) != 0:
                        temp_tensor_sum = temp_tensor_sum/(len(self.sharedops[operations[i][0]][operations[i][1]]))



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

    def get_rand_FCM(self, asum = 15, force = 10):    
        """
        Generate a symmeterized force constant matrix from an unsymmeterized matrix
        that has no unstable modes and also obeys the acoustic sum rule through an 
        iterative procedure

        Args:
            fcm (numpy array): unsymmeterized force constant matrix 
            operations: list of operation mappings for indexed sites in the force
                constant matrix
            sharedops: list of point operations shared by a pair of atomic sites
            force (float): maximum force constant
            asum (int): number of iterations to attempt to obey the acoustic sum
                rule

        Return:
            force constant matrix representing the force constant matrix

        """

        numsites = len(self.structure.sites)
        structure = pymatgen.io.phonopy.get_phonopy_structure(self.structure)
        pnstruc = Phonopy(structure, np.eye(3), np.eye(3))
        
        dyn = self.get_unstable_FCM(force)
        # dyn = self.get_asum_FCM(dyn, self.FCM_operations, self.sharedops, numiter = asum)
        dyn = self.get_stable_FCM(dyn)
        
        dyn = np.reshape(dyn, (numsites,3,numsites,3)).swapaxes(1,2)
        

        dynmass= np.zeros([len(self.structure),len(self.structure),3,3])
        masses = []
        for j in range(numsites):
            masses.append(self.structure.sites[j].specie.atomic_mass)
        dynmass = np.zeros([numsites,numsites,3,3])
        for m in range(numsites):
            for n in range(numsites):
                dynmass[m][n] = dyn[m][n]*np.sqrt(masses[m])*np.sqrt(masses[n])

        supercell = pnstruc.get_supercell()
        primitive = pnstruc.get_primitive()
    
        converter = dyntofc.DynmatToForceConstants(primitive, supercell)

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

        return fc
 
def get_piezo(BEC, IST, FCM, rcondy = 0.0001):
    numsites = len(BEC)
    temp_fcm = np.reshape(np.swapaxes(FCM,1,2),(numsites*3,numsites*3))
    
    piezo = np.zeros([3,3,3])
    eigs,vecs = np.linalg.eig(temp_fcm) 
    K = np.linalg.pinv(-temp_fcm, \
        rcond = np.abs(eigs[np.argsort(np.abs(eigs))[2]])/np.abs(eigs[np.argsort(np.abs(eigs))[-1]])+rcondy)
    
    K = np.reshape(K, (numsites,3,numsites,3)).swapaxes(1,2)
    return np.einsum("ikl,ijlm,jmno->kno",BEC,K,IST)*16.0216559424

def rand_piezo(struc, pointops, sharedops, BEC, IST, FCM, anumiter = 10):
    """
    Generate a random piezoelectric tensor based on a structure and corresponding
    symmetry

    Args:
        struc (pymatgen structure): structure whose symmetry operations the piezo tensor must obey
        pointops: list of point operations obeyed by a single atomic site
        sharedops: list of point operations shared by a pair of atomic sites
        BEC_IST_operations: list of point operations a site must obey or the transformation
            that brings one site's properties to anothers
        FCM_operations: list of point operations a pair of sites must obey or the transformation
            that brings one pair of site's properties to anothers

    Return:
        PiezoTensor object
    """
    numsites = len(struc.sites)
    bec = BornEffectiveCharge(struc, BEC, pointops)
    bec.get_BEC_IST_operations() 
    rand_BEC = bec.get_rand_BEC()

    ist = InternalStrainTensor(struc, IST, pointops)
    ist.get_BEC_IST_operations(BEC) 
    rand_IST = ist.get_rand_IST()

    fcm = ForceConstantMatrix(struc, FCM, pointops, sharedops)
    fcm.get_FCM_operations()
    rand_FCM = fcm.get_rand_FCM()

    P = get_piezo(rand_BEC, rand_IST, rand_FCM)*16.0216559424/struc.volume

    return (rand_BEC,rand_IST,rand_FCM,P)
