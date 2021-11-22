"""
Piezo sensitivity analysis module.
"""

import warnings

import numpy as np
from monty.dev import requires

import pymatgen.io.phonopy
from pymatgen.core.tensors import Tensor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga

try:
    from phonopy import Phonopy
    from phonopy.harmonic import dynmat_to_fc as dyntofc
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
        Create an BornEffectiveChargeTensor object defined by a
        structure, point operations of the structure's atomic sites.
        Note that the constructor uses __new__ rather than __init__
        according to the standard method ofsubclassing numpy ndarrays.

        Args:
            input_matrix (Nx3x3 array-like): the Nx3x3 array-like
                representing the born effective charge tensor
        """

        self.structure = structure
        self.bec = bec
        self.pointops = pointops
        self.BEC_operations = None
        if np.sum(self.bec) >= tol:
            warnings.warn("Input born effective charge tensor does not satisfy charge neutrality")

    def get_BEC_operations(self, eigtol=1e-05, opstol=1e-03):
        """
        Returns the symmetry operations which maps the tensors
        belonging to equivalent sites onto each other in the form
        [site index 1, site index 2, [Symmops mapping from site
        index 1 to site index 2]]


        Args:
            eigtol (float): tolerance for determining if two sites are
            related by symmetry
            opstol (float): tolerance for determining if a symmetry
            operation relates two sites

        Return:
            list of symmetry operations mapping equivalent sites and
            the indexes of those sites.
        """
        bec = self.bec
        struc = self.structure
        ops = sga(struc).get_symmetry_operations(cartesian=True)
        uniquepointops = []
        for op in ops:
            uniquepointops.append(op)

        for ops in self.pointops:
            for op in ops:
                if op not in uniquepointops:
                    uniquepointops.append(op)

        passed = []
        relations = []
        for site, val in enumerate(bec):
            unique = 1
            eig1, vecs1 = np.linalg.eig(val)
            index = np.argsort(eig1)
            neweig = np.real([eig1[index[0]], eig1[index[1]], eig1[index[2]]])
            for index, p in enumerate(passed):
                if np.allclose(neweig, p[1], atol=eigtol):
                    relations.append([site, index])
                    unique = 0
                    passed.append([site, p[0], neweig])
                    break
            if unique == 1:
                relations.append([site, site])
                passed.append([site, neweig])
        BEC_operations = []
        for atom, r in enumerate(relations):
            BEC_operations.append(r)
            BEC_operations[atom].append([])

            for op in uniquepointops:
                new = op.transform_tensor(self.bec[relations[atom][1]])

                # Check the matrix it references
                if np.allclose(new, self.bec[r[0]], atol=opstol):
                    BEC_operations[atom][2].append(op)

        self.BEC_operations = BEC_operations

    def get_rand_BEC(self, max_charge=1):
        """
        Generate a random born effective charge tensor which obeys a structure's
        symmetry and the acoustic sum rule

        Args:
            max_charge (float): maximum born effective charge value

        Return:
            np.array Born effective charge tensor
        """

        struc = self.structure
        symstruc = sga(struc)
        symstruc = symstruc.get_symmetrized_structure()

        l = len(struc)
        BEC = np.zeros((l, 3, 3))
        for atom, ops in enumerate(self.BEC_operations):
            if ops[0] == ops[1]:
                temp_tensor = Tensor(np.random.rand(3, 3) - 0.5)
                temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in self.pointops[atom]]) / len(
                    self.pointops[atom]
                )
                BEC[atom] = temp_tensor
            else:
                tempfcm = np.zeros([3, 3])
                for op in ops[2]:

                    tempfcm += op.transform_tensor(BEC[self.BEC_operations[atom][1]])
                BEC[ops[0]] = tempfcm
                if len(ops[2]) != 0:
                    BEC[ops[0]] = BEC[ops[0]] / len(ops[2])

        #     Enforce Acoustic Sum
        disp_charge = np.einsum("ijk->jk", BEC) / l
        add = np.zeros([l, 3, 3])

        for atom, ops in enumerate(self.BEC_operations):

            if ops[0] == ops[1]:
                temp_tensor = Tensor(disp_charge)
                temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in self.pointops[atom]]) / len(
                    self.pointops[atom]
                )
                add[ops[0]] = temp_tensor
            else:
                temp_tensor = np.zeros([3, 3])
                for op in ops[2]:

                    temp_tensor += op.transform_tensor(add[self.BEC_operations[atom][1]])

                add[ops[0]] = temp_tensor

                if len(ops) != 0:
                    add[ops[0]] = add[ops[0]] / len(ops[2])

        BEC = BEC - add

        return BEC * max_charge


class InternalStrainTensor:
    """
    This class describes the Nx3x3x3 internal tensor defined by a
    structure, point operations of the structure's atomic sites.
    """

    def __init__(self, structure, ist, pointops, tol=1e-3):
        """
        Create an InternalStrainTensor object.

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
            warnings.warn("Input internal strain tensor does not satisfy standard symmetries")

    def get_IST_operations(self, opstol=1e-03):
        """
        Returns the symmetry operations which maps the tensors
        belonging to equivalent sites onto each other in the form
        [site index 1, site index 2, [Symmops mapping from site
        index 1 to site index 2]]


        Args:
            opstol (float): tolerance for determining if a symmetry
            operation relates two sites

        Return:
            list of symmetry operations mapping equivalent sites and
            the indexes of those sites.
        """

        struc = self.structure
        ops = sga(struc).get_symmetry_operations(cartesian=True)
        uniquepointops = []
        for op in ops:
            uniquepointops.append(op)

        for ops in self.pointops:
            for op in ops:
                if op not in uniquepointops:
                    uniquepointops.append(op)

        IST_operations = []
        for atom in range(len(self.ist)):  # pylint: disable=C0200
            IST_operations.append([])
            for j in range(0, atom):
                for op in uniquepointops:
                    new = op.transform_tensor(self.ist[j])

                    # Check the matrix it references
                    if np.allclose(new, self.ist[atom], atol=opstol):
                        IST_operations[atom].append([j, op])

        self.IST_operations = IST_operations

    def get_rand_IST(self, max_force=1):
        """
        Generate a random internal strain tensor which obeys a structure's
        symmetry and the acoustic sum rule

        Args:
            max_force(float): maximum born effective charge value

        Return:
            InternalStrainTensor object
        """

        l = len(self.structure)
        IST = np.zeros((l, 3, 3, 3))
        for atom, ops in enumerate(self.IST_operations):
            temp_tensor = np.zeros([3, 3, 3])
            for op in ops:
                temp_tensor += op[1].transform_tensor(IST[op[0]])

            if len(ops) == 0:
                temp_tensor = Tensor(np.random.rand(3, 3, 3) - 0.5)
                for dim in range(3):
                    temp_tensor[dim] = (temp_tensor[dim] + temp_tensor[dim].T) / 2
                temp_tensor = sum([temp_tensor.transform(symm_op) for symm_op in self.pointops[atom]]) / len(
                    self.pointops[atom]
                )
            IST[atom] = temp_tensor
            if len(ops) != 0:
                IST[atom] = IST[atom] / len(ops)

        return IST * max_force


class ForceConstantMatrix:
    """
    This class describes the NxNx3x3 force constant matrix defined by a
    structure, point operations of the structure's atomic sites, and the
    shared symmetry operations between pairs of atomic sites.
    """

    def __init__(self, structure, fcm, pointops, sharedops, tol=1e-3):
        """
        Create an ForceConstantMatrix object.

        Args:
            input_matrix (NxNx3x3 array-like): the NxNx3x3 array-like
                representing the force constant matrix
        """

        self.structure = structure
        self.fcm = fcm
        self.pointops = pointops
        self.sharedops = sharedops
        self.FCM_operations = None

    def get_FCM_operations(self, eigtol=1e-05, opstol=1e-05):
        """
        Returns the symmetry operations which maps the tensors
        belonging to equivalent sites onto each other in the form
        [site index 1a, site index 1b, site index 2a, site index 2b,
        [Symmops mapping from site index 1a, 1b to site index 2a, 2b]]


        Args:
            eigtol (float): tolerance for determining if two sites are
            related by symmetry
            opstol (float): tolerance for determining if a symmetry
            operation relates two sites

        Return:
            list of symmetry operations mapping equivalent sites and
            the indexes of those sites.
        """
        struc = self.structure
        ops = sga(struc).get_symmetry_operations(cartesian=True)
        uniquepointops = []
        for op in ops:
            uniquepointops.append(op)

        for ops in self.pointops:
            for op in ops:
                if op not in uniquepointops:
                    uniquepointops.append(op)

        passed = []
        relations = []
        for atom1 in range(len(self.fcm)):  # pylint: disable=C0200
            for atom2 in range(atom1, len(self.fcm)):
                unique = 1
                eig1, vecs1 = np.linalg.eig(self.fcm[atom1][atom2])
                index = np.argsort(eig1)
                neweig = np.real([eig1[index[0]], eig1[index[1]], eig1[index[2]]])

                for entry, p in enumerate(passed):
                    if np.allclose(neweig, p[2], atol=eigtol):
                        relations.append([atom1, atom2, p[0], p[1]])
                        unique = 0
                        break
                if unique == 1:
                    relations.append([atom1, atom2, atom2, atom1])
                    passed.append([atom1, atom2, np.real(neweig)])
        FCM_operations = []
        for entry, r in enumerate(relations):
            FCM_operations.append(r)
            FCM_operations[entry].append([])

            good = 0
            for op in uniquepointops:
                new = op.transform_tensor(self.fcm[r[2]][r[3]])

                if np.allclose(new, self.fcm[r[0]][r[1]], atol=opstol):
                    FCM_operations[entry][4].append(op)
                    good = 1
            if r[0] == r[3] and r[1] == r[2]:
                good = 1
            if r[0] == r[2] and r[1] == r[3]:
                good = 1
            if good == 0:
                FCM_operations[entry] = [
                    r[0],
                    r[1],
                    r[3],
                    r[2],
                ]
                FCM_operations[entry].append([])
                for op in uniquepointops:
                    new = op.transform_tensor(self.fcm[r[2]][r[3]])
                    if np.allclose(
                        new.T,
                        self.fcm[r[0]][r[1]],
                        atol=opstol,
                    ):
                        FCM_operations[entry][4].append(op)

        self.FCM_operations = FCM_operations
        return FCM_operations

    def get_unstable_FCM(self, max_force=1):
        """
        Generate an unsymmeterized force constant matrix

        Args:
            max_charge (float): maximum born effective charge value

        Return:
            numpy array representing the force constant matrix
        """

        struc = self.structure
        operations = self.FCM_operations
        # set max force in reciprocal space
        numsites = len(struc.sites)
        D = (1 / max_force) * 2 * (np.ones([numsites * 3, numsites * 3]))
        for op in operations:
            same = 0
            transpose = 0
            if op[0] == op[1] and op[0] == op[2] and op[0] == op[3]:
                same = 1
            if op[0] == op[3] and op[1] == op[2]:
                transpose = 1
            if transpose == 0 and same == 0:
                D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] = np.zeros([3, 3])
                D[3 * op[1] : 3 * op[1] + 3, 3 * op[0] : 3 * op[0] + 3] = np.zeros([3, 3])

                for symop in op[4]:

                    tempfcm = D[3 * op[2] : 3 * op[2] + 3, 3 * op[3] : 3 * op[3] + 3]
                    tempfcm = symop.transform_tensor(tempfcm)
                    D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] += tempfcm

                if len(op[4]) != 0:
                    D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] = D[
                        3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3
                    ] / len(op[4])

                D[3 * op[1] : 3 * op[1] + 3, 3 * op[0] : 3 * op[0] + 3] = D[
                    3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3
                ].T
                continue

            temp_tensor = Tensor(np.random.rand(3, 3) - 0.5) * max_force

            temp_tensor_sum = sum([temp_tensor.transform(symm_op) for symm_op in self.sharedops[op[0]][op[1]]])
            temp_tensor_sum = temp_tensor_sum / (len(self.sharedops[op[0]][op[1]]))
            if op[0] != op[1]:
                for pair in range(len(op[4])):

                    temp_tensor2 = temp_tensor_sum.T
                    temp_tensor2 = op[4][pair].transform_tensor(temp_tensor2)
                    temp_tensor_sum = (temp_tensor_sum + temp_tensor2) / 2

            else:
                temp_tensor_sum = (temp_tensor_sum + temp_tensor_sum.T) / 2

            D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] = temp_tensor_sum
            D[3 * op[1] : 3 * op[1] + 3, 3 * op[0] : 3 * op[0] + 3] = temp_tensor_sum.T

        return D

    def get_symmetrized_FCM(self, unsymmetrized_fcm, max_force=1):
        """
        Generate a symmeterized force constant matrix from an unsymmeterized matrix

        Args:
            unsymmetrized_fcm (numpy array): unsymmeterized force constant matrix
            max_charge (float): maximum born effective charge value

        Return:
            3Nx3N numpy array representing the force constant matrix
        """

        operations = self.FCM_operations
        D = unsymmetrized_fcm
        for op in operations:
            same = 0
            transpose = 0
            if op[0] == op[1] and op[0] == operations[2] and op[0] == op[3]:
                same = 1
            if op[0] == op[3] and op[1] == op[2]:
                transpose = 1
            if transpose == 0 and same == 0:
                D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] = np.zeros([3, 3])

                for symop in op[4]:

                    tempfcm = D[3 * op[2] : 3 * op[2] + 3, 3 * op[3] : 3 * op[3] + 3]
                    tempfcm = symop.transform_tensor(tempfcm)

                    D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] += tempfcm

                if len(op[4]) != 0:
                    D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] = D[
                        3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3
                    ] / len(op[4])
                D[3 * op[1] : 3 * op[1] + 3, 3 * op[0] : 3 * op[0] + 3] = D[
                    3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3
                ].T
                continue

            temp_tensor = Tensor(D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3])
            temp_tensor_sum = sum([temp_tensor.transform(symm_op) for symm_op in self.sharedops[op[0]][op[1]]])
            if len(self.sharedops[op[0]][op[1]]) != 0:
                temp_tensor_sum = temp_tensor_sum / (len(self.sharedops[op[0]][op[1]]))

            # Apply the proper transformation if there is an equivalent already
            if op[0] != op[1]:

                for pair in range(len(op[4])):

                    temp_tensor2 = temp_tensor_sum.T
                    temp_tensor2 = op[4][pair].transform_tensor(temp_tensor2)
                    temp_tensor_sum = (temp_tensor_sum + temp_tensor2) / 2

            else:
                temp_tensor_sum = (temp_tensor_sum + temp_tensor_sum.T) / 2

            D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] = temp_tensor_sum
            D[3 * op[1] : 3 * op[1] + 3, 3 * op[0] : 3 * op[0] + 3] = temp_tensor_sum.T

        return D

    def get_stable_FCM(self, fcm, fcmasum=10):
        """
        Generate a symmeterized force constant matrix that obeys the objects symmetry
        constraints, has no unstable modes and also obeys the acoustic sum rule through an
        iterative procedure

        Args:
            fcm (numpy array): unsymmeterized force constant matrix
            fcmasum (int): number of iterations to attempt to obey the acoustic sum
                rule

        Return:
            3Nx3N numpy array representing the force constant matrix
        """

        check = 0
        count = 0
        while check == 0:
            # if resymmetrizing brings back unstable modes 20 times, the method breaks
            if count > 20:
                check = 1
                break

            eigs, vecs = np.linalg.eig(fcm)

            maxeig = np.max(-1 * eigs)
            eigsort = np.argsort(np.abs(eigs))
            for i in range(3, len(eigs)):
                if eigs[eigsort[i]] > 1e-06:
                    eigs[eigsort[i]] = -1 * maxeig * np.random.rand()
            diag = np.real(np.eye(len(fcm)) * eigs)

            fcm = np.real(np.matmul(np.matmul(vecs, diag), vecs.T))
            fcm = self.get_symmetrized_FCM(fcm)
            fcm = self.get_asum_FCM(fcm)
            eigs, vecs = np.linalg.eig(fcm)
            unstable_modes = 0
            eigsort = np.argsort(np.abs(eigs))
            for i in range(3, len(eigs)):
                if eigs[eigsort[i]] > 1e-06:
                    unstable_modes = 1
            if unstable_modes == 1:
                count = count + 1
                continue
            check = 1

        return fcm

    # acoustic sum

    def get_asum_FCM(self, fcm, numiter=15):
        """
        Generate a symmeterized force constant matrix that obeys the objects symmetry
        constraints and obeys the acoustic sum rule through an iterative procedure

        Args:
            fcm (numpy array): 3Nx3N unsymmeterized force constant matrix
            numiter (int): number of iterations to attempt to obey the acoustic sum
                rule

        Return:
            numpy array representing the force constant matrix
        """

        # set max force in reciprocal space
        operations = self.FCM_operations
        numsites = len(self.structure)

        D = np.ones([numsites * 3, numsites * 3])
        for num in range(numiter):
            X = np.real(fcm)

            # symmetry operations
            pastrow = 0
            total = np.zeros([3, 3])
            for col in range(numsites):
                total = total + X[0:3, col * 3 : col * 3 + 3]

            total = total / (numsites)
            for op in operations:
                same = 0
                transpose = 0
                if op[0] == op[1] and op[0] == op[2] and op[0] == op[3]:
                    same = 1
                if op[0] == op[3] and op[1] == op[2]:
                    transpose = 1
                if transpose == 0 and same == 0:
                    D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] = np.zeros([3, 3])

                    for symop in op[4]:

                        tempfcm = D[3 * op[2] : 3 * op[2] + 3, 3 * op[3] : 3 * op[3] + 3]
                        tempfcm = symop.transform_tensor(tempfcm)

                        D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] += tempfcm

                    if len(op[4]) != 0:
                        D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] = D[
                            3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3
                        ] / len(op[4])
                    D[3 * op[1] : 3 * op[1] + 3, 3 * op[0] : 3 * op[0] + 3] = D[
                        3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3
                    ].T
                    continue
                # Get the difference in the sum up to this point
                currrow = op[0]
                if currrow != pastrow:
                    total = np.zeros([3, 3])
                    for col in range(numsites):
                        total = total + X[currrow * 3 : currrow * 3 + 3, col * 3 : col * 3 + 3]
                    for col in range(currrow):
                        total = total - D[currrow * 3 : currrow * 3 + 3, col * 3 : col * 3 + 3]
                    total = total / (numsites - currrow)
                pastrow = currrow

                # Apply the point symmetry operations of the site
                temp_tensor = Tensor(total)
                temp_tensor_sum = sum([temp_tensor.transform(symm_op) for symm_op in self.sharedops[op[0]][op[1]]])

                if len(self.sharedops[op[0]][op[1]]) != 0:
                    temp_tensor_sum = temp_tensor_sum / (len(self.sharedops[op[0]][op[1]]))

                # Apply the proper transformation if there is an equivalent already
                if op[0] != op[1]:

                    for pair in range(len(op[4])):

                        temp_tensor2 = temp_tensor_sum.T
                        temp_tensor2 = op[4][pair].transform_tensor(temp_tensor2)
                        temp_tensor_sum = (temp_tensor_sum + temp_tensor2) / 2

                else:
                    temp_tensor_sum = (temp_tensor_sum + temp_tensor_sum.T) / 2

                D[3 * op[0] : 3 * op[0] + 3, 3 * op[1] : 3 * op[1] + 3] = temp_tensor_sum
                D[3 * op[1] : 3 * op[1] + 3, 3 * op[0] : 3 * op[0] + 3] = temp_tensor_sum.T
            fcm = fcm - D

        return fcm

    @requires(Phonopy, "phonopy not installed!")
    def get_rand_FCM(self, asum=15, force=10):
        """
        Generate a symmeterized force constant matrix from an unsymmeterized matrix
        that has no unstable modes and also obeys the acoustic sum rule through an
        iterative procedure

        Args:
            force (float): maximum force constant
            asum (int): number of iterations to attempt to obey the acoustic sum
                rule

        Return:
            NxNx3x3 np.array representing the force constant matrix

        """

        numsites = len(self.structure.sites)
        structure = pymatgen.io.phonopy.get_phonopy_structure(self.structure)
        pnstruc = Phonopy(structure, np.eye(3), np.eye(3))

        dyn = self.get_unstable_FCM(force)
        dyn = self.get_stable_FCM(dyn)

        dyn = np.reshape(dyn, (numsites, 3, numsites, 3)).swapaxes(1, 2)

        dynmass = np.zeros([len(self.structure), len(self.structure), 3, 3])
        masses = []
        for j in range(numsites):
            masses.append(self.structure.sites[j].specie.atomic_mass)
        dynmass = np.zeros([numsites, numsites, 3, 3])
        for m in range(numsites):
            for n in range(numsites):
                dynmass[m][n] = dyn[m][n] * np.sqrt(masses[m]) * np.sqrt(masses[n])

        supercell = pnstruc.get_supercell()
        primitive = pnstruc.get_primitive()

        converter = dyntofc.DynmatToForceConstants(primitive, supercell)

        dyn = np.reshape(np.swapaxes(dynmass, 1, 2), (numsites * 3, numsites * 3))

        converter.set_dynamical_matrices(dynmat=[dyn])

        converter.run()
        fc = converter.get_force_constants()

        return fc


def get_piezo(BEC, IST, FCM, rcond=0.0001):
    """
    Generate a random piezoelectric tensor based on a structure and corresponding
    symmetry

    Args:
        BEC (numpy array): Nx3x3 array representing the born effective charge tensor
        IST (numpy array): Nx3x3x3 array representing the internal strain tensor
        FCM (numpy array): NxNx3x3 array representing the born effective charge tensor
        rcondy (float): condition for excluding eigenvalues in the pseudoinverse

    Return:
        3x3x3 calculated Piezo tensor
    """

    numsites = len(BEC)
    temp_fcm = np.reshape(np.swapaxes(FCM, 1, 2), (numsites * 3, numsites * 3))

    eigs, vecs = np.linalg.eig(temp_fcm)
    K = np.linalg.pinv(
        -temp_fcm,
        rcond=np.abs(eigs[np.argsort(np.abs(eigs))[2]]) / np.abs(eigs[np.argsort(np.abs(eigs))[-1]]) + rcond,
    )

    K = np.reshape(K, (numsites, 3, numsites, 3)).swapaxes(1, 2)
    return np.einsum("ikl,ijlm,jmno->kno", BEC, K, IST) * 16.0216559424


@requires(Phonopy, "phonopy not installed!")
def rand_piezo(struc, pointops, sharedops, BEC, IST, FCM, anumiter=10):
    """
    Generate a random piezoelectric tensor based on a structure and corresponding
    symmetry

    Args:
        struc (pymatgen structure): structure whose symmetry operations the piezo tensor must obey
        pointops: list of point operations obeyed by a single atomic site
        sharedops: list of point operations shared by a pair of atomic sites
        BEC (numpy array): Nx3x3 array representing the born effective charge tensor
        IST (numpy array): Nx3x3x3 array representing the internal strain tensor
        FCM (numpy array): NxNx3x3 array representing the born effective charge tensor
        anumiter (int): number of iterations for acoustic sum rule convergence
    Return:
        list in the form of [Nx3x3 random born effective charge tenosr,
        Nx3x3x3 random internal strain tensor, NxNx3x3 random force constant matrix, 3x3x3 piezo tensor]
    """
    bec = BornEffectiveCharge(struc, BEC, pointops)
    bec.get_BEC_operations()
    rand_BEC = bec.get_rand_BEC()

    ist = InternalStrainTensor(struc, IST, pointops)
    ist.get_IST_operations()
    rand_IST = ist.get_rand_IST()

    fcm = ForceConstantMatrix(struc, FCM, pointops, sharedops)
    fcm.get_FCM_operations()
    rand_FCM = fcm.get_rand_FCM()

    P = get_piezo(rand_BEC, rand_IST, rand_FCM) * 16.0216559424 / struc.volume

    return (rand_BEC, rand_IST, rand_FCM, P)
