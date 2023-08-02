---
layout: default
title: pymatgen.analysis.piezo_sensitivity.md
nav_exclude: true
---

# pymatgen.analysis.piezo_sensitivity module

Piezo sensitivity analysis module.


### _class_ pymatgen.analysis.piezo_sensitivity.BornEffectiveCharge(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), bec, pointops, tol: float = 0.001)
Bases: `object`

This class describes the Nx3x3 born effective charge tensor.

Create an BornEffectiveChargeTensor object defined by a
structure, point operations of the structure’s atomic sites.
Note that the constructor uses __new__ rather than __init__
according to the standard method of subclassing numpy ndarrays.


* **Parameters**

    **input_matrix** (*Nx3x3 array-like*) – the Nx3x3 array-like
    representing the born effective charge tensor



#### get_BEC_operations(eigtol=1e-05, opstol=0.001)
Returns the symmetry operations which maps the tensors
belonging to equivalent sites onto each other in the form
[site index 1, site index 2, [Symmops mapping from site
index 1 to site index 2]].


* **Parameters**


    * **eigtol** (*float*) – tolerance for determining if two sites are


    * **symmetry** (*related by*) –


    * **opstol** (*float*) – tolerance for determining if a symmetry


    * **sites** (*operation relates two*) –



* **Returns**

    list of symmetry operations mapping equivalent sites and
    the indexes of those sites.



#### get_rand_BEC(max_charge=1)
Generate a random born effective charge tensor which obeys a structure’s
symmetry and the acoustic sum rule.


* **Parameters**

    **max_charge** (*float*) – maximum born effective charge value



* **Returns**

    np.array Born effective charge tensor



### _class_ pymatgen.analysis.piezo_sensitivity.ForceConstantMatrix(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), fcm, pointops, sharedops, tol: float = 0.001)
Bases: `object`

This class describes the NxNx3x3 force constant matrix defined by a
structure, point operations of the structure’s atomic sites, and the
shared symmetry operations between pairs of atomic sites.

Create an ForceConstantMatrix object.


* **Parameters**

    **input_matrix** (*NxNx3x3 array-like*) – the NxNx3x3 array-like
    representing the force constant matrix



#### get_FCM_operations(eigtol=1e-05, opstol=1e-05)
Returns the symmetry operations which maps the tensors
belonging to equivalent sites onto each other in the form
[site index 1a, site index 1b, site index 2a, site index 2b,
[Symmops mapping from site index 1a, 1b to site index 2a, 2b]].


* **Parameters**


    * **eigtol** (*float*) – tolerance for determining if two sites are


    * **symmetry** (*related by*) –


    * **opstol** (*float*) – tolerance for determining if a symmetry


    * **sites** (*operation relates two*) –



* **Returns**

    list of symmetry operations mapping equivalent sites and
    the indexes of those sites.



#### get_asum_FCM(fcm: ndarray, numiter: int = 15)
Generate a symmeterized force constant matrix that obeys the objects symmetry
constraints and obeys the acoustic sum rule through an iterative procedure.


* **Parameters**


    * **fcm** (*numpy array*) – 3Nx3N unsymmeterized force constant matrix


    * **numiter** (*int*) – number of iterations to attempt to obey the acoustic sum
    rule



* **Returns**

    numpy array representing the force constant matrix



#### get_rand_FCM(asum=15, force=10)
Generate a symmeterized force constant matrix from an unsymmeterized matrix
that has no unstable modes and also obeys the acoustic sum rule through an
iterative procedure.


* **Parameters**


    * **force** (*float*) – maximum force constant


    * **asum** (*int*) – number of iterations to attempt to obey the acoustic sum
    rule



* **Returns**

    NxNx3x3 np.array representing the force constant matrix



#### get_stable_FCM(fcm, fcmasum=10)
Generate a symmeterized force constant matrix that obeys the objects symmetry
constraints, has no unstable modes and also obeys the acoustic sum rule through an
iterative procedure.


* **Parameters**


    * **fcm** (*numpy array*) – unsymmeterized force constant matrix


    * **fcmasum** (*int*) – number of iterations to attempt to obey the acoustic sum
    rule



* **Returns**

    3Nx3N numpy array representing the force constant matrix



#### get_symmetrized_FCM(unsymmetrized_fcm, max_force=1)
Generate a symmeterized force constant matrix from an unsymmeterized matrix.


* **Parameters**


    * **unsymmetrized_fcm** (*numpy array*) – unsymmeterized force constant matrix


    * **max_charge** (*float*) – maximum born effective charge value



* **Returns**

    3Nx3N numpy array representing the force constant matrix



#### get_unstable_FCM(max_force=1)
Generate an unsymmeterized force constant matrix.


* **Parameters**

    **max_charge** (*float*) – maximum born effective charge value



* **Returns**

    numpy array representing the force constant matrix



### _class_ pymatgen.analysis.piezo_sensitivity.InternalStrainTensor(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), ist, pointops, tol: float = 0.001)
Bases: `object`

This class describes the Nx3x3x3 internal tensor defined by a
structure, point operations of the structure’s atomic sites.

Create an InternalStrainTensor object.


* **Parameters**

    **input_matrix** (*Nx3x3x3 array-like*) – the Nx3x3x3 array-like
    representing the internal strain tensor



#### get_IST_operations(opstol=0.001)
Returns the symmetry operations which maps the tensors
belonging to equivalent sites onto each other in the form
[site index 1, site index 2, [Symmops mapping from site
index 1 to site index 2]].


* **Parameters**


    * **opstol** (*float*) – tolerance for determining if a symmetry


    * **sites** (*operation relates two*) –



* **Returns**

    list of symmetry operations mapping equivalent sites and
    the indexes of those sites.



#### get_rand_IST(max_force=1)
Generate a random internal strain tensor which obeys a structure’s
symmetry and the acoustic sum rule.


* **Parameters**

    **max_force** (*float*) – maximum born effective charge value



* **Returns**

    InternalStrainTensor object



### pymatgen.analysis.piezo_sensitivity.get_piezo(BEC, IST, FCM, rcond=0.0001)
Generate a random piezoelectric tensor based on a structure and corresponding
symmetry.


* **Parameters**


    * **BEC** (*numpy array*) – Nx3x3 array representing the born effective charge tensor


    * **IST** (*numpy array*) – Nx3x3x3 array representing the internal strain tensor


    * **FCM** (*numpy array*) – NxNx3x3 array representing the born effective charge tensor


    * **rcondy** (*float*) – condition for excluding eigenvalues in the pseudoinverse



* **Returns**

    3x3x3 calculated Piezo tensor



### pymatgen.analysis.piezo_sensitivity.rand_piezo(struct, pointops, sharedops, BEC, IST, FCM, anumiter=10)
Generate a random piezoelectric tensor based on a structure and corresponding
symmetry.


* **Parameters**


    * **struct** (*pymatgen structure*) – structure whose symmetry operations the piezo tensor must obey


    * **pointops** – list of point operations obeyed by a single atomic site


    * **sharedops** – list of point operations shared by a pair of atomic sites


    * **BEC** (*numpy array*) – Nx3x3 array representing the born effective charge tensor


    * **IST** (*numpy array*) – Nx3x3x3 array representing the internal strain tensor


    * **FCM** (*numpy array*) – NxNx3x3 array representing the born effective charge tensor


    * **anumiter** (*int*) – number of iterations for acoustic sum rule convergence



* **Returns**

    list in the form of [Nx3x3 random born effective charge tenosr,
    Nx3x3x3 random internal strain tensor, NxNx3x3 random force constant matrix, 3x3x3 piezo tensor]