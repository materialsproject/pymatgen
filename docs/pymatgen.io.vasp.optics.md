---
layout: default
title: pymatgen.io.vasp.optics.md
nav_exclude: true
---

# pymatgen.io.vasp.optics module

Classes for parsing and manipulating VASP optical properties calculations.


### _class_ pymatgen.io.vasp.optics.DielectricFunctionCalculator(cder_real: NDArray, cder_imag: NDArray, eigs: NDArray, kweights: NDArray, nedos: int, deltae: float, ismear: int, sigma: float, efermi: float, cshift: float, ispin: int, volume: float)
Bases: `MSONable`

Class for postprocessing VASP optical properties calculations.

This objects helps load the different parameters from the vasprun.xml file but allows users to override
them as needed.

The standard vasprun.xml from an `LOPTICS=.True.` calculation already contains
the complex frequency dependent dielectric functions.  However you have no way to decompose
the different contributions.  Since the `WAVEDER` file is also written during an optical calculation,
you can reconstruct the dielectric functions purely in Python and have full control over contribution
from different bands and k-points.

VASP’s linear optics follow these steps:


    * Calculate the imaginary part


    * Perform symmetry operations (this is not implemented here)


    * Calculate the real part

Currently, this Calculator only works for `ISYM=0` calculations since we cannot guarantee that our
externally defined symmetry operations are the same as VASP’s. This can be fixed by printing the
symmetry operators into the vasprun.xml file. If this happens in future versions of VASP,
we can dramatically speed up the calculations here by considering only the irreducible kpoints.


#### _property_ cder()
Complex CDER from WAVEDER.


#### cder_imag(_: NDArra_ )

#### cder_real(_: NDArra_ )

#### cshift(_: floa_ )

#### deltae(_: floa_ )

#### efermi(_: floa_ )

#### eigs(_: NDArra_ )

#### _classmethod_ from_directory(directory: Path | str)
Construct a DielectricFunction from a directory containing vasprun.xml and WAVEDER files.


#### _classmethod_ from_vasp_objects(vrun: [Vasprun](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun), waveder: [Waveder](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder))
Construct a DielectricFunction from Vasprun, Kpoint, and Waveder objects.


* **Parameters**


    * **vrun** – Vasprun object


    * **kpoint** – Kpoint object


    * **waveder** – Waveder object



#### get_epsilon(idir: int, jdir: int, efermi: float | None = None, nedos: int | None = None, deltae: float | None = None, ismear: int | None = None, sigma: float | None = None, cshift: float | None = None, mask: NDArray | None = None)
Compute the frequency dependent dielectric function.


* **Parameters**


    * **idir** – First direction of the dielectric tensor


    * **jdir** – Second direction of the dielectric tensor


    * **efermi** – Fermi energy


    * **nedos** – Number of points in the DOS


    * **deltae** – Energy step in the DOS


    * **ismear** – Smearing method (only has 0:gaussian, >0:Methfessel-Paxton)


    * **sigma** – Smearing width


    * **cshift** – Complex shift used for Kramer-Kronig transformation


    * **mask** – Mask for the bands/kpoint/spin index to include in the calculation



#### ismear(_: in_ )

#### ispin(_: in_ )

#### kweights(_: NDArra_ )

#### nedos(_: in_ )

#### plot_weighted_transition_data(idir: int, jdir: int, mask: NDArray | None = None, min_val: float = 0.0)
Data for plotting the weight matrix elements as a scatter plot.

Since the computation of the final spectrum (especially the smearing part)
is still fairly expensive.  This function can be used to check the values
of some portion of the spectrum (defined by the mask).
In a sense, we are lookin at the imaginary part of the dielectric function
before the smearing is applied.


* **Parameters**


    * **idir** – First direction of the dielectric tensor.


    * **jdir** – Second direction of the dielectric tensor.


    * **mask** – Mask to apply to the CDER for the bands/kpoint/spin
    index to include in the calculation


    * **min_val** – Minimum value below this value the matrix element will not be shown.



#### sigma(_: floa_ )

#### volume(_: floa_ )

### pymatgen.io.vasp.optics.delta_func(x, ismear)
Replication of VASP’s delta function.


### pymatgen.io.vasp.optics.delta_methfessel_paxton(x, n)
D_n (x) = exp -x^2 \* sum_i=0^n A_i H_2i(x)
where H is a Hermite polynomial and
A_i = (-1)^i / ( i! 4^i sqrt(pi) ).


### pymatgen.io.vasp.optics.epsilon_imag(cder: NDArray, eigs: NDArray, kweights: ArrayLike, efermi: float, nedos: int, deltae: float, ismear: int, sigma: float, idir: int, jdir: int, mask: NDArray | None = None)
Replicate the EPSILON_IMAG function of VASP.


* **Parameters**


    * **cder** – The data written to the WAVEDER (nbands, nbands, nkpoints, nspin, diri, dirj)


    * **eigs** – The eigenvalues (nbands, nkpoints, nspin)


    * **kweights** – The kpoint weights (nkpoints)


    * **efermi** – The fermi energy


    * **nedos** – The sampling of the energy values


    * **deltae** – The energy grid spacing


    * **ismear** – The smearing parameter used by the `step_func`.


    * **sigma** – The width of the smearing


    * **idir** – The first direction of the dielectric tensor


    * **jdir** – The second direction of the dielectric tensor


    * **mask** – Mask for the bands/kpoint/spin index to include in the calculation



* **Returns**

    Array of size nedos with the imaginary part of the dielectric function.



* **Return type**

    np.array



### pymatgen.io.vasp.optics.get_delta(x0: float, sigma: float, nx: int, dx: float, ismear: int = 3)
Get the smeared delta function to be added to form the spectrum.

This replaces the SLOT function from VASP. Uses finite differences instead of
evaluating the delta function since the step function is more likely to have analytic form.


* **Parameters**


    * **x0** – The center of the dielectric function.


    * **sigma** – The width of the smearing


    * **nx** – The number of grid points in the output grid.


    * **dx** – The gridspacing of the output grid.


    * **ismear** – The smearing parameter used by the `step_func`.



* **Returns**

    Array of size nx with delta function on the desired outputgrid.



* **Return type**

    np.array



### pymatgen.io.vasp.optics.get_step(x0, sigma, nx, dx, ismear)
Get the smeared step function to be added to form the spectrum.

This replaces the SLOT function from VASP.


* **Parameters**


    * **x0** – The center of the dielectric function.


    * **sigma** – The width of the smearing


    * **nx** – The number of grid points in the output grid.


    * **dx** – The gridspacing of the output grid.


    * **ismear** – The smearing parameter used by the `step_func`.



* **Returns**

    Array of size nx with step function on the desired outputgrid.



* **Return type**

    np.array



### pymatgen.io.vasp.optics.kramers_kronig(eps: np.ndarray, nedos: int, deltae: float, cshift: float = 0.1)
Perform the Kramers-Kronig transformation.

Perform the Kramers-Kronig transformation exactly as VASP does it.
The input eps should be complex and the imaginary part of the dielectric function
should be stored as the real part of the complex input array.
The output should be the complex dielectric function.


* **Parameters**


    * **eps** – The dielectric function with the imaginary part stored as the real part and nothing in the imaginary part.


    * **nedos** – The sampling of the energy values


    * **deltae** – The energy grid spacing


    * **cshift** – The shift of the imaginary part of the dielectric function.



* **Returns**

    Array of size nedos with the complex dielectric function.



* **Return type**

    np.array



### pymatgen.io.vasp.optics.step_func(x, ismear)
Replication of VASP’s step function.


### pymatgen.io.vasp.optics.step_methfessel_paxton(x, n)
S_n (x) = (1 + erf x)/2 - exp -x^2 \* sum_i=1^n A_i H_{2i-1}(x)
where H is a Hermite polynomial and
A_i = (-1)^i / ( i! 4^i sqrt(pi) ).