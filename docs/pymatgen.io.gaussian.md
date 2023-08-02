---
layout: default
title: pymatgen.io.gaussian.md
nav_exclude: true
---

# pymatgen.io.gaussian module

This module implements input and output processing from Gaussian.


### _class_ pymatgen.io.gaussian.GaussianInput(mol, charge=None, spin_multiplicity=None, title=None, functional='HF', basis_set='6-31G(d)', route_parameters=None, input_parameters=None, link0_parameters=None, dieze_tag='#P', gen_basis=None)
Bases: `object`

An object representing a Gaussian input file.


* **Parameters**


    * **mol** – Input molecule. It can either be a Molecule object,
    a string giving the geometry in a format supported by Gaussian,
    or `None`. If the molecule is `None`, you will need to use
    read it in from a checkpoint. Consider adding `CHK` to the
    `link0_parameters`.


    * **charge** – Charge of the molecule. If None, charge on molecule is used.
    Defaults to None. This allows the input file to be set a
    charge independently from the molecule itself.
    If `mol` is not a Molecule object, then you must specify a charge.


    * **spin_multiplicity** – Spin multiplicity of molecule. Defaults to None,
    which means that the spin multiplicity is set to 1 if the
    molecule has no unpaired electrons and to 2 if there are
    unpaired electrons. If `mol` is not a Molecule object, then you

    > must specify the multiplicity



    * **title** – Title for run. Defaults to formula of molecule if None.


    * **functional** – Functional for run.


    * **basis_set** – Basis set for run.


    * **route_parameters** – Additional route parameters as a dict. For example,
    {‘SP’:””, “SCF”:”Tight”}


    * **input_parameters** – Additional input parameters for run as a dict. Used
    for example, in PCM calculations. E.g., {“EPS”:12}


    * **link0_parameters** – Link0 parameters as a dict. E.g., {“%mem”: “1000MW”}


    * **dieze_tag** – # preceding the route line. E.g. “#p”


    * **gen_basis** – allows a user-specified basis set to be used in a Gaussian
    calculation. If this is not None, the attribute `basis_set` will
    be set to “Gen”.



#### as_dict()

* **Returns**

    MSONable dict



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – dict



* **Returns**

    GaussianInput



#### _static_ from_file(filename)
Creates GaussianInput from a file.


* **Parameters**

    **filename** – Gaussian input filename



* **Returns**

    GaussianInput object



#### _static_ from_str(contents)
Creates GaussianInput from a string.


* **Parameters**

    **contents** – String representing an Gaussian input file.



* **Returns**

    GaussianInput object



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_cart_coords()
Return the Cartesian coordinates of the molecule.


#### get_zmatrix()
Returns a z-matrix representation of the molecule.


#### _property_ molecule()
Returns molecule associated with this GaussianInput.


#### to_str(cart_coords=False)
Return GaussianInput string.


* **Parameters**

    **cart_coords** (*bool*) – If True, return Cartesian coordinates instead of z-matrix.
    Defaults to False.



#### to_string(\*\*kwds)
to_string is deprecated!
Use to_str instead


#### write_file(filename, cart_coords=False)
Write the input string into a file.

Option: see __str__ method


### _class_ pymatgen.io.gaussian.GaussianOutput(filename)
Bases: `object`

Parser for Gaussian output files.

**NOTE**: Still in early beta.

Attributes:
.. attribute:: structures

> All structures from the calculation in the standard orientation. If the
> symmetry is not considered, the standard orientation is not printed out
> and the input orientation is used instead. Check the standard_orientation
> attribute.


#### structures_input_orientation()
All structures from the calculation in the input orientation or the
Z-matrix orientation (if an opt=z-matrix was requested).


#### opt_structures()
All optimized structures from the calculation in the standard orientation,
if the attribute ‘standard_orientation’ is True, otherwise in the input
or the Z-matrix orientation.


#### energies()
All energies from the calculation.


#### eigenvalues()
List of eigenvalues for the last geometry


#### MO_coefficients()
Matrix of MO coefficients for the last geometry


#### cart_forces()
All Cartesian forces from the calculation.


#### frequencies()
A list for each freq calculation and for each mode of a dict with
{

> > “frequency”: freq in cm-1,
> > “symmetry”: symmetry tag
> > “r_mass”: Reduce mass,
> > “f_constant”: force constant,
> > “IR_intensity”: IR Intensity,
> > “mode”: normal mode

> }

The normal mode is a 1D vector of dx, dy dz of each atom.


#### hessian()
Matrix of second derivatives of the energy with respect to cartesian
coordinates in the **input orientation** frame. Need #P in the
route section in order to be in the output.


#### properly_terminated()
True if run has properly terminated


#### is_pcm()
True if run is a PCM run.


#### is_spin()
True if it is an unrestricted run


#### stationary_type()
If it is a relaxation run, indicates whether it is a minimum (Minimum)
or a saddle point (“Saddle”).


#### corrections()
Thermochemical corrections if this run is a Freq run as a dict. Keys
are “Zero-point”, “Thermal”, “Enthalpy” and “Gibbs Free Energy”


#### functional()
Functional used in the run.


#### basis_set()
Basis set used in the run


#### route()
Additional route parameters as a dict. For example,

    {‘SP’:””, “SCF”:”Tight”}


#### dieze_tag()
# preceding the route line, e.g. “#P”


#### link0()
Link0 parameters as a dict. E.g., {“%mem”: “1000MW”}


#### charge()
Charge for structure


#### spin_multiplicity()
Spin multiplicity for structure


#### num_basis_func()
Number of basis functions in the run.


#### electrons()
number of alpha and beta electrons as (N alpha, N beta)


#### pcm()
PCM parameters and output if available.


#### errors()
error if not properly terminated (list to be completed in error_defs)


#### Mulliken_charges()
Mulliken atomic charges


#### eigenvectors()
Matrix of shape (num_basis_func, num_basis_func). Each column is an
eigenvectors and contains AO coefficients of an MO.

eigenvectors[Spin] = mat(num_basis_func, num_basis_func)


#### molecular_orbital()
MO development coefficients on AO in a more convenient array dict
for each atom and basis set label.

mo[Spin][OM j][atom i] = {AO_k: coeff, AO_k: coeff … }


#### atom_basis_labels()
Labels of AO for each atoms. These labels are those used in the output
of molecular orbital coefficients (POP=Full) and in the
molecular_orbital array dict.

atom_basis_labels[iatom] = [AO_k, AO_k, …]


#### resumes()
List of gaussian data resume given at the end of the output file before
the quotation. The resumes are given as string.


#### title()
Title of the gaussian run.


#### standard_orientation()
If True, the geometries stored in the structures are in the standard
orientation. Else, the geometries are in the input orientation.


#### bond_orders()
Dict of bond order values read in the output file such as:
{(0, 1): 0.8709, (1, 6): 1.234, …}

The keys are the atom indexes and the values are the Wiberg bond indexes
that are printed using pop=NBOREAD and $nbo bndidx $end.

Methods:
.. method:: to_input()

> Return a GaussianInput object using the last geometry and the same
> calculation parameters.


#### read_scan()
Read a potential energy surface from a gaussian scan calculation.


#### get_scan_plot()
Get a matplotlib plot of the potential energy surface


#### save_scan_plot()
Save a matplotlib plot of the potential energy surface to a file


* **Parameters**

    **filename** – Filename of Gaussian output file.



#### as_dict()
JSON-serializable dict representation.


#### _property_ final_energy()
Final energy in Gaussian output.


* **Type**

    return



#### _property_ final_structure()
Final structure in Gaussian output.


* **Type**

    return



#### get_scan_plot(coords=None)
Get a matplotlib plot of the potential energy surface.


* **Parameters**

    **coords** – internal coordinate name to use as abscissa.



#### get_spectre_plot(sigma=0.05, step=0.01)
Get a matplotlib plot of the UV-visible xas. Transitions are plotted
as vertical lines and as a sum of normal functions with sigma with. The
broadening is applied in energy and the xas is plotted as a function
of the wavelength.


* **Parameters**


    * **sigma** – Full width at half maximum in eV for normal functions.


    * **step** – bin interval in eV



* **Returns**

    {“energies”: values, “lambda”: values, “xas”: values}

        where values are lists of abscissa (energies, lamba) and
        the sum of gaussian functions (xas).

    A matplotlib plot.




* **Return type**

    A dict



#### read_excitation_energies()
Read a excitation energies after a TD-DFT calculation.


* **Returns**

    A list of tuple for each transition such as

        [(energie (eV), lambda (nm), oscillatory strength), … ]




* **Return type**

    A list



#### read_scan()
Read a potential energy surface from a gaussian scan calculation.


* **Returns**

    {“energies”: [ values ],

        ”coords”: {“d1”: [ values ], “A2”, [ values ], … }}

    ”energies” are the energies of all points of the potential energy
    surface. “coords” are the internal coordinates used to compute the
    potential energy surface and the internal coordinates optimized,
    labelled by their name as defined in the calculation.




* **Return type**

    A dict



#### save_scan_plot(filename='scan.pdf', img_format='pdf', coords=None)
Save matplotlib plot of the potential energy surface to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **coords** – internal coordinate name to use as abcissa.



#### save_spectre_plot(filename='spectre.pdf', img_format='pdf', sigma=0.05, step=0.01)
Save matplotlib plot of the spectre to a file.


* **Parameters**


    * **filename** – Filename to write to.


    * **img_format** – Image format to use. Defaults to EPS.


    * **sigma** – Full width at half maximum in eV for normal functions.


    * **step** – bin interval in eV



#### to_input(mol=None, charge=None, spin_multiplicity=None, title=None, functional=None, basis_set=None, route_parameters=None, input_parameters=None, link0_parameters=None, dieze_tag=None, cart_coords=False)
Create a new input object using by default the last geometry read in
the output file and with the same calculation parameters. Arguments
are the same as GaussianInput class.


* **Returns**

    the gaussian input object



* **Return type**

    gaunip (GaussianInput)



### pymatgen.io.gaussian.read_route_line(route)
read route line in gaussian input/output and return functional basis_set
and a dictionary of other route parameters.


* **Parameters**

    **route** (*str*) – the route line



* **Returns**

    the method (HF, PBE …)
    basis_set (str) : the basis set
    route (dict) : dictionary of parameters



* **Return type**

    functional (str)