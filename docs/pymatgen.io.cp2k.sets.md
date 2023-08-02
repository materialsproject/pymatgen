---
layout: default
title: pymatgen.io.cp2k.sets.md
nav_exclude: true
---

# pymatgen.io.cp2k.sets module

This module defines input sets for CP2K and is a work in progress. The structure/philosophy
of this module is based on the Vasp input sets in Pymatgen. These sets are meant to contain
tested parameters that will result in successful, reproducible, consistent calculations without
need for intervention 99% of the time. 99% of the time, you only need to provide a pymatgen
structure object and let the defaults take over from there.

The sets are intended to be very general, e.g. a set for geometry relaxation, and so most of the
time, if you have specific needs, you can simply specify them via the keyword argument
override_default_params (see Section.update() method). If you have the need to create a new input
set (say for a standardized high throughput calculation) then you can create a new child of the
Cp2kInputSet class.

In order to implement a new Set within the current code structure, follow this 3 step flow:


    1. Inherit from Cp2kInputSet or one of its children and call the super() constructor


    2. Create the new sections and insert them into self and its subsections as needed


    3. Call self.update(override_default_params) in order to allow user settings.


### _class_ pymatgen.io.cp2k.sets.CellOptSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for cell optimization relaxation.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** ([*Kpoints*](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints)) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



### _exception_ pymatgen.io.cp2k.sets.Cp2kValidationError(message)
Bases: `Exception`

Cp2k Validation Exception. Not exhausted. May raise validation
errors for features which actually do work if using a newer version
of cp2k.


#### CP2K_VERSION(_ = 'v2022.1_ )

### _class_ pymatgen.io.cp2k.sets.DftSet(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule), project_name: str = 'CP2K', basis_and_potential: dict | None = None, xc_functionals: list | str | None = None, multiplicity: int = 0, ot: bool = True, energy_gap: float = -1, qs_method: str = 'GPW', eps_default: float = 1e-12, eps_scf: float = 1e-06, max_scf: int | None = None, minimizer: str = 'DIIS', preconditioner: str = 'FULL_SINGLE_INVERSE', algorithm: str = 'IRAC', linesearch: str = '2PNT', rotation: bool = True, occupation_preconditioner: bool = False, cutoff: int | float | None = None, rel_cutoff: int = 50, ngrids: int = 5, progression_factor: int = 3, override_default_params: dict | None = None, wfn_restart_file_name: str | None = None, kpoints: VaspKpoints | None = None, smearing: bool = False, \*\*kwargs)
Bases: [`Cp2kInput`](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Cp2kInput)

Base for an input set using the Quickstep module (i.e. a DFT calculation). The DFT section is
pretty vast in CP2K, so this set hopes to make the DFT setup fairly simple. The provided
parameters are pretty conservative, and so they should not need to be changed very often.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** ([*Kpoints*](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints)) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



#### activate_epr(\*\*kwargs)
Calculate g-tensor. Requires localize. Suggested with GAPW.


#### activate_fast_minimization(on)
Method to modify the set to use fast SCF minimization.


#### activate_hybrid(hybrid_functional: str = 'PBE0', hf_fraction: float = 0.25, gga_x_fraction: float = 0.75, gga_c_fraction: float = 1, max_memory: int = 2000, cutoff_radius: float = 8.0, potential_type: str | None = None, omega: float = 0.11, scale_coulomb: float = 1, scale_gaussian: float = 1, scale_longrange: float = 1, admm: bool = True, admm_method: str = 'BASIS_PROJECTION', admm_purification_method: str = 'NONE', admm_exch_correction_func: str = 'DEFAULT', eps_schwarz: float = 1e-07, eps_schwarz_forces: float = 1e-06, screen_on_initial_p: bool = True, screen_p_forces: bool = True)
Basic set for activating hybrid DFT calculation using Auxiliary Density Matrix Method.

Note 1: When running ADMM with cp2k, memory is very important. If the memory requirements
exceed what is available (see max_memory), then CP2K will have to calculate the 4-electron
integrals for HFX during each step of the SCF cycle. ADMM provides a huge speed up by
making the memory requirements *feasible* to fit into RAM, which means you only need to
calculate the integrals once each SCF cycle. But, this only works if it fits into memory.
When setting up ADMM calculations, we recommend doing whatever is possible to fit all the
4EI into memory.

Note 2: This set is designed for reliable high-throughput calculations, NOT for extreme
accuracy. Please review the in-line comments in this method if you want more control.


* **Parameters**


    * **hybrid_functional** (*str*) – Type of hybrid functional. This set supports HSE (screened)
    and PBE0 (truncated). Default is PBE0, which converges easier in the GPW basis
    used by cp2k.


    * **hf_fraction** (*float*) – fraction of exact HF exchange energy to mix. Default: 0.25


    * **gga_x_fraction** (*float*) – fraction of gga exchange energy to retain. Default: 0.75


    * **gga_c_fraction** (*float*) – fraction of gga correlation energy to retain. Default: 1.0


    * **max_memory** (*int*) – Maximum memory available to each MPI process (in Mb) in the
    calculation. Most modern computing nodes will have ~2Gb per core, or 2048 Mb,
    but check for your specific system. This value should be as large as possible
    while still leaving some memory for the other parts of cp2k. Important: If
    this value is set larger than the memory limits, CP2K will likely seg-fault.
    Default: 2000


    * **cutoff_radius** (*float*) – for truncated hybrid functional (i.e. PBE0), this is the cutoff
    radius. The default is selected as that which generally gives convergence, but
    maybe too low (if you want very high accuracy) or too high (if you want a quick
    screening). Default: 8 angstroms


    * **potential_type** (*str*) – what interaction potential to use for HFX. Available in CP2K are
    COULOMB, GAUSSIAN, IDENTITY, LOGRANGE, MIX_CL, MIX_CL_TRUNC, MIX_LG, SHORTRANGE,
    and TRUNCATED. Default is None, and it will be set automatically depending on the
    named hybrid_functional that you use, but setting it to one of the acceptable
    values will constitute a user-override.


    * **omega** (*float*) – For HSE, this specifies the screening parameter. HSE06 sets this as
    0.2, which is the default.


    * **scale_coulomb** – Scale for the coulomb operator if using a range separated functional


    * **scale_gaussian** – Scale for the gaussian operator (if applicable)


    * **scale_longrange** – Scale for the coulomb operator if using a range separated functional


    * **admm** – Whether or not to use the auxiliary density matrix method for the exact
    HF exchange contribution. Highly recommended. Speed ups between 10x and 1000x are
    possible when compared to non ADMM hybrid calculations.


    * **admm_method** – Method for constructing the auxiliary basis


    * **admm_purification_method** – Method for purifying the auxiliary density matrix so as to
    preserve properties, such as idempotency. May lead to shifts in the
    eigenvalues.


    * **admm_exch_correction_func** – Which functional to use to calculate the exchange correction
    E_x(primary) - E_x(aux)


    * **eps_schwarz** – Screening threshold for HFX, in Ha. Contributions smaller than
    this will be screened. The smaller the value, the more accurate, but also the more
    costly. Default value is 1e-7. 1e-6 works in a large number of cases, but is
    quite aggressive, which can lead to convergence issues.


    * **eps_schwarz_forces** – Same as for eps_schwarz, but for screening contributions to
    forces. Convergence is not as sensitive with respect to eps_schwarz forces as
    compared to eps_schwarz, and so 1e-6 should be good default.


    * **screen_on_initial_p** – If an initial density matrix is provided, in the form of a
    CP2K wfn restart file, then this initial density will be used for screening. This
    is generally very computationally efficient, but, as with eps_schwarz, can lead to
    instabilities if the initial density matrix is poor.


    * **screen_p_forces** – Same as screen_on_initial_p, but for screening of forces.



#### activate_hyperfine()
Print the hyperfine coupling constants.


#### activate_localize(states='OCCUPIED', preconditioner='FULL_ALL', restart=False)
Activate calculation of the maximally localized wannier functions.


* **Parameters**


    * **states** – Which states to calculate. occupied, unoccupied, mixed states, or all states. At
    present, unoccupied orbitals are only implemented for GPW.


    * **preconditioner** – Preconditioner to use for optimize


    * **restart** – Initialize from the localization restart file



#### activate_motion(max_drift: float = 0.003, rms_drift: float = 0.0015, max_force: float = 0.00045, rms_force: float = 0.0003, max_iter: int = 200, optimizer: str = 'BFGS', trust_radius: float = 0.25, line_search: str = '2PNT', ensemble: str = 'NVE', temperature: float | int = 300, timestep: float | int = 0.5, nsteps: int = 3, thermostat: str = 'NOSE', nproc_rep: int = 1)
Turns on the motion section for GEO_OPT, CELL_OPT, etc. calculations.
Will turn on the printing subsections and also bind any constraints
to their respective atoms.


#### activate_nmr(\*\*kwargs)
Calculate nmr shifts. Requires localize. Suggested with GAPW.


#### activate_nonperiodic(solver='ANALYTIC')
Activates a calculation with non-periodic calculations by turning of PBC and
changing the poisson solver. Still requires a CELL to put the atoms.


#### activate_polar(\*\*kwargs)
Calculate polarizations (including raman).


#### activate_robust_minimization()
Method to modify the set to use more robust SCF minimization technique.


#### activate_spinspin(\*\*kwargs)
Calculate spin-spin coupling tensor. Requires localize.


#### activate_tddfpt(\*\*kwargs)
Activate TDDFPT for calculating excited states. Only works with GPW. Supports hfx.


#### activate_vdw_potential(dispersion_functional: str, potential_type: str)
Activate van der Waals dispersion corrections.


* **Parameters**


    * **dispersion_functional** – Type of dispersion functional.
    Options: pair_potential or non_local


    * **potential_type** – What type of potential to use, given a dispersion functional type
    Options: DFTD2, DFTD3, DFTD3(BJ), DRSLL, LMKLL, RVV10



#### activate_very_strict_minimization()
Method to modify the set to use very strict SCF minimization scheme
:return:


#### create_subsys(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure) | [Molecule](pymatgen.core.structure.md#pymatgen.core.structure.Molecule))
Create the structure for the input.


#### _static_ get_basis_and_potential(structure, basis_and_potential)
Get a dictionary of basis and potential info for constructing the input file.

data in basis_and_potential argument can be specified in several ways:

> Strategy 1: Element-specific info (takes precedence)

> >
> > 1. Provide a basis and potential object:

> > > el: {‘basis’: obj, ‘potential’: obj}


> > 2. Provide a hash of the object that matches the keys in the pmg configured cp2k data files.

> > > el: {‘basis’: hash, ‘potential’: hash}

> > 3. Provide the name of the basis and potential AND the basis_filenames and potential_filename
> > keywords specifying where to find these objects

> > > el: {

> > >     ‘basis’: name, ‘potential’: name, ‘basis_filenames’: [filenames],
> > >     ‘potential_filename’: filename

> > > }

> Strategy 2: global descriptors

> > In this case, any elements not present in the argument will be dealt with by searching the pmg
> > configured cp2k data files to find a objects matching your requirements.

> > >
> > > * functional: Find potential and basis that have been optimized for a specific functional like PBE.

> > >     Can be None if you do not require them to match.


> > > * basis_type: type of basis to search for (e.g. DZVP-MOLOPT).


> > > * aux_basis_type: type of basis to search for (e.g. pFIT). Some elements do not have all aux types

> > >     available. Use aux_basis_type that is universal to avoid issues, or avoid using this strategy.


> > > * potential_type: “Pseudopotential” or “All Electron”

> > **\*BE WARNED\*** CP2K data objects can have the same name, this will sort those and choose the first one
> > that matches.

Will raise an error if no basis/potential info can be found according to the input.


#### _static_ get_cutoff_from_basis(basis_sets, rel_cutoff)
Given a basis and a relative cutoff. Determine the ideal cutoff variable.


#### _static_ get_xc_functionals(xc_functionals: list | str | None = None)
Get XC functionals. If simplified names are provided in kwargs, they
will be expanded into their corresponding X and C names.


#### modify_dft_print_iters(iters, add_last='no')
Modify all DFT print iterations at once. Common use is to set iters to the max
number of iterations + 1 and then set add_last to numeric. This would have the
effect of printing only the first and last iteration, which might be useful for
speeding up/saving space on GEO_OPT or MD runs where you don’t need the intermediate
values.


* **Parameters**


    * **iters** (*int*) – print each “iters” iterations.


    * **add_last** (*str*) – Whether to explicitly include the last iteration, and how to mark it.
    numeric: mark last iteration with the iteration number
    symbolic: mark last iteration with the letter “l”
    no: do not explicitly include the last iteration



#### print_bandstructure(kpoints_line_density: int = 20)
Attaches a non-scf band structure calc the end of an SCF loop.

This requires a kpoint calculation, which is not always default in cp2k.


* **Parameters**

    **kpoints_line_density** – number of kpoints along each branch in line-mode calc.



#### print_dos(ndigits=6)
Activate printing of the overall DOS file.

Note: As of 2022.1, ndigits needs to be set to a sufficient value to ensure data is not lost.
Note: As of 2022.1, can only be used with a k-point calculation.


#### print_e_density(stride=(2, 2, 2))
Controls the printing of cube files with electronic density and, for UKS, the spin density.


#### print_forces()
Print out the forces and stress during calculation.


#### print_hirshfeld(on=True)
Activate or deactivate printing of Hirshfeld charges.


#### print_ldos(nlumo: int = -1)
Activate the printing of LDOS files, printing one for each atom kind by default.


* **Parameters**

    **nlumo** (*int*) – Number of virtual orbitals to be added to the MO set (-1=all).
    CAUTION: Setting this value to be higher than the number of states present may
    cause a Cholesky error.



#### print_mo()
Print molecular orbitals when running non-OT diagonalization.


#### print_mo_cubes(write_cube: bool = False, nlumo: int = -1, nhomo: int = -1)
Activate printing of molecular orbitals.


* **Parameters**


    * **write_cube** (*bool*) – whether to write cube file for the MOs instead of out file


    * **nlumo** (*int*) – Controls the number of lumos printed and dumped as a cube (-1=all)


    * **nhomo** (*int*) – Controls the number of homos printed and dumped as a cube (-1=all)



#### print_mulliken(on=False)
Activate or deactivate printing of Mulliken charges.


#### print_pdos(nlumo: int = -1)
Activate creation of the PDOS file.


* **Parameters**

    **nlumo** (*int*) – Number of virtual orbitals to be added to the MO set (-1=all).
    CAUTION: Setting this value to be higher than the number of states present may
    cause a Cholesky error.



#### print_v_hartree(stride=(2, 2, 2))
Controls the printing of a cube file with eletrostatic potential generated by the

    total density (electrons+ions). It is valid only for QS with GPW formalism.

Note that by convention the potential has opposite sign than the expected physical one.


#### set_charge(charge: int)
Set the overall charge of the simulation cell.


#### validate()
Implements a few checks for a valid input set.


#### write_basis_set_file(basis_sets, fn='BASIS')
Write the basis sets to a file.


#### write_potential_file(potentials, fn='POTENTIAL')
Write the potentials to a file.


### _class_ pymatgen.io.cp2k.sets.HybridCellOptSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for hybrid cell optimization relaxation.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** ([*Kpoints*](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints)) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



### _class_ pymatgen.io.cp2k.sets.HybridRelaxSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for hybrid geometry relaxation.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** ([*Kpoints*](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints)) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



### _class_ pymatgen.io.cp2k.sets.HybridStaticSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for static calculations.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** ([*Kpoints*](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints)) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



### _class_ pymatgen.io.cp2k.sets.RelaxSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for geometry relaxation.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** ([*Kpoints*](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints)) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.



### _class_ pymatgen.io.cp2k.sets.StaticSet(\*\*kwargs)
Bases: `DftSet`

Quick Constructor for static calculations.


* **Parameters**


    * **structure** – Pymatgen structure or molecule object


    * **ot** (*bool*) – Whether or not to use orbital transformation method for matrix
    diagonalization. OT is the flagship scf solver of CP2K, and will provide
    speed-ups for this part of the calculation, but the system must have a band gap
    for OT to be used (higher band-gap –> faster convergence).


    * **energy_gap** (*float*) – Estimate of energy gap for pre-conditioner. Default is -1, leaving
    it up to cp2k.


    * **eps_default** (*float*) – Replaces all EPS_XX Keywords in the DFT section value, ensuring
    an overall accuracy of at least this much.


    * **eps_scf** (*float*) – The convergence criteria for leaving the SCF loop. Default is 1e-6.
    Should ensure reasonable results, but is not applicable to all situations.

    > Note: eps_scf is *not* in units of energy, as in most DFT codes. For OT method,
    > it is the largest gradient of the energy with respect to changing any of the
    > molecular orbital coefficients. For diagonalization, it is the largest change
    > in the density matrix from the last step.



    * **max_scf** (*int*) – The max number of SCF cycles before terminating the solver. NOTE: With
    the OT solver, this corresponds to the max number of INNER scf loops, and then
    the outer loops are set with outer_max_scf, while with diagonalization it
    corresponds to the overall (INNER\*OUTER) number of SCF steps, with the
    inner loop limit set by


    * **minimizer** (*str*) – The minimization scheme. DIIS can be as much as 50% faster than the
    more robust conjugate gradient method, and so it is chosen as default. Switch to CG
    if dealing with a difficult system.


    * **preconditioner** (*str*) – Pre-conditioner for the OT method. FULL_SINGLE_INVERSE is very
    robust and compatible with non-integer occupations from IRAC+rotation. FULL_ALL is
    considered “best” but needs algorithm to be set to STRICT. Only change from these
    two when simulation cell gets to be VERY large, in which case FULL_KINETIC might be
    preferred.


    * **algorithm** (*str*) – Algorithm for the OT method. STRICT assumes that the orbitals are
    strictly orthogonal to each other, which works well for wide gap ionic systems,
    but can diverge for systems with small gaps, fractional occupations, and some
    other cases. IRAC (iterative refinement of the approximate congruency)
    transformation is not analytically correct and uses a truncated polynomial
    expansion, but is robust to the problems with STRICT, and so is the default.


    * **linesearch** (*str*) – Linesearch method for CG. 2PNT is the default, and is the fastest,
    but is not as robust as 3PNT. 2PNT is required as of cp2k v9.1 for compatibility
    with irac+rotation. This may be upgraded in the future. 3PNT can be good for wide
    gapped transition metal systems as an alternative.


    * **rotation** (*bool*) – Whether or not to allow for rotation of the orbitals in the OT method.
    This equates to allowing for fractional occupations in the calculation.


    * **occupation_preconditioner** (*bool*) – Whether or not to account for fractional occupations
    in the preconditioner. This method is not fully integrated as of cp2k v9.1 and is
    set to false by default.


    * **cutoff** (*int*) – Cutoff energy (in Ry) for the finest level of the multigrid. A high
    cutoff will allow you to have very accurate calculations PROVIDED that REL_CUTOFF
    is appropriate. By default cutoff is set to 0, leaving it up to the set.


    * **rel_cutoff** (*int*) – This cutoff decides how the Gaussians are mapped onto the different
    levels of the multigrid. If REL_CUTOFF is too low, then even if you have a high
    CUTOFF, all Gaussians will be mapped onto the coarsest level of the multi-grid,
    and thus the effective integration grid for the calculation may still be too
    coarse. By default 50Ry is chosen, which should be sufficient given the cutoff is
    large enough.


    * **ngrids** (*int*) – number of multi-grids to use. CP2K default is 4, but the molopt basis
    files recommend 5.


    * **progression_factor** (*int*) – Divisor of CUTOFF to get the cutoff for the next level of
    the multigrid.


    * **wfn_restart_file_name** (*str*) – RESTART file for the initial wavefunction guess.


    * **kpoints** ([*Kpoints*](pymatgen.io.cp2k.inputs.md#pymatgen.io.cp2k.inputs.Kpoints)) – kpoints object from pymatgen.io.vasp.inputs.Kpoints. By default,
    CP2K runs with gamma point only.


    * **smearing** (*bool*) – whether or not to activate smearing (should be done for systems
    containing no (or a very small) band gap.