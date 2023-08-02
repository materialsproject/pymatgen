---
layout: default
title: pymatgen.io.vasp.md
nav_exclude: true
---

# pymatgen.io.vasp package

This package implements modules for input and output to and from VASP. It
imports the key classes form both vasp_input and vasp_output to allow most
classes to be simply called as pymatgen.io.vasp.Incar for example, to retain
backwards compatibility.



* [pymatgen.io.vasp.help module](pymatgen.io.vasp.help.md)


    * [`VaspDoc`](pymatgen.io.vasp.help.md#pymatgen.io.vasp.help.VaspDoc)


        * [`VaspDoc.get_help()`](pymatgen.io.vasp.help.md#pymatgen.io.vasp.help.VaspDoc.get_help)


        * [`VaspDoc.get_incar_tags()`](pymatgen.io.vasp.help.md#pymatgen.io.vasp.help.VaspDoc.get_incar_tags)


        * [`VaspDoc.print_help()`](pymatgen.io.vasp.help.md#pymatgen.io.vasp.help.VaspDoc.print_help)


        * [`VaspDoc.print_jupyter_help()`](pymatgen.io.vasp.help.md#pymatgen.io.vasp.help.VaspDoc.print_jupyter_help)


* [pymatgen.io.vasp.inputs module](pymatgen.io.vasp.inputs.md)


    * [`BadIncarWarning`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.BadIncarWarning)


    * [`Incar`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar)


        * [`Incar.as_dict()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar.as_dict)


        * [`Incar.check_params()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar.check_params)


        * [`Incar.diff()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar.diff)


        * [`Incar.from_dict()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar.from_dict)


        * [`Incar.from_file()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar.from_file)


        * [`Incar.from_str()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar.from_str)


        * [`Incar.from_string()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar.from_string)


        * [`Incar.get_string()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar.get_string)


        * [`Incar.proc_val()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar.proc_val)


        * [`Incar.write_file()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Incar.write_file)


    * [`Kpoints`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints)


        * [`Kpoints.as_dict()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.as_dict)


        * [`Kpoints.automatic()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.automatic)


        * [`Kpoints.automatic_density()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.automatic_density)


        * [`Kpoints.automatic_density_by_lengths()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.automatic_density_by_lengths)


        * [`Kpoints.automatic_density_by_vol()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.automatic_density_by_vol)


        * [`Kpoints.automatic_gamma_density()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.automatic_gamma_density)


        * [`Kpoints.automatic_linemode()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.automatic_linemode)


        * [`Kpoints.from_dict()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.from_dict)


        * [`Kpoints.from_file()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.from_file)


        * [`Kpoints.from_str()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.from_str)


        * [`Kpoints.from_string()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.from_string)


        * [`Kpoints.gamma_automatic()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.gamma_automatic)


        * [`Kpoints.monkhorst_automatic()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.monkhorst_automatic)


        * [`Kpoints.style`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.style)


        * [`Kpoints.supported_modes`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.supported_modes)


        * [`Kpoints.write_file()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Kpoints.write_file)


    * [`KpointsSupportedModes`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.KpointsSupportedModes)


        * [`KpointsSupportedModes.Automatic`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Automatic)


        * [`KpointsSupportedModes.Cartesian`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Cartesian)


        * [`KpointsSupportedModes.Gamma`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Gamma)


        * [`KpointsSupportedModes.Line_mode`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Line_mode)


        * [`KpointsSupportedModes.Monkhorst`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Monkhorst)


        * [`KpointsSupportedModes.Reciprocal`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.Reciprocal)


        * [`KpointsSupportedModes.from_str()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.from_str)


        * [`KpointsSupportedModes.from_string()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.KpointsSupportedModes.from_string)


    * [`Orbital`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Orbital)


        * [`Orbital.E`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Orbital.E)


        * [`Orbital.j`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Orbital.j)


        * [`Orbital.l`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Orbital.l)


        * [`Orbital.n`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Orbital.n)


        * [`Orbital.occ`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Orbital.occ)


    * [`OrbitalDescription`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.OrbitalDescription)


        * [`OrbitalDescription.E`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.OrbitalDescription.E)


        * [`OrbitalDescription.Rcut`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.OrbitalDescription.Rcut)


        * [`OrbitalDescription.Rcut2`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.OrbitalDescription.Rcut2)


        * [`OrbitalDescription.Type`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.OrbitalDescription.Type)


        * [`OrbitalDescription.Type2`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.OrbitalDescription.Type2)


        * [`OrbitalDescription.l`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.OrbitalDescription.l)


    * [`Poscar`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar)


        * [`Poscar.structure`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.structure)


        * [`Poscar.comment`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.comment)


        * [`Poscar.true_names`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.true_names)


        * [`Poscar.selective_dynamics`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.selective_dynamics)


        * [`Poscar.velocities`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.velocities)


        * [`Poscar.predictor_corrector`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.predictor_corrector)


        * [`Poscar.predictor_corrector_preamble`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.predictor_corrector_preamble)


        * [`Poscar.temperature`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.temperature)


        * [`Poscar.as_dict()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.as_dict)


        * [`Poscar.from_dict()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.from_dict)


        * [`Poscar.from_file()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.from_file)


        * [`Poscar.from_str()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.from_str)


        * [`Poscar.from_string()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.from_string)


        * [`Poscar.get_string()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.get_string)


        * [`Poscar.natoms`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.natoms)


        * [`Poscar.predictor_corrector`](pymatgen.io.vasp.inputs.md#id0)


        * [`Poscar.selective_dynamics`](pymatgen.io.vasp.inputs.md#id1)


        * [`Poscar.set_temperature()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.set_temperature)


        * [`Poscar.site_symbols`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.site_symbols)


        * [`Poscar.velocities`](pymatgen.io.vasp.inputs.md#id2)


        * [`Poscar.write_file()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Poscar.write_file)


    * [`Potcar`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Potcar)


        * [`Potcar.FUNCTIONAL_CHOICES`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Potcar.FUNCTIONAL_CHOICES)


        * [`Potcar.as_dict()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Potcar.as_dict)


        * [`Potcar.from_dict()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Potcar.from_dict)


        * [`Potcar.from_file()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Potcar.from_file)


        * [`Potcar.set_symbols()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Potcar.set_symbols)


        * [`Potcar.spec`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Potcar.spec)


        * [`Potcar.symbols`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Potcar.symbols)


        * [`Potcar.write_file()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.Potcar.write_file)


    * [`PotcarSingle`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle)


        * [`PotcarSingle.data`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.data)


        * [`PotcarSingle.keywords`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.keywords)


        * [`PotcarSingle.atomic_no`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.atomic_no)


        * [`PotcarSingle.electron_configuration`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.electron_configuration)


        * [`PotcarSingle.element`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.element)


        * [`PotcarSingle.from_file()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.from_file)


        * [`PotcarSingle.from_symbol_and_functional()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.from_symbol_and_functional)


        * [`PotcarSingle.functional`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.functional)


        * [`PotcarSingle.functional_class`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.functional_class)


        * [`PotcarSingle.functional_dir`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.functional_dir)


        * [`PotcarSingle.functional_tags`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.functional_tags)


        * [`PotcarSingle.get_potcar_file_hash()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.get_potcar_file_hash)


        * [`PotcarSingle.get_potcar_hash()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.get_potcar_hash)


        * [`PotcarSingle.get_sha256_file_hash()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.get_sha256_file_hash)


        * [`PotcarSingle.identify_potcar()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.identify_potcar)


        * [`PotcarSingle.nelectrons`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.nelectrons)


        * [`PotcarSingle.parse_functions`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.parse_functions)


        * [`PotcarSingle.potential_type`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.potential_type)


        * [`PotcarSingle.symbol`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.symbol)


        * [`PotcarSingle.verify_potcar()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.verify_potcar)


        * [`PotcarSingle.write_file()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.PotcarSingle.write_file)


    * [`UnknownPotcarWarning`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.UnknownPotcarWarning)


    * [`VaspInput`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.VaspInput)


        * [`VaspInput.as_dict()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.VaspInput.as_dict)


        * [`VaspInput.from_dict()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.VaspInput.from_dict)


        * [`VaspInput.from_directory()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.VaspInput.from_directory)


        * [`VaspInput.run_vasp()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.VaspInput.run_vasp)


        * [`VaspInput.write_input()`](pymatgen.io.vasp.inputs.md#pymatgen.io.vasp.inputs.VaspInput.write_input)


* [pymatgen.io.vasp.optics module](pymatgen.io.vasp.optics.md)


    * [`DielectricFunctionCalculator`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator)


        * [`DielectricFunctionCalculator.cder`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.cder)


        * [`DielectricFunctionCalculator.cder_imag`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.cder_imag)


        * [`DielectricFunctionCalculator.cder_real`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.cder_real)


        * [`DielectricFunctionCalculator.cshift`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.cshift)


        * [`DielectricFunctionCalculator.deltae`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.deltae)


        * [`DielectricFunctionCalculator.efermi`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.efermi)


        * [`DielectricFunctionCalculator.eigs`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.eigs)


        * [`DielectricFunctionCalculator.from_directory()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.from_directory)


        * [`DielectricFunctionCalculator.from_vasp_objects()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.from_vasp_objects)


        * [`DielectricFunctionCalculator.get_epsilon()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.get_epsilon)


        * [`DielectricFunctionCalculator.ismear`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.ismear)


        * [`DielectricFunctionCalculator.ispin`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.ispin)


        * [`DielectricFunctionCalculator.kweights`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.kweights)


        * [`DielectricFunctionCalculator.nedos`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.nedos)


        * [`DielectricFunctionCalculator.plot_weighted_transition_data()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.plot_weighted_transition_data)


        * [`DielectricFunctionCalculator.sigma`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.sigma)


        * [`DielectricFunctionCalculator.volume`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.DielectricFunctionCalculator.volume)


    * [`delta_func()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.delta_func)


    * [`delta_methfessel_paxton()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.delta_methfessel_paxton)


    * [`epsilon_imag()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.epsilon_imag)


    * [`get_delta()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.get_delta)


    * [`get_step()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.get_step)


    * [`kramers_kronig()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.kramers_kronig)


    * [`step_func()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.step_func)


    * [`step_methfessel_paxton()`](pymatgen.io.vasp.optics.md#pymatgen.io.vasp.optics.step_methfessel_paxton)


* [pymatgen.io.vasp.outputs module](pymatgen.io.vasp.outputs.md)


    * [`BSVasprun`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.BSVasprun)


        * [`BSVasprun.as_dict()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.BSVasprun.as_dict)


    * [`Chgcar`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Chgcar)


        * [`Chgcar.from_file()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Chgcar.from_file)


        * [`Chgcar.net_magnetization`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Chgcar.net_magnetization)


    * [`Dynmat`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Dynmat)


        * [`Dynmat.data`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Dynmat.data)


        * [`Dynmat.get_phonon_frequencies()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Dynmat.get_phonon_frequencies)


        * [`Dynmat.masses`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Dynmat.masses)


        * [`Dynmat.natoms`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Dynmat.natoms)


        * [`Dynmat.ndisps`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Dynmat.ndisps)


        * [`Dynmat.nspecs`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Dynmat.nspecs)


    * [`Eigenval`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval)


        * [`Eigenval.filename`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval.filename)


        * [`Eigenval.occu_tol`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval.occu_tol)


        * [`Eigenval.ispin`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval.ispin)


        * [`Eigenval.nelect`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval.nelect)


        * [`Eigenval.nkpt`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval.nkpt)


        * [`Eigenval.nbands`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval.nbands)


        * [`Eigenval.kpoints`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval.kpoints)


        * [`Eigenval.kpoints_weights`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval.kpoints_weights)


        * [`Eigenval.eigenvalues`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval.eigenvalues)


        * [`Eigenval.eigenvalue_band_properties`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Eigenval.eigenvalue_band_properties)


    * [`Elfcar`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Elfcar)


        * [`Elfcar.from_file()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Elfcar.from_file)


        * [`Elfcar.get_alpha()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Elfcar.get_alpha)


    * [`Locpot`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Locpot)


        * [`Locpot.from_file()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Locpot.from_file)


    * [`Oszicar`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Oszicar)


        * [`Oszicar.electronic_steps`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Oszicar.electronic_steps)


        * [`Oszicar.all_energies`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Oszicar.all_energies)


        * [`Oszicar.as_dict()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Oszicar.as_dict)


        * [`Oszicar.final_energy`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Oszicar.final_energy)


    * [`Outcar`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar)


        * [`Outcar.magnetization`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.magnetization)


        * [`Outcar.chemical_shielding`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.chemical_shielding)


        * [`Outcar.unsym_cs_tensor`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.unsym_cs_tensor)


        * [`Outcar.cs_g0_contribution`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.cs_g0_contribution)


        * [`Outcar.cs_core_contribution`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.cs_core_contribution)


        * [`Outcar.efg`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.efg)


        * [`Outcar.charge`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.charge)


        * [`Outcar.is_stopped`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.is_stopped)


        * [`Outcar.run_stats`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.run_stats)


        * [`Outcar.elastic_tensor`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.elastic_tensor)


        * [`Outcar.drift`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.drift)


        * [`Outcar.ngf`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.ngf)


        * [`Outcar.efermi`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.efermi)


        * [`Outcar.filename`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.filename)


        * [`Outcar.final_energy`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.final_energy)


        * [`Outcar.final_energy_wo_entrp`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.final_energy_wo_entrp)


        * [`Outcar.final_fr_energy`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.final_fr_energy)


        * [`Outcar.has_onsite_density_matrices`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.has_onsite_density_matrices)


        * [`Outcar.lcalcpol`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.lcalcpol)


        * [`Outcar.lepsilon`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.lepsilon)


        * [`Outcar.nelect`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.nelect)


        * [`Outcar.spin`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.spin)


        * [`Outcar.total_mag`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.total_mag)


        * [`Outcar.as_dict()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.as_dict)


        * [`Outcar.read_avg_core_poten()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_avg_core_poten)


        * [`Outcar.read_chemical_shielding()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_chemical_shielding)


        * [`Outcar.read_core_state_eigen()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_core_state_eigen)


        * [`Outcar.read_corrections()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_corrections)


        * [`Outcar.read_cs_core_contribution()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_cs_core_contribution)


        * [`Outcar.read_cs_g0_contribution()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_cs_g0_contribution)


        * [`Outcar.read_cs_raw_symmetrized_tensors()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_cs_raw_symmetrized_tensors)


        * [`Outcar.read_elastic_tensor()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_elastic_tensor)


        * [`Outcar.read_electrostatic_potential()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_electrostatic_potential)


        * [`Outcar.read_fermi_contact_shift()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_fermi_contact_shift)


        * [`Outcar.read_freq_dielectric()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_freq_dielectric)


        * [`Outcar.read_igpar()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_igpar)


        * [`Outcar.read_internal_strain_tensor()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_internal_strain_tensor)


        * [`Outcar.read_lcalcpol()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_lcalcpol)


        * [`Outcar.read_lepsilon()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_lepsilon)


        * [`Outcar.read_lepsilon_ionic()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_lepsilon_ionic)


        * [`Outcar.read_neb()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_neb)


        * [`Outcar.read_nmr_efg()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_nmr_efg)


        * [`Outcar.read_nmr_efg_tensor()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_nmr_efg_tensor)


        * [`Outcar.read_onsite_density_matrices()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_onsite_density_matrices)


        * [`Outcar.read_pattern()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_pattern)


        * [`Outcar.read_piezo_tensor()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_piezo_tensor)


        * [`Outcar.read_pseudo_zval()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_pseudo_zval)


        * [`Outcar.read_table_pattern()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Outcar.read_table_pattern)


    * [`Procar`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Procar)


        * [`Procar.data`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Procar.data)


        * [`Procar.weights`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Procar.weights)


        * [`Procar.get_occupation()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Procar.get_occupation)


        * [`Procar.get_projection_on_elements()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Procar.get_projection_on_elements)


    * [`UnconvergedVASPWarning`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.UnconvergedVASPWarning)


    * [`VaspParseError`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.VaspParseError)


    * [`Vasprun`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun)


        * [`Vasprun.ionic_steps`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.ionic_steps)


        * [`Vasprun.tdos`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.tdos)


        * [`Vasprun.idos`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.idos)


        * [`Vasprun.pdos`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.pdos)


        * [`Vasprun.efermi`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.efermi)


        * [`Vasprun.eigenvalues`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.eigenvalues)


        * [`Vasprun.projected_eigenvalues`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.projected_eigenvalues)


        * [`Vasprun.projected_magnetisation`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.projected_magnetisation)


        * [`Vasprun.other_dielectric`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.other_dielectric)


        * [`Vasprun.nionic_steps`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.nionic_steps)


        * [`Vasprun.force_constants`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.force_constants)


        * [`Vasprun.normalmode_eigenvals`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.normalmode_eigenvals)


        * [`Vasprun.normalmode_eigenvecs`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.normalmode_eigenvecs)


        * [`Vasprun.md_data`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.md_data)


        * [`Vasprun.incar`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.incar)


        * [`Vasprun.parameters`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.parameters)


        * [`Vasprun.kpoints`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.kpoints)


        * [`Vasprun.actual_kpoints`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.actual_kpoints)


        * [`Vasprun.actual_kpoints_weights`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.actual_kpoints_weights)


        * [`Vasprun.atomic_symbols`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.atomic_symbols)


        * [`Vasprun.potcar_symbols`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.potcar_symbols)


        * [`Vasprun.as_dict()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.as_dict)


        * [`Vasprun.calculate_efermi()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.calculate_efermi)


        * [`Vasprun.complete_dos`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.complete_dos)


        * [`Vasprun.complete_dos_normalized`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.complete_dos_normalized)


        * [`Vasprun.converged`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.converged)


        * [`Vasprun.converged_electronic`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.converged_electronic)


        * [`Vasprun.converged_ionic`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.converged_ionic)


        * [`Vasprun.dielectric`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.dielectric)


        * [`Vasprun.eigenvalue_band_properties`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.eigenvalue_band_properties)


        * [`Vasprun.epsilon_ionic`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.epsilon_ionic)


        * [`Vasprun.epsilon_static`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.epsilon_static)


        * [`Vasprun.epsilon_static_wolfe`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.epsilon_static_wolfe)


        * [`Vasprun.final_energy`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.final_energy)


        * [`Vasprun.get_band_structure()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.get_band_structure)


        * [`Vasprun.get_computed_entry()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.get_computed_entry)


        * [`Vasprun.get_potcars()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.get_potcars)


        * [`Vasprun.get_trajectory()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.get_trajectory)


        * [`Vasprun.hubbards`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.hubbards)


        * [`Vasprun.is_hubbard`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.is_hubbard)


        * [`Vasprun.is_spin`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.is_spin)


        * [`Vasprun.optical_absorption_coeff`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.optical_absorption_coeff)


        * [`Vasprun.run_type`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.run_type)


        * [`Vasprun.structures`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.structures)


        * [`Vasprun.update_charge_from_potcar()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.update_charge_from_potcar)


        * [`Vasprun.update_potcar_spec()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Vasprun.update_potcar_spec)


    * [`VolumetricData`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.VolumetricData)


        * [`VolumetricData.parse_file()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.VolumetricData.parse_file)


        * [`VolumetricData.write_file()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.VolumetricData.write_file)


    * [`WSWQ`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.WSWQ)


        * [`WSWQ.nspin`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.WSWQ.nspin)


        * [`WSWQ.nkpoints`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.WSWQ.nkpoints)


        * [`WSWQ.nbands`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.WSWQ.nbands)


        * [`WSWQ.me_real`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.WSWQ.me_real)


        * [`WSWQ.me_imag`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.WSWQ.me_imag)


        * [`WSWQ.data`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.WSWQ.data)


        * [`WSWQ.from_file()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.WSWQ.from_file)


        * [`WSWQ.me_imag`](pymatgen.io.vasp.outputs.md#id0)


        * [`WSWQ.me_real`](pymatgen.io.vasp.outputs.md#id1)


        * [`WSWQ.nbands`](pymatgen.io.vasp.outputs.md#id2)


        * [`WSWQ.nkpoints`](pymatgen.io.vasp.outputs.md#id3)


        * [`WSWQ.nspin`](pymatgen.io.vasp.outputs.md#id4)


    * [`Wavecar`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar)


        * [`Wavecar.filename`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.filename)


        * [`Wavecar.vasp_type`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.vasp_type)


        * [`Wavecar.nk`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.nk)


        * [`Wavecar.nb`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.nb)


        * [`Wavecar.encut`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.encut)


        * [`Wavecar.efermi`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.efermi)


        * [`Wavecar.a`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.a)


        * [`Wavecar.b`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.b)


        * [`Wavecar.vol`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.vol)


        * [`Wavecar.kpoints`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.kpoints)


        * [`Wavecar.band_energy`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.band_energy)


        * [`Wavecar.Gpoints`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.Gpoints)


        * [`Wavecar.coeffs`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.coeffs)


        * [`Wavecar.evaluate_wavefunc()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.evaluate_wavefunc)


        * [`Wavecar.fft_mesh()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.fft_mesh)


        * [`Wavecar.get_parchg()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.get_parchg)


        * [`Wavecar.write_unks()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Wavecar.write_unks)


    * [`Waveder`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder)


        * [`Waveder.cder_real`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder.cder_real)


        * [`Waveder.cder_imag`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder.cder_imag)


        * [`Waveder.cder`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder.cder)


        * [`Waveder.cder_imag`](pymatgen.io.vasp.outputs.md#id5)


        * [`Waveder.cder_real`](pymatgen.io.vasp.outputs.md#id6)


        * [`Waveder.from_binary()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder.from_binary)


        * [`Waveder.from_formatted()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder.from_formatted)


        * [`Waveder.get_orbital_derivative_between_states()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder.get_orbital_derivative_between_states)


        * [`Waveder.nbands`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder.nbands)


        * [`Waveder.nkpoints`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder.nkpoints)


        * [`Waveder.nspin`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Waveder.nspin)


    * [`Xdatcar`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Xdatcar)


        * [`Xdatcar.structures`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Xdatcar.structures)


        * [`Xdatcar.comment`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Xdatcar.comment)


        * [`Xdatcar.concatenate()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Xdatcar.concatenate)


        * [`Xdatcar.get_string()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Xdatcar.get_string)


        * [`Xdatcar.natoms`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Xdatcar.natoms)


        * [`Xdatcar.site_symbols`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Xdatcar.site_symbols)


        * [`Xdatcar.write_file()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.Xdatcar.write_file)


    * [`get_adjusted_fermi_level()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.get_adjusted_fermi_level)


    * [`get_band_structure_from_vasp_multiple_branches()`](pymatgen.io.vasp.outputs.md#pymatgen.io.vasp.outputs.get_band_structure_from_vasp_multiple_branches)


* [pymatgen.io.vasp.sets module](pymatgen.io.vasp.sets.md)


    * [`BadInputSetWarning`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.BadInputSetWarning)


    * [`DictSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.DictSet)


        * [`DictSet.calculate_ng()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.DictSet.calculate_ng)


        * [`DictSet.estimate_nbands()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.DictSet.estimate_nbands)


        * [`DictSet.incar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.DictSet.incar)


        * [`DictSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.DictSet.kpoints)


        * [`DictSet.nelect`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.DictSet.nelect)


        * [`DictSet.poscar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.DictSet.poscar)


        * [`DictSet.potcar_functional`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.DictSet.potcar_functional)


        * [`DictSet.structure`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.DictSet.structure)


        * [`DictSet.write_input()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.DictSet.write_input)


    * [`LobsterSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.LobsterSet)


        * [`LobsterSet.CONFIG`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.LobsterSet.CONFIG)


    * [`MITMDSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MITMDSet)


        * [`MITMDSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MITMDSet.kpoints)


    * [`MITNEBSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MITNEBSet)


        * [`MITNEBSet.poscar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MITNEBSet.poscar)


        * [`MITNEBSet.poscars`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MITNEBSet.poscars)


        * [`MITNEBSet.write_input()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MITNEBSet.write_input)


    * [`MITRelaxSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MITRelaxSet)


        * [`MITRelaxSet.CONFIG`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MITRelaxSet.CONFIG)


    * [`MPAbsorptionSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPAbsorptionSet)


        * [`MPAbsorptionSet.SUPPORTED_MODES`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPAbsorptionSet.SUPPORTED_MODES)


        * [`MPAbsorptionSet.from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPAbsorptionSet.from_prev_calc)


        * [`MPAbsorptionSet.incar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPAbsorptionSet.incar)


        * [`MPAbsorptionSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPAbsorptionSet.kpoints)


        * [`MPAbsorptionSet.override_from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPAbsorptionSet.override_from_prev_calc)


    * [`MPHSEBSSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPHSEBSSet)


        * [`MPHSEBSSet.from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPHSEBSSet.from_prev_calc)


        * [`MPHSEBSSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPHSEBSSet.kpoints)


        * [`MPHSEBSSet.override_from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPHSEBSSet.override_from_prev_calc)


    * [`MPHSERelaxSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPHSERelaxSet)


        * [`MPHSERelaxSet.CONFIG`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPHSERelaxSet.CONFIG)


    * [`MPMDSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPMDSet)


        * [`MPMDSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPMDSet.kpoints)


    * [`MPMetalRelaxSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPMetalRelaxSet)


        * [`MPMetalRelaxSet.CONFIG`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPMetalRelaxSet.CONFIG)


    * [`MPNMRSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPNMRSet)


        * [`MPNMRSet.incar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPNMRSet.incar)


    * [`MPNonSCFSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPNonSCFSet)


        * [`MPNonSCFSet.from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPNonSCFSet.from_prev_calc)


        * [`MPNonSCFSet.incar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPNonSCFSet.incar)


        * [`MPNonSCFSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPNonSCFSet.kpoints)


        * [`MPNonSCFSet.override_from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPNonSCFSet.override_from_prev_calc)


    * [`MPRelaxSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPRelaxSet)


        * [`MPRelaxSet.CONFIG`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPRelaxSet.CONFIG)


    * [`MPSOCSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPSOCSet)


        * [`MPSOCSet.from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPSOCSet.from_prev_calc)


        * [`MPSOCSet.incar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPSOCSet.incar)


        * [`MPSOCSet.override_from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPSOCSet.override_from_prev_calc)


        * [`MPSOCSet.user_potcar_functional`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPSOCSet.user_potcar_functional)


    * [`MPScanRelaxSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPScanRelaxSet)


        * [`MPScanRelaxSet.CONFIG`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPScanRelaxSet.CONFIG)


    * [`MPScanStaticSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPScanStaticSet)


        * [`MPScanStaticSet.from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPScanStaticSet.from_prev_calc)


        * [`MPScanStaticSet.incar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPScanStaticSet.incar)


        * [`MPScanStaticSet.override_from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPScanStaticSet.override_from_prev_calc)


        * [`MPScanStaticSet.user_potcar_functional`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPScanStaticSet.user_potcar_functional)


    * [`MPStaticSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPStaticSet)


        * [`MPStaticSet.from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPStaticSet.from_prev_calc)


        * [`MPStaticSet.incar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPStaticSet.incar)


        * [`MPStaticSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPStaticSet.kpoints)


        * [`MPStaticSet.override_from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPStaticSet.override_from_prev_calc)


        * [`MPStaticSet.user_potcar_functional`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MPStaticSet.user_potcar_functional)


    * [`MVLElasticSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLElasticSet)


        * [`MVLElasticSet.user_potcar_functional`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLElasticSet.user_potcar_functional)


    * [`MVLGBSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGBSet)


        * [`MVLGBSet.incar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGBSet.incar)


        * [`MVLGBSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGBSet.kpoints)


        * [`MVLGBSet.user_potcar_functional`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGBSet.user_potcar_functional)


    * [`MVLGWSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGWSet)


        * [`MVLGWSet.CONFIG`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGWSet.CONFIG)


        * [`MVLGWSet.SUPPORTED_MODES`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGWSet.SUPPORTED_MODES)


        * [`MVLGWSet.from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGWSet.from_prev_calc)


        * [`MVLGWSet.incar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGWSet.incar)


        * [`MVLGWSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGWSet.kpoints)


        * [`MVLGWSet.override_from_prev_calc()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLGWSet.override_from_prev_calc)


    * [`MVLNPTMDSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLNPTMDSet)


        * [`MVLNPTMDSet.user_potcar_functional`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLNPTMDSet.user_potcar_functional)


    * [`MVLRelax52Set`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLRelax52Set)


        * [`MVLRelax52Set.CONFIG`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLRelax52Set.CONFIG)


    * [`MVLScanRelaxSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLScanRelaxSet)


        * [`MVLScanRelaxSet.user_potcar_functional`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLScanRelaxSet.user_potcar_functional)


    * [`MVLSlabSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLSlabSet)


        * [`MVLSlabSet.as_dict()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLSlabSet.as_dict)


        * [`MVLSlabSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLSlabSet.kpoints)


        * [`MVLSlabSet.user_potcar_functional`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.MVLSlabSet.user_potcar_functional)


    * [`VaspInputSet`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet)


        * [`VaspInputSet.as_dict()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet.as_dict)


        * [`VaspInputSet.get_vasp_input()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet.get_vasp_input)


        * [`VaspInputSet.incar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet.incar)


        * [`VaspInputSet.kpoints`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet.kpoints)


        * [`VaspInputSet.poscar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet.poscar)


        * [`VaspInputSet.potcar`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet.potcar)


        * [`VaspInputSet.potcar_symbols`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet.potcar_symbols)


        * [`VaspInputSet.write_input()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.VaspInputSet.write_input)


    * [`batch_write_input()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.batch_write_input)


    * [`get_structure_from_prev_run()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.get_structure_from_prev_run)


    * [`get_valid_magmom_struct()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.get_valid_magmom_struct)


    * [`get_vasprun_outcar()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.get_vasprun_outcar)


    * [`next_num_with_prime_factors()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.next_num_with_prime_factors)


    * [`primes_less_than()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.primes_less_than)


    * [`standardize_structure()`](pymatgen.io.vasp.sets.md#pymatgen.io.vasp.sets.standardize_structure)