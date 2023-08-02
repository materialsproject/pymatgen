---
layout: default
title: pymatgen.command_line.vampire_caller.md
nav_exclude: true
---

# pymatgen.command_line.vampire_caller module

This module implements an interface to the VAMPIRE code for atomistic
simulations of magnetic materials.

This module depends on a compiled vampire executable available in the path.
Please download at [https://vampire.york.ac.uk/download/](https://vampire.york.ac.uk/download/) and
follow the instructions to compile the executable.

If you use this module, please cite:

“Atomistic spin model simulations of magnetic nanomaterials.”
R. F. L. Evans, W. J. Fan, P. Chureemart, T. A. Ostler, M. O. A. Ellis
and R. W. Chantrell. J. Phys.: Condens. Matter 26, 103202 (2014)


### _class_ pymatgen.command_line.vampire_caller.VampireCaller(ordered_structures=None, energies=None, mc_box_size=4.0, equil_timesteps=2000, mc_timesteps=4000, save_inputs=False, hm=None, avg=True, user_input_settings=None)
Bases: `object`

Run Vampire on a material with magnetic ordering and exchange parameter information to compute the critical
temperature with classical Monte Carlo.

user_input_settings is a dictionary that can contain:
\* start_t (int): Start MC sim at this temp, defaults to 0 K.
\* end_t (int): End MC sim at this temp, defaults to 1500 K.
\* temp_increment (int): Temp step size, defaults to 25 K.


* **Parameters**


    * **ordered_structures** (*list*) – Structure objects with magmoms.


    * **energies** (*list*) – Energies of each relaxed magnetic structure.


    * **mc_box_size** (*float*) – x=y=z dimensions (nm) of MC simulation box


    * **equil_timesteps** (*int*) – number of MC steps for equilibrating


    * **mc_timesteps** (*int*) – number of MC steps for averaging


    * **save_inputs** (*bool*) – if True, save scratch dir of vampire input files


    * **hm** ([*HeisenbergModel*](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel)) – object already fit to low energy
    magnetic orderings.


    * **avg** (*bool*) – If True, simply use <J> exchange parameter estimate.
    If False, attempt to use NN, NNN, etc. interactions.


    * **user_input_settings** (*dict*) – optional commands for VAMPIRE Monte Carlo


    * **sgraph** ([*StructureGraph*](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph)) – Ground state graph.


    * **unique_site_ids** (*dict*) – Maps each site to its unique identifier


    * **nn_interactions** (*dict*) – {i: j} pairs of NN interactions
    between unique sites.


    * **ex_params** (*dict*) – Exchange parameter values (meV/atom)


    * **mft_t** (*float*) – Mean field theory estimate of critical T


    * **mat_name** (*str*) – Formula unit label for input files


    * **mat_id_dict** (*dict*) – Maps sites to material id # for vampire
    indexing.



#### _static_ parse_stdout(vamp_stdout, n_mats: int)
Parse stdout from Vampire.


* **Parameters**


    * **vamp_stdout** (*txt file*) – Vampire ‘output’ file.


    * **n_mats** (*int*) – Number of materials in Vampire simulation.



* **Returns**

    MSONable vampire output.
    critical_temp (float): Calculated critical temp.



* **Return type**

    parsed_out (DataFrame)



### _class_ pymatgen.command_line.vampire_caller.VampireOutput(parsed_out=None, nmats=None, critical_temp=None)
Bases: `MSONable`

This class processes results from a Vampire Monte Carlo simulation
and returns the critical temperature.


* **Parameters**


    * **parsed_out** (*json*) – json rep of parsed stdout DataFrame.


    * **nmats** (*int*) – Number of distinct materials (1 for each specie and up/down spin).


    * **critical_temp** (*float*) – Monte Carlo Tc result.