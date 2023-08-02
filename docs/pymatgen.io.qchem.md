---
layout: default
title: pymatgen.io.qchem.md
nav_exclude: true
---

# pymatgen.io.qchem package

This package implements modules for input and output to and from Qchem.



* [pymatgen.io.qchem.inputs module](pymatgen.io.qchem.inputs.md)


    * [`QCInput`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput)


        * [`QCInput.almo_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.almo_template)


        * [`QCInput.cdft_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.cdft_template)


        * [`QCInput.find_sections()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.find_sections)


        * [`QCInput.from_file()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.from_file)


        * [`QCInput.from_multi_jobs_file()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.from_multi_jobs_file)


        * [`QCInput.from_str()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.from_str)


        * [`QCInput.geom_opt_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.geom_opt_template)


        * [`QCInput.get_string()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.get_string)


        * [`QCInput.molecule_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.molecule_template)


        * [`QCInput.multi_job_string()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.multi_job_string)


        * [`QCInput.nbo_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.nbo_template)


        * [`QCInput.opt_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.opt_template)


        * [`QCInput.pcm_nonels_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.pcm_nonels_template)


        * [`QCInput.pcm_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.pcm_template)


        * [`QCInput.plots_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.plots_template)


        * [`QCInput.read_almo()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_almo)


        * [`QCInput.read_cdft()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_cdft)


        * [`QCInput.read_geom_opt()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_geom_opt)


        * [`QCInput.read_molecule()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_molecule)


        * [`QCInput.read_nbo()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_nbo)


        * [`QCInput.read_opt()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_opt)


        * [`QCInput.read_pcm()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_pcm)


        * [`QCInput.read_pcm_nonels()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_pcm_nonels)


        * [`QCInput.read_plots()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_plots)


        * [`QCInput.read_rem()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_rem)


        * [`QCInput.read_scan()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_scan)


        * [`QCInput.read_smx()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_smx)


        * [`QCInput.read_solvent()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_solvent)


        * [`QCInput.read_svp()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_svp)


        * [`QCInput.read_vdw()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.read_vdw)


        * [`QCInput.rem_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.rem_template)


        * [`QCInput.scan_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.scan_template)


        * [`QCInput.smx_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.smx_template)


        * [`QCInput.solvent_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.solvent_template)


        * [`QCInput.svp_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.svp_template)


        * [`QCInput.van_der_waals_template()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.van_der_waals_template)


        * [`QCInput.write_multi_job_file()`](pymatgen.io.qchem.inputs.md#pymatgen.io.qchem.inputs.QCInput.write_multi_job_file)


* [pymatgen.io.qchem.outputs module](pymatgen.io.qchem.outputs.md)


    * [`QCOutput`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.QCOutput)


        * [`QCOutput.as_dict()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.QCOutput.as_dict)


        * [`QCOutput.multiple_outputs_from_file()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.QCOutput.multiple_outputs_from_file)


    * [`check_for_structure_changes()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.check_for_structure_changes)


    * [`get_percentage()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.get_percentage)


    * [`jump_to_header()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.jump_to_header)


    * [`nbo_parser()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.nbo_parser)


    * [`parse_hybridization_character()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.parse_hybridization_character)


    * [`parse_hyperbonds()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.parse_hyperbonds)


    * [`parse_natural_populations()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.parse_natural_populations)


    * [`parse_perturbation_energy()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.parse_perturbation_energy)


    * [`z_int()`](pymatgen.io.qchem.outputs.md#pymatgen.io.qchem.outputs.z_int)


* [pymatgen.io.qchem.sets module](pymatgen.io.qchem.sets.md)


    * [`ForceSet`](pymatgen.io.qchem.sets.md#pymatgen.io.qchem.sets.ForceSet)


    * [`FreqSet`](pymatgen.io.qchem.sets.md#pymatgen.io.qchem.sets.FreqSet)


    * [`OptSet`](pymatgen.io.qchem.sets.md#pymatgen.io.qchem.sets.OptSet)


    * [`PESScanSet`](pymatgen.io.qchem.sets.md#pymatgen.io.qchem.sets.PESScanSet)


    * [`QChemDictSet`](pymatgen.io.qchem.sets.md#pymatgen.io.qchem.sets.QChemDictSet)


        * [`QChemDictSet.write()`](pymatgen.io.qchem.sets.md#pymatgen.io.qchem.sets.QChemDictSet.write)


    * [`SinglePointSet`](pymatgen.io.qchem.sets.md#pymatgen.io.qchem.sets.SinglePointSet)


    * [`TransitionStateSet`](pymatgen.io.qchem.sets.md#pymatgen.io.qchem.sets.TransitionStateSet)


* [pymatgen.io.qchem.utils module](pymatgen.io.qchem.utils.md)


    * [`lower_and_check_unique()`](pymatgen.io.qchem.utils.md#pymatgen.io.qchem.utils.lower_and_check_unique)


    * [`process_parsed_HESS()`](pymatgen.io.qchem.utils.md#pymatgen.io.qchem.utils.process_parsed_HESS)


    * [`process_parsed_coords()`](pymatgen.io.qchem.utils.md#pymatgen.io.qchem.utils.process_parsed_coords)


    * [`process_parsed_fock_matrix()`](pymatgen.io.qchem.utils.md#pymatgen.io.qchem.utils.process_parsed_fock_matrix)


    * [`read_matrix_pattern()`](pymatgen.io.qchem.utils.md#pymatgen.io.qchem.utils.read_matrix_pattern)


    * [`read_pattern()`](pymatgen.io.qchem.utils.md#pymatgen.io.qchem.utils.read_pattern)


    * [`read_table_pattern()`](pymatgen.io.qchem.utils.md#pymatgen.io.qchem.utils.read_table_pattern)