---
layout: default
title: pymatgen.command_line.md
nav_exclude: true
---

# pymatgen.command_line package

This package contains various command line wrappers to programs used in
pymatgen that do not have Python equivalents.



* [pymatgen.command_line.bader_caller module](pymatgen.command_line.bader_caller.md)


    * [`BaderAnalysis`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis)


        * [`BaderAnalysis.data`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.data)


        * [`BaderAnalysis.vacuum_volume`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.vacuum_volume)


        * [`BaderAnalysis.vacuum_charge`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.vacuum_charge)


        * [`BaderAnalysis.nelectrons`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.nelectrons)


        * [`BaderAnalysis.chgcar`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.chgcar)


        * [`BaderAnalysis.atomic_densities`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.atomic_densities)


        * [`BaderAnalysis.from_path()`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.from_path)


        * [`BaderAnalysis.get_charge()`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.get_charge)


        * [`BaderAnalysis.get_charge_decorated_structure()`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.get_charge_decorated_structure)


        * [`BaderAnalysis.get_charge_transfer()`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.get_charge_transfer)


        * [`BaderAnalysis.get_decorated_structure()`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.get_decorated_structure)


        * [`BaderAnalysis.get_oxidation_state_decorated_structure()`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.get_oxidation_state_decorated_structure)


        * [`BaderAnalysis.get_partial_charge()`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.get_partial_charge)


        * [`BaderAnalysis.summary`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.BaderAnalysis.summary)


    * [`bader_analysis_from_objects()`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.bader_analysis_from_objects)


    * [`bader_analysis_from_path()`](pymatgen.command_line.bader_caller.md#pymatgen.command_line.bader_caller.bader_analysis_from_path)


* [pymatgen.command_line.chargemol_caller module](pymatgen.command_line.chargemol_caller.md)


    * [`ChargemolAnalysis`](pymatgen.command_line.chargemol_caller.md#pymatgen.command_line.chargemol_caller.ChargemolAnalysis)


        * [`ChargemolAnalysis.get_bond_order()`](pymatgen.command_line.chargemol_caller.md#pymatgen.command_line.chargemol_caller.ChargemolAnalysis.get_bond_order)


        * [`ChargemolAnalysis.get_charge()`](pymatgen.command_line.chargemol_caller.md#pymatgen.command_line.chargemol_caller.ChargemolAnalysis.get_charge)


        * [`ChargemolAnalysis.get_charge_transfer()`](pymatgen.command_line.chargemol_caller.md#pymatgen.command_line.chargemol_caller.ChargemolAnalysis.get_charge_transfer)


        * [`ChargemolAnalysis.get_partial_charge()`](pymatgen.command_line.chargemol_caller.md#pymatgen.command_line.chargemol_caller.ChargemolAnalysis.get_partial_charge)


        * [`ChargemolAnalysis.get_property_decorated_structure()`](pymatgen.command_line.chargemol_caller.md#pymatgen.command_line.chargemol_caller.ChargemolAnalysis.get_property_decorated_structure)


        * [`ChargemolAnalysis.summary`](pymatgen.command_line.chargemol_caller.md#pymatgen.command_line.chargemol_caller.ChargemolAnalysis.summary)


* [pymatgen.command_line.critic2_caller module](pymatgen.command_line.critic2_caller.md)


    * [`Critic2Analysis`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.Critic2Analysis)


        * [`Critic2Analysis.get_critical_point_for_site()`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.Critic2Analysis.get_critical_point_for_site)


        * [`Critic2Analysis.get_volume_and_charge_for_site()`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.Critic2Analysis.get_volume_and_charge_for_site)


        * [`Critic2Analysis.structure_graph()`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.Critic2Analysis.structure_graph)


    * [`Critic2Caller`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.Critic2Caller)


        * [`Critic2Caller.from_chgcar()`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.Critic2Caller.from_chgcar)


        * [`Critic2Caller.from_path()`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.Critic2Caller.from_path)


    * [`CriticalPoint`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.CriticalPoint)


        * [`CriticalPoint.ellipticity`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.CriticalPoint.ellipticity)


        * [`CriticalPoint.laplacian`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.CriticalPoint.laplacian)


        * [`CriticalPoint.type`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.CriticalPoint.type)


    * [`CriticalPointType`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.CriticalPointType)


        * [`CriticalPointType.bond`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.CriticalPointType.bond)


        * [`CriticalPointType.cage`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.CriticalPointType.cage)


        * [`CriticalPointType.nnattr`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.CriticalPointType.nnattr)


        * [`CriticalPointType.nucleus`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.CriticalPointType.nucleus)


        * [`CriticalPointType.ring`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.CriticalPointType.ring)


    * [`get_filepath()`](pymatgen.command_line.critic2_caller.md#pymatgen.command_line.critic2_caller.get_filepath)


* [pymatgen.command_line.enumlib_caller module](pymatgen.command_line.enumlib_caller.md)


    * [`EnumError`](pymatgen.command_line.enumlib_caller.md#pymatgen.command_line.enumlib_caller.EnumError)


* [pymatgen.command_line.gulp_caller module](pymatgen.command_line.gulp_caller.md)


    * [`BuckinghamPotential`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.BuckinghamPotential)


    * [`GulpCaller`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpCaller)


        * [`GulpCaller.run()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpCaller.run)


    * [`GulpConvergenceError`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpConvergenceError)


    * [`GulpError`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpError)


    * [`GulpIO`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO)


        * [`GulpIO.buckingham_input()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO.buckingham_input)


        * [`GulpIO.buckingham_potential()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO.buckingham_potential)


        * [`GulpIO.get_energy()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO.get_energy)


        * [`GulpIO.get_relaxed_structure()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO.get_relaxed_structure)


        * [`GulpIO.keyword_line()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO.keyword_line)


        * [`GulpIO.library_line()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO.library_line)


        * [`GulpIO.specie_potential_lines()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO.specie_potential_lines)


        * [`GulpIO.structure_lines()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO.structure_lines)


        * [`GulpIO.tersoff_input()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO.tersoff_input)


        * [`GulpIO.tersoff_potential()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.GulpIO.tersoff_potential)


    * [`TersoffPotential`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.TersoffPotential)


    * [`get_energy_buckingham()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.get_energy_buckingham)


    * [`get_energy_relax_structure_buckingham()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.get_energy_relax_structure_buckingham)


    * [`get_energy_tersoff()`](pymatgen.command_line.gulp_caller.md#pymatgen.command_line.gulp_caller.get_energy_tersoff)


* [pymatgen.command_line.mcsqs_caller module](pymatgen.command_line.mcsqs_caller.md)


    * [`Sqs`](pymatgen.command_line.mcsqs_caller.md#pymatgen.command_line.mcsqs_caller.Sqs)


        * [`Sqs.allsqs`](pymatgen.command_line.mcsqs_caller.md#pymatgen.command_line.mcsqs_caller.Sqs.allsqs)


        * [`Sqs.bestsqs`](pymatgen.command_line.mcsqs_caller.md#pymatgen.command_line.mcsqs_caller.Sqs.bestsqs)


        * [`Sqs.clusters`](pymatgen.command_line.mcsqs_caller.md#pymatgen.command_line.mcsqs_caller.Sqs.clusters)


        * [`Sqs.directory`](pymatgen.command_line.mcsqs_caller.md#pymatgen.command_line.mcsqs_caller.Sqs.directory)


        * [`Sqs.objective_function`](pymatgen.command_line.mcsqs_caller.md#pymatgen.command_line.mcsqs_caller.Sqs.objective_function)


    * [`run_mcsqs()`](pymatgen.command_line.mcsqs_caller.md#pymatgen.command_line.mcsqs_caller.run_mcsqs)


* [pymatgen.command_line.vampire_caller module](pymatgen.command_line.vampire_caller.md)


    * [`VampireCaller`](pymatgen.command_line.vampire_caller.md#pymatgen.command_line.vampire_caller.VampireCaller)


        * [`VampireCaller.parse_stdout()`](pymatgen.command_line.vampire_caller.md#pymatgen.command_line.vampire_caller.VampireCaller.parse_stdout)


    * [`VampireOutput`](pymatgen.command_line.vampire_caller.md#pymatgen.command_line.vampire_caller.VampireOutput)