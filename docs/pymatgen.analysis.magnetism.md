---
layout: default
title: pymatgen.analysis.magnetism.md
nav_exclude: true
---

# pymatgen.analysis.magnetism package

Package for analysis of magnetic structures.



* [pymatgen.analysis.magnetism.analyzer module](pymatgen.analysis.magnetism.analyzer.md)


    * [`CollinearMagneticStructureAnalyzer`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer)


        * [`CollinearMagneticStructureAnalyzer.get_exchange_group_info()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_exchange_group_info)


        * [`CollinearMagneticStructureAnalyzer.get_ferromagnetic_structure()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_ferromagnetic_structure)


        * [`CollinearMagneticStructureAnalyzer.get_nonmagnetic_structure()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_nonmagnetic_structure)


        * [`CollinearMagneticStructureAnalyzer.get_structure_with_only_magnetic_atoms()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_structure_with_only_magnetic_atoms)


        * [`CollinearMagneticStructureAnalyzer.get_structure_with_spin()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_structure_with_spin)


        * [`CollinearMagneticStructureAnalyzer.is_magnetic`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.is_magnetic)


        * [`CollinearMagneticStructureAnalyzer.magmoms`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.magmoms)


        * [`CollinearMagneticStructureAnalyzer.magnetic_species_and_magmoms`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.magnetic_species_and_magmoms)


        * [`CollinearMagneticStructureAnalyzer.matches_ordering()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.matches_ordering)


        * [`CollinearMagneticStructureAnalyzer.number_of_magnetic_sites`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.number_of_magnetic_sites)


        * [`CollinearMagneticStructureAnalyzer.number_of_unique_magnetic_sites()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.number_of_unique_magnetic_sites)


        * [`CollinearMagneticStructureAnalyzer.ordering`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.ordering)


        * [`CollinearMagneticStructureAnalyzer.types_of_magnetic_specie`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.types_of_magnetic_specie)


        * [`CollinearMagneticStructureAnalyzer.types_of_magnetic_species`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.types_of_magnetic_species)


    * [`MagneticDeformation`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation)


        * [`MagneticDeformation.deformation`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation.deformation)


        * [`MagneticDeformation.type`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation.type)


    * [`MagneticStructureEnumerator`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator)


        * [`MagneticStructureEnumerator.available_strategies`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator.available_strategies)


    * [`Ordering`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering)


        * [`Ordering.AFM`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering.AFM)


        * [`Ordering.FM`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering.FM)


        * [`Ordering.FiM`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering.FiM)


        * [`Ordering.NM`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering.NM)


        * [`Ordering.Unknown`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering.Unknown)


    * [`OverwriteMagmomMode`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode)


        * [`OverwriteMagmomMode.none`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.none)


        * [`OverwriteMagmomMode.normalize`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.normalize)


        * [`OverwriteMagmomMode.replace_all`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.replace_all)


        * [`OverwriteMagmomMode.respect_sign`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.respect_sign)


        * [`OverwriteMagmomMode.respect_zero`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.respect_zero)


    * [`magnetic_deformation()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.magnetic_deformation)


* [pymatgen.analysis.magnetism.heisenberg module](pymatgen.analysis.magnetism.heisenberg.md)


    * [`HeisenbergMapper`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper)


        * [`HeisenbergMapper.estimate_exchange()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.estimate_exchange)


        * [`HeisenbergMapper.get_exchange()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_exchange)


        * [`HeisenbergMapper.get_heisenberg_model()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_heisenberg_model)


        * [`HeisenbergMapper.get_interaction_graph()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_interaction_graph)


        * [`HeisenbergMapper.get_low_energy_orderings()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_low_energy_orderings)


        * [`HeisenbergMapper.get_mft_temperature()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_mft_temperature)


    * [`HeisenbergModel`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel)


        * [`HeisenbergModel.as_dict()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel.as_dict)


        * [`HeisenbergModel.from_dict()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel.from_dict)


    * [`HeisenbergScreener`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener)


        * [`HeisenbergScreener.screened_structures`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener.screened_structures)


        * [`HeisenbergScreener.screened_energies`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener.screened_energies)


* [pymatgen.analysis.magnetism.jahnteller module](pymatgen.analysis.magnetism.jahnteller.md)


    * [`JahnTellerAnalyzer`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer)


        * [`JahnTellerAnalyzer.get_analysis()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_analysis)


        * [`JahnTellerAnalyzer.get_analysis_and_structure()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_analysis_and_structure)


        * [`JahnTellerAnalyzer.get_magnitude_of_effect_from_species()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_magnitude_of_effect_from_species)


        * [`JahnTellerAnalyzer.get_magnitude_of_effect_from_spin_config()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_magnitude_of_effect_from_spin_config)


        * [`JahnTellerAnalyzer.is_jahn_teller_active()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.is_jahn_teller_active)


        * [`JahnTellerAnalyzer.mu_so()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.mu_so)


        * [`JahnTellerAnalyzer.tag_structure()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.tag_structure)