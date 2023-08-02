---
layout: default
title: pymatgen.phonon.md
nav_exclude: true
---

# pymatgen.phonon package

Phonon DOS and bandstructure analysis package.



* [pymatgen.phonon.bandstructure module](pymatgen.phonon.bandstructure.md)


    * [`PhononBandStructure`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure)


        * [`PhononBandStructure.as_dict()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure.as_dict)


        * [`PhononBandStructure.asr_breaking()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure.asr_breaking)


        * [`PhononBandStructure.from_dict()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure.from_dict)


        * [`PhononBandStructure.get_nac_eigendisplacements_along_dir()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure.get_nac_eigendisplacements_along_dir)


        * [`PhononBandStructure.get_nac_frequencies_along_dir()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure.get_nac_frequencies_along_dir)


        * [`PhononBandStructure.has_eigendisplacements`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure.has_eigendisplacements)


        * [`PhononBandStructure.has_imaginary_freq()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure.has_imaginary_freq)


        * [`PhononBandStructure.has_nac`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure.has_nac)


        * [`PhononBandStructure.min_freq()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructure.min_freq)


    * [`PhononBandStructureSymmLine`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructureSymmLine)


        * [`PhononBandStructureSymmLine.as_dict()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructureSymmLine.as_dict)


        * [`PhononBandStructureSymmLine.as_phononwebsite()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructureSymmLine.as_phononwebsite)


        * [`PhononBandStructureSymmLine.band_reorder()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructureSymmLine.band_reorder)


        * [`PhononBandStructureSymmLine.from_dict()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructureSymmLine.from_dict)


        * [`PhononBandStructureSymmLine.get_branch()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructureSymmLine.get_branch)


        * [`PhononBandStructureSymmLine.get_equivalent_qpoints()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructureSymmLine.get_equivalent_qpoints)


        * [`PhononBandStructureSymmLine.write_phononwebsite()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.PhononBandStructureSymmLine.write_phononwebsite)


    * [`eigenvectors_from_displacements()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.eigenvectors_from_displacements)


    * [`estimate_band_connection()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.estimate_band_connection)


    * [`get_reasonable_repetitions()`](pymatgen.phonon.bandstructure.md#pymatgen.phonon.bandstructure.get_reasonable_repetitions)


* [pymatgen.phonon.dos module](pymatgen.phonon.dos.md)


    * [`CompletePhononDos`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.CompletePhononDos)


        * [`CompletePhononDos.pdos`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.CompletePhononDos.pdos)


        * [`CompletePhononDos.as_dict()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.CompletePhononDos.as_dict)


        * [`CompletePhononDos.from_dict()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.CompletePhononDos.from_dict)


        * [`CompletePhononDos.get_element_dos()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.CompletePhononDos.get_element_dos)


        * [`CompletePhononDos.get_site_dos()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.CompletePhononDos.get_site_dos)


    * [`PhononDos`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos)


        * [`PhononDos.as_dict()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos.as_dict)


        * [`PhononDos.cv()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos.cv)


        * [`PhononDos.entropy()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos.entropy)


        * [`PhononDos.from_dict()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos.from_dict)


        * [`PhononDos.get_interpolated_value()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos.get_interpolated_value)


        * [`PhononDos.get_smeared_densities()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos.get_smeared_densities)


        * [`PhononDos.helmholtz_free_energy()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos.helmholtz_free_energy)


        * [`PhononDos.ind_zero_freq()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos.ind_zero_freq)


        * [`PhononDos.internal_energy()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos.internal_energy)


        * [`PhononDos.zero_point_energy()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.PhononDos.zero_point_energy)


    * [`coth()`](pymatgen.phonon.dos.md#pymatgen.phonon.dos.coth)


* [pymatgen.phonon.gruneisen module](pymatgen.phonon.gruneisen.md)


    * [`GruneisenParameter`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenParameter)


        * [`GruneisenParameter.acoustic_debye_temp`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenParameter.acoustic_debye_temp)


        * [`GruneisenParameter.average_gruneisen()`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenParameter.average_gruneisen)


        * [`GruneisenParameter.debye_temp_limit`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenParameter.debye_temp_limit)


        * [`GruneisenParameter.debye_temp_phonopy()`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenParameter.debye_temp_phonopy)


        * [`GruneisenParameter.phdos`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenParameter.phdos)


        * [`GruneisenParameter.tdos`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenParameter.tdos)


        * [`GruneisenParameter.thermal_conductivity_slack()`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenParameter.thermal_conductivity_slack)


    * [`GruneisenPhononBandStructure`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenPhononBandStructure)


        * [`GruneisenPhononBandStructure.as_dict()`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenPhononBandStructure.as_dict)


        * [`GruneisenPhononBandStructure.from_dict()`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenPhononBandStructure.from_dict)


    * [`GruneisenPhononBandStructureSymmLine`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenPhononBandStructureSymmLine)


        * [`GruneisenPhononBandStructureSymmLine.from_dict()`](pymatgen.phonon.gruneisen.md#pymatgen.phonon.gruneisen.GruneisenPhononBandStructureSymmLine.from_dict)


* [pymatgen.phonon.ir_spectra module](pymatgen.phonon.ir_spectra.md)


    * [`IRDielectricTensor`](pymatgen.phonon.ir_spectra.md#pymatgen.phonon.ir_spectra.IRDielectricTensor)


        * [`IRDielectricTensor.as_dict()`](pymatgen.phonon.ir_spectra.md#pymatgen.phonon.ir_spectra.IRDielectricTensor.as_dict)


        * [`IRDielectricTensor.from_dict()`](pymatgen.phonon.ir_spectra.md#pymatgen.phonon.ir_spectra.IRDielectricTensor.from_dict)


        * [`IRDielectricTensor.get_ir_spectra()`](pymatgen.phonon.ir_spectra.md#pymatgen.phonon.ir_spectra.IRDielectricTensor.get_ir_spectra)


        * [`IRDielectricTensor.get_plotter()`](pymatgen.phonon.ir_spectra.md#pymatgen.phonon.ir_spectra.IRDielectricTensor.get_plotter)


        * [`IRDielectricTensor.get_spectrum()`](pymatgen.phonon.ir_spectra.md#pymatgen.phonon.ir_spectra.IRDielectricTensor.get_spectrum)


        * [`IRDielectricTensor.max_phfreq`](pymatgen.phonon.ir_spectra.md#pymatgen.phonon.ir_spectra.IRDielectricTensor.max_phfreq)


        * [`IRDielectricTensor.nph_freqs`](pymatgen.phonon.ir_spectra.md#pymatgen.phonon.ir_spectra.IRDielectricTensor.nph_freqs)


        * [`IRDielectricTensor.plot()`](pymatgen.phonon.ir_spectra.md#pymatgen.phonon.ir_spectra.IRDielectricTensor.plot)


        * [`IRDielectricTensor.write_json()`](pymatgen.phonon.ir_spectra.md#pymatgen.phonon.ir_spectra.IRDielectricTensor.write_json)


* [pymatgen.phonon.plotter module](pymatgen.phonon.plotter.md)


    * [`FreqUnits`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.FreqUnits)


        * [`FreqUnits.factor`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.FreqUnits.factor)


        * [`FreqUnits.label`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.FreqUnits.label)


    * [`GruneisenPhononBSPlotter`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.GruneisenPhononBSPlotter)


        * [`GruneisenPhononBSPlotter.bs_plot_data()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.GruneisenPhononBSPlotter.bs_plot_data)


        * [`GruneisenPhononBSPlotter.get_plot_gs()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.GruneisenPhononBSPlotter.get_plot_gs)


        * [`GruneisenPhononBSPlotter.plot_compare_gs()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.GruneisenPhononBSPlotter.plot_compare_gs)


        * [`GruneisenPhononBSPlotter.save_plot_gs()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.GruneisenPhononBSPlotter.save_plot_gs)


        * [`GruneisenPhononBSPlotter.show_gs()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.GruneisenPhononBSPlotter.show_gs)


    * [`GruneisenPlotter`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.GruneisenPlotter)


        * [`GruneisenPlotter.get_plot()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.GruneisenPlotter.get_plot)


        * [`GruneisenPlotter.save_plot()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.GruneisenPlotter.save_plot)


        * [`GruneisenPlotter.show()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.GruneisenPlotter.show)


    * [`PhononBSPlotter`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononBSPlotter)


        * [`PhononBSPlotter.bs_plot_data()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononBSPlotter.bs_plot_data)


        * [`PhononBSPlotter.get_plot()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononBSPlotter.get_plot)


        * [`PhononBSPlotter.get_proj_plot()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononBSPlotter.get_proj_plot)


        * [`PhononBSPlotter.get_ticks()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononBSPlotter.get_ticks)


        * [`PhononBSPlotter.plot_brillouin()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononBSPlotter.plot_brillouin)


        * [`PhononBSPlotter.plot_compare()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononBSPlotter.plot_compare)


        * [`PhononBSPlotter.save_plot()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononBSPlotter.save_plot)


        * [`PhononBSPlotter.show()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononBSPlotter.show)


        * [`PhononBSPlotter.show_proj()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononBSPlotter.show_proj)


    * [`PhononDosPlotter`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononDosPlotter)


        * [`PhononDosPlotter.add_dos()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononDosPlotter.add_dos)


        * [`PhononDosPlotter.add_dos_dict()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononDosPlotter.add_dos_dict)


        * [`PhononDosPlotter.get_dos_dict()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononDosPlotter.get_dos_dict)


        * [`PhononDosPlotter.get_plot()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononDosPlotter.get_plot)


        * [`PhononDosPlotter.save_plot()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononDosPlotter.save_plot)


        * [`PhononDosPlotter.show()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.PhononDosPlotter.show)


    * [`ThermoPlotter`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.ThermoPlotter)


        * [`ThermoPlotter.plot_cv()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.ThermoPlotter.plot_cv)


        * [`ThermoPlotter.plot_entropy()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.ThermoPlotter.plot_entropy)


        * [`ThermoPlotter.plot_helmholtz_free_energy()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.ThermoPlotter.plot_helmholtz_free_energy)


        * [`ThermoPlotter.plot_internal_energy()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.ThermoPlotter.plot_internal_energy)


        * [`ThermoPlotter.plot_thermodynamic_properties()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.ThermoPlotter.plot_thermodynamic_properties)


    * [`freq_units()`](pymatgen.phonon.plotter.md#pymatgen.phonon.plotter.freq_units)


* [pymatgen.phonon.thermal_displacements module](pymatgen.phonon.thermal_displacements.md)


    * [`ThermalDisplacementMatrices`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices)


        * [`ThermalDisplacementMatrices.B`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.B)


        * [`ThermalDisplacementMatrices.U1U2U3`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.U1U2U3)


        * [`ThermalDisplacementMatrices.Ucif`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.Ucif)


        * [`ThermalDisplacementMatrices.Ustar`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.Ustar)


        * [`ThermalDisplacementMatrices.beta`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.beta)


        * [`ThermalDisplacementMatrices.compute_directionality_quality_criterion()`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.compute_directionality_quality_criterion)


        * [`ThermalDisplacementMatrices.from_Ucif()`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.from_Ucif)


        * [`ThermalDisplacementMatrices.from_cif_P1()`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.from_cif_P1)


        * [`ThermalDisplacementMatrices.from_structure_with_site_properties_Ucif()`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.from_structure_with_site_properties_Ucif)


        * [`ThermalDisplacementMatrices.get_full_matrix()`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.get_full_matrix)


        * [`ThermalDisplacementMatrices.get_reduced_matrix()`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.get_reduced_matrix)


        * [`ThermalDisplacementMatrices.ratio_prolate`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.ratio_prolate)


        * [`ThermalDisplacementMatrices.to_structure_with_site_properties_Ucif()`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.to_structure_with_site_properties_Ucif)


        * [`ThermalDisplacementMatrices.visualize_directionality_quality_criterion()`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.visualize_directionality_quality_criterion)


        * [`ThermalDisplacementMatrices.write_cif()`](pymatgen.phonon.thermal_displacements.md#pymatgen.phonon.thermal_displacements.ThermalDisplacementMatrices.write_cif)