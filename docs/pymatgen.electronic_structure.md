---
layout: default
title: pymatgen.electronic_structure.md
nav_exclude: true
---

# pymatgen.electronic_structure package

This package contains electronic structure related tools and analyses.



* [pymatgen.electronic_structure.bandstructure module](pymatgen.electronic_structure.bandstructure.md)


    * [`BandStructure`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure)


        * [`BandStructure.lattice_rec`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.lattice_rec)


        * [`BandStructure.efermi`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.efermi)


        * [`BandStructure.is_spin_polarized`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.is_spin_polarized)


        * [`BandStructure.bands`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.bands)


        * [`BandStructure.nb_bands`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.nb_bands)


        * [`BandStructure.structure`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.structure)


        * [`BandStructure.projections`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.projections)


        * [`BandStructure.as_dict()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.as_dict)


        * [`BandStructure.from_dict()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.from_dict)


        * [`BandStructure.from_old_dict()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.from_old_dict)


        * [`BandStructure.get_band_gap()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.get_band_gap)


        * [`BandStructure.get_cbm()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.get_cbm)


        * [`BandStructure.get_direct_band_gap()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.get_direct_band_gap)


        * [`BandStructure.get_direct_band_gap_dict()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.get_direct_band_gap_dict)


        * [`BandStructure.get_kpoint_degeneracy()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.get_kpoint_degeneracy)


        * [`BandStructure.get_projection_on_elements()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.get_projection_on_elements)


        * [`BandStructure.get_projections_on_elements_and_orbitals()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.get_projections_on_elements_and_orbitals)


        * [`BandStructure.get_sym_eq_kpoints()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.get_sym_eq_kpoints)


        * [`BandStructure.get_vbm()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.get_vbm)


        * [`BandStructure.is_metal()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructure.is_metal)


    * [`BandStructureSymmLine`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine)


        * [`BandStructureSymmLine.apply_scissor()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine.apply_scissor)


        * [`BandStructureSymmLine.as_dict()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine.as_dict)


        * [`BandStructureSymmLine.get_branch()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine.get_branch)


        * [`BandStructureSymmLine.get_equivalent_kpoints()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.BandStructureSymmLine.get_equivalent_kpoints)


    * [`Kpoint`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint)


        * [`Kpoint.a`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint.a)


        * [`Kpoint.as_dict()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint.as_dict)


        * [`Kpoint.b`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint.b)


        * [`Kpoint.c`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint.c)


        * [`Kpoint.cart_coords`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint.cart_coords)


        * [`Kpoint.frac_coords`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint.frac_coords)


        * [`Kpoint.from_dict()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint.from_dict)


        * [`Kpoint.label`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint.label)


        * [`Kpoint.lattice`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.Kpoint.lattice)


    * [`LobsterBandStructureSymmLine`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.LobsterBandStructureSymmLine)


        * [`LobsterBandStructureSymmLine.as_dict()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.LobsterBandStructureSymmLine.as_dict)


        * [`LobsterBandStructureSymmLine.from_dict()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.LobsterBandStructureSymmLine.from_dict)


        * [`LobsterBandStructureSymmLine.from_old_dict()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.LobsterBandStructureSymmLine.from_old_dict)


        * [`LobsterBandStructureSymmLine.get_projection_on_elements()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.LobsterBandStructureSymmLine.get_projection_on_elements)


        * [`LobsterBandStructureSymmLine.get_projections_on_elements_and_orbitals()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.LobsterBandStructureSymmLine.get_projections_on_elements_and_orbitals)


    * [`get_reconstructed_band_structure()`](pymatgen.electronic_structure.bandstructure.md#pymatgen.electronic_structure.bandstructure.get_reconstructed_band_structure)


* [pymatgen.electronic_structure.boltztrap module](pymatgen.electronic_structure.boltztrap.md)


    * [`BoltztrapAnalyzer`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer)


        * [`BoltztrapAnalyzer.as_dict()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.as_dict)


        * [`BoltztrapAnalyzer.check_acc_bzt_bands()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.check_acc_bzt_bands)


        * [`BoltztrapAnalyzer.from_dict()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.from_dict)


        * [`BoltztrapAnalyzer.from_files()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.from_files)


        * [`BoltztrapAnalyzer.get_average_eff_mass()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_average_eff_mass)


        * [`BoltztrapAnalyzer.get_carrier_concentration()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_carrier_concentration)


        * [`BoltztrapAnalyzer.get_complete_dos()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_complete_dos)


        * [`BoltztrapAnalyzer.get_complexity_factor()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_complexity_factor)


        * [`BoltztrapAnalyzer.get_conductivity()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_conductivity)


        * [`BoltztrapAnalyzer.get_extreme()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_extreme)


        * [`BoltztrapAnalyzer.get_hall_carrier_concentration()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_hall_carrier_concentration)


        * [`BoltztrapAnalyzer.get_mu_bounds()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_mu_bounds)


        * [`BoltztrapAnalyzer.get_power_factor()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_power_factor)


        * [`BoltztrapAnalyzer.get_seebeck()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_seebeck)


        * [`BoltztrapAnalyzer.get_seebeck_eff_mass()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_seebeck_eff_mass)


        * [`BoltztrapAnalyzer.get_symm_bands()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_symm_bands)


        * [`BoltztrapAnalyzer.get_thermal_conductivity()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_thermal_conductivity)


        * [`BoltztrapAnalyzer.get_zt()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.get_zt)


        * [`BoltztrapAnalyzer.parse_cond_and_hall()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.parse_cond_and_hall)


        * [`BoltztrapAnalyzer.parse_intrans()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.parse_intrans)


        * [`BoltztrapAnalyzer.parse_outputtrans()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.parse_outputtrans)


        * [`BoltztrapAnalyzer.parse_struct()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.parse_struct)


        * [`BoltztrapAnalyzer.parse_transdos()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapAnalyzer.parse_transdos)


    * [`BoltztrapError`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapError)


    * [`BoltztrapRunner`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner)


        * [`BoltztrapRunner.as_dict()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner.as_dict)


        * [`BoltztrapRunner.bs`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner.bs)


        * [`BoltztrapRunner.nelec`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner.nelec)


        * [`BoltztrapRunner.run()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner.run)


        * [`BoltztrapRunner.write_def()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner.write_def)


        * [`BoltztrapRunner.write_energy()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner.write_energy)


        * [`BoltztrapRunner.write_input()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner.write_input)


        * [`BoltztrapRunner.write_intrans()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner.write_intrans)


        * [`BoltztrapRunner.write_proj()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner.write_proj)


        * [`BoltztrapRunner.write_struct()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.BoltztrapRunner.write_struct)


    * [`compare_sym_bands()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.compare_sym_bands)


    * [`eta_from_seebeck()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.eta_from_seebeck)


    * [`read_cube_file()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.read_cube_file)


    * [`seebeck_eff_mass_from_carr()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.seebeck_eff_mass_from_carr)


    * [`seebeck_eff_mass_from_seebeck_carr()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.seebeck_eff_mass_from_seebeck_carr)


    * [`seebeck_spb()`](pymatgen.electronic_structure.boltztrap.md#pymatgen.electronic_structure.boltztrap.seebeck_spb)


* [pymatgen.electronic_structure.boltztrap2 module](pymatgen.electronic_structure.boltztrap2.md)


    * [`BandstructureLoader`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BandstructureLoader)


        * [`BandstructureLoader.bandana()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BandstructureLoader.bandana)


        * [`BandstructureLoader.get_lattvec()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BandstructureLoader.get_lattvec)


        * [`BandstructureLoader.get_volume()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BandstructureLoader.get_volume)


        * [`BandstructureLoader.set_upper_lower_bands()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BandstructureLoader.set_upper_lower_bands)


    * [`BztInterpolator`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztInterpolator)


        * [`BztInterpolator.get_band_structure()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztInterpolator.get_band_structure)


        * [`BztInterpolator.get_dos()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztInterpolator.get_dos)


        * [`BztInterpolator.get_partial_doses()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztInterpolator.get_partial_doses)


        * [`BztInterpolator.load()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztInterpolator.load)


        * [`BztInterpolator.save()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztInterpolator.save)


    * [`BztPlotter`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztPlotter)


        * [`BztPlotter.plot_bands()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztPlotter.plot_bands)


        * [`BztPlotter.plot_dos()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztPlotter.plot_dos)


        * [`BztPlotter.plot_props()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztPlotter.plot_props)


    * [`BztTransportProperties`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztTransportProperties)


        * [`BztTransportProperties.compute_properties_doping()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztTransportProperties.compute_properties_doping)


        * [`BztTransportProperties.load()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztTransportProperties.load)


        * [`BztTransportProperties.save()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.BztTransportProperties.save)


    * [`VasprunBSLoader`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.VasprunBSLoader)


        * [`VasprunBSLoader.bandana()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.VasprunBSLoader.bandana)


        * [`VasprunBSLoader.from_file()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.VasprunBSLoader.from_file)


        * [`VasprunBSLoader.get_lattvec()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.VasprunBSLoader.get_lattvec)


        * [`VasprunBSLoader.get_volume()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.VasprunBSLoader.get_volume)


    * [`VasprunLoader`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.VasprunLoader)


        * [`VasprunLoader.bandana()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.VasprunLoader.bandana)


        * [`VasprunLoader.from_file()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.VasprunLoader.from_file)


        * [`VasprunLoader.get_lattvec()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.VasprunLoader.get_lattvec)


        * [`VasprunLoader.get_volume()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.VasprunLoader.get_volume)


    * [`merge_up_down_doses()`](pymatgen.electronic_structure.boltztrap2.md#pymatgen.electronic_structure.boltztrap2.merge_up_down_doses)


* [pymatgen.electronic_structure.cohp module](pymatgen.electronic_structure.cohp.md)


    * [`Cohp`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.Cohp)


        * [`Cohp.as_dict()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.Cohp.as_dict)


        * [`Cohp.from_dict()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.Cohp.from_dict)


        * [`Cohp.get_cohp()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.Cohp.get_cohp)


        * [`Cohp.get_icohp()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.Cohp.get_icohp)


        * [`Cohp.get_interpolated_value()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.Cohp.get_interpolated_value)


        * [`Cohp.has_antibnd_states_below_efermi()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.Cohp.has_antibnd_states_below_efermi)


    * [`CompleteCohp`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.CompleteCohp)


        * [`CompleteCohp.as_dict()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.CompleteCohp.as_dict)


        * [`CompleteCohp.from_dict()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.CompleteCohp.from_dict)


        * [`CompleteCohp.from_file()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.CompleteCohp.from_file)


        * [`CompleteCohp.get_cohp_by_label()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.CompleteCohp.get_cohp_by_label)


        * [`CompleteCohp.get_orbital_resolved_cohp()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.CompleteCohp.get_orbital_resolved_cohp)


        * [`CompleteCohp.get_summed_cohp_by_label_and_orbital_list()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.CompleteCohp.get_summed_cohp_by_label_and_orbital_list)


        * [`CompleteCohp.get_summed_cohp_by_label_list()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.CompleteCohp.get_summed_cohp_by_label_list)


    * [`IcohpCollection`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpCollection)


        * [`IcohpCollection.are_coops`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpCollection.are_coops)


        * [`IcohpCollection.are_cobis`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpCollection.are_cobis)


        * [`IcohpCollection.is_spin_polarized`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpCollection.is_spin_polarized)


        * [`IcohpCollection.are_cobis`](pymatgen.electronic_structure.cohp.md#id0)


        * [`IcohpCollection.are_coops`](pymatgen.electronic_structure.cohp.md#id1)


        * [`IcohpCollection.extremum_icohpvalue()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpCollection.extremum_icohpvalue)


        * [`IcohpCollection.get_icohp_by_label()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpCollection.get_icohp_by_label)


        * [`IcohpCollection.get_icohp_dict_by_bondlengths()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpCollection.get_icohp_dict_by_bondlengths)


        * [`IcohpCollection.get_icohp_dict_of_site()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpCollection.get_icohp_dict_of_site)


        * [`IcohpCollection.get_summed_icohp_by_label_list()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpCollection.get_summed_icohp_by_label_list)


        * [`IcohpCollection.is_spin_polarized`](pymatgen.electronic_structure.cohp.md#id2)


    * [`IcohpValue`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpValue)


        * [`IcohpValue.num_bonds`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpValue.num_bonds)


        * [`IcohpValue.are_coops`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpValue.are_coops)


        * [`IcohpValue.are_cobis`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpValue.are_cobis)


        * [`IcohpValue.icohp`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpValue.icohp)


        * [`IcohpValue.are_cobis`](pymatgen.electronic_structure.cohp.md#id3)


        * [`IcohpValue.are_coops`](pymatgen.electronic_structure.cohp.md#id4)


        * [`IcohpValue.icohp`](pymatgen.electronic_structure.cohp.md#id5)


        * [`IcohpValue.icohpvalue()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpValue.icohpvalue)


        * [`IcohpValue.icohpvalue_orbital()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpValue.icohpvalue_orbital)


        * [`IcohpValue.is_spin_polarized`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpValue.is_spin_polarized)


        * [`IcohpValue.num_bonds`](pymatgen.electronic_structure.cohp.md#id6)


        * [`IcohpValue.summed_icohp`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpValue.summed_icohp)


        * [`IcohpValue.summed_orbital_icohp`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.IcohpValue.summed_orbital_icohp)


    * [`get_integrated_cohp_in_energy_range()`](pymatgen.electronic_structure.cohp.md#pymatgen.electronic_structure.cohp.get_integrated_cohp_in_energy_range)


* [pymatgen.electronic_structure.core module](pymatgen.electronic_structure.core.md)


    * [`Magmom`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom)


        * [`Magmom.are_collinear()`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.are_collinear)


        * [`Magmom.from_global_moment_and_saxis()`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.from_global_moment_and_saxis)


        * [`Magmom.from_moment_relative_to_crystal_axes()`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.from_moment_relative_to_crystal_axes)


        * [`Magmom.get_00t_magmom_with_xyz_saxis()`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.get_00t_magmom_with_xyz_saxis)


        * [`Magmom.get_consistent_set_and_saxis()`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.get_consistent_set_and_saxis)


        * [`Magmom.get_moment()`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.get_moment)


        * [`Magmom.get_moment_relative_to_crystal_axes()`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.get_moment_relative_to_crystal_axes)


        * [`Magmom.get_suggested_saxis()`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.get_suggested_saxis)


        * [`Magmom.get_xyz_magmom_with_001_saxis()`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.get_xyz_magmom_with_001_saxis)


        * [`Magmom.global_moment`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.global_moment)


        * [`Magmom.have_consistent_saxis()`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.have_consistent_saxis)


        * [`Magmom.projection`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Magmom.projection)


    * [`Orbital`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital)


        * [`Orbital.dx2`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.dx2)


        * [`Orbital.dxy`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.dxy)


        * [`Orbital.dxz`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.dxz)


        * [`Orbital.dyz`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.dyz)


        * [`Orbital.dz2`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.dz2)


        * [`Orbital.f0`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.f0)


        * [`Orbital.f1`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.f1)


        * [`Orbital.f2`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.f2)


        * [`Orbital.f3`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.f3)


        * [`Orbital.f_1`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.f_1)


        * [`Orbital.f_2`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.f_2)


        * [`Orbital.f_3`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.f_3)


        * [`Orbital.orbital_type`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.orbital_type)


        * [`Orbital.px`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.px)


        * [`Orbital.py`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.py)


        * [`Orbital.pz`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.pz)


        * [`Orbital.s`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Orbital.s)


    * [`OrbitalType`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType)


        * [`OrbitalType.d`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType.d)


        * [`OrbitalType.f`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType.f)


        * [`OrbitalType.p`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType.p)


        * [`OrbitalType.s`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.OrbitalType.s)


    * [`Spin`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin)


        * [`Spin.down`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin.down)


        * [`Spin.up`](pymatgen.electronic_structure.core.md#pymatgen.electronic_structure.core.Spin.up)


* [pymatgen.electronic_structure.dos module](pymatgen.electronic_structure.dos.md)


    * [`CompleteDos`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos)


        * [`CompleteDos.structure`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.structure)


        * [`CompleteDos.pdos`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.pdos)


        * [`CompleteDos.as_dict()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.as_dict)


        * [`CompleteDos.fp_to_dict()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.fp_to_dict)


        * [`CompleteDos.from_dict()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.from_dict)


        * [`CompleteDos.get_band_center()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_band_center)


        * [`CompleteDos.get_band_filling()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_band_filling)


        * [`CompleteDos.get_band_kurtosis()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_band_kurtosis)


        * [`CompleteDos.get_band_skewness()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_band_skewness)


        * [`CompleteDos.get_band_width()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_band_width)


        * [`CompleteDos.get_dos_fp()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_dos_fp)


        * [`CompleteDos.get_dos_fp_similarity()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_dos_fp_similarity)


        * [`CompleteDos.get_element_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_element_dos)


        * [`CompleteDos.get_element_spd_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_element_spd_dos)


        * [`CompleteDos.get_hilbert_transform()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_hilbert_transform)


        * [`CompleteDos.get_n_moment()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_n_moment)


        * [`CompleteDos.get_normalized()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_normalized)


        * [`CompleteDos.get_site_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_site_dos)


        * [`CompleteDos.get_site_orbital_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_site_orbital_dos)


        * [`CompleteDos.get_site_spd_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_site_spd_dos)


        * [`CompleteDos.get_site_t2g_eg_resolved_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_site_t2g_eg_resolved_dos)


        * [`CompleteDos.get_spd_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_spd_dos)


        * [`CompleteDos.get_upper_band_edge()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.get_upper_band_edge)


        * [`CompleteDos.spin_polarization`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.CompleteDos.spin_polarization)


    * [`DOS`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.DOS)


        * [`DOS.XLABEL`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.DOS.XLABEL)


        * [`DOS.YLABEL`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.DOS.YLABEL)


        * [`DOS.get_cbm_vbm()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.DOS.get_cbm_vbm)


        * [`DOS.get_gap()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.DOS.get_gap)


        * [`DOS.get_interpolated_gap()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.DOS.get_interpolated_gap)


    * [`Dos`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos)


        * [`Dos.as_dict()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos.as_dict)


        * [`Dos.from_dict()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos.from_dict)


        * [`Dos.get_cbm_vbm()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos.get_cbm_vbm)


        * [`Dos.get_densities()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos.get_densities)


        * [`Dos.get_gap()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos.get_gap)


        * [`Dos.get_interpolated_gap()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos.get_interpolated_gap)


        * [`Dos.get_interpolated_value()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos.get_interpolated_value)


        * [`Dos.get_smeared_densities()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.Dos.get_smeared_densities)


    * [`FermiDos`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.FermiDos)


        * [`FermiDos.as_dict()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.FermiDos.as_dict)


        * [`FermiDos.from_dict()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.FermiDos.from_dict)


        * [`FermiDos.get_doping()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.FermiDos.get_doping)


        * [`FermiDos.get_fermi()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.FermiDos.get_fermi)


        * [`FermiDos.get_fermi_interextrapolated()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.FermiDos.get_fermi_interextrapolated)


    * [`LobsterCompleteDos`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.LobsterCompleteDos)


        * [`LobsterCompleteDos.from_dict()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.LobsterCompleteDos.from_dict)


        * [`LobsterCompleteDos.get_element_spd_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.LobsterCompleteDos.get_element_spd_dos)


        * [`LobsterCompleteDos.get_site_orbital_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.LobsterCompleteDos.get_site_orbital_dos)


        * [`LobsterCompleteDos.get_site_t2g_eg_resolved_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.LobsterCompleteDos.get_site_t2g_eg_resolved_dos)


        * [`LobsterCompleteDos.get_spd_dos()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.LobsterCompleteDos.get_spd_dos)


    * [`add_densities()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.add_densities)


    * [`f0()`](pymatgen.electronic_structure.dos.md#pymatgen.electronic_structure.dos.f0)


* [pymatgen.electronic_structure.plotter module](pymatgen.electronic_structure.plotter.md)


    * [`BSDOSPlotter`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSDOSPlotter)


        * [`BSDOSPlotter.get_plot()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSDOSPlotter.get_plot)


    * [`BSPlotter`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotter)


        * [`BSPlotter.add_bs()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotter.add_bs)


        * [`BSPlotter.bs_plot_data()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotter.bs_plot_data)


        * [`BSPlotter.get_plot()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotter.get_plot)


        * [`BSPlotter.get_ticks()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotter.get_ticks)


        * [`BSPlotter.get_ticks_old()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotter.get_ticks_old)


        * [`BSPlotter.plot_brillouin()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotter.plot_brillouin)


        * [`BSPlotter.plot_compare()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotter.plot_compare)


        * [`BSPlotter.save_plot()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotter.save_plot)


        * [`BSPlotter.show()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotter.show)


    * [`BSPlotterProjected`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotterProjected)


        * [`BSPlotterProjected.get_elt_projected_plots()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotterProjected.get_elt_projected_plots)


        * [`BSPlotterProjected.get_elt_projected_plots_color()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotterProjected.get_elt_projected_plots_color)


        * [`BSPlotterProjected.get_projected_plots_dots()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotterProjected.get_projected_plots_dots)


        * [`BSPlotterProjected.get_projected_plots_dots_patom_pmorb()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BSPlotterProjected.get_projected_plots_dots_patom_pmorb)


    * [`BoltztrapPlotter`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter)


        * [`BoltztrapPlotter.plot_carriers()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_carriers)


        * [`BoltztrapPlotter.plot_complexity_factor_mu()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_complexity_factor_mu)


        * [`BoltztrapPlotter.plot_conductivity_dop()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_conductivity_dop)


        * [`BoltztrapPlotter.plot_conductivity_mu()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_conductivity_mu)


        * [`BoltztrapPlotter.plot_conductivity_temp()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_conductivity_temp)


        * [`BoltztrapPlotter.plot_dos()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_dos)


        * [`BoltztrapPlotter.plot_eff_mass_dop()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_eff_mass_dop)


        * [`BoltztrapPlotter.plot_eff_mass_temp()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_eff_mass_temp)


        * [`BoltztrapPlotter.plot_hall_carriers()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_hall_carriers)


        * [`BoltztrapPlotter.plot_power_factor_dop()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_power_factor_dop)


        * [`BoltztrapPlotter.plot_power_factor_mu()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_power_factor_mu)


        * [`BoltztrapPlotter.plot_power_factor_temp()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_power_factor_temp)


        * [`BoltztrapPlotter.plot_seebeck_dop()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_seebeck_dop)


        * [`BoltztrapPlotter.plot_seebeck_eff_mass_mu()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_seebeck_eff_mass_mu)


        * [`BoltztrapPlotter.plot_seebeck_mu()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_seebeck_mu)


        * [`BoltztrapPlotter.plot_seebeck_temp()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_seebeck_temp)


        * [`BoltztrapPlotter.plot_zt_dop()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_zt_dop)


        * [`BoltztrapPlotter.plot_zt_mu()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_zt_mu)


        * [`BoltztrapPlotter.plot_zt_temp()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.BoltztrapPlotter.plot_zt_temp)


    * [`CohpPlotter`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.CohpPlotter)


        * [`CohpPlotter.add_cohp()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.CohpPlotter.add_cohp)


        * [`CohpPlotter.add_cohp_dict()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.CohpPlotter.add_cohp_dict)


        * [`CohpPlotter.get_cohp_dict()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.CohpPlotter.get_cohp_dict)


        * [`CohpPlotter.get_plot()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.CohpPlotter.get_plot)


        * [`CohpPlotter.save_plot()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.CohpPlotter.save_plot)


        * [`CohpPlotter.show()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.CohpPlotter.show)


    * [`DosPlotter`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.DosPlotter)


        * [`DosPlotter.add_dos()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.DosPlotter.add_dos)


        * [`DosPlotter.add_dos_dict()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.DosPlotter.add_dos_dict)


        * [`DosPlotter.get_dos_dict()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.DosPlotter.get_dos_dict)


        * [`DosPlotter.get_plot()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.DosPlotter.get_plot)


        * [`DosPlotter.save_plot()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.DosPlotter.save_plot)


        * [`DosPlotter.show()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.DosPlotter.show)


    * [`fold_point()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.fold_point)


    * [`plot_brillouin_zone()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.plot_brillouin_zone)


    * [`plot_brillouin_zone_from_kpath()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.plot_brillouin_zone_from_kpath)


    * [`plot_ellipsoid()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.plot_ellipsoid)


    * [`plot_fermi_surface()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.plot_fermi_surface)


    * [`plot_labels()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.plot_labels)


    * [`plot_lattice_vectors()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.plot_lattice_vectors)


    * [`plot_path()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.plot_path)


    * [`plot_points()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.plot_points)


    * [`plot_wigner_seitz()`](pymatgen.electronic_structure.plotter.md#pymatgen.electronic_structure.plotter.plot_wigner_seitz)