---
layout: default
title: pymatgen.io.lobster.md
nav_exclude: true
---

# pymatgen.io.lobster package

This package implements modules for input and output to and from Lobster. It
imports the key classes form both lobster.inputs and lobster_outputs to allow most
classes to be simply called as pymatgen.io.lobster.Lobsterin for example, to retain
backwards compatibility.



* [pymatgen.io.lobster.inputs module](pymatgen.io.lobster.inputs.md)


    * [`Lobsterin`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin)


        * [`Lobsterin.AVAILABLEKEYWORDS`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.AVAILABLEKEYWORDS)


        * [`Lobsterin.BOOLEAN_KEYWORDS`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.BOOLEAN_KEYWORDS)


        * [`Lobsterin.FLOAT_KEYWORDS`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.FLOAT_KEYWORDS)


        * [`Lobsterin.LISTKEYWORDS`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.LISTKEYWORDS)


        * [`Lobsterin.STRING_KEYWORDS`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.STRING_KEYWORDS)


        * [`Lobsterin.as_dict()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.as_dict)


        * [`Lobsterin.diff()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.diff)


        * [`Lobsterin.from_dict()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.from_dict)


        * [`Lobsterin.from_file()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.from_file)


        * [`Lobsterin.get_all_possible_basis_functions()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.get_all_possible_basis_functions)


        * [`Lobsterin.get_basis()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.get_basis)


        * [`Lobsterin.standard_calculations_from_vasp_files()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.standard_calculations_from_vasp_files)


        * [`Lobsterin.write_INCAR()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.write_INCAR)


        * [`Lobsterin.write_KPOINTS()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.write_KPOINTS)


        * [`Lobsterin.write_POSCAR_with_standard_primitive()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.write_POSCAR_with_standard_primitive)


        * [`Lobsterin.write_lobsterin()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.Lobsterin.write_lobsterin)


    * [`get_all_possible_basis_combinations()`](pymatgen.io.lobster.inputs.md#pymatgen.io.lobster.inputs.get_all_possible_basis_combinations)


* [pymatgen.io.lobster.lobsterenv module](pymatgen.io.lobster.lobsterenv.md)


    * [`ICOHPNeighborsInfo`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo)


        * [`ICOHPNeighborsInfo.atoms`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.atoms)


        * [`ICOHPNeighborsInfo.central_isites`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.central_isites)


        * [`ICOHPNeighborsInfo.labels`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.labels)


        * [`ICOHPNeighborsInfo.list_icohps`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.list_icohps)


        * [`ICOHPNeighborsInfo.n_bonds`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.n_bonds)


        * [`ICOHPNeighborsInfo.total_icohp`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.ICOHPNeighborsInfo.total_icohp)


    * [`LobsterLightStructureEnvironments`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterLightStructureEnvironments)


        * [`LobsterLightStructureEnvironments.as_dict()`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterLightStructureEnvironments.as_dict)


        * [`LobsterLightStructureEnvironments.from_Lobster()`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterLightStructureEnvironments.from_Lobster)


        * [`LobsterLightStructureEnvironments.uniquely_determines_coordination_environments`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterLightStructureEnvironments.uniquely_determines_coordination_environments)


    * [`LobsterNeighbors`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors)


        * [`LobsterNeighbors.anion_types`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.anion_types)


        * [`LobsterNeighbors.get_anion_types()`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_anion_types)


        * [`LobsterNeighbors.get_info_cohps_to_neighbors()`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_info_cohps_to_neighbors)


        * [`LobsterNeighbors.get_info_icohps_between_neighbors()`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_info_icohps_between_neighbors)


        * [`LobsterNeighbors.get_info_icohps_to_neighbors()`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_info_icohps_to_neighbors)


        * [`LobsterNeighbors.get_light_structure_environment()`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_light_structure_environment)


        * [`LobsterNeighbors.get_nn_info()`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.get_nn_info)


        * [`LobsterNeighbors.molecules_allowed`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.molecules_allowed)


        * [`LobsterNeighbors.plot_cohps_of_neighbors()`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.plot_cohps_of_neighbors)


        * [`LobsterNeighbors.structures_allowed`](pymatgen.io.lobster.lobsterenv.md#pymatgen.io.lobster.lobsterenv.LobsterNeighbors.structures_allowed)


* [pymatgen.io.lobster.outputs module](pymatgen.io.lobster.outputs.md)


    * [`Bandoverlaps`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Bandoverlaps)


        * [`Bandoverlaps.has_good_quality_check_occupied_bands()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Bandoverlaps.has_good_quality_check_occupied_bands)


        * [`Bandoverlaps.has_good_quality_maxDeviation()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Bandoverlaps.has_good_quality_maxDeviation)


    * [`Charge`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Charge)


        * [`Charge.get_structure_with_charges()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Charge.get_structure_with_charges)


    * [`Cohpcar`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Cohpcar)


    * [`Doscar`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Doscar)


        * [`Doscar.completedos`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Doscar.completedos)


        * [`Doscar.pdos`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Doscar.pdos)


        * [`Doscar.tdos`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Doscar.tdos)


        * [`Doscar.energies`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Doscar.energies)


        * [`Doscar.tdensities`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Doscar.tdensities)


        * [`Doscar.energies`](pymatgen.io.lobster.outputs.md#id0)


        * [`Doscar.energies`](pymatgen.io.lobster.outputs.md#id1)


        * [`Doscar.is_spin_polarized`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Doscar.is_spin_polarized)


        * [`Doscar.completedos`](pymatgen.io.lobster.outputs.md#id2)


        * [`Doscar.energies`](pymatgen.io.lobster.outputs.md#id3)


        * [`Doscar.is_spin_polarized`](pymatgen.io.lobster.outputs.md#id4)


        * [`Doscar.itdensities`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Doscar.itdensities)


        * [`Doscar.pdos`](pymatgen.io.lobster.outputs.md#id5)


        * [`Doscar.tdensities`](pymatgen.io.lobster.outputs.md#id6)


        * [`Doscar.tdos`](pymatgen.io.lobster.outputs.md#id7)


    * [`Fatband`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Fatband)


        * [`Fatband.get_bandstructure()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Fatband.get_bandstructure)


    * [`Grosspop`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Grosspop)


        * [`Grosspop.get_structure_with_total_grosspop()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Grosspop.get_structure_with_total_grosspop)


    * [`Icohplist`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Icohplist)


        * [`Icohplist.icohpcollection`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Icohplist.icohpcollection)


        * [`Icohplist.icohplist`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Icohplist.icohplist)


    * [`Lobsterout`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Lobsterout)


        * [`Lobsterout.get_doc()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Lobsterout.get_doc)


    * [`MadelungEnergies`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.MadelungEnergies)


    * [`SitePotential`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.SitePotential)


        * [`SitePotential.get_structure_with_site_potentials()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.SitePotential.get_structure_with_site_potentials)


    * [`Wavefunction`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Wavefunction)


        * [`Wavefunction.get_volumetricdata_density()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Wavefunction.get_volumetricdata_density)


        * [`Wavefunction.get_volumetricdata_imaginary()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Wavefunction.get_volumetricdata_imaginary)


        * [`Wavefunction.get_volumetricdata_real()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Wavefunction.get_volumetricdata_real)


        * [`Wavefunction.set_volumetric_data()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Wavefunction.set_volumetric_data)


        * [`Wavefunction.write_file()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.Wavefunction.write_file)


    * [`get_orb_from_str()`](pymatgen.io.lobster.outputs.md#pymatgen.io.lobster.outputs.get_orb_from_str)