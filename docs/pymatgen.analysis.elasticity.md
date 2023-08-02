---
layout: default
title: pymatgen.analysis.elasticity.md
nav_exclude: true
---

# pymatgen.analysis.elasticity package

Package for analyzing elastic tensors and properties.



* [pymatgen.analysis.elasticity.elastic module](pymatgen.analysis.elasticity.elastic.md)


    * [`ComplianceTensor`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ComplianceTensor)


    * [`ElasticTensor`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor)


        * [`ElasticTensor.cahill_thermalcond()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.cahill_thermalcond)


        * [`ElasticTensor.clarke_thermalcond()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.clarke_thermalcond)


        * [`ElasticTensor.compliance_tensor`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.compliance_tensor)


        * [`ElasticTensor.debye_temperature()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.debye_temperature)


        * [`ElasticTensor.directional_elastic_mod()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.directional_elastic_mod)


        * [`ElasticTensor.directional_poisson_ratio()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.directional_poisson_ratio)


        * [`ElasticTensor.from_independent_strains()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.from_independent_strains)


        * [`ElasticTensor.from_pseudoinverse()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.from_pseudoinverse)


        * [`ElasticTensor.g_reuss`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.g_reuss)


        * [`ElasticTensor.g_voigt`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.g_voigt)


        * [`ElasticTensor.g_vrh`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.g_vrh)


        * [`ElasticTensor.get_structure_property_dict()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.get_structure_property_dict)


        * [`ElasticTensor.green_kristoffel()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.green_kristoffel)


        * [`ElasticTensor.homogeneous_poisson`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.homogeneous_poisson)


        * [`ElasticTensor.k_reuss`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.k_reuss)


        * [`ElasticTensor.k_voigt`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.k_voigt)


        * [`ElasticTensor.k_vrh`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.k_vrh)


        * [`ElasticTensor.long_v()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.long_v)


        * [`ElasticTensor.property_dict`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.property_dict)


        * [`ElasticTensor.snyder_ac()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.snyder_ac)


        * [`ElasticTensor.snyder_opt()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.snyder_opt)


        * [`ElasticTensor.snyder_total()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.snyder_total)


        * [`ElasticTensor.trans_v()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.trans_v)


        * [`ElasticTensor.universal_anisotropy`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.universal_anisotropy)


        * [`ElasticTensor.y_mod`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.y_mod)


    * [`ElasticTensorExpansion`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion)


        * [`ElasticTensorExpansion.calculate_stress()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.calculate_stress)


        * [`ElasticTensorExpansion.energy_density()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.energy_density)


        * [`ElasticTensorExpansion.from_diff_fit()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.from_diff_fit)


        * [`ElasticTensorExpansion.get_compliance_expansion()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_compliance_expansion)


        * [`ElasticTensorExpansion.get_effective_ecs()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_effective_ecs)


        * [`ElasticTensorExpansion.get_ggt()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_ggt)


        * [`ElasticTensorExpansion.get_gruneisen_parameter()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_gruneisen_parameter)


        * [`ElasticTensorExpansion.get_heat_capacity()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_heat_capacity)


        * [`ElasticTensorExpansion.get_stability_criteria()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_stability_criteria)


        * [`ElasticTensorExpansion.get_strain_from_stress()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_strain_from_stress)


        * [`ElasticTensorExpansion.get_symmetric_wallace_tensor()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_symmetric_wallace_tensor)


        * [`ElasticTensorExpansion.get_tgt()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_tgt)


        * [`ElasticTensorExpansion.get_wallace_tensor()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_wallace_tensor)


        * [`ElasticTensorExpansion.get_yield_stress()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_yield_stress)


        * [`ElasticTensorExpansion.omega()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.omega)


        * [`ElasticTensorExpansion.order`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.order)


        * [`ElasticTensorExpansion.thermal_expansion_coeff()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.thermal_expansion_coeff)


    * [`NthOrderElasticTensor`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor)


        * [`NthOrderElasticTensor.GPa_to_eV_A3`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.GPa_to_eV_A3)


        * [`NthOrderElasticTensor.calculate_stress()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.calculate_stress)


        * [`NthOrderElasticTensor.energy_density()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.energy_density)


        * [`NthOrderElasticTensor.from_diff_fit()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.from_diff_fit)


        * [`NthOrderElasticTensor.order`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.order)


        * [`NthOrderElasticTensor.symbol`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.symbol)


    * [`diff_fit()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.diff_fit)


    * [`find_eq_stress()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.find_eq_stress)


    * [`generate_pseudo()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.generate_pseudo)


    * [`get_diff_coeff()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.get_diff_coeff)


    * [`get_strain_state_dict()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.get_strain_state_dict)


    * [`get_symbol_list()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.get_symbol_list)


    * [`raise_error_if_unphysical()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.raise_error_if_unphysical)


    * [`subs()`](pymatgen.analysis.elasticity.elastic.md#pymatgen.analysis.elasticity.elastic.subs)


* [pymatgen.analysis.elasticity.strain module](pymatgen.analysis.elasticity.strain.md)


    * [`Deformation`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Deformation)


        * [`Deformation.apply_to_structure()`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Deformation.apply_to_structure)


        * [`Deformation.from_index_amount()`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Deformation.from_index_amount)


        * [`Deformation.get_perturbed_indices()`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Deformation.get_perturbed_indices)


        * [`Deformation.green_lagrange_strain`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Deformation.green_lagrange_strain)


        * [`Deformation.is_independent()`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Deformation.is_independent)


        * [`Deformation.symbol`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Deformation.symbol)


    * [`DeformedStructureSet`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.DeformedStructureSet)


    * [`Strain`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Strain)


        * [`Strain.from_deformation()`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Strain.from_deformation)


        * [`Strain.from_index_amount()`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Strain.from_index_amount)


        * [`Strain.get_deformation_matrix()`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Strain.get_deformation_matrix)


        * [`Strain.symbol`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Strain.symbol)


        * [`Strain.von_mises_strain`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.Strain.von_mises_strain)


    * [`convert_strain_to_deformation()`](pymatgen.analysis.elasticity.strain.md#pymatgen.analysis.elasticity.strain.convert_strain_to_deformation)


* [pymatgen.analysis.elasticity.stress module](pymatgen.analysis.elasticity.stress.md)


    * [`Stress`](pymatgen.analysis.elasticity.stress.md#pymatgen.analysis.elasticity.stress.Stress)


        * [`Stress.dev_principal_invariants`](pymatgen.analysis.elasticity.stress.md#pymatgen.analysis.elasticity.stress.Stress.dev_principal_invariants)


        * [`Stress.deviator_stress`](pymatgen.analysis.elasticity.stress.md#pymatgen.analysis.elasticity.stress.Stress.deviator_stress)


        * [`Stress.mean_stress`](pymatgen.analysis.elasticity.stress.md#pymatgen.analysis.elasticity.stress.Stress.mean_stress)


        * [`Stress.piola_kirchoff_1()`](pymatgen.analysis.elasticity.stress.md#pymatgen.analysis.elasticity.stress.Stress.piola_kirchoff_1)


        * [`Stress.piola_kirchoff_2()`](pymatgen.analysis.elasticity.stress.md#pymatgen.analysis.elasticity.stress.Stress.piola_kirchoff_2)


        * [`Stress.symbol`](pymatgen.analysis.elasticity.stress.md#pymatgen.analysis.elasticity.stress.Stress.symbol)


        * [`Stress.von_mises`](pymatgen.analysis.elasticity.stress.md#pymatgen.analysis.elasticity.stress.Stress.von_mises)