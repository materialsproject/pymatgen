---
layout: default
title: pymatgen.analysis.interfaces.md
nav_exclude: true
---

# pymatgen.analysis.interfaces package

Module that implements various algorithms related to interface construction and analysis.



* [pymatgen.analysis.interfaces.coherent_interfaces module](pymatgen.analysis.interfaces.coherent_interfaces.md)


    * [`CoherentInterfaceBuilder`](pymatgen.analysis.interfaces.coherent_interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.CoherentInterfaceBuilder)


        * [`CoherentInterfaceBuilder.get_interfaces()`](pymatgen.analysis.interfaces.coherent_interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.CoherentInterfaceBuilder.get_interfaces)


    * [`from_2d_to_3d()`](pymatgen.analysis.interfaces.coherent_interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.from_2d_to_3d)


    * [`get_2d_transform()`](pymatgen.analysis.interfaces.coherent_interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.get_2d_transform)


    * [`get_rot_3d_for_2d()`](pymatgen.analysis.interfaces.coherent_interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.get_rot_3d_for_2d)


* [pymatgen.analysis.interfaces.substrate_analyzer module](pymatgen.analysis.interfaces.substrate_analyzer.md)


    * [`SubstrateAnalyzer`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateAnalyzer)


        * [`SubstrateAnalyzer.calculate()`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateAnalyzer.calculate)


        * [`SubstrateAnalyzer.generate_surface_vectors()`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateAnalyzer.generate_surface_vectors)


    * [`SubstrateMatch`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch)


        * [`SubstrateMatch.elastic_energy`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.elastic_energy)


        * [`SubstrateMatch.film_miller`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.film_miller)


        * [`SubstrateMatch.from_zsl()`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.from_zsl)


        * [`SubstrateMatch.ground_state_energy`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.ground_state_energy)


        * [`SubstrateMatch.strain`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.strain)


        * [`SubstrateMatch.substrate_miller`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.substrate_miller)


        * [`SubstrateMatch.total_energy`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.total_energy)


        * [`SubstrateMatch.von_mises_strain`](pymatgen.analysis.interfaces.substrate_analyzer.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.von_mises_strain)


* [pymatgen.analysis.interfaces.zsl module](pymatgen.analysis.interfaces.zsl.md)


    * [`ZSLGenerator`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLGenerator)


        * [`ZSLGenerator.generate_sl_transformation_sets()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLGenerator.generate_sl_transformation_sets)


        * [`ZSLGenerator.get_equiv_transformations()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLGenerator.get_equiv_transformations)


    * [`ZSLMatch`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch)


        * [`ZSLMatch.film_sl_vectors`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.film_sl_vectors)


        * [`ZSLMatch.film_transformation`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.film_transformation)


        * [`ZSLMatch.film_vectors`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.film_vectors)


        * [`ZSLMatch.match_area`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.match_area)


        * [`ZSLMatch.match_transformation`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.match_transformation)


        * [`ZSLMatch.substrate_sl_vectors`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.substrate_sl_vectors)


        * [`ZSLMatch.substrate_transformation`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.substrate_transformation)


        * [`ZSLMatch.substrate_vectors`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.substrate_vectors)


    * [`fast_norm()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.fast_norm)


    * [`gen_sl_transform_matrices()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.gen_sl_transform_matrices)


    * [`get_factors()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.get_factors)


    * [`is_same_vectors()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.is_same_vectors)


    * [`reduce_vectors()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.reduce_vectors)


    * [`rel_angle()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.rel_angle)


    * [`rel_strain()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.rel_strain)


    * [`vec_angle()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.vec_angle)


    * [`vec_area()`](pymatgen.analysis.interfaces.zsl.md#pymatgen.analysis.interfaces.zsl.vec_area)