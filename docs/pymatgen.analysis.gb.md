---
layout: default
title: pymatgen.analysis.gb.md
nav_exclude: true
---

# pymatgen.analysis.gb package

This package implements various grain boundary analyses.



* [pymatgen.analysis.gb.grain module](pymatgen.analysis.gb.grain.md)


    * [`GrainBoundary`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundary)


        * [`GrainBoundary.as_dict()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundary.as_dict)


        * [`GrainBoundary.bottom_grain`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundary.bottom_grain)


        * [`GrainBoundary.coincidents`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundary.coincidents)


        * [`GrainBoundary.copy()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundary.copy)


        * [`GrainBoundary.from_dict()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundary.from_dict)


        * [`GrainBoundary.get_sorted_structure()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundary.get_sorted_structure)


        * [`GrainBoundary.sigma`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundary.sigma)


        * [`GrainBoundary.sigma_from_site_prop`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundary.sigma_from_site_prop)


        * [`GrainBoundary.top_grain`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundary.top_grain)


    * [`GrainBoundaryGenerator`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator)


        * [`GrainBoundaryGenerator.enum_possible_plane_cubic()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_possible_plane_cubic)


        * [`GrainBoundaryGenerator.enum_sigma_cubic()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_sigma_cubic)


        * [`GrainBoundaryGenerator.enum_sigma_hex()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_sigma_hex)


        * [`GrainBoundaryGenerator.enum_sigma_ort()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_sigma_ort)


        * [`GrainBoundaryGenerator.enum_sigma_rho()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_sigma_rho)


        * [`GrainBoundaryGenerator.enum_sigma_tet()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_sigma_tet)


        * [`GrainBoundaryGenerator.gb_from_parameters()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.gb_from_parameters)


        * [`GrainBoundaryGenerator.get_ratio()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.get_ratio)


        * [`GrainBoundaryGenerator.get_rotation_angle_from_sigma()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.get_rotation_angle_from_sigma)


        * [`GrainBoundaryGenerator.get_trans_mat()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.get_trans_mat)


        * [`GrainBoundaryGenerator.reduce_mat()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.reduce_mat)


        * [`GrainBoundaryGenerator.slab_from_csl()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.slab_from_csl)


        * [`GrainBoundaryGenerator.vec_to_surface()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.vec_to_surface)


    * [`fix_pbc()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.fix_pbc)


    * [`symm_group_cubic()`](pymatgen.analysis.gb.grain.md#pymatgen.analysis.gb.grain.symm_group_cubic)