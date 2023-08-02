---
layout: default
title: pymatgen.analysis.chemenv.utils.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.utils package

Utility package for chemenv.



* [pymatgen.analysis.chemenv.utils.chemenv_config module](pymatgen.analysis.chemenv.utils.chemenv_config.md)


    * [`ChemEnvConfig`](pymatgen.analysis.chemenv.utils.chemenv_config.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig)


        * [`ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS`](pymatgen.analysis.chemenv.utils.chemenv_config.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS)


        * [`ChemEnvConfig.auto_load()`](pymatgen.analysis.chemenv.utils.chemenv_config.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.auto_load)


        * [`ChemEnvConfig.has_materials_project_access`](pymatgen.analysis.chemenv.utils.chemenv_config.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.has_materials_project_access)


        * [`ChemEnvConfig.package_options_description()`](pymatgen.analysis.chemenv.utils.chemenv_config.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.package_options_description)


        * [`ChemEnvConfig.save()`](pymatgen.analysis.chemenv.utils.chemenv_config.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.save)


        * [`ChemEnvConfig.setup()`](pymatgen.analysis.chemenv.utils.chemenv_config.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.setup)


        * [`ChemEnvConfig.setup_package_options()`](pymatgen.analysis.chemenv.utils.chemenv_config.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.setup_package_options)


* [pymatgen.analysis.chemenv.utils.chemenv_errors module](pymatgen.analysis.chemenv.utils.chemenv_errors.md)


    * [`AbstractChemenvError`](pymatgen.analysis.chemenv.utils.chemenv_errors.md#pymatgen.analysis.chemenv.utils.chemenv_errors.AbstractChemenvError)


    * [`ChemenvError`](pymatgen.analysis.chemenv.utils.chemenv_errors.md#pymatgen.analysis.chemenv.utils.chemenv_errors.ChemenvError)


    * [`EquivalentSiteSearchError`](pymatgen.analysis.chemenv.utils.chemenv_errors.md#pymatgen.analysis.chemenv.utils.chemenv_errors.EquivalentSiteSearchError)


    * [`NeighborsNotComputedChemenvError`](pymatgen.analysis.chemenv.utils.chemenv_errors.md#pymatgen.analysis.chemenv.utils.chemenv_errors.NeighborsNotComputedChemenvError)


    * [`SolidAngleError`](pymatgen.analysis.chemenv.utils.chemenv_errors.md#pymatgen.analysis.chemenv.utils.chemenv_errors.SolidAngleError)


* [pymatgen.analysis.chemenv.utils.coordination_geometry_utils module](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md)


    * [`Plane`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane)


        * [`Plane.TEST_2D_POINTS`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.TEST_2D_POINTS)


        * [`Plane.a`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.a)


        * [`Plane.abcd`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.abcd)


        * [`Plane.b`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.b)


        * [`Plane.c`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.c)


        * [`Plane.coefficients`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.coefficients)


        * [`Plane.crosses_origin`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.crosses_origin)


        * [`Plane.d`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.d)


        * [`Plane.distance_to_origin`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.distance_to_origin)


        * [`Plane.distance_to_point()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.distance_to_point)


        * [`Plane.distances()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.distances)


        * [`Plane.distances_indices_groups()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.distances_indices_groups)


        * [`Plane.distances_indices_sorted()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.distances_indices_sorted)


        * [`Plane.fit_error()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.fit_error)


        * [`Plane.fit_least_square_distance_error()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.fit_least_square_distance_error)


        * [`Plane.fit_maximum_distance_error()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.fit_maximum_distance_error)


        * [`Plane.from_2points_and_origin()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_2points_and_origin)


        * [`Plane.from_3points()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_3points)


        * [`Plane.from_coefficients()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_coefficients)


        * [`Plane.from_npoints()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_npoints)


        * [`Plane.from_npoints_least_square_distance()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_npoints_least_square_distance)


        * [`Plane.from_npoints_maximum_distance()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_npoints_maximum_distance)


        * [`Plane.indices_separate()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.indices_separate)


        * [`Plane.init_3points()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.init_3points)


        * [`Plane.is_in_list()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.is_in_list)


        * [`Plane.is_in_plane()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.is_in_plane)


        * [`Plane.is_same_plane_as()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.is_same_plane_as)


        * [`Plane.orthonormal_vectors()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.orthonormal_vectors)


        * [`Plane.perpendicular_bisector()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.perpendicular_bisector)


        * [`Plane.project_and_to2dim()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.project_and_to2dim)


        * [`Plane.project_and_to2dim_ordered_indices()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.project_and_to2dim_ordered_indices)


        * [`Plane.projectionpoints()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.projectionpoints)


    * [`anticlockwise_sort()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.anticlockwise_sort)


    * [`anticlockwise_sort_indices()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.anticlockwise_sort_indices)


    * [`changebasis()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.changebasis)


    * [`collinear()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.collinear)


    * [`diamond_functions()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.diamond_functions)


    * [`function_comparison()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.function_comparison)


    * [`get_lower_and_upper_f()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.get_lower_and_upper_f)


    * [`is_anion_cation_bond()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.is_anion_cation_bond)


    * [`matrixTimesVector()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.matrixTimesVector)


    * [`my_solid_angle()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.my_solid_angle)


    * [`quarter_ellipsis_functions()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.quarter_ellipsis_functions)


    * [`rectangle_surface_intersection()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.rectangle_surface_intersection)


    * [`rotateCoords()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.rotateCoords)


    * [`rotateCoordsOpt()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.rotateCoordsOpt)


    * [`separation_in_list()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.separation_in_list)


    * [`sort_separation()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.sort_separation)


    * [`sort_separation_tuple()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.sort_separation_tuple)


    * [`spline_functions()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.spline_functions)


    * [`vectorsToMatrix()`](pymatgen.analysis.chemenv.utils.coordination_geometry_utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.vectorsToMatrix)


* [pymatgen.analysis.chemenv.utils.defs_utils module](pymatgen.analysis.chemenv.utils.defs_utils.md)


    * [`AdditionalConditions`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions)


        * [`AdditionalConditions.ALL`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ALL)


        * [`AdditionalConditions.CONDITION_DESCRIPTION`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.CONDITION_DESCRIPTION)


        * [`AdditionalConditions.NONE`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.NONE)


        * [`AdditionalConditions.NO_AC`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.NO_AC)


        * [`AdditionalConditions.NO_ADDITIONAL_CONDITION`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.NO_ADDITIONAL_CONDITION)


        * [`AdditionalConditions.NO_E2SEB`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.NO_E2SEB)


        * [`AdditionalConditions.NO_ELEMENT_TO_SAME_ELEMENT_BONDS`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.NO_ELEMENT_TO_SAME_ELEMENT_BONDS)


        * [`AdditionalConditions.ONLY_ACB`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_ACB)


        * [`AdditionalConditions.ONLY_ACB_AND_NO_E2SEB`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_ACB_AND_NO_E2SEB)


        * [`AdditionalConditions.ONLY_ANION_CATION_BONDS`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_ANION_CATION_BONDS)


        * [`AdditionalConditions.ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS)


        * [`AdditionalConditions.ONLY_E2OB`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_E2OB)


        * [`AdditionalConditions.ONLY_ELEMENT_TO_OXYGEN_BONDS`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_ELEMENT_TO_OXYGEN_BONDS)


        * [`AdditionalConditions.check_condition()`](pymatgen.analysis.chemenv.utils.defs_utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.check_condition)


* [pymatgen.analysis.chemenv.utils.func_utils module](pymatgen.analysis.chemenv.utils.func_utils.md)


    * [`AbstractRatioFunction`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction)


        * [`AbstractRatioFunction.ALLOWED_FUNCTIONS`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction.ALLOWED_FUNCTIONS)


        * [`AbstractRatioFunction.evaluate()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction.evaluate)


        * [`AbstractRatioFunction.from_dict()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction.from_dict)


        * [`AbstractRatioFunction.setup_parameters()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction.setup_parameters)


    * [`CSMFiniteRatioFunction`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction)


        * [`CSMFiniteRatioFunction.ALLOWED_FUNCTIONS`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.ALLOWED_FUNCTIONS)


        * [`CSMFiniteRatioFunction.fractions()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.fractions)


        * [`CSMFiniteRatioFunction.mean_estimator()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.mean_estimator)


        * [`CSMFiniteRatioFunction.power2_decreasing_exp()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.power2_decreasing_exp)


        * [`CSMFiniteRatioFunction.ratios()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.ratios)


        * [`CSMFiniteRatioFunction.smootherstep()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.smootherstep)


        * [`CSMFiniteRatioFunction.smoothstep()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.smoothstep)


    * [`CSMInfiniteRatioFunction`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction)


        * [`CSMInfiniteRatioFunction.ALLOWED_FUNCTIONS`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.ALLOWED_FUNCTIONS)


        * [`CSMInfiniteRatioFunction.fractions()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.fractions)


        * [`CSMInfiniteRatioFunction.mean_estimator()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.mean_estimator)


        * [`CSMInfiniteRatioFunction.power2_inverse_decreasing()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.power2_inverse_decreasing)


        * [`CSMInfiniteRatioFunction.power2_inverse_power2_decreasing()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.power2_inverse_power2_decreasing)


        * [`CSMInfiniteRatioFunction.ratios()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.ratios)


    * [`DeltaCSMRatioFunction`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.DeltaCSMRatioFunction)


        * [`DeltaCSMRatioFunction.ALLOWED_FUNCTIONS`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.DeltaCSMRatioFunction.ALLOWED_FUNCTIONS)


        * [`DeltaCSMRatioFunction.smootherstep()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.DeltaCSMRatioFunction.smootherstep)


    * [`RatioFunction`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction)


        * [`RatioFunction.ALLOWED_FUNCTIONS`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.ALLOWED_FUNCTIONS)


        * [`RatioFunction.inverse_smootherstep()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.inverse_smootherstep)


        * [`RatioFunction.inverse_smoothstep()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.inverse_smoothstep)


        * [`RatioFunction.power2_decreasing_exp()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.power2_decreasing_exp)


        * [`RatioFunction.power2_inverse_decreasing()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.power2_inverse_decreasing)


        * [`RatioFunction.power2_inverse_power2_decreasing()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.power2_inverse_power2_decreasing)


        * [`RatioFunction.smootherstep()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.smootherstep)


        * [`RatioFunction.smoothstep()`](pymatgen.analysis.chemenv.utils.func_utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.smoothstep)


* [pymatgen.analysis.chemenv.utils.graph_utils module](pymatgen.analysis.chemenv.utils.graph_utils.md)


    * [`MultiGraphCycle`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.MultiGraphCycle)


        * [`MultiGraphCycle.order()`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.MultiGraphCycle.order)


        * [`MultiGraphCycle.validate()`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.MultiGraphCycle.validate)


    * [`SimpleGraphCycle`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle)


        * [`SimpleGraphCycle.as_dict()`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle.as_dict)


        * [`SimpleGraphCycle.from_dict()`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle.from_dict)


        * [`SimpleGraphCycle.from_edges()`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle.from_edges)


        * [`SimpleGraphCycle.order()`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle.order)


        * [`SimpleGraphCycle.validate()`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle.validate)


    * [`get_all_elementary_cycles()`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.get_all_elementary_cycles)


    * [`get_all_simple_paths_edges()`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.get_all_simple_paths_edges)


    * [`get_delta()`](pymatgen.analysis.chemenv.utils.graph_utils.md#pymatgen.analysis.chemenv.utils.graph_utils.get_delta)


* [pymatgen.analysis.chemenv.utils.math_utils module](pymatgen.analysis.chemenv.utils.math_utils.md)


    * [`cosinus_step()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.cosinus_step)


    * [`divisors()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.divisors)


    * [`get_center_of_arc()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.get_center_of_arc)


    * [`get_linearly_independent_vectors()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.get_linearly_independent_vectors)


    * [`normal_cdf_step()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.normal_cdf_step)


    * [`power2_decreasing_exp()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.power2_decreasing_exp)


    * [`power2_inverse_decreasing()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.power2_inverse_decreasing)


    * [`power2_inverse_power2_decreasing()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.power2_inverse_power2_decreasing)


    * [`power2_inverse_powern_decreasing()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.power2_inverse_powern_decreasing)


    * [`power2_tangent_decreasing()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.power2_tangent_decreasing)


    * [`power3_step()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.power3_step)


    * [`powern_decreasing()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.powern_decreasing)


    * [`powern_parts_step()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.powern_parts_step)


    * [`prime_factors()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.prime_factors)


    * [`scale_and_clamp()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.scale_and_clamp)


    * [`smootherstep()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.smootherstep)


    * [`smoothstep()`](pymatgen.analysis.chemenv.utils.math_utils.md#pymatgen.analysis.chemenv.utils.math_utils.smoothstep)


* [pymatgen.analysis.chemenv.utils.scripts_utils module](pymatgen.analysis.chemenv.utils.scripts_utils.md)


    * [`compute_environments()`](pymatgen.analysis.chemenv.utils.scripts_utils.md#pymatgen.analysis.chemenv.utils.scripts_utils.compute_environments)


    * [`draw_cg()`](pymatgen.analysis.chemenv.utils.scripts_utils.md#pymatgen.analysis.chemenv.utils.scripts_utils.draw_cg)


    * [`visualize()`](pymatgen.analysis.chemenv.utils.scripts_utils.md#pymatgen.analysis.chemenv.utils.scripts_utils.visualize)