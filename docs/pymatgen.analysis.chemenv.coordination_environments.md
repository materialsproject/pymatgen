---
layout: default
title: pymatgen.analysis.chemenv.coordination_environments.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.coordination_environments package

Package for analyzing coordination environments.

## Subpackages


* [pymatgen.analysis.chemenv.coordination_environments.coordination_geometries_files package](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries_files.md)




* [pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies module](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md)


    * [`AbstractChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy)


        * [`AbstractChemenvStrategy.AC`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.AC)


        * [`AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE)


        * [`AbstractChemenvStrategy.STRATEGY_DESCRIPTION`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.STRATEGY_DESCRIPTION)


        * [`AbstractChemenvStrategy.STRATEGY_INFO_FIELDS`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.STRATEGY_INFO_FIELDS)


        * [`AbstractChemenvStrategy.STRATEGY_OPTIONS`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.STRATEGY_OPTIONS)


        * [`AbstractChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.as_dict)


        * [`AbstractChemenvStrategy.equivalent_site_index_and_transform()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.equivalent_site_index_and_transform)


        * [`AbstractChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.from_dict)


        * [`AbstractChemenvStrategy.get_site_ce_fractions_and_neighbors()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.get_site_ce_fractions_and_neighbors)


        * [`AbstractChemenvStrategy.get_site_coordination_environment()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.get_site_coordination_environment)


        * [`AbstractChemenvStrategy.get_site_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.get_site_coordination_environments)


        * [`AbstractChemenvStrategy.get_site_coordination_environments_fractions()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.get_site_coordination_environments_fractions)


        * [`AbstractChemenvStrategy.get_site_neighbors()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.get_site_neighbors)


        * [`AbstractChemenvStrategy.prepare_symmetries()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.prepare_symmetries)


        * [`AbstractChemenvStrategy.set_option()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.set_option)


        * [`AbstractChemenvStrategy.set_structure_environments()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.set_structure_environments)


        * [`AbstractChemenvStrategy.setup_options()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.setup_options)


        * [`AbstractChemenvStrategy.symmetry_measure_type`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.symmetry_measure_type)


        * [`AbstractChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.uniquely_determines_coordination_environments)


    * [`AdditionalConditionInt`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt)


        * [`AdditionalConditionInt.allowed_values`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt.allowed_values)


        * [`AdditionalConditionInt.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt.as_dict)


        * [`AdditionalConditionInt.description`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt.description)


        * [`AdditionalConditionInt.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt.from_dict)


        * [`AdditionalConditionInt.integer`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt.integer)


    * [`AngleCutoffFloat`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat)


        * [`AngleCutoffFloat.allowed_values`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat.allowed_values)


        * [`AngleCutoffFloat.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat.as_dict)


        * [`AngleCutoffFloat.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat.from_dict)


    * [`AngleNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight)


        * [`AngleNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.SHORT_NAME)


        * [`AngleNbSetWeight.angle_sum()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.angle_sum)


        * [`AngleNbSetWeight.angle_sumn()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.angle_sumn)


        * [`AngleNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.as_dict)


        * [`AngleNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.from_dict)


        * [`AngleNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.weight)


    * [`AnglePlateauNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight)


        * [`AnglePlateauNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight.SHORT_NAME)


        * [`AnglePlateauNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight.as_dict)


        * [`AnglePlateauNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight.from_dict)


        * [`AnglePlateauNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight.weight)


    * [`CNBiasNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight)


        * [`CNBiasNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.SHORT_NAME)


        * [`CNBiasNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.as_dict)


        * [`CNBiasNbSetWeight.explicit()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.explicit)


        * [`CNBiasNbSetWeight.from_description()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.from_description)


        * [`CNBiasNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.from_dict)


        * [`CNBiasNbSetWeight.geometrically_equidistant()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.geometrically_equidistant)


        * [`CNBiasNbSetWeight.linearly_equidistant()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.linearly_equidistant)


        * [`CNBiasNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.weight)


    * [`CSMFloat`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat)


        * [`CSMFloat.allowed_values`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat.allowed_values)


        * [`CSMFloat.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat.as_dict)


        * [`CSMFloat.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat.from_dict)


    * [`DeltaCSMNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight)


        * [`DeltaCSMNbSetWeight.DEFAULT_EFFECTIVE_CSM_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.DEFAULT_EFFECTIVE_CSM_ESTIMATOR)


        * [`DeltaCSMNbSetWeight.DEFAULT_SYMMETRY_MEASURE_TYPE`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.DEFAULT_SYMMETRY_MEASURE_TYPE)


        * [`DeltaCSMNbSetWeight.DEFAULT_WEIGHT_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.DEFAULT_WEIGHT_ESTIMATOR)


        * [`DeltaCSMNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.SHORT_NAME)


        * [`DeltaCSMNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.as_dict)


        * [`DeltaCSMNbSetWeight.delta_cn_specifics()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.delta_cn_specifics)


        * [`DeltaCSMNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.from_dict)


        * [`DeltaCSMNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.weight)


    * [`DeltaDistanceNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight)


        * [`DeltaDistanceNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight.SHORT_NAME)


        * [`DeltaDistanceNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight.as_dict)


        * [`DeltaDistanceNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight.from_dict)


        * [`DeltaDistanceNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight.weight)


    * [`DistanceAngleAreaNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight)


        * [`DistanceAngleAreaNbSetWeight.AC`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.AC)


        * [`DistanceAngleAreaNbSetWeight.DEFAULT_SURFACE_DEFINITION`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.DEFAULT_SURFACE_DEFINITION)


        * [`DistanceAngleAreaNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.SHORT_NAME)


        * [`DistanceAngleAreaNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.as_dict)


        * [`DistanceAngleAreaNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.from_dict)


        * [`DistanceAngleAreaNbSetWeight.rectangle_crosses_area()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.rectangle_crosses_area)


        * [`DistanceAngleAreaNbSetWeight.w_area_has_intersection()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.w_area_has_intersection)


        * [`DistanceAngleAreaNbSetWeight.w_area_intersection_nbsfh_fbs_onb0()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.w_area_intersection_nbsfh_fbs_onb0)


        * [`DistanceAngleAreaNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.weight)


    * [`DistanceCutoffFloat`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat)


        * [`DistanceCutoffFloat.allowed_values`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat.allowed_values)


        * [`DistanceCutoffFloat.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat.as_dict)


        * [`DistanceCutoffFloat.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat.from_dict)


    * [`DistanceNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight)


        * [`DistanceNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight.SHORT_NAME)


        * [`DistanceNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight.as_dict)


        * [`DistanceNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight.from_dict)


        * [`DistanceNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight.weight)


    * [`DistancePlateauNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight)


        * [`DistancePlateauNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight.SHORT_NAME)


        * [`DistancePlateauNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight.as_dict)


        * [`DistancePlateauNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight.from_dict)


        * [`DistancePlateauNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight.weight)


    * [`MultiWeightsChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy)


        * [`MultiWeightsChemenvStrategy.DEFAULT_CE_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.DEFAULT_CE_ESTIMATOR)


        * [`MultiWeightsChemenvStrategy.STRATEGY_DESCRIPTION`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.STRATEGY_DESCRIPTION)


        * [`MultiWeightsChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.as_dict)


        * [`MultiWeightsChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.from_dict)


        * [`MultiWeightsChemenvStrategy.stats_article_weights_parameters()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.stats_article_weights_parameters)


        * [`MultiWeightsChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.uniquely_determines_coordination_environments)


    * [`NbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NbSetWeight)


        * [`NbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NbSetWeight.as_dict)


        * [`NbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NbSetWeight.weight)


    * [`NormalizedAngleDistanceNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight)


        * [`NormalizedAngleDistanceNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.SHORT_NAME)


        * [`NormalizedAngleDistanceNbSetWeight.ang()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.ang)


        * [`NormalizedAngleDistanceNbSetWeight.anginvdist()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.anginvdist)


        * [`NormalizedAngleDistanceNbSetWeight.anginvndist()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.anginvndist)


        * [`NormalizedAngleDistanceNbSetWeight.angn()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.angn)


        * [`NormalizedAngleDistanceNbSetWeight.angninvdist()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.angninvdist)


        * [`NormalizedAngleDistanceNbSetWeight.angninvndist()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.angninvndist)


        * [`NormalizedAngleDistanceNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.as_dict)


        * [`NormalizedAngleDistanceNbSetWeight.aweight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.aweight)


        * [`NormalizedAngleDistanceNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.from_dict)


        * [`NormalizedAngleDistanceNbSetWeight.gweight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.gweight)


        * [`NormalizedAngleDistanceNbSetWeight.invdist()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.invdist)


        * [`NormalizedAngleDistanceNbSetWeight.invndist()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.invndist)


        * [`NormalizedAngleDistanceNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.weight)


    * [`SelfCSMNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight)


        * [`SelfCSMNbSetWeight.DEFAULT_EFFECTIVE_CSM_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.DEFAULT_EFFECTIVE_CSM_ESTIMATOR)


        * [`SelfCSMNbSetWeight.DEFAULT_SYMMETRY_MEASURE_TYPE`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.DEFAULT_SYMMETRY_MEASURE_TYPE)


        * [`SelfCSMNbSetWeight.DEFAULT_WEIGHT_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.DEFAULT_WEIGHT_ESTIMATOR)


        * [`SelfCSMNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.SHORT_NAME)


        * [`SelfCSMNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.as_dict)


        * [`SelfCSMNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.from_dict)


        * [`SelfCSMNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.weight)


    * [`SimpleAbundanceChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy)


        * [`SimpleAbundanceChemenvStrategy.DEFAULT_ADDITIONAL_CONDITION`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.DEFAULT_ADDITIONAL_CONDITION)


        * [`SimpleAbundanceChemenvStrategy.DEFAULT_MAX_DIST`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.DEFAULT_MAX_DIST)


        * [`SimpleAbundanceChemenvStrategy.STRATEGY_DESCRIPTION`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.STRATEGY_DESCRIPTION)


        * [`SimpleAbundanceChemenvStrategy.STRATEGY_OPTIONS`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.STRATEGY_OPTIONS)


        * [`SimpleAbundanceChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.as_dict)


        * [`SimpleAbundanceChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.from_dict)


        * [`SimpleAbundanceChemenvStrategy.get_site_coordination_environment()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.get_site_coordination_environment)


        * [`SimpleAbundanceChemenvStrategy.get_site_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.get_site_coordination_environments)


        * [`SimpleAbundanceChemenvStrategy.get_site_neighbors()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.get_site_neighbors)


        * [`SimpleAbundanceChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.uniquely_determines_coordination_environments)


    * [`SimplestChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy)


        * [`SimplestChemenvStrategy.DEFAULT_ADDITIONAL_CONDITION`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.DEFAULT_ADDITIONAL_CONDITION)


        * [`SimplestChemenvStrategy.DEFAULT_ANGLE_CUTOFF`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.DEFAULT_ANGLE_CUTOFF)


        * [`SimplestChemenvStrategy.DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF)


        * [`SimplestChemenvStrategy.DEFAULT_DISTANCE_CUTOFF`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.DEFAULT_DISTANCE_CUTOFF)


        * [`SimplestChemenvStrategy.STRATEGY_DESCRIPTION`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.STRATEGY_DESCRIPTION)


        * [`SimplestChemenvStrategy.STRATEGY_OPTIONS`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.STRATEGY_OPTIONS)


        * [`SimplestChemenvStrategy.add_strategy_visualization_to_subplot()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.add_strategy_visualization_to_subplot)


        * [`SimplestChemenvStrategy.additional_condition`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.additional_condition)


        * [`SimplestChemenvStrategy.angle_cutoff`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.angle_cutoff)


        * [`SimplestChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.as_dict)


        * [`SimplestChemenvStrategy.continuous_symmetry_measure_cutoff`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.continuous_symmetry_measure_cutoff)


        * [`SimplestChemenvStrategy.distance_cutoff`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.distance_cutoff)


        * [`SimplestChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.from_dict)


        * [`SimplestChemenvStrategy.get_site_coordination_environment()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.get_site_coordination_environment)


        * [`SimplestChemenvStrategy.get_site_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.get_site_coordination_environments)


        * [`SimplestChemenvStrategy.get_site_coordination_environments_fractions()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.get_site_coordination_environments_fractions)


        * [`SimplestChemenvStrategy.get_site_neighbors()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.get_site_neighbors)


        * [`SimplestChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.uniquely_determines_coordination_environments)


    * [`StrategyOption`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.StrategyOption)


        * [`StrategyOption.allowed_values`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.StrategyOption.allowed_values)


        * [`StrategyOption.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.StrategyOption.as_dict)


    * [`TargettedPenaltiedAbundanceChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy)


        * [`TargettedPenaltiedAbundanceChemenvStrategy.DEFAULT_TARGET_ENVIRONMENTS`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy.DEFAULT_TARGET_ENVIRONMENTS)


        * [`TargettedPenaltiedAbundanceChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy.as_dict)


        * [`TargettedPenaltiedAbundanceChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy.from_dict)


        * [`TargettedPenaltiedAbundanceChemenvStrategy.get_site_coordination_environment()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy.get_site_coordination_environment)


        * [`TargettedPenaltiedAbundanceChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy.uniquely_determines_coordination_environments)


    * [`WeightedNbSetChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy)


        * [`WeightedNbSetChemenvStrategy.DEFAULT_CE_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.DEFAULT_CE_ESTIMATOR)


        * [`WeightedNbSetChemenvStrategy.STRATEGY_DESCRIPTION`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.STRATEGY_DESCRIPTION)


        * [`WeightedNbSetChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.as_dict)


        * [`WeightedNbSetChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.from_dict)


        * [`WeightedNbSetChemenvStrategy.get_site_coordination_environment()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.get_site_coordination_environment)


        * [`WeightedNbSetChemenvStrategy.get_site_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.get_site_coordination_environments)


        * [`WeightedNbSetChemenvStrategy.get_site_coordination_environments_fractions()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.get_site_coordination_environments_fractions)


        * [`WeightedNbSetChemenvStrategy.get_site_neighbors()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.get_site_neighbors)


        * [`WeightedNbSetChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.uniquely_determines_coordination_environments)


    * [`get_effective_csm()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.get_effective_csm)


    * [`set_info()`](pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.set_info)


* [pymatgen.analysis.chemenv.coordination_environments.coordination_geometries module](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md)


    * [`AbstractChemenvAlgorithm`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AbstractChemenvAlgorithm)


        * [`AbstractChemenvAlgorithm.algorithm_type`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AbstractChemenvAlgorithm.algorithm_type)


        * [`AbstractChemenvAlgorithm.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AbstractChemenvAlgorithm.as_dict)


    * [`AllCoordinationGeometries`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries)


        * [`AllCoordinationGeometries.get_geometries()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_geometries)


        * [`AllCoordinationGeometries.get_geometry_from_IUCr_symbol()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_geometry_from_IUCr_symbol)


        * [`AllCoordinationGeometries.get_geometry_from_IUPAC_symbol()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_geometry_from_IUPAC_symbol)


        * [`AllCoordinationGeometries.get_geometry_from_mp_symbol()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_geometry_from_mp_symbol)


        * [`AllCoordinationGeometries.get_geometry_from_name()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_geometry_from_name)


        * [`AllCoordinationGeometries.get_implemented_geometries()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_implemented_geometries)


        * [`AllCoordinationGeometries.get_not_implemented_geometries()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_not_implemented_geometries)


        * [`AllCoordinationGeometries.get_symbol_cn_mapping()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_symbol_cn_mapping)


        * [`AllCoordinationGeometries.get_symbol_name_mapping()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_symbol_name_mapping)


        * [`AllCoordinationGeometries.is_a_valid_coordination_geometry()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.is_a_valid_coordination_geometry)


        * [`AllCoordinationGeometries.pretty_print()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.pretty_print)


    * [`CoordinationGeometry`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry)


        * [`CoordinationGeometry.CSM_SKIP_SEPARATION_PLANE_ALGO`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.CSM_SKIP_SEPARATION_PLANE_ALGO)


        * [`CoordinationGeometry.IUCr_symbol`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.IUCr_symbol)


        * [`CoordinationGeometry.IUCr_symbol_str`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.IUCr_symbol_str)


        * [`CoordinationGeometry.IUPAC_symbol`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.IUPAC_symbol)


        * [`CoordinationGeometry.IUPAC_symbol_str`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.IUPAC_symbol_str)


        * [`CoordinationGeometry.NeighborsSetsHints`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints)


            * [`CoordinationGeometry.NeighborsSetsHints.ALLOWED_HINTS_TYPES`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.ALLOWED_HINTS_TYPES)


            * [`CoordinationGeometry.NeighborsSetsHints.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.as_dict)


            * [`CoordinationGeometry.NeighborsSetsHints.double_cap_hints()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.double_cap_hints)


            * [`CoordinationGeometry.NeighborsSetsHints.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.from_dict)


            * [`CoordinationGeometry.NeighborsSetsHints.hints()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.hints)


            * [`CoordinationGeometry.NeighborsSetsHints.single_cap_hints()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.single_cap_hints)


            * [`CoordinationGeometry.NeighborsSetsHints.triple_cap_hints()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.triple_cap_hints)


        * [`CoordinationGeometry.algorithms`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.algorithms)


        * [`CoordinationGeometry.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.as_dict)


        * [`CoordinationGeometry.ce_symbol`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.ce_symbol)


        * [`CoordinationGeometry.coordination_number`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.coordination_number)


        * [`CoordinationGeometry.distfactor_max`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.distfactor_max)


        * [`CoordinationGeometry.edges()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.edges)


        * [`CoordinationGeometry.faces()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.faces)


        * [`CoordinationGeometry.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.from_dict)


        * [`CoordinationGeometry.get_central_site()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.get_central_site)


        * [`CoordinationGeometry.get_coordination_number()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.get_coordination_number)


        * [`CoordinationGeometry.get_name()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.get_name)


        * [`CoordinationGeometry.get_pmeshes()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.get_pmeshes)


        * [`CoordinationGeometry.is_implemented()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.is_implemented)


        * [`CoordinationGeometry.mp_symbol`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.mp_symbol)


        * [`CoordinationGeometry.number_of_permutations`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.number_of_permutations)


        * [`CoordinationGeometry.pauling_stability_ratio`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.pauling_stability_ratio)


        * [`CoordinationGeometry.ref_permutation()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.ref_permutation)


        * [`CoordinationGeometry.set_permutations_safe_override()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.set_permutations_safe_override)


        * [`CoordinationGeometry.solid_angles()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.solid_angles)


    * [`ExplicitPermutationsAlgorithm`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.ExplicitPermutationsAlgorithm)


        * [`ExplicitPermutationsAlgorithm.as_dict`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.ExplicitPermutationsAlgorithm.as_dict)


        * [`ExplicitPermutationsAlgorithm.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.ExplicitPermutationsAlgorithm.from_dict)


        * [`ExplicitPermutationsAlgorithm.permutations`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.ExplicitPermutationsAlgorithm.permutations)


    * [`SeparationPlane`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane)


        * [`SeparationPlane.argsorted_ref_separation_perm`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.argsorted_ref_separation_perm)


        * [`SeparationPlane.as_dict`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.as_dict)


        * [`SeparationPlane.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.from_dict)


        * [`SeparationPlane.permutations`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.permutations)


        * [`SeparationPlane.ref_separation_perm`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.ref_separation_perm)


        * [`SeparationPlane.safe_separation_permutations()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.safe_separation_permutations)


* [pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder module](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md)


    * [`AbstractGeometry`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry)


        * [`AbstractGeometry.cn`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.cn)


        * [`AbstractGeometry.coordination_number`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.coordination_number)


        * [`AbstractGeometry.from_cg()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.from_cg)


        * [`AbstractGeometry.points_wcs_csc()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wcs_csc)


        * [`AbstractGeometry.points_wcs_ctwcc()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wcs_ctwcc)


        * [`AbstractGeometry.points_wcs_ctwocc()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wcs_ctwocc)


        * [`AbstractGeometry.points_wocs_csc()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wocs_csc)


        * [`AbstractGeometry.points_wocs_ctwcc()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wocs_ctwcc)


        * [`AbstractGeometry.points_wocs_ctwocc()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wocs_ctwocc)


    * [`LocalGeometryFinder`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder)


        * [`LocalGeometryFinder.BVA_DISTANCE_SCALE_FACTORS`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.BVA_DISTANCE_SCALE_FACTORS)


        * [`LocalGeometryFinder.DEFAULT_BVA_DISTANCE_SCALE_FACTOR`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.DEFAULT_BVA_DISTANCE_SCALE_FACTOR)


        * [`LocalGeometryFinder.DEFAULT_SPG_ANALYZER_OPTIONS`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.DEFAULT_SPG_ANALYZER_OPTIONS)


        * [`LocalGeometryFinder.DEFAULT_STRATEGY`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.DEFAULT_STRATEGY)


        * [`LocalGeometryFinder.PRESETS`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.PRESETS)


        * [`LocalGeometryFinder.STRUCTURE_REFINEMENT_NONE`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.STRUCTURE_REFINEMENT_NONE)


        * [`LocalGeometryFinder.STRUCTURE_REFINEMENT_REFINED`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.STRUCTURE_REFINEMENT_REFINED)


        * [`LocalGeometryFinder.STRUCTURE_REFINEMENT_SYMMETRIZED`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.STRUCTURE_REFINEMENT_SYMMETRIZED)


        * [`LocalGeometryFinder.compute_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.compute_coordination_environments)


        * [`LocalGeometryFinder.compute_structure_environments()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.compute_structure_environments)


        * [`LocalGeometryFinder.coordination_geometry_symmetry_measures()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures)


        * [`LocalGeometryFinder.coordination_geometry_symmetry_measures_fallback_random()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures_fallback_random)


        * [`LocalGeometryFinder.coordination_geometry_symmetry_measures_separation_plane()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures_separation_plane)


        * [`LocalGeometryFinder.coordination_geometry_symmetry_measures_separation_plane_optim()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures_separation_plane_optim)


        * [`LocalGeometryFinder.coordination_geometry_symmetry_measures_sepplane_optim()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures_sepplane_optim)


        * [`LocalGeometryFinder.coordination_geometry_symmetry_measures_standard()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures_standard)


        * [`LocalGeometryFinder.get_coordination_symmetry_measures()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.get_coordination_symmetry_measures)


        * [`LocalGeometryFinder.get_coordination_symmetry_measures_optim()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.get_coordination_symmetry_measures_optim)


        * [`LocalGeometryFinder.get_structure()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.get_structure)


        * [`LocalGeometryFinder.set_structure()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.set_structure)


        * [`LocalGeometryFinder.setup_explicit_indices_local_geometry()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_explicit_indices_local_geometry)


        * [`LocalGeometryFinder.setup_local_geometry()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_local_geometry)


        * [`LocalGeometryFinder.setup_ordered_indices_local_geometry()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_ordered_indices_local_geometry)


        * [`LocalGeometryFinder.setup_parameter()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_parameter)


        * [`LocalGeometryFinder.setup_parameters()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_parameters)


        * [`LocalGeometryFinder.setup_random_indices_local_geometry()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_random_indices_local_geometry)


        * [`LocalGeometryFinder.setup_random_structure()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_random_structure)


        * [`LocalGeometryFinder.setup_structure()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_structure)


        * [`LocalGeometryFinder.setup_test_perfect_environment()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_test_perfect_environment)


        * [`LocalGeometryFinder.update_nb_set_environments()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.update_nb_set_environments)


    * [`find_rotation()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.find_rotation)


    * [`find_scaling_factor()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.find_scaling_factor)


    * [`symmetry_measure()`](pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.symmetry_measure)


* [pymatgen.analysis.chemenv.coordination_environments.structure_environments module](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md)


    * [`ChemicalEnvironments`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments)


        * [`ChemicalEnvironments.add_coord_geom()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.add_coord_geom)


        * [`ChemicalEnvironments.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.as_dict)


        * [`ChemicalEnvironments.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.from_dict)


        * [`ChemicalEnvironments.is_close_to()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.is_close_to)


        * [`ChemicalEnvironments.minimum_geometries()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.minimum_geometries)


        * [`ChemicalEnvironments.minimum_geometry()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.minimum_geometry)


    * [`LightStructureEnvironments`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments)


        * [`LightStructureEnvironments.DEFAULT_STATISTICS_FIELDS`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.DEFAULT_STATISTICS_FIELDS)


        * [`LightStructureEnvironments.DELTA_MAX_OXIDATION_STATE`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.DELTA_MAX_OXIDATION_STATE)


        * [`LightStructureEnvironments.NeighborsSet`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet)


            * [`LightStructureEnvironments.NeighborsSet.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.as_dict)


            * [`LightStructureEnvironments.NeighborsSet.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.from_dict)


            * [`LightStructureEnvironments.NeighborsSet.neighb_coords`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.neighb_coords)


            * [`LightStructureEnvironments.NeighborsSet.neighb_indices_and_images`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.neighb_indices_and_images)


            * [`LightStructureEnvironments.NeighborsSet.neighb_sites`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.neighb_sites)


            * [`LightStructureEnvironments.NeighborsSet.neighb_sites_and_indices`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.neighb_sites_and_indices)


        * [`LightStructureEnvironments.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.as_dict)


        * [`LightStructureEnvironments.clear_environments()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.clear_environments)


        * [`LightStructureEnvironments.contains_only_one_anion()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.contains_only_one_anion)


        * [`LightStructureEnvironments.contains_only_one_anion_atom()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.contains_only_one_anion_atom)


        * [`LightStructureEnvironments.environments_identified()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.environments_identified)


        * [`LightStructureEnvironments.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.from_dict)


        * [`LightStructureEnvironments.from_structure_environments()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.from_structure_environments)


        * [`LightStructureEnvironments.get_site_info_for_specie_allces()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.get_site_info_for_specie_allces)


        * [`LightStructureEnvironments.get_site_info_for_specie_ce()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.get_site_info_for_specie_ce)


        * [`LightStructureEnvironments.get_statistics()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.get_statistics)


        * [`LightStructureEnvironments.setup_statistic_lists()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.setup_statistic_lists)


        * [`LightStructureEnvironments.site_contains_environment()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.site_contains_environment)


        * [`LightStructureEnvironments.site_has_clear_environment()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.site_has_clear_environment)


        * [`LightStructureEnvironments.structure_contains_atom_environment()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.structure_contains_atom_environment)


        * [`LightStructureEnvironments.structure_has_clear_environments()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.structure_has_clear_environments)


        * [`LightStructureEnvironments.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.uniquely_determines_coordination_environments)


    * [`StructureEnvironments`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments)


        * [`StructureEnvironments.AC`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.AC)


        * [`StructureEnvironments.NeighborsSet`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet)


            * [`StructureEnvironments.NeighborsSet.add_source()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.add_source)


            * [`StructureEnvironments.NeighborsSet.angle_plateau()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.angle_plateau)


            * [`StructureEnvironments.NeighborsSet.angles`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.angles)


            * [`StructureEnvironments.NeighborsSet.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.as_dict)


            * [`StructureEnvironments.NeighborsSet.coords`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.coords)


            * [`StructureEnvironments.NeighborsSet.distance_plateau()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.distance_plateau)


            * [`StructureEnvironments.NeighborsSet.distances`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.distances)


            * [`StructureEnvironments.NeighborsSet.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.from_dict)


            * [`StructureEnvironments.NeighborsSet.get_neighb_voronoi_indices()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.get_neighb_voronoi_indices)


            * [`StructureEnvironments.NeighborsSet.info`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.info)


            * [`StructureEnvironments.NeighborsSet.neighb_coords`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.neighb_coords)


            * [`StructureEnvironments.NeighborsSet.neighb_coordsOpt`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.neighb_coordsOpt)


            * [`StructureEnvironments.NeighborsSet.neighb_sites`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.neighb_sites)


            * [`StructureEnvironments.NeighborsSet.neighb_sites_and_indices`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.neighb_sites_and_indices)


            * [`StructureEnvironments.NeighborsSet.normalized_angles`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.normalized_angles)


            * [`StructureEnvironments.NeighborsSet.normalized_distances`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.normalized_distances)


            * [`StructureEnvironments.NeighborsSet.source`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.source)


            * [`StructureEnvironments.NeighborsSet.voronoi_grid_surface_points()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.voronoi_grid_surface_points)


        * [`StructureEnvironments.add_neighbors_set()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.add_neighbors_set)


        * [`StructureEnvironments.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.as_dict)


        * [`StructureEnvironments.differences_wrt()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.differences_wrt)


        * [`StructureEnvironments.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.from_dict)


        * [`StructureEnvironments.get_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.get_coordination_environments)


        * [`StructureEnvironments.get_csm()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.get_csm)


        * [`StructureEnvironments.get_csm_and_maps()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.get_csm_and_maps)


        * [`StructureEnvironments.get_csms()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.get_csms)


        * [`StructureEnvironments.get_environments_figure()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.get_environments_figure)


        * [`StructureEnvironments.init_neighbors_sets()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.init_neighbors_sets)


        * [`StructureEnvironments.plot_csm_and_maps()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.plot_csm_and_maps)


        * [`StructureEnvironments.plot_environments()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.plot_environments)


        * [`StructureEnvironments.save_environments_figure()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.save_environments_figure)


        * [`StructureEnvironments.update_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.update_coordination_environments)


        * [`StructureEnvironments.update_site_info()`](pymatgen.analysis.chemenv.coordination_environments.structure_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.update_site_info)


* [pymatgen.analysis.chemenv.coordination_environments.voronoi module](pymatgen.analysis.chemenv.coordination_environments.voronoi.md)


    * [`DetailedVoronoiContainer`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer)


        * [`DetailedVoronoiContainer.AC`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.AC)


        * [`DetailedVoronoiContainer.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.as_dict)


        * [`DetailedVoronoiContainer.default_normalized_angle_tolerance`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.default_normalized_angle_tolerance)


        * [`DetailedVoronoiContainer.default_normalized_distance_tolerance`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.default_normalized_distance_tolerance)


        * [`DetailedVoronoiContainer.default_voronoi_cutoff`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.default_voronoi_cutoff)


        * [`DetailedVoronoiContainer.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.from_dict)


        * [`DetailedVoronoiContainer.get_rdf_figure()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.get_rdf_figure)


        * [`DetailedVoronoiContainer.get_sadf_figure()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.get_sadf_figure)


        * [`DetailedVoronoiContainer.is_close_to()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.is_close_to)


        * [`DetailedVoronoiContainer.maps_and_surfaces()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.maps_and_surfaces)


        * [`DetailedVoronoiContainer.maps_and_surfaces_bounded()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.maps_and_surfaces_bounded)


        * [`DetailedVoronoiContainer.neighbors()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.neighbors)


        * [`DetailedVoronoiContainer.neighbors_surfaces()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.neighbors_surfaces)


        * [`DetailedVoronoiContainer.neighbors_surfaces_bounded()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.neighbors_surfaces_bounded)


        * [`DetailedVoronoiContainer.setup_neighbors_distances_and_angles()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.setup_neighbors_distances_and_angles)


        * [`DetailedVoronoiContainer.setup_voronoi_list()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.setup_voronoi_list)


        * [`DetailedVoronoiContainer.to_bson_voronoi_list2()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.to_bson_voronoi_list2)


        * [`DetailedVoronoiContainer.voronoi_parameters_bounds_and_limits()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.voronoi_parameters_bounds_and_limits)


    * [`from_bson_voronoi_list2()`](pymatgen.analysis.chemenv.coordination_environments.voronoi.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.from_bson_voronoi_list2)