---
layout: default
title: pymatgen.analysis.chemenv.coordination_environments.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.chemenv.coordination_environments package

Package for analyzing coordination environments.

## Subpackages


* [pymatgen.analysis.chemenv.coordination_environments.coordination_geometries_files package](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries_files.md)



## pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies module

This module provides so-called “strategies” to determine the coordination environments of an atom in a structure.
Some strategies can favour larger or smaller environments. Some strategies uniquely identifies the environments while
some others can identify the environment as a “mix” of several environments, each of which is assigned with a given
fraction. The choice of the strategy depends on the purpose of the user.


### _class_ AbstractChemenvStrategy(structure_environments=None, symmetry_measure_type='csm_wcs_ctwcc')
Bases: `MSONable`

Class used to define a Chemenv strategy for the neighbors and coordination environment to be applied to a
StructureEnvironments object.

Abstract constructor for the all chemenv strategies.
:param structure_environments: StructureEnvironments object containing all the information on the

> coordination of the sites in a structure.


#### AC(_ = <pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions object_ )

#### DEFAULT_SYMMETRY_MEASURE_TYPE(_ = 'csm_wcs_ctwcc_ )

#### STRATEGY_DESCRIPTION(_: str | Non_ _ = Non_ )

#### STRATEGY_INFO_FIELDS(_: ClassVar[list_ _ = [_ )

#### STRATEGY_OPTIONS(_: ClassVar[dict[str, dict]_ _ = {_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract_ as_dict()
Bson-serializable dict representation of the SimplestChemenvStrategy object.


* **Returns**

    Bson-serializable dict representation of the SimplestChemenvStrategy object.



#### equivalent_site_index_and_transform(psite)
Get the equivalent site and corresponding symmetry+translation transformations.


* **Parameters**

    **psite** – Periodic site.



* **Returns**

    Equivalent site in the unit cell, translations and symmetry transformation.



#### _classmethod_ from_dict(d)
Reconstructs the SimpleAbundanceChemenvStrategy object from a dict representation of the
SimpleAbundanceChemenvStrategy object created using the as_dict method.
:param d: dict representation of the SimpleAbundanceChemenvStrategy object


* **Returns**

    StructureEnvironments object.



#### get_site_ce_fractions_and_neighbors(site, full_ce_info=False, strategy_info=False)
Applies the strategy to the structure_environments object in order to get coordination environments, their
fraction, csm, geometry_info, and neighbors
:param site: Site for which the above information is sought


* **Returns**

    The list of neighbors of the site. For complex strategies, where one allows multiple solutions, this


can return a list of list of neighbors.


#### _abstract_ get_site_coordination_environment(site)
Applies the strategy to the structure_environments object in order to define the coordination environment of
a given site.
:param site: Site for which the coordination environment is looked for


* **Returns**

    The coordination environment of the site. For complex strategies, where one allows multiple
    solutions, this can return a list of coordination environments for the site.



#### _abstract_ get_site_coordination_environments(site)
Applies the strategy to the structure_environments object in order to define the coordination environment of
a given site.
:param site: Site for which the coordination environment is looked for


* **Returns**

    The coordination environment of the site. For complex strategies, where one allows multiple
    solutions, this can return a list of coordination environments for the site.



#### _abstract_ get_site_coordination_environments_fractions(site, isite=None, dequivsite=None, dthissite=None, mysym=None, ordered=True, min_fraction=0, return_maps=True, return_strategy_dict_info=False)
Applies the strategy to the structure_environments object in order to define the coordination environment of
a given site.
:param site: Site for which the coordination environment is looked for


* **Returns**

    The coordination environment of the site. For complex strategies, where one allows multiple
    solutions, this can return a list of coordination environments for the site.



#### _abstract_ get_site_neighbors(site)
Applies the strategy to the structure_environments object in order to get the neighbors of a given site.
:param site: Site for which the neighbors are looked for
:param structure_environments: StructureEnvironments object containing all the information needed to get the

> neighbors of the site


* **Returns**

    The list of neighbors of the site. For complex strategies, where one allows multiple solutions, this
    can return a list of list of neighbors.



#### prepare_symmetries()
Prepare the symmetries for the structure contained in the structure environments.


#### set_option(option_name, option_value)
Set up a given option for this strategy.


* **Parameters**


    * **option_name** – Name of the option.


    * **option_value** – Value for this option.



#### set_structure_environments(structure_environments)
Set the structure environments to this strategy.


* **Parameters**

    **structure_environments** – StructureEnvironments object.



#### setup_options(all_options_dict)
Set up options for this strategy based on a dict.


* **Parameters**

    **all_options_dict** – Dict of option_name->option_value.



#### _property_ symmetry_measure_type()
Type of symmetry measure.


#### _property_ uniquely_determines_coordination_environments()
Returns True if the strategy leads to a unique coordination environment.


### _class_ AdditionalConditionInt(integer)
Bases: `int`, `StrategyOption`

Integer representing an additional condition in a strategy.

Special int representing additional conditions.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### allowed_values(_: str | Non_ _ = "Integer amongst :\\n - 0 for 'No additional condition'\\n - 1 for 'Only anion-cation bonds'\\n - 2 for 'No element-element bonds (same elements)'\\n - 3 for 'Only anion-cation bonds and no element-element bonds (same elements)'\\n - 4 for 'Only element-oxygen bonds'\\n_ )

#### as_dict()
MSONable dict.


#### description(_ = 'Only element-oxygen bonds_ )

#### _classmethod_ from_dict(dct)
Initialize additional condition from dict.


* **Parameters**

    **d** – Dict representation of the additional condition.



#### integer(_ = _ )

### _class_ AngleCutoffFloat(cutoff)
Bases: `float`, `StrategyOption`

Angle cutoff in a strategy.

Special float that should be between 0 and 1.


* **Parameters**

    **cutoff** – Angle cutoff.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### allowed_values(_: str | Non_ _ = 'Real number between 0 and 1_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(d)
Initialize angle cutoff from dict.


* **Parameters**

    **d** – Dict representation of the angle cutoff.



### _class_ AngleNbSetWeight(aa=1)
Bases: `NbSetWeight`

Weight of neighbors set based on the angle.

Initialize AngleNbSetWeight estimator.


* **Parameters**

    **aa** – Exponent of the angle for the estimator.



#### SHORT_NAME(_ = 'AngleWeight_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ angle_sum(nb_set)
Sum of all angles in a neighbors set.


* **Parameters**

    **nb_set** – Neighbors set.



* **Returns**

    Sum of solid angles for the neighbors set.



#### angle_sumn(nb_set)
Sum of all angles to a given power in a neighbors set.


* **Parameters**

    **nb_set** – Neighbors set.



* **Returns**

    Sum of solid angles to the power aa for the neighbors set.



#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dct)
Construct AngleNbSetWeight from dict representation.


#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ AnglePlateauNbSetWeight(angle_function=None, weight_function=None)
Bases: `NbSetWeight`

Weight of neighbors set based on the angle.

Initialize AnglePlateauNbSetWeight.


* **Parameters**


    * **angle_function** – Angle function to use.


    * **weight_function** – Ratio function to use.



#### SHORT_NAME(_ = 'AnglePlateauWeight_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dct)
Initialize from dict.


* **Parameters**

    **dct** – Dict representation of AnglePlateauNbSetWeight.



* **Returns**

    AnglePlateauNbSetWeight.



#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ CNBiasNbSetWeight(cn_weights, initialization_options)
Bases: `NbSetWeight`

Weight of neighbors set based on specific biases towards specific coordination numbers.

Initialize CNBiasNbSetWeight.


* **Parameters**


    * **cn_weights** – Weights for each coordination.


    * **initialization_options** – Options for initialization.



#### SHORT_NAME(_ = 'CNBiasWeight_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict.


#### _classmethod_ explicit(cn_weights)
Initializes weights explicitly for each coordination.


* **Parameters**

    **cn_weights** – Weights for each coordination.



* **Returns**

    CNBiasNbSetWeight.



#### _classmethod_ from_description(dct)
Initializes weights from description.


* **Parameters**

    **dct** – Dictionary description.



* **Returns**

    CNBiasNbSetWeight.



#### _classmethod_ from_dict(dct)
Initialize from dict.


* **Parameters**

    **dct** – Dict representation of CNBiasNbSetWeight.



* **Returns**

    CNBiasNbSetWeight.



#### _classmethod_ geometrically_equidistant(weight_cn1, weight_cn13)
Initializes geometrically equidistant weights for each coordination.


* **Parameters**


    * **weight_cn1** – Weight of coordination 1.


    * **weight_cn13** – Weight of coordination 13.



* **Returns**

    CNBiasNbSetWeight.



#### _classmethod_ linearly_equidistant(weight_cn1, weight_cn13)
Initializes linearly equidistant weights for each coordination.


* **Parameters**


    * **weight_cn1** – Weight of coordination 1.


    * **weight_cn13** – Weight of coordination 13.



* **Returns**

    CNBiasNbSetWeight.



#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ CSMFloat(cutoff)
Bases: `float`, `StrategyOption`

Real number representing a Continuous Symmetry Measure.

Special float that should be between 0 and 100.


* **Parameters**

    **cutoff** – CSM.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### allowed_values(_: str | Non_ _ = 'Real number between 0 and 100_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dct)
Initialize CSM from dict.


* **Parameters**

    **d** – Dict representation of the CSM.



### _class_ DeltaCSMNbSetWeight(effective_csm_estimator={'function': 'power2_inverse_decreasing', 'options': {'max_csm': 8.0}}, weight_estimator={'function': 'smootherstep', 'options': {'delta_csm_max': 3.0, 'delta_csm_min': 0.5}}, delta_cn_weight_estimators=None, symmetry_measure_type='csm_wcs_ctwcc')
Bases: `NbSetWeight`

Weight of neighbors set based on the differences of CSM.

Initialize SelfCSMNbSetWeight.


* **Parameters**


    * **effective_csm_estimator** – Ratio function used for the effective CSM (comparison between neighbors sets).


    * **weight_estimator** – Weight estimator within a given neighbors set.


    * **delta_cn_weight_estimators** – Specific weight estimators for specific cn


    * **symmetry_measure_type** – Type of symmetry measure to be used.



#### DEFAULT_EFFECTIVE_CSM_ESTIMATOR(_ = {'function': 'power2_inverse_decreasing', 'options': {'max_csm': 8.0}_ )

#### DEFAULT_SYMMETRY_MEASURE_TYPE(_ = 'csm_wcs_ctwcc_ )

#### DEFAULT_WEIGHT_ESTIMATOR(_ = {'function': 'smootherstep', 'options': {'delta_csm_max': 3.0, 'delta_csm_min': 0.5}_ )

#### SHORT_NAME(_ = 'DeltaCSMWeight_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict.


#### _classmethod_ delta_cn_specifics(delta_csm_mins=None, delta_csm_maxs=None, function='smootherstep', symmetry_measure_type='csm_wcs_ctwcc', effective_csm_estimator={'function': 'power2_inverse_decreasing', 'options': {'max_csm': 8.0}})
Initializes DeltaCSMNbSetWeight from specific coordination number differences.


* **Parameters**


    * **delta_csm_mins** – Minimums for each coordination number.


    * **delta_csm_maxs** – Maximums for each coordination number.


    * **function** – Ratio function used.


    * **symmetry_measure_type** – Type of symmetry measure to be used.


    * **effective_csm_estimator** – Ratio function used for the effective CSM (comparison between neighbors sets).



* **Returns**

    DeltaCSMNbSetWeight.



#### _classmethod_ from_dict(dct)
Initialize from dict.


* **Parameters**

    **dct** – Dict representation of DeltaCSMNbSetWeight.



* **Returns**

    DeltaCSMNbSetWeight.



#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ DeltaDistanceNbSetWeight(weight_function=None, nbs_source='voronoi')
Bases: `NbSetWeight`

Weight of neighbors set based on the difference of distances.

Initialize DeltaDistanceNbSetWeight.


* **Parameters**


    * **weight_function** – Ratio function to use.


    * **nbs_source** – Source of the neighbors.



#### SHORT_NAME(_ = 'DeltaDistanceNbSetWeight_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dct)
Initialize from dict.


* **Parameters**

    **dct** – Dict representation of DeltaDistanceNbSetWeight.



* **Returns**

    DeltaDistanceNbSetWeight.



#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ DistanceAngleAreaNbSetWeight(weight_type='has_intersection', surface_definition={'angle_bounds': {'lower': 0.1, 'upper': 0.8}, 'distance_bounds': {'lower': 1.2, 'upper': 1.8}, 'type': 'standard_elliptic'}, nb_sets_from_hints='fallback_to_source', other_nb_sets='0_weight', additional_condition=1, smoothstep_distance=None, smoothstep_angle=None)
Bases: `NbSetWeight`

Weight of neighbors set based on the area in the distance-angle space.

Initialize CNBiasNbSetWeight.


* **Parameters**


    * **weight_type** – Type of weight.


    * **surface_definition** – Definition of the surface.


    * **nb_sets_from_hints** – How to deal with neighbors sets obtained from “hints”.


    * **other_nb_sets** – What to do with other neighbors sets.


    * **additional_condition** – Additional condition to be used.


    * **smoothstep_distance** – Smoothstep distance.


    * **smoothstep_angle** – Smoothstep angle.



#### AC(_ = <pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions object_ )

#### DEFAULT_SURFACE_DEFINITION(_ = {'angle_bounds': {'lower': 0.1, 'upper': 0.8}, 'distance_bounds': {'lower': 1.2, 'upper': 1.8}, 'type': 'standard_elliptic'_ )

#### SHORT_NAME(_ = 'DistAngleAreaWeight_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dct)
Initialize from dict.


* **Parameters**

    **dct** – Dict representation of DistanceAngleAreaNbSetWeight.



* **Returns**

    DistanceAngleAreaNbSetWeight.



#### rectangle_crosses_area(d1, d2, a1, a2)
Whether a given rectangle crosses the area defined by the upper and lower curves.


* **Parameters**


    * **d1** – lower d.


    * **d2** – upper d.


    * **a1** – lower a.


    * **a2** – upper a.



#### w_area_has_intersection(nb_set, structure_environments, cn_map, additional_info)
Get intersection of the neighbors set area with the surface.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments.


    * **cn_map** – Mapping index of the neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Area intersection between neighbors set and surface.



#### w_area_intersection_nbsfh_fbs_onb0(nb_set, structure_environments, cn_map, additional_info)
Get intersection of the neighbors set area with the surface.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments.


    * **cn_map** – Mapping index of the neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Area intersection between neighbors set and surface.



#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ DistanceCutoffFloat(cutoff)
Bases: `float`, `StrategyOption`

Distance cutoff in a strategy.

Special float that should be between 1 and infinity.


* **Parameters**

    **cutoff** – Distance cutoff.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### allowed_values(_: str | Non_ _ = 'Real number between 1 and +infinity_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(d)
Initialize distance cutoff from dict.


* **Parameters**

    **d** – Dict representation of the distance cutoff.



### _class_ DistanceNbSetWeight(weight_function=None, nbs_source='voronoi')
Bases: `NbSetWeight`

Weight of neighbors set based on the distance.

Initialize DistanceNbSetWeight.


* **Parameters**


    * **weight_function** – Ratio function to use.


    * **nbs_source** – Source of the neighbors.



#### SHORT_NAME(_ = 'DistanceNbSetWeight_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSOnable dict.


#### _classmethod_ from_dict(dct)
Initialize from dict.


* **Parameters**

    **dct** – Dict representation of DistanceNbSetWeight.



* **Returns**

    DistanceNbSetWeight.



#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ DistancePlateauNbSetWeight(distance_function=None, weight_function=None)
Bases: `NbSetWeight`

Weight of neighbors set based on the distance.

Initialize DistancePlateauNbSetWeight.


* **Parameters**


    * **distance_function** – Distance function to use.


    * **weight_function** – Ratio function to use.



#### SHORT_NAME(_ = 'DistancePlateauWeight_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dct)
Initialize from dict.


* **Parameters**

    **dct** – Dict representation of DistancePlateauNbSetWeight.



* **Returns**

    DistancePlateauNbSetWeight.



#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ MultiWeightsChemenvStrategy(structure_environments=None, additional_condition=1, symmetry_measure_type='csm_wcs_ctwcc', dist_ang_area_weight=None, self_csm_weight=None, delta_csm_weight=None, cn_bias_weight=None, angle_weight=None, normalized_angle_distance_weight=None, ce_estimator={'function': 'power2_inverse_power2_decreasing', 'options': {'max_csm': 8.0}})
Bases: `WeightedNbSetChemenvStrategy`

MultiWeightsChemenvStrategy.

Constructor for the MultiWeightsChemenvStrategy.
:param structure_environments: StructureEnvironments object containing all the information on the

> coordination of the sites in a structure.


#### DEFAULT_CE_ESTIMATOR(_ = {'function': 'power2_inverse_power2_decreasing', 'options': {'max_csm': 8.0}_ )

#### STRATEGY_DESCRIPTION(_: str | Non_ _ = '    Multi Weights ChemenvStrategy_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()

* **Returns**

    Bson-serializable dict representation of the MultiWeightsChemenvStrategy object.



#### _classmethod_ from_dict(d)
Reconstructs the MultiWeightsChemenvStrategy object from a dict representation of the
MultipleAbundanceChemenvStrategy object created using the as_dict method.
:param d: dict representation of the MultiWeightsChemenvStrategy object


* **Returns**

    MultiWeightsChemenvStrategy object.



#### _classmethod_ stats_article_weights_parameters()
Initialize strategy used in the statistics article.


#### _property_ uniquely_determines_coordination_environments()
Whether this strategy uniquely determines coordination environments.


### _class_ NbSetWeight()
Bases: `MSONable`

Abstract object for neighbors sets weights estimations.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract_ as_dict()
A JSON-serializable dict representation of this neighbors set weight.


#### _abstract_ weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ NormalizedAngleDistanceNbSetWeight(average_type, aa, bb)
Bases: `NbSetWeight`

Weight of neighbors set based on the normalized angle/distance.

Initialize NormalizedAngleDistanceNbSetWeight.


* **Parameters**


    * **average_type** – Average function.


    * **aa** – Exponent for the angle values.


    * **bb** – Exponent for the distance values.



#### SHORT_NAME(_ = 'NormAngleDistWeight_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ ang(nb_set)
Angle weight.


* **Parameters**

    **nb_set** – Neighbors set.



* **Returns**

    List of angle weights.



#### _static_ anginvdist(nb_set)
Angle/distance weight.


* **Parameters**

    **nb_set** – Neighbors set.



* **Returns**

    List of angle/distance weights.



#### anginvndist(nb_set)
Angle/power distance weight.


* **Parameters**

    **nb_set** – Neighbors set.



* **Returns**

    List of angle/power distance weights.



#### angn(nb_set)
Power angle weight.


* **Parameters**

    **nb_set** – Neighbors set.



* **Returns**

    List of power angle weights.



#### angninvdist(nb_set)
Power angle/distance weight.


* **Parameters**

    **nb_set** – Neighbors set.



* **Returns**

    List of power angle/distance weights.



#### angninvndist(nb_set)
Power angle/power distance weight.


* **Parameters**

    **nb_set** – Neighbors set.



* **Returns**

    List of power angle/power distance weights.



#### as_dict()
MSONable dict.


#### _static_ aweight(fda_list)
Standard mean of the weights.


* **Parameters**

    **fda_list** – List of estimator weights for each neighbor.



* **Returns**

    Standard mean of the weights.



#### _classmethod_ from_dict(dct)
Initialize from dict.


* **Parameters**

    **dct** – Dict representation of NormalizedAngleDistanceNbSetWeight.



* **Returns**

    NormalizedAngleDistanceNbSetWeight.



#### _static_ gweight(fda_list)
Geometric mean of the weights.


* **Parameters**

    **fda_list** – List of estimator weights for each neighbor.



* **Returns**

    Geometric mean of the weights.



#### _static_ invdist(nb_set)
Inverse distance weight.


* **Parameters**

    **nb_set** – Neighbors set.



* **Returns**

    List of inverse distances.



#### invndist(nb_set)
Inverse power distance weight.


* **Parameters**

    **nb_set** – Neighbors set.



* **Returns**

    List of inverse power distances.



#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ SelfCSMNbSetWeight(effective_csm_estimator={'function': 'power2_inverse_decreasing', 'options': {'max_csm': 8.0}}, weight_estimator={'function': 'power2_decreasing_exp', 'options': {'alpha': 1, 'max_csm': 8.0}}, symmetry_measure_type='csm_wcs_ctwcc')
Bases: `NbSetWeight`

Weight of neighbors set based on the Self CSM.

Initialize SelfCSMNbSetWeight.


* **Parameters**


    * **effective_csm_estimator** – Ratio function used for the effective CSM (comparison between neighbors sets).


    * **weight_estimator** – Weight estimator within a given neighbors set.


    * **symmetry_measure_type** – Type of symmetry measure to be used.



#### DEFAULT_EFFECTIVE_CSM_ESTIMATOR(_ = {'function': 'power2_inverse_decreasing', 'options': {'max_csm': 8.0}_ )

#### DEFAULT_SYMMETRY_MEASURE_TYPE(_ = 'csm_wcs_ctwcc_ )

#### DEFAULT_WEIGHT_ESTIMATOR(_ = {'function': 'power2_decreasing_exp', 'options': {'alpha': 1, 'max_csm': 8.0}_ )

#### SHORT_NAME(_ = 'SelfCSMWeight_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dct)
Initialize from dict.


* **Parameters**

    **dct** – Dict representation of SelfCSMNbSetWeight.



* **Returns**

    SelfCSMNbSetWeight.



#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ SimpleAbundanceChemenvStrategy(structure_environments=None, additional_condition=1, symmetry_measure_type='csm_wcs_ctwcc')
Bases: `AbstractChemenvStrategy`

Simple ChemenvStrategy using the neighbors that are the most “abundant” in the grid of angle and distance
parameters for the definition of neighbors in the Voronoi approach.
The coordination environment is then given as the one with the lowest continuous symmetry measure.

Constructor for the SimpleAbundanceChemenvStrategy.
:param structure_environments: StructureEnvironments object containing all the information on the

> coordination of the sites in a structure.


#### DEFAULT_ADDITIONAL_CONDITION(_ = _ )

#### DEFAULT_MAX_DIST(_ = 2._ )

#### STRATEGY_DESCRIPTION(_: str | Non_ _ = '    Simple Abundance ChemenvStrategy using the most "abundant" neighbors map \\n    for the definition of neighbors in the Voronoi approach. \\n    The coordination environment is then given as the one with the \\n    lowest continuous symmetry measure._ )

#### STRATEGY_OPTIONS(_: ClassVar[dict[str, dict]_ _ = {'additional_condition': {'default': 1, 'internal': '_additional_condition', 'type': <class 'pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt'>}, 'surface_calculation_type': {}_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _get_map(isite)

#### _get_maps_surfaces(isite, surface_calculation_type=None)

#### as_dict()
Bson-serializable dict representation of the SimpleAbundanceChemenvStrategy object.


* **Returns**

    Bson-serializable dict representation of the SimpleAbundanceChemenvStrategy object.



#### _classmethod_ from_dict(d)
Reconstructs the SimpleAbundanceChemenvStrategy object from a dict representation of the
SimpleAbundanceChemenvStrategy object created using the as_dict method.
:param d: dict representation of the SimpleAbundanceChemenvStrategy object


* **Returns**

    StructureEnvironments object.



#### get_site_coordination_environment(site, isite=None, dequivsite=None, dthissite=None, mysym=None, return_map=False)
Get the coordination environment of a given site.


* **Parameters**


    * **site** – Site for which coordination environment is needed.


    * **isite** – Index of the site.


    * **dequivsite** – Translation of the equivalent site.


    * **dthissite** – Translation of this site.


    * **mysym** – Symmetry to be applied.


    * **return_map** – Whether to return cn_map (identifies the NeighborsSet used).



* **Returns**

    Coordination environment of site.



#### get_site_coordination_environments(site, isite=None, dequivsite=None, dthissite=None, mysym=None, return_maps=False)
Get the coordination environments of a given site.


* **Parameters**


    * **site** – Site for which coordination environment is needed.


    * **isite** – Index of the site.


    * **dequivsite** – Translation of the equivalent site.


    * **dthissite** – Translation of this site.


    * **mysym** – Symmetry to be applied.


    * **return_maps** – Whether to return cn_maps (identifies all the NeighborsSet used).



* **Returns**

    List of coordination environment.



#### get_site_neighbors(site)
Get the neighbors of a given site with this strategy.


* **Parameters**

    **site** – Periodic site.



* **Returns**

    List of neighbors of site.



#### _property_ uniquely_determines_coordination_environments()
Whether this strategy uniquely determines coordination environments.


### _class_ SimplestChemenvStrategy(structure_environments=None, distance_cutoff=1.4, angle_cutoff=0.3, additional_condition=1, continuous_symmetry_measure_cutoff=10, symmetry_measure_type='csm_wcs_ctwcc')
Bases: `AbstractChemenvStrategy`

Simplest ChemenvStrategy using fixed angle and distance parameters for the definition of neighbors in the
Voronoi approach. The coordination environment is then given as the one with the lowest continuous symmetry measure.

Constructor for this SimplestChemenvStrategy.
:param distance_cutoff: Distance cutoff used
:param angle_cutoff: Angle cutoff used.


#### DEFAULT_ADDITIONAL_CONDITION(_ = _ )

#### DEFAULT_ANGLE_CUTOFF(_ = 0._ )

#### DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF(_ = 1_ )

#### DEFAULT_DISTANCE_CUTOFF(_ = 1._ )

#### STRATEGY_DESCRIPTION(_: str | Non_ _ = '    Simplest ChemenvStrategy using fixed angle and distance parameters \\n    for the definition of neighbors in the Voronoi approach. \\n    The coordination environment is then given as the one with the \\n    lowest continuous symmetry measure._ )

#### STRATEGY_OPTIONS(_: ClassVar[dict[str, dict]_ _ = {'additional_condition': {'default': 1, 'internal': '_additional_condition', 'type': <class 'pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt'>}, 'angle_cutoff': {'default': 0.3, 'internal': '_angle_cutoff', 'type': <class 'pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat'>}, 'continuous_symmetry_measure_cutoff': {'default': 10, 'internal': '_continuous_symmetry_measure_cutoff', 'type': <class 'pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat'>}, 'distance_cutoff': {'default': 1.4, 'internal': '_distance_cutoff', 'type': <class 'pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat'>}_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### add_strategy_visualization_to_subplot(subplot, visualization_options=None, plot_type=None)
Add a visual of the strategy on a distance-angle plot.


* **Parameters**


    * **subplot** – Axes object onto the visual should be added.


    * **visualization_options** – Options for the visual.


    * **plot_type** – Type of distance-angle plot.



#### _property_ additional_condition()
Additional condition for this strategy.


#### _property_ angle_cutoff()
Angle cutoff used.


#### as_dict()
Bson-serializable dict representation of the SimplestChemenvStrategy object.


* **Returns**

    Bson-serializable dict representation of the SimplestChemenvStrategy object.



#### _property_ continuous_symmetry_measure_cutoff()
CSM cutoff used.


#### _property_ distance_cutoff()
Distance cutoff used.


#### _classmethod_ from_dict(d)
Reconstructs the SimplestChemenvStrategy object from a dict representation of the SimplestChemenvStrategy object
created using the as_dict method.
:param d: dict representation of the SimplestChemenvStrategy object


* **Returns**

    StructureEnvironments object.



#### get_site_coordination_environment(site, isite=None, dequivsite=None, dthissite=None, mysym=None, return_map=False)
Get the coordination environment of a given site.


* **Parameters**


    * **site** – Site for which coordination environment is needed.


    * **isite** – Index of the site.


    * **dequivsite** – Translation of the equivalent site.


    * **dthissite** – Translation of this site.


    * **mysym** – Symmetry to be applied.


    * **return_map** – Whether to return cn_map (identifies the NeighborsSet used).



* **Returns**

    Coordination environment of site.



#### get_site_coordination_environments(site, isite=None, dequivsite=None, dthissite=None, mysym=None, return_maps=False)
Get the coordination environments of a given site.


* **Parameters**


    * **site** – Site for which coordination environment is needed.


    * **isite** – Index of the site.


    * **dequivsite** – Translation of the equivalent site.


    * **dthissite** – Translation of this site.


    * **mysym** – Symmetry to be applied.


    * **return_maps** – Whether to return cn_maps (identifies all the NeighborsSet used).



* **Returns**

    List of coordination environment.



#### get_site_coordination_environments_fractions(site, isite=None, dequivsite=None, dthissite=None, mysym=None, ordered=True, min_fraction=0, return_maps=True, return_strategy_dict_info=False)
Get the coordination environments of a given site and additional information.


* **Parameters**


    * **site** – Site for which coordination environment is needed.


    * **isite** – Index of the site.


    * **dequivsite** – Translation of the equivalent site.


    * **dthissite** – Translation of this site.


    * **mysym** – Symmetry to be applied.


    * **ordered** – Whether to order the list by fractions.


    * **min_fraction** – Minimum fraction to include in the list


    * **return_maps** – Whether to return cn_maps (identifies all the NeighborsSet used).


    * **return_strategy_dict_info** – Whether to add the info about the strategy used.



* **Returns**

    List of Dict with coordination environment, fraction and additional info.



#### get_site_neighbors(site, isite=None, dequivsite=None, dthissite=None, mysym=None)
Get the neighbors of a given site.


* **Parameters**


    * **site** – Site for which neighbors are needed.


    * **isite** – Index of the site.


    * **dequivsite** – Translation of the equivalent site.


    * **dthissite** – Translation of this site.


    * **mysym** – Symmetry to be applied.



* **Returns**

    List of coordinated neighbors of site.



#### _property_ uniquely_determines_coordination_environments()
Whether this strategy uniquely determines coordination environments.


### _class_ StrategyOption()
Bases: `MSONable`

Abstract class for the options of the chemenv strategies.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### allowed_values(_: str | Non_ _ = Non_ )

#### _abstract_ as_dict()
A JSON-serializable dict representation of this strategy option.


### _class_ TargettedPenaltiedAbundanceChemenvStrategy(structure_environments=None, truncate_dist_ang=True, additional_condition=1, max_nabundant=5, target_environments=('O:6',), target_penalty_type='max_csm', max_csm=5.0, symmetry_measure_type='csm_wcs_ctwcc')
Bases: `SimpleAbundanceChemenvStrategy`

Simple ChemenvStrategy using the neighbors that are the most “abundant” in the grid of angle and distance
parameters for the definition of neighbors in the Voronoi approach, with a bias for a given list of target
environments. This can be useful in the case of, e.g. connectivity search of some given environment.
The coordination environment is then given as the one with the lowest continuous symmetry measure.

Initializes strategy.

Not yet implemented.
:param structure_environments:
:param truncate_dist_ang:
:param additional_condition:
:param max_nabundant:
:param target_environments:
:param target_penalty_type:
:param max_csm:
:param symmetry_measure_type:


#### DEFAULT_TARGET_ENVIRONMENTS(_ = ('O:6',_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### _get_map(isite)

#### as_dict()
Bson-serializable dict representation of the TargettedPenaltiedAbundanceChemenvStrategy object.


* **Returns**

    Bson-serializable dict representation of the TargettedPenaltiedAbundanceChemenvStrategy object.



#### _classmethod_ from_dict(d)
Reconstructs the TargettedPenaltiedAbundanceChemenvStrategy object from a dict representation of the
TargettedPenaltiedAbundanceChemenvStrategy object created using the as_dict method.
:param d: dict representation of the TargettedPenaltiedAbundanceChemenvStrategy object


* **Returns**

    TargettedPenaltiedAbundanceChemenvStrategy object.



#### get_site_coordination_environment(site, isite=None, dequivsite=None, dthissite=None, mysym=None, return_map=False)
Get the coordination environment of a given site.


* **Parameters**


    * **site** – Site for which coordination environment is needed.


    * **isite** – Index of the site.


    * **dequivsite** – Translation of the equivalent site.


    * **dthissite** – Translation of this site.


    * **mysym** – Symmetry to be applied.


    * **return_map** – Whether to return cn_map (identifies the NeighborsSet used).



* **Returns**

    Coordination environment of site.



#### _property_ uniquely_determines_coordination_environments()
Whether this strategy uniquely determines coordination environments.


### _class_ WeightedNbSetChemenvStrategy(structure_environments=None, additional_condition=1, symmetry_measure_type='csm_wcs_ctwcc', nb_set_weights=None, ce_estimator={'function': 'power2_inverse_power2_decreasing', 'options': {'max_csm': 8.0}})
Bases: `AbstractChemenvStrategy`

WeightedNbSetChemenvStrategy.

Constructor for the WeightedNbSetChemenvStrategy.
:param structure_environments: StructureEnvironments object containing all the information on the

> coordination of the sites in a structure.


#### DEFAULT_CE_ESTIMATOR(_ = {'function': 'power2_inverse_power2_decreasing', 'options': {'max_csm': 8.0}_ )

#### STRATEGY_DESCRIPTION(_: str | Non_ _ = '    WeightedNbSetChemenvStrategy_ )

#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
Bson-serializable dict representation of the WeightedNbSetChemenvStrategy object.


* **Returns**

    Bson-serializable dict representation of the WeightedNbSetChemenvStrategy object.



#### _classmethod_ from_dict(d)
Reconstructs the WeightedNbSetChemenvStrategy object from a dict representation of the
WeightedNbSetChemenvStrategy object created using the as_dict method.
:param d: dict representation of the WeightedNbSetChemenvStrategy object


* **Returns**

    WeightedNbSetChemenvStrategy object.



#### get_site_coordination_environment(site)
Get the coordination environment of a given site.

Not implemented for this strategy


#### get_site_coordination_environments(site, isite=None, dequivsite=None, dthissite=None, mysym=None, return_maps=False)
Get the coordination environments of a given site.


* **Parameters**


    * **site** – Site for which coordination environment is needed.


    * **isite** – Index of the site.


    * **dequivsite** – Translation of the equivalent site.


    * **dthissite** – Translation of this site.


    * **mysym** – Symmetry to be applied.


    * **return_maps** – Whether to return cn_maps (identifies all the NeighborsSet used).



* **Returns**

    List of coordination environment.



#### get_site_coordination_environments_fractions(site, isite=None, dequivsite=None, dthissite=None, mysym=None, ordered=True, min_fraction=0, return_maps=True, return_strategy_dict_info=False, return_all=False)
Get the coordination environments of a given site and additional information.


* **Parameters**


    * **site** – Site for which coordination environment is needed.


    * **isite** – Index of the site.


    * **dequivsite** – Translation of the equivalent site.


    * **dthissite** – Translation of this site.


    * **mysym** – Symmetry to be applied.


    * **ordered** – Whether to order the list by fractions.


    * **min_fraction** – Minimum fraction to include in the list


    * **return_maps** – Whether to return cn_maps (identifies all the NeighborsSet used).


    * **return_strategy_dict_info** – Whether to add the info about the strategy used.



* **Returns**

    List of Dict with coordination environment, fraction and additional info.



#### get_site_neighbors(site)
Get the neighbors of a given site.

Not implemented for this strategy.


#### _property_ uniquely_determines_coordination_environments()
Whether this strategy uniquely determines coordination environments.


### get_effective_csm(nb_set, cn_map, structure_environments, additional_info, symmetry_measure_type, max_effective_csm, effective_csm_estimator_ratio_function)
Get the effective continuous symmetry measure of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **cn_map** – Mapping index of this neighbors set.


    * **structure_environments** – Structure environments.


    * **additional_info** – Additional information for the neighbors set.


    * **symmetry_measure_type** – Type of symmetry measure to be used in the effective CSM.


    * **max_effective_csm** – Max CSM to use for the effective CSM calculation.


    * **effective_csm_estimator_ratio_function** – Ratio function to use to compute effective CSM.



* **Returns**

    Effective CSM of a given Neighbors set.



### set_info(additional_info, field, isite, cn_map, value)
Set additional information for the weights.


* **Parameters**


    * **additional_info** – Additional information.


    * **field** – Type of additional information.


    * **isite** – Index of site to add info.


    * **cn_map** – Mapping index of the neighbors set.


    * **value** – Value of this additional information.



* **Returns**

    None


## pymatgen.analysis.chemenv.coordination_environments.coordination_geometries module

This module contains the class describing the coordination geometries that can exist in a given structure. These
“model” coordination geometries are described in the following articles :

>
> * Pure Appl. Chem., Vol. 79, No. 10, pp. 1779–1799, 2007.


> * Acta Cryst. A, Vol. 46, No. 1, pp. 1–11, 1990.

The module also contains descriptors of part of these geometries (plane of separation, …) that are used in the
identification algorithms.


### _class_ AbstractChemenvAlgorithm(algorithm_type)
Bases: `MSONable`

Base class used to define a Chemenv algorithm used to identify the correct permutation for the computation
of the Continuous Symmetry Measure.

Base constructor for ChemenvAlgorithm.


* **Parameters**

    **algorithm_type** (*str*) – Type of algorithm.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ algorithm_type()
Return the type of algorithm.


* **Returns**

    Type of the algorithm



* **Return type**

    str



#### _abstract_ as_dict()
A JSON-serializable dict representation of the algorithm.


### _class_ AllCoordinationGeometries(permutations_safe_override=False, only_symbols=None)
Bases: `dict`

Class used to store all the reference “coordination geometries” (list with instances of the CoordinationGeometry
classes).

Initializes the list of Coordination Geometries.


* **Parameters**


    * **permutations_safe_override** – Whether to use safe permutations.


    * **only_symbols** – Whether to restrict the list of environments to be identified.



#### get_geometries(coordination=None, returned='cg')
Returns a list of coordination geometries with the given coordination number.


* **Parameters**


    * **coordination** – The coordination number of which the list of coordination geometries are returned.


    * **returned** – Type of objects in the list.



#### get_geometry_from_IUCr_symbol(IUCr_symbol)
Returns the coordination geometry of the given IUCr symbol.


* **Parameters**

    **IUCr_symbol** – The IUCr symbol of the coordination geometry.



#### get_geometry_from_IUPAC_symbol(IUPAC_symbol)
Returns the coordination geometry of the given IUPAC symbol.


* **Parameters**

    **IUPAC_symbol** – The IUPAC symbol of the coordination geometry.



#### get_geometry_from_mp_symbol(mp_symbol)
Returns the coordination geometry of the given mp_symbol.


* **Parameters**

    **mp_symbol** – The mp_symbol of the coordination geometry.



#### get_geometry_from_name(name)
Returns the coordination geometry of the given name.


* **Parameters**

    **name** – The name of the coordination geometry.



#### get_implemented_geometries(coordination=None, returned='cg', include_deactivated=False)
Returns a list of the implemented coordination geometries with the given coordination number.


* **Parameters**


    * **coordination** – The coordination number of which the list of implemented coordination geometries
    are returned.


    * **returned** – Type of objects in the list.


    * **include_deactivated** – Whether to include CoordinationGeometry that are deactivated.



#### get_not_implemented_geometries(coordination=None, returned='mp_symbol')
Returns a list of the implemented coordination geometries with the given coordination number.


* **Parameters**


    * **coordination** – The coordination number of which the list of implemented coordination geometries
    are returned.


    * **returned** – Type of objects in the list.



#### get_symbol_cn_mapping(coordination=None)
Return a dictionary mapping the symbol of a CoordinationGeometry to its coordination.


* **Parameters**

    **coordination** – Whether to restrict the dictionary to a given coordination.



* **Returns**

    map of symbol of a CoordinationGeometry to its coordination.



* **Return type**

    dict



#### get_symbol_name_mapping(coordination=None)
Return a dictionary mapping the symbol of a CoordinationGeometry to its name.


* **Parameters**

    **coordination** – Whether to restrict the dictionary to a given coordination.



* **Returns**

    map symbol of a CoordinationGeometry to its name.



* **Return type**

    dict



#### is_a_valid_coordination_geometry(mp_symbol=None, IUPAC_symbol=None, IUCr_symbol=None, name=None, cn=None)
Checks whether a given coordination geometry is valid (exists) and whether the parameters are coherent with
each other.


* **Parameters**


    * **mp_symbol** – The mp_symbol of the coordination geometry.


    * **IUPAC_symbol** – The IUPAC_symbol of the coordination geometry.


    * **IUCr_symbol** – The IUCr_symbol of the coordination geometry.


    * **name** – The name of the coordination geometry.


    * **cn** – The coordination of the coordination geometry.



#### pretty_print(type='implemented_geometries', maxcn=8, additional_info=None)
Return a string with a list of the Coordination Geometries.


* **Parameters**


    * **type** – Type of string to be returned (all_geometries, all_geometries_latex_images, all_geometries_latex,
    implemented_geometries).


    * **maxcn** – Maximum coordination.


    * **additional_info** – Whether to add some additional info for each coordination geometry.



* **Returns**

    description of the list of coordination geometries.



* **Return type**

    str



### _class_ CoordinationGeometry(mp_symbol, name, alternative_names=None, IUPAC_symbol=None, IUCr_symbol=None, coordination=None, central_site=None, points=None, solid_angles=None, permutations_safe_override=False, deactivate=False, faces=None, edges=None, algorithms=None, equivalent_indices=None, neighbors_sets_hints=None)
Bases: `object`

Class used to store the ideal representation of a chemical environment or “coordination geometry”.

Initializes one “coordination geometry” according to [Pure Appl. Chem., Vol. 79, No. 10, pp. 1779–1799, 2007]
and [Acta Cryst. A, Vol. 46, No. 1, pp. 1–11, 1990].


* **Parameters**


    * **mp_symbol** – Symbol used internally for the coordination geometry.


    * **name** – Name of the coordination geometry.


    * **alternative_names** – Alternative names for this coordination geometry.


    * **IUPAC_symbol** – The IUPAC symbol of this coordination geometry.


    * **IUCr_symbol** – The IUCr symbol of this coordination geometry.


    * **coordination** – The coordination number of this coordination geometry (number of neighboring atoms).


    * **central_site** – The coordinates of the central site of this coordination geometry.


    * **points** – The list of the coordinates of all the points of this coordination geometry.


    * **solid_angles** – The list of solid angles for each neighbor in this coordination geometry.


    * **permutations_safe_override** – Computes all the permutations if set to True (overrides the plane separation
    algorithms or any other algorithm, for testing purposes)


    * **deactivate** – Whether to deactivate this coordination geometry


    * **faces** – List of the faces with their vertices given in a clockwise or anticlockwise order, for drawing
    purposes.


    * **edges** – List of edges, for drawing purposes.


    * **algorithms** – Algorithms used to identify this coordination geometry.


    * **equivalent_indices** – The equivalent sets of indices in this coordination geometry (can be used to skip
    equivalent permutations that have already been performed).


    * **neighbors_sets_hints** – Neighbors sets hints for this coordination geometry.



#### CSM_SKIP_SEPARATION_PLANE_ALGO(_ = 10._ )

#### _property_ IUCr_symbol()
Returns the IUCr symbol of this coordination geometry.


#### _property_ IUCr_symbol_str()
Returns a string representation of the IUCr symbol of this coordination geometry.


#### _property_ IUPAC_symbol()
Returns the IUPAC symbol of this coordination geometry.


#### _property_ IUPAC_symbol_str()
Returns a string representation of the IUPAC symbol of this coordination geometry.


#### _class_ NeighborsSetsHints(hints_type, options)
Bases: `object`

Class used to describe neighbors sets hints.

This allows to possibly get a lower coordination from a capped-like model polyhedron.

Constructor for this NeighborsSetsHints.


* **Parameters**


    * **hints_type** – type of hint (single, double or triple cap)


    * **options** – options for the “hinting”, e.g. the maximum csm value beyond which no additional
    neighbors set could be found from a “cap hint”.



#### ALLOWED_HINTS_TYPES(_ = ('single_cap', 'double_cap', 'triple_cap'_ )

#### as_dict()
A JSON-serializable dict representation of this NeighborsSetsHints.


#### double_cap_hints(hints_info)
Return hints for an additional neighbors set, i.e. the voronoi indices that
constitute this new neighbors set, in case of a “Double cap” hint.


* **Parameters**

    **hints_info** – Info needed to build new “hinted” neighbors set.



* **Returns**

    Voronoi indices of the new “hinted” neighbors set.



* **Return type**

    list[int]



#### _classmethod_ from_dict(dct)
Reconstructs the NeighborsSetsHints from its JSON-serializable dict representation.


#### hints(hints_info)
Return hints for an additional neighbors set, i.e. the voronoi indices that
constitute this new neighbors set.


* **Parameters**

    **hints_info** – Info needed to build new “hinted” neighbors set.



* **Returns**

    Voronoi indices of the new “hinted” neighbors set.



* **Return type**

    list[int]



#### single_cap_hints(hints_info)
Return hints for an additional neighbors set, i.e. the voronoi indices that
constitute this new neighbors set, in case of a “Single cap” hint.


* **Parameters**

    **hints_info** – Info needed to build new “hinted” neighbors set.



* **Returns**

    Voronoi indices of the new “hinted” neighbors set.



* **Return type**

    list[int]



#### triple_cap_hints(hints_info)
Return hints for an additional neighbors set, i.e. the voronoi indices that
constitute this new neighbors set, in case of a “Triple cap” hint.


* **Parameters**

    **hints_info** – Info needed to build new “hinted” neighbors set.



* **Returns**

    Voronoi indices of the new “hinted” neighbors set.



* **Return type**

    list[int]



#### _property_ algorithms()
Returns the list of algorithms that are used to identify this coordination geometry.


#### as_dict()
A JSON-serializable dict representation of this CoordinationGeometry.


#### _property_ ce_symbol()
Returns the symbol of this coordination geometry.


#### _property_ coordination_number()
Returns the coordination number of this coordination geometry.


#### _property_ distfactor_max()
The maximum distfactor for the perfect CoordinationGeometry (usually 1.0 for symmetric polyhedrons).


#### edges(sites, permutation=None, input='sites')
Returns the list of edges of this coordination geometry. Each edge is given as a
list of its end vertices coordinates.


#### faces(sites, permutation=None)
Returns the list of faces of this coordination geometry. Each face is given as a
list of its vertices coordinates.


#### _classmethod_ from_dict(dct)
Reconstructs the CoordinationGeometry from its JSON-serializable dict representation.


* **Parameters**

    **dct** – a JSON-serializable dict representation of a CoordinationGeometry.



* **Returns**

    CoordinationGeometry



#### get_central_site()
Returns the central site of this coordination geometry.


#### get_coordination_number()
Returns the coordination number of this coordination geometry.


#### get_name()
Returns the name of this coordination geometry.


#### get_pmeshes(sites, permutation=None)
Returns the pmesh strings used for jmol to show this geometry.


#### is_implemented()
Returns True if this coordination geometry is implemented.


#### _property_ mp_symbol()
Returns the MP symbol of this coordination geometry.


#### _property_ number_of_permutations()
Returns the number of permutations of this coordination geometry.


#### _property_ pauling_stability_ratio()
Returns the theoretical Pauling stability ratio (rC/rA) for this environment.


#### ref_permutation(permutation)
Returns the reference permutation for a set of equivalent permutations.

Can be useful to skip permutations that have already been performed.


* **Parameters**

    **permutation** – Current permutation



* **Returns**

    Reference permutation of the perfect CoordinationGeometry.



* **Return type**

    Permutation



#### solid_angles(permutation=None)
Returns the list of “perfect” solid angles Each edge is given as a
list of its end vertices coordinates.


### _class_ ExplicitPermutationsAlgorithm(permutations)
Bases: `AbstractChemenvAlgorithm`

Class representing the algorithm doing the explicit permutations for the calculation of
the Continuous Symmetry Measure.

Initializes a separation plane for a given perfect coordination geometry.


* **Parameters**

    **permutations** – Permutations used for this algorithm.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ as_dict()
Returns:
dict: JSON-serializable representation of this ExplicitPermutationsAlgorithm


#### _classmethod_ from_dict(dct)
Reconstruct ExplicitPermutationsAlgorithm from its JSON-serializable dict representation.


#### _property_ permutations()
Return the permutations to be performed for this algorithm.


* **Returns**

    Permutations to be performed.



* **Return type**

    list



### _class_ SeparationPlane(plane_points, mirror_plane=False, ordered_plane=False, point_groups=None, ordered_point_groups=None, explicit_permutations=None, minimum_number_of_points=None, explicit_optimized_permutations=None, multiplicity=None, other_plane_points=None)
Bases: `AbstractChemenvAlgorithm`

Class representing the algorithm using separation planes for the calculation of
the Continuous Symmetry Measure.

Initializes a separation plane for a given perfect coordination geometry.


* **Parameters**


    * **plane_points** – Indices of the points that are in the plane in the perfect structure (and should be
    found in the defective one as well).


    * **mirror_plane** – True if the separation plane is a mirror plane, in which case there is a correspondence
    of the points in each point_group (can reduce the number of permutations).


    * **ordered_plane** – True if the order of the points in the plane can be taken into account to reduce the
    number of permutations.


    * **point_groups** – Indices of the points in the two groups of points separated by the plane.


    * **ordered_point_groups** – Whether the order of the points in each group of points can be taken into account to
    reduce the number of permutations.


    * **explicit_permutations** – Explicit permutations to be performed in this separation plane algorithm.


    * **minimum_number_of_points** – Minimum number of points needed to initialize a separation plane
    for this algorithm.


    * **explicit_optimized_permutations** – Optimized set of explicit permutations to be performed in this
    separation plane algorithm.


    * **multiplicity** – Number of such planes in the model geometry.


    * **other_plane_points** – Indices of the points that are in the plane in the perfect structure for the other
    planes. The multiplicity should be equal to the length of this list + 1 (“main” separation plane +
    the other ones).



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _property_ argsorted_ref_separation_perm()
“Arg sorted” ordered indices of the separation plane.

This is used in the identification of the final permutation to be used.


* **Returns**

    “arg sorted” ordered indices of the separation plane.



* **Return type**

    list[int]



#### _property_ as_dict()
Return the JSON-serializable dict representation of this SeparationPlane algorithm.


* **Returns**

    JSON-serializable representation of this SeparationPlane algorithm.



* **Return type**

    dict



#### _classmethod_ from_dict(dct)
Reconstructs the SeparationPlane algorithm from its JSON-serializable dict representation.


* **Parameters**

    **dct** – a JSON-serializable dict representation of an SeparationPlane algorithm.



* **Returns**

    algorithm object



* **Return type**

    SeparationPlane



#### _property_ permutations()
Permutations used for this separation plane algorithm.


* **Returns**

    to be performed.



* **Return type**

    list[Permutations]



#### _property_ ref_separation_perm()
Ordered indices of the separation plane.

### Examples

For a separation plane of type 2|4|3, with plane_points indices [0, 3, 5, 8] and
point_groups indices [1, 4] and [2, 7, 6], the list of ordered indices is :
[0, 3, 5, 8, 1, 4, 2, 7, 6].


* **Returns**

    of ordered indices of this separation plane.



* **Return type**

    list[int]



#### safe_separation_permutations(ordered_plane=False, ordered_point_groups=None, add_opposite=False)
Simple and safe permutations for this separation plane.

This is not meant to be used in production. Default configuration for ChemEnv does not use this method.


* **Parameters**


    * **ordered_plane** – Whether the order of the points in the plane can be used to reduce the
    number of permutations.


    * **ordered_point_groups** – Whether the order of the points in each point group can be used to reduce the
    number of permutations.


    * **add_opposite** – Whether to add the permutations from the second group before the first group as well.


Returns

    list[int]: safe permutations.

## pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder module

This module contains the main object used to identify the coordination environments in a given structure.
If you use this module, please cite:
David Waroquiers, Xavier Gonze, Gian-Marco Rignanese, Cathrin Welker-Nieuwoudt, Frank Rosowski,
Michael Goebel, Stephan Schenk, Peter Degelmann, Rute Andre, Robert Glaum, and Geoffroy Hautier,
“Statistical analysis of coordination environments in oxides”,
Chem. Mater., 2017, 29 (19), pp 8346-8360,
DOI: 10.1021/acs.chemmater.7b02766
D. Waroquiers, J. George, M. Horton, S. Schenk, K. A. Persson, G.-M. Rignanese, X. Gonze, G. Hautier
“ChemEnv: a fast and robust coordination environment identification tool”,
Acta Cryst. B 2020, 76, pp 683-695,
DOI: 10.1107/S2052520620007994.


### _class_ AbstractGeometry(central_site=None, bare_coords=None, centering_type='standard', include_central_site_in_centroid=False, optimization=None)
Bases: `object`

Class used to describe a geometry (perfect or distorted).

Constructor for the abstract geometry
:param central_site: Coordinates of the central site
:param bare_coords: Coordinates of the neighbors of the central site
:param centering_type: How to center the abstract geometry
:param include_central_site_in_centroid: When the centering is on the centroid, the central site is included

> if this parameter is set to True.


* **Raise**

    ValueError if the parameters are not consistent.



#### _property_ cn()
Coordination number


#### _property_ coordination_number()
Coordination number


#### _classmethod_ from_cg(cg, centering_type='standard', include_central_site_in_centroid=False)

* **Parameters**


    * **cg** –


    * **centering_type** –


    * **include_central_site_in_centroid** –



#### points_wcs_csc(permutation=None)

* **Parameters**

    **permutation** –



#### points_wcs_ctwcc(permutation=None)

* **Parameters**

    **permutation** –



#### points_wcs_ctwocc(permutation=None)

* **Parameters**

    **permutation** –



#### points_wocs_csc(permutation=None)

* **Parameters**

    **permutation** –



#### points_wocs_ctwcc(permutation=None)

* **Parameters**

    **permutation** –



#### points_wocs_ctwocc(permutation=None)

* **Parameters**

    **permutation** –



### _class_ LocalGeometryFinder(permutations_safe_override: bool = False, plane_ordering_override: bool = True, plane_safe_permutations: bool = False, only_symbols=None)
Bases: `object`

Main class used to find the local environments in a structure.


* **Parameters**


    * **permutations_safe_override** – If set to True, all permutations are tested (very time-consuming for large


    * **numbers!****)** (*coordination*) –


    * **plane_ordering_override** – If set to False, the ordering of the points in the plane is disabled


    * **plane_safe_permutations** – Whether to use safe permutations.


    * **only_symbols** – Whether to restrict the list of environments to be identified.



#### BVA_DISTANCE_SCALE_FACTORS(_ = {'GGA_relaxed': 1.015, 'LDA_relaxed': 0.995, 'experimental': 1.0_ )

#### DEFAULT_BVA_DISTANCE_SCALE_FACTOR(_ = 1._ )

#### DEFAULT_SPG_ANALYZER_OPTIONS(_ = {'angle_tolerance': 5, 'symprec': 0.001_ )

#### DEFAULT_STRATEGY(_ = <pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy object_ )

#### PRESETS(_ = {'DEFAULT': {'maximum_distance_factor': 2.0, 'minimum_angle_factor': 0.05, 'optimization': 2, 'voronoi_normalized_angle_tolerance': 0.03, 'voronoi_normalized_distance_tolerance': 0.05}_ )

#### STRUCTURE_REFINEMENT_NONE(_ = 'none_ )

#### STRUCTURE_REFINEMENT_REFINED(_ = 'refined_ )

#### STRUCTURE_REFINEMENT_SYMMETRIZED(_ = 'symmetrized_ )

#### _cg_csm_separation_plane(coordination_geometry, sep_plane, local_plane, plane_separations, dist_tolerances=None, testing=False, tested_permutations=False, points_perfect=None)

#### _cg_csm_separation_plane_optim1(coordination_geometry, sepplane, local_plane, points_perfect=None, separation_indices=None)

#### _cg_csm_separation_plane_optim2(coordination_geometry, sepplane, local_plane, points_perfect=None, separation_indices=None)

#### _update_results_all_csms(result_dict, permutations, imin, geometry)

#### compute_coordination_environments(structure, indices=None, only_cations=True, strategy=<pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy object>, valences='bond-valence-analysis', initial_structure_environments=None)

* **Parameters**


    * **structure** –


    * **indices** –


    * **only_cations** –


    * **strategy** –


    * **valences** –


    * **initial_structure_environments** –



#### compute_structure_environments(excluded_atoms=None, only_atoms=None, only_cations=True, only_indices=None, maximum_distance_factor=2.0, minimum_angle_factor=0.05, max_cn=None, min_cn=None, only_symbols=None, valences='undefined', additional_conditions=None, info=None, timelimit=None, initial_structure_environments=None, get_from_hints=False, voronoi_normalized_distance_tolerance=0.05, voronoi_normalized_angle_tolerance=0.03, voronoi_distance_cutoff=None, recompute=None, optimization=2)
Computes and returns the StructureEnvironments object containing all the information about the coordination
environments in the structure
:param excluded_atoms: Atoms for which the coordination geometries does not have to be identified
:param only_atoms: If not set to None, atoms for which the coordination geometries have to be identified
:param only_cations: If set to True, will only compute environments for cations
:param only_indices: If not set to None, will only compute environments the atoms of the given indices
:param maximum_distance_factor: If not set to None, neighbors beyond

> maximum_distance_factor\*closest_neighbor_distance are not considered


* **Parameters**


    * **minimum_angle_factor** – If not set to None, neighbors for which the angle is lower than
    minimum_angle_factor\*largest_angle_neighbor are not considered


    * **max_cn** – maximum coordination number to be considered


    * **min_cn** – minimum coordination number to be considered


    * **only_symbols** – if not set to None, consider only coordination environments with the given symbols


    * **valences** – valences of the atoms


    * **additional_conditions** – additional conditions to be considered in the bonds (example : only bonds
    between cation and anion


    * **info** – additional info about the calculation


    * **timelimit** – time limit (in secs) after which the calculation of the StructureEnvironments object stops


    * **initial_structure_environments** – initial StructureEnvironments object (most probably incomplete)


    * **get_from_hints** – whether to add neighbors sets from “hints” (e.g. capped environment => test the
    neighbors without the cap)


    * **voronoi_normalized_distance_tolerance** – tolerance for the normalized distance used to distinguish
    neighbors sets


    * **voronoi_normalized_angle_tolerance** – tolerance for the normalized angle used to distinguish
    neighbors sets


    * **voronoi_distance_cutoff** – determines distance of considered neighbors. Especially important to increase it
    for molecules in a box.


    * **recompute** – whether to recompute the sites already computed (when initial_structure_environments
    is not None)


    * **optimization** – optimization algorithm



* **Returns**

    The StructureEnvironments object containing all the information about the coordination
    environments in the structure.



#### coordination_geometry_symmetry_measures(coordination_geometry, tested_permutations=False, points_perfect=None, optimization=None)
Returns the symmetry measures of a given coordination_geometry for a set of
permutations depending on the permutation setup. Depending on the parameters of
the LocalGeometryFinder and on the coordination geometry, different methods are called.


* **Parameters**

    **coordination_geometry** – Coordination geometry for which the symmetry measures are looked for



* **Raises**

    **NotImplementedError** – if the permutation_setup does not exist



* **Returns**

    the symmetry measures of a given coordination_geometry for a set of permutations



#### coordination_geometry_symmetry_measures_fallback_random(coordination_geometry, NRANDOM=10, points_perfect=None)
Returns the symmetry measures for a random set of permutations for the coordination geometry
“coordination_geometry”. Fallback implementation for the plane separation algorithms measures
of each permutation
:param coordination_geometry: The coordination geometry to be investigated
:param NRANDOM: Number of random permutations to be tested


* **Returns**

    The symmetry measures for the given coordination geometry for each permutation investigated.



#### coordination_geometry_symmetry_measures_separation_plane(coordination_geometry, separation_plane_algo, testing=False, tested_permutations=False, points_perfect=None)
Returns the symmetry measures of the given coordination geometry “coordination_geometry” using separation
facets to reduce the complexity of the system. Caller to the refined 2POINTS, 3POINTS and other …
:param coordination_geometry: The coordination geometry to be investigated


* **Returns**

    The symmetry measures for the given coordination geometry for each plane and permutation investigated.



#### coordination_geometry_symmetry_measures_separation_plane_optim(coordination_geometry, separation_plane_algo, points_perfect=None, nb_set=None, optimization=None)
Returns the symmetry measures of the given coordination geometry “coordination_geometry” using separation
facets to reduce the complexity of the system. Caller to the refined 2POINTS, 3POINTS and other …


* **Parameters**


    * **coordination_geometry** – The coordination geometry to be investigated.


    * **separation_plane_algo** – Separation Plane algorithm used.


    * **points_perfect** – Points corresponding to the perfect geometry.


    * **nb_set** – Neighbor set for this set of points. (used to store already computed separation planes)


    * **optimization** – Optimization level (1 or 2).



* **Returns**

    Continuous symmetry measures for the given coordination geometry for each plane and permutation

        investigated, corresponding permutations, corresponding algorithms,
        corresponding mappings from local to perfect environment and corresponding mappings
        from perfect to local environment.




* **Return type**

    tuple



#### coordination_geometry_symmetry_measures_sepplane_optim(coordination_geometry, points_perfect=None, nb_set=None, optimization=None)
Returns the symmetry measures of a given coordination_geometry for a set of
permutations depending on the permutation setup. Depending on the parameters of
the LocalGeometryFinder and on the coordination geometry, different methods are called.


* **Parameters**

    **coordination_geometry** – Coordination geometry for which the symmetry measures are looked for



* **Raises**

    **NotImplementedError** – if the permutation_setup does not exist



* **Returns**

    the symmetry measures of a given coordination_geometry for a set of permutations



#### coordination_geometry_symmetry_measures_standard(coordination_geometry, algo, points_perfect=None, optimization=None)
Returns the symmetry measures for a set of permutations (whose setup depends on the coordination geometry)
for the coordination geometry “coordination_geometry”. Standard implementation looking for the symmetry
measures of each permutation
:param coordination_geometry: The coordination geometry to be investigated


* **Returns**

    The symmetry measures for the given coordination geometry for each permutation investigated.



#### get_coordination_symmetry_measures(only_minimum=True, all_csms=True, optimization=None)
Returns the continuous symmetry measures of the current local geometry in a dictionary.


* **Returns**

    the continuous symmetry measures of the current local geometry in a dictionary.



#### get_coordination_symmetry_measures_optim(only_minimum=True, all_csms=True, nb_set=None, optimization=None)
Returns the continuous symmetry measures of the current local geometry in a dictionary.


* **Returns**

    the continuous symmetry measures of the current local geometry in a dictionary.



#### get_structure()
Returns the pymatgen Structure that has been setup for the identification of geometries (the initial one
might have been refined/symmetrized using the SpaceGroupAnalyzer).


* **Returns**

    The pymatgen Structure that has been setup for the identification of geometries (the initial one


might have been refined/symmetrized using the SpaceGroupAnalyzer).


#### set_structure(lattice: [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice), species, coords, coords_are_cartesian)
Sets up the pymatgen structure for which the coordination geometries have to be identified starting from the
lattice, the species and the coordinates
:param lattice: The lattice of the structure
:param species: The species on the sites
:param coords: The coordinates of the sites
:param coords_are_cartesian: If set to True, the coordinates are given in Cartesian coordinates.


#### setup_explicit_indices_local_geometry(explicit_indices)
Sets up explicit indices for the local geometry, for testing purposes
:param explicit_indices: explicit indices for the neighbors (set of numbers
from 0 to CN-1 in a given order).


#### setup_local_geometry(isite, coords, optimization=None)
Sets up the AbstractGeometry for the local geometry of site with index isite.
:param isite: Index of the site for which the local geometry has to be set up
:param coords: The coordinates of the (local) neighbors.


#### setup_ordered_indices_local_geometry(coordination)
Sets up ordered indices for the local geometry, for testing purposes
:param coordination: coordination of the local geometry.


#### setup_parameter(parameter, value)
Setup of one specific parameter to the given value. The other parameters are unchanged. See setup_parameters
method for the list of possible parameters
:param parameter: Parameter to setup/update
:param value: Value of the parameter.


#### setup_parameters(centering_type='standard', include_central_site_in_centroid=False, bva_distance_scale_factor=None, structure_refinement='refined', spg_analyzer_options=None)
Setup of the parameters for the coordination geometry finder. A reference point for the geometries has to be
chosen. This can be the centroid of the structure (including or excluding the atom for which the coordination
geometry is looked for) or the atom itself. In the ‘standard’ centering_type, the reference point is the central
atom for coordination numbers 1, 2, 3 and 4 and the centroid for coordination numbers > 4.
:param centering_type: Type of the reference point (centering) ‘standard’, ‘centroid’ or ‘central_site’
:param include_central_site_in_centroid: In case centering_type is ‘centroid’, the central site is included if

> this value is set to True.


* **Parameters**


    * **bva_distance_scale_factor** – Scaling factor for the bond valence analyzer (this might be different whether
    the structure is an experimental one, an LDA or a GGA relaxed one, or any other relaxation scheme (where
    under- or over-estimation of bond lengths is known).


    * **structure_refinement** – Refinement of the structure. Can be “none”, “refined” or “symmetrized”.


    * **spg_analyzer_options** – Options for the SpaceGroupAnalyzer (dictionary specifying “symprec”
    and “angle_tolerance”. See pymatgen’s SpaceGroupAnalyzer for more information.



#### setup_random_indices_local_geometry(coordination)
Sets up random indices for the local geometry, for testing purposes
:param coordination: coordination of the local geometry.


#### setup_random_structure(coordination)
Sets up a purely random structure with a given coordination.
:param coordination: coordination number for the random structure.


#### setup_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Sets up the structure for which the coordination geometries have to be identified. The structure is analyzed
with the space group analyzer and a refined structure is used
:param structure: A pymatgen Structure.


#### setup_test_perfect_environment(symbol, randomness=False, max_random_dist=0.1, symbol_type='mp_symbol', indices='RANDOM', random_translation='NONE', random_rotation='NONE', random_scale='NONE', points=None)

* **Parameters**


    * **symbol** –


    * **randomness** –


    * **max_random_dist** –


    * **symbol_type** –


    * **indices** –


    * **random_translation** –


    * **random_rotation** –


    * **random_scale** –


    * **points** –



#### update_nb_set_environments(se, isite, cn, inb_set, nb_set, recompute=False, optimization=None)

* **Parameters**


    * **se** –


    * **isite** –


    * **cn** –


    * **inb_set** –


    * **nb_set** –


    * **recompute** –


    * **optimization** –



### find_rotation(points_distorted, points_perfect)
This finds the rotation matrix that aligns the (distorted) set of points “points_distorted” with respect to the
(perfect) set of points “points_perfect” in a least-square sense.
:param points_distorted: List of points describing a given (distorted) polyhedron for which the rotation that

> aligns these points in a least-square sense to the set of perfect points “points_perfect”


* **Parameters**

    **points_perfect** – List of “perfect” points describing a given model polyhedron.



* **Returns**

    The rotation matrix.



### find_scaling_factor(points_distorted, points_perfect, rot)
This finds the scaling factor between the (distorted) set of points “points_distorted” and the
(perfect) set of points “points_perfect” in a least-square sense.
:param points_distorted: List of points describing a given (distorted) polyhedron for which the scaling factor has

> to be obtained.


* **Parameters**


    * **points_perfect** – List of “perfect” points describing a given model polyhedron.


    * **rot** – The rotation matrix



* **Returns**

    The scaling factor between the two structures and the rotated set of (distorted) points.



### symmetry_measure(points_distorted, points_perfect)
Computes the continuous symmetry measure of the (distorted) set of points “points_distorted” with respect to the
(perfect) set of points “points_perfect”.
:param points_distorted: List of points describing a given (distorted) polyhedron for which the symmetry measure

> has to be computed with respect to the model polyhedron described by the list of points
> “points_perfect”.


* **Parameters**

    **points_perfect** – List of “perfect” points describing a given model polyhedron.



* **Returns**

    The continuous symmetry measure of the distorted polyhedron with respect to the perfect polyhedron.


## pymatgen.analysis.chemenv.coordination_environments.structure_environments module

This module contains objects that are used to describe the environments in a structure. The most detailed object
(StructureEnvironments) contains a very thorough analysis of the environments of a given atom but is difficult to
used as such. The LightStructureEnvironments object is a lighter version that is obtained by applying a “strategy”
on the StructureEnvironments object. Basically, the LightStructureEnvironments provides the coordination environment(s)
and possibly some fraction corresponding to these.


### _class_ ChemicalEnvironments(coord_geoms=None)
Bases: `MSONable`

Class used to store all the information about the chemical environment of a given site for a given list of
coordinated neighbors (internally called “cn_map”).

Initializes the ChemicalEnvironments object containing all the information about the chemical
environment of a given site.


* **Parameters**

    **coord_geoms** – coordination geometries to be added to the chemical environment.



#### add_coord_geom(mp_symbol, symmetry_measure, algo='UNKNOWN', permutation=None, override=False, local2perfect_map=None, perfect2local_map=None, detailed_voronoi_index=None, other_symmetry_measures=None, rotation_matrix=None, scaling_factor=None)
Adds a coordination geometry to the ChemicalEnvironments object.


* **Parameters**


    * **mp_symbol** – Symbol of the coordination geometry added.


    * **symmetry_measure** – Symmetry measure of the coordination geometry added.


    * **algo** – Algorithm used for the search of the coordination geometry added.


    * **permutation** – Permutation of the neighbors that leads to the csm stored.


    * **override** – If set to True, the coordination geometry will override the existent one if present.


    * **local2perfect_map** – Mapping of the local indices to the perfect indices.


    * **perfect2local_map** – Mapping of the perfect indices to the local indices.


    * **detailed_voronoi_index** – Index in the voronoi containing the neighbors set.


    * **other_symmetry_measures** – Other symmetry measure of the coordination geometry added (with/without the
    central atom, centered on the central atom or on the centroid with/without the central atom).


    * **rotation_matrix** – Rotation matrix mapping the local geometry to the perfect geometry.


    * **scaling_factor** – Scaling factor mapping the local geometry to the perfect geometry.



* **Raises**

    **ChemenvError if the coordination geometry is already added and override is set to False** –



#### as_dict()
Returns a dictionary representation of the ChemicalEnvironments object.


* **Returns**

    A dictionary representation of the ChemicalEnvironments object.



#### _classmethod_ from_dict(d)
Reconstructs the ChemicalEnvironments object from a dict representation of the ChemicalEnvironments created
using the as_dict method.


* **Parameters**

    **d** – dict representation of the ChemicalEnvironments object.



* **Returns**

    ChemicalEnvironments object.



#### is_close_to(other, rtol=0.0, atol=1e-08)
Whether this ChemicalEnvironments object is close to another one.


* **Parameters**


    * **other** – Another ChemicalEnvironments object.


    * **rtol** – Relative tolerance for the comparison of Continuous Symmetry Measures.


    * **atol** – Absolute tolerance for the comparison of Continuous Symmetry Measures.



* **Returns**

    True if the two ChemicalEnvironments objects are close to each other.



* **Return type**

    bool



#### minimum_geometries(n=None, symmetry_measure_type=None, max_csm=None)
Returns a list of geometries with increasing continuous symmetry measure in this ChemicalEnvironments object.


* **Parameters**

    **n** – Number of geometries to be included in the list.



* **Returns**

    List of geometries with increasing continuous symmetry measure in this ChemicalEnvironments object.



* **Raises**

    **ValueError if no coordination geometry is found in this ChemicalEnvironments object.** –



#### minimum_geometry(symmetry_measure_type=None, max_csm=None)
Returns the geometry with the minimum continuous symmetry measure of this ChemicalEnvironments.


* **Returns**

    tuple (symbol, csm) with symbol being the geometry with the minimum continuous symmetry measure and
    csm being the continuous symmetry measure associated to it.



* **Raises**

    **ValueError if no coordination geometry is found in this ChemicalEnvironments object.** –



### _class_ LightStructureEnvironments(strategy, coordination_environments=None, all_nbs_sites=None, neighbors_sets=None, structure=None, valences=None, valences_origin=None)
Bases: `MSONable`

Class used to store the chemical environments of a given structure obtained from a given ChemenvStrategy. Currently,
only strategies leading to the determination of a unique environment for each site is allowed
This class does not store all the information contained in the StructureEnvironments object, only the coordination
environment found.

Constructor for the LightStructureEnvironments object.


* **Parameters**


    * **strategy** – ChemEnv strategy used to get the environments.


    * **coordination_environments** – The coordination environments identified.


    * **all_nbs_sites** – All the possible neighbors for each site in the structure.


    * **neighbors_sets** – The neighbors sets of each site in the structure.


    * **structure** – The structure.


    * **valences** – The valences used to get the environments (if needed).


    * **valences_origin** – How the valences were obtained (e.g. from the Bond-valence analysis or from the original
    structure).



#### DEFAULT_STATISTICS_FIELDS(_ = ('anion_list', 'anion_atom_list', 'cation_list', 'cation_atom_list', 'neutral_list', 'neutral_atom_list', 'atom_coordination_environments_present', 'ion_coordination_environments_present', 'fraction_atom_coordination_environments_present', 'fraction_ion_coordination_environments_present', 'coordination_environments_atom_present', 'coordination_environments_ion_present'_ )

#### DELTA_MAX_OXIDATION_STATE(_ = 0._ )

#### _class_ NeighborsSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), isite, all_nbs_sites, all_nbs_sites_indices)
Bases: `object`

Class used to store a given set of neighbors of a given site (based on a list of sites, the voronoi
container is not part of the LightStructureEnvironments object).

Constructor for NeighborsSet.


* **Parameters**


    * **structure** – Structure object.


    * **isite** – Index of the site for which neighbors are stored in this NeighborsSet.


    * **all_nbs_sites** – All the possible neighbors for this site.


    * **all_nbs_sites_indices** – Indices of the sites in all_nbs_sites that make up this NeighborsSet.



#### as_dict()
A JSON-serializable dict representation of the NeighborsSet.


#### _classmethod_ from_dict(dct, structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), all_nbs_sites)
Reconstructs the NeighborsSet algorithm from its JSON-serializable dict representation, together with
the structure and all the possible neighbors sites.

As an inner (nested) class, the NeighborsSet is not supposed to be used anywhere else that inside the
LightStructureEnvironments. The from_dict method is thus using the structure and all_nbs_sites when
reconstructing itself. These two are both in the LightStructureEnvironments object.


* **Parameters**


    * **dct** – a JSON-serializable dict representation of a NeighborsSet.


    * **structure** – The structure.


    * **all_nbs_sites** – The list of all the possible neighbors for a given site.



* **Returns**

    NeighborsSet



#### _property_ neighb_coords()
Coordinates of neighbors for this NeighborsSet.


#### _property_ neighb_indices_and_images(_: list[dict[str, int]_ )
List of indices and images with respect to the original unit cell sites for this NeighborsSet.


#### _property_ neighb_sites()
Neighbors for this NeighborsSet as pymatgen Sites.


#### _property_ neighb_sites_and_indices()
List of neighbors for this NeighborsSet as pymatgen Sites and their index in the original structure.


#### as_dict()

* **Returns**

    Bson-serializable representation of the LightStructureEnvironments object.



* **Return type**

    dict



#### clear_environments(conditions=None)
Get the clear environments in the structure.


* **Parameters**

    **conditions** – Conditions to be checked for an environment to be “clear”.



* **Returns**

    Clear environments in this structure.



* **Return type**

    list



#### contains_only_one_anion(anion)
Whether this LightStructureEnvironments concerns a structure with only one given anion type.


* **Parameters**

    **anion** – Anion (e.g. O2-, …).



* **Returns**

    True if this LightStructureEnvironments concerns a structure with only one given anion.



* **Return type**

    bool



#### contains_only_one_anion_atom(anion_atom)
Whether this LightStructureEnvironments concerns a structure with only one given anion atom type.


* **Parameters**

    **anion_atom** – Anion (e.g. O, …). The structure could contain O2- and O- though.



* **Returns**

    True if this LightStructureEnvironments concerns a structure with only one given anion_atom.



* **Return type**

    bool



#### environments_identified()
Return the set of environments identified in this structure.


* **Returns**

    environments identified in this structure.



* **Return type**

    set



#### _classmethod_ from_dict(d)
Reconstructs the LightStructureEnvironments object from a dict representation of the
LightStructureEnvironments created using the as_dict method.


* **Parameters**

    **d** – dict representation of the LightStructureEnvironments object.



* **Returns**

    LightStructureEnvironments object.



#### _classmethod_ from_structure_environments(strategy, structure_environments, valences=None, valences_origin=None)
Construct a LightStructureEnvironments object from a strategy and a StructureEnvironments object.


* **Parameters**


    * **strategy** – ChemEnv strategy used.


    * **structure_environments** – StructureEnvironments object from which to construct the LightStructureEnvironments.


    * **valences** – The valences of each site in the structure.


    * **valences_origin** – How the valences were obtained (e.g. from the Bond-valence analysis or from the original
    structure).



* **Returns**

    LightStructureEnvironments



#### get_site_info_for_specie_allces(specie, min_fraction=0)
Get list of indices that have the given specie.


* **Parameters**


    * **specie** – Species to get.


    * **min_fraction** – Minimum fraction of the coordination environment.



* **Returns**

    with the list of coordination environments for the given species, the indices of the sites

        in which they appear, their fractions and continuous symmetry measures.




* **Return type**

    dict



#### get_site_info_for_specie_ce(specie, ce_symbol)
Get list of indices that have the given specie with a given Coordination environment.


* **Parameters**


    * **specie** – Species to get.


    * **ce_symbol** – Symbol of the coordination environment to get.



* **Returns**

    Keys are ‘isites’, ‘fractions’, ‘csms’ which contain list of indices in the structure

        that have the given specie in the given environment, their fraction and continuous
        symmetry measures.




* **Return type**

    dict



#### get_statistics(statistics_fields=('anion_list', 'anion_atom_list', 'cation_list', 'cation_atom_list', 'neutral_list', 'neutral_atom_list', 'atom_coordination_environments_present', 'ion_coordination_environments_present', 'fraction_atom_coordination_environments_present', 'fraction_ion_coordination_environments_present', 'coordination_environments_atom_present', 'coordination_environments_ion_present'), bson_compatible=False)
Get the statistics of environments for this structure.


* **Parameters**


    * **statistics_fields** – Which statistics to get.


    * **bson_compatible** – Whether to make the dictionary BSON-compatible.



* **Returns**

    with the requested statistics.



* **Return type**

    dict



#### setup_statistic_lists()
Set up the statistics of environments for this LightStructureEnvironments.


#### site_contains_environment(isite, ce_symbol)
Whether a given site contains a given coordination environment.


* **Parameters**


    * **isite** – Index of the site.


    * **ce_symbol** – Symbol of the coordination environment.



* **Returns**

    True if the site contains the given coordination environment.



* **Return type**

    bool



#### site_has_clear_environment(isite, conditions=None)
Whether a given site has a “clear” environments.

A “clear” environment is somewhat arbitrary. You can pass (multiple) conditions, e.g. the environment should
have a continuous symmetry measure lower than this, a fraction higher than that, …


* **Parameters**


    * **isite** – Index of the site.


    * **conditions** – Conditions to be checked for an environment to be “clear”.



* **Returns**

    True if the site has a clear environment.



* **Return type**

    bool



#### structure_contains_atom_environment(atom_symbol, ce_symbol)
Checks whether the structure contains a given atom in a given environment.


* **Parameters**


    * **atom_symbol** – Symbol of the atom.


    * **ce_symbol** – Symbol of the coordination environment.



* **Returns**

    True if the coordination environment is found for the given atom.



* **Return type**

    bool



#### structure_has_clear_environments(conditions=None, skip_none=True, skip_empty=False)
Whether all sites in a structure have “clear” environments.


* **Parameters**


    * **conditions** – Conditions to be checked for an environment to be “clear”.


    * **skip_none** – Whether to skip sites for which no environments have been computed.


    * **skip_empty** – Whether to skip sites for which no environments could be found.



* **Returns**

    True if all the sites in the structure have clear environments.



* **Return type**

    bool



#### _property_ uniquely_determines_coordination_environments()
True if the coordination environments are uniquely determined.


### _class_ StructureEnvironments(voronoi, valences, sites_map, equivalent_sites, ce_list, structure, neighbors_sets=None, info=None)
Bases: `MSONable`

Class used to store the chemical environments of a given structure.

Constructor for the StructureEnvironments object.


* **Parameters**


    * **voronoi** – VoronoiContainer object for the structure.


    * **valences** – Valences provided.


    * **sites_map** – Mapping of equivalent sites to the unequivalent sites that have been computed.


    * **equivalent_sites** – List of list of equivalent sites of the structure.


    * **ce_list** – List of chemical environments.


    * **structure** – Structure object.


    * **neighbors_sets** – List of neighbors sets.


    * **info** – Additional information for this StructureEnvironments object.



#### AC(_ = <pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions object_ )

#### _class_ NeighborsSet(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), isite, detailed_voronoi, site_voronoi_indices, sources=None)
Bases: `object`

Class used to store a given set of neighbors of a given site (based on the detailed_voronoi).

Constructor for NeighborsSet.


* **Parameters**


    * **structure** – Structure object.


    * **isite** – Index of the site for which neighbors are stored in this NeighborsSet.


    * **detailed_voronoi** – Corresponding DetailedVoronoiContainer object containing all the possible
    neighbors of the give site.


    * **site_voronoi_indices** – Indices of the voronoi sites in the DetailedVoronoiContainer object that
    make up this NeighborsSet.


    * **sources** – Sources for this NeighborsSet, i.e. how this NeighborsSet was generated.



#### add_source(source)
Add a source to this NeighborsSet.


* **Parameters**

    **source** – Information about the generation of this NeighborsSet.



#### angle_plateau()
Returns the angles plateau’s for this NeighborsSet.


#### _property_ angles()
Angles for each neighbor in this NeighborsSet.


#### as_dict()
A JSON-serializable dict representation of the NeighborsSet.


#### _property_ coords()
Coordinates of the current central atom and its neighbors for this NeighborsSet.


#### distance_plateau()
Returns the distances plateau’s for this NeighborsSet.


#### _property_ distances()
Distances to each neighbor in this NeighborsSet.


#### _classmethod_ from_dict(dct, structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), detailed_voronoi)
Reconstructs the NeighborsSet algorithm from its JSON-serializable dict representation, together with
the structure and the DetailedVoronoiContainer.

As an inner (nested) class, the NeighborsSet is not supposed to be used anywhere else that inside the
StructureEnvironments. The from_dict method is thus using the structure and  detailed_voronoi when
reconstructing itself. These two are both in the StructureEnvironments object.


* **Parameters**


    * **dct** – a JSON-serializable dict representation of a NeighborsSet.


    * **structure** – The structure.


    * **detailed_voronoi** – The Voronoi object containing all the neighboring atoms from which the subset of
    neighbors for this NeighborsSet is extracted.



* **Returns**

    NeighborsSet



#### get_neighb_voronoi_indices(permutation)
Get indices in the detailed_voronoi corresponding to the current permutation.


* **Parameters**

    **permutation** – Current permutation for which the indices in the detailed_voronoi are needed.



* **Returns**

    indices in the detailed_voronoi.



* **Return type**

    list[int]



#### _property_ info()
Summarized information about this NeighborsSet.


#### _property_ neighb_coords()
Coordinates of neighbors for this NeighborsSet.


#### _property_ neighb_coordsOpt()
Optimized access to the coordinates of neighbors for this NeighborsSet.


#### _property_ neighb_sites()
Neighbors for this NeighborsSet as pymatgen Sites.


#### _property_ neighb_sites_and_indices()
List of neighbors for this NeighborsSet as pymatgen Sites and their index in the original structure.


#### _property_ normalized_angles()
Normalized angles for each neighbor in this NeighborsSet.


#### _property_ normalized_distances()
Normalized distances to each neighbor in this NeighborsSet.


#### _property_ source()
Returns the source of this NeighborsSet (how it was generated, e.g. from which Voronoi
cutoffs, or from hints).


#### voronoi_grid_surface_points(additional_condition=1, other_origins='DO_NOTHING')
Get the surface points in the Voronoi grid for this neighbor from the sources.
The general shape of the points should look like a staircase such as in the following figure :

> ^

0.0|

    B—-C|    ||    |a  |      k    D——-E
n  |      |            |
g  |      |            |
l  |      |            |
e  |      j            F—-n———G

> |                           ||                           |A—-g——-h—-i———H1.0+————————————————->

    1.0              distance              2.0   ->+Inf


* **Parameters**


    * **additional_condition** – Additional condition for the neighbors.


    * **other_origins** – What to do with sources that do not come from the Voronoi grid (e.g. “from hints”).



#### add_neighbors_set(isite, nb_set)
Adds a neighbor set to the list of neighbors sets for this site.


* **Parameters**


    * **isite** – Index of the site under consideration.


    * **nb_set** – NeighborsSet to be added.



#### as_dict()
Bson-serializable dict representation of the StructureEnvironments object.


* **Returns**

    Bson-serializable dict representation of the StructureEnvironments object.



#### differences_wrt(other)
Return differences found in the current StructureEnvironments with respect to another StructureEnvironments.


* **Parameters**

    **other** – A StructureEnvironments object.



* **Returns**

    List of differences between the two StructureEnvironments objects.



#### _classmethod_ from_dict(d)
Reconstructs the StructureEnvironments object from a dict representation of the StructureEnvironments created
using the as_dict method.


* **Parameters**

    **d** – dict representation of the StructureEnvironments object.



* **Returns**

    StructureEnvironments object.



#### get_coordination_environments(isite, cn, nb_set)
Get the ChemicalEnvironments for a given site, coordination and neighbors set.


* **Parameters**


    * **isite** – Index of the site for which the ChemicalEnvironments is looked for.


    * **cn** – Coordination for which the ChemicalEnvironments is looked for.


    * **nb_set** – Neighbors set for which the ChemicalEnvironments is looked for.



* **Returns**

    ChemicalEnvironments



#### get_csm(isite, mp_symbol)
Get the continuous symmetry measure for a given site in the given coordination environment.


* **Parameters**


    * **isite** – Index of the site.


    * **mp_symbol** – Symbol of the coordination environment for which we want the continuous symmetry measure.



* **Returns**

    Continuous symmetry measure of the given site in the given environment.



#### get_csm_and_maps(isite, max_csm=8.0, figsize=None, symmetry_measure_type=None)
Plotting of the coordination numbers of a given site for all the distfactor/angfactor parameters. If the
chemical environments are given, a color map is added to the plot, with the lowest continuous symmetry measure
as the value for the color of that distfactor/angfactor set.


* **Parameters**


    * **isite** – Index of the site for which the plot has to be done.


    * **max_csm** – Maximum continuous symmetry measure to be shown.


    * **figsize** – Size of the figure.


    * **symmetry_measure_type** – Type of continuous symmetry measure to be used.



* **Returns**

    Matplotlib figure and axes representing the CSM and maps.



#### get_csms(isite, mp_symbol)
Returns the continuous symmetry measure(s) of site with index isite with respect to the
perfect coordination environment with mp_symbol. For some environments, a given mp_symbol might not
be available (if there is no voronoi parameters leading to a number of neighbors corresponding to
the coordination number of environment mp_symbol). For some environments, a given mp_symbol might
lead to more than one csm (when two or more different voronoi parameters lead to different neighbors
but with same number of neighbors).


* **Parameters**


    * **isite** – Index of the site.


    * **mp_symbol** – MP symbol of the perfect environment for which the csm has to be given.



* **Returns**

    for site isite with respect to geometry mp_symbol



* **Return type**

    list[CSM]



#### get_environments_figure(isite, plot_type=None, title='Coordination numbers', max_dist=2.0, colormap=None, figsize=None, strategy=None)
Plotting of the coordination environments of a given site for all the distfactor/angfactor regions. The
chemical environments with the lowest continuous symmetry measure is shown for each distfactor/angfactor
region as the value for the color of that distfactor/angfactor region (using a colormap).


* **Parameters**


    * **isite** – Index of the site for which the plot has to be done.


    * **plot_type** – How to plot the coordinations.


    * **title** – Title for the figure.


    * **max_dist** – Maximum distance to be plotted when the plotting of the distance is set to ‘initial_normalized’
    or ‘initial_real’ (Warning: this is not the same meaning in both cases! In the first case, the
    closest atom lies at a “normalized” distance of 1.0 so that 2.0 means refers to this normalized
    distance while in the second case, the real distance is used).


    * **colormap** – Color map to be used for the continuous symmetry measure.


    * **figsize** – Size of the figure.


    * **strategy** – Whether to plot information about one of the Chemenv Strategies.



* **Returns**

    matplotlib figure and axes representing the environments.



* **Return type**

    tuple[plt.Figure, plt.Axes]



#### init_neighbors_sets(isite, additional_conditions=None, valences=None)
Initialize the list of neighbors sets for the current site.


* **Parameters**


    * **isite** – Index of the site under consideration.


    * **additional_conditions** – Additional conditions to be used for the initialization of the list of
    neighbors sets, e.g. “Only anion-cation bonds”, …


    * **valences** – List of valences for each site in the structure (needed if an additional condition based on the
    valence is used, e.g. only anion-cation bonds).



#### plot_csm_and_maps(isite, max_csm=8.0)
Plotting of the coordination numbers of a given site for all the distfactor/angfactor parameters. If the
chemical environments are given, a color map is added to the plot, with the lowest continuous symmetry measure
as the value for the color of that distfactor/angfactor set.


* **Parameters**


    * **isite** – Index of the site for which the plot has to be done


    * **max_csm** – Maximum continuous symmetry measure to be shown.



#### plot_environments(isite, plot_type=None, title='Coordination numbers', max_dist=2.0, figsize=None, strategy=None)
Plotting of the coordination numbers of a given site for all the distfactor/angfactor parameters. If the
chemical environments are given, a color map is added to the plot, with the lowest continuous symmetry measure
as the value for the color of that distfactor/angfactor set.


* **Parameters**


    * **isite** – Index of the site for which the plot has to be done.


    * **plot_type** – How to plot the coordinations.


    * **title** – Title for the figure.


    * **max_dist** – Maximum distance to be plotted when the plotting of the distance is set to ‘initial_normalized’
    or ‘initial_real’ (Warning: this is not the same meaning in both cases! In the first case, the
    closest atom lies at a “normalized” distance of 1.0 so that 2.0 means refers to this normalized
    distance while in the second case, the real distance is used).


    * **figsize** – Size of the figure.


    * **strategy** – Whether to plot information about one of the Chemenv Strategies.



#### save_environments_figure(isite, imagename='image.png', plot_type=None, title='Coordination numbers', max_dist=2.0, figsize=None)
Saves the environments figure to a given file.


* **Parameters**


    * **isite** – Index of the site for which the plot has to be done.


    * **imagename** – Name of the file to which the figure has to be saved.


    * **plot_type** – How to plot the coordinations.


    * **title** – Title for the figure.


    * **max_dist** – Maximum distance to be plotted when the plotting of the distance is set to ‘initial_normalized’
    or ‘initial_real’ (Warning: this is not the same meaning in both cases! In the first case, the
    closest atom lies at a “normalized” distance of 1.0 so that 2.0 means refers to this normalized
    distance while in the second case, the real distance is used).


    * **figsize** – Size of the figure.



#### update_coordination_environments(isite, cn, nb_set, ce)
Updates the coordination environment for this site, coordination and neighbor set.


* **Parameters**


    * **isite** – Index of the site to be updated.


    * **cn** – Coordination to be updated.


    * **nb_set** – Neighbors set to be updated.


    * **ce** – ChemicalEnvironments object for this neighbors set.



#### update_site_info(isite, info_dict)
Update information about this site.


* **Parameters**


    * **isite** – Index of the site for which info has to be updated.


    * **info_dict** – Dictionary of information to be added for this site.


## pymatgen.analysis.chemenv.coordination_environments.voronoi module

This module contains the object used to describe the possible bonded atoms based on a Voronoi analysis.


### _class_ DetailedVoronoiContainer(structure=None, voronoi_list2=None, voronoi_cutoff=10.0, isites=None, normalized_distance_tolerance=1e-05, normalized_angle_tolerance=0.001, additional_conditions=None, valences=None, maximum_distance_factor=None, minimum_angle_factor=None)
Bases: `MSONable`

Class used to store the full Voronoi of a given structure.

Constructor for the VoronoiContainer object. Either a structure is given, in which case the Voronoi is
computed, or the different components of the VoronoiContainer are given (used in the from_dict method).


* **Parameters**


    * **structure** – Structure for which the Voronoi is computed.


    * **voronoi_list2** – List of voronoi polyhedrons for each site.


    * **voronoi_cutoff** – cutoff used for the voronoi.


    * **isites** – indices of sites for which the Voronoi has to be computed.


    * **normalized_distance_tolerance** – Tolerance for two normalized distances to be considered equal.


    * **normalized_angle_tolerance** – Tolerance for two normalized angles to be considered equal.


    * **additional_conditions** – Additional conditions to be used.


    * **valences** – Valences of all the sites in the structure (used when additional conditions require it).


    * **maximum_distance_factor** – The maximum distance factor to be considered.


    * **minimum_angle_factor** – The minimum angle factor to be considered.



* **Raises**

    **RuntimeError if the Voronoi cannot be constructed.** –



#### AC(_ = <pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions object_ )

#### _static_ _get_vertices_dist_ang_indices(parameter_indices_list)

#### _precompute_additional_conditions(ivoronoi, voronoi, valences)

#### _precompute_angle_conditions(ivoronoi, voronoi)

#### _precompute_distance_conditions(ivoronoi, voronoi)

#### as_dict()
Bson-serializable dict representation of the VoronoiContainer.


* **Returns**

    dictionary that is BSON-encodable.



#### default_normalized_angle_tolerance(_ = 0.00_ )

#### default_normalized_distance_tolerance(_ = 1e-0_ )

#### default_voronoi_cutoff(_ = 10._ )

#### _classmethod_ from_dict(dct)
Reconstructs the VoronoiContainer object from a dict representation of the VoronoiContainer created using
the as_dict method.


* **Parameters**

    **dct** – dict representation of the VoronoiContainer object.



* **Returns**

    VoronoiContainer object.



#### get_rdf_figure(isite, normalized=True, figsize=None, step_function=None)
Get the Radial Distribution Figure for a given site.


* **Parameters**


    * **isite** – Index of the site.


    * **normalized** – Whether to normalize distances.


    * **figsize** – Size of the figure.


    * **step_function** – Type of step function to be used for the RDF.



* **Returns**

    Matplotlib figure.



* **Return type**

    plt.figure



#### get_sadf_figure(isite, normalized=True, figsize=None, step_function=None)
Get the Solid Angle Distribution Figure for a given site.


* **Parameters**


    * **isite** – Index of the site.


    * **normalized** – Whether to normalize angles.


    * **figsize** – Size of the figure.


    * **step_function** – Type of step function to be used for the SADF.



* **Returns**

    matplotlib figure.



* **Return type**

    plt.figure



#### is_close_to(other, rtol=0.0, atol=1e-08)
Whether two DetailedVoronoiContainer objects are close to each other.


* **Parameters**


    * **other** – Another DetailedVoronoiContainer to be compared with.


    * **rtol** – Relative tolerance to compare values.


    * **atol** – Absolute tolerance to compare values.



* **Returns**

    True if the two DetailedVoronoiContainer are close to each other.



* **Return type**

    bool



#### maps_and_surfaces(isite, surface_calculation_type=None, max_dist=2.0, additional_conditions=None)
Get the different surfaces and their cn_map corresponding to the different distance-angle cutoffs
for a given site.


* **Parameters**


    * **isite** – Index of the site


    * **surface_calculation_type** – How to compute the surface.


    * **max_dist** – The maximum distance factor to be considered.


    * **additional_conditions** – If additional conditions have to be considered.



* **Returns**

    Surfaces and cn_map’s for each distance-angle cutoff.



#### maps_and_surfaces_bounded(isite, surface_calculation_options=None, additional_conditions=None)
Get the different surfaces (using boundaries) and their cn_map corresponding to the different
distance-angle cutoffs for a given site.


* **Parameters**


    * **isite** – Index of the site


    * **surface_calculation_options** – Options for the boundaries.


    * **additional_conditions** – If additional conditions have to be considered.



* **Returns**

    Surfaces and cn_map’s for each distance-angle cutoff.



#### neighbors(isite, distfactor, angfactor, additional_condition=None)
Get the neighbors of a given site corresponding to a given distance and angle factor.


* **Parameters**


    * **isite** – Index of the site.


    * **distfactor** – Distance factor.


    * **angfactor** – Angle factor.


    * **additional_condition** – Additional condition to be used (currently not implemented).



* **Returns**

    List of neighbors of the given site for the given distance and angle factors.



#### neighbors_surfaces(isite, surface_calculation_type=None, max_dist=2.0)
Get the different surfaces corresponding to the different distance-angle cutoffs for a given site.


* **Parameters**


    * **isite** – Index of the site


    * **surface_calculation_type** – How to compute the surface.


    * **max_dist** – The maximum distance factor to be considered.



* **Returns**

    Surfaces for each distance-angle cutoff.



#### neighbors_surfaces_bounded(isite, surface_calculation_options=None)
Get the different surfaces (using boundaries) corresponding to the different distance-angle cutoffs
for a given site.


* **Parameters**


    * **isite** – Index of the site.


    * **surface_calculation_options** – Options for the boundaries.



* **Returns**

    Surfaces for each distance-angle cutoff.



#### setup_neighbors_distances_and_angles(indices)
Initializes the angle and distance separations.


* **Parameters**

    **indices** – Indices of the sites for which the Voronoi is needed.



#### setup_voronoi_list(indices, voronoi_cutoff)
Set up of the voronoi list of neighbors by calling qhull.


* **Parameters**


    * **indices** – indices of the sites for which the Voronoi is needed.


    * **voronoi_cutoff** – Voronoi cutoff for the search of neighbors.



* **Raises**

    **RuntimeError** – If an infinite vertex is found in the voronoi construction.



#### to_bson_voronoi_list2()
Transforms the voronoi_list into a vlist + bson_nb_voro_list, that are BSON-encodable.


* **Returns**

    [vlist, bson_nb_voro_list], to be used in the as_dict method.



#### voronoi_parameters_bounds_and_limits(isite, plot_type, max_dist)
Get the different boundaries and limits of the distance and angle factors for the given site.


* **Parameters**


    * **isite** – Index of the site.


    * **plot_type** – Types of distance/angle parameters to get.


    * **max_dist** – Maximum distance factor.



* **Returns**

    Distance and angle bounds and limits.



### from_bson_voronoi_list2(bson_nb_voro_list2, structure)
Returns the voronoi_list needed for the VoronoiContainer object from a bson-encoded voronoi_list.


* **Parameters**


    * **bson_nb_voro_list2** – List of periodic sites involved in the Voronoi.


    * **structure** – Structure object.



* **Returns**

    The voronoi_list needed for the VoronoiContainer (with PeriodicSites as keys of the dictionary - not
    allowed in the BSON format).