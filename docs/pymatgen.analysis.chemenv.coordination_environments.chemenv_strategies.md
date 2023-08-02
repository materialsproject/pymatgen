---
layout: default
title: pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies module

This module provides so-called “strategies” to determine the coordination environments of an atom in a structure.
Some strategies can favour larger or smaller environments. Some strategies uniquely identifies the environments while
some others can identify the environment as a “mix” of several environments, each of which is assigned with a given
fraction. The choice of the strategy depends on the purpose of the user.


### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy(structure_environments=None, symmetry_measure_type='csm_wcs_ctwcc')
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

#### _abstract_ as_dict()
Bson-serializable dict representation of the SimplestChemenvStrategy object.
:return: Bson-serializable dict representation of the SimplestChemenvStrategy object.


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
:return: StructureEnvironments object.


#### get_site_ce_fractions_and_neighbors(site, full_ce_info=False, strategy_info=False)
Applies the strategy to the structure_environments object in order to get coordination environments, their
fraction, csm, geometry_info, and neighbors
:param site: Site for which the above information is seeked
:return: The list of neighbors of the site. For complex strategies, where one allows multiple solutions, this
can return a list of list of neighbors.


#### _abstract_ get_site_coordination_environment(site)
Applies the strategy to the structure_environments object in order to define the coordination environment of
a given site.
:param site: Site for which the coordination environment is looked for
:return: The coordination environment of the site. For complex strategies, where one allows multiple

> solutions, this can return a list of coordination environments for the site.


#### _abstract_ get_site_coordination_environments(site)
Applies the strategy to the structure_environments object in order to define the coordination environment of
a given site.
:param site: Site for which the coordination environment is looked for
:return: The coordination environment of the site. For complex strategies, where one allows multiple

> solutions, this can return a list of coordination environments for the site.


#### _abstract_ get_site_coordination_environments_fractions(site, isite=None, dequivsite=None, dthissite=None, mysym=None, ordered=True, min_fraction=0, return_maps=True, return_strategy_dict_info=False)
Applies the strategy to the structure_environments object in order to define the coordination environment of
a given site.
:param site: Site for which the coordination environment is looked for
:return: The coordination environment of the site. For complex strategies, where one allows multiple

> solutions, this can return a list of coordination environments for the site.


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



* **Returns**

    None



#### set_structure_environments(structure_environments)
Set the structure environments to this strategy.


* **Parameters**

    **structure_environments** – StructureEnvironments object.



* **Returns**

    None



#### setup_options(all_options_dict)
Set up options for this strategy based on a dict.


* **Parameters**

    **all_options_dict** – Dict of option_name->option_value.



* **Returns**

    None



#### _property_ symmetry_measure_type()
Type of symmetry measure.


#### _property_ uniquely_determines_coordination_environments()
Returns True if the strategy leads to a unique coordination environment, False otherwise.
:return: True if the strategy leads to a unique coordination environment, False otherwise.


### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt(integer)
Bases: `int`, `StrategyOption`

Integer representing an additional condition in a strategy.

Special int representing additional conditions.


#### allowed_values(_: str | Non_ _ = "Integer amongst :\\n - 0 for 'No additional condition'\\n - 1 for 'Only anion-cation bonds'\\n - 2 for 'No element-element bonds (same elements)'\\n - 3 for 'Only anion-cation bonds and no element-element bonds (same elements)'\\n - 4 for 'Only element-oxygen bonds'\\n_ )

#### as_dict()
MSONable dict.


#### description(_ = 'Only element-oxygen bonds_ )

#### _classmethod_ from_dict(dct)
Initialize additional condition from dict.


* **Parameters**

    **d** – Dict representation of the additional condition.



#### integer(_ = _ )

### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat(myfloat)
Bases: `float`, `StrategyOption`

Angle cutoff in a strategy.

Special float that should be between 0 and 1.


* **Parameters**

    **myfloat** – Angle cutoff.



#### allowed_values(_: str | Non_ _ = 'Real number between 0 and 1_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(d)
Initialize angle cutoff from dict.


* **Parameters**

    **d** – Dict representation of the angle cutoff.



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight(aa=1)
Bases: `NbSetWeight`

Weight of neighbors set based on the angle.

Initialize AngleNbSetWeight estimator.


* **Parameters**

    **aa** – Exponent of the angle for the estimator.



#### SHORT_NAME(_ = 'AngleWeight_ )

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


#### _classmethod_ from_dict(dd)
From dict
:param dd:
:return:


#### weight(nb_set, structure_environments, cn_map=None, additional_info=None)
Get the weight of a given neighbors set.


* **Parameters**


    * **nb_set** – Neighbors set.


    * **structure_environments** – Structure environments used to estimate weight.


    * **cn_map** – Mapping index for this neighbors set.


    * **additional_info** – Additional information.



* **Returns**

    Weight of the neighbors set.



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight(angle_function=None, weight_function=None)
Bases: `NbSetWeight`

Weight of neighbors set based on the angle.

Initialize AnglePlateauNbSetWeight.


* **Parameters**


    * **angle_function** – Angle function to use.


    * **weight_function** – Ratio function to use.



#### SHORT_NAME(_ = 'AnglePlateauWeight_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dd)
Initialize from dict.


* **Parameters**

    **dd** – Dict representation of AnglePlateauNbSetWeight.



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



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight(cn_weights, initialization_options)
Bases: `NbSetWeight`

Weight of neighbors set based on specific biases towards specific coordination numbers.

Initialize CNBiasNbSetWeight.


* **Parameters**


    * **cn_weights** – Weights for each coordination.


    * **initialization_options** – Options for initialization.



#### SHORT_NAME(_ = 'CNBiasWeight_ )

#### as_dict()
MSONable dict.


#### _classmethod_ explicit(cn_weights)
Initializes weights explicitly for each coordination.


* **Parameters**

    **cn_weights** – Weights for each coordination.



* **Returns**

    CNBiasNbSetWeight.



#### _classmethod_ from_description(dd)
Initializes weights from description.


* **Parameters**

    **dd** – Dictionary description.



* **Returns**

    CNBiasNbSetWeight.



#### _classmethod_ from_dict(dd)
Initialize from dict.


* **Parameters**

    **dd** – Dict representation of CNBiasNbSetWeight.



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



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat(myfloat)
Bases: `float`, `StrategyOption`

Real number representing a Continuous Symmetry Measure.

Special float that should be between 0 and 100.


* **Parameters**

    **myfloat** – CSM.



#### allowed_values(_: str | Non_ _ = 'Real number between 0 and 100_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dct)
Initialize CSM from dict.


* **Parameters**

    **d** – Dict representation of the CSM.



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight(effective_csm_estimator={'function': 'power2_inverse_decreasing', 'options': {'max_csm': 8.0}}, weight_estimator={'function': 'smootherstep', 'options': {'delta_csm_max': 3.0, 'delta_csm_min': 0.5}}, delta_cn_weight_estimators=None, symmetry_measure_type='csm_wcs_ctwcc')
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

#### as_dict()
MSONable dict.
:return:


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



#### _classmethod_ from_dict(dd)
Initialize from dict.


* **Parameters**

    **dd** – Dict representation of DeltaCSMNbSetWeight.



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



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight(weight_function=None, nbs_source='voronoi')
Bases: `NbSetWeight`

Weight of neighbors set based on the difference of distances.

Initialize DeltaDistanceNbSetWeight.


* **Parameters**


    * **weight_function** – Ratio function to use.


    * **nbs_source** – Source of the neighbors.



#### SHORT_NAME(_ = 'DeltaDistanceNbSetWeight_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dd)
Initialize from dict.


* **Parameters**

    **dd** – Dict representation of DeltaDistanceNbSetWeight.



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



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight(weight_type='has_intersection', surface_definition={'angle_bounds': {'lower': 0.1, 'upper': 0.8}, 'distance_bounds': {'lower': 1.2, 'upper': 1.8}, 'type': 'standard_elliptic'}, nb_sets_from_hints='fallback_to_source', other_nb_sets='0_weight', additional_condition=1, smoothstep_distance=None, smoothstep_angle=None)
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

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dd)
Initialize from dict.


* **Parameters**

    **dd** – Dict representation of DistanceAngleAreaNbSetWeight.



* **Returns**

    DistanceAngleAreaNbSetWeight.



#### rectangle_crosses_area(d1, d2, a1, a2)
Whether a given rectangle crosses the area defined by the upper and lower curves.


* **Parameters**


    * **d1** – lower d.


    * **d2** – upper d.


    * **a1** – lower a.


    * **a2** – upper a.



* **Returns**




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



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat(myfloat)
Bases: `float`, `StrategyOption`

Distance cutoff in a strategy.

Special float that should be between 1 and infinity.


* **Parameters**

    **myfloat** – Distance cutoff.



#### allowed_values(_: str | Non_ _ = 'Real number between 1 and +infinity_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(d)
Initialize distance cutoff from dict.


* **Parameters**

    **d** – Dict representation of the distance cutoff.



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight(weight_function=None, nbs_source='voronoi')
Bases: `NbSetWeight`

Weight of neighbors set based on the distance.

Initialize DistanceNbSetWeight.


* **Parameters**


    * **weight_function** – Ratio function to use.


    * **nbs_source** – Source of the neighbors.



#### SHORT_NAME(_ = 'DistanceNbSetWeight_ )

#### as_dict()
MSOnable dict.


#### _classmethod_ from_dict(dd)
Initialize from dict.


* **Parameters**

    **dd** – Dict representation of DistanceNbSetWeight.



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



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight(distance_function=None, weight_function=None)
Bases: `NbSetWeight`

Weight of neighbors set based on the distance.

Initialize DistancePlateauNbSetWeight.


* **Parameters**


    * **distance_function** – Distance function to use.


    * **weight_function** – Ratio function to use.



#### SHORT_NAME(_ = 'DistancePlateauWeight_ )

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dd)
Initialize from dict.


* **Parameters**

    **dd** – Dict representation of DistancePlateauNbSetWeight.



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



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy(structure_environments=None, additional_condition=1, symmetry_measure_type='csm_wcs_ctwcc', dist_ang_area_weight=None, self_csm_weight=None, delta_csm_weight=None, cn_bias_weight=None, angle_weight=None, normalized_angle_distance_weight=None, ce_estimator={'function': 'power2_inverse_power2_decreasing', 'options': {'max_csm': 8.0}})
Bases: `WeightedNbSetChemenvStrategy`

MultiWeightsChemenvStrategy.

Constructor for the MultiWeightsChemenvStrategy.
:param structure_environments: StructureEnvironments object containing all the information on the

> coordination of the sites in a structure.


#### DEFAULT_CE_ESTIMATOR(_ = {'function': 'power2_inverse_power2_decreasing', 'options': {'max_csm': 8.0}_ )

#### STRATEGY_DESCRIPTION(_: str | Non_ _ = '    Multi Weights ChemenvStrategy_ )

#### as_dict()

* **Returns**

    Bson-serializable dict representation of the MultiWeightsChemenvStrategy object.



#### _classmethod_ from_dict(d)
Reconstructs the MultiWeightsChemenvStrategy object from a dict representation of the
MultipleAbundanceChemenvStrategy object created using the as_dict method.
:param d: dict representation of the MultiWeightsChemenvStrategy object
:return: MultiWeightsChemenvStrategy object.


#### _classmethod_ stats_article_weights_parameters()
Initialize strategy used in the statistics article.


#### _property_ uniquely_determines_coordination_environments()
Whether this strategy uniquely determines coordination environments.


### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NbSetWeight()
Bases: `MSONable`

Abstract object for neighbors sets weights estimations.


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



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight(average_type, aa, bb)
Bases: `NbSetWeight`

Weight of neighbors set based on the normalized angle/distance.

Initialize NormalizedAngleDistanceNbSetWeight.


* **Parameters**


    * **average_type** – Average function.


    * **aa** – Exponent for the angle values.


    * **bb** – Exponent for the distance values.



#### SHORT_NAME(_ = 'NormAngleDistWeight_ )

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



#### _classmethod_ from_dict(dd)
Initialize from dict.


* **Parameters**

    **dd** – Dict representation of NormalizedAngleDistanceNbSetWeight.



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



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight(effective_csm_estimator={'function': 'power2_inverse_decreasing', 'options': {'max_csm': 8.0}}, weight_estimator={'function': 'power2_decreasing_exp', 'options': {'alpha': 1, 'max_csm': 8.0}}, symmetry_measure_type='csm_wcs_ctwcc')
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

#### as_dict()
MSONable dict.


#### _classmethod_ from_dict(dd)
Initialize from dict.


* **Parameters**

    **dd** – Dict representation of SelfCSMNbSetWeight.



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



### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy(structure_environments=None, additional_condition=1, symmetry_measure_type='csm_wcs_ctwcc')
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

#### as_dict()
Bson-serializable dict representation of the SimpleAbundanceChemenvStrategy object.
:return: Bson-serializable dict representation of the SimpleAbundanceChemenvStrategy object.


#### _classmethod_ from_dict(d)
Reconstructs the SimpleAbundanceChemenvStrategy object from a dict representation of the
SimpleAbundanceChemenvStrategy object created using the as_dict method.
:param d: dict representation of the SimpleAbundanceChemenvStrategy object
:return: StructureEnvironments object.


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


### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy(structure_environments=None, distance_cutoff=1.4, angle_cutoff=0.3, additional_condition=1, continuous_symmetry_measure_cutoff=10, symmetry_measure_type='csm_wcs_ctwcc')
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

#### add_strategy_visualization_to_subplot(subplot, visualization_options=None, plot_type=None)
Add a visual of the strategy on a distance-angle plot.


* **Parameters**


    * **subplot** – Axes object onto the visual should be added.


    * **visualization_options** – Options for the visual.


    * **plot_type** – Type of distance-angle plot.



* **Returns**

    None



#### _property_ additional_condition()
Additional condition for this strategy.


#### _property_ angle_cutoff()
Angle cutoff used.


#### as_dict()
Bson-serializable dict representation of the SimplestChemenvStrategy object.
:return: Bson-serializable dict representation of the SimplestChemenvStrategy object.


#### _property_ continuous_symmetry_measure_cutoff()
CSM cutoff used.


#### _property_ distance_cutoff()
Distance cutoff used.


#### _classmethod_ from_dict(d)
Reconstructs the SimplestChemenvStrategy object from a dict representation of the SimplestChemenvStrategy object
created using the as_dict method.
:param d: dict representation of the SimplestChemenvStrategy object
:return: StructureEnvironments object.


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


### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.StrategyOption()
Bases: `MSONable`

Abstract class for the options of the chemenv strategies.


#### allowed_values(_: str | Non_ _ = Non_ )

#### _abstract_ as_dict()
A JSON-serializable dict representation of this strategy option.


### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy(structure_environments=None, truncate_dist_ang=True, additional_condition=1, max_nabundant=5, target_environments=('O:6',), target_penalty_type='max_csm', max_csm=5.0, symmetry_measure_type='csm_wcs_ctwcc')
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

#### as_dict()
Bson-serializable dict representation of the TargettedPenaltiedAbundanceChemenvStrategy object.
:return: Bson-serializable dict representation of the TargettedPenaltiedAbundanceChemenvStrategy object.


#### _classmethod_ from_dict(d)
Reconstructs the TargettedPenaltiedAbundanceChemenvStrategy object from a dict representation of the
TargettedPenaltiedAbundanceChemenvStrategy object created using the as_dict method.
:param d: dict representation of the TargettedPenaltiedAbundanceChemenvStrategy object
:return: TargettedPenaltiedAbundanceChemenvStrategy object.


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


### _class_ pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy(structure_environments=None, additional_condition=1, symmetry_measure_type='csm_wcs_ctwcc', nb_set_weights=None, ce_estimator={'function': 'power2_inverse_power2_decreasing', 'options': {'max_csm': 8.0}})
Bases: `AbstractChemenvStrategy`

WeightedNbSetChemenvStrategy.

Constructor for the WeightedNbSetChemenvStrategy.
:param structure_environments: StructureEnvironments object containing all the information on the

> coordination of the sites in a structure.


#### DEFAULT_CE_ESTIMATOR(_ = {'function': 'power2_inverse_power2_decreasing', 'options': {'max_csm': 8.0}_ )

#### STRATEGY_DESCRIPTION(_: str | Non_ _ = '    WeightedNbSetChemenvStrategy_ )

#### as_dict()
Bson-serializable dict representation of the WeightedNbSetChemenvStrategy object.
:return: Bson-serializable dict representation of the WeightedNbSetChemenvStrategy object.


#### _classmethod_ from_dict(d)
Reconstructs the WeightedNbSetChemenvStrategy object from a dict representation of the
WeightedNbSetChemenvStrategy object created using the as_dict method.
:param d: dict representation of the WeightedNbSetChemenvStrategy object
:return: WeightedNbSetChemenvStrategy object.


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


### pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.get_effective_csm(nb_set, cn_map, structure_environments, additional_info, symmetry_measure_type, max_effective_csm, effective_csm_estimator_ratio_function)
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



### pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.set_info(additional_info, field, isite, cn_map, value)
Set additional information for the weights.


* **Parameters**


    * **additional_info** – Additional information.


    * **field** – Type of additional information.


    * **isite** – Index of site to add info.


    * **cn_map** – Mapping index of the neighbors set.


    * **value** – Value of this additional information.



* **Returns**

    None