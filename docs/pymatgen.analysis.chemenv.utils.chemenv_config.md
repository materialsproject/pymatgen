---
layout: default
title: pymatgen.analysis.chemenv.utils.chemenv_config.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.utils.chemenv_config module

This module contains the classes for configuration of the chemenv package.


### _class_ pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig(package_options=None)
Bases: `object`

Class used to store the configuration of the chemenv package :


    * Materials project access


    * ICSD database access


    * Default options (strategies, …).


* **Parameters**

    **package_options** –



#### DEFAULT_PACKAGE_OPTIONS(_ = {'default_max_distance_factor': 1.5, 'default_strategy': {'strategy': 'SimplestChemenvStrategy', 'strategy_options': {'additional_condition': 1, 'angle_cutoff': 0.3, 'continuous_symmetry_measure_cutoff': 10, 'distance_cutoff': 1.4}}_ )

#### _classmethod_ auto_load(root_dir=None)
Autoload options.
:param root_dir:


#### _property_ has_materials_project_access()
Whether MP access is enabled.
:return:


#### package_options_description()
Describe package options.


#### save(root_dir=None)
Save the options.
:param root_dir:


#### setup()
Setup the class.


#### setup_package_options()
Setup the package options.