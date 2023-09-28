---
layout: default
title: pymatgen.analysis.chemenv.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis.chemenv package

Package for analyzing chemical environments.

## Subpackages


* [pymatgen.analysis.chemenv.connectivity package](pymatgen.analysis.chemenv.connectivity.md)




    * [pymatgen.analysis.chemenv.connectivity.connected_components module](pymatgen.analysis.chemenv.connectivity.md#module-pymatgen.analysis.chemenv.connectivity.connected_components)


        * [`ConnectedComponent`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent)


            * [`ConnectedComponent._edgedictkey_to_edgekey()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent._edgedictkey_to_edgekey)


            * [`ConnectedComponent._edgekey_to_edgedictkey()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent._edgekey_to_edgedictkey)


            * [`ConnectedComponent._order_periodicity_vectors()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent._order_periodicity_vectors)


            * [`ConnectedComponent._order_vectors()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent._order_vectors)


            * [`ConnectedComponent._retuplify_edgedata()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent._retuplify_edgedata)


            * [`ConnectedComponent.as_dict()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.as_dict)


            * [`ConnectedComponent.compute_periodicity()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.compute_periodicity)


            * [`ConnectedComponent.compute_periodicity_all_simple_paths_algorithm()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.compute_periodicity_all_simple_paths_algorithm)


            * [`ConnectedComponent.compute_periodicity_cycle_basis()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.compute_periodicity_cycle_basis)


            * [`ConnectedComponent.coordination_sequence()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.coordination_sequence)


            * [`ConnectedComponent.description()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.description)


            * [`ConnectedComponent.elastic_centered_graph()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.elastic_centered_graph)


            * [`ConnectedComponent.from_dict()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.from_dict)


            * [`ConnectedComponent.from_graph()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.from_graph)


            * [`ConnectedComponent.graph`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.graph)


            * [`ConnectedComponent.is_0d`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.is_0d)


            * [`ConnectedComponent.is_1d`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.is_1d)


            * [`ConnectedComponent.is_2d`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.is_2d)


            * [`ConnectedComponent.is_3d`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.is_3d)


            * [`ConnectedComponent.is_periodic`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.is_periodic)


            * [`ConnectedComponent.make_supergraph()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.make_supergraph)


            * [`ConnectedComponent.periodicity`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.periodicity)


            * [`ConnectedComponent.periodicity_vectors`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.periodicity_vectors)


            * [`ConnectedComponent.show_graph()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.show_graph)


        * [`draw_network()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.draw_network)


        * [`make_supergraph()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connected_components.make_supergraph)


    * [pymatgen.analysis.chemenv.connectivity.connectivity_finder module](pymatgen.analysis.chemenv.connectivity.md#module-pymatgen.analysis.chemenv.connectivity.connectivity_finder)


        * [`ConnectivityFinder`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connectivity_finder.ConnectivityFinder)


            * [`ConnectivityFinder.get_structure_connectivity()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connectivity_finder.ConnectivityFinder.get_structure_connectivity)


            * [`ConnectivityFinder.setup_parameters()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.connectivity_finder.ConnectivityFinder.setup_parameters)


    * [pymatgen.analysis.chemenv.connectivity.environment_nodes module](pymatgen.analysis.chemenv.connectivity.md#module-pymatgen.analysis.chemenv.connectivity.environment_nodes)


        * [`AbstractEnvironmentNode`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode)


            * [`AbstractEnvironmentNode.ATOM`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.ATOM)


            * [`AbstractEnvironmentNode.CE_NNBCES_NBCES_LIGANDS`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.CE_NNBCES_NBCES_LIGANDS)


            * [`AbstractEnvironmentNode.COORDINATION_ENVIRONMENT`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.COORDINATION_ENVIRONMENT)


            * [`AbstractEnvironmentNode.DEFAULT_EXTENSIONS`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.DEFAULT_EXTENSIONS)


            * [`AbstractEnvironmentNode.LIGANDS_ARRANGEMENT`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.LIGANDS_ARRANGEMENT)


            * [`AbstractEnvironmentNode.NEIGHBORING_CES`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NEIGHBORING_CES)


            * [`AbstractEnvironmentNode.NEIGHBORING_COORDINATION_ENVIRONMENTS`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NEIGHBORING_COORDINATION_ENVIRONMENTS)


            * [`AbstractEnvironmentNode.NEIGHBORS_LIGANDS_ARRANGEMENT`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NEIGHBORS_LIGANDS_ARRANGEMENT)


            * [`AbstractEnvironmentNode.NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_CE`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_CE)


            * [`AbstractEnvironmentNode.NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_COORDINATION_ENVIRONMENT`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_COORDINATION_ENVIRONMENT)


            * [`AbstractEnvironmentNode.NUMBER_OF_NEIGHBORING_CES`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NUMBER_OF_NEIGHBORING_CES)


            * [`AbstractEnvironmentNode.NUMBER_OF_NEIGHBORING_COORDINATION_ENVIRONMENTS`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NUMBER_OF_NEIGHBORING_COORDINATION_ENVIRONMENTS)


            * [`AbstractEnvironmentNode.atom_symbol`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.atom_symbol)


            * [`AbstractEnvironmentNode.ce`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.ce)


            * [`AbstractEnvironmentNode.ce_symbol`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.ce_symbol)


            * [`AbstractEnvironmentNode.coordination_environment`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.coordination_environment)


            * [`AbstractEnvironmentNode.everything_equal()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.everything_equal)


            * [`AbstractEnvironmentNode.isite`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.isite)


            * [`AbstractEnvironmentNode.mp_symbol`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.mp_symbol)


        * [`EnvironmentNode`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.EnvironmentNode)


            * [`EnvironmentNode.coordination_environment`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.EnvironmentNode.coordination_environment)


            * [`EnvironmentNode.everything_equal()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.EnvironmentNode.everything_equal)


        * [`get_environment_node()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.get_environment_node)


    * [pymatgen.analysis.chemenv.connectivity.structure_connectivity module](pymatgen.analysis.chemenv.connectivity.md#module-pymatgen.analysis.chemenv.connectivity.structure_connectivity)


        * [`StructureConnectivity`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity)


            * [`StructureConnectivity.add_bonds()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.add_bonds)


            * [`StructureConnectivity.add_sites()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.add_sites)


            * [`StructureConnectivity.as_dict()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.as_dict)


            * [`StructureConnectivity.environment_subgraph()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.environment_subgraph)


            * [`StructureConnectivity.from_dict()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.from_dict)


            * [`StructureConnectivity.get_connected_components()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.get_connected_components)


            * [`StructureConnectivity.print_links()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.print_links)


            * [`StructureConnectivity.setup_atom_environment_subgraph()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.setup_atom_environment_subgraph)


            * [`StructureConnectivity.setup_atom_environments_subgraph()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.setup_atom_environments_subgraph)


            * [`StructureConnectivity.setup_connectivity_description()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.setup_connectivity_description)


            * [`StructureConnectivity.setup_environment_subgraph()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.setup_environment_subgraph)


            * [`StructureConnectivity.setup_environments_subgraph()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.setup_environments_subgraph)


        * [`get_delta_image()`](pymatgen.analysis.chemenv.connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.get_delta_image)


* [pymatgen.analysis.chemenv.coordination_environments package](pymatgen.analysis.chemenv.coordination_environments.md)


    * [Subpackages](pymatgen.analysis.chemenv.coordination_environments.md#subpackages)


        * [pymatgen.analysis.chemenv.coordination_environments.coordination_geometries_files package](pymatgen.analysis.chemenv.coordination_environments.coordination_geometries_files.md)




    * [pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies module](pymatgen.analysis.chemenv.coordination_environments.md#module-pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies)


        * [`AbstractChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy)


            * [`AbstractChemenvStrategy.AC`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.AC)


            * [`AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.DEFAULT_SYMMETRY_MEASURE_TYPE)


            * [`AbstractChemenvStrategy.STRATEGY_DESCRIPTION`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.STRATEGY_DESCRIPTION)


            * [`AbstractChemenvStrategy.STRATEGY_INFO_FIELDS`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.STRATEGY_INFO_FIELDS)


            * [`AbstractChemenvStrategy.STRATEGY_OPTIONS`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.STRATEGY_OPTIONS)


            * [`AbstractChemenvStrategy._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy._abc_impl)


            * [`AbstractChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.as_dict)


            * [`AbstractChemenvStrategy.equivalent_site_index_and_transform()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.equivalent_site_index_and_transform)


            * [`AbstractChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.from_dict)


            * [`AbstractChemenvStrategy.get_site_ce_fractions_and_neighbors()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.get_site_ce_fractions_and_neighbors)


            * [`AbstractChemenvStrategy.get_site_coordination_environment()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.get_site_coordination_environment)


            * [`AbstractChemenvStrategy.get_site_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.get_site_coordination_environments)


            * [`AbstractChemenvStrategy.get_site_coordination_environments_fractions()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.get_site_coordination_environments_fractions)


            * [`AbstractChemenvStrategy.get_site_neighbors()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.get_site_neighbors)


            * [`AbstractChemenvStrategy.prepare_symmetries()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.prepare_symmetries)


            * [`AbstractChemenvStrategy.set_option()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.set_option)


            * [`AbstractChemenvStrategy.set_structure_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.set_structure_environments)


            * [`AbstractChemenvStrategy.setup_options()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.setup_options)


            * [`AbstractChemenvStrategy.symmetry_measure_type`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.symmetry_measure_type)


            * [`AbstractChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AbstractChemenvStrategy.uniquely_determines_coordination_environments)


        * [`AdditionalConditionInt`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt)


            * [`AdditionalConditionInt._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt._abc_impl)


            * [`AdditionalConditionInt.allowed_values`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt.allowed_values)


            * [`AdditionalConditionInt.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt.as_dict)


            * [`AdditionalConditionInt.description`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt.description)


            * [`AdditionalConditionInt.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt.from_dict)


            * [`AdditionalConditionInt.integer`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AdditionalConditionInt.integer)


        * [`AngleCutoffFloat`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat)


            * [`AngleCutoffFloat._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat._abc_impl)


            * [`AngleCutoffFloat.allowed_values`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat.allowed_values)


            * [`AngleCutoffFloat.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat.as_dict)


            * [`AngleCutoffFloat.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleCutoffFloat.from_dict)


        * [`AngleNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight)


            * [`AngleNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.SHORT_NAME)


            * [`AngleNbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight._abc_impl)


            * [`AngleNbSetWeight.angle_sum()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.angle_sum)


            * [`AngleNbSetWeight.angle_sumn()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.angle_sumn)


            * [`AngleNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.as_dict)


            * [`AngleNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.from_dict)


            * [`AngleNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AngleNbSetWeight.weight)


        * [`AnglePlateauNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight)


            * [`AnglePlateauNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight.SHORT_NAME)


            * [`AnglePlateauNbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight._abc_impl)


            * [`AnglePlateauNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight.as_dict)


            * [`AnglePlateauNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight.from_dict)


            * [`AnglePlateauNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.AnglePlateauNbSetWeight.weight)


        * [`CNBiasNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight)


            * [`CNBiasNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.SHORT_NAME)


            * [`CNBiasNbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight._abc_impl)


            * [`CNBiasNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.as_dict)


            * [`CNBiasNbSetWeight.explicit()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.explicit)


            * [`CNBiasNbSetWeight.from_description()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.from_description)


            * [`CNBiasNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.from_dict)


            * [`CNBiasNbSetWeight.geometrically_equidistant()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.geometrically_equidistant)


            * [`CNBiasNbSetWeight.linearly_equidistant()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.linearly_equidistant)


            * [`CNBiasNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CNBiasNbSetWeight.weight)


        * [`CSMFloat`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat)


            * [`CSMFloat._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat._abc_impl)


            * [`CSMFloat.allowed_values`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat.allowed_values)


            * [`CSMFloat.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat.as_dict)


            * [`CSMFloat.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.CSMFloat.from_dict)


        * [`DeltaCSMNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight)


            * [`DeltaCSMNbSetWeight.DEFAULT_EFFECTIVE_CSM_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.DEFAULT_EFFECTIVE_CSM_ESTIMATOR)


            * [`DeltaCSMNbSetWeight.DEFAULT_SYMMETRY_MEASURE_TYPE`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.DEFAULT_SYMMETRY_MEASURE_TYPE)


            * [`DeltaCSMNbSetWeight.DEFAULT_WEIGHT_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.DEFAULT_WEIGHT_ESTIMATOR)


            * [`DeltaCSMNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.SHORT_NAME)


            * [`DeltaCSMNbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight._abc_impl)


            * [`DeltaCSMNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.as_dict)


            * [`DeltaCSMNbSetWeight.delta_cn_specifics()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.delta_cn_specifics)


            * [`DeltaCSMNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.from_dict)


            * [`DeltaCSMNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaCSMNbSetWeight.weight)


        * [`DeltaDistanceNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight)


            * [`DeltaDistanceNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight.SHORT_NAME)


            * [`DeltaDistanceNbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight._abc_impl)


            * [`DeltaDistanceNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight.as_dict)


            * [`DeltaDistanceNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight.from_dict)


            * [`DeltaDistanceNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DeltaDistanceNbSetWeight.weight)


        * [`DistanceAngleAreaNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight)


            * [`DistanceAngleAreaNbSetWeight.AC`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.AC)


            * [`DistanceAngleAreaNbSetWeight.DEFAULT_SURFACE_DEFINITION`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.DEFAULT_SURFACE_DEFINITION)


            * [`DistanceAngleAreaNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.SHORT_NAME)


            * [`DistanceAngleAreaNbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight._abc_impl)


            * [`DistanceAngleAreaNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.as_dict)


            * [`DistanceAngleAreaNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.from_dict)


            * [`DistanceAngleAreaNbSetWeight.rectangle_crosses_area()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.rectangle_crosses_area)


            * [`DistanceAngleAreaNbSetWeight.w_area_has_intersection()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.w_area_has_intersection)


            * [`DistanceAngleAreaNbSetWeight.w_area_intersection_nbsfh_fbs_onb0()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.w_area_intersection_nbsfh_fbs_onb0)


            * [`DistanceAngleAreaNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceAngleAreaNbSetWeight.weight)


        * [`DistanceCutoffFloat`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat)


            * [`DistanceCutoffFloat._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat._abc_impl)


            * [`DistanceCutoffFloat.allowed_values`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat.allowed_values)


            * [`DistanceCutoffFloat.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat.as_dict)


            * [`DistanceCutoffFloat.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceCutoffFloat.from_dict)


        * [`DistanceNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight)


            * [`DistanceNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight.SHORT_NAME)


            * [`DistanceNbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight._abc_impl)


            * [`DistanceNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight.as_dict)


            * [`DistanceNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight.from_dict)


            * [`DistanceNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistanceNbSetWeight.weight)


        * [`DistancePlateauNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight)


            * [`DistancePlateauNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight.SHORT_NAME)


            * [`DistancePlateauNbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight._abc_impl)


            * [`DistancePlateauNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight.as_dict)


            * [`DistancePlateauNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight.from_dict)


            * [`DistancePlateauNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.DistancePlateauNbSetWeight.weight)


        * [`MultiWeightsChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy)


            * [`MultiWeightsChemenvStrategy.DEFAULT_CE_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.DEFAULT_CE_ESTIMATOR)


            * [`MultiWeightsChemenvStrategy.STRATEGY_DESCRIPTION`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.STRATEGY_DESCRIPTION)


            * [`MultiWeightsChemenvStrategy._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy._abc_impl)


            * [`MultiWeightsChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.as_dict)


            * [`MultiWeightsChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.from_dict)


            * [`MultiWeightsChemenvStrategy.stats_article_weights_parameters()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.stats_article_weights_parameters)


            * [`MultiWeightsChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.MultiWeightsChemenvStrategy.uniquely_determines_coordination_environments)


        * [`NbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NbSetWeight)


            * [`NbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NbSetWeight._abc_impl)


            * [`NbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NbSetWeight.as_dict)


            * [`NbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NbSetWeight.weight)


        * [`NormalizedAngleDistanceNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight)


            * [`NormalizedAngleDistanceNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.SHORT_NAME)


            * [`NormalizedAngleDistanceNbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight._abc_impl)


            * [`NormalizedAngleDistanceNbSetWeight.ang()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.ang)


            * [`NormalizedAngleDistanceNbSetWeight.anginvdist()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.anginvdist)


            * [`NormalizedAngleDistanceNbSetWeight.anginvndist()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.anginvndist)


            * [`NormalizedAngleDistanceNbSetWeight.angn()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.angn)


            * [`NormalizedAngleDistanceNbSetWeight.angninvdist()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.angninvdist)


            * [`NormalizedAngleDistanceNbSetWeight.angninvndist()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.angninvndist)


            * [`NormalizedAngleDistanceNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.as_dict)


            * [`NormalizedAngleDistanceNbSetWeight.aweight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.aweight)


            * [`NormalizedAngleDistanceNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.from_dict)


            * [`NormalizedAngleDistanceNbSetWeight.gweight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.gweight)


            * [`NormalizedAngleDistanceNbSetWeight.invdist()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.invdist)


            * [`NormalizedAngleDistanceNbSetWeight.invndist()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.invndist)


            * [`NormalizedAngleDistanceNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.NormalizedAngleDistanceNbSetWeight.weight)


        * [`SelfCSMNbSetWeight`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight)


            * [`SelfCSMNbSetWeight.DEFAULT_EFFECTIVE_CSM_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.DEFAULT_EFFECTIVE_CSM_ESTIMATOR)


            * [`SelfCSMNbSetWeight.DEFAULT_SYMMETRY_MEASURE_TYPE`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.DEFAULT_SYMMETRY_MEASURE_TYPE)


            * [`SelfCSMNbSetWeight.DEFAULT_WEIGHT_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.DEFAULT_WEIGHT_ESTIMATOR)


            * [`SelfCSMNbSetWeight.SHORT_NAME`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.SHORT_NAME)


            * [`SelfCSMNbSetWeight._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight._abc_impl)


            * [`SelfCSMNbSetWeight.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.as_dict)


            * [`SelfCSMNbSetWeight.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.from_dict)


            * [`SelfCSMNbSetWeight.weight()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SelfCSMNbSetWeight.weight)


        * [`SimpleAbundanceChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy)


            * [`SimpleAbundanceChemenvStrategy.DEFAULT_ADDITIONAL_CONDITION`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.DEFAULT_ADDITIONAL_CONDITION)


            * [`SimpleAbundanceChemenvStrategy.DEFAULT_MAX_DIST`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.DEFAULT_MAX_DIST)


            * [`SimpleAbundanceChemenvStrategy.STRATEGY_DESCRIPTION`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.STRATEGY_DESCRIPTION)


            * [`SimpleAbundanceChemenvStrategy.STRATEGY_OPTIONS`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.STRATEGY_OPTIONS)


            * [`SimpleAbundanceChemenvStrategy._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy._abc_impl)


            * [`SimpleAbundanceChemenvStrategy._get_map()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy._get_map)


            * [`SimpleAbundanceChemenvStrategy._get_maps_surfaces()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy._get_maps_surfaces)


            * [`SimpleAbundanceChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.as_dict)


            * [`SimpleAbundanceChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.from_dict)


            * [`SimpleAbundanceChemenvStrategy.get_site_coordination_environment()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.get_site_coordination_environment)


            * [`SimpleAbundanceChemenvStrategy.get_site_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.get_site_coordination_environments)


            * [`SimpleAbundanceChemenvStrategy.get_site_neighbors()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.get_site_neighbors)


            * [`SimpleAbundanceChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimpleAbundanceChemenvStrategy.uniquely_determines_coordination_environments)


        * [`SimplestChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy)


            * [`SimplestChemenvStrategy.DEFAULT_ADDITIONAL_CONDITION`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.DEFAULT_ADDITIONAL_CONDITION)


            * [`SimplestChemenvStrategy.DEFAULT_ANGLE_CUTOFF`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.DEFAULT_ANGLE_CUTOFF)


            * [`SimplestChemenvStrategy.DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.DEFAULT_CONTINUOUS_SYMMETRY_MEASURE_CUTOFF)


            * [`SimplestChemenvStrategy.DEFAULT_DISTANCE_CUTOFF`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.DEFAULT_DISTANCE_CUTOFF)


            * [`SimplestChemenvStrategy.STRATEGY_DESCRIPTION`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.STRATEGY_DESCRIPTION)


            * [`SimplestChemenvStrategy.STRATEGY_OPTIONS`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.STRATEGY_OPTIONS)


            * [`SimplestChemenvStrategy._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy._abc_impl)


            * [`SimplestChemenvStrategy.add_strategy_visualization_to_subplot()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.add_strategy_visualization_to_subplot)


            * [`SimplestChemenvStrategy.additional_condition`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.additional_condition)


            * [`SimplestChemenvStrategy.angle_cutoff`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.angle_cutoff)


            * [`SimplestChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.as_dict)


            * [`SimplestChemenvStrategy.continuous_symmetry_measure_cutoff`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.continuous_symmetry_measure_cutoff)


            * [`SimplestChemenvStrategy.distance_cutoff`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.distance_cutoff)


            * [`SimplestChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.from_dict)


            * [`SimplestChemenvStrategy.get_site_coordination_environment()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.get_site_coordination_environment)


            * [`SimplestChemenvStrategy.get_site_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.get_site_coordination_environments)


            * [`SimplestChemenvStrategy.get_site_coordination_environments_fractions()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.get_site_coordination_environments_fractions)


            * [`SimplestChemenvStrategy.get_site_neighbors()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.get_site_neighbors)


            * [`SimplestChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.SimplestChemenvStrategy.uniquely_determines_coordination_environments)


        * [`StrategyOption`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.StrategyOption)


            * [`StrategyOption._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.StrategyOption._abc_impl)


            * [`StrategyOption.allowed_values`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.StrategyOption.allowed_values)


            * [`StrategyOption.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.StrategyOption.as_dict)


        * [`TargettedPenaltiedAbundanceChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy)


            * [`TargettedPenaltiedAbundanceChemenvStrategy.DEFAULT_TARGET_ENVIRONMENTS`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy.DEFAULT_TARGET_ENVIRONMENTS)


            * [`TargettedPenaltiedAbundanceChemenvStrategy._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy._abc_impl)


            * [`TargettedPenaltiedAbundanceChemenvStrategy._get_map()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy._get_map)


            * [`TargettedPenaltiedAbundanceChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy.as_dict)


            * [`TargettedPenaltiedAbundanceChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy.from_dict)


            * [`TargettedPenaltiedAbundanceChemenvStrategy.get_site_coordination_environment()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy.get_site_coordination_environment)


            * [`TargettedPenaltiedAbundanceChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.TargettedPenaltiedAbundanceChemenvStrategy.uniquely_determines_coordination_environments)


        * [`WeightedNbSetChemenvStrategy`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy)


            * [`WeightedNbSetChemenvStrategy.DEFAULT_CE_ESTIMATOR`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.DEFAULT_CE_ESTIMATOR)


            * [`WeightedNbSetChemenvStrategy.STRATEGY_DESCRIPTION`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.STRATEGY_DESCRIPTION)


            * [`WeightedNbSetChemenvStrategy._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy._abc_impl)


            * [`WeightedNbSetChemenvStrategy.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.as_dict)


            * [`WeightedNbSetChemenvStrategy.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.from_dict)


            * [`WeightedNbSetChemenvStrategy.get_site_coordination_environment()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.get_site_coordination_environment)


            * [`WeightedNbSetChemenvStrategy.get_site_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.get_site_coordination_environments)


            * [`WeightedNbSetChemenvStrategy.get_site_coordination_environments_fractions()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.get_site_coordination_environments_fractions)


            * [`WeightedNbSetChemenvStrategy.get_site_neighbors()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.get_site_neighbors)


            * [`WeightedNbSetChemenvStrategy.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.WeightedNbSetChemenvStrategy.uniquely_determines_coordination_environments)


        * [`get_effective_csm()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.get_effective_csm)


        * [`set_info()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies.set_info)


    * [pymatgen.analysis.chemenv.coordination_environments.coordination_geometries module](pymatgen.analysis.chemenv.coordination_environments.md#module-pymatgen.analysis.chemenv.coordination_environments.coordination_geometries)


        * [`AbstractChemenvAlgorithm`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AbstractChemenvAlgorithm)


            * [`AbstractChemenvAlgorithm._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AbstractChemenvAlgorithm._abc_impl)


            * [`AbstractChemenvAlgorithm.algorithm_type`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AbstractChemenvAlgorithm.algorithm_type)


            * [`AbstractChemenvAlgorithm.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AbstractChemenvAlgorithm.as_dict)


        * [`AllCoordinationGeometries`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries)


            * [`AllCoordinationGeometries.get_geometries()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_geometries)


            * [`AllCoordinationGeometries.get_geometry_from_IUCr_symbol()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_geometry_from_IUCr_symbol)


            * [`AllCoordinationGeometries.get_geometry_from_IUPAC_symbol()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_geometry_from_IUPAC_symbol)


            * [`AllCoordinationGeometries.get_geometry_from_mp_symbol()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_geometry_from_mp_symbol)


            * [`AllCoordinationGeometries.get_geometry_from_name()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_geometry_from_name)


            * [`AllCoordinationGeometries.get_implemented_geometries()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_implemented_geometries)


            * [`AllCoordinationGeometries.get_not_implemented_geometries()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_not_implemented_geometries)


            * [`AllCoordinationGeometries.get_symbol_cn_mapping()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_symbol_cn_mapping)


            * [`AllCoordinationGeometries.get_symbol_name_mapping()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.get_symbol_name_mapping)


            * [`AllCoordinationGeometries.is_a_valid_coordination_geometry()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.is_a_valid_coordination_geometry)


            * [`AllCoordinationGeometries.pretty_print()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.AllCoordinationGeometries.pretty_print)


        * [`CoordinationGeometry`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry)


            * [`CoordinationGeometry.CSM_SKIP_SEPARATION_PLANE_ALGO`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.CSM_SKIP_SEPARATION_PLANE_ALGO)


            * [`CoordinationGeometry.IUCr_symbol`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.IUCr_symbol)


            * [`CoordinationGeometry.IUCr_symbol_str`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.IUCr_symbol_str)


            * [`CoordinationGeometry.IUPAC_symbol`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.IUPAC_symbol)


            * [`CoordinationGeometry.IUPAC_symbol_str`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.IUPAC_symbol_str)


            * [`CoordinationGeometry.NeighborsSetsHints`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints)


                * [`CoordinationGeometry.NeighborsSetsHints.ALLOWED_HINTS_TYPES`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.ALLOWED_HINTS_TYPES)


                * [`CoordinationGeometry.NeighborsSetsHints.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.as_dict)


                * [`CoordinationGeometry.NeighborsSetsHints.double_cap_hints()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.double_cap_hints)


                * [`CoordinationGeometry.NeighborsSetsHints.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.from_dict)


                * [`CoordinationGeometry.NeighborsSetsHints.hints()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.hints)


                * [`CoordinationGeometry.NeighborsSetsHints.single_cap_hints()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.single_cap_hints)


                * [`CoordinationGeometry.NeighborsSetsHints.triple_cap_hints()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.NeighborsSetsHints.triple_cap_hints)


            * [`CoordinationGeometry.algorithms`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.algorithms)


            * [`CoordinationGeometry.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.as_dict)


            * [`CoordinationGeometry.ce_symbol`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.ce_symbol)


            * [`CoordinationGeometry.coordination_number`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.coordination_number)


            * [`CoordinationGeometry.distfactor_max`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.distfactor_max)


            * [`CoordinationGeometry.edges()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.edges)


            * [`CoordinationGeometry.faces()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.faces)


            * [`CoordinationGeometry.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.from_dict)


            * [`CoordinationGeometry.get_central_site()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.get_central_site)


            * [`CoordinationGeometry.get_coordination_number()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.get_coordination_number)


            * [`CoordinationGeometry.get_name()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.get_name)


            * [`CoordinationGeometry.get_pmeshes()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.get_pmeshes)


            * [`CoordinationGeometry.is_implemented()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.is_implemented)


            * [`CoordinationGeometry.mp_symbol`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.mp_symbol)


            * [`CoordinationGeometry.number_of_permutations`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.number_of_permutations)


            * [`CoordinationGeometry.pauling_stability_ratio`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.pauling_stability_ratio)


            * [`CoordinationGeometry.ref_permutation()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.ref_permutation)


            * [`CoordinationGeometry.solid_angles()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.CoordinationGeometry.solid_angles)


        * [`ExplicitPermutationsAlgorithm`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.ExplicitPermutationsAlgorithm)


            * [`ExplicitPermutationsAlgorithm._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.ExplicitPermutationsAlgorithm._abc_impl)


            * [`ExplicitPermutationsAlgorithm.as_dict`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.ExplicitPermutationsAlgorithm.as_dict)


            * [`ExplicitPermutationsAlgorithm.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.ExplicitPermutationsAlgorithm.from_dict)


            * [`ExplicitPermutationsAlgorithm.permutations`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.ExplicitPermutationsAlgorithm.permutations)


        * [`SeparationPlane`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane)


            * [`SeparationPlane._abc_impl`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane._abc_impl)


            * [`SeparationPlane.argsorted_ref_separation_perm`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.argsorted_ref_separation_perm)


            * [`SeparationPlane.as_dict`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.as_dict)


            * [`SeparationPlane.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.from_dict)


            * [`SeparationPlane.permutations`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.permutations)


            * [`SeparationPlane.ref_separation_perm`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.ref_separation_perm)


            * [`SeparationPlane.safe_separation_permutations()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometries.SeparationPlane.safe_separation_permutations)


    * [pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder module](pymatgen.analysis.chemenv.coordination_environments.md#module-pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder)


        * [`AbstractGeometry`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry)


            * [`AbstractGeometry.cn`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.cn)


            * [`AbstractGeometry.coordination_number`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.coordination_number)


            * [`AbstractGeometry.from_cg()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.from_cg)


            * [`AbstractGeometry.points_wcs_csc()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wcs_csc)


            * [`AbstractGeometry.points_wcs_ctwcc()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wcs_ctwcc)


            * [`AbstractGeometry.points_wcs_ctwocc()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wcs_ctwocc)


            * [`AbstractGeometry.points_wocs_csc()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wocs_csc)


            * [`AbstractGeometry.points_wocs_ctwcc()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wocs_ctwcc)


            * [`AbstractGeometry.points_wocs_ctwocc()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.AbstractGeometry.points_wocs_ctwocc)


        * [`LocalGeometryFinder`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder)


            * [`LocalGeometryFinder.BVA_DISTANCE_SCALE_FACTORS`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.BVA_DISTANCE_SCALE_FACTORS)


            * [`LocalGeometryFinder.DEFAULT_BVA_DISTANCE_SCALE_FACTOR`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.DEFAULT_BVA_DISTANCE_SCALE_FACTOR)


            * [`LocalGeometryFinder.DEFAULT_SPG_ANALYZER_OPTIONS`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.DEFAULT_SPG_ANALYZER_OPTIONS)


            * [`LocalGeometryFinder.DEFAULT_STRATEGY`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.DEFAULT_STRATEGY)


            * [`LocalGeometryFinder.PRESETS`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.PRESETS)


            * [`LocalGeometryFinder.STRUCTURE_REFINEMENT_NONE`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.STRUCTURE_REFINEMENT_NONE)


            * [`LocalGeometryFinder.STRUCTURE_REFINEMENT_REFINED`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.STRUCTURE_REFINEMENT_REFINED)


            * [`LocalGeometryFinder.STRUCTURE_REFINEMENT_SYMMETRIZED`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.STRUCTURE_REFINEMENT_SYMMETRIZED)


            * [`LocalGeometryFinder._cg_csm_separation_plane()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder._cg_csm_separation_plane)


            * [`LocalGeometryFinder._cg_csm_separation_plane_optim1()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder._cg_csm_separation_plane_optim1)


            * [`LocalGeometryFinder._cg_csm_separation_plane_optim2()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder._cg_csm_separation_plane_optim2)


            * [`LocalGeometryFinder._update_results_all_csms()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder._update_results_all_csms)


            * [`LocalGeometryFinder.compute_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.compute_coordination_environments)


            * [`LocalGeometryFinder.compute_structure_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.compute_structure_environments)


            * [`LocalGeometryFinder.coordination_geometry_symmetry_measures()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures)


            * [`LocalGeometryFinder.coordination_geometry_symmetry_measures_fallback_random()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures_fallback_random)


            * [`LocalGeometryFinder.coordination_geometry_symmetry_measures_separation_plane()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures_separation_plane)


            * [`LocalGeometryFinder.coordination_geometry_symmetry_measures_separation_plane_optim()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures_separation_plane_optim)


            * [`LocalGeometryFinder.coordination_geometry_symmetry_measures_sepplane_optim()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures_sepplane_optim)


            * [`LocalGeometryFinder.coordination_geometry_symmetry_measures_standard()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.coordination_geometry_symmetry_measures_standard)


            * [`LocalGeometryFinder.get_coordination_symmetry_measures()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.get_coordination_symmetry_measures)


            * [`LocalGeometryFinder.get_coordination_symmetry_measures_optim()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.get_coordination_symmetry_measures_optim)


            * [`LocalGeometryFinder.get_structure()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.get_structure)


            * [`LocalGeometryFinder.set_structure()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.set_structure)


            * [`LocalGeometryFinder.setup_explicit_indices_local_geometry()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_explicit_indices_local_geometry)


            * [`LocalGeometryFinder.setup_local_geometry()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_local_geometry)


            * [`LocalGeometryFinder.setup_ordered_indices_local_geometry()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_ordered_indices_local_geometry)


            * [`LocalGeometryFinder.setup_parameter()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_parameter)


            * [`LocalGeometryFinder.setup_parameters()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_parameters)


            * [`LocalGeometryFinder.setup_random_indices_local_geometry()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_random_indices_local_geometry)


            * [`LocalGeometryFinder.setup_random_structure()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_random_structure)


            * [`LocalGeometryFinder.setup_structure()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_structure)


            * [`LocalGeometryFinder.setup_test_perfect_environment()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.setup_test_perfect_environment)


            * [`LocalGeometryFinder.update_nb_set_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.LocalGeometryFinder.update_nb_set_environments)


        * [`find_rotation()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.find_rotation)


        * [`find_scaling_factor()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.find_scaling_factor)


        * [`symmetry_measure()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder.symmetry_measure)


    * [pymatgen.analysis.chemenv.coordination_environments.structure_environments module](pymatgen.analysis.chemenv.coordination_environments.md#module-pymatgen.analysis.chemenv.coordination_environments.structure_environments)


        * [`ChemicalEnvironments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments)


            * [`ChemicalEnvironments.add_coord_geom()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.add_coord_geom)


            * [`ChemicalEnvironments.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.as_dict)


            * [`ChemicalEnvironments.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.from_dict)


            * [`ChemicalEnvironments.is_close_to()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.is_close_to)


            * [`ChemicalEnvironments.minimum_geometries()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.minimum_geometries)


            * [`ChemicalEnvironments.minimum_geometry()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.ChemicalEnvironments.minimum_geometry)


        * [`LightStructureEnvironments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments)


            * [`LightStructureEnvironments.DEFAULT_STATISTICS_FIELDS`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.DEFAULT_STATISTICS_FIELDS)


            * [`LightStructureEnvironments.DELTA_MAX_OXIDATION_STATE`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.DELTA_MAX_OXIDATION_STATE)


            * [`LightStructureEnvironments.NeighborsSet`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet)


                * [`LightStructureEnvironments.NeighborsSet.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.as_dict)


                * [`LightStructureEnvironments.NeighborsSet.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.from_dict)


                * [`LightStructureEnvironments.NeighborsSet.neighb_coords`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.neighb_coords)


                * [`LightStructureEnvironments.NeighborsSet.neighb_indices_and_images`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.neighb_indices_and_images)


                * [`LightStructureEnvironments.NeighborsSet.neighb_sites`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.neighb_sites)


                * [`LightStructureEnvironments.NeighborsSet.neighb_sites_and_indices`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.NeighborsSet.neighb_sites_and_indices)


            * [`LightStructureEnvironments.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.as_dict)


            * [`LightStructureEnvironments.clear_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.clear_environments)


            * [`LightStructureEnvironments.contains_only_one_anion()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.contains_only_one_anion)


            * [`LightStructureEnvironments.contains_only_one_anion_atom()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.contains_only_one_anion_atom)


            * [`LightStructureEnvironments.environments_identified()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.environments_identified)


            * [`LightStructureEnvironments.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.from_dict)


            * [`LightStructureEnvironments.from_structure_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.from_structure_environments)


            * [`LightStructureEnvironments.get_site_info_for_specie_allces()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.get_site_info_for_specie_allces)


            * [`LightStructureEnvironments.get_site_info_for_specie_ce()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.get_site_info_for_specie_ce)


            * [`LightStructureEnvironments.get_statistics()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.get_statistics)


            * [`LightStructureEnvironments.setup_statistic_lists()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.setup_statistic_lists)


            * [`LightStructureEnvironments.site_contains_environment()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.site_contains_environment)


            * [`LightStructureEnvironments.site_has_clear_environment()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.site_has_clear_environment)


            * [`LightStructureEnvironments.structure_contains_atom_environment()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.structure_contains_atom_environment)


            * [`LightStructureEnvironments.structure_has_clear_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.structure_has_clear_environments)


            * [`LightStructureEnvironments.uniquely_determines_coordination_environments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.LightStructureEnvironments.uniquely_determines_coordination_environments)


        * [`StructureEnvironments`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments)


            * [`StructureEnvironments.AC`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.AC)


            * [`StructureEnvironments.NeighborsSet`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet)


                * [`StructureEnvironments.NeighborsSet.add_source()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.add_source)


                * [`StructureEnvironments.NeighborsSet.angle_plateau()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.angle_plateau)


                * [`StructureEnvironments.NeighborsSet.angles`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.angles)


                * [`StructureEnvironments.NeighborsSet.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.as_dict)


                * [`StructureEnvironments.NeighborsSet.coords`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.coords)


                * [`StructureEnvironments.NeighborsSet.distance_plateau()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.distance_plateau)


                * [`StructureEnvironments.NeighborsSet.distances`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.distances)


                * [`StructureEnvironments.NeighborsSet.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.from_dict)


                * [`StructureEnvironments.NeighborsSet.get_neighb_voronoi_indices()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.get_neighb_voronoi_indices)


                * [`StructureEnvironments.NeighborsSet.info`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.info)


                * [`StructureEnvironments.NeighborsSet.neighb_coords`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.neighb_coords)


                * [`StructureEnvironments.NeighborsSet.neighb_coordsOpt`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.neighb_coordsOpt)


                * [`StructureEnvironments.NeighborsSet.neighb_sites`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.neighb_sites)


                * [`StructureEnvironments.NeighborsSet.neighb_sites_and_indices`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.neighb_sites_and_indices)


                * [`StructureEnvironments.NeighborsSet.normalized_angles`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.normalized_angles)


                * [`StructureEnvironments.NeighborsSet.normalized_distances`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.normalized_distances)


                * [`StructureEnvironments.NeighborsSet.source`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.source)


                * [`StructureEnvironments.NeighborsSet.voronoi_grid_surface_points()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.NeighborsSet.voronoi_grid_surface_points)


            * [`StructureEnvironments.add_neighbors_set()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.add_neighbors_set)


            * [`StructureEnvironments.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.as_dict)


            * [`StructureEnvironments.differences_wrt()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.differences_wrt)


            * [`StructureEnvironments.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.from_dict)


            * [`StructureEnvironments.get_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.get_coordination_environments)


            * [`StructureEnvironments.get_csm()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.get_csm)


            * [`StructureEnvironments.get_csm_and_maps()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.get_csm_and_maps)


            * [`StructureEnvironments.get_csms()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.get_csms)


            * [`StructureEnvironments.get_environments_figure()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.get_environments_figure)


            * [`StructureEnvironments.init_neighbors_sets()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.init_neighbors_sets)


            * [`StructureEnvironments.plot_csm_and_maps()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.plot_csm_and_maps)


            * [`StructureEnvironments.plot_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.plot_environments)


            * [`StructureEnvironments.save_environments_figure()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.save_environments_figure)


            * [`StructureEnvironments.update_coordination_environments()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.update_coordination_environments)


            * [`StructureEnvironments.update_site_info()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.structure_environments.StructureEnvironments.update_site_info)


    * [pymatgen.analysis.chemenv.coordination_environments.voronoi module](pymatgen.analysis.chemenv.coordination_environments.md#module-pymatgen.analysis.chemenv.coordination_environments.voronoi)


        * [`DetailedVoronoiContainer`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer)


            * [`DetailedVoronoiContainer.AC`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.AC)


            * [`DetailedVoronoiContainer._get_vertices_dist_ang_indices()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer._get_vertices_dist_ang_indices)


            * [`DetailedVoronoiContainer._precompute_additional_conditions()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer._precompute_additional_conditions)


            * [`DetailedVoronoiContainer._precompute_angle_conditions()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer._precompute_angle_conditions)


            * [`DetailedVoronoiContainer._precompute_distance_conditions()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer._precompute_distance_conditions)


            * [`DetailedVoronoiContainer.as_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.as_dict)


            * [`DetailedVoronoiContainer.default_normalized_angle_tolerance`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.default_normalized_angle_tolerance)


            * [`DetailedVoronoiContainer.default_normalized_distance_tolerance`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.default_normalized_distance_tolerance)


            * [`DetailedVoronoiContainer.default_voronoi_cutoff`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.default_voronoi_cutoff)


            * [`DetailedVoronoiContainer.from_dict()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.from_dict)


            * [`DetailedVoronoiContainer.get_rdf_figure()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.get_rdf_figure)


            * [`DetailedVoronoiContainer.get_sadf_figure()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.get_sadf_figure)


            * [`DetailedVoronoiContainer.is_close_to()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.is_close_to)


            * [`DetailedVoronoiContainer.maps_and_surfaces()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.maps_and_surfaces)


            * [`DetailedVoronoiContainer.maps_and_surfaces_bounded()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.maps_and_surfaces_bounded)


            * [`DetailedVoronoiContainer.neighbors()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.neighbors)


            * [`DetailedVoronoiContainer.neighbors_surfaces()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.neighbors_surfaces)


            * [`DetailedVoronoiContainer.neighbors_surfaces_bounded()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.neighbors_surfaces_bounded)


            * [`DetailedVoronoiContainer.setup_neighbors_distances_and_angles()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.setup_neighbors_distances_and_angles)


            * [`DetailedVoronoiContainer.setup_voronoi_list()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.setup_voronoi_list)


            * [`DetailedVoronoiContainer.to_bson_voronoi_list2()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.to_bson_voronoi_list2)


            * [`DetailedVoronoiContainer.voronoi_parameters_bounds_and_limits()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.DetailedVoronoiContainer.voronoi_parameters_bounds_and_limits)


        * [`from_bson_voronoi_list2()`](pymatgen.analysis.chemenv.coordination_environments.md#pymatgen.analysis.chemenv.coordination_environments.voronoi.from_bson_voronoi_list2)


* [pymatgen.analysis.chemenv.utils package](pymatgen.analysis.chemenv.utils.md)




    * [pymatgen.analysis.chemenv.utils.chemenv_config module](pymatgen.analysis.chemenv.utils.md#module-pymatgen.analysis.chemenv.utils.chemenv_config)


        * [`ChemEnvConfig`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig)


            * [`ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.DEFAULT_PACKAGE_OPTIONS)


            * [`ChemEnvConfig.auto_load()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.auto_load)


            * [`ChemEnvConfig.has_materials_project_access`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.has_materials_project_access)


            * [`ChemEnvConfig.package_options_description()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.package_options_description)


            * [`ChemEnvConfig.save()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.save)


            * [`ChemEnvConfig.setup()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.setup)


            * [`ChemEnvConfig.setup_package_options()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_config.ChemEnvConfig.setup_package_options)


    * [pymatgen.analysis.chemenv.utils.chemenv_errors module](pymatgen.analysis.chemenv.utils.md#module-pymatgen.analysis.chemenv.utils.chemenv_errors)


        * [`AbstractChemenvError`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_errors.AbstractChemenvError)


        * [`ChemenvError`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_errors.ChemenvError)


        * [`EquivalentSiteSearchError`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_errors.EquivalentSiteSearchError)


        * [`NeighborsNotComputedChemenvError`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_errors.NeighborsNotComputedChemenvError)


        * [`SolidAngleError`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.chemenv_errors.SolidAngleError)


    * [pymatgen.analysis.chemenv.utils.coordination_geometry_utils module](pymatgen.analysis.chemenv.utils.md#module-pymatgen.analysis.chemenv.utils.coordination_geometry_utils)


        * [`Plane`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane)


            * [`Plane.TEST_2D_POINTS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.TEST_2D_POINTS)


            * [`Plane.a`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.a)


            * [`Plane.abcd`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.abcd)


            * [`Plane.b`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.b)


            * [`Plane.c`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.c)


            * [`Plane.coefficients`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.coefficients)


            * [`Plane.crosses_origin`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.crosses_origin)


            * [`Plane.d`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.d)


            * [`Plane.distance_to_origin`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.distance_to_origin)


            * [`Plane.distance_to_point()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.distance_to_point)


            * [`Plane.distances()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.distances)


            * [`Plane.distances_indices_groups()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.distances_indices_groups)


            * [`Plane.distances_indices_sorted()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.distances_indices_sorted)


            * [`Plane.fit_error()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.fit_error)


            * [`Plane.fit_least_square_distance_error()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.fit_least_square_distance_error)


            * [`Plane.fit_maximum_distance_error()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.fit_maximum_distance_error)


            * [`Plane.from_2points_and_origin()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_2points_and_origin)


            * [`Plane.from_3points()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_3points)


            * [`Plane.from_coefficients()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_coefficients)


            * [`Plane.from_npoints()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_npoints)


            * [`Plane.from_npoints_least_square_distance()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_npoints_least_square_distance)


            * [`Plane.from_npoints_maximum_distance()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.from_npoints_maximum_distance)


            * [`Plane.indices_separate()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.indices_separate)


            * [`Plane.init_3points()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.init_3points)


            * [`Plane.is_in_list()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.is_in_list)


            * [`Plane.is_in_plane()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.is_in_plane)


            * [`Plane.is_same_plane_as()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.is_same_plane_as)


            * [`Plane.orthonormal_vectors()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.orthonormal_vectors)


            * [`Plane.perpendicular_bisector()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.perpendicular_bisector)


            * [`Plane.project_and_to2dim()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.project_and_to2dim)


            * [`Plane.project_and_to2dim_ordered_indices()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.project_and_to2dim_ordered_indices)


            * [`Plane.projectionpoints()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.Plane.projectionpoints)


        * [`anticlockwise_sort()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.anticlockwise_sort)


        * [`anticlockwise_sort_indices()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.anticlockwise_sort_indices)


        * [`changebasis()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.changebasis)


        * [`collinear()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.collinear)


        * [`diamond_functions()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.diamond_functions)


        * [`function_comparison()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.function_comparison)


        * [`get_lower_and_upper_f()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.get_lower_and_upper_f)


        * [`is_anion_cation_bond()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.is_anion_cation_bond)


        * [`matrixTimesVector()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.matrixTimesVector)


        * [`quarter_ellipsis_functions()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.quarter_ellipsis_functions)


        * [`rectangle_surface_intersection()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.rectangle_surface_intersection)


        * [`rotateCoords()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.rotateCoords)


        * [`rotateCoordsOpt()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.rotateCoordsOpt)


        * [`separation_in_list()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.separation_in_list)


        * [`solid_angle()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.solid_angle)


        * [`sort_separation()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.sort_separation)


        * [`sort_separation_tuple()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.sort_separation_tuple)


        * [`spline_functions()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.spline_functions)


        * [`vectorsToMatrix()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.coordination_geometry_utils.vectorsToMatrix)


    * [pymatgen.analysis.chemenv.utils.defs_utils module](pymatgen.analysis.chemenv.utils.md#module-pymatgen.analysis.chemenv.utils.defs_utils)


        * [`AdditionalConditions`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions)


            * [`AdditionalConditions.ALL`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ALL)


            * [`AdditionalConditions.CONDITION_DESCRIPTION`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.CONDITION_DESCRIPTION)


            * [`AdditionalConditions.NONE`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.NONE)


            * [`AdditionalConditions.NO_AC`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.NO_AC)


            * [`AdditionalConditions.NO_ADDITIONAL_CONDITION`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.NO_ADDITIONAL_CONDITION)


            * [`AdditionalConditions.NO_E2SEB`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.NO_E2SEB)


            * [`AdditionalConditions.NO_ELEMENT_TO_SAME_ELEMENT_BONDS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.NO_ELEMENT_TO_SAME_ELEMENT_BONDS)


            * [`AdditionalConditions.ONLY_ACB`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_ACB)


            * [`AdditionalConditions.ONLY_ACB_AND_NO_E2SEB`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_ACB_AND_NO_E2SEB)


            * [`AdditionalConditions.ONLY_ANION_CATION_BONDS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_ANION_CATION_BONDS)


            * [`AdditionalConditions.ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS)


            * [`AdditionalConditions.ONLY_E2OB`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_E2OB)


            * [`AdditionalConditions.ONLY_ELEMENT_TO_OXYGEN_BONDS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.ONLY_ELEMENT_TO_OXYGEN_BONDS)


            * [`AdditionalConditions.check_condition()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.defs_utils.AdditionalConditions.check_condition)


    * [pymatgen.analysis.chemenv.utils.func_utils module](pymatgen.analysis.chemenv.utils.md#module-pymatgen.analysis.chemenv.utils.func_utils)


        * [`AbstractRatioFunction`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction)


            * [`AbstractRatioFunction.ALLOWED_FUNCTIONS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction.ALLOWED_FUNCTIONS)


            * [`AbstractRatioFunction.evaluate()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction.evaluate)


            * [`AbstractRatioFunction.from_dict()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction.from_dict)


            * [`AbstractRatioFunction.setup_parameters()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.AbstractRatioFunction.setup_parameters)


        * [`CSMFiniteRatioFunction`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction)


            * [`CSMFiniteRatioFunction.ALLOWED_FUNCTIONS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.ALLOWED_FUNCTIONS)


            * [`CSMFiniteRatioFunction.fractions()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.fractions)


            * [`CSMFiniteRatioFunction.mean_estimator()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.mean_estimator)


            * [`CSMFiniteRatioFunction.power2_decreasing_exp()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.power2_decreasing_exp)


            * [`CSMFiniteRatioFunction.ratios()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.ratios)


            * [`CSMFiniteRatioFunction.smootherstep()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.smootherstep)


            * [`CSMFiniteRatioFunction.smoothstep()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMFiniteRatioFunction.smoothstep)


        * [`CSMInfiniteRatioFunction`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction)


            * [`CSMInfiniteRatioFunction.ALLOWED_FUNCTIONS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.ALLOWED_FUNCTIONS)


            * [`CSMInfiniteRatioFunction.fractions()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.fractions)


            * [`CSMInfiniteRatioFunction.mean_estimator()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.mean_estimator)


            * [`CSMInfiniteRatioFunction.power2_inverse_decreasing()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.power2_inverse_decreasing)


            * [`CSMInfiniteRatioFunction.power2_inverse_power2_decreasing()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.power2_inverse_power2_decreasing)


            * [`CSMInfiniteRatioFunction.ratios()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.CSMInfiniteRatioFunction.ratios)


        * [`DeltaCSMRatioFunction`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.DeltaCSMRatioFunction)


            * [`DeltaCSMRatioFunction.ALLOWED_FUNCTIONS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.DeltaCSMRatioFunction.ALLOWED_FUNCTIONS)


            * [`DeltaCSMRatioFunction.smootherstep()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.DeltaCSMRatioFunction.smootherstep)


        * [`RatioFunction`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction)


            * [`RatioFunction.ALLOWED_FUNCTIONS`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.ALLOWED_FUNCTIONS)


            * [`RatioFunction.inverse_smootherstep()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.inverse_smootherstep)


            * [`RatioFunction.inverse_smoothstep()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.inverse_smoothstep)


            * [`RatioFunction.power2_decreasing_exp()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.power2_decreasing_exp)


            * [`RatioFunction.power2_inverse_decreasing()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.power2_inverse_decreasing)


            * [`RatioFunction.power2_inverse_power2_decreasing()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.power2_inverse_power2_decreasing)


            * [`RatioFunction.smootherstep()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.smootherstep)


            * [`RatioFunction.smoothstep()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.func_utils.RatioFunction.smoothstep)


    * [pymatgen.analysis.chemenv.utils.graph_utils module](pymatgen.analysis.chemenv.utils.md#module-pymatgen.analysis.chemenv.utils.graph_utils)


        * [`MultiGraphCycle`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.MultiGraphCycle)


            * [`MultiGraphCycle._is_valid()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.MultiGraphCycle._is_valid)


            * [`MultiGraphCycle.order()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.MultiGraphCycle.order)


            * [`MultiGraphCycle.validate()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.MultiGraphCycle.validate)


        * [`SimpleGraphCycle`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle)


            * [`SimpleGraphCycle._is_valid()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle._is_valid)


            * [`SimpleGraphCycle.as_dict()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle.as_dict)


            * [`SimpleGraphCycle.from_dict()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle.from_dict)


            * [`SimpleGraphCycle.from_edges()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle.from_edges)


            * [`SimpleGraphCycle.order()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle.order)


            * [`SimpleGraphCycle.validate()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.SimpleGraphCycle.validate)


        * [`_c2index_isreverse()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils._c2index_isreverse)


        * [`get_all_elementary_cycles()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.get_all_elementary_cycles)


        * [`get_all_simple_paths_edges()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.get_all_simple_paths_edges)


        * [`get_delta()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.graph_utils.get_delta)


    * [pymatgen.analysis.chemenv.utils.math_utils module](pymatgen.analysis.chemenv.utils.md#module-pymatgen.analysis.chemenv.utils.math_utils)


        * [`_append_es2sequences()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils._append_es2sequences)


        * [`_cartesian_product()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils._cartesian_product)


        * [`_factor_generator()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils._factor_generator)


        * [`cosinus_step()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.cosinus_step)


        * [`divisors()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.divisors)


        * [`get_center_of_arc()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.get_center_of_arc)


        * [`get_linearly_independent_vectors()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.get_linearly_independent_vectors)


        * [`normal_cdf_step()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.normal_cdf_step)


        * [`power2_decreasing_exp()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.power2_decreasing_exp)


        * [`power2_inverse_decreasing()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.power2_inverse_decreasing)


        * [`power2_inverse_power2_decreasing()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.power2_inverse_power2_decreasing)


        * [`power2_inverse_powern_decreasing()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.power2_inverse_powern_decreasing)


        * [`power2_tangent_decreasing()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.power2_tangent_decreasing)


        * [`power3_step()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.power3_step)


        * [`powern_decreasing()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.powern_decreasing)


        * [`powern_parts_step()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.powern_parts_step)


        * [`prime_factors()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.prime_factors)


        * [`scale_and_clamp()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.scale_and_clamp)


        * [`smootherstep()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.smootherstep)


        * [`smoothstep()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.math_utils.smoothstep)


    * [pymatgen.analysis.chemenv.utils.scripts_utils module](pymatgen.analysis.chemenv.utils.md#module-pymatgen.analysis.chemenv.utils.scripts_utils)


        * [`compute_environments()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.scripts_utils.compute_environments)


        * [`draw_cg()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.scripts_utils.draw_cg)


        * [`visualize()`](pymatgen.analysis.chemenv.utils.md#pymatgen.analysis.chemenv.utils.scripts_utils.visualize)