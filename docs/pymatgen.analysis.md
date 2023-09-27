---
layout: default
title: pymatgen.analysis.md
nav_exclude: true
---

1. TOC
{:toc}

# pymatgen.analysis namespace

## Subpackages


* [pymatgen.analysis.chemenv package](pymatgen.analysis.chemenv.md)


    * [Subpackages](pymatgen.analysis.chemenv.md#subpackages)


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


* [pymatgen.analysis.diffraction package](pymatgen.analysis.diffraction.md)




    * [pymatgen.analysis.diffraction.core module](pymatgen.analysis.diffraction.md#module-pymatgen.analysis.diffraction.core)


        * [`AbstractDiffractionPatternCalculator`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator)


            * [`AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL)


            * [`AbstractDiffractionPatternCalculator.TWO_THETA_TOL`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.TWO_THETA_TOL)


            * [`AbstractDiffractionPatternCalculator._abc_impl`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator._abc_impl)


            * [`AbstractDiffractionPatternCalculator.get_pattern()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.get_pattern)


            * [`AbstractDiffractionPatternCalculator.get_plot()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.get_plot)


            * [`AbstractDiffractionPatternCalculator.plot_structures()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.plot_structures)


            * [`AbstractDiffractionPatternCalculator.show_plot()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.show_plot)


        * [`DiffractionPattern`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.DiffractionPattern)


            * [`DiffractionPattern.XLABEL`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.DiffractionPattern.XLABEL)


            * [`DiffractionPattern.YLABEL`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.DiffractionPattern.YLABEL)


        * [`get_unique_families()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.core.get_unique_families)


    * [pymatgen.analysis.diffraction.neutron module](pymatgen.analysis.diffraction.md#module-pymatgen.analysis.diffraction.neutron)


        * [`NDCalculator`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.neutron.NDCalculator)


            * [`NDCalculator._abc_impl`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.neutron.NDCalculator._abc_impl)


            * [`NDCalculator.get_pattern()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.neutron.NDCalculator.get_pattern)


    * [pymatgen.analysis.diffraction.tem module](pymatgen.analysis.diffraction.md#module-pymatgen.analysis.diffraction.tem)


        * [`TEMCalculator`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator)


            * [`TEMCalculator._abc_impl`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator._abc_impl)


            * [`TEMCalculator.bragg_angles()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.bragg_angles)


            * [`TEMCalculator.cell_intensity()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.cell_intensity)


            * [`TEMCalculator.cell_scattering_factors()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.cell_scattering_factors)


            * [`TEMCalculator.electron_scattering_factors()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.electron_scattering_factors)


            * [`TEMCalculator.generate_points()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.generate_points)


            * [`TEMCalculator.get_first_point()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_first_point)


            * [`TEMCalculator.get_interplanar_angle()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_interplanar_angle)


            * [`TEMCalculator.get_interplanar_spacings()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_interplanar_spacings)


            * [`TEMCalculator.get_pattern()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_pattern)


            * [`TEMCalculator.get_plot_2d()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_plot_2d)


            * [`TEMCalculator.get_plot_2d_concise()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_plot_2d_concise)


            * [`TEMCalculator.get_plot_coeffs()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_plot_coeffs)


            * [`TEMCalculator.get_positions()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_positions)


            * [`TEMCalculator.get_s2()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_s2)


            * [`TEMCalculator.is_parallel()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.is_parallel)


            * [`TEMCalculator.normalized_cell_intensity()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.normalized_cell_intensity)


            * [`TEMCalculator.tem_dots()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.tem_dots)


            * [`TEMCalculator.wavelength_rel()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.wavelength_rel)


            * [`TEMCalculator.x_ray_factors()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.x_ray_factors)


            * [`TEMCalculator.zone_axis_filter()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.tem.TEMCalculator.zone_axis_filter)


    * [pymatgen.analysis.diffraction.xrd module](pymatgen.analysis.diffraction.md#module-pymatgen.analysis.diffraction.xrd)


        * [`XRDCalculator`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.xrd.XRDCalculator)


            * [`XRDCalculator.AVAILABLE_RADIATION`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.xrd.XRDCalculator.AVAILABLE_RADIATION)


            * [`XRDCalculator._abc_impl`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.xrd.XRDCalculator._abc_impl)


            * [`XRDCalculator.get_pattern()`](pymatgen.analysis.diffraction.md#pymatgen.analysis.diffraction.xrd.XRDCalculator.get_pattern)


* [pymatgen.analysis.elasticity package](pymatgen.analysis.elasticity.md)




    * [pymatgen.analysis.elasticity.elastic module](pymatgen.analysis.elasticity.md#module-pymatgen.analysis.elasticity.elastic)


        * [`ComplianceTensor`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ComplianceTensor)


        * [`ElasticTensor`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor)


            * [`ElasticTensor.cahill_thermalcond()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.cahill_thermalcond)


            * [`ElasticTensor.clarke_thermalcond()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.clarke_thermalcond)


            * [`ElasticTensor.compliance_tensor`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.compliance_tensor)


            * [`ElasticTensor.debye_temperature()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.debye_temperature)


            * [`ElasticTensor.directional_elastic_mod()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.directional_elastic_mod)


            * [`ElasticTensor.directional_poisson_ratio()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.directional_poisson_ratio)


            * [`ElasticTensor.from_independent_strains()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.from_independent_strains)


            * [`ElasticTensor.from_pseudoinverse()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.from_pseudoinverse)


            * [`ElasticTensor.g_reuss`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.g_reuss)


            * [`ElasticTensor.g_voigt`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.g_voigt)


            * [`ElasticTensor.g_vrh`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.g_vrh)


            * [`ElasticTensor.get_structure_property_dict()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.get_structure_property_dict)


            * [`ElasticTensor.green_kristoffel()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.green_kristoffel)


            * [`ElasticTensor.homogeneous_poisson`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.homogeneous_poisson)


            * [`ElasticTensor.k_reuss`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.k_reuss)


            * [`ElasticTensor.k_voigt`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.k_voigt)


            * [`ElasticTensor.k_vrh`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.k_vrh)


            * [`ElasticTensor.long_v()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.long_v)


            * [`ElasticTensor.property_dict`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.property_dict)


            * [`ElasticTensor.snyder_ac()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.snyder_ac)


            * [`ElasticTensor.snyder_opt()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.snyder_opt)


            * [`ElasticTensor.snyder_total()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.snyder_total)


            * [`ElasticTensor.trans_v()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.trans_v)


            * [`ElasticTensor.universal_anisotropy`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.universal_anisotropy)


            * [`ElasticTensor.y_mod`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensor.y_mod)


        * [`ElasticTensorExpansion`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion)


            * [`ElasticTensorExpansion._abc_impl`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion._abc_impl)


            * [`ElasticTensorExpansion.calculate_stress()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.calculate_stress)


            * [`ElasticTensorExpansion.energy_density()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.energy_density)


            * [`ElasticTensorExpansion.from_diff_fit()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.from_diff_fit)


            * [`ElasticTensorExpansion.get_compliance_expansion()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_compliance_expansion)


            * [`ElasticTensorExpansion.get_effective_ecs()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_effective_ecs)


            * [`ElasticTensorExpansion.get_ggt()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_ggt)


            * [`ElasticTensorExpansion.get_gruneisen_parameter()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_gruneisen_parameter)


            * [`ElasticTensorExpansion.get_heat_capacity()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_heat_capacity)


            * [`ElasticTensorExpansion.get_stability_criteria()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_stability_criteria)


            * [`ElasticTensorExpansion.get_strain_from_stress()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_strain_from_stress)


            * [`ElasticTensorExpansion.get_symmetric_wallace_tensor()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_symmetric_wallace_tensor)


            * [`ElasticTensorExpansion.get_tgt()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_tgt)


            * [`ElasticTensorExpansion.get_wallace_tensor()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_wallace_tensor)


            * [`ElasticTensorExpansion.get_yield_stress()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.get_yield_stress)


            * [`ElasticTensorExpansion.omega()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.omega)


            * [`ElasticTensorExpansion.order`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.order)


            * [`ElasticTensorExpansion.thermal_expansion_coeff()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.ElasticTensorExpansion.thermal_expansion_coeff)


        * [`NthOrderElasticTensor`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor)


            * [`NthOrderElasticTensor.GPa_to_eV_A3`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.GPa_to_eV_A3)


            * [`NthOrderElasticTensor.calculate_stress()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.calculate_stress)


            * [`NthOrderElasticTensor.energy_density()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.energy_density)


            * [`NthOrderElasticTensor.from_diff_fit()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.from_diff_fit)


            * [`NthOrderElasticTensor.order`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.order)


            * [`NthOrderElasticTensor.symbol`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.NthOrderElasticTensor.symbol)


        * [`diff_fit()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.diff_fit)


        * [`find_eq_stress()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.find_eq_stress)


        * [`generate_pseudo()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.generate_pseudo)


        * [`get_diff_coeff()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.get_diff_coeff)


        * [`get_strain_state_dict()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.get_strain_state_dict)


        * [`get_symbol_list()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.get_symbol_list)


        * [`raise_if_unphysical()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.raise_if_unphysical)


        * [`subs()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.elastic.subs)


    * [pymatgen.analysis.elasticity.strain module](pymatgen.analysis.elasticity.md#module-pymatgen.analysis.elasticity.strain)


        * [`Deformation`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Deformation)


            * [`Deformation.apply_to_structure()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Deformation.apply_to_structure)


            * [`Deformation.from_index_amount()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Deformation.from_index_amount)


            * [`Deformation.get_perturbed_indices()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Deformation.get_perturbed_indices)


            * [`Deformation.green_lagrange_strain`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Deformation.green_lagrange_strain)


            * [`Deformation.is_independent()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Deformation.is_independent)


            * [`Deformation.symbol`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Deformation.symbol)


        * [`DeformedStructureSet`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.DeformedStructureSet)


            * [`DeformedStructureSet._abc_impl`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.DeformedStructureSet._abc_impl)


        * [`Strain`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Strain)


            * [`Strain.from_deformation()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Strain.from_deformation)


            * [`Strain.from_index_amount()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Strain.from_index_amount)


            * [`Strain.get_deformation_matrix()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Strain.get_deformation_matrix)


            * [`Strain.symbol`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Strain.symbol)


            * [`Strain.von_mises_strain`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.Strain.von_mises_strain)


        * [`convert_strain_to_deformation()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.strain.convert_strain_to_deformation)


    * [pymatgen.analysis.elasticity.stress module](pymatgen.analysis.elasticity.md#module-pymatgen.analysis.elasticity.stress)


        * [`Stress`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.stress.Stress)


            * [`Stress.dev_principal_invariants`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.stress.Stress.dev_principal_invariants)


            * [`Stress.deviator_stress`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.stress.Stress.deviator_stress)


            * [`Stress.mean_stress`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.stress.Stress.mean_stress)


            * [`Stress.piola_kirchoff_1()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.stress.Stress.piola_kirchoff_1)


            * [`Stress.piola_kirchoff_2()`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.stress.Stress.piola_kirchoff_2)


            * [`Stress.symbol`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.stress.Stress.symbol)


            * [`Stress.von_mises`](pymatgen.analysis.elasticity.md#pymatgen.analysis.elasticity.stress.Stress.von_mises)


* [pymatgen.analysis.ferroelectricity package](pymatgen.analysis.ferroelectricity.md)




    * [pymatgen.analysis.ferroelectricity.polarization module](pymatgen.analysis.ferroelectricity.md#module-pymatgen.analysis.ferroelectricity.polarization)


        * [`EnergyTrend`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.EnergyTrend)


            * [`EnergyTrend.endpoints_minima()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.EnergyTrend.endpoints_minima)


            * [`EnergyTrend.max_spline_jump()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.EnergyTrend.max_spline_jump)


            * [`EnergyTrend.smoothness()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.EnergyTrend.smoothness)


            * [`EnergyTrend.spline()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.EnergyTrend.spline)


        * [`Polarization`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.Polarization)


            * [`Polarization.from_outcars_and_structures()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.from_outcars_and_structures)


            * [`Polarization.get_lattice_quanta()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.get_lattice_quanta)


            * [`Polarization.get_pelecs_and_pions()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.get_pelecs_and_pions)


            * [`Polarization.get_polarization_change()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.get_polarization_change)


            * [`Polarization.get_polarization_change_norm()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.get_polarization_change_norm)


            * [`Polarization.get_same_branch_polarization_data()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.get_same_branch_polarization_data)


            * [`Polarization.max_spline_jumps()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.max_spline_jumps)


            * [`Polarization.same_branch_splines()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.same_branch_splines)


            * [`Polarization.smoothness()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.smoothness)


        * [`PolarizationLattice`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.PolarizationLattice)


            * [`PolarizationLattice._abc_impl`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.PolarizationLattice._abc_impl)


            * [`PolarizationLattice._properties`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.PolarizationLattice._properties)


            * [`PolarizationLattice.get_nearest_site()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.PolarizationLattice.get_nearest_site)


        * [`calc_ionic()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.calc_ionic)


        * [`get_total_ionic_dipole()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.get_total_ionic_dipole)


        * [`zval_dict_from_potcar()`](pymatgen.analysis.ferroelectricity.md#pymatgen.analysis.ferroelectricity.polarization.zval_dict_from_potcar)


* [pymatgen.analysis.gb package](pymatgen.analysis.gb.md)




    * [pymatgen.analysis.gb.grain module](pymatgen.analysis.gb.md#module-pymatgen.analysis.gb.grain)


        * [`GrainBoundary`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary)


            * [`GrainBoundary._abc_impl`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary._abc_impl)


            * [`GrainBoundary._properties`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary._properties)


            * [`GrainBoundary.as_dict()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary.as_dict)


            * [`GrainBoundary.bottom_grain`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary.bottom_grain)


            * [`GrainBoundary.coincidents`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary.coincidents)


            * [`GrainBoundary.copy()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary.copy)


            * [`GrainBoundary.from_dict()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary.from_dict)


            * [`GrainBoundary.get_sorted_structure()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary.get_sorted_structure)


            * [`GrainBoundary.sigma`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary.sigma)


            * [`GrainBoundary.sigma_from_site_prop`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary.sigma_from_site_prop)


            * [`GrainBoundary.top_grain`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundary.top_grain)


        * [`GrainBoundaryGenerator`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator)


            * [`GrainBoundaryGenerator.enum_possible_plane_cubic()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_possible_plane_cubic)


            * [`GrainBoundaryGenerator.enum_sigma_cubic()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_sigma_cubic)


            * [`GrainBoundaryGenerator.enum_sigma_hex()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_sigma_hex)


            * [`GrainBoundaryGenerator.enum_sigma_ort()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_sigma_ort)


            * [`GrainBoundaryGenerator.enum_sigma_rho()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_sigma_rho)


            * [`GrainBoundaryGenerator.enum_sigma_tet()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.enum_sigma_tet)


            * [`GrainBoundaryGenerator.gb_from_parameters()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.gb_from_parameters)


            * [`GrainBoundaryGenerator.get_ratio()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.get_ratio)


            * [`GrainBoundaryGenerator.get_rotation_angle_from_sigma()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.get_rotation_angle_from_sigma)


            * [`GrainBoundaryGenerator.get_trans_mat()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.get_trans_mat)


            * [`GrainBoundaryGenerator.reduce_mat()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.reduce_mat)


            * [`GrainBoundaryGenerator.slab_from_csl()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.slab_from_csl)


            * [`GrainBoundaryGenerator.vec_to_surface()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.GrainBoundaryGenerator.vec_to_surface)


        * [`fix_pbc()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.fix_pbc)


        * [`symm_group_cubic()`](pymatgen.analysis.gb.md#pymatgen.analysis.gb.grain.symm_group_cubic)


* [pymatgen.analysis.interfaces package](pymatgen.analysis.interfaces.md)




    * [pymatgen.analysis.interfaces.coherent_interfaces module](pymatgen.analysis.interfaces.md#module-pymatgen.analysis.interfaces.coherent_interfaces)


        * [`CoherentInterfaceBuilder`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.CoherentInterfaceBuilder)


            * [`CoherentInterfaceBuilder._find_matches()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.CoherentInterfaceBuilder._find_matches)


            * [`CoherentInterfaceBuilder._find_terminations()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.CoherentInterfaceBuilder._find_terminations)


            * [`CoherentInterfaceBuilder.get_interfaces()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.CoherentInterfaceBuilder.get_interfaces)


        * [`from_2d_to_3d()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.from_2d_to_3d)


        * [`get_2d_transform()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.get_2d_transform)


        * [`get_rot_3d_for_2d()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.coherent_interfaces.get_rot_3d_for_2d)


    * [pymatgen.analysis.interfaces.substrate_analyzer module](pymatgen.analysis.interfaces.md#module-pymatgen.analysis.interfaces.substrate_analyzer)


        * [`SubstrateAnalyzer`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateAnalyzer)


            * [`SubstrateAnalyzer.calculate()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateAnalyzer.calculate)


            * [`SubstrateAnalyzer.generate_surface_vectors()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateAnalyzer.generate_surface_vectors)


        * [`SubstrateMatch`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch)


            * [`SubstrateMatch.elastic_energy`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.elastic_energy)


            * [`SubstrateMatch.film_miller`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.film_miller)


            * [`SubstrateMatch.from_zsl()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.from_zsl)


            * [`SubstrateMatch.ground_state_energy`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.ground_state_energy)


            * [`SubstrateMatch.strain`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.strain)


            * [`SubstrateMatch.substrate_miller`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.substrate_miller)


            * [`SubstrateMatch.total_energy`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.total_energy)


            * [`SubstrateMatch.von_mises_strain`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.substrate_analyzer.SubstrateMatch.von_mises_strain)


    * [pymatgen.analysis.interfaces.zsl module](pymatgen.analysis.interfaces.md#module-pymatgen.analysis.interfaces.zsl)


        * [`ZSLGenerator`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLGenerator)


            * [`ZSLGenerator.generate_sl_transformation_sets()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLGenerator.generate_sl_transformation_sets)


            * [`ZSLGenerator.get_equiv_transformations()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLGenerator.get_equiv_transformations)


        * [`ZSLMatch`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLMatch)


            * [`ZSLMatch.film_sl_vectors`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.film_sl_vectors)


            * [`ZSLMatch.film_transformation`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.film_transformation)


            * [`ZSLMatch.film_vectors`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.film_vectors)


            * [`ZSLMatch.match_area`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.match_area)


            * [`ZSLMatch.match_transformation`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.match_transformation)


            * [`ZSLMatch.substrate_sl_vectors`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.substrate_sl_vectors)


            * [`ZSLMatch.substrate_transformation`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.substrate_transformation)


            * [`ZSLMatch.substrate_vectors`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.ZSLMatch.substrate_vectors)


        * [`_bidirectional_same_vectors()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl._bidirectional_same_vectors)


        * [`_unidirectional_is_same_vectors()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl._unidirectional_is_same_vectors)


        * [`fast_norm()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.fast_norm)


        * [`gen_sl_transform_matrices()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.gen_sl_transform_matrices)


        * [`get_factors()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.get_factors)


        * [`is_same_vectors()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.is_same_vectors)


        * [`reduce_vectors()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.reduce_vectors)


        * [`rel_angle()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.rel_angle)


        * [`rel_strain()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.rel_strain)


        * [`vec_angle()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.vec_angle)


        * [`vec_area()`](pymatgen.analysis.interfaces.md#pymatgen.analysis.interfaces.zsl.vec_area)


* [pymatgen.analysis.magnetism package](pymatgen.analysis.magnetism.md)




    * [pymatgen.analysis.magnetism.analyzer module](pymatgen.analysis.magnetism.md#module-pymatgen.analysis.magnetism.analyzer)


        * [`CollinearMagneticStructureAnalyzer`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer)


            * [`CollinearMagneticStructureAnalyzer._round_magmoms()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer._round_magmoms)


            * [`CollinearMagneticStructureAnalyzer.get_exchange_group_info()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_exchange_group_info)


            * [`CollinearMagneticStructureAnalyzer.get_ferromagnetic_structure()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_ferromagnetic_structure)


            * [`CollinearMagneticStructureAnalyzer.get_nonmagnetic_structure()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_nonmagnetic_structure)


            * [`CollinearMagneticStructureAnalyzer.get_structure_with_only_magnetic_atoms()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_structure_with_only_magnetic_atoms)


            * [`CollinearMagneticStructureAnalyzer.get_structure_with_spin()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_structure_with_spin)


            * [`CollinearMagneticStructureAnalyzer.is_magnetic`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.is_magnetic)


            * [`CollinearMagneticStructureAnalyzer.magmoms`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.magmoms)


            * [`CollinearMagneticStructureAnalyzer.magnetic_species_and_magmoms`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.magnetic_species_and_magmoms)


            * [`CollinearMagneticStructureAnalyzer.matches_ordering()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.matches_ordering)


            * [`CollinearMagneticStructureAnalyzer.number_of_magnetic_sites`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.number_of_magnetic_sites)


            * [`CollinearMagneticStructureAnalyzer.number_of_unique_magnetic_sites()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.number_of_unique_magnetic_sites)


            * [`CollinearMagneticStructureAnalyzer.ordering`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.ordering)


            * [`CollinearMagneticStructureAnalyzer.types_of_magnetic_specie`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.types_of_magnetic_specie)


            * [`CollinearMagneticStructureAnalyzer.types_of_magnetic_species`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.types_of_magnetic_species)


        * [`MagneticDeformation`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation)


            * [`MagneticDeformation._asdict()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation._asdict)


            * [`MagneticDeformation._field_defaults`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation._field_defaults)


            * [`MagneticDeformation._fields`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation._fields)


            * [`MagneticDeformation._make()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation._make)


            * [`MagneticDeformation._replace()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation._replace)


            * [`MagneticDeformation.deformation`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation.deformation)


            * [`MagneticDeformation.type`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation.type)


        * [`MagneticStructureEnumerator`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator)


            * [`MagneticStructureEnumerator._generate_ordered_structures()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator._generate_ordered_structures)


            * [`MagneticStructureEnumerator._generate_transformations()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator._generate_transformations)


            * [`MagneticStructureEnumerator._sanitize_input_structure()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator._sanitize_input_structure)


            * [`MagneticStructureEnumerator.available_strategies`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator.available_strategies)


        * [`Ordering`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.Ordering)


            * [`Ordering.AFM`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.Ordering.AFM)


            * [`Ordering.FM`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.Ordering.FM)


            * [`Ordering.FiM`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.Ordering.FiM)


            * [`Ordering.NM`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.Ordering.NM)


            * [`Ordering.Unknown`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.Ordering.Unknown)


        * [`OverwriteMagmomMode`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode)


            * [`OverwriteMagmomMode.none`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.none)


            * [`OverwriteMagmomMode.normalize`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.normalize)


            * [`OverwriteMagmomMode.replace_all`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.replace_all)


            * [`OverwriteMagmomMode.respect_sign`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.respect_sign)


            * [`OverwriteMagmomMode.respect_zero`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.respect_zero)


        * [`magnetic_deformation()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.analyzer.magnetic_deformation)


    * [pymatgen.analysis.magnetism.heisenberg module](pymatgen.analysis.magnetism.md#module-pymatgen.analysis.magnetism.heisenberg)


        * [`HeisenbergMapper`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper)


            * [`HeisenbergMapper._get_exchange_df()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper._get_exchange_df)


            * [`HeisenbergMapper._get_graphs()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper._get_graphs)


            * [`HeisenbergMapper._get_j_exc()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper._get_j_exc)


            * [`HeisenbergMapper._get_nn_dict()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper._get_nn_dict)


            * [`HeisenbergMapper._get_unique_sites()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper._get_unique_sites)


            * [`HeisenbergMapper.estimate_exchange()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.estimate_exchange)


            * [`HeisenbergMapper.get_exchange()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_exchange)


            * [`HeisenbergMapper.get_heisenberg_model()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_heisenberg_model)


            * [`HeisenbergMapper.get_interaction_graph()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_interaction_graph)


            * [`HeisenbergMapper.get_low_energy_orderings()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_low_energy_orderings)


            * [`HeisenbergMapper.get_mft_temperature()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_mft_temperature)


        * [`HeisenbergModel`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel)


            * [`HeisenbergModel._get_j_exc()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel._get_j_exc)


            * [`HeisenbergModel.as_dict()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel.as_dict)


            * [`HeisenbergModel.from_dict()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel.from_dict)


        * [`HeisenbergScreener`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener)


            * [`HeisenbergScreener.screened_structures`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener.screened_structures)


            * [`HeisenbergScreener.screened_energies`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener.screened_energies)


            * [`HeisenbergScreener._do_cleanup()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener._do_cleanup)


            * [`HeisenbergScreener._do_screen()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener._do_screen)


    * [pymatgen.analysis.magnetism.jahnteller module](pymatgen.analysis.magnetism.md#module-pymatgen.analysis.magnetism.jahnteller)


        * [`JahnTellerAnalyzer`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer)


            * [`JahnTellerAnalyzer._estimate_spin_state()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer._estimate_spin_state)


            * [`JahnTellerAnalyzer._get_number_of_d_electrons()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer._get_number_of_d_electrons)


            * [`JahnTellerAnalyzer.get_analysis()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_analysis)


            * [`JahnTellerAnalyzer.get_analysis_and_structure()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_analysis_and_structure)


            * [`JahnTellerAnalyzer.get_magnitude_of_effect_from_species()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_magnitude_of_effect_from_species)


            * [`JahnTellerAnalyzer.get_magnitude_of_effect_from_spin_config()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_magnitude_of_effect_from_spin_config)


            * [`JahnTellerAnalyzer.is_jahn_teller_active()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.is_jahn_teller_active)


            * [`JahnTellerAnalyzer.mu_so()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.mu_so)


            * [`JahnTellerAnalyzer.tag_structure()`](pymatgen.analysis.magnetism.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.tag_structure)


* [pymatgen.analysis.solar package](pymatgen.analysis.solar.md)




    * [pymatgen.analysis.solar.slme module](pymatgen.analysis.solar.md#module-pymatgen.analysis.solar.slme)


        * [`absorption_coefficient()`](pymatgen.analysis.solar.md#pymatgen.analysis.solar.slme.absorption_coefficient)


        * [`get_dir_indir_gap()`](pymatgen.analysis.solar.md#pymatgen.analysis.solar.slme.get_dir_indir_gap)


        * [`optics()`](pymatgen.analysis.solar.md#pymatgen.analysis.solar.slme.optics)


        * [`parse_dielectric_data()`](pymatgen.analysis.solar.md#pymatgen.analysis.solar.slme.parse_dielectric_data)


        * [`slme()`](pymatgen.analysis.solar.md#pymatgen.analysis.solar.slme.slme)


        * [`to_matrix()`](pymatgen.analysis.solar.md#pymatgen.analysis.solar.slme.to_matrix)


* [pymatgen.analysis.structure_prediction package](pymatgen.analysis.structure_prediction.md)




    * [pymatgen.analysis.structure_prediction.dopant_predictor module](pymatgen.analysis.structure_prediction.md#module-pymatgen.analysis.structure_prediction.dopant_predictor)


        * [`_get_dopants()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.dopant_predictor._get_dopants)


        * [`_int_to_roman()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.dopant_predictor._int_to_roman)


        * [`_shannon_radii_from_cn()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.dopant_predictor._shannon_radii_from_cn)


        * [`get_dopants_from_shannon_radii()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.dopant_predictor.get_dopants_from_shannon_radii)


        * [`get_dopants_from_substitution_probabilities()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.dopant_predictor.get_dopants_from_substitution_probabilities)


    * [pymatgen.analysis.structure_prediction.substitution_probability module](pymatgen.analysis.structure_prediction.md#module-pymatgen.analysis.structure_prediction.substitution_probability)


        * [`SubstitutionPredictor`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionPredictor)


            * [`SubstitutionPredictor.composition_prediction()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionPredictor.composition_prediction)


            * [`SubstitutionPredictor.list_prediction()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionPredictor.list_prediction)


        * [`SubstitutionProbability`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionProbability)


    * [pymatgen.analysis.structure_prediction.substitutor module](pymatgen.analysis.structure_prediction.md#module-pymatgen.analysis.structure_prediction.substitutor)


        * [`Substitutor`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor)


            * [`Substitutor._is_charge_balanced()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor._is_charge_balanced)


            * [`Substitutor._is_from_chemical_system()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor._is_from_chemical_system)


            * [`Substitutor.as_dict()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.as_dict)


            * [`Substitutor.from_dict()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.from_dict)


            * [`Substitutor.get_allowed_species()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.get_allowed_species)


            * [`Substitutor.pred_from_comp()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.pred_from_comp)


            * [`Substitutor.pred_from_list()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.pred_from_list)


            * [`Substitutor.pred_from_structures()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.pred_from_structures)


    * [pymatgen.analysis.structure_prediction.volume_predictor module](pymatgen.analysis.structure_prediction.md#module-pymatgen.analysis.structure_prediction.volume_predictor)


        * [`DLSVolumePredictor`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.volume_predictor.DLSVolumePredictor)


            * [`DLSVolumePredictor.get_predicted_structure()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.volume_predictor.DLSVolumePredictor.get_predicted_structure)


            * [`DLSVolumePredictor.predict()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.volume_predictor.DLSVolumePredictor.predict)


        * [`RLSVolumePredictor`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.volume_predictor.RLSVolumePredictor)


            * [`RLSVolumePredictor.get_predicted_structure()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.volume_predictor.RLSVolumePredictor.get_predicted_structure)


            * [`RLSVolumePredictor.predict()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.volume_predictor.RLSVolumePredictor.predict)


        * [`_is_ox()`](pymatgen.analysis.structure_prediction.md#pymatgen.analysis.structure_prediction.volume_predictor._is_ox)


* [pymatgen.analysis.topological package](pymatgen.analysis.topological.md)




    * [pymatgen.analysis.topological.spillage module](pymatgen.analysis.topological.md#module-pymatgen.analysis.topological.spillage)


        * [`SOCSpillage`](pymatgen.analysis.topological.md#pymatgen.analysis.topological.spillage.SOCSpillage)


            * [`SOCSpillage.isclose()`](pymatgen.analysis.topological.md#pymatgen.analysis.topological.spillage.SOCSpillage.isclose)


            * [`SOCSpillage.orth()`](pymatgen.analysis.topological.md#pymatgen.analysis.topological.spillage.SOCSpillage.orth)


            * [`SOCSpillage.overlap_so_spinpol()`](pymatgen.analysis.topological.md#pymatgen.analysis.topological.spillage.SOCSpillage.overlap_so_spinpol)


* [pymatgen.analysis.xas package](pymatgen.analysis.xas.md)




    * [pymatgen.analysis.xas.spectrum module](pymatgen.analysis.xas.md#module-pymatgen.analysis.xas.spectrum)


        * [`XAS`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.XAS)


            * [`XAS.x`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.XAS.x)


            * [`XAS.y`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.XAS.y)


            * [`XAS.absorbing_element`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.XAS.absorbing_element)


            * [`XAS.edge`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.XAS.edge)


            * [`XAS.spectrum_type`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.XAS.spectrum_type)


            * [`XAS.absorbing_index`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.XAS.absorbing_index)


            * [`XAS.XLABEL`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.XAS.XLABEL)


            * [`XAS.YLABEL`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.XAS.YLABEL)


            * [`XAS.stitch()`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.XAS.stitch)


        * [`site_weighted_spectrum()`](pymatgen.analysis.xas.md#pymatgen.analysis.xas.spectrum.site_weighted_spectrum)



## pymatgen.analysis.adsorption module

This module provides classes used to enumerate surface sites and to find
adsorption sites on slabs.


### _class_ AdsorbateSiteFinder(slab, selective_dynamics: bool = False, height: float = 0.9, mi_vec: ArrayLike | None = None)
Bases: `object`

This class finds adsorbate sites on slabs and generates adsorbate
structures according to user-defined criteria.

The algorithm for finding sites is essentially as follows:


    1. Determine “surface sites” by finding those within

        a height threshold along the miller index of the
        highest site


    2. Create a network of surface sites using the Delaunay

        triangulation of the surface sites


    3. Assign on-top, bridge, and hollow adsorption sites

        at the nodes, edges, and face centers of the Del.
        Triangulation


    4. Generate structures from a molecule positioned at

        these sites

Create an AdsorbateSiteFinder object.


* **Parameters**


    * **slab** ([*Slab*](pymatgen.core.md#pymatgen.core.surface.Slab)) – slab object for which to find adsorbate sites


    * **selective_dynamics** (*bool*) – flag for whether to assign
    non-surface sites as fixed for selective dynamics


    * **height** (*float*) – height criteria for selection of surface sites


    * **mi_vec** (*3-D array-like*) – vector corresponding to the vector
    concurrent with the miller index, this enables use with
    slabs that have been reoriented, but the miller vector
    must be supplied manually



#### add_adsorbate(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), ads_coord, repeat=None, translate=True, reorient=True)
Adds an adsorbate at a particular coordinate. Adsorbate represented
by a Molecule object and is translated to (0, 0, 0) if translate is
True, or positioned relative to the input adsorbate coordinate if
translate is False.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – molecule object representing the adsorbate


    * **ads_coord** (*array*) – coordinate of adsorbate position


    * **repeat** (*3-tuple** or **list*) – input for making a supercell of slab
    prior to placing the adsorbate


    * **translate** (*bool*) – flag on whether to translate the molecule so
    that its CoM is at the origin prior to adding it to the surface


    * **reorient** (*bool*) – flag on whether to reorient the molecule to
    have its z-axis concurrent with miller index



#### adsorb_both_surfaces(molecule, repeat=None, min_lw=5.0, translate=True, reorient=True, find_args=None)
Function that generates all adsorption structures for a given
molecular adsorbate on both surfaces of a slab. This is useful for
calculating surface energy where both surfaces need to be equivalent or
if we want to calculate nonpolar systems.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – molecule corresponding to adsorbate


    * **repeat** (*3-tuple** or **list*) – repeat argument for supercell generation


    * **min_lw** (*float*) – minimum length and width of the slab, only used
    if repeat is None


    * **reorient** (*bool*) – flag on whether or not to reorient adsorbate
    along the miller index


    * **find_args** (*dict*) – dictionary of arguments to be passed to the
    call to self.find_adsorption_sites, e.g. {“distance”:2.0}



#### _classmethod_ assign_selective_dynamics(slab)
Helper function to assign selective dynamics site_properties based
on surface, subsurface site properties.


* **Parameters**

    **slab** ([*Slab*](pymatgen.core.md#pymatgen.core.surface.Slab)) – slab for which to assign selective dynamics



#### assign_site_properties(slab, height=0.9)
Assigns site properties.


#### _classmethod_ ensemble_center(site_list, indices, cartesian=True)
Finds the center of an ensemble of sites selected from a list of
sites. Helper method for the find_adsorption_sites algorithm.


* **Parameters**


    * **site_list** (*list** of **sites*) – list of sites


    * **indices** (*list** of **ints*) – list of ints from which to select
    sites from site list


    * **cartesian** (*bool*) – whether to get average fractional or
    Cartesian coordinate



#### find_adsorption_sites(distance=2.0, put_inside=True, symm_reduce=0.01, near_reduce=0.01, positions=('ontop', 'bridge', 'hollow'), no_obtuse_hollow=True)
Finds surface sites according to the above algorithm. Returns a list
of corresponding Cartesian coordinates.


* **Parameters**


    * **distance** (*float*) – distance from the coordinating ensemble
    of atoms along the miller index for the site (i. e.
    the distance from the slab itself)


    * **put_inside** (*bool*) – whether to put the site inside the cell


    * **symm_reduce** (*float*) – symm reduction threshold


    * **near_reduce** (*float*) – near reduction threshold


    * **positions** (*list*) – which positions to include in the site finding
    “ontop”: sites on top of surface sites
    “bridge”: sites at edges between surface sites in Delaunay

    > triangulation of surface sites in the miller plane

    ”hollow”: sites at centers of Delaunay triangulation faces
    “subsurface”: subsurface positions projected into miller plane



    * **no_obtuse_hollow** (*bool*) – flag to indicate whether to include
    obtuse triangular ensembles in hollow sites



#### find_surface_sites_by_height(slab, height=0.9, xy_tol=0.05)
This method finds surface sites by determining which sites are
within a threshold value in height from the topmost site in a list of
sites.


* **Parameters**


    * **site_list** (*list*) – list of sites from which to select surface sites


    * **height** (*float*) – threshold in angstroms of distance from topmost
    site in slab along the slab c-vector to include in surface
    site determination


    * **xy_tol** (*float*) – if supplied, will remove any sites which are
    within a certain distance in the miller plane.



* **Returns**

    list of sites selected to be within a threshold of the highest



#### _classmethod_ from_bulk_and_miller(structure, miller_index, min_slab_size=8.0, min_vacuum_size=10.0, max_normal_search=None, center_slab=True, selective_dynamics=False, undercoord_threshold=0.09)
This method constructs the adsorbate site finder from a bulk
structure and a miller index, which allows the surface sites to be
determined from the difference in bulk and slab coordination, as
opposed to the height threshold.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure from which slab
    input to the ASF is constructed


    * **miller_index** (*3-tuple** or **list*) – miller index to be used


    * **min_slab_size** (*float*) – min slab size for slab generation


    * **min_vacuum_size** (*float*) – min vacuum size for slab generation


    * **max_normal_search** (*int*) – max normal search for slab generation


    * **center_slab** (*bool*) – whether to center slab in slab generation


    * **dynamics** (*selective*) – whether to assign surface sites
    to selective dynamics


    * **undercoord_threshold** (*float*) – threshold of “undercoordation”
    to use for the assignment of surface sites. Default is
    0.1, for which surface sites will be designated if they
    are 10% less coordinated than their bulk counterpart



#### generate_adsorption_structures(molecule, repeat=None, min_lw=5.0, translate=True, reorient=True, find_args=None)
Function that generates all adsorption structures for a given
molecular adsorbate. Can take repeat argument or minimum length/width
of precursor slab as an input.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – molecule corresponding to adsorbate


    * **repeat** (*3-tuple** or **list*) – repeat argument for supercell generation


    * **min_lw** (*float*) – minimum length and width of the slab, only used
    if repeat is None


    * **translate** (*bool*) – flag on whether to translate the molecule so
    that its CoM is at the origin prior to adding it to the surface


    * **reorient** (*bool*) – flag on whether or not to reorient adsorbate
    along the miller index


    * **find_args** (*dict*) – dictionary of arguments to be passed to the
    call to self.find_adsorption_sites, e.g. {“distance”:2.0}



#### generate_substitution_structures(atom, target_species=None, sub_both_sides=False, range_tol=0.01, dist_from_surf=0)
Function that performs substitution-type doping on the surface and
returns all possible configurations where one dopant is substituted per
surface. Can substitute one surface or both.


* **Parameters**


    * **atom** (*str*) – atom corresponding to substitutional dopant


    * **sub_both_sides** (*bool*) – If true, substitute an equivalent
    site on the other surface


    * **target_species** (*list*) – List of specific species to substitute


    * **range_tol** (*float*) – Find viable substitution sites at a specific
    distance from the surface +- this tolerance


    * **dist_from_surf** (*float*) – Distance from the surface to find viable
    substitution sites, defaults to 0 to substitute at the surface



#### get_extended_surface_mesh(repeat=(5, 5, 1))
Gets an extended surface mesh for to use for adsorption site finding
by constructing supercell of surface sites.


* **Parameters**

    **repeat** (*3-tuple*) – repeat for getting extended surface mesh



#### near_reduce(coords_set, threshold=0.0001)
Prunes coordinate set for coordinates that are within threshold.


* **Parameters**


    * **coords_set** (*Nx3 array-like*) – list or array of coordinates


    * **threshold** (*float*) – threshold value for distance



#### subsurface_sites()
Convenience method to return list of subsurface sites.


#### _property_ surface_sites()
Convenience method to return a list of surface sites.


#### symm_reduce(coords_set, threshold=1e-06)
Reduces the set of adsorbate sites by finding removing symmetrically
equivalent duplicates.


* **Parameters**


    * **coords_set** – coordinate set in Cartesian coordinates


    * **threshold** – tolerance for distance equivalence, used
    as input to in_coord_list_pbc for dupl. checking



### get_mi_vec(slab)
Convenience function which returns the unit vector aligned with the
miller index.


### get_rot(slab)
Gets the transformation to rotate the z axis into the miller index.


### plot_slab(slab, ax, scale=0.8, repeat=5, window=1.5, draw_unit_cell=True, decay=0.2, adsorption_sites=True, inverse=False)
Function that helps visualize the slab in a 2-D plot, for convenient
viewing of output of AdsorbateSiteFinder.


* **Parameters**


    * **slab** (*slab*) – Slab object to be visualized


    * **ax** (*axes*) – matplotlib axes with which to visualize


    * **scale** (*float*) – radius scaling for sites


    * **repeat** (*int*) – number of repeating unit cells to visualize


    * **window** (*float*) – window for setting the axes limits, is essentially
    a fraction of the unit cell limits


    * **draw_unit_cell** (*bool*) – flag indicating whether or not to draw cell


    * **decay** (*float*) – how the alpha-value decays along the z-axis


    * **inverse** (*bool*) – invert z axis to plot opposite surface



### put_coord_inside(lattice, cart_coordinate)
Converts a Cartesian coordinate such that it is inside the unit cell.


### reorient_z(structure)
Reorients a structure such that the z axis is concurrent with the normal
to the A-B plane.

## pymatgen.analysis.bond_dissociation module

Module for BondDissociationEnergies.


### _class_ BondDissociationEnergies(molecule_entry: dict[str, str | dict[str, str | int]], fragment_entries: list[dict[str, str | dict[str, str | int]]], allow_additional_charge_separation: bool = False, multibreak: bool = False)
Bases: `MSONable`

Standard constructor for bond dissociation energies. All bonds in the principle molecule are
looped through and their dissociation energies are calculated given the energies of the resulting
fragments, or, in the case of a ring bond, from the energy of the molecule obtained from breaking
the bond and opening the ring. This class should only be called after the energies of the optimized
principle molecule and all relevant optimized fragments have been determined, either from quantum
chemistry or elsewhere. It was written to provide the analysis after running an Atomate fragmentation
workflow.

Note that the entries passed by the user must have the following keys: formula_pretty, initial_molecule,
final_molecule. If a PCM is present, all entries should also have a pcm_dielectric key.


* **Parameters**


    * **molecule_entry** (*dict*) – Entry for the principle molecule. Should have the keys mentioned above.


    * **fragment_entries** (*list** of **dicts*) – List of fragment entries. Each should have the keys mentioned above.


    * **allow_additional_charge_separation** (*bool*) – If True, consider larger than normal charge separation
    among fragments. Defaults to False. See the definition of self.expected_charges below for more
    specific information.


    * **multibreak** (*bool*) – If True, additionally attempt to break pairs of bonds. Defaults to False.



#### build_new_entry(frags, bonds)
Simple function to format a bond dissociation entry that will eventually be returned to the user.


* **Parameters**


    * **frags** –


    * **bonds** –



#### filter_fragment_entries(fragment_entries)
Filter the fragment entries.


* **Parameters**

    **fragment_entries** –



#### fragment_and_process(bonds)
Fragment and process bonds.


* **Parameters**

    **bonds** – Bonds to process.



#### search_fragment_entries(frag)
Search all fragment entries for those isomorphic to the given fragment.
We distinguish between entries where both initial and final molgraphs are isomorphic to the
given fragment (entries) vs those where only the initial molgraph is isomorphic to the given
fragment (initial_entries) vs those where only the final molgraph is isomorphic (final_entries).


* **Parameters**

    **frag** – Fragment


## pymatgen.analysis.bond_valence module

This module implements classes to perform bond valence analyses.


### _class_ BVAnalyzer(symm_tol=0.1, max_radius=4, max_permutations=100000, distance_scale_factor=1.015, charge_neutrality_tolerance=1e-05, forbidden_species=None)
Bases: `object`

This class implements a maximum a posteriori (MAP) estimation method to
determine oxidation states in a structure. The algorithm is as follows:
1) The bond valence sum of all symmetrically distinct sites in a structure
is calculated using the element-based parameters in M. O’Keefe, & N. Brese,
JACS, 1991, 113(9), 3226-3229. doi:10.1021/ja00009a002.
2) The posterior probabilities of all oxidation states is then calculated
using: P(oxi_state/BV) = K \* P(BV/oxi_state) \* P(oxi_state), where K is
a constant factor for each element. P(BV/oxi_state) is calculated as a
Gaussian with mean and std deviation determined from an analysis of
the ICSD. The posterior P(oxi_state) is determined from a frequency
analysis of the ICSD.
3) The oxidation states are then ranked in order of decreasing probability
and the oxidation state combination that result in a charge neutral cell
is selected.

Initializes the BV analyzer, with useful defaults.


* **Parameters**


    * **symm_tol** – Symmetry tolerance used to determine which sites are
    symmetrically equivalent. Set to 0 to turn off symmetry.


    * **max_radius** – Maximum radius in Angstrom used to find nearest neighbors.


    * **max_permutations** – The maximum number of permutations of oxidation states to test.


    * **distance_scale_factor** – A scale factor to be applied. This is useful for scaling
    distances, esp in the case of calculation-relaxed structures
    which may tend to under (GGA) or over bind (LDA). The default
    of 1.015 works for GGA. For experimental structure, set this to
    1.


    * **charge_neutrality_tolerance** – Tolerance on the charge neutrality when unordered structures
    are at stake.


    * **forbidden_species** – List of species that are forbidden (example : [“O-”] cannot be
    used) It is used when e.g. someone knows that some oxidation
    state cannot occur for some atom in a structure or list of
    structures.



#### CHARGE_NEUTRALITY_TOLERANCE(_ = 1e-0_ )

#### _calc_site_probabilities(site, nn)

#### _calc_site_probabilities_unordered(site, nn)

#### get_oxi_state_decorated_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Get an oxidation state decorated structure. This currently works only
for ordered structures only.


* **Parameters**

    **structure** – Structure to analyze



* **Returns**

    A modified structure that is oxidation state decorated.



* **Raises**

    **ValueError if the valences cannot be determined.** –



#### get_valences(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Returns a list of valences for each site in the structure.


* **Parameters**

    **structure** – Structure to analyze



* **Returns**

    A list of valences for each site in the structure (for an ordered structure),
    e.g., [1, 1, -2] or a list of lists with the valences for each fractional
    element of each site in the structure (for an unordered structure), e.g., [[2,
    4], [3], [-2], [-2], [-2]]



* **Raises**

    **A ValueError if the valences cannot be determined.** –



### add_oxidation_state_by_site_fraction(structure, oxidation_states)
Add oxidation states to a structure by fractional site.


* **Parameters**

    **oxidation_states** (*list*) – List of list of oxidation states for each
    site fraction for each site.
    E.g., [[2, 4], [3], [-2], [-2], [-2]]



### calculate_bv_sum(site, nn_list, scale_factor=1.0)
Calculates the BV sum of a site.


* **Parameters**


    * **site** ([*PeriodicSite*](pymatgen.core.md#pymatgen.core.sites.PeriodicSite)) – The central site to calculate the bond valence


    * **nn_list** (*[*[*Neighbor*](pymatgen.core.md#pymatgen.core.structure.Neighbor)*]*) – A list of namedtuple Neighbors having “distance”
    and “site” attributes


    * **scale_factor** (*float*) – A scale factor to be applied. This is useful for
    scaling distance, esp in the case of calculation-relaxed structures
    which may tend to under (GGA) or over bind (LDA).



### calculate_bv_sum_unordered(site, nn_list, scale_factor=1)
Calculates the BV sum of a site for unordered structures.


* **Parameters**


    * **site** ([*PeriodicSite*](pymatgen.core.md#pymatgen.core.sites.PeriodicSite)) – The central site to calculate the bond valence


    * **nn_list** (*[*[*Neighbor*](pymatgen.core.md#pymatgen.core.structure.Neighbor)*]*) – A list of namedtuple Neighbors having “distance”
    and “site” attributes


    * **scale_factor** (*float*) – A scale factor to be applied. This is useful for
    scaling distance, esp in the case of calculation-relaxed structures
    which may tend to under (GGA) or over bind (LDA).



### get_z_ordered_elmap(comp)
Arbitrary ordered element map on the elements/species of a composition of a
given site in an unordered structure. Returns a list of tuples (
element_or_specie: occupation) in the arbitrary order.

The arbitrary order is based on the Z of the element and the smallest
fractional occupations first.
Example : {“Ni3+”: 0.2, “Ni4+”: 0.2, “Cr3+”: 0.15, “Zn2+”: 0.34,
“Cr4+”: 0.11} will yield the species in the following order :
Cr4+, Cr3+, Ni3+, Ni4+, Zn2+ … or
Cr4+, Cr3+, Ni4+, Ni3+, Zn2+

## pymatgen.analysis.chempot_diagram module

This module implements the construction and plotting of chemical potential diagrams
from a list of entries within a chemical system containing 2 or more elements. The
chemical potential diagram is the mathematical dual to the traditional compositional
phase diagram.

For more information, please cite/reference the paper below:

> Todd, P. K., McDermott, M. J., Rom, C. L., Corrao, A. A., Denney, J. J., Dwaraknath,
> S. S.,  Khalifah, P. G., Persson, K. A., & Neilson, J. R. (2021). Selectivity in
> Yttrium Manganese Oxide Synthesis via Local Chemical Potentials in Hyperdimensional
> Phase Space. Journal of the American Chemical Society, 143(37), 15185-15194.
> [https://doi.org/10.1021/jacs.1c06229](https://doi.org/10.1021/jacs.1c06229)

Please also consider referencing the original 1999 paper by H. Yokokawa,
who outlined many of its possible uses:

> Yokokawa, H. “Generalized chemical potential diagram and its applications to
> chemical reactions at interfaces between dissimilar materials.” JPE 20,
> 258 (1999). [https://doi.org/10.1361/105497199770335794](https://doi.org/10.1361/105497199770335794)


### _class_ ChemicalPotentialDiagram(entries: list[PDEntry], limits: dict[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element), tuple[float, float]] | None = None, default_min_limit: float = -50.0, formal_chempots: bool = True)
Bases: `MSONable`

The chemical potential diagram is the mathematical dual to the compositional
phase diagram. To create the diagram, convex minimization is
performed in energy (E) vs. chemical potential (μ) space by taking the lower convex
envelope of hyperplanes. Accordingly, “points” on the compositional phase diagram
become N-dimensional convex polytopes (domains) in chemical potential space.

For more information on this specific implementation of the algorithm,
please cite/reference the paper below:

> Todd, P. K., McDermott, M. J., Rom, C. L., Corrao, A. A., Denney, J. J., Dwaraknath,
> S. S.,  Khalifah, P. G., Persson, K. A., & Neilson, J. R. (2021). Selectivity in
> Yttrium Manganese Oxide Synthesis via Local Chemical Potentials in Hyperdimensional
> Phase Space. Journal of the American Chemical Society, 143(37), 15185-15194.
> [https://doi.org/10.1021/jacs.1c06229](https://doi.org/10.1021/jacs.1c06229)


* **Parameters**


    * **entries** (*list**[**PDEntry**]*) – PDEntry-like objects containing a composition and
    energy. Must contain elemental references and be suitable for typical
    phase diagram construction. Entries must be within a chemical system
    of with 2+ elements.


    * **limits** (*dict**[*[*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)*, **float**] **| **None*) – Bounds of elemental chemical potentials (min, max),
    which are used to construct the border hyperplanes used in the HalfSpaceIntersection
    algorithm; these constrain the space over which the domains are calculated and also
    determine the size of the plotted diagram. Any elemental limits not specified are
    covered in the default_min_limit argument. e.g., {Element(“Li”): [-12.0, 0.0], …}


    * **default_min_limit** (*float*) – Default minimum chemical potential limit (i.e.,
    lower bound) for unspecified elements within the “limits” argument.


    * **formal_chempots** (*bool*) – Whether to plot the formal (‘reference’) chemical potentials
    (i.e. μ_X - μ_X^0) or the absolute DFT reference energies (i.e. μ_X(DFT)).
    Default is True (i.e. plot formal chemical potentials).



#### _static_ _get_2d_domain_lines(draw_domains)
Returns a list of Scatter objects tracing the domain lines on a
2-dimensional chemical potential diagram.


#### _get_2d_plot(elements: list[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)], label_stable: bool | None, element_padding: float | None)
Returns a Plotly figure for a 2-dimensional chemical potential diagram.


#### _static_ _get_3d_domain_lines(domains: dict[str, list[[Simplex](pymatgen.util.md#pymatgen.util.coord.Simplex)] | None])
Returns a list of Scatter3d objects tracing the domain lines on a
3-dimensional chemical potential diagram.


#### _static_ _get_3d_domain_simplexes_and_ann_loc(points_3d: ndarray)
Returns a list of Simplex objects and coordinates of annotation for one
domain in a 3-d chemical potential diagram. Uses PCA to project domain
into 2-dimensional space so that ConvexHull can be used to identify the
bounding polygon.


#### _static_ _get_3d_formula_lines(draw_domains: dict[str, np.ndarray], formula_colors: list[str] | None)
Returns a list of Scatter3d objects defining the bounding polyhedra.


#### _static_ _get_3d_formula_meshes(draw_domains: dict[str, np.ndarray], formula_colors: list[str] | None)
Returns a list of Mesh3d objects for the domains specified by the
user (i.e., draw_domains).


#### _get_3d_plot(elements: list[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)], label_stable: bool | None, formulas_to_draw: list[str] | None, draw_formula_meshes: bool | None, draw_formula_lines: bool | None, formula_colors: list[str] | None, element_padding: float | None)
Returns a Plotly figure for a 3-dimensional chemical potential diagram.


#### _static_ _get_annotation(ann_loc: np.ndarray, formula: str)
Returns a Plotly annotation dict given a formula and location.


#### _static_ _get_axis_layout_dict(elements: list[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)])
Returns a Plotly layout dict for either 2-d or 3-d axes.


#### _get_border_hyperplanes()
Returns an array of the bounding hyperplanes given by elemental limits.


#### _get_domains()
Returns a dictionary of domains as {formula: np.ndarray}.


#### _get_hyperplanes_and_entries()
Returns both the array of hyperplanes, as well as a list of the minimum
entries.


#### _static_ _get_min_entries_and_el_refs(entries: list[PDEntry])
Returns a list of the minimum-energy entries at each composition and the
entries corresponding to the elemental references.


#### _static_ _get_new_limits_from_padding(domains: dict[str, ndarray], elem_indices: list[int], element_padding: float, default_min_limit: float)
Gets new minimum limits for each element by subtracting specified padding
from the minimum for each axis found in any of the domains.


#### _property_ border_hyperplanes(_: ndarra_ )
Returns bordering hyperplanes.


#### _property_ chemical_system(_: st_ )
Returns the chemical system (A-B-C-…) of diagram object.


#### _property_ domains(_: dict[str, ndarray_ )
Mapping of formulas to array of domain boundary points.


#### _property_ el_refs(_: dict[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element), PDEntry_ )
Returns a dictionary of elements and reference entries.


#### _property_ entry_dict(_: dict[str, [ComputedEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)_ )
Mapping between reduced formula and ComputedEntry.


#### get_plot(elements: list[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element) | str] | None = None, label_stable: bool | None = True, formulas_to_draw: list[str] | None = None, draw_formula_meshes: bool | None = True, draw_formula_lines: bool | None = True, formula_colors: list[str] = ['rgb(27,158,119)', 'rgb(217,95,2)', 'rgb(117,112,179)', 'rgb(231,41,138)', 'rgb(102,166,30)', 'rgb(230,171,2)', 'rgb(166,118,29)', 'rgb(102,102,102)'], element_padding: float | None = 1.0)
Plot the 2-dimensional or 3-dimensional chemical potential diagram using an
interactive Plotly interface.

Elemental axes can be specified; if none provided, will automatically default
to first 2-3 elements within the “elements” attribute.

In 3D, this method also allows for plotting of lower-dimensional “slices” of
hyperdimensional polytopes (e.g., the LiMnO2 domain within a Y-Mn-O diagram).
This allows for visualization of some of the phase boundaries that can only
be seen fully in high dimensional space; see the “formulas_to_draw” argument.


* **Parameters**


    * **elements** – list of elements to use as axes in the diagram. If None,
    automatically defaults to the first 2 or elements within the
    object’s “elements” attribute.


    * **label_stable** – whether to label stable phases by their reduced
    formulas. Defaults to True.


    * **formulas_to_draw** – for 3-dimensional diagrams, an optional list of
    formulas to plot on the diagram; if these are from a different
    chemical system a 3-d polyhedron “slice” will be plotted. Defaults to None.


    * **draw_formula_meshes** – whether to draw a colored mesh for the
    optionally specified formulas_to_draw. Defaults to True.


    * **draw_formula_lines** – whether to draw bounding lines for the
    optionally specified formulas_to_draw. Defaults to True.


    * **formula_colors** – a list of colors to use in the plotting of the optionally
    specified formulas_to-draw. Defaults to the Plotly Dark2 color scheme.


    * **element_padding** – if provided, automatically adjusts chemical potential axis
    limits of the plot such that elemental domains have the specified padding
    (in eV/atom), helping provide visual clarity. Defaults to 1.0.



* **Returns**

    A Plotly Figure object



#### _property_ hyperplane_entries(_: list[PDEntry_ )
Returns list of entries corresponding to hyperplanes.


#### _property_ hyperplanes(_: ndarra_ )
Returns array of hyperplane data.


#### _property_ lims(_: ndarra_ )
Returns array of limits used in constructing hyperplanes.


### _renormalize_entry(entry: PDEntry, renormalization_energy_per_atom: float)
Regenerate the input entry with an energy per atom decreased by renormalization_energy_per_atom.


### get_2d_orthonormal_vector(line_pts: ndarray)
Calculates a vector that is orthonormal to a line given by a set of points. Used
for determining the location of an annotation on a 2-d chemical potential diagram.


* **Parameters**

    **line_pts** – a 2x2 array in the form of [[x0, y0], [x1, y1]] giving the
    coordinates of a line



* **Returns**

    A length-2 vector that is orthonormal to the line.



* **Return type**

    np.ndarray



### get_centroid_2d(vertices: ndarray)
A bare-bones implementation of the formula for calculating the centroid of a 2D
polygon. Useful for calculating the location of an annotation on a chemical
potential domain within a 3D chemical potential diagram.

**NOTE**: vertices must be ordered circumferentially!


* **Parameters**

    **vertices** – array of 2-d coordinates corresponding to a polygon, ordered
    circumferentially



* **Returns**

    Array giving 2-d centroid coordinates



### simple_pca(data: ndarray, k: int = 2)
A bare-bones implementation of principal component analysis (PCA) used in the
ChemicalPotentialDiagram class for plotting.


* **Parameters**


    * **data** – array of observations


    * **k** – Number of principal components returned



* **Returns**

    Tuple of projected data, eigenvalues, eigenvectors


## pymatgen.analysis.cost module

This module is used to estimate the cost of various compounds. Costs are taken
from the a CostDB instance, for example a CSV file via CostDBCSV.
For compounds with no cost listed, a Phase Diagram style convex hull
optimization is performed to determine a set of compositions that can be mixed
to give the desired compound with lowest total cost.


### _class_ CostAnalyzer(costdb)
Bases: `object`

Given a CostDB, figures out the minimum cost solutions via convex hull.


* **Parameters**

    **(****)** (*costdb*) – Cost database.



#### get_cost_per_kg(comp)
Get best estimate of minimum cost/kg based on known data.


* **Parameters**

    **comp** – Composition as a pymatgen.core.structure.Composition



* **Returns**

    float of cost/kg



#### get_cost_per_mol(comp)
Get best estimate of minimum cost/mol based on known data.


* **Parameters**

    **comp** – Composition as a pymatgen.core.structure.Composition



* **Returns**

    float of cost/mol



#### get_lowest_decomposition(composition)
Get the decomposition leading to lowest cost.


* **Parameters**

    **composition** – Composition as a pymatgen.core.structure.Composition



* **Returns**

    amount}



* **Return type**

    Decomposition as a dict of {Entry



### _class_ CostDB()
Bases: `object`

Abstract class for representing a Cost database.
Can be extended, e.g. for file-based or REST-based databases.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract_ get_entries(chemsys)
For a given chemical system, return an array of CostEntries.


* **Parameters**

    **chemsys** – array of Elements defining the chemical system.



* **Returns**

    array of CostEntries



### _class_ CostDBCSV(filename)
Bases: `CostDB`

Read a CSV file to get costs
Format is formula,cost_per_kg,name,BibTeX.


* **Parameters**

    **filename** (*str*) – Filename of cost database.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### get_entries(chemsys)
For a given chemical system, return an array of CostEntries.


* **Parameters**

    **chemsys** – array of Elements defining the chemical system.



* **Returns**

    array of CostEntries



### _class_ CostEntry(composition, cost, name, reference)
Bases: `PDEntry`

Extends PDEntry to include a BibTeX reference and include language about cost.


* **Parameters**


    * **composition** – Composition as a pymatgen.core.structure.Composition


    * **cost** – Cost (per mol, NOT per kg) of the full Composition


    * **name** – Optional parameter to name the entry. Defaults to the reduced
    chemical formula as in PDEntry.


    * **reference** – Reference data as BiBTeX string.



#### _abc_impl(_ = <_abc._abc_data object_ )
## pymatgen.analysis.dimensionality module

This module provides functions to get the dimensionality of a structure.

A number of different algorithms are implemented. These are based on the
following publications:

get_dimensionality_larsen:


    * P. M. Larsen, M. Pandey, M. Strange, K. W. Jacobsen. Definition of a
    scoring parameter to identify low-dimensional materials components.
    Phys. Rev. Materials 3, 034003 (2019).

get_dimensionality_cheon:


    * Cheon, G.; Duerloo, K.-A. N.; Sendek, A. D.; Porter, C.; Chen, Y.; Reed,
    E. J. Data Mining for New Two- and One-Dimensional Weakly Bonded Solids
    and Lattice-Commensurate Heterostructures. Nano Lett. 2017.

get_dimensionality_gorai:


    * Gorai, P., Toberer, E. & Stevanovic, V. Computational Identification of
    Promising Thermoelectric Materials Among Known Quasi-2D Binary Compounds.
    J. Mater. Chem. A 2, 4136 (2016).


### calculate_dimensionality_of_site(bonded_structure, site_index, inc_vertices=False)
Calculates the dimensionality of the component containing the given site.

Implements directly the modified breadth-first-search algorithm described in
Algorithm 1 of:

P. M. Larsen, M. Pandey, M. Strange, K. W. Jacobsen. Definition of a
scoring parameter to identify low-dimensional materials components.
Phys. Rev. Materials 3, 034003 (2019).


* **Parameters**


    * **bonded_structure** (*StructureGraph*) – A structure with bonds, represented
    as a pymatgen structure graph. For example, generated using the
    CrystalNN.get_bonded_structure() method.


    * **site_index** (*int*) – The index of a site in the component of interest.


    * **inc_vertices** (*bool**, **optional*) – Whether to return the vertices (site
    images) of the component.



* **Returns**

    If inc_vertices is False, the dimensionality of the
    component will be returned as an int. If inc_vertices is true, the
    function will return a tuple of (dimensionality, vertices), where
    vertices is a list of tuples. E.g. [(0, 0, 0), (1, 1, 1)].



* **Return type**

    (int or tuple)



### find_clusters(struct, connected_matrix)
Finds bonded clusters of atoms in the structure with periodic boundary
conditions.

If there are atoms that are not bonded to anything, returns [0,1,0]. (For
faster computation time)

Author: Gowoon Cheon
Email: [gcheon@stanford.edu](mailto:gcheon@stanford.edu)


* **Parameters**


    * **struct** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **connected_matrix** – Must be made from the same structure with
    find_connected_atoms() function.



* **Returns**

    the size of the largest cluster in the crystal structure
    min_cluster: the size of the smallest cluster in the crystal structure
    clusters: list of bonded clusters found here, clusters are formatted as
    sets of indices of atoms



* **Return type**

    max_cluster



### find_connected_atoms(struct, tolerance=0.45, ldict=None)
Finds bonded atoms and returns a adjacency matrix of bonded atoms.

Author: “Gowoon Cheon”
Email: “[gcheon@stanford.edu](mailto:gcheon@stanford.edu)”


* **Parameters**


    * **struct** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **tolerance** – length in angstroms used in finding bonded atoms. Two atoms
    are considered bonded if (radius of atom 1) + (radius of atom 2) +
    (tolerance) < (distance between atoms 1 and 2). Default
    value = 0.45, the value used by JMol and Cheon et al.


    * **ldict** – dictionary of bond lengths used in finding bonded atoms. Values
    from JMol are used as default



* **Returns**

    A numpy array of shape (number of atoms, number of atoms);
    If any image of atom j is bonded to atom i with periodic boundary
    conditions, the matrix element [atom i, atom j] is 1.



* **Return type**

    (np.ndarray)



### get_dimensionality_cheon(structure_raw, tolerance=0.45, ldict=None, standardize=True, larger_cell=False)
Algorithm for finding the dimensions of connected subunits in a structure.
This method finds the dimensionality of the material even when the material
is not layered along low-index planes, or does not have flat
layers/molecular wires.

Author: “Gowoon Cheon”
Email: “[gcheon@stanford.edu](mailto:gcheon@stanford.edu)”

See details at :

Cheon, G.; Duerloo, K.-A. N.; Sendek, A. D.; Porter, C.; Chen, Y.; Reed,
E. J. Data Mining for New Two- and One-Dimensional Weakly Bonded Solids and
Lattice-Commensurate Heterostructures. Nano Lett. 2017.


* **Parameters**


    * **structure_raw** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – A pymatgen Structure object.


    * **tolerance** (*float*) – length in angstroms used in finding bonded atoms.
    Two atoms are considered bonded if (radius of atom 1) + (radius of
    atom 2) + (tolerance) < (distance between atoms 1 and 2). Default
    value = 0.45, the value used by JMol and Cheon et al.


    * **ldict** (*dict*) – dictionary of bond lengths used in finding bonded atoms.
    Values from JMol are used as default


    * **standardize** – works with conventional standard structures if True. It is
    recommended to keep this as True.


    * **larger_cell** – tests with 3x3x3 supercell instead of 2x2x2. Testing with
    2x2x2 supercell is faster but misclassifies rare interpenetrated 3D

    > structures. Testing with a larger cell circumvents this problem




* **Returns**

    dimension of the largest cluster as a string. If there are ions
    or molecules it returns ‘intercalated ion/molecule’



* **Return type**

    (str)



### get_dimensionality_gorai(structure, max_hkl=2, el_radius_updates=None, min_slab_size=5, min_vacuum_size=5, standardize=True, bonds=None)
This method returns whether a structure is 3D, 2D (layered), or 1D (linear
chains or molecules) according to the algorithm published in Gorai, P.,
Toberer, E. & Stevanovic, V. Computational Identification of Promising
Thermoelectric Materials Among Known Quasi-2D Binary Compounds. J. Mater.
Chem. A 2, 4136 (2016).

Note that a 1D structure detection might indicate problems in the bonding
algorithm, particularly for ionic crystals (e.g., NaCl)

Users can change the behavior of bonds detection by passing either
el_radius_updates to update atomic radii for auto-detection of max bond
distances, or bonds to explicitly specify max bond distances for atom pairs.
Note that if you pass both, el_radius_updates are ignored.


* **Parameters**


    * **structure** – (Structure) structure to analyze dimensionality for


    * **max_hkl** – (int) max index of planes to look for layers


    * **el_radius_updates** – (dict) symbol->float to update atomic radii


    * **min_slab_size** – (float) internal surface construction parameter


    * **min_vacuum_size** – (float) internal surface construction parameter


    * **standardize** (*bool*) – whether to standardize the structure before
    analysis. Set to False only if you already have the structure in a
    convention where layers / chains will be along low <hkl> indexes.


    * **bonds** (*{**(**specie1**, **specie2*) – max_bond_dist}: bonds are
    specified as a dict of tuples: float of specie1, specie2
    and the max bonding distance. For example, PO4 groups may be
    defined as {(“P”, “O”): 3}.



* **Returns**

    the dimensionality of the structure - 1 (molecules/chains),

        2 (layered), or 3 (3D)




* **Return type**

    int



### get_dimensionality_larsen(bonded_structure)
Gets the dimensionality of a bonded structure.

The dimensionality of the structure is the highest dimensionality of all
structure components. This method is very robust and can handle
many tricky structures, regardless of structure type or improper connections
due to periodic boundary conditions.

Requires a StructureGraph object as input. This can be generated using one
of the NearNeighbor classes. For example, using the CrystalNN class:

```default
bonded_structure = CrystalNN().get_bonded_structure(structure)
```

Based on the modified breadth-first-search algorithm described in:

P. M. Larsen, M. Pandey, M. Strange, K. W. Jacobsen. Definition of a
scoring parameter to identify low-dimensional materials components.
Phys. Rev. Materials 3, 034003 (2019).


* **Parameters**

    **bonded_structure** (*StructureGraph*) – A structure with bonds, represented
    as a pymatgen structure graph. For example, generated using the
    CrystalNN.get_bonded_structure() method.



* **Returns**

    The dimensionality of the structure.



* **Return type**

    (int)



### get_structure_components(bonded_structure, inc_orientation=False, inc_site_ids=False, inc_molecule_graph=False)
Gets information on the components in a bonded structure.

Correctly determines the dimensionality of all structures, regardless of
structure type or improper connections due to periodic boundary conditions.

Requires a StructureGraph object as input. This can be generated using one
of the NearNeighbor classes. For example, using the CrystalNN class:

```default
bonded_structure = CrystalNN().get_bonded_structure(structure)
```

Based on the modified breadth-first-search algorithm described in:

P. M. Larsen, M. Pandey, M. Strange, K. W. Jacobsen. Definition of a
scoring parameter to identify low-dimensional materials components.
Phys. Rev. Materials 3, 034003 (2019).


* **Parameters**


    * **bonded_structure** (*StructureGraph*) – A structure with bonds, represented
    as a pymatgen structure graph. For example, generated using the
    CrystalNN.get_bonded_structure() method.


    * **inc_orientation** (*bool**, **optional*) – Whether to include the orientation
    of the structure component. For surfaces, the miller index is given,
    for one-dimensional structures, the direction of the chain is given.


    * **inc_site_ids** (*bool**, **optional*) – Whether to include the site indices
    of the sites in the structure component.


    * **inc_molecule_graph** (*bool**, **optional*) – Whether to include MoleculeGraph
    objects for zero-dimensional components.



* **Returns**

    Information on the components in a structure as a list
    of dictionaries with the keys:


    * ”structure_graph”: A pymatgen StructureGraph object for the

        component.


    * ”dimensionality”: The dimensionality of the structure component as an

        int.


    * ”orientation”: If inc_orientation is True, the orientation of the

        component as a tuple. E.g. (1, 1, 1)


    * ”site_ids”: If inc_site_ids is True, the site indices of the

        sites in the component as a tuple.


    * ”molecule_graph”: If inc_molecule_graph is True, the site a

        MoleculeGraph object for zero-dimensional components.




* **Return type**

    (list of dict)



### zero_d_graph_to_molecule_graph(bonded_structure, graph)
Converts a zero-dimensional networkx Graph object into a MoleculeGraph.

Implements a similar breadth-first search to that in
calculate_dimensionality_of_site().


* **Parameters**


    * **bonded_structure** (*StructureGraph*) – A structure with bonds, represented
    as a pymatgen structure graph. For example, generated using the
    CrystalNN.get_bonded_structure() method.


    * **graph** (*nx.Graph*) – A networkx Graph object for the component of
    interest.



* **Returns**

    A MoleculeGraph object of the component.



* **Return type**

    (MoleculeGraph)


## pymatgen.analysis.disorder module

This module provides various methods to analyze order/disorder in materials.


### get_warren_cowley_parameters(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), r: float, dr: float)
Warren-Crowley parameters.


* **Parameters**


    * **structure** – Pymatgen Structure.


    * **r** – Radius


    * **dr** – Shell width



* **Returns**

    -1.0, …}



* **Return type**

    Warren-Crowley parameters in the form of a dict, e.g., {(Element Mo, Element W)


## pymatgen.analysis.energy_models module

This module implements a EnergyModel abstract class and some basic
implementations. Basically, an EnergyModel is any model that returns an
“energy” for any given structure.


### _class_ EnergyModel()
Bases: `MSONable`

Abstract structure filter class.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict representation.



* **Returns**

    EnergyModel



#### _abstract_ get_energy(structure)

* **Parameters**

    **structure** – Structure



* **Returns**

    Energy value



### _class_ EwaldElectrostaticModel(real_space_cut=None, recip_space_cut=None, eta=None, acc_factor=8.0)
Bases: `EnergyModel`

Wrapper around EwaldSum to calculate the electrostatic energy.

Initializes the model. Args have the same definitions as in
pymatgen.analysis.ewald.EwaldSummation.


* **Parameters**


    * **real_space_cut** (*float*) – Real space cutoff radius dictating how
    many terms are used in the real space sum. Defaults to None,
    which means determine automatically using the formula given
    in gulp 3.1 documentation.


    * **recip_space_cut** (*float*) – Reciprocal space cutoff radius.
    Defaults to None, which means determine automatically using
    the formula given in gulp 3.1 documentation.


    * **eta** (*float*) – Screening parameter. Defaults to None, which means
    determine automatically.


    * **acc_factor** (*float*) – No. of significant figures each sum is
    converged to.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict


#### get_energy(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** – Structure



* **Returns**

    Energy value



### _class_ IsingModel(j, max_radius)
Bases: `EnergyModel`

A very simple Ising model, with r^2 decay.


* **Parameters**


    * **j** (*float*) – The interaction parameter. E = J \* spin1 \* spin2.


    * **radius** (*float*) – max_radius for the interaction.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict


#### get_energy(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** – Structure



* **Returns**

    Energy value



### _class_ NsitesModel()
Bases: `EnergyModel`

Sets the energy to the number of sites. More sites => higher “energy”.
Used to rank structures from smallest number of sites to largest number
of sites after enumeration.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict


#### get_energy(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** – Structure



* **Returns**

    Energy value



### _class_ SymmetryModel(symprec: float = 0.1, angle_tolerance=5)
Bases: `EnergyModel`

Sets the energy to the negative of the spacegroup number. Higher symmetry =>
lower “energy”.

Args have same meaning as in pymatgen.symmetry.SpacegroupAnalyzer.


* **Parameters**


    * **symprec** (*float*) – Symmetry tolerance. Defaults to 0.1.


    * **angle_tolerance** (*float*) – Tolerance for angles. Defaults to 5 degrees.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
MSONable dict


#### get_energy(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** – Structure



* **Returns**

    Energy value


## pymatgen.analysis.eos module

This module implements various equation of states.

Note: Most of the code were initially adapted from ASE and deltafactor by
@gmatteo but has since undergone major refactoring.


### _class_ Birch(volumes, energies)
Bases: `EOSBase`

Birch EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _func(volume, params)
From Intermetallic compounds: Principles and Practice, Vol. I:
Principles Chapter 9 pages 195-210 by M. Mehl. B. Klein,
D. Papaconstantopoulos.
case where n=0.


### _class_ BirchMurnaghan(volumes, energies)
Bases: `EOSBase`

BirchMurnaghan EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _func(volume, params)
BirchMurnaghan equation from PRB 70, 224107.


### _class_ DeltaFactor(volumes, energies)
Bases: `PolynomialEOS`

Fitting a polynomial EOS using delta factor.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _func(volume, params)
The equation of state function. This must be implemented by all classes
that derive from this abstract class.


* **Parameters**

    **volume** (*float/numpy.array*) – params (list/tuple): values for the parameters other than the

        volume used by the eos.




#### _set_params()
Overridden to account for the fact the fit with volume\*\*(2/3) instead
of volume.


#### fit(order=3)
Overridden since this eos works with volume\*\*(2/3) instead of volume.


### _class_ EOS(eos_name='murnaghan')
Bases: `object`

Convenient wrapper. Retained in its original state to ensure backward
compatibility.

Fit equation of state for bulk systems.

The following equations are supported:

```default
murnaghan: PRB 28, 5480 (1983)

birch: Intermetallic compounds: Principles and Practice, Vol I:
    Principles. pages 195-210

birch_murnaghan: PRB 70, 224107

pourier_tarantola: PRB 70, 224107

vinet: PRB 70, 224107

deltafactor

numerical_eos: 10.1103/PhysRevB.90.174107.
```

Usage:

```default
eos = EOS(eos_name='murnaghan')
eos_fit = eos.fit(volumes, energies)
eos_fit.plot()
```


* **Parameters**

    **eos_name** (*str*) – Type of EOS to fit.



#### MODELS(_ = {'birch': <class 'pymatgen.analysis.eos.Birch'>, 'birch_murnaghan': <class 'pymatgen.analysis.eos.BirchMurnaghan'>, 'deltafactor': <class 'pymatgen.analysis.eos.DeltaFactor'>, 'murnaghan': <class 'pymatgen.analysis.eos.Murnaghan'>, 'numerical_eos': <class 'pymatgen.analysis.eos.NumericalEOS'>, 'pourier_tarantola': <class 'pymatgen.analysis.eos.PourierTarantola'>, 'vinet': <class 'pymatgen.analysis.eos.Vinet'>_ )

#### fit(volumes, energies)
Fit energies as function of volumes.


* **Parameters**


    * **volumes** (*list/np.array*) –


    * **energies** (*list/np.array*) –



* **Returns**

    EOSBase object



* **Return type**

    EOSBase



### _class_ EOSBase(volumes, energies)
Bases: `object`

Abstract class that must be subclassed by all equation of state
implementations.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract_ _func(volume, params)
The equation of state function. This must be implemented by all classes
that derive from this abstract class.


* **Parameters**

    **volume** (*float/numpy.array*) – params (list/tuple): values for the parameters other than the

        volume used by the eos.




#### _initial_guess()
Quadratic fit to get an initial guess for the parameters.


* **Returns**

    (e0, b0, b1, v0)



* **Return type**

    tuple



#### _property_ b0()
Returns the bulk modulus.
Note: the units for the bulk modulus: unit of energy/unit of volume^3.


#### _property_ b0_GPa()
Returns the bulk modulus in GPa.
Note: This assumes that the energy and volumes are in eV and Ang^3

> respectively.


#### _property_ b1()
Returns the derivative of bulk modulus wrt pressure(dimensionless).


#### _property_ e0()
Returns the min energy.


#### fit()
Do the fitting. Does least square fitting. If you want to use custom
fitting, must override this.


#### func(volume)
The equation of state function with the parameters other than volume set
to the ones obtained from fitting.


* **Parameters**

    **volume** (*list/numpy.array*) –



* **Returns**

    numpy.array



#### plot(width=8, height=None, ax: plt.Axes = None, dpi=None, \*\*kwargs)
Plot the equation of state.


* **Parameters**


    * **width** (*float*) – Width of plot in inches. Defaults to 8in.


    * **height** (*float*) – Height of plot in inches. Defaults to width \*
    golden ratio.


    * **ax** (*plt.Axes*) – If supplied, changes will be made to the existing Axes.
    Otherwise, new Axes will be created.


    * **dpi** –


    * **kwargs** (*dict*) – additional args fed to pyplot.plot.
    supported keys: style, color, text, label



* **Returns**

    The matplotlib axes.



* **Return type**

    plt.Axes



#### plot_ax(ax: plt.Axes = None, fontsize=12, \*\*kwargs)
Plot the equation of state on axis ax.


* **Parameters**


    * **ax** – matplotlib Axes or None if a new figure should be created.


    * **fontsize** – Legend fontsize.


    * **color** (*str*) – plot color.


    * **label** (*str*) – Plot label


    * **text** (*str*) – Legend text (options)



* **Returns**

    matplotlib figure.



* **Return type**

    plt.Figure


Keyword arguments controlling the display of the figure:

| kwargs

 | Meaning

 |
| ------ | ------- |
| title

  | Title of the plot (Default: None).

 |
| show

   | True to show the figure (default: True).

 |
| savefig

 | “abc.png” or “abc.eps” to save the figure to a file.

 |
| size_kwargs

 | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

 |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                        |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### _property_ results()
Returns a summary dict.


* **Returns**

    dict



#### _property_ v0()
Returns the minimum or the reference volume in Ang^3.


### _exception_ EOSError()
Bases: `Exception`

Error class for EOS fitting.


### _class_ Murnaghan(volumes, energies)
Bases: `EOSBase`

Murnaghan EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _func(volume, params)
From PRB 28,5480 (1983).


### _class_ NumericalEOS(volumes, energies)
Bases: `PolynomialEOS`

A numerical EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### fit(min_ndata_factor=3, max_poly_order_factor=5, min_poly_order=2)
Fit the input data to the ‘numerical eos’, the equation of state employed
in the quasiharmonic Debye model described in the paper:
10.1103/PhysRevB.90.174107.

credits: Cormac Toher


* **Parameters**


    * **min_ndata_factor** (*int*) – parameter that controls the minimum number
    of data points that will be used for fitting.
    minimum number of data points =

    > total data points-2\*min_ndata_factor



    * **max_poly_order_factor** (*int*) – parameter that limits the max order
    of the polynomial used for fitting.
    max_poly_order = number of data points used for fitting -

    > max_poly_order_factor



    * **min_poly_order** (*int*) – minimum order of the polynomial to be
    considered for fitting.



### _class_ PolynomialEOS(volumes, energies)
Bases: `EOSBase`

Derives from EOSBase. Polynomial based equations of states must subclass
this.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _func(volume, params)
The equation of state function. This must be implemented by all classes
that derive from this abstract class.


* **Parameters**

    **volume** (*float/numpy.array*) – params (list/tuple): values for the parameters other than the

        volume used by the eos.




#### _set_params()
Use the fit polynomial to compute the parameter e0, b0, b1 and v0
and set to the _params attribute.


#### fit(order)
Do polynomial fitting and set the parameters. Uses numpy polyfit.


* **Parameters**

    **order** (*int*) – order of the fit polynomial



### _class_ PourierTarantola(volumes, energies)
Bases: `EOSBase`

PourierTarantola EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _func(volume, params)
Pourier-Tarantola equation from PRB 70, 224107.


### _class_ Vinet(volumes, energies)
Bases: `EOSBase`

Vinet EOS.


* **Parameters**


    * **volumes** (*list/numpy.array*) – volumes in Ang^3


    * **energies** (*list/numpy.array*) – energy in eV.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _func(volume, params)
Vinet equation from PRB 70, 224107.

## pymatgen.analysis.ewald module

This module provides classes for calculating the Ewald sum of a structure.


### _class_ EwaldMinimizer(matrix, m_list, num_to_return=1, algo=0)
Bases: `object`

This class determines the manipulations that will minimize an Ewald matrix,
given a list of possible manipulations. This class does not perform the
manipulations on a structure, but will return the list of manipulations
that should be done on one to produce the minimal structure. It returns the
manipulations for the n lowest energy orderings. This class should be used
to perform fractional species substitution or fractional species removal to
produce a new structure. These manipulations create large numbers of
candidate structures, and this class can be used to pick out those with the
lowest Ewald sum.

An alternative (possibly more intuitive) interface to this class is the
order disordered structure transformation.

Author - Will Richards


* **Parameters**


    * **matrix** – A matrix of the Ewald sum interaction energies. This is stored
    in the class as a diagonally symmetric array and so
    self._matrix will not be the same as the input matrix.


    * **m_list** – list of manipulations. each item is of the form
    (multiplication fraction, number_of_indices, indices, species)
    These are sorted such that the first manipulation contains the
    most permutations. this is actually evaluated last in the
    recursion since I’m using pop.


    * **num_to_return** – The minimizer will find the number_returned lowest
    energy structures. This is likely to return a number of duplicate
    structures so it may be necessary to overestimate and then
    remove the duplicates later. (duplicate checking in this
    process is extremely expensive).



#### ALGO_BEST_FIRST(_ = _ )
Slowly increases the speed (with the cost of decreasing
accuracy) as the minimizer runs. Attempts to limit the run time to
approximately 30 minutes.


* **Type**

    ALGO_TIME_LIMIT



#### ALGO_COMPLETE(_ = _ )

#### ALGO_FAST(_ = _ )

#### ALGO_TIME_LIMIT(_ = _ )

#### _recurse(matrix, m_list, indices, output_m_list=None)
This method recursively finds the minimal permutations using a binary
tree search strategy.


* **Parameters**


    * **matrix** – The current matrix (with some permutations already
    performed).


    * **m_list** – The list of permutations still to be performed


    * **indices** – Set of indices which haven’t had a permutation
    performed on them.



#### add_m_list(matrix_sum, m_list)
This adds an m_list to the output_lists and updates the current
minimum if the list is full.


#### best_case(matrix, m_list, indices_left)
Computes a best case given a matrix and manipulation list.


* **Parameters**


    * **matrix** – the current matrix (with some permutations already
    performed)


    * **m_list** – [(multiplication fraction, number_of_indices, indices,
    species)] describing the manipulation


    * **indices** – Set of indices which haven’t had a permutation
    performed on them.



#### _property_ best_m_list()
Best m_list found.


* **Type**

    Returns



#### _classmethod_ get_next_index(matrix, manipulation, indices_left)
Returns an index that should have the most negative effect on the
matrix sum.


#### minimize_matrix()
This method finds and returns the permutations that produce the lowest
Ewald sum calls recursive function to iterate through permutations.


#### _property_ minimized_sum()
Minimized sum.


* **Type**

    Returns



#### _property_ output_lists()
output lists.


* **Type**

    Returns



### _class_ EwaldSummation(structure, real_space_cut=None, recip_space_cut=None, eta=None, acc_factor=12.0, w=0.7071067811865475, compute_forces=False)
Bases: `MSONable`

Calculates the electrostatic energy of a periodic array of charges using
the Ewald technique.

Ref:

    Ewald summation techniques in perspective: a survey
    Abdulnour Y. Toukmaji and John A. Board Jr.
    DOI: 10.1016/0010-4655(96)00016-1
    URL: [http://www.ee.duke.edu/~ayt/ewaldpaper/ewaldpaper.html](http://www.ee.duke.edu/~ayt/ewaldpaper/ewaldpaper.html)

This matrix can be used to do fast calculations of Ewald sums after species
removal.

E = E_recip + E_real + E_point

Atomic units used in the code, then converted to eV.

Initializes and calculates the Ewald sum. Default convergence
parameters have been specified, but you can override them if you wish.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure that must have proper
    Species on all sites, i.e. Element with oxidation state. Use
    Structure.add_oxidation_state… for example.


    * **real_space_cut** (*float*) – Real space cutoff radius dictating how
    many terms are used in the real space sum. Defaults to None,
    which means determine automatically using the formula given
    in gulp 3.1 documentation.


    * **recip_space_cut** (*float*) – Reciprocal space cutoff radius.
    Defaults to None, which means determine automatically using
    the formula given in gulp 3.1 documentation.


    * **eta** (*float*) – The screening parameter. Defaults to None, which means
    determine automatically.


    * **acc_factor** (*float*) – No. of significant figures each sum is
    converged to.


    * **w** (*float*) – Weight parameter, w, has been included that represents
    the relative computational expense of calculating a term in
    real and reciprocal space. Default of 0.7 reproduces result
    similar to GULP 4.2. This has little effect on the total
    energy, but may influence speed of computation in large
    systems. Note that this parameter is used only when the
    cutoffs are set to None.


    * **compute_forces** (*bool*) – Whether to compute forces. False by
    default since it is usually not needed.



#### CONV_FACT(_ = 14.3996454784256_ )

#### _calc_ewald_terms()
Calculates and sets all Ewald terms (point, real and reciprocal).


#### _calc_real_and_point()
Determines the self energy -(eta/pi)\*\*(1/2) \* sum_{i=1}^{N} q_i\*\*2.


#### _calc_recip()
Perform the reciprocal space summation. Calculates the quantity
E_recip = 1/(2PiV) sum_{G < Gmax} exp(-(G.G/4/eta))/(G.G) S(G)S(-G)
where
S(G) = sum_{k=1,N} q_k exp(-i G.r_k)
S(G)S(-G) =

```
|
```

S(G)|\*\*2.

This method is heavily vectorized to utilize numpy’s C backend for
speed.


#### as_dict(verbosity: int = 0)
Json-serialization dict representation of EwaldSummation.


* **Parameters**

    **verbosity** (*int*) – Verbosity level. Default of 0 only includes the
    matrix representation. Set to 1 for more details.



#### compute_partial_energy(removed_indices)
Gives total Ewald energy for certain sites being removed, i.e. zeroed
out.


#### compute_sub_structure(sub_structure, tol: float = 0.001)
Gives total Ewald energy for an sub structure in the same
lattice. The sub_structure must be a subset of the original
structure, with possible different charges.


* **Parameters**


    * **substructure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Substructure to compute Ewald sum for.


    * **tol** (*float*) – Tolerance for site matching in fractional coordinates.



* **Returns**

    Ewald sum of substructure.



#### _property_ eta()
eta value used in Ewald summation.


* **Type**

    Returns



#### _property_ forces()
The forces on each site as a Nx3 matrix. Each row corresponds to a
site.


#### _classmethod_ from_dict(d: dict[str, Any], fmt: str | None = None, \*\*kwargs)
Create an EwaldSummation instance from JSON-serialized dictionary.


* **Parameters**


    * **d** (*dict*) – Dictionary representation


    * **fmt** (*str**, **optional*) – Unused. Defaults to None.



* **Returns**

    class instance



* **Return type**

    EwaldSummation



#### get_site_energy(site_index)
Compute the energy for a single site in the structure.


* **Parameters**

    **site_index** (*int*) – Index of site


Returns:
(float) - Energy of that site


#### _property_ point_energy()
The point energy.


#### _property_ point_energy_matrix()
The point space matrix. A diagonal matrix with the point terms for each
site in the diagonal elements.


#### _property_ real_space_energy()
The real space energy.


#### _property_ real_space_energy_matrix()
The real space energy matrix. Each matrix element (i, j) corresponds to
the interaction energy between site i and site j in real space.


#### _property_ reciprocal_space_energy()
The reciprocal space energy.


#### _property_ reciprocal_space_energy_matrix()
The reciprocal space energy matrix. Each matrix element (i, j)
corresponds to the interaction energy between site i and site j in
reciprocal space.


#### _property_ total_energy()
The total energy.


#### _property_ total_energy_matrix()
The total energy matrix. Each matrix element (i, j) corresponds to the
total interaction energy between site i and site j.

Note that this does not include the charged-cell energy, which is only important
when the simulation cell is not charge balanced.


### compute_average_oxidation_state(site)
Calculates the average oxidation state of a site.


* **Parameters**

    **site** – Site to compute average oxidation state



* **Returns**

    Average oxidation state of site.


## pymatgen.analysis.excitation module

This module defines an excitation spectrum class.


### _class_ ExcitationSpectrum(x, y)
Bases: [`Spectrum`](pymatgen.core.md#pymatgen.core.spectrum.Spectrum)

Basic excitation spectrum object.


#### x()
The sequence of energies.


* **Type**

    Sequence[float]



#### y()
The sequence of mu(E).


* **Type**

    Sequence[float]



* **Parameters**


    * **x** – A sequence of x-ray energies in eV


    * **y** – A sequence of intensity values.



#### XLABEL(_ = 'Energy (eV)_ )

#### YLABEL(_ = 'Intensity_ )
## pymatgen.analysis.fragmenter module

Perform fragmentation of molecules.


### _class_ Fragmenter(molecule, edges=None, depth=1, open_rings=False, use_metal_edge_extender=False, opt_steps=10000, prev_unique_frag_dict=None, assume_previous_thoroughness=True)
Bases: `MSONable`

Molecule fragmenter class.

Standard constructor for molecule fragmentation.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – The molecule to fragment.


    * **edges** (*list*) – List of index pairs that define graph edges, aka molecule bonds. If not set,
    edges will be determined with OpenBabel. Defaults to None.


    * **depth** (*int*) – The number of levels of iterative fragmentation to perform, where each level
    will include fragments obtained by breaking one bond of a fragment one level up.
    Defaults to 1. However, if set to 0, instead all possible fragments are generated
    using an alternative, non-iterative scheme.


    * **open_rings** (*bool*) – Whether or not to open any rings encountered during fragmentation.
    Defaults to False. If true, any bond that fails to yield disconnected graphs when
    broken is instead removed and the entire structure is optimized with OpenBabel in
    order to obtain a good initial guess for an opened geometry that can then be put
    back into QChem to be optimized without the ring just reforming.


    * **use_metal_edge_extender** (*bool*) – Whether or not to attempt to add additional edges from
    O, N, F, or Cl to any Li or Mg atoms present that OpenBabel may have missed. Defaults
    to False. Most important for ionic bonding. Note that additional metal edges may yield
    new “rings” (e.g. -C-O-Li-O- in LiEC) that will not play nicely with ring opening.


    * **opt_steps** (*int*) – Number of optimization steps when opening rings. Defaults to 10000.


    * **prev_unique_frag_dict** (*dict*) – A dictionary of previously identified unique fragments.
    Defaults to None. Typically only used when trying to find the set of unique fragments
    that come from multiple molecules.


    * **assume_previous_thoroughness** (*bool*) – Whether or not to assume that a molecule / fragment
    provided in prev_unique_frag_dict has all of its unique subfragments also provided in
    prev_unique_frag_dict. Defaults to True. This is an essential optimization when trying
    to find the set of unique fragments that come from multiple molecules if all of those
    molecules are being fully iteratively fragmented. However, if you’re passing a
    prev_unique_frag_dict which includes a molecule and its fragments that were generated
    at insufficient depth to find all possible subfragments to a fragmentation calculation
    of a different molecule that you aim to find all possible subfragments of and which has
    common subfragments with the previous molecule, this optimization will cause you to
    miss some unique subfragments.



#### _fragment_one_level(old_frag_dict)
Perform one step of iterative fragmentation on a list of molecule graphs. Loop through the graphs,
then loop through each graph’s edges and attempt to remove that edge in order to obtain two
disconnected subgraphs, aka two new fragments. If successful, check to see if the new fragments
are already present in self.unique_fragments, and append them if not. If unsuccessful, we know
that edge belongs to a ring. If we are opening rings, do so with that bond, and then again
check if the resulting fragment is present in self.unique_fragments and add it if it is not.


#### _open_all_rings()
Having already generated all unique fragments that did not require ring opening,
now we want to also obtain fragments that do require opening. We achieve this by
looping through all unique fragments and opening each bond present in any ring
we find. We also temporarily add the principle molecule graph to self.unique_fragments
so that its rings are opened as well.


### open_ring(mol_graph, bond, opt_steps)
Function to actually open a ring using OpenBabel’s local opt. Given a molecule
graph and a bond, convert the molecule graph into an OpenBabel molecule, remove
the given bond, perform the local opt with the number of steps determined by
self.steps, and then convert the resulting structure back into a molecule graph
to be returned.

## pymatgen.analysis.functional_groups module

Determine functional groups present in a Molecule.


### _class_ FunctionalGroupExtractor(molecule, optimize=False)
Bases: `object`

This class is used to algorithmically parse a molecule (represented by an
instance of pymatgen.analysis.graphs.MoleculeGraph) and determine arbitrary
functional groups.

Instantiation method for FunctionalGroupExtractor.


* **Parameters**


    * **molecule** – Either a filename, a pymatgen.core.structure.Molecule
    object, or a pymatgen.analysis.graphs.MoleculeGraph object.


    * **optimize** – Default False. If True, then the input molecule will be
    modified, adding Hydrogens, performing a simple conformer search,
    etc.



#### categorize_functional_groups(groups)
Determine classes of functional groups present in a set.


* **Parameters**

    **groups** – Set of functional groups.



* **Returns**

    dict containing representations of the groups, the indices of
    where the group occurs in the MoleculeGraph, and how many of each
    type of group there is.



#### get_all_functional_groups(elements=None, func_groups=None, catch_basic=True)
Identify all functional groups (or all within a certain subset) in the
molecule, combining the methods described above.


* **Parameters**


    * **elements** – List of elements that will qualify a carbon as special
    (if only certain functional groups are of interest).
    Default None.


    * **func_groups** – List of strs representing the functional groups of
    interest. Default to None, meaning that all of the functional groups
    defined in this function will be sought.


    * **catch_basic** – bool. If True, use get_basic_functional_groups and
    other methods



* **Returns**

    list of sets of ints, representing groups of connected atoms



#### get_basic_functional_groups(func_groups=None)
Identify functional groups that cannot be identified by the Ertl method
of get_special_carbon and get_heteroatoms, such as benzene rings, methyl
groups, and ethyl groups.

TODO: Think of other functional groups that are important enough to be
added (ex: do we need ethyl, butyl, propyl?)


* **Parameters**

    **func_groups** – List of strs representing the functional groups of
    interest. Default to None, meaning that all of the functional groups
    defined in this function will be sought.



* **Returns**

    list of sets of ints, representing groups of connected atoms



#### get_heteroatoms(elements=None)
Identify non-H, non-C atoms in the MoleculeGraph, returning a list of
their node indices.


* **Parameters**

    **elements** – List of elements to identify (if only certain
    functional groups are of interest).



* **Returns**

    set of ints representing node indices



#### get_special_carbon(elements=None)
Identify Carbon atoms in the MoleculeGraph that fit the characteristics
defined Ertl (2017), returning a list of their node indices.

The conditions for marking carbon atoms are (quoted from Ertl):

    “- atoms connected by non-aromatic double or triple bond to any
    heteroatom
    - atoms in nonaromatic carbon-carbon double or triple bonds
    - acetal carbons, i.e. sp3 carbons connected to two or more oxygens,
    nitrogens or sulfurs; these O, N or S atoms must have only single bonds
    - all atoms in oxirane, aziridine and thiirane rings”


* **Parameters**

    **elements** – List of elements that will qualify a carbon as special
    (if only certain functional groups are of interest).
    Default None.



* **Returns**

    set of ints representing node indices



#### link_marked_atoms(atoms)
Take a list of marked “interesting” atoms (heteroatoms, special carbons)
and attempt to connect them, returning a list of disjoint groups of
special atoms (and their connected hydrogens).


* **Parameters**

    **atoms** – set of marked “interesting” atoms, presumably identified
    using other functions in this class.



* **Returns**

    list of sets of ints, representing groups of connected atoms


## pymatgen.analysis.graphs module

Module for graph representations of crystals and molecules.


### _class_ ConnectedSite(site, jimage, index, weight, dist)
Bases: `tuple`

Create new instance of ConnectedSite(site, jimage, index, weight, dist)


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('site', 'jimage', 'index', 'weight', 'dist'_ )

#### _classmethod_ _make(iterable)
Make a new ConnectedSite object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new ConnectedSite object replacing specified fields with new values


#### dist()
Alias for field number 4


#### index()
Alias for field number 2


#### jimage()
Alias for field number 1


#### site()
Alias for field number 0


#### weight()
Alias for field number 3


### _exception_ MolGraphSplitError()
Bases: `Exception`

Raised when a molecule graph is failed to split into two disconnected
subgraphs.


### _class_ MoleculeGraph(molecule, graph_data=None)
Bases: `MSONable`

This is a class for annotating a Molecule with
bond information, stored in the form of a graph. A “bond” does
not necessarily have to be a chemical bond, but can store any
kind of information that connects two Sites.

If constructing this class manually, use the with_empty_graph
method or with_local_env_strategy method (using an algorithm
provided by the local_env module, such as O’Keeffe).

This class that contains connection information:
relationships between sites represented by a Graph structure,
and an associated structure object.

This class uses the NetworkX package to store and operate
on the graph itself, but contains a lot of helper methods
to make associating a graph with a given molecule easier.

Use cases for this include storing bonding information,
NMR J-couplings, Heisenberg exchange parameters, etc.


* **Parameters**


    * **molecule** – Molecule object


    * **graph_data** – dict containing graph information in
    dict format (not intended to be constructed manually,
    see as_dict method for format)



#### _classmethod_ _edges_to_string(g)

#### add_edge(from_index, to_index, weight=None, warn_duplicates=True, edge_properties=None)
Add edge to graph.

Since physically a ‘bond’ (or other connection
between sites) doesn’t have a direction, from_index,
from_jimage can be swapped with to_index, to_jimage.

However, images will always be shifted so that
from_index < to_index and from_jimage becomes (0, 0, 0).


* **Parameters**


    * **from_index** – index of site connecting from


    * **to_index** – index of site connecting to


    * **(****float****)** (*weight*) – e.g. bond length


    * **(****bool****)** (*warn_duplicates*) – if True, will warn if
    trying to add duplicate edges (duplicate edges will not
    be added in either case)


    * **(****dict****)** (*edge_properties*) – any other information to
    store on graph edges, similar to Structure’s site_properties



#### alter_edge(from_index, to_index, new_weight=None, new_edge_properties=None)
Alters either the weight or the edge_properties of
an edge in the MoleculeGraph.


* **Parameters**


    * **from_index** – int


    * **to_index** – int


    * **new_weight** – alter_edge does not require
    that weight be altered. As such, by default, this
    is None. If weight is to be changed, it should be a
    float.


    * **new_edge_properties** – alter_edge does not require
    that edge_properties be altered. As such, by default,
    this is None. If any edge properties are to be changed,
    it should be a dictionary of edge properties to be changed.



#### as_dict()
As in pymatgen.core.Molecule except
with using to_dict_of_dicts from NetworkX
to store graph information.


#### break_edge(from_index, to_index, allow_reverse=False)
Remove an edge from the MoleculeGraph.


* **Parameters**


    * **from_index** – int


    * **to_index** – int


    * **allow_reverse** – If allow_reverse is True, then break_edge will
    attempt to break both (from_index, to_index) and, failing that,
    will attempt to break (to_index, from_index).



#### build_unique_fragments()
Find all possible fragment combinations of the MoleculeGraphs (in other
words, all connected induced subgraphs).


#### diff(other, strict=True)
Compares two MoleculeGraphs. Returns dict with
keys ‘self’, ‘other’, ‘both’ with edges that are
present in only one MoleculeGraph (‘self’ and
‘other’), and edges that are present in both.

The Jaccard distance is a simple measure of the
dissimilarity between two MoleculeGraphs (ignoring
edge weights), and is defined by 1 - (size of the
intersection / size of the union) of the sets of
edges. This is returned with key ‘dist’.

Important note: all node indices are in terms
of the MoleculeGraph this method is called
from, not the ‘other’ MoleculeGraph: there
is no guarantee the node indices will be the
same if the underlying Molecules are ordered
differently.


* **Parameters**


    * **other** – MoleculeGraph


    * **strict** – if False, will compare bonds
    from different Molecules, with node indices
    replaced by Species strings, will not count
    number of occurrences of bonds



#### draw_graph_to_file(filename='graph', diff=None, hide_unconnected_nodes=False, hide_image_edges=True, edge_colors=False, node_labels=False, weight_labels=False, image_labels=False, color_scheme='VESTA', keep_dot=False, algo='fdp')
Draws graph using GraphViz.

The networkx graph object itself can also be drawn
with networkx’s in-built graph drawing methods, but
note that this might give misleading results for
multigraphs (edges are super-imposed on each other).

If visualization is difficult to interpret,
hide_image_edges can help, especially in larger
graphs.


* **Parameters**


    * **filename** – filename to output, will detect filetype
    from extension (any graphviz filetype supported, such as
    pdf or png)


    * **(****StructureGraph****)** (*diff*) – an additional graph to
    compare with, will color edges red that do not exist in diff
    and edges green that are in diff graph but not in the
    reference graph


    * **hide_unconnected_nodes** – if True, hide unconnected
    nodes


    * **hide_image_edges** – if True, do not draw edges that
    go through periodic boundaries


    * **(****bool****)** (*keep_dot*) – if True, use node colors to
    color edges


    * **(****bool****)** – if True, label nodes with
    species and site index


    * **(****bool****)** – if True, label edges with
    weights


    * **(****bool****)** – if True, label edges with
    their periodic images (usually only used for debugging,
    edges to periodic images always appear as dashed lines)


    * **(****str****)** (*color_scheme*) – “VESTA” or “JMOL”


    * **(****bool****)** – keep GraphViz .dot file for later
    visualization


    * **algo** – any graphviz algo, “neato” (for simple graphs)
    or “fdp” (for more crowded graphs) usually give good outputs



#### _property_ edge_weight_name()
Name of the edge weight property of graph


#### _property_ edge_weight_unit()
Units of the edge weight property of graph


#### find_rings(including=None)
Find ring structures in the MoleculeGraph.


* **Parameters**


    * **including** (*list**[**int**]*) – list of site indices. If including is not None, then find_rings


    * **default** (*will only return those rings including the specified sites. By*) –


    * **parameter** (*this*) –


    * **None** (*is*) –


    * **returned.** (*and all rings will be*) –



* **Returns**

    Each entry will be a ring (cycle, in graph theory terms)

        including the index found in the Molecule. If there is no cycle including an index, the
        value will be an empty list.




* **Return type**

    list[list[tuple[int, int]]]



#### _classmethod_ from_dict(dct)
As in pymatgen.core.Molecule except
restoring graphs using from_dict_of_dicts
from NetworkX to restore graph information.


#### get_connected_sites(n)
Returns a named tuple of neighbors of site n:
periodic_site, jimage, index, weight.
Index is the index of the corresponding site
in the original structure, weight can be
None if not defined.
:param n: index of Site in Molecule
:param jimage: lattice vector of site


* **Returns**

    list of ConnectedSite tuples,
    sorted by closest first.



#### get_coordination_of_site(n)
Returns the number of neighbors of site n.
In graph terms, simply returns degree
of node corresponding to site n.
:param n: index of site
:return (int):


#### get_disconnected_fragments(return_index_map: bool = False)
Determine if the MoleculeGraph is connected. If it is not, separate the
MoleculeGraph into different MoleculeGraphs, where each resulting MoleculeGraph is
a disconnected subgraph of the original. Currently, this function naively assigns
the charge of the total molecule to a single submolecule. A later effort will be
to actually accurately assign charge.


* **Parameters**

    **return_index_map** (*bool*) – If True, return a dictionary that maps the
    new indices to the original indices. Defaults to False.


NOTE: This function does not modify the original MoleculeGraph. It creates a copy,
modifies that, and returns two or more new MoleculeGraph objects.


* **Returns**

    Each MoleculeGraph is a disconnected subgraph of the original MoleculeGraph.



* **Return type**

    list[MoleculeGraph]



#### insert_node(idx, species, coords, validate_proximity=False, site_properties=None, edges=None)
A wrapper around Molecule.insert(), which also incorporates the new
site into the MoleculeGraph.


* **Parameters**


    * **idx** – Index at which to insert the new site


    * **species** – Species for the new site


    * **coords** – 3x1 array representing coordinates of the new site


    * **validate_proximity** – For Molecule.insert(); if True (default
    False), distance will be checked to ensure that site can be safely
    added.


    * **site_properties** – Site properties for Molecule


    * **edges** – List of dicts representing edges to be added to the
    MoleculeGraph. These edges must include the index of the new site i,
    and all indices used for these edges should reflect the
    MoleculeGraph AFTER the insertion, NOT before. Each dict should at
    least have a “to_index” and “from_index” key, and can also have a
    “weight” and a “properties” key.



#### isomorphic_to(other)
Checks if the graphs of two MoleculeGraphs are isomorphic to one
another. In order to prevent problems with misdirected edges, both
graphs are converted into undirected nx.Graph objects.


* **Parameters**

    **other** – MoleculeGraph object to be compared.



* **Returns**

    bool



#### _property_ name()
Name of graph


#### remove_nodes(indices)
A wrapper for Molecule.remove_sites().


* **Parameters**

    **indices** – list of indices in the current Molecule (and graph) to
    be removed.



#### replace_group(index, func_grp, strategy, bond_order=1, graph_dict=None, strategy_params=None)
Builds off of Molecule.substitute and MoleculeGraph.substitute_group
to replace a functional group in self.molecule with a functional group.
This method also amends self.graph to incorporate the new functional
group.

TODO: Figure out how to replace into a ring structure.


* **Parameters**


    * **index** – Index of atom to substitute.


    * **func_grp** – Substituent molecule. There are three options:


        1. Providing an actual molecule as the input. The first atom
    must be a DummySpecies X, indicating the position of
    nearest neighbor. The second atom must be the next
    nearest atom. For example, for a methyl group
    substitution, func_grp should be X-CH3, where X is the
    first site and C is the second site. What the code will
    do is to remove the index site, and connect the nearest
    neighbor to the C atom in CH3. The X-C bond indicates the
    directionality to connect the atoms.


        2. A string name. The molecule will be obtained from the
    relevant template in func_groups.json.


        3. A MoleculeGraph object.



    * **strategy** – Class from pymatgen.analysis.local_env.


    * **bond_order** – A specified bond order to calculate the bond
    length between the attached functional group and the nearest
    neighbor site. Defaults to 1.


    * **graph_dict** – Dictionary representing the bonds of the functional
    group (format: {(u, v): props}, where props is a dictionary of
    properties, including weight. If None, then the algorithm
    will attempt to automatically determine bonds using one of
    a list of strategies defined in pymatgen.analysis.local_env.


    * **strategy_params** – dictionary of keyword arguments for strategy.
    If None, default parameters will be used.



#### set_node_attributes()
Replicates molecule site properties (specie, coords, etc.) in the
MoleculeGraph.


#### sort(key: Callable[[[Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule)], float] | None = None, reverse: bool = False)
Same as Molecule.sort(). Also remaps nodes in graph.


* **Parameters**


    * **key** (*callable**, **optional*) – Sort key. Defaults to None.


    * **reverse** (*bool**, **optional*) – Reverse sort order. Defaults to False.



#### split_molecule_subgraphs(bonds, allow_reverse=False, alterations=None)
Split MoleculeGraph into two or more MoleculeGraphs by
breaking a set of bonds. This function uses
MoleculeGraph.break_edge repeatedly to create
disjoint graphs (two or more separate molecules).
This function does not only alter the graph
information, but also changes the underlying
Molecules.
If the bonds parameter does not include sufficient
bonds to separate two molecule fragments, then this
function will fail.
Currently, this function naively assigns the charge
of the total molecule to a single submolecule. A
later effort will be to actually accurately assign
charge.
NOTE: This function does not modify the original
MoleculeGraph. It creates a copy, modifies that, and
returns two or more new MoleculeGraph objects.
:param bonds: list of tuples (from_index, to_index)

> representing bonds to be broken to split the MoleculeGraph.


* **Parameters**


    * **alterations** – a dict {(from_index, to_index): alt},
    where alt is a dictionary including weight and/or edge
    properties to be changed following the split.


    * **allow_reverse** – If allow_reverse is True, then break_edge will
    attempt to break both (from_index, to_index) and, failing that,
    will attempt to break (to_index, from_index).



* **Returns**

    list of MoleculeGraphs.



#### substitute_group(index, func_grp, strategy, bond_order=1, graph_dict=None, strategy_params=None)
Builds off of Molecule.substitute to replace an atom in self.molecule
with a functional group. This method also amends self.graph to
incorporate the new functional group.

NOTE: using a MoleculeGraph will generally produce a different graph
compared with using a Molecule or str (when not using graph_dict).


* **Parameters**


    * **index** – Index of atom to substitute.


    * **func_grp** – Substituent molecule. There are three options:


        1. Providing an actual molecule as the input. The first atom

        must be a DummySpecies X, indicating the position of
        nearest neighbor. The second atom must be the next
        nearest atom. For example, for a methyl group
        substitution, func_grp should be X-CH3, where X is the
        first site and C is the second site. What the code will
        do is to remove the index site, and connect the nearest
        neighbor to the C atom in CH3. The X-C bond indicates the
        directionality to connect the atoms.


        2. A string name. The molecule will be obtained from the

        relevant template in func_groups.json.


        3. A MoleculeGraph object.



    * **strategy** – Class from pymatgen.analysis.local_env.


    * **bond_order** – A specified bond order to calculate the bond
    length between the attached functional group and the nearest
    neighbor site. Defaults to 1.


    * **graph_dict** – Dictionary representing the bonds of the functional
    group (format: {(u, v): props}, where props is a dictionary of
    properties, including weight. If None, then the algorithm
    will attempt to automatically determine bonds using one of
    a list of strategies defined in pymatgen.analysis.local_env.


    * **strategy_params** – dictionary of keyword arguments for strategy.
    If None, default parameters will be used.



#### _static_ with_edges(molecule: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), edges: dict[tuple[int, int], dict])
Constructor for MoleculeGraph, using pre-existing or pre-defined edges
with optional edge parameters.


* **Parameters**


    * **molecule** – Molecule object


    * **edges** – dict representing the bonds of the functional
    group (format: {(u, v): props}, where props is a dictionary of
    properties, including weight. Props should be None if no
    additional properties are to be specified.



* **Returns**

    mg, a MoleculeGraph



#### _classmethod_ with_empty_graph(molecule, name='bonds', edge_weight_name=None, edge_weight_units=None)
Constructor for MoleculeGraph, returns a MoleculeGraph
object with an empty graph (no edges, only nodes defined
that correspond to Sites in Molecule).


* **Parameters**


    * **(****Molecule****)** (*molecule*) –


    * **(****str****)** (*edge_weight_units*) – name of graph, e.g. “bonds”


    * **(****str****)** – name of edge weights,
    e.g. “bond_length” or “exchange_constant”


    * **(****str****)** – name of edge weight units
    e.g. “Å” or “eV”



* **Return (MoleculeGraph)**



#### _static_ with_local_env_strategy(molecule, strategy)
Constructor for MoleculeGraph, using a strategy
from pymatgen.analysis.local_env.


* **Parameters**


    * **molecule** – Molecule object


    * **strategy** – an instance of a
    pymatgen.analysis.local_env.NearNeighbors object



* **Returns**

    mg, a MoleculeGraph



### _class_ StructureGraph(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), graph_data=None)
Bases: `MSONable`

This is a class for annotating a Structure with bond information, stored in the form
of a graph. A “bond” does not necessarily have to be a chemical bond, but can store
any kind of information that connects two Sites.

If constructing this class manually, use the with_empty_graph method or
with_local_env_strategy method (using an algorithm provided by the local_env
module, such as O’Keeffe).
This class that contains connection information: relationships between sites
represented by a Graph structure, and an associated structure object.

StructureGraph uses the NetworkX package to store and operate on the graph itself, but
contains a lot of helper methods to make associating a graph with a given
crystallographic structure easier.
Use cases for this include storing bonding information, NMR J-couplings,
Heisenberg exchange parameters, etc.
For periodic graphs, class stores information on the graph edges of what lattice
image the edge belongs to.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object to be analyzed.


    * **graph_data** (*dict*) – Dictionary containing graph information. Not intended to be
    constructed manually see as_dict method for format.



#### _classmethod_ _edges_to_string(g)

#### add_edge(from_index, to_index, from_jimage=(0, 0, 0), to_jimage=None, weight=None, warn_duplicates=True, edge_properties=None)
Add edge to graph.

Since physically a ‘bond’ (or other connection
between sites) doesn’t have a direction, from_index,
from_jimage can be swapped with to_index, to_jimage.

However, images will always be shifted so that
from_index < to_index and from_jimage becomes (0, 0, 0).


* **Parameters**


    * **from_index** – index of site connecting from


    * **to_index** – index of site connecting to


    * **ints****)** (*to_jimage** (**tuple of*) – lattice vector of periodic
    image, e.g. (1, 0, 0) for periodic image in +x direction


    * **ints****)** – lattice vector of image


    * **(****float****)** (*weight*) – e.g. bond length


    * **(****bool****)** (*warn_duplicates*) – if True, will warn if
    trying to add duplicate edges (duplicate edges will not
    be added in either case)


    * **(****dict****)** (*edge_properties*) – any other information to
    store on graph edges, similar to Structure’s site_properties



#### alter_edge(from_index, to_index, to_jimage=None, new_weight=None, new_edge_properties=None)
Alters either the weight or the edge_properties of
an edge in the StructureGraph.


* **Parameters**


    * **from_index** – int


    * **to_index** – int


    * **to_jimage** – tuple


    * **new_weight** – alter_edge does not require
    that weight be altered. As such, by default, this
    is None. If weight is to be changed, it should be a
    float.


    * **new_edge_properties** – alter_edge does not require
    that edge_properties be altered. As such, by default,
    this is None. If any edge properties are to be changed,
    it should be a dictionary of edge properties to be changed.



#### as_dict()
As in pymatgen.core.Structure except
with using to_dict_of_dicts from NetworkX
to store graph information.


#### break_edge(from_index, to_index, to_jimage=None, allow_reverse=False)
Remove an edge from the StructureGraph. If no image is given, this method will fail.


* **Parameters**


    * **from_index** – int


    * **to_index** – int


    * **to_jimage** – tuple


    * **allow_reverse** – If allow_reverse is True, then break_edge will
    attempt to break both (from_index, to_index) and, failing that,
    will attempt to break (to_index, from_index).



#### diff(other, strict=True)
Compares two StructureGraphs. Returns dict with
keys ‘self’, ‘other’, ‘both’ with edges that are
present in only one StructureGraph (‘self’ and
‘other’), and edges that are present in both.

The Jaccard distance is a simple measure of the
dissimilarity between two StructureGraphs (ignoring
edge weights), and is defined by 1 - (size of the
intersection / size of the union) of the sets of
edges. This is returned with key ‘dist’.

Important note: all node indices are in terms
of the StructureGraph this method is called
from, not the ‘other’ StructureGraph: there
is no guarantee the node indices will be the
same if the underlying Structures are ordered
differently.


* **Parameters**


    * **other** – StructureGraph


    * **strict** – if False, will compare bonds
    from different Structures, with node indices
    replaced by Species strings, will not count
    number of occurrences of bonds



#### draw_graph_to_file(filename='graph', diff=None, hide_unconnected_nodes=False, hide_image_edges=True, edge_colors=False, node_labels=False, weight_labels=False, image_labels=False, color_scheme='VESTA', keep_dot=False, algo='fdp')
Draws graph using GraphViz.

The networkx graph object itself can also be drawn
with networkx’s in-built graph drawing methods, but
note that this might give misleading results for
multigraphs (edges are super-imposed on each other).

If visualization is difficult to interpret,
hide_image_edges can help, especially in larger
graphs.


* **Parameters**


    * **filename** – filename to output, will detect filetype
    from extension (any graphviz filetype supported, such as
    pdf or png)


    * **(****StructureGraph****)** (*diff*) – an additional graph to
    compare with, will color edges red that do not exist in diff
    and edges green that are in diff graph but not in the
    reference graph


    * **hide_unconnected_nodes** – if True, hide unconnected
    nodes


    * **hide_image_edges** – if True, do not draw edges that
    go through periodic boundaries


    * **(****bool****)** (*keep_dot*) – if True, use node colors to
    color edges


    * **(****bool****)** – if True, label nodes with
    species and site index


    * **(****bool****)** – if True, label edges with
    weights


    * **(****bool****)** – if True, label edges with
    their periodic images (usually only used for debugging,
    edges to periodic images always appear as dashed lines)


    * **(****str****)** (*color_scheme*) – “VESTA” or “JMOL”


    * **(****bool****)** – keep GraphViz .dot file for later
    visualization


    * **algo** – any graphviz algo, “neato” (for simple graphs)
    or “fdp” (for more crowded graphs) usually give good outputs



#### _property_ edge_weight_name()
Name of the edge weight property of graph


#### _property_ edge_weight_unit()
Units of the edge weight property of graph


#### _classmethod_ from_dict(d)
As in pymatgen.core.Structure except
restoring graphs using from_dict_of_dicts
from NetworkX to restore graph information.


#### get_connected_sites(n, jimage=(0, 0, 0))
Returns a named tuple of neighbors of site n:
periodic_site, jimage, index, weight.
Index is the index of the corresponding site
in the original structure, weight can be
None if not defined.
:param n: index of Site in Structure
:param jimage: lattice vector of site


* **Returns**

    list of ConnectedSite tuples,
    sorted by closest first.



#### get_coordination_of_site(n)
Returns the number of neighbors of site n.
In graph terms, simply returns degree
of node corresponding to site n.
:param n: index of site
:return (int):


#### get_subgraphs_as_molecules(use_weights=False)
Retrieve subgraphs as molecules, useful for extracting
molecules from periodic crystals.

Will only return unique molecules, not any duplicates
present in the crystal (a duplicate defined as an
isomorphic subgraph).


* **Parameters**

    **(****bool****)** (*use_weights*) – If True, only treat subgraphs
    as isomorphic if edges have the same weights. Typically,
    this means molecules will need to have the same bond
    lengths to be defined as duplicates, otherwise bond
    lengths can differ. This is a fairly robust approach,
    but will treat e.g. enantiomers as being duplicates.



* **Returns**

    list of unique Molecules in Structure



#### insert_node(idx, species, coords, coords_are_cartesian=False, validate_proximity=False, site_properties=None, edges=None)
A wrapper around Molecule.insert(), which also incorporates the new
site into the MoleculeGraph.


* **Parameters**


    * **idx** – Index at which to insert the new site


    * **species** – Species for the new site


    * **coords** – 3x1 array representing coordinates of the new site


    * **coords_are_cartesian** – Whether coordinates are cartesian.
    Defaults to False.


    * **validate_proximity** – For Molecule.insert(); if True (default
    False), distance will be checked to ensure that site can be safely
    added.


    * **site_properties** – Site properties for Molecule


    * **edges** – List of dicts representing edges to be added to the
    MoleculeGraph. These edges must include the index of the new site i,
    and all indices used for these edges should reflect the
    MoleculeGraph AFTER the insertion, NOT before. Each dict should at
    least have a “to_index” and “from_index” key, and can also have a
    “weight” and a “properties” key.



#### _property_ name()
Name of graph


#### remove_nodes(indices)
A wrapper for Molecule.remove_sites().


* **Parameters**

    **indices** – list of indices in the current Molecule (and graph) to
    be removed.



#### set_node_attributes()
Gives each node a “specie” and a “coords” attribute, updated with the
current species and coordinates.


#### sort(key=None, reverse=False)
Same as Structure.sort(). Also remaps nodes in graph.


* **Parameters**


    * **key** – key to sort by


    * **reverse** – reverse sort order



#### substitute_group(index, func_grp, strategy, bond_order=1, graph_dict=None, strategy_params=None)
Builds off of Structure.substitute to replace an atom in self.structure
with a functional group. This method also amends self.graph to
incorporate the new functional group.

NOTE: Care must be taken to ensure that the functional group that is
substituted will not place atoms to close to each other, or violate the
dimensions of the Lattice.


* **Parameters**


    * **index** – Index of atom to substitute.


    * **func_grp** – Substituent molecule. There are two options:


        1. Providing an actual Molecule as the input. The first atom

        must be a DummySpecies X, indicating the position of
        nearest neighbor. The second atom must be the next
        nearest atom. For example, for a methyl group
        substitution, func_grp should be X-CH3, where X is the
        first site and C is the second site. What the code will
        do is to remove the index site, and connect the nearest
        neighbor to the C atom in CH3. The X-C bond indicates the
        directionality to connect the atoms.


        2. A string name. The molecule will be obtained from the

        relevant template in func_groups.json.



    * **strategy** – Class from pymatgen.analysis.local_env.


    * **bond_order** – A specified bond order to calculate the bond
    length between the attached functional group and the nearest
    neighbor site. Defaults to 1.


    * **graph_dict** – Dictionary representing the bonds of the functional
    group (format: {(u, v): props}, where props is a dictionary of
    properties, including weight. If None, then the algorithm
    will attempt to automatically determine bonds using one of
    a list of strategies defined in pymatgen.analysis.local_env.


    * **strategy_params** – dictionary of keyword arguments for strategy.
    If None, default parameters will be used.



#### _property_ types_and_weights_of_connections()
Extract a dictionary summarizing the types and weights
of edges in the graph.


* **Returns**

    A dictionary with keys specifying the
    species involved in a connection in alphabetical order
    (e.g. string ‘Fe-O’) and values which are a list of
    weights for those connections (e.g. bond lengths).



#### types_of_coordination_environments(anonymous=False)
Extract information on the different co-ordination environments
present in the graph.


* **Parameters**

    **anonymous** – if anonymous, will replace specie names
    with A, B, C, etc.



* **Returns**

    a list of co-ordination environments,
    e.g. [‘Mo-S(6)’, ‘S-Mo(3)’]



#### _property_ weight_statistics()
Extract a statistical summary of edge weights present in
the graph.


* **Returns**

    A dict with an ‘all_weights’ list, ‘minimum’,
    ‘maximum’, ‘median’, ‘mean’, ‘std_dev’



#### _static_ with_edges(structure, edges)
Constructor for MoleculeGraph, using pre-existing or pre-defined edges
with optional edge parameters.


* **Parameters**


    * **molecule** – Molecule object


    * **edges** – dict representing the bonds of the functional
    group (format: {(from_index, to_index, from_image, to_image): props},
    where props is a dictionary of properties, including weight.
    Props should be None if no additional properties are to be
    specified.



* **Returns**

    sg, a StructureGraph



#### _classmethod_ with_empty_graph(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), name: str = 'bonds', edge_weight_name: str | None = None, edge_weight_units: str | None = None)
Constructor for an empty StructureGraph, i.e. no edges, containing only nodes corresponding
to sites in Structure.


* **Parameters**


    * **structure** – A pymatgen Structure object.


    * **name** – Name of the graph, e.g. “bonds”.


    * **edge_weight_name** – Name of the edge weights, e.g. “bond_length” or “exchange_constant”.


    * **edge_weight_units** – Name of the edge weight units, e.g. “Å” or “eV”.



* **Returns**

    an empty graph with no edges, only nodes defined

        that correspond to sites in Structure.




* **Return type**

    StructureGraph



#### _static_ with_local_env_strategy(structure, strategy, weights=False, edge_properties=False)
Constructor for StructureGraph, using a strategy
from pymatgen.analysis.local_env.


* **Parameters**


    * **structure** – Structure object


    * **strategy** – an instance of a
    pymatgen.analysis.local_env.NearNeighbors object


    * **weights** – if True, use weights from local_env class
    (consult relevant class for their meaning)


    * **edge_properties** – if True, edge_properties from neighbors will be used



### _compare(g1, g2, i1, i2)
Helper function called by isomorphic to ensure comparison of node identities.


### _igraph_from_nxgraph(graph)
Helper function that converts a networkx graph object into an igraph graph object.


### _isomorphic(frag1: Graph, frag2: Graph)
Internal function to check if two graph objects are isomorphic, using igraph if
if is available and networkx if it is not.

## pymatgen.analysis.hhi module

This module is used to estimate the Herfindahl-Hirschman Index, or HHI, of
chemical compounds. The HHI is a measure of how geographically confined or
dispersed the elements comprising a compound are. A low HHI is desirable
because it means the component elements are geographically dispersed.

Data/strategy from “Data-Driven Review of Thermoelectric Materials:
Performance and Resource Considerations” by Gaultois et al., published
in Chemistry of Materials (2013).

## pymatgen.analysis.interface module

This module provides classes to store, generate, and manipulate material interfaces.

## pymatgen.analysis.interface_reactions module

This module provides a class to predict and analyze interfacial reactions between two
solids, with or without an open element (e.g., flowing O2).


### _class_ GrandPotentialInterfacialReactivity(c1: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition), c2: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition), grand_pd: GrandPotentialPhaseDiagram, pd_non_grand: PhaseDiagram, include_no_mixing_energy: bool = False, norm: bool = True, use_hull_energy: bool = True)
Bases: `InterfacialReactivity`

Extends upon InterfacialReactivity to allow for modelling possible reactions
at the interface between two solids in the presence of an open element. The
thermodynamics of the open system are provided by the user via the
GrandPotentialPhaseDiagram class.


* **Parameters**


    * **c1** – Reactant 1 composition


    * **c2** – Reactant 2 composition


    * **grand_pd** – Grand potential phase diagram object built from all elements in
    composition c1 and c2.


    * **include_no_mixing_energy** – No_mixing_energy for a reactant is the
    opposite number of its energy above grand potential convex hull. In
    cases where reactions involve elements reservoir, this param
    determines whether no_mixing_energy of reactants will be included
    in the final reaction energy calculation. By definition, if pd is
    not a GrandPotentialPhaseDiagram object, this param is False.


    * **pd_non_grand** – PhaseDiagram object but not
    GrandPotentialPhaseDiagram object built from elements in c1 and c2.


    * **norm** – Whether or not the total number of atoms in composition
    of reactant will be normalized to 1.


    * **use_hull_energy** – Whether or not use the convex hull energy for
    a given composition for reaction energy calculation. If false,
    the energy of ground state structure will be used instead.
    Note that in case when ground state can not be found for a
    composition, convex hull energy will be used associated with a
    warning message.



#### _get_grand_potential(composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))
Computes the grand potential Phi at a given composition and
chemical potential(s).


* **Parameters**

    **composition** – Composition object.



* **Returns**

    Grand potential at a given composition at chemical potential(s).



#### _get_reactants(x: float)
Returns a list of relevant reactant compositions given an x coordinate.


#### get_no_mixing_energy()
Generates the opposite number of energy above grand potential
convex hull for both reactants.


* **Returns**

    [(reactant1, no_mixing_energy1),(reactant2,no_mixing_energy2)].



### _class_ InterfacialReactivity(c1: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition), c2: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition), pd: PhaseDiagram, norm: bool = True, use_hull_energy: bool = False, \*\*kwargs)
Bases: `MSONable`

Class for modeling an interface between two solids and its possible reactions.
The two reactants are provided as Composition objects (c1 and c2), along with the
relevant compositional PhaseDiagram object. Possible reactions are calculated by
finding all points along a tie-line between c1 and c2 where there is a “kink” in
the phase diagram; i.e. a point or facet of the phase diagram.

Please consider citing one or both of the following papers if you use this code
in your own work.

### References

Richards, W. D., Miara, L. J., Wang, Y., Kim, J. C., &amp; Ceder, G. (2015).
Interface stability in solid-state batteries. Chemistry of Materials, 28(1),
266-273. [https://doi.org/10.1021/acs.chemmater.5b04082](https://doi.org/10.1021/acs.chemmater.5b04082)

Xiao, Y., Wang, Y., Bo, S.-H., Kim, J. C., Miara, L. J., &amp; Ceder, G. (2019).
Understanding interface stability in solid-state batteries.
Nature Reviews Materials, 5(2), 105-126.
[https://doi.org/10.1038/s41578-019-0157-5](https://doi.org/10.1038/s41578-019-0157-5)


* **Parameters**


    * **c1** – Reactant 1 composition


    * **c2** – Reactant 2 composition


    * **pd** – Phase diagram object built from all elements in composition c1 and c2.


    * **norm** – Whether or not the total number of atoms in composition
    of reactant will be normalized to 1.


    * **use_hull_energy** – Whether or not use the convex hull energy for
    a given composition for reaction energy calculation. If false,
    the energy of ground state structure will be used instead.
    Note that in case when ground state can not be found for a
    composition, convex hull energy will be used associated with a
    warning message.



#### EV_TO_KJ_PER_MOL(_ = 96.485_ )

#### _static_ _convert(x: float, factor1: float, factor2: float)
Converts mixing ratio x in comp1 - comp2 tie line to that in
c1 - c2 tie line.


* **Parameters**


    * **x** – Mixing ratio x in comp1 - comp2 tie line, a float
    between 0 and 1.


    * **factor1** – Compositional ratio between composition c1 and
    processed composition comp1. E.g., factor for
    Composition(‘SiO2’) and Composition(‘O’) is 2.0.


    * **factor2** – Compositional ratio between composition c2 and
    processed composition comp2.



* **Returns**

    Mixing ratio in c1 - c2 tie line, a float between 0 and 1.



#### _get_elmt_amt_in_rxn(rxn: Reaction)
Computes total number of atoms in a reaction formula for elements
not in external reservoir. This method is used in the calculation
of reaction energy per mol of reaction formula.


* **Parameters**

    **rxn** – a Reaction object.



* **Returns**

    Total number of atoms for non_reservoir elements.



#### _get_energy(x)
Computes reaction energy in eV/atom at mixing ratio x : (1-x) for
self.comp1 : self.comp2.


* **Parameters**

    **x** (*float*) – Mixing ratio x of reactants, a float between 0 and 1.



* **Returns**

    Reaction energy.



#### _static_ _get_entry_energy(pd: PhaseDiagram, composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))
Finds the lowest entry energy for entries matching the composition.
Entries with non-negative formation energies are excluded. If no
entry is found, use the convex hull energy for the composition.


* **Parameters**


    * **pd** – Phase diagram object


    * **composition** – Composition object that the target entry should match



* **Returns**

    The lowest entry energy among entries matching the composition.



#### _get_matplotlib_figure()
Returns a matplotlib figure of reaction kinks diagram.


#### _get_original_composition_ratio(reaction)
Returns the molar mixing ratio between the reactants with ORIGINAL (
instead of processed) compositions for a reaction.


* **Parameters**

    **reaction** (*Reaction*) – Reaction object that contains the original
    reactant compositions.



* **Returns**

    The molar mixing ratio between the original reactant
    compositions for a reaction.



#### _static_ _get_plotly_annotations(x: list[float], y: list[float], reactions: list[Reaction])
Returns dictionary of annotations for the Plotly figure layout.


#### _get_plotly_figure()
Returns a Plotly figure of reaction kinks diagram.


#### _get_reactants(x: float)
Returns a list of relevant reactant compositions given an x coordinate.


#### _get_reaction(x: float)
Generates balanced reaction at mixing ratio x : (1-x) for
self.comp1 : self.comp2.


* **Parameters**

    **x** (*float*) – Mixing ratio x of reactants, a float between 0 and 1.



* **Returns**

    Reaction object.



#### _get_xaxis_title(latex: bool = True)
Returns the formatted title of the x axis (using either html/latex).


#### _static_ _reverse_convert(x: float, factor1: float, factor2: float)
Converts mixing ratio x in c1 - c2 tie line to that in
comp1 - comp2 tie line.


* **Parameters**


    * **x** – Mixing ratio x in c1 - c2 tie line, a float between
    0 and 1.


    * **factor1** – Compositional ratio between composition c1 and
    processed composition comp1. E.g., factor for
    Composition(‘SiO2’) and Composition(‘O’) is 2.


    * **factor2** – Compositional ratio between composition c2 and
    processed composition comp2.



* **Returns**

    Mixing ratio in comp1 - comp2 tie line, a float between 0 and 1.



#### _classmethod_ get_chempot_correction(element: str, temp: float, pres: float)
Get the normalized correction term Δμ for chemical potential of a gas
phase consisting of element at given temperature and pressure,
referenced to that in the standard state (T_std = 298.15 K,
T_std = 1 bar). The gas phase is limited to be one of O2, N2, Cl2,
F2, H2. Calculation formula can be found in the documentation of
Materials Project website.


* **Parameters**


    * **element** – The string representing the element.


    * **temp** – The temperature of the gas phase in Kelvin.


    * **pres** – The pressure of the gas phase in Pa.



* **Returns**

    The correction of chemical potential in eV/atom of the gas
    phase at given temperature and pressure.



#### get_critical_original_kink_ratio()
Returns a list of molar mixing ratio for each kink between ORIGINAL
(instead of processed) reactant compositions. This is the
same list as mixing ratio obtained from get_kinks method
if self.norm = False.


* **Returns**

    A list of floats representing molar mixing ratios between
    the original reactant compositions for each kink.



#### get_dataframe()
Returns a pandas DataFrame representation of the data produced by the
get_kinks() method.


#### get_kinks()
Finds all the kinks in mixing ratio where reaction products changes
along the tie-line of composition self.c1 and composition self.c2.


* **Returns**

    (index, mixing ratio, reaction energy in eV/atom, Reaction object, reaction
    energy per mol of formula in kJ/mol).



* **Return type**

    List object of tuples, each of which contains 5 elements



#### _property_ labels()
Returns a dictionary containing kink information:
{index: ‘x= mixing_ratio energy= reaction_energy reaction_equation’}.
E.g., {1: ‘x= 0 energy = 0 Mn -> Mn’,

> 2: ‘x= 0.5 energy = -15 O2 + Mn -> MnO2’,
> 3: ‘x= 1 energy = 0 O2 -> O2’}.


#### _property_ minimum()
Finds the minimum reaction energy E_min and corresponding
mixing ratio x_min.


* **Returns**

    Tuple (x_min, E_min).



#### plot(backend: Literal['plotly', 'matplotlib'] = 'plotly')
Plots reaction energy as a function of mixing ratio x in self.c1 - self.c2
tie line.


* **Parameters**

    **backend** (*"plotly"** | **"matplotlib"*) – Plotting library used to create the plot. Defaults to
    “plotly” but can also be “matplotlib”.



* **Returns**

    Plot of reaction energies as a function of mixing ratio



#### _property_ products()
List of formulas of potential products. E.g., [‘Li’,’O2’,’Mn’].

## pymatgen.analysis.local_env module

This module provides classes to perform analyses of
the local environments (e.g., finding near neighbors)
of single sites in molecules and structures.


### _class_ BrunnerNN_real(tol: float = 0.0001, cutoff=8.0)
Bases: `NearNeighbors`

Determine coordination number using Brunner’s algorithm which counts the
atoms that are within the largest gap in differences in real space
interatomic distances. This algorithm uses Brunner’s method of
largest gap in interatomic distances.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for bond determination
    (default: 1E-4).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for near-neighbor
    atoms. Defaults to 8.0.



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ BrunnerNN_reciprocal(tol: float = 0.0001, cutoff=8.0)
Bases: `NearNeighbors`

Determine coordination number using Brunner’s algorithm which counts the
atoms that are within the largest gap in differences in real space
interatomic distances. This algorithm uses Brunner’s method of
largest reciprocal gap in interatomic distances.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for bond determination
    (default: 1E-4).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for near-neighbor
    atoms. Defaults to 8.0.



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ BrunnerNN_relative(tol: float = 0.0001, cutoff=8.0)
Bases: `NearNeighbors`

Determine coordination number using Brunner’s algorithm which counts the
atoms that are within the largest gap in differences in real space
interatomic distances. This algorithm uses Brunner’s method of
of largest relative gap in interatomic distances.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for bond determination
    (default: 1E-4).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for near-neighbor
    atoms. Defaults to 8.0.



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ CovalentBondNN(tol: float = 0.2, order=True)
Bases: `NearNeighbors`

Determine near-neighbor sites and bond orders using built-in
pymatgen.Molecule CovalentBond functionality.

NOTE: This strategy is only appropriate for molecules, and not for
structures.


* **Parameters**


    * **tol** (*float*) – Tolerance for covalent bond checking.


    * **order** (*bool*) – If True (default), this class will compute bond
    orders. If False, bond lengths will be computed.



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_bonded_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), decorate: bool = False)
Obtain a MoleculeGraph object using this NearNeighbor class.


* **Parameters**


    * **structure** – Molecule object.


    * **decorate** (*bool*) – whether to annotate site properties


    * **by** (*with order parameters using neighbors determined*) –


    * **class** (*this NearNeighbor*) –



* **Returns**

    object from pymatgen.analysis.graphs



* **Return type**

    MoleculeGraph



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites and weights (orders) of bonds for a given
atom.


* **Parameters**


    * **structure** – input Molecule.


    * **n** – index of site for which to determine near neighbors.



* **Returns**

    [dict] representing a neighboring site and the type of
    bond present between site n and the neighboring site.



#### get_nn_shell_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), site_idx, shell)
Get a certain nearest neighbor shell for a certain site.

Determines all non-backtracking paths through the neighbor network
computed by get_nn_info. The weight is determined by multiplying
the weight of the neighbor at each hop through the network. For
example, a 2nd-nearest-neighbor that has a weight of 1 from its
1st-nearest-neighbor and weight 0.5 from the original site will
be assigned a weight of 0.5.

As this calculation may involve computing the nearest neighbors of
atoms multiple times, the calculation starts by computing all of the
neighbor info and then calling _get_nn_shell_info. If you are likely
to call this method for more than one site, consider calling get_all_nn
first and then calling this protected method yourself.


* **Parameters**


    * **structure** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Input structure


    * **site_idx** (*int*) – index of site for which to determine neighbor
    information.


    * **shell** (*int*) – Which neighbor shell to retrieve (1 == 1st NN shell)



* **Returns**

    list of dictionaries. Each entry in the list is information about

        a certain neighbor in the structure, in the same format as
        get_nn_info.




#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ Critic2NN()
Bases: `NearNeighbors`

Performs a topological analysis using critic2 to obtain
neighbor information, using a sum of atomic charge
densities. If an actual charge density is available
(e.g. from a VASP CHGCAR), see Critic2Caller directly
instead.

Init for Critic2NN.


#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_bonded_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), decorate: bool = False)

* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **decorate** (*bool**, **optional*) – Whether to decorate the structure. Defaults to False.



* **Returns**

    Bonded structure



* **Return type**

    StructureGraph



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ CrystalNN(weighted_cn=False, cation_anion=False, distance_cutoffs=(0.5, 1), x_diff_weight=3.0, porous_adjustment=True, search_cutoff=7, fingerprint_length=None)
Bases: `NearNeighbors`

This is a custom near-neighbor method intended for use in all kinds of periodic structures
(metals, minerals, porous structures, etc). It is based on a Voronoi algorithm and uses the
solid angle weights to determine the probability of various coordination environments. The
algorithm can also modify probability using smooth distance cutoffs as well as Pauling
electronegativity differences. The output can either be the most probable coordination
environment or a weighted list of coordination environments.

Initialize CrystalNN with desired parameters. Default parameters assume
“chemical bond” type behavior is desired. For geometric neighbor
finding (e.g., structural framework), set (i) distance_cutoffs=None,
(ii) x_diff_weight=0 and (optionally) (iii) porous_adjustment=False
which will disregard the atomic identities and perform best for a purely
geometric match.


* **Parameters**


    * **weighted_cn** – (bool) if set to True, will return fractional weights
    for each potential near neighbor.


    * **cation_anion** – (bool) if set True, will restrict bonding targets to
    sites with opposite or zero charge. Requires an oxidation states
    on all sites in the structure.


    * **distance_cutoffs** – ([float, float]) - if not None, penalizes neighbor
    distances greater than sum of covalent radii plus
    distance_cutoffs[0]. Distances greater than covalent radii sum
    plus distance_cutoffs[1] are enforced to have zero weight.


    * **x_diff_weight** – (float) - if multiple types of neighbor elements are
    possible, this sets preferences for targets with higher
    electronegativity difference.


    * **porous_adjustment** – (bool) - if True, readjusts Voronoi weights to
    better describe layered / porous structures


    * **search_cutoff** – (float) cutoff in Angstroms for initial neighbor
    search; this will be adjusted if needed internally


    * **fingerprint_length** – (int) if a fixed_length CN “fingerprint” is
    desired from get_nn_data(), set this parameter



#### _class_ NNData(all_nninfo, cn_weights, cn_nninfo)
Bases: `tuple`

Create new instance of NNData(all_nninfo, cn_weights, cn_nninfo)


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('all_nninfo', 'cn_weights', 'cn_nninfo'_ )

#### _classmethod_ _make(iterable)
Make a new NNData object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new NNData object replacing specified fields with new values


#### all_nninfo()
Alias for field number 0


#### cn_nninfo()
Alias for field number 2


#### cn_weights()
Alias for field number 1


#### _static_ _semicircle_integral(dist_bins, idx)
An internal method to get an integral between two bounds of a unit
semicircle. Used in algorithm to determine bond probabilities.


* **Parameters**


    * **dist_bins** – (float) list of all possible bond weights


    * **idx** – (float) index of starting bond weight



* **Returns**

    (float) integral of portion of unit semicircle



#### get_cn(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int, \*\*kwargs)
Get coordination number, CN, of site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine CN.


    * **use_weights** (*bool*) – flag indicating whether (True)
    to use weights for computing the coordination number
    or not (False, default: each coordinated site has equal
    weight).


    * **on_disorder** (*'take_majority_strict'** | **'take_majority_drop'** | **'take_max_species'** | **'error'*) – What to do when encountering a disordered structure. ‘error’ will raise ValueError.
    ‘take_majority_strict’ will use the majority specie on each site and raise
    ValueError if no majority exists. ‘take_max_species’ will use the first max specie
    on each site. For {{Fe: 0.4, O: 0.4, C: 0.2}}, ‘error’ and ‘take_majority_strict’
    will raise ValueError, while ‘take_majority_drop’ ignores this site altogether and
    ‘take_max_species’ will use Fe as the site specie.



* **Returns**

    coordination number.



* **Return type**

    cn (int or float)



#### get_cn_dict(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int, use_weights: bool = False, \*\*kwargs)
Get coordination number, CN, of each element bonded to site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure


    * **n** (*int*) – index of site for which to determine CN.


    * **use_weights** (*bool*) – flag indicating whether (True)
    to use weights for computing the coordination number
    or not (False, default: each coordinated site has equal
    weight).



* **Returns**

    dictionary of CN of each element bonded to site



* **Return type**

    cn (dict)



#### get_nn_data(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int, length=None)
The main logic of the method to compute near neighbor.


* **Parameters**


    * **structure** – (Structure) enclosing structure object


    * **n** – (int) index of target site to get NN info for


    * **length** – (int) if set, will return a fixed range of CN numbers



* **Returns**


    * all near neighbor sites with weights


    * a dict of CN -> weight


    * a dict of CN -> associated near neighbor sites




* **Return type**

    a namedtuple (NNData) object that contains



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor information.


* **Parameters**


    * **structure** – (Structure) pymatgen Structure


    * **n** – (int) index of target site



* **Returns**

    each dictionary provides information

        about a single near neighbor, where key ‘site’ gives access to the
        corresponding Site object, ‘image’ gives the image location, and
        ‘weight’ provides the weight that a given near-neighbor site contributes
        to the coordination number (1 or smaller), ‘site_index’ gives index of
        the corresponding site in the original structure.




* **Return type**

    siw (list[dict])



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



#### _static_ transform_to_length(nn_data, length)
Given NNData, transforms data to the specified fingerprint length


* **Parameters**


    * **nn_data** – (NNData)


    * **length** – (int) desired length of NNData.



### _class_ CutOffDictNN(cut_off_dict=None)
Bases: `NearNeighbors`

A basic NN class using a dictionary of fixed cut-off distances.
Only pairs of elements listed in the cut-off dictionary are considered
during construction of the neighbor lists.

Omit passing a dictionary for a Null/Empty NN class.


* **Parameters**


    * **cut_off_dict** (*dict**[**str**, **float**]*) – a dictionary


    * **distances** (*of cut-off*) – 2.0} for


    * **{** (*e.g.*) – 2.0} for


    * **Angstroms.** (*a maximum Fe-O bond length** of **2*) –


    * **listed** (*Bonds will only be created between pairs*) –


    * **dictionary.** (*in the cut-off*) –


    * **decorated** (*If your structure is oxidation state*) –


:param :
:param the cut-off distances will have to explicitly include:
:param the oxidation state: 2.0}.
:type the oxidation state: ‘Fe2+’, ‘O2-’
:param e.g. {: 2.0}.
:type e.g. {: ‘Fe2+’, ‘O2-’


#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### _classmethod_ from_preset(preset)
Initialize a CutOffDictNN according to a preset set of cutoffs.


* **Parameters**

    **preset** (*str*) – A preset name. The list of supported presets are:
    - “vesta_2019”: The distance cutoffs used by the VESTA

    > visualisation program.




* **Returns**

    A CutOffDictNN using the preset cut-off dictionary.



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ EconNN(tol: float = 0.2, cutoff: float = 10.0, cation_anion: bool = False, use_fictive_radius: bool = False)
Bases: `NearNeighbors`

Determines the average effective coordination number for each cation in a
given structure using Hoppe’s algorithm.

This method follows the procedure outlined in:

Hoppe, Rudolf. “Effective coordination numbers (ECoN) and mean fictive ionic
radii (MEFIR).” Zeitschrift für Kristallographie-Crystalline Materials
150.1-4 (1979): 23-52.


* **Parameters**


    * **tol** – Tolerance parameter for bond determination.


    * **cutoff** – Cutoff radius in Angstrom to look for near-neighbor atoms.


    * **cation_anion** – If set to True, will restrict bonding targets to
    sites with opposite or zero charge. Requires an oxidation states
    on all sites in the structure.


    * **use_fictive_radius** – Whether to use the fictive radius in the
    EcoN calculation. If False, the bond distance will be used.



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ IsayevNN(tol: float = 0.25, targets: [Element](pymatgen.core.md#pymatgen.core.periodic_table.Element) | list[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)] | None = None, cutoff: float = 13.0, allow_pathological: bool = False, extra_nn_info: bool = True, compute_adj_neighbors: bool = True)
Bases: `VoronoiNN`

Uses the algorithm defined in 10.1038/ncomms15679.

Sites are considered neighbors if (i) they share a Voronoi facet and (ii) the
bond distance is less than the sum of the Cordero covalent radii + 0.25 Å.


* **Parameters**


    * **tol** – Tolerance in Å for bond distances that are considered coordinated.


    * **targets** – Target element(s).


    * **cutoff** – Cutoff radius in Angstrom to look for near-neighbor atoms.


    * **allow_pathological** – Whether to allow infinite vertices in Voronoi
    coordination.


    * **extra_nn_info** – Add all polyhedron info to get_nn_info.


    * **compute_adj_neighbors** – Whether to compute which neighbors are adjacent. Turn
    off for faster performance.



#### _filter_nns(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int, nns: dict[str, Any])
Extract and filter the NN info into the format needed by NearestNeighbors.


* **Parameters**


    * **structure** – The structure.


    * **n** – The central site index.


    * **nns** – Nearest neighbor information for the structure.



* **Returns**

    See get_nn_info for the format of the returned data.



#### get_all_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.



* **Returns**

    List of near neighbor information for each site. See get_nn_info for the
    format of the data for each site.



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor site information.

Gets the associated image locations and weights of the site with index n
in structure using Voronoi decomposition and distance cutoff.


* **Parameters**


    * **structure** – Input structure.


    * **n** – Index of site for which to determine near-neighbor sites.



* **Returns**

    List of dicts containing the near-neighbor information. Each dict has the
    keys:


    * ”site”: The near-neighbor site.


    * ”image”: The periodic image of the near-neighbor site.


    * ”weight”: The face weight of the Voronoi decomposition.


    * ”site_index”: The index of the near-neighbor site in the original
    structure.




### _class_ JmolNN(tol: float = 0.45, min_bond_distance: float = 0.4, el_radius_updates: dict[SpeciesLike, float] | None = None)
Bases: `NearNeighbors`

Determine near-neighbor sites and coordination number using an emulation
of Jmol’s default autoBond() algorithm. This version of the algorithm
does not take into account any information regarding known charge
states.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for bond determination
    Defaults to 0.56.


    * **min_bond_distance** (*float*) – minimum bond distance for consideration
    Defaults to 0.4.


    * **el_radius_updates** – (dict) symbol->float to override default atomic
    radii table values.



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_max_bond_distance(el1_sym, el2_sym)
Use Jmol algorithm to determine bond length from atomic parameters


* **Parameters**


    * **el1_sym** – (str) symbol of atom 1


    * **el2_sym** – (str) symbol of atom 2.



* **Returns**

    max bond length



* **Return type**

    float



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n using the bond identification
algorithm underlying Jmol.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near
    neighbors.



* **Returns**

    tuples, each one

        of which represents a neighbor site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ LocalStructOrderParams(types, parameters=None, cutoff=-10.0)
Bases: `object`

This class permits the calculation of various types of local
structure order parameters.


* **Parameters**


    * **types** (*[**string**]*) – list of strings representing the types of
    order parameters to be calculated. Note that multiple
    mentions of the same type may occur. Currently available
    types recognize following environments:

    > ”cn”: simple coordination number—normalized

    >     if desired;

    > ”sgl_bd”: single bonds;
    > “bent”: bent (angular) coordinations

    > > (Zimmermann & Jain, in progress, 2017);

    > ”T”: T-shape coordinations;
    > “see_saw_rect”: see saw-like coordinations;
    > “tet”: tetrahedra

    > > (Zimmermann et al., submitted, 2017);

    > ”oct”: octahedra

    >     (Zimmermann et al., submitted, 2017);

    > ”bcc”: body-centered cubic environments (Peters,


    >         1. Chem. Phys., 131, 244103, 2009);

    > ”tri_plan”: trigonal planar environments;
    > “sq_plan”: square planar environments;
    > “pent_plan”: pentagonal planar environments;
    > “tri_pyr”: trigonal pyramids (coordinated atom is in

    > > the center of the basal plane);

    > ”sq_pyr”: square pyramids;
    > “pent_pyr”: pentagonal pyramids;
    > “hex_pyr”: hexagonal pyramids;
    > “tri_bipyr”: trigonal bipyramids;
    > “sq_bipyr”: square bipyramids;
    > “pent_bipyr”: pentagonal bipyramids;
    > “hex_bipyr”: hexagonal bipyramids;
    > “cuboct”: cuboctahedra;
    > “q2”: motif-unspecific bond orientational order

    > > parameter (BOOP) of weight l=2 (Steinhardt
    > > et al., Phys. Rev. B, 28, 784-805, 1983);

    > ”q4”: BOOP of weight l=4;
    > “q6”: BOOP of weight l=6.
    > “reg_tri”: regular triangle with varying height

    > > to basal plane;

    > ”sq”: square coordination (cf., “reg_tri”);
    > “oct_legacy”: original Peters-style OP recognizing

    > > octahedral coordination environments
    > > (Zimmermann et al., J. Am. Chem. Soc.,
    > > 137, 13352-13361, 2015) that can, however,
    > > produce small negative values sometimes.

    > ”sq_pyr_legacy”: square pyramids (legacy);



    * **parameters** (*[**dict**]*) – list of dictionaries
    that store float-type parameters associated with the
    definitions of the different order parameters
    (length of list = number of OPs). If an entry
    is None, default values are used that are read from
    the op_params.yaml file. With few exceptions, 9 different
    parameters are used across all OPs:

    > ”norm”: normalizing constant (used in “cn”

    >     (default value: 1)).

    > ”TA”: target angle (TA) in fraction of 180 degrees

    >     (“bent” (1), “tet” (0.6081734479693927),
    >     “tri_plan” (0.66666666667), “pent_plan” (0.6),
    >     “sq_pyr_legacy” (0.5)).

    > ”IGW_TA”: inverse Gaussian width (IGW) for penalizing

    >     angles away from the target angle in inverse
    >     fractions of 180 degrees to (“bent” and “tet” (15),
    >     “tri_plan” (13.5), “pent_plan” (18),
    >     “sq_pyr_legacy” (30)).

    > ”IGW_EP”: IGW for penalizing angles away from the

    >     equatorial plane (EP) at 90 degrees (“T”, “see_saw_rect”,
    >     “oct”, “sq_plan”, “tri_pyr”, “sq_pyr”, “pent_pyr”,
    >     “hex_pyr”, “tri_bipyr”, “sq_bipyr”, “pent_bipyr”,
    >     “hex_bipyr”, and “oct_legacy” (18)).

    > ”fac_AA”: factor applied to azimuth angle (AA) in cosine

    >     term (“T”, “tri_plan”, and “sq_plan” (1), “tet”,
    >     “tri_pyr”, and “tri_bipyr” (1.5), “oct”, “sq_pyr”,
    >     “sq_bipyr”, and “oct_legacy” (2), “pent_pyr”
    >     and “pent_bipyr” (2.5), “hex_pyr” and
    >     “hex_bipyr” (3)).

    > ”exp_cos_AA”: exponent applied to cosine term of AA

    >     (“T”, “tet”, “oct”, “tri_plan”, “sq_plan”,
    >     “tri_pyr”, “sq_pyr”, “pent_pyr”, “hex_pyr”,
    >     “tri_bipyr”, “sq_bipyr”, “pent_bipyr”, “hex_bipyr”,
    >     and “oct_legacy” (2)).

    > ”min_SPP”: smallest angle (in radians) to consider

    >     a neighbor to be
    >     at South pole position (“see_saw_rect”, “oct”, “bcc”,
    >     “sq_plan”, “tri_bipyr”, “sq_bipyr”, “pent_bipyr”,
    >     “hex_bipyr”, “cuboct”, and “oct_legacy”
    >     (2.792526803190927)).

    > ”IGW_SPP”: IGW for penalizing angles away from South

    >     pole position (“see_saw_rect”, “oct”, “bcc”, “sq_plan”,
    >     “tri_bipyr”, “sq_bipyr”, “pent_bipyr”, “hex_bipyr”,
    >     “cuboct”, and “oct_legacy” (15)).

    > ”w_SPP”: weight for South pole position relative to

    >     equatorial positions (“see_saw_rect” and “sq_plan” (1),
    >     “cuboct” (1.8), “tri_bipyr” (2), “oct”,
    >     “sq_bipyr”, and “oct_legacy” (3), “pent_bipyr” (4),
    >     “hex_bipyr” (5), “bcc” (6)).



    * **cutoff** (*float*) – Cutoff radius to determine which nearest
    neighbors are supposed to contribute to the order
    parameters. If the value is negative the neighboring
    sites found by distance and cutoff radius are further
    pruned using the get_nn method from the
    VoronoiNN class.



#### \__supported_types(_ = ('cn', 'sgl_bd', 'bent', 'tri_plan', 'tri_plan_max', 'reg_tri', 'sq_plan', 'sq_plan_max', 'pent_plan', 'pent_plan_max', 'sq', 'tet', 'tet_max', 'tri_pyr', 'sq_pyr', 'sq_pyr_legacy', 'tri_bipyr', 'sq_bipyr', 'oct', 'oct_legacy', 'pent_pyr', 'hex_pyr', 'pent_bipyr', 'hex_bipyr', 'T', 'cuboct', 'cuboct_max', 'see_saw_rect', 'bcc', 'q2', 'q4', 'q6', 'oct_max', 'hex_plan_max', 'sq_face_cap_trig_pris'_ )

#### compute_trigonometric_terms(thetas, phis)
Computes trigonometric terms that are required to
calculate bond orientational order parameters using
internal variables.


* **Parameters**


    * **thetas** (*[**float**]*) – polar angles of all neighbors in radians.


    * **phis** (*[**float**]*) – azimuth angles of all neighbors in radians.
    The list of
    azimuth angles of all neighbors in radians. The list of
    azimuth angles is expected to have the same size as the
    list of polar angles; otherwise, a ValueError is raised.
    Also, the two lists of angles have to be coherent in
    order. That is, it is expected that the order in the list
    of azimuth angles corresponds to a distinct sequence of
    neighbors. And, this sequence has to equal the sequence
    of neighbors in the list of polar angles.



#### get_order_parameters(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int, indices_neighs: list[int] | None = None, tol: float = 0.0, target_spec=None)
Compute all order parameters of site n.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site in input structure,
    for which OPs are to be
    calculated. Note that we do not use the sites iterator
    here, but directly access sites via struct[index].


    * **indices_neighs** (*list**[**int**]*) – list of indices of those neighbors
    in Structure object
    structure that are to be considered for OP computation.
    This optional argument overwrites the way neighbors are
    to be determined as defined in the constructor (i.e.,
    Voronoi coordination finder via negative cutoff radius
    vs constant cutoff radius if cutoff was positive).
    We do not use information about the underlying
    structure lattice if the neighbor indices are explicitly
    provided. This has two important consequences. First,
    the input Structure object can, in fact, be a
    simple list of Site objects. Second, no nearest images
    of neighbors are determined when providing an index list.
    Note furthermore that this neighbor
    determination type ignores the optional target_spec
    argument.


    * **tol** (*float*) – threshold of weight
    (= solid angle / maximal solid angle)
    to determine if a particular pair is
    considered neighbors; this is relevant only in the case
    when Voronoi polyhedra are used to determine coordination


    * **target_spec** ([*Species*](pymatgen.core.md#pymatgen.core.periodic_table.Species)) – target species to be considered
    when calculating the order
    parameters of site n; None includes all species of input
    structure.



* **Returns**

    representing order parameters. Should it not be
    possible to compute a given OP for a conceptual reason, the
    corresponding entry is None instead of a float. For Steinhardt
    et al.’s bond orientational OPs and the other geometric OPs
    (“tet”, “oct”, “bcc”, etc.),
    this can happen if there is a single
    neighbor around site n in the structure because that
    does not permit calculation of angles between multiple
    neighbors.



* **Return type**

    [floats]



#### get_parameters(index)
Returns list of floats that represents
the parameters associated
with calculation of the order
parameter that was defined at the index provided.
Attention: the parameters do not need to equal those originally
inputted because of processing out of efficiency reasons.


* **Parameters**

    **index** (*int*) – index of order parameter for which associated parameters
    are to be returned.



* **Returns**

    parameters of a given OP.



* **Return type**

    [float]



#### get_q2(thetas=None, phis=None)
Calculates the value of the bond orientational order parameter of
weight l=2. If the function is called with non-empty lists of
polar and azimuthal angles the corresponding trigonometric terms
are computed afresh. Otherwise, it is expected that the
compute_trigonometric_terms function has been just called.


* **Parameters**


    * **thetas** (*[**float**]*) – polar angles of all neighbors in radians.


    * **phis** (*[**float**]*) – azimuth angles of all neighbors in radians.



* **Returns**

    bond orientational order parameter of weight l=2

        corresponding to the input angles thetas and phis.




* **Return type**

    float



#### get_q4(thetas=None, phis=None)
Calculates the value of the bond orientational order parameter of
weight l=4. If the function is called with non-empty lists of
polar and azimuthal angles the corresponding trigonometric terms
are computed afresh. Otherwise, it is expected that the
compute_trigonometric_terms function has been just called.


* **Parameters**


    * **thetas** (*[**float**]*) – polar angles of all neighbors in radians.


    * **phis** (*[**float**]*) – azimuth angles of all neighbors in radians.



* **Returns**

    bond orientational order parameter of weight l=4

        corresponding to the input angles thetas and phis.




* **Return type**

    float



#### get_q6(thetas=None, phis=None)
Calculates the value of the bond orientational order parameter of
weight l=6. If the function is called with non-empty lists of
polar and azimuthal angles the corresponding trigonometric terms
are computed afresh. Otherwise, it is expected that the
compute_trigonometric_terms function has been just called.


* **Parameters**


    * **thetas** (*[**float**]*) – polar angles of all neighbors in radians.


    * **phis** (*[**float**]*) – azimuth angles of all neighbors in radians.



* **Returns**

    bond orientational order parameter of weight l=6

        corresponding to the input angles thetas and phis.




* **Return type**

    float



#### get_type(index)
Return type of order parameter at the index provided and
represented by a short string.


* **Parameters**

    **index** (*int*) – index of order parameter for which type is
    to be returned.



* **Returns**

    OP type.



* **Return type**

    str



#### _property_ last_nneigh()
Returns:
int: the number of neighbors encountered during the most

> recent order parameter calculation. A value of -1 indicates
> that no such calculation has yet been performed for this instance.


#### _property_ num_ops()
Returns:
int: the number of different order parameters that are targeted to be calculated.


### _class_ MinimumDistanceNN(tol: float = 0.1, cutoff=10, get_all_sites=False)
Bases: `NearNeighbors`

Determine near-neighbor sites and coordination number using the
nearest neighbor(s) at distance, d_min, plus all neighbors
within a distance (1 + tol) \* d_min, where tol is a
(relative) distance tolerance parameter.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for neighbor identification
    (default: 0.1).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for trial
    near-neighbor sites (default: 10).


    * **get_all_sites** (*bool*) – If this is set to True then the neighbor
    sites are only determined by the cutoff radius, tol is ignored.



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n using the closest neighbor
distance-based method.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near
    neighbors.



* **Returns**

    dicts with (Site, array, float) each one of which represents a

        neighbor site, its image location, and its weight.




* **Return type**

    siw (list[dict])



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ MinimumOKeeffeNN(tol: float = 0.1, cutoff=10)
Bases: `NearNeighbors`

Determine near-neighbor sites and coordination number using the
neighbor(s) at closest relative distance, d_min_OKeffee, plus some
relative tolerance, where bond valence parameters from O’Keeffe’s
bond valence method (J. Am. Chem. Soc. 1991, 3226-3229) are used
to calculate relative distances.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for neighbor identification
    (default: 0.1).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for trial
    near-neighbor sites (default: 10).



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n using the closest relative
neighbor distance-based method with O’Keeffe parameters.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near
    neighbors.



* **Returns**

    tuples, each one

        of which represents a neighbor site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ MinimumVIRENN(tol: float = 0.1, cutoff=10)
Bases: `NearNeighbors`

Determine near-neighbor sites and coordination number using the
neighbor(s) at closest relative distance, d_min_VIRE, plus some
relative tolerance, where atom radii from the
ValenceIonicRadiusEvaluator (VIRE) are used
to calculate relative distances.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for neighbor identification
    (default: 0.1).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for trial
    near-neighbor sites (default: 10).



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n using the closest relative
neighbor distance-based method with VIRE atomic/ionic radii.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near
    neighbors.



* **Returns**

    tuples, each one

        of which represents a neighbor site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ NearNeighbors()
Bases: `object`

Base class to determine near neighbors that typically include nearest
neighbors and others that are within some tolerable distance.


#### _static_ _get_image(structure, site)
Private convenience method for get_nn_info,
gives lattice image from provided PeriodicSite and Structure.

Image is defined as displacement from original site in structure to a given site.
i.e. if structure has a site at (-0.1, 1.0, 0.3), then (0.9, 0, 2.3) -> jimage = (1, -1, 2).
Note that this method takes O(number of sites) due to searching an original site.


* **Parameters**


    * **structure** – Structure Object


    * **site** – PeriodicSite Object



* **Returns**

    ((int)\*3) Lattice image



* **Return type**

    image



#### _get_nn_shell_info(structure, all_nn_info, site_idx, shell, _previous_steps=frozenset({}), _cur_image=(0, 0, 0))
Private method for computing the neighbor shell information.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) –


    * **all_nn_info** (*[**[**dict**]**]*) –


    * **site_idx** (*int*) – information.


    * **shell** (*int**) **- Which neighbor shell to retrieve** (**1 == 1st NN shell*) –


    * **_previous_steps** (*{**(**site_idx**, **image}*) – Set of
    sites that have already been traversed.


    * **_cur_image** (*tuple*) –



* **Returns**

    list of dictionaries. Each entry in the list is information about

        a certain neighbor in the structure, in the same format as
        get_nn_info. Does not update the site positions




#### _static_ _get_original_site(structure, site)
Private convenience method for get_nn_info,
gives original site index from ProvidedPeriodicSite.


#### _property_ extend_structure_molecules(_: boo_ )
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_all_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Get a listing of all neighbors for all sites in a structure.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure



* **Returns**

    List of NN site information for each site in the structure. Each

        entry has the same format as get_nn_info




#### get_bonded_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), decorate: bool = False, weights: bool = True, edge_properties: bool = False, on_disorder: on_disorder_options = 'take_majority_strict')
Obtain a StructureGraph object using this NearNeighbor
class. Requires the optional dependency networkx
(pip install networkx).


* **Parameters**


    * **structure** – Structure object.


    * **decorate** (*bool*) – whether to annotate site properties with order parameters using neighbors
    determined by this NearNeighbor class


    * **weights** (*bool*) – whether to include edge weights from NearNeighbor class in StructureGraph


    * **edge_properties** (*bool*) –


    * **on_disorder** (*'take_majority_strict'** | **'take_majority_drop'** | **'take_max_species'** | **'error'*) – What to do when encountering a disordered structure. ‘error’ will raise ValueError.
    ‘take_majority_strict’ will use the majority specie on each site and raise
    ValueError if no majority exists. ‘take_max_species’ will use the first max specie
    on each site. For {{Fe: 0.4, O: 0.4, C: 0.2}}, ‘error’ and ‘take_majority_strict’
    will raise ValueError, while ‘take_majority_drop’ ignores this site altogether and
    ‘take_max_species’ will use Fe as the site specie.



* **Returns**

    object from pymatgen.analysis.graphs



* **Return type**

    StructureGraph



#### get_cn(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int, use_weights: bool = False, on_disorder: Literal['take_majority_strict', 'take_majority_drop', 'take_max_species', 'error'] = 'take_majority_strict')
Get coordination number, CN, of site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine CN.


    * **use_weights** (*bool*) – flag indicating whether (True) to use weights for computing the coordination
    number or not (False, default: each coordinated site has equal weight).


    * **on_disorder** (*'take_majority_strict'** | **'take_majority_drop'** | **'take_max_species'** | **'error'*) – What to do when encountering a disordered structure. ‘error’ will raise ValueError.
    ‘take_majority_strict’ will use the majority specie on each site and raise
    ValueError if no majority exists. ‘take_max_species’ will use the first max specie
    on each site. For {{Fe: 0.4, O: 0.4, C: 0.2}}, ‘error’ and ‘take_majority_strict’
    will raise ValueError, while ‘take_majority_drop’ ignores this site altogether and
    ‘take_max_species’ will use Fe as the site specie.



* **Returns**

    coordination number.



* **Return type**

    cn (int or float)



#### get_cn_dict(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int, use_weights: bool = False)
Get coordination number, CN, of each element bonded to site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure


    * **n** (*int*) – index of site for which to determine CN.


    * **use_weights** (*bool*) – flag indicating whether (True)
    to use weights for computing the coordination number
    or not (False, default: each coordinated site has equal
    weight).



* **Returns**

    dictionary of CN of each element bonded to site



* **Return type**

    cn (dict)



#### get_local_order_parameters(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Calculate those local structure order parameters for
the given site whose ideal CN corresponds to the
underlying motif (e.g., CN=4, then calculate the
square planar, tetrahedral, see-saw-like,
rectangular see-saw-like order parameters).


* **Parameters**


    * **structure** – Structure object


    * **n** (*int*) – site index.


Returns (dict[str, float]):

    A dict of order parameters (values) and the
    underlying motif type (keys; for example, tetrahedral).


#### get_nn(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get near neighbors of site with index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site in structure for which to determine
    neighbors.



* **Returns**

    near neighbors.



* **Return type**

    sites (list of Site objects)



#### get_nn_images(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get image location of all near neighbors of site with index n in
structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine the image
    location of near neighbors.



* **Returns**

    image locations of

        near neighbors.




* **Return type**

    images (list of 3D integer array)



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    information.



* **Returns**

    each dictionary provides information

        about a single near neighbor, where key ‘site’ gives access to the
        corresponding Site object, ‘image’ gives the image location, and
        ‘weight’ provides the weight that a given near-neighbor site contributes
        to the coordination number (1 or smaller), ‘site_index’ gives index of
        the corresponding site in the original structure.




* **Return type**

    siw (list[dict])



#### get_nn_shell_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), site_idx, shell)
Get a certain nearest neighbor shell for a certain site.

Determines all non-backtracking paths through the neighbor network
computed by get_nn_info. The weight is determined by multiplying
the weight of the neighbor at each hop through the network. For
example, a 2nd-nearest-neighbor that has a weight of 1 from its
1st-nearest-neighbor and weight 0.5 from the original site will
be assigned a weight of 0.5.

As this calculation may involve computing the nearest neighbors of
atoms multiple times, the calculation starts by computing all of the
neighbor info and then calling _get_nn_shell_info. If you are likely
to call this method for more than one site, consider calling get_all_nn
first and then calling this protected method yourself.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **site_idx** (*int*) – index of site for which to determine neighbor
    information.


    * **shell** (*int*) – Which neighbor shell to retrieve (1 == 1st NN shell)



* **Returns**

    list of dictionaries. Each entry in the list is information about

        a certain neighbor in the structure, in the same format as
        get_nn_info.




#### get_weights_of_nn_sites(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get weight associated with each near neighbor of site with
index n in structure.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine the weights.



* **Returns**

    near-neighbor weights.



* **Return type**

    weights (list of floats)



#### _property_ molecules_allowed(_: boo_ )
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed(_: boo_ )
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ OpenBabelNN(order=True)
Bases: `NearNeighbors`

Determine near-neighbor sites and bond orders using OpenBabel API.

NOTE: This strategy is only appropriate for molecules, and not for
structures.


* **Parameters**


    * **order** (*bool*) – True if bond order should be returned as a weight, False


    * **weight.** (*if bond length should be used as a*) –



#### _property_ extend_structure_molecules()
Do Molecules need to be converted to Structures to use
this NearNeighbors class? Note: this property is not defined for classes
for which molecules_allowed is False.


* **Type**

    Boolean property



#### get_bonded_structure(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), decorate: bool = False)
Obtain a MoleculeGraph object using this NearNeighbor
class. Requires the optional dependency networkx
(pip install networkx).


* **Parameters**


    * **structure** – Molecule object.


    * **decorate** (*bool*) – whether to annotate site properties


    * **by** (*with order parameters using neighbors determined*) –


    * **class** (*this NearNeighbor*) –



* **Returns**

    object from pymatgen.analysis.graphs



* **Return type**

    MoleculeGraph



#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites and weights (orders) of bonds for a given
atom.


* **Parameters**


    * **structure** – Molecule object.


    * **n** – index of site for which to determine near neighbors.



* **Returns**

    representing a neighboring site and the type of
    bond present between site n and the neighboring site.



* **Return type**

    (dict)



#### get_nn_shell_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), site_idx, shell)
Get a certain nearest neighbor shell for a certain site.

Determines all non-backtracking paths through the neighbor network
computed by get_nn_info. The weight is determined by multiplying
the weight of the neighbor at each hop through the network. For
example, a 2nd-nearest-neighbor that has a weight of 1 from its
1st-nearest-neighbor and weight 0.5 from the original site will
be assigned a weight of 0.5.

As this calculation may involve computing the nearest neighbors of
atoms multiple times, the calculation starts by computing all of the
neighbor info and then calling _get_nn_shell_info. If you are likely
to call this method for more than one site, consider calling get_all_nn
first and then calling this protected method yourself.


* **Parameters**


    * **structure** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Input structure


    * **site_idx** (*int*) – index of site for which to determine neighbor
    information.


    * **shell** (*int*) – Which neighbor shell to retrieve (1 == 1st NN shell)



* **Returns**

    list of dictionaries. Each entry in the list is information about

        a certain neighbor in the structure, in the same format as
        get_nn_info.




#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _class_ ValenceIonicRadiusEvaluator(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Bases: `object`

Computes site valences and ionic radii for a structure using bond valence
analyzer.


* **Parameters**

    **structure** – pymatgen.core.structure.Structure.



#### _get_ionic_radii()
Computes ionic radii of elements for all sites in the structure.
If valence is zero, atomic radius is used.


#### _get_valences()
Computes ionic valences of elements for all sites in the structure.


#### _property_ radii()
List of ionic radii of elements in the order of sites.


#### _property_ structure()
Returns oxidation state decorated structure.


#### _property_ valences()
List of oxidation states of elements in the order of sites.


### _class_ VoronoiNN(tol=0, targets=None, cutoff=13.0, allow_pathological=False, weight='solid_angle', extra_nn_info=True, compute_adj_neighbors=True)
Bases: `NearNeighbors`

Uses a Voronoi algorithm to determine near neighbors for each site in a
structure.


* **Parameters**


    * **tol** (*float*) – tolerance parameter for near-neighbor finding. Faces that are
    smaller than tol fraction of the largest face are not included in the
    tessellation. (default: 0).


    * **targets** ([*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)* or **list** of **Elements*) – target element(s).


    * **cutoff** (*float*) – cutoff radius in Angstrom to look for near-neighbor
    atoms. Defaults to 13.0.


    * **allow_pathological** (*bool*) – whether to allow infinite vertices in
    determination of Voronoi coordination.


    * **weight** (*string*) – available in get_voronoi_polyhedra)


    * **extra_nn_info** (*bool*) –


    * **compute_adj_neighbors** (*bool*) – adjacent. Turn off for faster performance.



#### _extract_cell_info(site_idx, sites, targets, voro, compute_adj_neighbors=False)
Get the information about a certain atom from the results of a tessellation.


* **Parameters**


    * **site_idx** (*int*) –


    * **sites** (*[*[*Site*](pymatgen.core.md#pymatgen.core.sites.Site)*]*) –


    * **targets** (*[*[*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)*]*) –


    * **qvoronoi** (*voro - Output of*) –


    * **compute_adj_neighbors** (*boolean*) –



* **Returns**

    A dict of sites sharing a common Voronoi facet. Key is facet id

        (not useful) and values are dictionaries containing statistics
        about the facet:

        >
        > * site: Pymatgen site


        > * solid_angle - Solid angle subtended by face


        > * angle_normalized - Solid angle normalized such that the

        >     faces with the largest


        > * area - Area of the facet


        > * face_dist - Distance between site n and the facet


        > * volume - Volume of Voronoi cell for this face


        > * n_verts - Number of vertices on the facet


        > * adj_neighbors - Facet id’s for the adjacent neighbors




#### _extract_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), nns)
Given Voronoi NNs, extract the NN info in the form needed by NearestNeighbors.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure being evaluated


    * **nns** (*[**dicts**]*) – Nearest neighbor information for a structure



* **Returns**

    See nn_info



* **Return type**

    (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### get_all_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))

* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.



* **Returns**

    All nn info for all sites.



#### get_all_voronoi_polyhedra(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Get the Voronoi polyhedra for all site in a simulation cell.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure to be evaluated



* **Returns**

    A dict of sites sharing a common Voronoi facet with the site
    n mapped to a directory containing statistics about the facet:

    >
    > * solid_angle - Solid angle subtended by face


    > * angle_normalized - Solid angle normalized such that the

    >     faces with the largest


    > * area - Area of the facet


    > * face_dist - Distance between site n and the facet


    > * volume - Volume of Voronoi cell for this face


    > * n_verts - Number of vertices on the facet




#### get_nn_info(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Get all near-neighbor sites as well as the associated image locations
and weights of the site with index n in structure
using Voronoi decomposition.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site for which to determine near-neighbor
    sites.



* **Returns**

    tuples, each one

        of which represents a coordinated site, its image location,
        and its weight.




* **Return type**

    siw (list of tuples ([Site](pymatgen.core.md#pymatgen.core.sites.Site), array, float))



#### get_voronoi_polyhedra(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n: int)
Gives a weighted polyhedra around a site.

See ref: A Proposed Rigorous Definition of Coordination Number,
M. O’Keeffe, Acta Cryst. (1979). A35, 772-775


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure for which to evaluate the
    coordination environment.


    * **n** (*int*) – site index.



* **Returns**

    A dict of sites sharing a common Voronoi facet with the site
    n mapped to a directory containing statistics about the facet:

    >
    > * solid_angle - Solid angle subtended by face


    > * angle_normalized - Solid angle normalized such that the

    >     faces with the largest


    > * area - Area of the facet


    > * face_dist - Distance between site n and the facet


    > * volume - Volume of Voronoi cell for this face


    > * n_verts - Number of vertices on the facet




#### _property_ molecules_allowed()
can this NearNeighbors class be used with Molecule
objects?


* **Type**

    Boolean property



#### _property_ structures_allowed()
can this NearNeighbors class be used with Structure
objects?


* **Type**

    Boolean property



### _get_default_radius(site)
An internal method to get a “default” covalent/element radius.


* **Parameters**

    **site** – (Site)



* **Returns**

    Covalent radius of element on site, or Atomic radius if unavailable



### _get_elements(site)
Get the list of elements for a Site.


* **Parameters**

    **site** ([*Site*](pymatgen.core.md#pymatgen.core.sites.Site)) – Site to assess



* **Returns**

    List of elements



* **Return type**

    [[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)]



### _get_fictive_ionic_radius(site: [Site](pymatgen.core.md#pymatgen.core.sites.Site), neighbor: [PeriodicNeighbor](pymatgen.core.md#pymatgen.core.structure.PeriodicNeighbor))
Get fictive ionic radius.

Follows equation 1 of:

Hoppe, Rudolf. “Effective coordination numbers (ECoN) and mean fictive ionic
radii (MEFIR).” Zeitschrift für Kristallographie-Crystalline Materials
150.1-4 (1979): 23-52.


* **Parameters**


    * **site** – The central site.


    * **site.** (*neighbor neighboring*) –



* **Returns**

    Hoppe’s fictive ionic radius.



### _get_mean_fictive_ionic_radius(fictive_ionic_radii: list[float], minimum_fir: float | None = None)
Returns the mean fictive ionic radius.

Follows equation 2:

Hoppe, Rudolf. “Effective coordination numbers (ECoN) and mean fictive ionic
radii (MEFIR).” Zeitschrift für Kristallographie-Crystalline Materials
150.1-4 (1979): 23-52.


* **Parameters**


    * **fictive_ionic_radii** – List of fictive ionic radii for a center site
    and its neighbors.


    * **minimum_fir** – Minimum fictive ionic radius to use.



* **Returns**

    Hoppe’s mean fictive ionic radius.



### _get_radius(site)
An internal method to get the expected radius for a site with
oxidation state.


* **Parameters**

    **site** – (Site)



* **Returns**

    ionic, covalent, or atomic.
    Returns 0 if no oxidation state or appropriate radius is found.



* **Return type**

    Oxidation-state dependent radius



### _get_vire(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure) | [IStructure](pymatgen.core.md#pymatgen.core.structure.IStructure))
Get the ValenceIonicRadiusEvaluator object for an structure taking
advantage of caching.


* **Parameters**

    **structure** – A structure.



* **Returns**

    Output of ValenceIonicRadiusEvaluator(structure)



### _get_vire_istructure(structure: [IStructure](pymatgen.core.md#pymatgen.core.structure.IStructure))
Get the ValenceIonicRadiusEvaluator object for an immutable structure
taking advantage of caching.


* **Parameters**

    **structure** – A structure.



* **Returns**

    Output of ValenceIonicRadiusEvaluator(structure)



### _handle_disorder(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), on_disorder: Literal['take_majority_strict', 'take_majority_drop', 'take_max_species', 'error'])
What to do in bonding and coordination number analysis if a site is disordered.


### _is_in_targets(site, targets)
Test whether a site contains elements in the target list.


* **Parameters**


    * **site** ([*Site*](pymatgen.core.md#pymatgen.core.sites.Site)) – Site to assess


    * **targets** (*[*[*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)*]*) –



* **Returns**

    (boolean) Whether this site contains a certain list of elements



### get_neighbors_of_site_with_index(struct, n, approach='min_dist', delta=0.1, cutoff=10)
Returns the neighbors of a given site using a specific neighbor-finding
method.


* **Parameters**


    * **struct** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site in Structure object for which motif type
    is to be determined.


    * **approach** (*str*) – type of neighbor-finding approach, where
    “min_dist” will use the MinimumDistanceNN class,
    “voronoi” the VoronoiNN class, “min_OKeeffe” the
    MinimumOKeeffe class, and “min_VIRE” the MinimumVIRENN class.


    * **delta** (*float*) – tolerance involved in neighbor finding.


    * **cutoff** (*float*) – (large) radius to find tentative neighbors.



* **Returns**

    neighbor sites.



### get_okeeffe_distance_prediction(el1, el2)
Returns an estimate of the bond valence parameter (bond length) using
the derived parameters from ‘Atoms Sizes and Bond Lengths in Molecules
and Crystals’ (O’Keeffe & Brese, 1991). The estimate is based on two
experimental parameters: r and c. The value for r  is based off radius,
while c is (usually) the Allred-Rochow electronegativity. Values used
are *not* generated from pymatgen, and are found in
‘okeeffe_params.json’.


* **Parameters**


    * **el1** ([*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)) – two Element objects


    * **el2** ([*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)) – two Element objects



* **Returns**

    a float value of the predicted bond length



### get_okeeffe_params(el_symbol)
Returns the elemental parameters related to atom size and
electronegativity which are used for estimating bond-valence
parameters (bond length) of pairs of atoms on the basis of data
provided in ‘Atoms Sizes and Bond Lengths in Molecules and Crystals’
(O’Keeffe & Brese, 1991).


* **Parameters**

    **el_symbol** (*str*) – element symbol.



* **Returns**

    atom-size (‘r’) and electronegativity-related (‘c’) parameter.



* **Return type**

    (dict)



### gramschmidt(vin, uin)
Returns that part of the first input vector
that is orthogonal to the second input vector.
The output vector is not normalized.


* **Parameters**


    * **vin** (*numpy array*) – first input vector


    * **uin** (*numpy array*) – second input vector



### metal_edge_extender(mol_graph, cutoff: float = 2.5, metals: list | tuple | None = ('Li', 'Mg', 'Ca', 'Zn', 'B', 'Al'), coordinators: list | tuple = ('O', 'N', 'F', 'S', 'Cl'))
Function to identify and add missed coordinate bond edges for metals.


* **Parameters**


    * **mol_graph** – pymatgen.analysis.graphs.MoleculeGraph object


    * **cutoff** – cutoff in Angstrom. Metal-coordinator sites that are closer
    together than this value will be considered coordination bonds.
    If the MoleculeGraph contains a metal, but no coordination bonds are found
    with the chosen cutoff, the cutoff will be increased by 1 Angstrom
    and another attempt will be made to identify coordination bonds.


    * **metals** – Species considered metals for the purpose of identifying
    missed coordinate bond edges. The set {“Li”, “Mg”, “Ca”, “Zn”, “B”, “Al”}
    (default) corresponds to the settings used in the LIBE dataset.
    Alternatively, set to None to cause any Species classified as a metal
    by Specie.is_metal to be considered a metal.


    * **coordinators** – Possible coordinating species to consider when identifying
    missed coordinate bonds. The default set {“O”, “N”, “F”, “S”, “Cl”} was
    used in the LIBE dataset.



* **Returns**

    pymatgen.analysis.graphs.MoleculeGraph object with additional

        metal bonds (if any found) added




* **Return type**

    mol_graph



### oxygen_edge_extender(mol_graph: MoleculeGraph)
Identify and add missed O-C or O-H bonds. This is particularly
important when oxygen is forming three bonds, e.g. in H3O+ or XOH2+.
See [https://github.com/materialsproject/pymatgen/pull/2903](https://github.com/materialsproject/pymatgen/pull/2903) for details.


* **Parameters**

    **mol_graph** (*MoleculeGraph*) – molecule graph to extend



* **Returns**

    object with additional O-C or O-H bonds added (if any found)



* **Return type**

    MoleculeGraph



### site_is_of_motif_type(struct, n, approach='min_dist', delta=0.1, cutoff=10, thresh=None)
Returns the motif type of the site with index n in structure struct;
currently featuring “tetrahedral”, “octahedral”, “bcc”, and “cp”
(close-packed: fcc and hcp) as well as “square pyramidal” and
“trigonal bipyramidal”. If the site is not recognized,
“unrecognized” is returned. If a site should be assigned to two
different motifs, “multiple assignments” is returned.


* **Parameters**


    * **struct** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **n** (*int*) – index of site in Structure object for which motif type
    is to be determined.


    * **approach** (*str*) – type of neighbor-finding approach, where
    “min_dist” will use the MinimumDistanceNN class,
    “voronoi” the VoronoiNN class, “min_OKeeffe” the
    MinimumOKeeffe class, and “min_VIRE” the MinimumVIRENN class.


    * **delta** (*float*) – tolerance involved in neighbor finding.


    * **cutoff** (*float*) – (large) radius to find tentative neighbors.


    * **thresh** (*dict*) – thresholds for motif criteria (currently, required
    keys and their default values are “qtet”: 0.5,
    “qoct”: 0.5, “qbcc”: 0.5, “q6”: 0.4).



* **Returns**

    motif type



* **Return type**

    str



### solid_angle(center, coords)
Helper method to calculate the solid angle of a set of coords from the
center.


* **Parameters**


    * **center** (*3x1 array*) – Center to measure solid angle from.


    * **coords** (*Nx3 array*) – List of coords to determine solid angle.



* **Returns**

    The solid angle.



### vol_tetra(vt1, vt2, vt3, vt4)
Calculate the volume of a tetrahedron, given the four vertices of vt1,
vt2, vt3 and vt4.


* **Parameters**


    * **vt1** (*array-like*) – coordinates of vertex 1.


    * **vt2** (*array-like*) – coordinates of vertex 2.


    * **vt3** (*array-like*) – coordinates of vertex 3.


    * **vt4** (*array-like*) – coordinates of vertex 4.



* **Returns**

    volume of the tetrahedron.



* **Return type**

    (float)


## pymatgen.analysis.molecule_matcher module

This module provides classes to perform fitting of molecule with arbitrary
atom orders.
This module is supposed to perform exact comparisons without the atom order
correspondence prerequisite, while molecule_structure_comparator is supposed
to do rough comparisons with the atom order correspondence prerequisite.

The implementation is based on an excellent python package called rmsd that
you can find at [https://github.com/charnley/rmsd](https://github.com/charnley/rmsd).


### _class_ AbstractMolAtomMapper()
Bases: `MSONable`

Abstract molecular atom order mapping class. A mapping will be able to
find the uniform atom order of two molecules that can pair the
geometrically equivalent atoms.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _classmethod_ from_dict(d)

* **Parameters**

    **(****)** (*d*) – Dict.



* **Returns**

    AbstractMolAtomMapper



#### _abstract_ get_molecule_hash(mol)
Defines a hash for molecules. This allows molecules to be grouped
efficiently for comparison.


* **Parameters**

    **mol** – The molecule. OpenBabel OBMol or pymatgen Molecule object



* **Returns**

    A hashable object. Examples can be string formulas, etc.



#### _abstract_ uniform_labels(mol1, mol2)
Pair the geometrically equivalent atoms of the molecules.


* **Parameters**


    * **mol1** – First molecule. OpenBabel OBMol or pymatgen Molecule object.


    * **mol2** – Second molecule. OpenBabel OBMol or pymatgen Molecule object.



* **Returns**

    (list1, list2) if uniform atom order is found. list1 and list2
    are for mol1 and mol2, respectively. Their length equal
    to the number of atoms. They represents the uniform atom order
    of the two molecules. The value of each element is the original
    atom index in mol1 or mol2 of the current atom in uniform atom
    order.
    (None, None) if unform atom is not available.



### _class_ BruteForceOrderMatcher(target: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Bases: `KabschMatcher`

Finding the best match between molecules by selecting molecule order
with the smallest RMSD from all the possible order combinations.

### Notes

When aligning molecules, the atoms of the two molecules **must** have same number
of atoms from the same species.

Constructor of the matcher object.


* **Parameters**

    **target** – a Molecule object used as a target during the alignment



#### fit(p: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), ignore_warning=False)
Order, rotate and transform p molecule according to the best match.

A ValueError will be raised when the total number of possible combinations
become unfeasible (more than a million combinations).


* **Parameters**


    * **p** – a Molecule object what will be matched with the target one.


    * **ignore_warning** – ignoring error when the number of combination is too large



* **Returns**

    Rotated and translated of the p Molecule object
    rmsd: Root-mean-square-deviation between p_prime and the target



* **Return type**

    p_prime



#### match(mol: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), ignore_warning: bool = False)
Similar as KabschMatcher.match but this method also finds the order of
atoms which belongs to the best match.

A ValueError will be raised when the total number of possible combinations
become unfeasible (more than a million combination).


* **Parameters**


    * **mol** – a Molecule object what will be matched with the target one.


    * **ignore_warning** – ignoring error when the number of combination is too large



* **Returns**

    The indices of atoms
    U: 3x3 rotation matrix
    V: Translation vector
    rmsd: Root mean squared deviation between P and Q



* **Return type**

    inds



#### _static_ permutations(atoms)
Generates all the possible permutations of atom order. To achieve better
performance all the cases where the atoms are different has been ignored.


### _class_ GeneticOrderMatcher(target: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), threshold: float)
Bases: `KabschMatcher`

This method was inspired by genetic algorithms and tries to match molecules
based on their already matched fragments.

It uses the fact that when two molecule is matching their sub-structures have to match as well.
The main idea here is that in each iteration (generation) we can check the match of all possible
fragments and ignore those which are not feasible.

Although in the worst case this method has N! complexity (same as the brute force one),
in practice it performs much faster because many of the combination can be eliminated
during the fragment matching.

### Notes

This method very robust and returns with all the possible orders.

There is a well known weakness/corner case: The case when there is
a outlier with large deviation with a small index might be ignored.
This happens due to the nature of the average function
used to calculate the RMSD for the fragments.

When aligning molecules, the atoms of the two molecules **must** have the
same number of atoms from the same species.

Constructor of the matcher object.


* **Parameters**


    * **target** – a Molecule object used as a target during the alignment


    * **threshold** – value used to match fragments and prune configuration



#### fit(p: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Order, rotate and transform all of the matched p molecule
according to the given threshold.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    p_prime: Rotated and translated of the p Molecule object
    rmsd: Root-mean-square-deviation between p_prime and the target



* **Return type**

    Array of the possible matches where the elements are



#### match(p: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Similar as KabschMatcher.match but this method also finds all of the
possible atomic orders according to the threshold.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    inds: The indices of atoms
    U: 3x3 rotation matrix
    V: Translation vector
    rmsd: Root mean squared deviation between P and Q



* **Return type**

    Array of the possible matches where the elements are



#### permutations(p: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Generates all of possible permutations of atom order according the threshold.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    Array of index arrays



### _class_ HungarianOrderMatcher(target: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Bases: `KabschMatcher`

This method pre-aligns the molecules based on their principal inertia
axis and then re-orders the input atom list using the Hungarian method.

### Notes

This method cannot guarantee the best match but is very fast.

When aligning molecules, the atoms of the two molecules **must** have same number
of atoms from the same species.

Constructor of the matcher object.


* **Parameters**

    **target** – a Molecule object used as a target during the alignment



#### fit(p: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Order, rotate and transform p molecule according to the best match.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    Rotated and translated of the p Molecule object
    rmsd: Root-mean-square-deviation between p_prime and the target



* **Return type**

    p_prime



#### _static_ get_principal_axis(coords, weights)
Get the molecule’s principal axis.


* **Parameters**


    * **coords** – coordinates of atoms


    * **weights** – the weight use for calculating the inertia tensor



* **Returns**

    Array of dim 3 containing the principal axis



#### match(p: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Similar as KabschMatcher.match but this method also finds the order of
atoms which belongs to the best match.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    The indices of atoms
    U: 3x3 rotation matrix
    V: Translation vector
    rmsd: Root mean squared deviation between P and Q



* **Return type**

    inds



#### _static_ permutations(p_atoms, p_centroid, p_weights, q_atoms, q_centroid, q_weights)
Generates two possible permutations of atom order. This method uses the principle component
of the inertia tensor to prealign the molecules and hungarian method to determine the order.
There are always two possible permutation depending on the way to pre-aligning the molecules.


* **Parameters**


    * **p_atoms** – atom numbers


    * **p_centroid** – array of atom positions


    * **p_weights** – array of atom weights


    * **q_atoms** – atom numbers


    * **q_centroid** – array of atom positions


    * **q_weights** – array of atom weights



* **Yields**

    *perm_inds* – array of atoms’ order



#### _static_ rotation_matrix_vectors(v1, v2)
Returns the rotation matrix that rotates v1 onto v2 using
Rodrigues’ rotation formula.

See more: [https://math.stackexchange.com/a/476311](https://math.stackexchange.com/a/476311)


* **Parameters**


    * **v1** – initial vector


    * **v2** – target vector



* **Returns**

    3x3 rotation matrix



### _class_ InchiMolAtomMapper(angle_tolerance=10.0)
Bases: `AbstractMolAtomMapper`

Pair atoms by inchi labels.


* **Parameters**

    **angle_tolerance** (*float*) – Angle threshold to assume linear molecule. In degrees.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### _static_ _align_heavy_atoms(mol1, mol2, vmol1, vmol2, ilabel1, ilabel2, eq_atoms)
Align the label of topologically identical atoms of second molecule
towards first molecule.


* **Parameters**


    * **mol1** – First molecule. OpenBabel OBMol object


    * **mol2** – Second molecule. OpenBabel OBMol object


    * **vmol1** – First virtual molecule constructed by centroids. OpenBabel
    OBMol object


    * **vmol2** – First virtual molecule constructed by centroids. OpenBabel
    OBMol object


    * **ilabel1** – inchi label map of the first molecule


    * **ilabel2** – inchi label map of the second molecule


    * **eq_atoms** – equivalent atom labels



* **Returns**

    corrected inchi labels of heavy atoms of the second molecule



#### _static_ _align_hydrogen_atoms(mol1, mol2, heavy_indices1, heavy_indices2)
Align the label of topologically identical atoms of second molecule
towards first molecule.


* **Parameters**


    * **mol1** – First molecule. OpenBabel OBMol object


    * **mol2** – Second molecule. OpenBabel OBMol object


    * **heavy_indices1** – inchi label map of the first molecule


    * **heavy_indices2** – label map of the second molecule



* **Returns**

    corrected label map of all atoms of the second molecule



#### _static_ _get_elements(mol, label)
The elements of the atoms in the specified order.


* **Parameters**


    * **mol** – The molecule. OpenBabel OBMol object.


    * **label** – The atom indices. List of integers.



* **Returns**

    Elements. List of integers.



#### _static_ _group_centroid(mol, ilabels, group_atoms)
Calculate the centroids of a group atoms indexed by the labels of inchi.


* **Parameters**


    * **mol** – The molecule. OpenBabel OBMol object


    * **ilabel** – inchi label map



* **Returns**

    Centroid. Tuple (x, y, z)



#### _static_ _inchi_labels(mol)
Get the inchi canonical labels of the heavy atoms in the molecule.


* **Parameters**

    **mol** – The molecule. OpenBabel OBMol object



* **Returns**

    The label mappings. List of tuple of canonical label,
    original label
    List of equivalent atoms.



#### _is_molecule_linear(mol)
Is the molecule a linear one.


* **Parameters**

    **mol** – The molecule. OpenBabel OBMol object.



* **Returns**

    Boolean value.



#### _virtual_molecule(mol, ilabels, eq_atoms)
Create a virtual molecule by unique atoms, the centroids of the
equivalent atoms.


* **Parameters**


    * **mol** – The molecule. OpenBabel OBMol object


    * **ilabels** – inchi label map


    * **eq_atoms** – equivalent atom labels


    * **farthest_group_idx** – The equivalent atom group index in which
    there is the farthest atom to the centroid



* **Returns**

    The virtual molecule



#### as_dict()

* **Returns**

    MSONable dict.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict Representation.



* **Returns**

    InchiMolAtomMapper



#### get_molecule_hash(mol)
Return inchi as molecular hash.


#### uniform_labels(mol1, mol2)

* **Parameters**


    * **mol1** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Molecule 1


    * **mol2** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Molecule 2.



* **Returns**

    Labels



### _class_ IsomorphismMolAtomMapper()
Bases: `AbstractMolAtomMapper`

Pair atoms by isomorphism permutations in the OpenBabel::OBAlign class.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()

* **Returns**

    Jsonable dict.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    IsomorphismMolAtomMapper



#### get_molecule_hash(mol)
Return inchi as molecular hash.


#### uniform_labels(mol1, mol2)
Pair the geometrically equivalent atoms of the molecules.
Calculate RMSD on all possible isomorphism mappings and return mapping
with the least RMSD.


* **Parameters**


    * **mol1** – First molecule. OpenBabel OBMol or pymatgen Molecule object.


    * **mol2** – Second molecule. OpenBabel OBMol or pymatgen Molecule object.



* **Returns**

    (list1, list2) if uniform atom order is found. list1 and list2
    are for mol1 and mol2, respectively. Their length equal
    to the number of atoms. They represents the uniform atom order
    of the two molecules. The value of each element is the original
    atom index in mol1 or mol2 of the current atom in uniform atom
    order.
    (None, None) if unform atom is not available.



### _class_ KabschMatcher(target: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Bases: `MSONable`

Molecule matcher using Kabsch algorithm.

The Kabsch algorithm capable aligning two molecules by finding the parameters
(translation, rotation) which minimize the root-mean-square-deviation (RMSD) of
two molecules which are topologically (atom types, geometry) similar two each other.

### Notes

When aligning molecules, the atoms of the two molecules **must** be in the same
order for the results to be sensible.

Constructor of the matcher object.


* **Parameters**

    **target** – a Molecule object used as a target during the alignment



#### fit(p: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Rotate and transform p molecule according to the best match.


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    Rotated and translated of the p Molecule object
    rmsd: Root-mean-square-deviation between p_prime and the target



* **Return type**

    p_prime



#### _static_ kabsch(P: ndarray, Q: ndarray)
The Kabsch algorithm is a method for calculating the optimal rotation matrix
that minimizes the root mean squared deviation (RMSD) between two paired sets of points
P and Q, centered around the their centroid.

For more info see:
- [http://en.wikipedia.org/wiki/Kabsch_algorithm](http://en.wikipedia.org/wiki/Kabsch_algorithm) and
- [https://cnx.org/contents/HV-RsdwL@23/Molecular-Distance-Measures](https://cnx.org/contents/HV-RsdwL@23/Molecular-Distance-Measures)


* **Parameters**


    * **P** – Nx3 matrix, where N is the number of points.


    * **Q** – Nx3 matrix, where N is the number of points.



* **Returns**

    3x3 rotation matrix



* **Return type**

    U



#### match(p: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule))
Using the Kabsch algorithm the alignment of two molecules (P, Q)
happens in three steps:
- translate the P and Q into their centroid
- compute of the optimal rotation matrix (U) using Kabsch algorithm
- compute the translation (V) and rmsd.

The function returns the rotation matrix (U), translation vector (V),
and RMSD between Q and P’, where P’ is:

> P’ = P \* U + V


* **Parameters**

    **p** – a Molecule object what will be matched with the target one.



* **Returns**

    Rotation matrix (D,D)
    V: Translation vector (D)
    RMSD : Root mean squared deviation between P and Q



* **Return type**

    U



### _class_ MoleculeMatcher(tolerance: float = 0.01, mapper=None)
Bases: `MSONable`

Class to match molecules and identify whether molecules are the same.


* **Parameters**


    * **tolerance** (*float*) – RMSD difference threshold whether two molecules are
    different


    * **mapper** (*AbstractMolAtomMapper*) – MolAtomMapper object that is able to map the atoms of two
    molecule to uniform order.



#### _static_ _calc_rms(mol1, mol2, clabel1, clabel2)
Calculate the RMSD.


* **Parameters**


    * **mol1** – The first molecule. OpenBabel OBMol or pymatgen Molecule
    object


    * **mol2** – The second molecule. OpenBabel OBMol or pymatgen Molecule
    object


    * **clabel1** – The atom indices that can reorder the first molecule to
    uniform atom order


    * **clabel1** – The atom indices that can reorder the second molecule to
    uniform atom order



* **Returns**

    The RMSD.



#### as_dict()

* **Returns**

    MSONable dict.



#### fit(mol1, mol2)
Fit two molecules.


* **Parameters**


    * **mol1** – First molecule. OpenBabel OBMol or pymatgen Molecule object


    * **mol2** – Second molecule. OpenBabel OBMol or pymatgen Molecule object



* **Returns**

    A boolean value indicates whether two molecules are the same.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    MoleculeMatcher



#### get_rmsd(mol1, mol2)
Get RMSD between two molecule with arbitrary atom order.


* **Returns**

    RMSD if topology of the two molecules are the same
    Infinite if  the topology is different



#### group_molecules(mol_list)
Group molecules by structural equality.


* **Parameters**

    **mol_list** – List of OpenBabel OBMol or pymatgen objects



* **Returns**

    A list of lists of matched molecules
    Assumption: if s1=s2 and s2=s3, then s1=s3
    This may not be true for small tolerances.


## pymatgen.analysis.molecule_structure_comparator module

This module provides classes to comparison the structures of the two
molecule. As long as the two molecule have the same bond connection tables,
the molecules are deemed to be same. The atom in the two molecule must be
paired accordingly.
This module is supposed to perform rough comparisons with the atom order
correspondence prerequisite, while molecule_matcher is supposed to do exact
comparisons without the atom order correspondence prerequisite.


### _class_ CovalentRadius()
Bases: `object`

Covalent radius of the elements.

Beatriz C. et al. Dalton Trans. 2008, 2832-2838. [https://doi.org/10.1039/b801115j](https://doi.org/10.1039/b801115j)


#### radius(_ = {'Ac': 2.15, 'Ag': 1.45, 'Al': 1.21, 'Am': 1.8, 'Ar': 1.06, 'As': 1.19, 'At': 1.5, 'Au': 1.36, 'B': 0.84, 'Ba': 2.15, 'Be': 0.96, 'Bi': 1.48, 'Br': 1.2, 'C': 0.73, 'Ca': 1.76, 'Cd': 1.44, 'Ce': 2.04, 'Cl': 1.02, 'Cm': 1.69, 'Co': 1.38, 'Cr': 1.39, 'Cs': 2.44, 'Cu': 1.32, 'Dy': 1.92, 'Er': 1.89, 'Eu': 1.98, 'F': 0.57, 'Fe': 1.42, 'Fr': 2.6, 'Ga': 1.22, 'Gd': 1.96, 'Ge': 1.2, 'H': 0.31, 'He': 0.28, 'Hf': 1.75, 'Hg': 1.32, 'Ho': 1.92, 'I': 1.39, 'In': 1.42, 'Ir': 1.41, 'K': 2.03, 'Kr': 1.16, 'La': 2.07, 'Li': 1.28, 'Lu': 1.87, 'Mg': 1.41, 'Mn': 1.5, 'Mo': 1.54, 'N': 0.71, 'Na': 1.66, 'Nb': 1.64, 'Nd': 2.01, 'Ne': 0.58, 'Ni': 1.24, 'Np': 1.9, 'O': 0.66, 'Os': 1.44, 'P': 1.07, 'Pa': 2, 'Pb': 1.46, 'Pd': 1.39, 'Pm': 1.99, 'Po': 1.4, 'Pr': 2.03, 'Pt': 1.36, 'Pu': 1.87, 'Ra': 2.21, 'Rb': 2.2, 'Re': 1.51, 'Rh': 1.42, 'Rn': 1.5, 'Ru': 1.46, 'S': 1.05, 'Sb': 1.39, 'Sc': 1.7, 'Se': 1.2, 'Si': 1.11, 'Sm': 1.98, 'Sn': 1.39, 'Sr': 1.95, 'Ta': 1.7, 'Tb': 1.94, 'Tc': 1.47, 'Te': 1.38, 'Th': 2.06, 'Ti': 1.6, 'Tl': 1.45, 'Tm': 1.9, 'U': 1.96, 'V': 1.53, 'W': 1.62, 'Xe': 1.4, 'Y': 1.9, 'Yb': 1.87, 'Zn': 1.22, 'Zr': 1.75_ )

### _class_ MoleculeStructureComparator(bond_length_cap=0.3, covalent_radius={'Ac': 2.15, 'Ag': 1.45, 'Al': 1.21, 'Am': 1.8, 'Ar': 1.06, 'As': 1.19, 'At': 1.5, 'Au': 1.36, 'B': 0.84, 'Ba': 2.15, 'Be': 0.96, 'Bi': 1.48, 'Br': 1.2, 'C': 0.73, 'Ca': 1.76, 'Cd': 1.44, 'Ce': 2.04, 'Cl': 1.02, 'Cm': 1.69, 'Co': 1.38, 'Cr': 1.39, 'Cs': 2.44, 'Cu': 1.32, 'Dy': 1.92, 'Er': 1.89, 'Eu': 1.98, 'F': 0.57, 'Fe': 1.42, 'Fr': 2.6, 'Ga': 1.22, 'Gd': 1.96, 'Ge': 1.2, 'H': 0.31, 'He': 0.28, 'Hf': 1.75, 'Hg': 1.32, 'Ho': 1.92, 'I': 1.39, 'In': 1.42, 'Ir': 1.41, 'K': 2.03, 'Kr': 1.16, 'La': 2.07, 'Li': 1.28, 'Lu': 1.87, 'Mg': 1.41, 'Mn': 1.5, 'Mo': 1.54, 'N': 0.71, 'Na': 1.66, 'Nb': 1.64, 'Nd': 2.01, 'Ne': 0.58, 'Ni': 1.24, 'Np': 1.9, 'O': 0.66, 'Os': 1.44, 'P': 1.07, 'Pa': 2, 'Pb': 1.46, 'Pd': 1.39, 'Pm': 1.99, 'Po': 1.4, 'Pr': 2.03, 'Pt': 1.36, 'Pu': 1.87, 'Ra': 2.21, 'Rb': 2.2, 'Re': 1.51, 'Rh': 1.42, 'Rn': 1.5, 'Ru': 1.46, 'S': 1.05, 'Sb': 1.39, 'Sc': 1.7, 'Se': 1.2, 'Si': 1.11, 'Sm': 1.98, 'Sn': 1.39, 'Sr': 1.95, 'Ta': 1.7, 'Tb': 1.94, 'Tc': 1.47, 'Te': 1.38, 'Th': 2.06, 'Ti': 1.6, 'Tl': 1.45, 'Tm': 1.9, 'U': 1.96, 'V': 1.53, 'W': 1.62, 'Xe': 1.4, 'Y': 1.9, 'Yb': 1.87, 'Zn': 1.22, 'Zr': 1.75}, priority_bonds=(), priority_cap=0.8, ignore_ionic_bond=True, bond_13_cap=0.05)
Bases: `MSONable`

Class to check whether the connection tables of the two molecules are the
same. The atom in the two molecule must be paired accordingly.


* **Parameters**


    * **bond_length_cap** – The ratio of the elongation of the bond to be
    acknowledged. If the distance between two atoms is less than (
    empirical covalent bond length) X (1 + bond_length_cap), the bond
    between the two atoms will be acknowledged.


    * **covalent_radius** – The covalent radius of the atoms.
    dict (element symbol -> radius)


    * **priority_bonds** – The bonds that are known to be existed in the initial
    molecule. Such bonds will be acknowledged in a loose criteria.
    The index should start from 0.


    * **priority_cap** – The ratio of the elongation of the bond to be
    acknowledged for the priority bonds.



#### _get_bonds(mol)
Find all the bond in a molcule.


* **Parameters**

    **mol** – the molecule. pymatgen Molecule object



* **Returns**

    List of tuple. Each tuple correspond to a bond represented by the
    id of the two end atoms.



#### are_equal(mol1, mol2)
Compare the bond table of the two molecules.


* **Parameters**


    * **mol1** – first molecule. pymatgen Molecule object.


    * **mol2** – second molecules. pymatgen Molecule object.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    MoleculeStructureComparator



#### _static_ get_13_bonds(priority_bonds)

* **Parameters**

    **(****)** (*priority_bonds*) –


Returns:


#### halogen_list(_ = ('F', 'Cl', 'Br', 'I'_ )

#### ionic_element_list(_ = ('Na', 'Mg', 'Al', 'Sc', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr'_ )
## pymatgen.analysis.nmr module

A module for NMR analysis.


### _class_ ChemicalShielding(cs_matrix, vscale=None)
Bases: [`SquareTensor`](pymatgen.core.md#pymatgen.core.tensors.SquareTensor)

This class extends the SquareTensor to perform extra analysis unique to
NMR Chemical shielding tensors.

Three notations to describe chemical shielding tensor (RK Harris; Magn. Resonance
Chem. 2008, 46, 582-598; DOI: 10.1002/mrc.2225) are supported.

Authors: Shyam Dwaraknath, Xiaohui Qu

Create a Chemical Shielding tensor.
Note that the constructor uses __new__
rather than __init__ according to the standard method of
subclassing numpy ndarrays.


* **Parameters**


    * **cs_matrix** (*1x3** or **3x3 array-like*) – the 3x3 array-like
    representing the chemical shielding tensor
    or a 1x3 array of the primary sigma values corresponding
    to the principal axis system


    * **vscale** (*6x1 array-like*) – 6x1 array-like scaling the
    Voigt-notation vector with the tensor entries



#### _class_ HaeberlenNotation(sigma_iso, delta_sigma_iso, zeta, eta)
Bases: `tuple`

Create new instance of HaeberlenNotation(sigma_iso, delta_sigma_iso, zeta, eta)


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('sigma_iso', 'delta_sigma_iso', 'zeta', 'eta'_ )

#### _classmethod_ _make(iterable)
Make a new HaeberlenNotation object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new HaeberlenNotation object replacing specified fields with new values


#### delta_sigma_iso()
Alias for field number 1


#### eta()
Alias for field number 3


#### sigma_iso()
Alias for field number 0


#### zeta()
Alias for field number 2


#### _class_ MarylandNotation(sigma_iso, omega, kappa)
Bases: `tuple`

Create new instance of MarylandNotation(sigma_iso, omega, kappa)


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('sigma_iso', 'omega', 'kappa'_ )

#### _classmethod_ _make(iterable)
Make a new MarylandNotation object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new MarylandNotation object replacing specified fields with new values


#### kappa()
Alias for field number 2


#### omega()
Alias for field number 1


#### sigma_iso()
Alias for field number 0


#### _class_ MehringNotation(sigma_iso, sigma_11, sigma_22, sigma_33)
Bases: `tuple`

Create new instance of MehringNotation(sigma_iso, sigma_11, sigma_22, sigma_33)


#### _asdict()
Return a new dict which maps field names to their values.


#### _field_defaults(_ = {_ )

#### _fields(_ = ('sigma_iso', 'sigma_11', 'sigma_22', 'sigma_33'_ )

#### _classmethod_ _make(iterable)
Make a new MehringNotation object from a sequence or iterable


#### _replace(\*\*kwds)
Return a new MehringNotation object replacing specified fields with new values


#### sigma_11()
Alias for field number 1


#### sigma_22()
Alias for field number 2


#### sigma_33()
Alias for field number 3


#### sigma_iso()
Alias for field number 0


#### _classmethod_ from_maryland_notation(sigma_iso, omega, kappa)
Initialize from Maryland notation.


* **Parameters**


    * **(****)** (*kappa*) –


    * **(****)** –


    * **(****)** –



* **Returns**

    ChemicalShielding



#### _property_ haeberlen_values()
the Chemical shielding tensor in Haeberlen Notation.


* **Type**

    Returns



#### _property_ maryland_values()
the Chemical shielding tensor in Maryland Notation.


* **Type**

    Returns



#### _property_ mehring_values()
the Chemical shielding tensor in Mehring Notation.


* **Type**

    Returns



#### _property_ principal_axis_system()
Returns a chemical shielding tensor aligned to the principle axis system
so that only the 3 diagonal components are non-zero.


### _class_ ElectricFieldGradient(efg_matrix, vscale=None)
Bases: [`SquareTensor`](pymatgen.core.md#pymatgen.core.tensors.SquareTensor)

This class extends the SquareTensor to perform extra analysis unique to
NMR Electric Field Gradient tensors in units of V/Angstrom^2.

Authors: Shyam Dwaraknath, Xiaohui Qu

Create a Chemical Shielding tensor.
Note that the constructor uses __new__
rather than __init__ according to the standard method of
subclassing numpy ndarrays.


* **Parameters**


    * **efg_matrix** (*1x3** or **3x3 array-like*) – the 3x3 array-like
    representing the electric field tensor
    or a 1x3 array of the primary values corresponding
    to the principal axis system


    * **vscale** (*6x1 array-like*) – 6x1 array-like scaling the
    voigt-notation vector with the tensor entries



#### _property_ V_xx()
First diagonal element.


* **Type**

    Returns



#### _property_ V_yy()
Second diagonal element.


* **Type**

    Returns



#### _property_ V_zz()
Third diagonal element.


* **Type**

    Returns



#### _property_ asymmetry()
Asymmetry of the electric field tensor defined as:
(V_yy - V_xx)/V_zz.


#### coupling_constant(specie)
Computes the coupling constant C_q as defined in:

    Wasylishen R E, Ashbrook S E, Wimperis S. NMR of quadrupolar nuclei
    in solid materials[M]. John Wiley & Sons, 2012. (Chapter 3.2).

C_q for a specific atom type for this electric field tensor:

    > C_q=e\*Q\*V_zz/h

    h: Planck’s constant
    Q: nuclear electric quadrupole moment in mb (millibarn
    e: elementary proton charge


* **Parameters**

    **specie** – flexible input to specify the species at this site.
    Can take a isotope or element string, Species object,
    or Site object



* **Returns**

    the coupling constant as a FloatWithUnit in MHz



#### _property_ principal_axis_system()
Returns a electric field gradient tensor aligned to the principle axis system so that only the 3 diagonal
components are non-zero.

## pymatgen.analysis.phase_diagram module

This module defines tools to generate and analyze phase diagrams.


### _class_ CompoundPhaseDiagram(entries, terminal_compositions, normalize_terminal_compositions=True)
Bases: `PhaseDiagram`

Generates phase diagrams from compounds as terminations instead of
elements.

Initializes a CompoundPhaseDiagram.


* **Parameters**


    * **entries** (*[**PDEntry**]*) – Sequence of input entries. For example,
    if you want a Li2O-P2O5 phase diagram, you might have all
    Li-P-O entries as an input.


    * **terminal_compositions** (*list**[*[*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)*]*) – Terminal compositions of
    phase space. In the Li2O-P2O5 example, these will be the
    Li2O and P2O5 compositions.


    * **normalize_terminal_compositions** (*bool*) – Whether to normalize the
    terminal compositions to a per atom basis. If normalized,
    the energy above hulls will be consistent
    for comparison across systems. Non-normalized terminals are
    more intuitive in terms of compositional breakdowns.



#### amount_tol(_ = 1e-0_ )

#### as_dict()

* **Returns**

    MSONable dictionary representation of CompoundPhaseDiagram.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – dictionary representation of CompoundPhaseDiagram.



* **Returns**

    CompoundPhaseDiagram



#### transform_entries(entries, terminal_compositions)
Method to transform all entries to the composition coordinate in the
terminal compositions. If the entry does not fall within the space
defined by the terminal compositions, they are excluded. For example,
Li3PO4 is mapped into a Li2O:1.5, P2O5:0.5 composition. The terminal
compositions are represented by DummySpecies.


* **Parameters**


    * **entries** – Sequence of all input entries


    * **terminal_compositions** – Terminal compositions of phase space.



* **Returns**

    Sequence of TransformedPDEntries falling within the phase space.



### _class_ GrandPotPDEntry(entry, chempots, name=None)
Bases: `PDEntry`

A grand potential pd entry object encompassing all relevant data for phase
diagrams. Chemical potentials are given as a element-chemical potential
dict.


* **Parameters**


    * **entry** – A PDEntry-like object.


    * **chempots** – Chemical potential specification as {Element: float}.


    * **name** – Optional parameter to name the entry. Defaults to the reduced
    chemical formula of the original entry.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()

* **Returns**

    MSONable dictionary representation of GrandPotPDEntry.



#### _property_ chemical_energy()
The chemical energy term mu\*N in the grand potential.


* **Returns**

    The chemical energy term mu\*N in the grand potential



#### _property_ composition(_: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition_ )
The composition after removing free species.


* **Returns**

    Composition



#### _property_ energy()
Returns:
The grand potential energy.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – dictionary representation of GrandPotPDEntry.



* **Returns**

    GrandPotPDEntry



### _class_ GrandPotentialPhaseDiagram(entries, chempots, elements=None, \*, computed_data=None)
Bases: `PhaseDiagram`

A class representing a Grand potential phase diagram. Grand potential phase
diagrams are essentially phase diagrams that are open to one or more
components. To construct such phase diagrams, the relevant free energy is
the grand potential, which can be written as the Legendre transform of the
Gibbs free energy as follows.

Grand potential = G - u_X N_X

The algorithm is based on the work in the following papers:


1. S. P. Ong, L. Wang, B. Kang, and G. Ceder, Li-Fe-P-O2 Phase Diagram from
First Principles Calculations. Chem. Mater., 2008, 20(5), 1798-1807.
doi:10.1021/cm702327g


2. S. P. Ong, A. Jain, G. Hautier, B. Kang, G. Ceder, Thermal stabilities
of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first
principles calculations. Electrochem. Comm., 2010, 12(3), 427-430.
doi:10.1016/j.elecom.2010.01.010

Standard constructor for grand potential phase diagram.


* **Parameters**


    * **entries** (*[**PDEntry**]*) – A list of PDEntry-like objects having an
    energy, energy_per_atom and composition.


    * **(****{Element** (*chempots*) – float}): Specify the chemical potentials
    of the open elements.


    * **elements** (*[*[*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)*]*) – Optional list of elements in the phase
    diagram. If set to None, the elements are determined from
    the entries themselves.


    * **computed_data** (*dict*) – A dict containing pre-computed data. This allows
    PhaseDiagram object to be reconstituted without performing the
    expensive convex hull computation. The dict is the output from the
    PhaseDiagram._compute() method and is stored in PhaseDiagram.computed_data
    when generated for the first time.



#### as_dict()

* **Returns**

    MSONable dictionary representation of GrandPotentialPhaseDiagram.



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – dictionary representation of GrandPotentialPhaseDiagram.



* **Returns**

    GrandPotentialPhaseDiagram



### _class_ PDEntry(composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition), energy: float, name: str | None = None, attribute: object = None)
Bases: [`Entry`](pymatgen.entries.md#pymatgen.entries.Entry)

An object encompassing all relevant data for phase diagrams.


#### composition()
The composition associated with the PDEntry.


* **Type**

    [Composition](pymatgen.core.md#pymatgen.core.composition.Composition)



#### energy()
The energy associated with the entry.


* **Type**

    float



#### name()
A name for the entry. This is the string shown in the phase diagrams.
By default, this is the reduced formula for the composition, but can be
set to some other string for display purposes.


* **Type**

    str



#### attribute()
A arbitrary attribute. Can be used to specify that the
entry is a newly found compound, or to specify a particular label for
the entry, etc. An attribute can be anything but must be MSONable.


* **Type**

    MSONable



* **Parameters**


    * **composition** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Composition


    * **energy** (*float*) – Energy for composition.


    * **name** (*str*) – Optional parameter to name the entry. Defaults
    to the reduced chemical formula.


    * **attribute** – Optional attribute of the entry. Must be MSONable.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()

* **Returns**

    MSONable dictionary representation of PDEntry.



#### _property_ energy(_: floa_ )
Returns:
the energy of the entry.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – dictionary representation of PDEntry.



* **Returns**

    PDEntry



### _class_ PDPlotter(phasediagram: PhaseDiagram, show_unstable: float = 0.2, backend: Literal['plotly', 'matplotlib'] = 'plotly', ternary_style: Literal['2d', '3d'] = '2d', \*\*plotkwargs)
Bases: `object`

A plotting class for compositional phase diagrams.

To use, initialize this class with a PhaseDiagram object containing 1-4 components
and call get_plot() or show().


* **Parameters**


    * **phasediagram** (*PhaseDiagram*) – PhaseDiagram object (must be 1-4 components).


    * **show_unstable** (*float*) – Whether unstable (above the hull) phases will be
    plotted. If a number > 0 is entered, all phases with
    e_hull < show_unstable (eV/atom) will be shown.


    * **backend** (*"plotly"** | **"matplotlib"*) – Python package to use for plotting.
    Defaults to “plotly”.


    * **ternary_style** (*"2d"** | **"3d"*) – Ternary phase diagrams are typically plotted in
    two-dimensions (2d), but can be plotted in three dimensions (3d) to visualize
    the depth of the hull. This argument only applies when backend=”plotly”.
    Defaults to “2d”.


    * **\*\*plotkwargs** (*dict*) – Keyword args passed to matplotlib.pyplot.plot (only
    applies when backend=”matplotlib”). Can be used to customize markers
    etc. If not set, the default is:

    > {

    >     “markerfacecolor”: “#4daf4a”,
    >     “markersize”: 10,
    >     “linewidth”: 3

    > }.




#### _create_plotly_element_annotations()
Creates terminal element annotations for Plotly phase diagrams. This method does
not apply to ternary_2d plots.

Functionality is included for phase diagrams with non-elemental endmembers
(as is true for grand potential phase diagrams).


* **Returns**

    List of annotation dicts.



#### _create_plotly_figure_layout(label_stable=True)
Creates layout for plotly phase diagram figure and updates with
figure annotations.


* **Parameters**

    **label_stable** (*bool*) – Whether to label stable compounds



* **Returns**

    Dictionary with Plotly figure layout settings.



#### _create_plotly_fill()
Creates shaded mesh traces for coloring the hull.

For tenrary_3d plots, the color shading is based on formation energy.


* **Returns**

    go.Mesh3d plot



#### _create_plotly_lines()
Create Plotly scatter plots containing line traces of phase diagram facets.


* **Returns**

    Either a go.Scatter (binary), go.Scatterternary (ternary_2d), or
    go.Scatter3d plot (ternary_3d, quaternary)



#### _create_plotly_markers(highlight_entries=None, label_uncertainties=False)
Creates stable and unstable marker plots for overlaying on the phase diagram.


* **Returns**

    Tuple of Plotly go.Scatter (unary, binary), go.Scatterternary(ternary_2d),
    or go.Scatter3d (ternary_3d, quaternary) objects in order:
    (stable markers, unstable markers)



#### _create_plotly_stable_labels(label_stable=True)
Creates a (hidable) scatter trace containing labels of stable phases.
Contains some functionality for creating sensible label positions. This method
does not apply to 2D ternary plots (stable labels are turned off).


* **Returns**

    go.Scatter (or go.Scatter3d) plot



#### _create_plotly_ternary_support_lines()
Creates support lines which aid in seeing the ternary hull in three
dimensions.


* **Returns**

    go.Scatter3d plot of support lines for ternary phase diagram.



#### _create_plotly_uncertainty_shading(stable_marker_plot)
Creates shaded uncertainty region for stable entries. Currently only works
for binary (dim=2) phase diagrams.


* **Parameters**


    * **stable_marker_plot** – go.Scatter object with stable markers and their


    * **bars.** (*error*) –



* **Returns**

    Plotly go.Scatter object with uncertainty window shading.



#### _get_matplotlib_2d_plot(label_stable=True, label_unstable=True, ordering=None, energy_colormap=None, vmin_mev=-60.0, vmax_mev=60.0, show_colorbar=True, process_attributes=False, ax: plt.Axes = None)
Shows the plot using matplotlib.

Imports are done within the function as matplotlib is no longer the default.


#### _get_matplotlib_3d_plot(label_stable=True, ax: plt.Axes = None)
Shows the plot using matplotlib.


* **Parameters**


    * **label_stable** (*bool*) – Whether to label stable compounds.


    * **ax** (*plt.Axes*) – An existing axes object (optional). If not provided, a new one will be created.



* **Returns**

    The axes object with the plot.



* **Return type**

    plt.Axes



#### get_chempot_range_map_plot(elements, referenced=True)
Returns a plot of the chemical potential range _map. Currently works
only for 3-component PDs.

Note: this functionality is now included in the ChemicalPotentialDiagram
class (pymatgen.analysis.chempot_diagram).


* **Parameters**


    * **elements** – Sequence of elements to be considered as independent
    variables. E.g., if you want to show the stability ranges of
    all Li-Co-O phases wrt to uLi and uO, you will supply
    [Element(“Li”), Element(“O”)]


    * **referenced** – if True, gives the results with a reference being the
    energy of the elemental phase. If False, gives absolute values.



* **Returns**

    matplotlib axes object.



* **Return type**

    plt.Axes



#### get_contour_pd_plot()
Plot a contour phase diagram plot, where phase triangles are colored
according to degree of instability by interpolation. Currently only
works for 3-component phase diagrams.


* **Returns**

    A matplotlib plot object.



#### get_plot(label_stable: bool = True, label_unstable: bool = True, ordering: Sequence[str] | None = None, energy_colormap=None, process_attributes: bool = False, ax: plt.Axes = None, label_uncertainties: bool = False, fill: bool = True, highlight_entries: Collection[PDEntry] | None = None)

* **Parameters**


    * **label_stable** – Whether to label stable compounds.


    * **label_unstable** – Whether to label unstable compounds.


    * **ordering** – Ordering of vertices, given as a list [‘Up’,
    ‘Left’,’Right’] (matplotlib only).


    * **energy_colormap** – Colormap for coloring energy (matplotlib only).


    * **process_attributes** – Whether to process the attributes (matplotlib only).


    * **ax** – Existing matplotlib Axes object if plotting multiple phase diagrams
    (matplotlib only).


    * **label_uncertainties** – Whether to add error bars to the hull.
    For binaries, this also shades the hull with the uncertainty window.
    (plotly only).


    * **fill** – Whether to shade the hull. For ternary_2d and quaternary plots, this
    colors facets arbitrarily for visual clarity. For ternary_3d plots, this
    shades the hull by formation energy (plotly only).


    * **highlight_entries** – Entries to highlight in the plot (plotly only). This will
    create a new marker trace that is separate from the other entries.



* **Returns**

    Plotly figure or matplotlib axes object depending on backend.



* **Return type**

    go.Figure | plt.Axes



#### _property_ pd_plot_data()
Plotting data for phase diagram. Cached for repetitive calls.

2-comp - Full hull with energies
3/4-comp - Projection into 2D or 3D Gibbs triangles


* **Returns**


    * lines is a list of list of coordinates for lines in the PD.


    * stable_entries is a dict of {coordinates

        in the phase diagram. (Each coordinate can only have one
        stable phase)


    * unstable_entries is a dict of {entry: coordinates} for all unstable

        nodes in the phase diagram.




* **Return type**

    A tuple containing three objects (lines, stable_entries, unstable_entries)



#### plot_chempot_range_map(elements, referenced=True)
Plot the chemical potential range _map using matplotlib. Currently works only for
3-component PDs. This shows the plot but does not return it.

Note: this functionality is now included in the ChemicalPotentialDiagram
class (pymatgen.analysis.chempot_diagram).


* **Parameters**


    * **elements** – Sequence of elements to be considered as independent
    variables. E.g., if you want to show the stability ranges of
    all Li-Co-O phases wrt to uLi and uO, you will supply
    [Element(“Li”), Element(“O”)]


    * **referenced** – if True, gives the results with a reference being the
    energy of the elemental phase. If False, gives absolute values.



#### plot_element_profile(element, comp, show_label_index=None, xlim=5)
Draw the element profile plot for a composition varying different
chemical potential of an element.

X value is the negative value of the chemical potential reference to
elemental chemical potential. For example, if choose Element(“Li”),
X= -(µLi-µLi0), which corresponds to the voltage versus metal anode.
Y values represent for the number of element uptake in this composition
(unit: per atom). All reactions are printed to help choosing the
profile steps you want to show label in the plot.


* **Parameters**


    * **element** ([*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)) – An element of which the chemical potential is
    considered. It also must be in the phase diagram.


    * **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – A composition.


    * **show_label_index** (*list** of **integers*) – The labels for reaction products
    you want to show in the plot. Default to None (not showing any
    annotation for reaction products). For the profile steps you want
    to show the labels, just add it to the show_label_index. The
    profile step counts from zero. For example, you can set
    show_label_index=[0, 2, 5] to label profile step 0,2,5.


    * **xlim** (*float*) – The max x value. x value is from 0 to xlim. Default to
    5 eV.



* **Returns**

    Plot of element profile evolution by varying the chemical potential
    of an element.



#### show(\*args, \*\*kwargs)
Draw the phase diagram with the provided arguments and display it. This shows
the figure but does not return it.


* **Parameters**


    * **\*args** – Passed to get_plot.


    * **\*\*kwargs** – Passed to get_plot.



#### write_image(stream: str | StringIO, image_format: str = 'svg', \*\*kwargs)
Directly save the plot to a file. This is a wrapper for calling plt.savefig() or
fig.write_image(), depending on the backend. For more customization, it is
recommended to call those methods directly.


* **Parameters**


    * **stream** (*str** | **StringIO*) – Filename or StringIO stream.


    * **image_format** (*str*) – Can be any supported image format for the plotting backend.
    Defaults to ‘svg’ (vector graphics).


    * **\*\*kwargs** – Optinoal kwargs passed to the get_plot function.



### _class_ PatchedPhaseDiagram(entries: Sequence[PDEntry] | set[PDEntry], elements: Sequence[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)] | None = None, keep_all_spaces: bool = False, verbose: bool = False)
Bases: `PhaseDiagram`

Computing the Convex Hull of a large set of data in multiple dimensions is
highly expensive. This class acts to breakdown large chemical spaces into
smaller chemical spaces which can be computed much more quickly due to having
both reduced dimensionality and data set sizes.


### subspaces ({str()
{Element, }}): Dictionary of the sets of elements for each of the
PhaseDiagrams within the PatchedPhaseDiagram.


### pds ({str()
PhaseDiagram}): Dictionary of PhaseDiagrams within the
PatchedPhaseDiagram.


#### all_entries()
All entries provided for Phase Diagram construction.
Note that this does not mean that all these entries are actually used in
the phase diagram. For example, this includes the positive formation energy
entries that are filtered out before Phase Diagram construction.


* **Type**

    list[PDEntry]



#### min_entries()
List of the  lowest energy entries for each composition
in the data provided for Phase Diagram construction.


* **Type**

    list[PDEntry]



#### el_refs()
List of elemental references for the phase diagrams.
These are entries corresponding to the lowest energy element entries for
simple compositional phase diagrams.


* **Type**

    list[PDEntry]



#### elements()
List of elements in the phase diagram.


* **Type**

    list[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)]



* **Parameters**


    * **entries** (*list**[**PDEntry**]*) – A list of PDEntry-like objects having an
    energy, energy_per_atom and composition.


    * **elements** (*list**[*[*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)*]**, **optional*) – Optional list of elements in the phase
    diagram. If set to None, the elements are determined from
    the entries themselves and are sorted alphabetically.
    If specified, element ordering (e.g. for pd coordinates)
    is preserved.


    * **keep_all_spaces** (*bool*) – Boolean control on whether to keep chemical spaces
    that are subspaces of other spaces.


    * **verbose** (*bool*) – Whether to show progress bar during convex hull construction.



#### _get_all_facets_and_simplexes()
Not Implemented - See PhaseDiagram.


#### _get_facet_and_simplex()
Not Implemented - See PhaseDiagram.


#### _get_facet_chempots()
Not Implemented - See PhaseDiagram.


#### _get_pd_patch_for_space(space: frozenset[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)])

* **Parameters**

    **space** (*frozenset**[*[*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)*]*) – chemical space of the form A-B-X.



* **Returns**

    space, PhaseDiagram for the given chemical space



#### _get_simplex_intersections()
Not Implemented - See PhaseDiagram.


#### as_dict()

* **Returns**

    MSONable dictionary representation of PatchedPhaseDiagram.



* **Return type**

    dict[str, Any]



#### _classmethod_ from_dict(dct)

* **Parameters**

    **d** (*dict*) – dictionary representation of PatchedPhaseDiagram.



* **Returns**

    PatchedPhaseDiagram



#### get_all_chempots()
Not Implemented - See PhaseDiagram.


#### get_chempot_range_map()
Not Implemented - See PhaseDiagram.


#### get_chempot_range_stability_phase()
Not Implemented - See PhaseDiagram.


#### get_composition_chempots()
Not Implemented - See PhaseDiagram.


#### get_critical_compositions()
Not Implemented - See PhaseDiagram.


#### get_decomp_and_e_above_hull(entry: PDEntry, allow_negative: bool = False, check_stable: bool = False, on_error: Literal['raise', 'warn', 'ignore'] = 'raise')
Same as method on parent class PhaseDiagram except check_stable defaults to False
for speed. See [https://github.com/materialsproject/pymatgen/issues/2840](https://github.com/materialsproject/pymatgen/issues/2840) for details.


#### get_decomposition(comp: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))
See PhaseDiagram.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – A composition



* **Returns**

    amount} where amount
    is the amount of the fractional composition.



* **Return type**

    Decomposition as a dict of {PDEntry



#### get_element_profile()
Not Implemented - See PhaseDiagram.


#### get_equilibrium_reaction_energy(entry: [Entry](pymatgen.entries.md#pymatgen.entries.Entry))
See PhaseDiagram.

NOTE this is only approximately the same as the what we would get
from PhaseDiagram as we make use of the slsqp approach inside
get_phase_separation_energy().


* **Parameters**

    **entry** (*PDEntry*) – A PDEntry like object



* **Returns**

    Equilibrium reaction energy of entry. Stable entries should have
    equilibrium reaction energy <= 0. The energy is given per atom.



#### get_pd_for_entry(entry: [Entry](pymatgen.entries.md#pymatgen.entries.Entry) | [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))
Get the possible phase diagrams for an entry.


* **Parameters**

    **entry** (*PDEntry** | *[*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – A PDEntry or Composition-like object



* **Returns**

    phase diagram that the entry is part of



* **Return type**

    PhaseDiagram



#### get_transition_chempots()
Not Implemented - See PhaseDiagram.


#### getmu_vertices_stability_phase()
Not Implemented - See PhaseDiagram.


### _class_ PhaseDiagram(entries: Sequence[PDEntry] | set[PDEntry], elements: Sequence[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)] = (), \*, computed_data: dict[str, Any] | None = None)
Bases: `MSONable`

Simple phase diagram class taking in elements and entries as inputs.
The algorithm is based on the work in the following papers:


1.
    1.
        1. Ong, L. Wang, B. Kang, and G. Ceder, Li-Fe-P-O2 Phase Diagram from

> First Principles Calculations. Chem. Mater., 2008, 20(5), 1798-1807.
> doi:10.1021/cm702327g


2.
    1.
        1. Ong, A. Jain, G. Hautier, B. Kang, G. Ceder, Thermal stabilities

> of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first
> principles calculations. Electrochem. Comm., 2010, 12(3), 427-430.
> doi:10.1016/j.elecom.2010.01.010


#### dim()
The dimensionality of the phase diagram.


* **Type**

    int



#### elements()
Elements in the phase diagram.


#### el_refs()
List of elemental references for the phase diagrams. These are
entries corresponding to the lowest energy element entries for simple
compositional phase diagrams.


#### all_entries()
All entries provided for Phase Diagram construction. Note that this
does not mean that all these entries are actually used in the phase
diagram. For example, this includes the positive formation energy
entries that are filtered out before Phase Diagram construction.


#### qhull_entries()
Actual entries used in convex hull. Excludes all positive formation
energy entries.


#### qhull_data()
Data used in the convex hull operation. This is essentially a matrix of
composition data and energy per atom values created from qhull_entries.


#### facets()
Facets of the phase diagram in the form of  [[1,2,3],[4,5,6]…].
For a ternary, it is the indices (references to qhull_entries and
qhull_data) for the vertices of the phase triangles. Similarly
extended to higher D simplices for higher dimensions.


#### simplices()
The simplices of the phase diagram as a list of np.ndarray, i.e.,
the list of stable compositional coordinates in the phase diagram.


* **Parameters**


    * **entries** (*list**[**PDEntry**]*) – A list of PDEntry-like objects having an
    energy, energy_per_atom and composition.


    * **elements** (*list**[*[*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)*]*) – Optional list of elements in the phase
    diagram. If set to None, the elements are determined from
    the entries themselves and are sorted alphabetically.
    If specified, element ordering (e.g. for pd coordinates)
    is preserved.


    * **computed_data** (*dict*) – A dict containing pre-computed data. This allows
    PhaseDiagram object to be reconstituted without performing the
    expensive convex hull computation. The dict is the output from the
    PhaseDiagram._compute() method and is stored in PhaseDiagram.computed_data
    when generated for the first time.



#### _compute()

#### _get_all_facets_and_simplexes(comp)
Get all facets that a composition falls into.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – A composition



#### _get_facet_and_simplex(comp: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))
Get any facet that a composition falls into. Cached so successive
calls at same composition are fast.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – A composition



#### _get_facet_chempots(facet)
Calculates the chemical potentials for each element within a facet.


* **Parameters**

    **facet** – Facet of the phase diagram.



* **Returns**

    chempot} for all elements in the phase diagram.



* **Return type**

    {element



#### _get_simplex_intersections(c1, c2)
Returns coordinates of the intersection of the tie line between two compositions
and the simplexes of the PhaseDiagram.


* **Parameters**


    * **c1** – Reduced dimension coordinates of first composition


    * **c2** – Reduced dimension coordinates of second composition



* **Returns**

    Array of the intersections between the tie line and the simplexes of
    the PhaseDiagram



#### _get_stable_entries_in_space(space)

* **Parameters**

    **space** (*set**[*[*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)*]*) – set of Element objects.



* **Returns**

    stable entries in the space.



* **Return type**

    list[[Entry](pymatgen.entries.md#pymatgen.entries.Entry)]



#### _property_ all_entries_hulldata()
Returns:
The actual ndarray used to construct the convex hull.


#### as_dict()

* **Returns**

    MSONable dictionary representation of PhaseDiagram.



#### formation_energy_tol(_ = 1e-1_ )

#### _classmethod_ from_dict(dct: dict[str, Any])

* **Parameters**

    **d** (*dict*) – dictionary representation of PhaseDiagram.



* **Returns**

    PhaseDiagram



#### get_all_chempots(comp)
Get chemical potentials at a given composition.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Composition



* **Returns**

    Chemical potentials.



#### get_chempot_range_map(elements: Sequence[[Element](pymatgen.core.md#pymatgen.core.periodic_table.Element)], referenced: bool = True, joggle: bool = True)
Returns a chemical potential range map for each stable entry.


* **Parameters**


    * **elements** – Sequence of elements to be considered as independent variables.
    E.g., if you want to show the stability ranges
    of all Li-Co-O phases with respect to mu_Li and mu_O, you will supply
    [Element(“Li”), Element(“O”)]


    * **referenced** – If True, gives the results with a reference being the
    energy of the elemental phase. If False, gives absolute values.


    * **joggle** (*bool*) – Whether to joggle the input to avoid precision
    errors.



* **Returns**

    [simplices]}. The list of
    simplices are the sides of the N-1 dim polytope bounding the
    allowable chemical potential range of each entry.



* **Return type**

    Returns a dict of the form {entry



#### get_chempot_range_stability_phase(target_comp, open_elt)
Returns a set of chemical potentials corresponding to the max and min
chemical potential of the open element for a given composition. It is
quite common to have for instance a ternary oxide (e.g., ABO3) for
which you want to know what are the A and B chemical potential leading
to the highest and lowest oxygen chemical potential (reducing and
oxidizing conditions). This is useful for defect computations.


* **Parameters**


    * **target_comp** – A Composition object


    * **open_elt** – Element that you want to constrain to be max or min



* **Returns**

    (mu_min, mu_max)}: Chemical potentials are given in

        ”absolute” values (i.e., not referenced to 0)




* **Return type**

    {Element



#### get_composition_chempots(comp)
Get the chemical potentials for all elements at a given composition.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Composition



* **Returns**

    Dictionary of chemical potentials.



#### get_critical_compositions(comp1, comp2)
Get the critical compositions along the tieline between two
compositions. I.e. where the decomposition products change.
The endpoints are also returned.


* **Parameters**


    * **comp1** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – First composition to define the tieline


    * **comp2** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Second composition to define the tieline



* **Returns**

    list of critical compositions. All are of

        the form x \* comp1 + (1-x) \* comp2




* **Return type**

    [([Composition](pymatgen.core.md#pymatgen.core.composition.Composition))]



#### get_decomp_and_e_above_hull(entry: PDEntry, allow_negative: bool = False, check_stable: bool = True, on_error: Literal['raise', 'warn', 'ignore'] = 'raise')
Provides the decomposition and energy above convex hull for an entry.
Due to caching, can be much faster if entries with the same composition
are processed together.


* **Parameters**


    * **entry** (*PDEntry*) – A PDEntry like object


    * **allow_negative** (*bool*) – Whether to allow negative e_above_hulls. Used to
    calculate equilibrium reaction energies. Defaults to False.


    * **check_stable** (*bool*) – Whether to first check whether an entry is stable.
    In normal circumstances, this is the faster option since checking for
    stable entries is relatively fast. However, if you have a huge proportion
    of unstable entries, then this check can slow things down. You should then
    set this to False.


    * **on_error** (*'raise'** | **'warn'** | **'ignore'*) – What to do if no valid decomposition was
    found. ‘raise’ will throw ValueError. ‘warn’ will print return (None, None).
    ‘ignore’ just returns (None, None). Defaults to ‘raise’.



* **Raises**

    **ValueError** – If no valid decomposition exists in this phase diagram for given entry.



* **Returns**

    (decomp, energy_above_hull). The decomposition is provided

        as a dict of {PDEntry: amount} where amount is the amount of the
        fractional composition. Stable entries should have energy above
        convex hull of 0. The energy is given per atom.




#### get_decomp_and_hull_energy_per_atom(comp: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))

* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Input composition.



* **Returns**

    Energy of lowest energy equilibrium at desired composition per atom



#### get_decomp_and_phase_separation_energy(entry: PDEntry, space_limit: int = 200, stable_only: bool = False, tols: Sequence[float] = (1e-08,), maxiter: int = 1000, \*\*kwargs: Any)
Provides the combination of entries in the PhaseDiagram that gives the
lowest formation enthalpy with the same composition as the given entry
excluding entries with the same composition and the energy difference
per atom between the given entry and the energy of the combination found.

For unstable entries that are not polymorphs of stable entries (or completely
novel entries) this is simply the energy above (or below) the convex hull.

For entries with the same composition as one of the stable entries in the
phase diagram setting stable_only to False (Default) allows for entries
not previously on the convex hull to be considered in the combination.
In this case the energy returned is what is referred to as the decomposition
enthalpy in:


1. Bartel, C., Trewartha, A., Wang, Q., Dunn, A., Jain, A., Ceder, G.,

    A critical examination of compound stability predictions from
    machine-learned formation energies, npj Computational Materials 6, 97 (2020)

For stable entries setting stable_only to True returns the same energy
as get_equilibrium_reaction_energy. This function is based on a constrained
optimization rather than recalculation of the convex hull making it
algorithmically cheaper. However, if tol is too loose there is potential
for this algorithm to converge to a different solution.


* **Parameters**


    * **entry** (*PDEntry*) – A PDEntry like object.


    * **space_limit** (*int*) – The maximum number of competing entries to consider
    before calculating a second convex hull to reducing the complexity
    of the optimization.


    * **stable_only** (*bool*) – Only use stable materials as competing entries.


    * **tols** (*list**[**float**]*) – Tolerances for convergence of the SLSQP optimization
    when finding the equilibrium reaction. Tighter tolerances tested first.


    * **maxiter** (*int*) – The maximum number of iterations of the SLSQP optimizer
    when finding the equilibrium reaction.


    * **\*\*kwargs** – Passed to get_decomp_and_e_above_hull.



* **Returns**

    (decomp, energy). The decomposition  is given as a dict of {PDEntry, amount}
    for all entries in the decomp reaction where amount is the amount of the
    fractional composition. The phase separation energy is given per atom.



#### get_decomposition(comp: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))
Provides the decomposition at a particular composition.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – A composition



* **Returns**

    amount} where amount
    is the amount of the fractional composition.



* **Return type**

    Decomposition as a dict of {PDEntry



#### get_e_above_hull(entry: PDEntry, \*\*kwargs: Any)
Provides the energy above convex hull for an entry.


* **Parameters**


    * **entry** (*PDEntry*) – A PDEntry like object.


    * **\*\*kwargs** – Passed to get_decomp_and_e_above_hull().



* **Returns**

    Energy above convex hull of entry. Stable entries should have

        energy above hull of 0. The energy is given per atom.




* **Return type**

    float | None



#### get_element_profile(element, comp, comp_tol=1e-05)
Provides the element evolution data for a composition. For example, can be used
to analyze Li conversion voltages by varying mu_Li and looking at the phases
formed. Also can be used to analyze O2 evolution by varying mu_O2.


* **Parameters**


    * **element** – An element. Must be in the phase diagram.


    * **comp** – A Composition


    * **comp_tol** – The tolerance to use when calculating decompositions.
    Phases with amounts less than this tolerance are excluded.
    Defaults to 1e-5.



* **Returns**

    [ {‘chempot’: -10.487582010000001, ‘evolution’: -2.0,
    ‘reaction’: Reaction Object], …]



* **Return type**

    Evolution data as a list of dictionaries of the following format



#### get_equilibrium_reaction_energy(entry: PDEntry)
Provides the reaction energy of a stable entry from the neighboring
equilibrium stable entries (also known as the inverse distance to
hull).


* **Parameters**

    **entry** (*PDEntry*) – A PDEntry like object



* **Returns**

    Equilibrium reaction energy of entry. Stable entries should have

        equilibrium reaction energy <= 0. The energy is given per atom.




* **Return type**

    float | None



#### get_form_energy(entry: PDEntry)
Returns the formation energy for an entry (NOT normalized) from the
elemental references.


* **Parameters**

    **entry** (*PDEntry*) – A PDEntry-like object.



* **Returns**

    Formation energy from the elemental references.



* **Return type**

    float



#### get_form_energy_per_atom(entry: PDEntry)
Returns the formation energy per atom for an entry from the
elemental references.


* **Parameters**

    **entry** (*PDEntry*) – An PDEntry-like object



* **Returns**

    Formation energy **per atom** from the elemental references.



#### get_hull_energy(comp: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))

* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Input composition.



* **Returns**

    Energy of lowest energy equilibrium at desired composition. Not

        normalized by atoms, i.e. E(Li4O2) = 2 \* E(Li2O)




#### get_hull_energy_per_atom(comp: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition), \*\*kwargs)

* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Input composition.



* **Returns**

    Energy of lowest energy equilibrium at desired composition.



#### get_phase_separation_energy(entry, \*\*kwargs)
Provides the energy to the convex hull for the given entry. For stable entries
already in the phase diagram the algorithm provides the phase separation energy
which is referred to as the decomposition enthalpy in:


1. Bartel, C., Trewartha, A., Wang, Q., Dunn, A., Jain, A., Ceder, G.,

    A critical examination of compound stability predictions from
    machine-learned formation energies, npj Computational Materials 6, 97 (2020)


* **Parameters**


    * **entry** (*PDEntry*) – A PDEntry like object


    * **\*\*kwargs** – Keyword args passed to get_decomp_and_decomp_energy
    space_limit (int): The maximum number of competing entries to consider.
    stable_only (bool): Only use stable materials as competing entries
    tol (float): The tolerance for convergence of the SLSQP optimization

    > when finding the equilibrium reaction.

    maxiter (int): The maximum number of iterations of the SLSQP optimizer

        when finding the equilibrium reaction.




* **Returns**

    phase separation energy per atom of entry. Stable entries should have
    energies <= 0, Stable elemental entries should have energies = 0 and
    unstable entries should have energies > 0. Entries that have the same
    composition as a stable energy may have positive or negative phase
    separation energies depending on their own energy.



#### get_plot(show_unstable: float = 0.2, backend: Literal['plotly', 'matplotlib'] = 'plotly', ternary_style: Literal['2d', '3d'] = '2d', label_stable: bool = True, label_unstable: bool = True, ordering: Sequence[str] | None = None, energy_colormap=None, process_attributes: bool = False, ax: plt.Axes = None, label_uncertainties: bool = False, fill: bool = True, \*\*kwargs)
Convenient wrapper for PDPlotter. Initializes a PDPlotter object and calls
get_plot() with provided combined arguments.

Plotting is only supported for phase diagrams with <=4 elements (unary,
binary, ternary, or quaternary systems).


* **Parameters**


    * **show_unstable** (*float*) – Whether unstable (above the hull) phases will be
    plotted. If a number > 0 is entered, all phases with
    e_hull < show_unstable (eV/atom) will be shown.


    * **backend** (*"plotly"** | **"matplotlib"*) – Python package to use for plotting.
    Defaults to “plotly”.


    * **ternary_style** (*"2d"** | **"3d"*) – Ternary phase diagrams are typically plotted in
    two-dimensions (2d), but can be plotted in three dimensions (3d) to visualize
    the depth of the hull. This argument only applies when backend=”plotly”.
    Defaults to “2d”.


    * **label_stable** – Whether to label stable compounds.


    * **label_unstable** – Whether to label unstable compounds.


    * **ordering** – Ordering of vertices (matplotlib backend only).


    * **energy_colormap** – Colormap for coloring energy (matplotlib backend only).


    * **process_attributes** – Whether to process the attributes (matplotlib
    backend only).


    * **ax** – Existing Axes object if plotting multiple phase diagrams (matplotlib backend only).


    * **label_uncertainties** – Whether to add error bars to the hull (plotly
    backend only). For binaries, this also shades the hull with the
    uncertainty window.


    * **fill** – Whether to shade the hull. For ternary_2d and quaternary plots, this
    colors facets arbitrarily for visual clarity. For ternary_3d plots, this
    shades the hull by formation energy (plotly backend only).


    * **\*\*kwargs** (*dict*) – Keyword args passed to PDPlotter.get_plot(). Can be used to customize markers
    etc. If not set, the default is { “markerfacecolor”: “#4daf4a”, “markersize”: 10, “linewidth”: 3 }



#### get_reference_energy(comp: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))
Sum of elemental reference energies over all elements in a composition.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Input composition.



* **Returns**

    Reference energy



* **Return type**

    float



#### get_reference_energy_per_atom(comp: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))
Sum of elemental reference energies over all elements in a composition.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Input composition.



* **Returns**

    Reference energy per atom



* **Return type**

    float



#### get_transition_chempots(element)
Get the critical chemical potentials for an element in the Phase
Diagram.


* **Parameters**

    **element** – An element. Has to be in the PD in the first place.



* **Returns**

    A sorted sequence of critical chemical potentials, from less
    negative to more negative.



#### getmu_vertices_stability_phase(target_comp, dep_elt, tol_en=0.01)
Returns a set of chemical potentials corresponding to the vertices of
the simplex in the chemical potential phase diagram.
The simplex is built using all elements in the target_composition
except dep_elt.
The chemical potential of dep_elt is computed from the target
composition energy.
This method is useful to get the limiting conditions for
defects computations for instance.


* **Parameters**


    * **target_comp** – A Composition object


    * **dep_elt** – the element for which the chemical potential is computed
    from the energy of the stable phase at the target composition


    * **tol_en** – a tolerance on the energy to set



* **Returns**

    mu}]: An array of conditions on simplex vertices for
    which each element has a chemical potential set to a given
    value. “absolute” values (i.e., not referenced to element energies)



* **Return type**

    [{Element



#### numerical_tol(_ = 1e-0_ )

#### pd_coords(comp: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))
The phase diagram is generated in a reduced dimensional space
(n_elements - 1). This function returns the coordinates in that space.
These coordinates are compatible with the stored simplex objects.


* **Parameters**

    **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – A composition



* **Returns**

    The coordinates for a given composition in the PhaseDiagram’s basis



#### _property_ stable_entries(_: set[[Entry](pymatgen.entries.md#pymatgen.entries.Entry)_ )
Returns:
set[Entry]: of stable entries in the phase diagram.


#### _property_ unstable_entries(_: set[[Entry](pymatgen.entries.md#pymatgen.entries.Entry)_ )
Returns:
set[Entry]: unstable entries in the phase diagram. Includes positive formation energy entries.


### _exception_ PhaseDiagramError()
Bases: `Exception`

An exception class for Phase Diagram generation.


### _class_ ReactionDiagram(entry1, entry2, all_entries, tol: float = 0.0001, float_fmt='%.4f')
Bases: `object`

Analyzes the possible reactions between a pair of compounds, e.g.,
an electrolyte and an electrode.


* **Parameters**


    * **entry1** ([*ComputedEntry*](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)) – Entry for 1st component. Note that
    corrections, if any, must already be pre-applied. This is to
    give flexibility for different kinds of corrections, e.g.,
    if a particular entry is fitted to an experimental data (such
    as EC molecule).


    * **entry2** ([*ComputedEntry*](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)) – Entry for 2nd component. Note that
    corrections must already be pre-applied. This is to
    give flexibility for different kinds of corrections, e.g.,
    if a particular entry is fitted to an experimental data (such
    as EC molecule).


    * **all_entries** (*[*[*ComputedEntry*](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)*]*) – All other entries to be
    considered in the analysis. Note that corrections, if any,
    must already be pre-applied.


    * **tol** (*float*) – Tolerance to be used to determine validity of reaction.


    * **float_fmt** (*str*) – Formatting string to be applied to all floats.
    Determines number of decimal places in reaction string.



#### get_compound_pd()
Get the CompoundPhaseDiagram object, which can then be used for
plotting.


* **Returns**

    CompoundPhaseDiagram



### _class_ TransformedPDEntry(entry, sp_mapping, name=None)
Bases: `PDEntry`

This class represents a TransformedPDEntry, which allows for a PDEntry to be
transformed to a different composition coordinate space. It is used in the
construction of phase diagrams that do not have elements as the terminal
compositions.


* **Parameters**


    * **entry** (*PDEntry*) – Original entry to be transformed.


    * **(****{Composition** (*sp_mapping*) – DummySpecies}): dictionary mapping Terminal Compositions to Dummy Species.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### amount_tol(_ = 1e-0_ )

#### as_dict()

* **Returns**

    MSONable dictionary representation of TransformedPDEntry.



#### _property_ composition(_: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition_ )
The composition in the dummy species space.


* **Returns**

    Composition



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – dictionary representation of TransformedPDEntry.



* **Returns**

    TransformedPDEntry



### _exception_ TransformedPDEntryError()
Bases: `Exception`

An exception class for TransformedPDEntry.


### _get_slsqp_decomp(comp, competing_entries, tols=(1e-08,), maxiter=1000)
Finds the amounts of competing compositions that minimize the energy of a
given composition.

The algorithm is based on the work in the following paper:


1. Bartel, C., Trewartha, A., Wang, Q., Dunn, A., Jain, A., Ceder, G.,

    A critical examination of compound stability predictions from
    machine-learned formation energies, npj Computational Materials 6, 97 (2020)


* **Parameters**


    * **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – A Composition to analyze


    * **competing_entries** (*[**PDEntry**]*) – List of entries to consider for decomposition


    * **tols** (*list*) – tolerances to try for SLSQP convergence. Issues observed for
    tol > 1e-7 in the fractional composition (default 1e-8)


    * **maxiter** (*int*) – maximum number of SLSQP iterations



* **Returns**

    amount} where amount
    is the amount of the fractional composition.



* **Return type**

    decomposition as a dict of {PDEntry



### get_facets(qhull_data: ArrayLike, joggle: bool = False)
Get the simplex facets for the Convex hull.


* **Parameters**


    * **qhull_data** (*np.ndarray*) – The data from which to construct the convex
    hull as a Nxd array (N being number of data points and d being the
    dimension)


    * **joggle** (*bool*) – Whether to joggle the input to avoid precision
    errors.



* **Returns**

    List of simplices of the Convex Hull.



### order_phase_diagram(lines, stable_entries, unstable_entries, ordering)
Orders the entries (their coordinates) in a phase diagram plot according
to the user specified ordering.
Ordering should be given as [‘Up’, ‘Left’, ‘Right’], where Up,
Left and Right are the names of the entries in the upper, left and right
corners of the triangle respectively.


* **Parameters**


    * **lines** – list of list of coordinates for lines in the PD.


    * **stable_entries** – {coordinate : entry} for each stable node in the
    phase diagram. (Each coordinate can only have one stable phase)


    * **unstable_entries** – {entry: coordinates} for all unstable nodes in the
    phase diagram.


    * **ordering** – Ordering of the phase diagram, given as a list [‘Up’,
    ‘Left’,’Right’]



* **Returns**


    * newlines is a list of list of coordinates for lines in the PD.


    * newstable_entries is a {coordinate : entry} for each stable node

    in the phase diagram. (Each coordinate can only have one
    stable phase)
    - newunstable_entries is a {entry: coordinates} for all unstable
    nodes in the phase diagram.




* **Return type**

    (newlines, newstable_entries, newunstable_entries)



### tet_coord(coord)
Convert a 3D coordinate into a tetrahedron based coordinate system for a
prettier phase diagram.


* **Parameters**

    **coord** – coordinate used in the convex hull computation.



* **Returns**

    coordinates in a tetrahedron-based coordinate system.



### triangular_coord(coord)
Convert a 2D coordinate into a triangle-based coordinate system for a
prettier phase diagram.


* **Parameters**

    **coord** – coordinate used in the convex hull computation.



* **Returns**

    coordinates in a triangular-based coordinate system.



### uniquelines(q)
Given all the facets, convert it into a set of unique lines. Specifically
used for converting convex hull facets into line pairs of coordinates.


* **Parameters**

    **q** – A 2-dim sequence, where each row represents a facet. E.g.,
    [[1,2,3],[3,6,7],…]



* **Returns**

    A set of tuple of lines. E.g., ((1,2), (1,3), (2,3), ….)



* **Return type**

    setoflines


## pymatgen.analysis.piezo module

This module provides classes for the Piezoelectric tensor.


### _class_ PiezoTensor(input_array, tol: float = 0.001)
Bases: [`Tensor`](pymatgen.core.md#pymatgen.core.tensors.Tensor)

This class describes the 3x6 piezo tensor in Voigt-notation.

Create an PiezoTensor object. The constructor throws an error if
the shape of the input_matrix argument is not 3x3x3, i. e. in true
tensor notation. Note that the constructor uses __new__ rather than
__init__ according to the standard method of subclassing numpy
ndarrays.


* **Parameters**

    **input_matrix** (*3x3x3 array-like*) – the 3x6 array-like
    representing the piezo tensor



#### _classmethod_ from_vasp_voigt(input_vasp_array)

* **Parameters**

    **input_vasp_array** (*nd.array*) – Voigt form of tensor.



* **Returns**

    PiezoTensor


## pymatgen.analysis.piezo_sensitivity module

Piezo sensitivity analysis module.


### _class_ BornEffectiveCharge(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), bec, pointops, tol: float = 0.001)
Bases: `object`

This class describes the Nx3x3 born effective charge tensor.

Create an BornEffectiveChargeTensor object defined by a
structure, point operations of the structure’s atomic sites.
Note that the constructor uses __new__ rather than __init__
according to the standard method of subclassing numpy ndarrays.


* **Parameters**

    **input_matrix** (*Nx3x3 array-like*) – the Nx3x3 array-like
    representing the born effective charge tensor



#### get_BEC_operations(eigtol=1e-05, opstol=0.001)
Returns the symmetry operations which maps the tensors
belonging to equivalent sites onto each other in the form
[site index 1, site index 2, [Symmops mapping from site
index 1 to site index 2]].


* **Parameters**


    * **eigtol** (*float*) – tolerance for determining if two sites are


    * **symmetry** (*related by*) –


    * **opstol** (*float*) – tolerance for determining if a symmetry


    * **sites** (*operation relates two*) –



* **Returns**

    list of symmetry operations mapping equivalent sites and
    the indexes of those sites.



#### get_rand_BEC(max_charge=1)
Generate a random born effective charge tensor which obeys a structure’s
symmetry and the acoustic sum rule.


* **Parameters**

    **max_charge** (*float*) – maximum born effective charge value



* **Returns**

    np.array Born effective charge tensor



### _class_ ForceConstantMatrix(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), fcm, pointops, sharedops, tol: float = 0.001)
Bases: `object`

This class describes the NxNx3x3 force constant matrix defined by a
structure, point operations of the structure’s atomic sites, and the
shared symmetry operations between pairs of atomic sites.

Create an ForceConstantMatrix object.


* **Parameters**

    **input_matrix** (*NxNx3x3 array-like*) – the NxNx3x3 array-like
    representing the force constant matrix



#### get_FCM_operations(eigtol=1e-05, opstol=1e-05)
Returns the symmetry operations which maps the tensors
belonging to equivalent sites onto each other in the form
[site index 1a, site index 1b, site index 2a, site index 2b,
[Symmops mapping from site index 1a, 1b to site index 2a, 2b]].


* **Parameters**


    * **eigtol** (*float*) – tolerance for determining if two sites are


    * **symmetry** (*related by*) –


    * **opstol** (*float*) – tolerance for determining if a symmetry


    * **sites** (*operation relates two*) –



* **Returns**

    list of symmetry operations mapping equivalent sites and
    the indexes of those sites.



#### get_asum_FCM(fcm: ndarray, numiter: int = 15)
Generate a symmeterized force constant matrix that obeys the objects symmetry
constraints and obeys the acoustic sum rule through an iterative procedure.


* **Parameters**


    * **fcm** (*numpy array*) – 3Nx3N unsymmeterized force constant matrix


    * **numiter** (*int*) – number of iterations to attempt to obey the acoustic sum
    rule



* **Returns**

    numpy array representing the force constant matrix



#### get_rand_FCM(asum=15, force=10)
Generate a symmeterized force constant matrix from an unsymmeterized matrix
that has no unstable modes and also obeys the acoustic sum rule through an
iterative procedure.


* **Parameters**


    * **force** (*float*) – maximum force constant


    * **asum** (*int*) – number of iterations to attempt to obey the acoustic sum
    rule



* **Returns**

    NxNx3x3 np.array representing the force constant matrix



#### get_stable_FCM(fcm, fcmasum=10)
Generate a symmeterized force constant matrix that obeys the objects symmetry
constraints, has no unstable modes and also obeys the acoustic sum rule through an
iterative procedure.


* **Parameters**


    * **fcm** (*numpy array*) – unsymmeterized force constant matrix


    * **fcmasum** (*int*) – number of iterations to attempt to obey the acoustic sum
    rule



* **Returns**

    3Nx3N numpy array representing the force constant matrix



#### get_symmetrized_FCM(unsymmetrized_fcm, max_force=1)
Generate a symmeterized force constant matrix from an unsymmeterized matrix.


* **Parameters**


    * **unsymmetrized_fcm** (*numpy array*) – unsymmeterized force constant matrix


    * **max_charge** (*float*) – maximum born effective charge value



* **Returns**

    3Nx3N numpy array representing the force constant matrix



#### get_unstable_FCM(max_force=1)
Generate an unsymmeterized force constant matrix.


* **Parameters**

    **max_charge** (*float*) – maximum born effective charge value



* **Returns**

    numpy array representing the force constant matrix



### _class_ InternalStrainTensor(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), ist, pointops, tol: float = 0.001)
Bases: `object`

This class describes the Nx3x3x3 internal tensor defined by a
structure, point operations of the structure’s atomic sites.

Create an InternalStrainTensor object.


* **Parameters**

    **input_matrix** (*Nx3x3x3 array-like*) – the Nx3x3x3 array-like
    representing the internal strain tensor



#### get_IST_operations(opstol=0.001)
Returns the symmetry operations which maps the tensors
belonging to equivalent sites onto each other in the form
[site index 1, site index 2, [Symmops mapping from site
index 1 to site index 2]].


* **Parameters**


    * **opstol** (*float*) – tolerance for determining if a symmetry


    * **sites** (*operation relates two*) –



* **Returns**

    list of symmetry operations mapping equivalent sites and
    the indexes of those sites.



#### get_rand_IST(max_force=1)
Generate a random internal strain tensor which obeys a structure’s
symmetry and the acoustic sum rule.


* **Parameters**

    **max_force** (*float*) – maximum born effective charge value



* **Returns**

    InternalStrainTensor object



### get_piezo(BEC, IST, FCM, rcond=0.0001)
Generate a random piezoelectric tensor based on a structure and corresponding
symmetry.


* **Parameters**


    * **BEC** (*numpy array*) – Nx3x3 array representing the born effective charge tensor


    * **IST** (*numpy array*) – Nx3x3x3 array representing the internal strain tensor


    * **FCM** (*numpy array*) – NxNx3x3 array representing the born effective charge tensor


    * **rcondy** (*float*) – condition for excluding eigenvalues in the pseudoinverse



* **Returns**

    3x3x3 calculated Piezo tensor



### rand_piezo(struct, pointops, sharedops, BEC, IST, FCM, anumiter=10)
Generate a random piezoelectric tensor based on a structure and corresponding
symmetry.


* **Parameters**


    * **struct** (*pymatgen structure*) – structure whose symmetry operations the piezo tensor must obey


    * **pointops** – list of point operations obeyed by a single atomic site


    * **sharedops** – list of point operations shared by a pair of atomic sites


    * **BEC** (*numpy array*) – Nx3x3 array representing the born effective charge tensor


    * **IST** (*numpy array*) – Nx3x3x3 array representing the internal strain tensor


    * **FCM** (*numpy array*) – NxNx3x3 array representing the born effective charge tensor


    * **anumiter** (*int*) – number of iterations for acoustic sum rule convergence



* **Returns**

    list in the form of [Nx3x3 random born effective charge tenosr,
    Nx3x3x3 random internal strain tensor, NxNx3x3 random force constant matrix, 3x3x3 piezo tensor]


## pymatgen.analysis.pourbaix_diagram module

This module is intended to be used to compute Pourbaix diagrams of arbitrary compositions
and formation energies.


### _class_ IonEntry(ion: [Ion](pymatgen.core.md#pymatgen.core.ion.Ion), energy: float, name: str | None = None, attribute=None)
Bases: `PDEntry`

Object similar to PDEntry, but contains an Ion object instead of a
Composition object.


#### name()
A name for the entry. This is the string shown in the phase diagrams.
By default, this is the reduced formula for the composition, but can be
set to some other string for display purposes.


* **Type**

    str



* **Parameters**


    * **ion** – Ion object


    * **energy** – Energy for composition.


    * **name** – Optional parameter to name the entry. Defaults to the
    chemical formula.


    * **attribute** – Optional attribute of the entry, e.g., band gap.



#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
Creates a dict of composition, energy, and ion name.


#### _classmethod_ from_dict(d)
Returns an IonEntry object from a dict.


### _class_ MultiEntry(entry_list, weights=None)
Bases: `PourbaixEntry`

PourbaixEntry-like object for constructing multi-elemental Pourbaix diagrams.

Initializes a MultiEntry.


* **Parameters**


    * **entry_list** (*[**PourbaixEntry**]*) – List of component PourbaixEntries


    * **weights** (*[**float**]*) – Weights associated with each entry. Default is None



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(dct)

* **Parameters**

    **dct** (*dict*) – Dict representation.



* **Returns**

    MultiEntry



#### _property_ name()
MultiEntry name, i. e. the name of each entry joined by ‘ + ‘.


### _class_ PourbaixDiagram(entries: list[PourbaixEntry] | list[MultiEntry], comp_dict: dict[str, float] | None = None, conc_dict: dict[str, float] | None = None, filter_solids: bool = True, nproc: int | None = None)
Bases: `MSONable`

Class to create a Pourbaix diagram from entries.


* **Parameters**


    * **entries** (*[**PourbaixEntry**] or **[**MultiEntry**]*) – Entries list
    containing Solids and Ions or a list of MultiEntries


    * **comp_dict** (*dict**[**str**, **float**]*) – Dictionary of compositions,
    defaults to equal parts of each elements


    * **conc_dict** (*dict**[**str**, **float**]*) – Dictionary of ion concentrations,
    defaults to 1e-6 for each element


    * **filter_solids** (*bool*) – applying this filter to a Pourbaix
    diagram ensures all included solid phases are filtered by
    stability on the compositional phase diagram. Defaults to True.
    The practical consequence of this is that highly oxidized or reduced
    phases that might show up in experiments due to kinetic limitations
    on oxygen/hydrogen evolution won’t appear in the diagram, but they are
    not actually “stable” (and are frequently overstabilized from DFT errors).
    Hence, including only the stable solid phases generally leads to the
    most accurate Pourbaix diagrams.


    * **nproc** (*int*) – number of processes to generate multi-entries with
    in parallel. Defaults to None (serial processing).



#### _convert_entries_to_points(pourbaix_entries)

* **Parameters**

    **pourbaix_entries** (*[**PourbaixEntry**]*) – list of Pourbaix entries
    to process into vectors in nph-nphi-composition space.



* **Returns**

    list of vectors, [[nph, nphi, e0, x1, x2, …, xn-1]]
    corresponding to each entry in nph-nphi-composition space



#### _generate_multielement_entries(entries, nproc=None)
Create entries for multi-element Pourbaix construction.

This works by finding all possible linear combinations
of entries that can result in the specified composition
from the initialized comp_dict.


* **Parameters**


    * **entries** (*[**PourbaixEntries**]*) – list of Pourbaix entries
    to process into MultiEntries


    * **nproc** (*int*) – number of processes to be used in parallel
    treatment of entry combos



#### _get_hull_in_nph_nphi_space(entries)
Generates convex hull of Pourbaix diagram entries in composition,
npH, and nphi space. This enables filtering of multi-entries
such that only compositionally stable combinations of entries
are included.


* **Parameters**

    **entries** (*[**PourbaixEntry**]*) – list of PourbaixEntries to construct
    the convex hull



* **Returns**

    PourbaixEntry list and stable

        facets corresponding to that list




* **Return type**

    tuple[list[PourbaixEntry], list[[Simplex](pymatgen.util.md#pymatgen.util.coord.Simplex)]]



#### _preprocess_pourbaix_entries(entries, nproc=None)
Generates multi-entries for Pourbaix diagram.


* **Parameters**


    * **entries** (*[**PourbaixEntry**]*) – list of PourbaixEntries to preprocess
    into MultiEntries


    * **nproc** (*int*) – number of processes to be used in parallel
    treatment of entry combos



* **Returns**

    ([MultiEntry]) list of stable MultiEntry candidates



#### _property_ all_entries()
Return all entries used to generate the Pourbaix diagram.


#### as_dict()

* **Returns**

    MSONable dict.



#### find_stable_entry(pH, V)
Finds stable entry at a pH,V condition


* **Parameters**


    * **pH** (*float*) – pH to find stable entry


    * **V** (*float*) – V to find stable entry.


Returns:


#### _classmethod_ from_dict(d)

* **Parameters**

    **(****)** (*d*) – Dict representation.



* **Returns**

    PourbaixDiagram



#### get_decomposition_energy(entry, pH, V)
Finds decomposition to most stable entries in eV/atom,
supports vectorized inputs for pH and V.


* **Parameters**


    * **entry** (*PourbaixEntry*) – PourbaixEntry corresponding to
    compound to find the decomposition for


    * **pH** (*float**, **[**float**]*) – pH at which to find the decomposition


    * **V** (*float**, **[**float**]*) – voltage at which to find the decomposition



* **Returns**

    Decomposition energy for the entry, i. e. the energy above

        the “Pourbaix hull” in eV/atom at the given conditions




#### get_hull_energy(pH, V)
Gets the minimum energy of the Pourbaix “basin” that is formed
from the stable Pourbaix planes. Vectorized.


* **Parameters**


    * **pH** (*float** or **[**float**]*) – pH at which to find the hull energy


    * **V** (*float** or **[**float**]*) – V at which to find the hull energy



* **Returns**

    (float or [float]) minimum Pourbaix energy at conditions



#### _static_ get_pourbaix_domains(pourbaix_entries, limits=None)
Returns a set of Pourbaix stable domains (i. e. polygons) in
pH-V space from a list of pourbaix_entries.

This function works by using scipy’s HalfspaceIntersection
function to construct all of the 2-D polygons that form the
boundaries of the planes corresponding to individual entry
gibbs free energies as a function of pH and V. Hyperplanes
of the form a\*pH + b\*V + 1 - g(0, 0) are constructed and
supplied to HalfspaceIntersection, which then finds the
boundaries of each Pourbaix region using the intersection
points.


* **Parameters**


    * **pourbaix_entries** (*[**PourbaixEntry**]*) – Pourbaix entries
    with which to construct stable Pourbaix domains


    * **limits** (*[**[**float**]**]*) – limits in which to do the pourbaix
    analysis



* **Returns**

    [boundary_points]}.
    The list of boundary points are the sides of the N-1
    dim polytope bounding the allowable ph-V range of each entry.



* **Return type**

    Returns a dict of the form {entry



#### get_stable_entry(pH, V)
Gets the stable entry at a given pH, V condition.


* **Parameters**


    * **pH** (*float*) – pH at a given condition


    * **V** (*float*) – V at a given condition



* **Returns**

    Pourbaix or multi-entry

        corresponding ot the minimum energy entry at a given
        pH, V condition




* **Return type**

    (PourbaixEntry or MultiEntry)



#### _static_ process_multientry(entry_list, prod_comp, coeff_threshold=0.0001)
Static method for finding a multientry based on
a list of entries and a product composition.
Essentially checks to see if a valid aqueous
reaction exists between the entries and the
product composition and returns a MultiEntry
with weights according to the coefficients if so.


* **Parameters**


    * **entry_list** (*[*[*Entry*](pymatgen.entries.md#pymatgen.entries.Entry)*]*) – list of entries from which to
    create a MultiEntry


    * **prod_comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – composition constraint for setting
    weights of MultiEntry


    * **coeff_threshold** (*float*) – threshold of stoichiometric
    coefficients to filter, if weights are lower than
    this value, the entry is not returned



#### _property_ stable_entries()
Returns the stable entries in the Pourbaix diagram.


#### _property_ unprocessed_entries()
Return unprocessed entries.


#### _property_ unstable_entries()
Returns all unstable entries in the Pourbaix diagram.


### _class_ PourbaixEntry(entry, entry_id=None, concentration=1e-06)
Bases: `MSONable`, [`Stringify`](pymatgen.util.md#pymatgen.util.string.Stringify)

An object encompassing all data relevant to a solid or ion
in a Pourbaix diagram. Each bulk solid/ion has an energy
g of the form: e = e0 + 0.0591 log10(conc) - nO mu_H2O
+ (nH - 2nO) pH + phi (-nH + 2nO + q).

Note that the energies corresponding to the input entries
should be formation energies with respect to hydrogen and
oxygen gas in order for the Pourbaix diagram formalism to
work. This may be changed to be more flexible in the future.


* **Parameters**


    * **entry** (*ComputedEntry/ComputedStructureEntry/PDEntry/IonEntry*) – An
    entry object


    * **(****)** (*concentration*) –


    * **(****)** –



#### as_dict()
Returns dict which contains Pourbaix Entry data.
Note that the pH, voltage, H2O factors are always calculated when
constructing a PourbaixEntry object.


#### _property_ composition()
Returns composition.


#### _property_ conc_term()
Returns the concentration contribution to the free energy,
and should only be present when there are ions in the entry.


#### _property_ elements()
Returns elements in the entry.


#### _property_ energy()
total energy of the Pourbaix
entry (at pH, V = 0 vs. SHE).


* **Type**

    Returns (float)



#### energy_at_conditions(pH, V)
Get free energy for a given pH and V.


* **Parameters**


    * **pH** (*float*) – pH at which to evaluate free energy


    * **V** (*float*) – voltage at which to evaluate free energy



* **Returns**

    free energy at conditions



#### _property_ energy_per_atom()
energy per atom of the Pourbaix entry.

Returns (float): energy per atom


#### _classmethod_ from_dict(d)
Invokes a PourbaixEntry from a dictionary.


#### get_element_fraction(element)
Gets the elemental fraction of a given non-OH element.


* **Parameters**

    **element** ([*Element*](pymatgen.core.md#pymatgen.core.periodic_table.Element)* or **str*) – string or element corresponding
    to element to get from composition



* **Returns**

    fraction of element / sum(all non-OH elements)



#### _property_ nH2O()
Get the number of H2O.


#### _property_ nPhi()
Get the number of electrons.


#### _property_ name()
Get the name for entry.


#### _property_ normalization_factor()
Sum of number of atoms minus the number of H and O in composition.


#### _property_ normalized_energy()
Returns:
energy normalized by number of non H or O atoms, e. g.
for Zn2O6, energy / 2 or for AgTe3(OH)3, energy / 4.


#### normalized_energy_at_conditions(pH, V)
Energy at an electrochemical condition, compatible with
numpy arrays for pH/V input.


* **Parameters**


    * **pH** (*float*) – pH at condition


    * **V** (*float*) – applied potential at condition



* **Returns**

    energy normalized by number of non-O/H atoms at condition



#### _property_ npH()
Get the number of H.


#### _property_ num_atoms()
Return number of atoms in current formula. Useful for normalization.


#### to_pretty_string()
A pretty string representation.


### _class_ PourbaixPlotter(pourbaix_diagram)
Bases: `object`

A plotter class for phase diagrams.


* **Parameters**

    **pourbaix_diagram** (*PourbaixDiagram*) – A PourbaixDiagram object.



#### domain_vertices(entry)
Returns the vertices of the Pourbaix domain.


* **Parameters**

    **entry** – Entry for which domain vertices are desired



* **Returns**

    list of vertices



#### get_pourbaix_plot(limits: tuple[float, float] | None = None, title: str = '', label_domains: bool = True, label_fontsize: int = 20, show_water_lines: bool = True, show_neutral_axes: bool = True, ax: plt.Axes = None)
Plot Pourbaix diagram.


* **Parameters**


    * **limits** – 2D list containing limits of the Pourbaix diagram
    of the form [[xlo, xhi], [ylo, yhi]]


    * **title** (*str*) – Title to display on plot


    * **label_domains** (*bool*) – whether to label Pourbaix domains


    * **label_fontsize** – font size for domain labels


    * **show_water_lines** – whether to show dashed lines indicating the region
    of water stability.


    * **lines** (*show_neutral_axes; whether to show dashed horizontal and vertical*) – at 0 V and pH 7, respectively.


    * **ax** (*Axes*) – Matplotlib Axes instance for plotting



* **Returns**

    matplotlib Axes object with Pourbaix diagram



* **Return type**

    Axes



#### plot_entry_stability(entry: Any, pH_range: tuple[float, float] | None = None, pH_resolution: int = 100, V_range: tuple[float, float] | None = None, V_resolution: int = 100, e_hull_max: float = 1, cmap: str = 'RdYlBu_r', ax: plt.Axes | None = None, \*\*kwargs: Any)
Plots the stability of an entry in the Pourbaix diagram.


* **Parameters**


    * **entry** (*Any*) – The entry to plot stability for.


    * **pH_range** (*tuple**[**float**, **float**]**, **optional*) – pH range for the plot. Defaults to [-2, 16].


    * **pH_resolution** (*int**, **optional*) – pH resolution. Defaults to 100.


    * **V_range** (*tuple**[**float**, **float**]**, **optional*) – Voltage range for the plot. Defaults to [-3, 3].


    * **V_resolution** (*int**, **optional*) – Voltage resolution. Defaults to 100.


    * **e_hull_max** (*float**, **optional*) – Maximum energy above the hull. Defaults to 1.


    * **cmap** (*str**, **optional*) – Colormap for the plot. Defaults to “RdYlBu_r”.


    * **ax** (*Axes**, **optional*) – Existing matplotlib Axes object for plotting. Defaults to None.


    * **\*\*kwargs** (*Any*) – Additional keyword arguments passed to get_pourbaix_plot.



* **Returns**

    Matplotlib Axes object with the plotted stability.



* **Return type**

    plt.Axes



#### show(\*args, \*\*kwargs)
Shows the Pourbaix plot.


* **Parameters**


    * **\*args** – args to get_pourbaix_plot


    * **\*\*kwargs** – kwargs to get_pourbaix_plot



### generate_entry_label(entry)
Generates a label for the Pourbaix plotter.


* **Parameters**

    **entry** (*PourbaixEntry** or **MultiEntry*) – entry to get a label for



### ion_or_solid_comp_object(formula)
Returns either an ion object or composition object given
a formula.


* **Parameters**

    **formula** – String formula. Eg. of ion: NaOH(aq), Na[+];
    Eg. of solid: Fe2O3(s), Fe(s), Na2O



* **Returns**

    Composition/Ion object


## pymatgen.analysis.prototypes module

This module is intended to match crystal structures against known crystallographic “prototype”
structures.

In this module, the AflowPrototypeMatcher uses the AFLOW LIBRARY OF CRYSTALLOGRAPHIC PROTOTYPES.
If using this particular class, please cite their publication appropriately:

Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
The AFLOW library of crystallographic prototypes: part 1.
Computational Materials Science, 136, S1-S828.
[https://doi.org/10.1016/j.commatsci.2017.01.017](https://doi.org/10.1016/j.commatsci.2017.01.017)


### _class_ AflowPrototypeMatcher(initial_ltol=0.2, initial_stol=0.3, initial_angle_tol=5)
Bases: `object`

This class will match structures to their crystal prototypes, and will
attempt to group species together to match structures derived from
prototypes (e.g. an A_xB_1-x_C from a binary prototype), and will
give these the names the “-like” suffix.

This class uses data from the AFLOW LIBRARY OF CRYSTALLOGRAPHIC PROTOTYPES.
If using this class, please cite their publication appropriately:

Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo, S. (2017).
The AFLOW library of crystallographic prototypes: part 1.
Computational Materials Science, 136, S1-S828.
[https://doi.org/10.1016/j.commatsci.2017.01.017](https://doi.org/10.1016/j.commatsci.2017.01.017)

Tolerances as defined in StructureMatcher. Tolerances will be
gradually decreased until only a single match is found (if possible).


* **Parameters**


    * **initial_ltol** – fractional length tolerance


    * **initial_stol** – site tolerance


    * **initial_angle_tol** – angle tolerance



#### _static_ _match_prototype(structure_matcher, structure)

#### _match_single_prototype(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))

#### get_prototypes(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Get prototype(s) structures for a given input structure. If you use this method in
your work, please cite the appropriate AFLOW publication:

Mehl, M. J., Hicks, D., Toher, C., Levy, O., Hanson, R. M., Hart, G., & Curtarolo,
S. (2017). The AFLOW library of crystallographic prototypes: part 1. Computational
Materials Science, 136, S1-S828. [https://doi.org/10.1016/j.commatsci.2017.01.017](https://doi.org/10.1016/j.commatsci.2017.01.017)


* **Parameters**

    **structure** – structure to match


Returns (list): A list of dicts with keys ‘snl’ for the matched prototype and

    ‘tags’, a dict of tags (‘mineral’, ‘strukturbericht’ and ‘aflow’) of that
    prototype. This should be a list containing just a single entry, but it is
    possible a material can match multiple prototypes.

## pymatgen.analysis.quasiharmonic module

This module implements the Quasi-harmonic Debye approximation that can
be used to compute thermal properties.

See the following papers for more info:

> [https://doi.org/10.1016/j.comphy.2003.12.001](https://doi.org/10.1016/j.comphy.2003.12.001) (2004)
> [https://doi.org/10.1103/PhysRevB.90.174107](https://doi.org/10.1103/PhysRevB.90.174107) (2014)


### _class_ QuasiharmonicDebyeApprox(energies, volumes, structure, t_min=300.0, t_step=100, t_max=300.0, eos='vinet', pressure=0.0, poisson=0.25, use_mie_gruneisen=False, anharmonic_contribution=False)
Bases: `object`

Quasiharmonic approximation.


* **Parameters**


    * **energies** (*list*) – list of DFT energies in eV


    * **volumes** (*list*) – list of volumes in Ang^3


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – pymatgen structure object


    * **t_min** (*float*) – min temperature


    * **t_step** (*float*) – temperature step


    * **t_max** (*float*) – max temperature


    * **eos** (*str*) – equation of state used for fitting the energies and the
    volumes.
    options supported by pymatgen: “quadratic”, “murnaghan”, “birch”,

    > ”birch_murnaghan”, “pourier_tarantola”, “vinet”,
    > “deltafactor”, “numerical_eos”



    * **pressure** (*float*) – in GPa, optional.


    * **poisson** (*float*) – poisson ratio.


    * **use_mie_gruneisen** (*bool*) – whether or not to use the mie-gruneisen
    formulation to compute the gruneisen parameter.
    The default is the slater-gamma formulation.


    * **anharmonic_contribution** (*bool*) – whether or not to consider the anharmonic
    contribution to the Debye temperature. Cannot be used with
    use_mie_gruneisen. Defaults to False.



#### _static_ debye_integral(y)
Debye integral. Eq(5) in  doi.org/10.1016/j.comphy.2003.12.001.


* **Parameters**

    **y** (*float*) – debye temperature/T, upper limit



* **Returns**

    unitless



* **Return type**

    float



#### debye_temperature(volume)
Calculates the debye temperature.
Eq(6) in doi.org/10.1016/j.comphy.2003.12.001. Thanks to Joey.

Eq(6) above is equivalent to Eq(3) in doi.org/10.1103/PhysRevB.37.790
which does not consider anharmonic effects. Eq(20) in the same paper
and Eq(18) in doi.org/10.1016/j.commatsci.2009.12.006 both consider
anharmonic contributions to the Debye temperature through the Gruneisen
parameter at 0K (Gruneisen constant).

The anharmonic contribution is toggled by setting the anharmonic_contribution
to True or False in the QuasiharmonicDebyeApprox constructor.


* **Parameters**

    **volume** (*float*) – in Ang^3



* **Returns**

    debye temperature in K



* **Return type**

    float



#### get_summary_dict()
Returns a dict with a summary of the computed properties.


#### gruneisen_parameter(temperature, volume)
Slater-gamma formulation(the default):

    gruneisen parameter = - d log(theta)/ d log(V) = - (1/6 + 0.5 d log(B)/ d log(V))

        = - (1/6 + 0.5 V/B dB/dV), where dB/dV = d^2E/dV^2 + V \* d^3E/dV^3.

Mie-gruneisen formulation:

    Eq(31) in doi.org/10.1016/j.comphy.2003.12.001
    Eq(7) in Blanco et. al. Joumal of Molecular Structure (Theochem)

    > 368 (1996) 245-255

    Also se J.P. Poirier, Introduction to the Physics of the Earth’s

        Interior, 2nd ed. (Cambridge University Press, Cambridge,
        2000) Eq(3.53)


* **Parameters**


    * **temperature** (*float*) – temperature in K


    * **volume** (*float*) – in Ang^3



* **Returns**

    unitless



* **Return type**

    float



#### optimize_gibbs_free_energy()
Evaluate the Gibbs free energy as a function of V, T and P i.e
G(V, T, P), minimize G(V, T, P) wrt V for each T and store the
optimum values.

Note: The data points for which the equation of state fitting fails

    are skipped.


#### optimizer(temperature)
Evaluate G(V, T, P) at the given temperature(and pressure) and
minimize it wrt V.


1. Compute the  vibrational Helmholtz free energy, A_vib.


2. Compute the Gibbs free energy as a function of volume, temperature

    and pressure, G(V,T,P).


3. Perform an equation of state fit to get the functional form of

    Gibbs free energy:G(V, T, P).


4. Finally G(V, P, T) is minimized with respect to V.


* **Parameters**

    **temperature** (*float*) – temperature in K



* **Returns**

    G_opt(V_opt, T, P) in eV and V_opt in Ang^3.



* **Return type**

    float, float



#### thermal_conductivity(temperature, volume)
Eq(17) in 10.1103/PhysRevB.90.174107.


* **Parameters**


    * **temperature** (*float*) – temperature in K


    * **volume** (*float*) – in Ang^3



* **Returns**

    thermal conductivity in W/K/m



* **Return type**

    float



#### vibrational_free_energy(temperature, volume)
Vibrational Helmholtz free energy, A_vib(V, T).
Eq(4) in doi.org/10.1016/j.comphy.2003.12.001.


* **Parameters**


    * **temperature** (*float*) – temperature in K


    * **volume** (*float*) –



* **Returns**

    vibrational free energy in eV



* **Return type**

    float



#### vibrational_internal_energy(temperature, volume)
Vibrational internal energy, U_vib(V, T).
Eq(4) in doi.org/10.1016/j.comphy.2003.12.001.


* **Parameters**


    * **temperature** (*float*) – temperature in K


    * **volume** (*float*) – in Ang^3



* **Returns**

    vibrational internal energy in eV



* **Return type**

    float


## pymatgen.analysis.quasirrho module

A module to calculate free energies using the Quasi-Rigid Rotor Harmonic
Oscillator approximation. Modified from a script by Steven Wheeler.
See: Grimme, S. Chem. Eur. J. 2012, 18, 9955


### _class_ QuasiRRHO(mol: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), frequencies: list[float], energy: float, mult: int, sigma_r: float = 1, temp: float = 298.15, press: float = 101317, v0: float = 100)
Bases: `object`

Class to calculate thermochemistry using Grimme’s Quasi-RRHO approximation.
All outputs are in atomic units, e.g. energy outputs are in Hartrees.
Citation: Grimme, S. Chemistry - A European Journal 18, 9955-9964 (2012).


#### temp()
Temperature [K]


* **Type**

    float



#### press()
Pressure [Pa]


* **Type**

    float



#### v0()
Cutoff frequency for Quasi-RRHO method [1/cm]


* **Type**

    float



#### entropy_quasiRRHO()
Quasi-RRHO entropy [Ha/K]


* **Type**

    float



#### entropy_ho()
Total entropy calculated with a harmonic
oscillator approximation for the vibrational entropy [Ha/K]


* **Type**

    float



#### h_corrected()
Thermal correction to the enthalpy [Ha]


* **Type**

    float



#### free_energy_quasiRRHO()
Quasi-RRHO free energy [Ha]


* **Type**

    float



#### free_energy_ho()
Free energy calculated without the Quasi-RRHO
method, i.e. with a harmonic oscillator approximation for the
vibrational entropy [Ha]


* **Type**

    float



* **Parameters**


    * **mol** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Pymatgen molecule


    * **frequencies** (*list*) – List of frequencies (float) [cm^-1]


    * **energy** (*float*) – Electronic energy [Ha]


    * **mult** (*int*) – Spin multiplicity


    * **sigma_r** (*int*) – Rotational symmetry number. Defaults to 1.


    * **temp** (*float*) – Temperature [K]. Defaults to 298.15.


    * **press** (*float*) – Pressure [Pa]. Defaults to 101_317.


    * **v0** (*float*) – Cutoff frequency for Quasi-RRHO method [cm^-1]. Defaults to 100.



#### _get_quasirrho_thermo(mol: [Molecule](pymatgen.core.md#pymatgen.core.structure.Molecule), mult: int, sigma_r: int, frequencies: list[float], elec_energy: float)
Calculate Quasi-RRHO thermochemistry


* **Parameters**


    * **mol** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Pymatgen molecule


    * **mult** (*int*) – Spin multiplicity


    * **sigma_r** (*int*) – Rotational symmetry number


    * **frequencies** (*list*) – List of frequencies [cm^-1]


    * **elec_energy** (*float*) – Electronic energy [Ha]



#### _classmethod_ from_gaussian_output(output: [GaussianOutput](pymatgen.io.md#pymatgen.io.gaussian.GaussianOutput), \*\*kwargs)

* **Parameters**

    **output** ([*GaussianOutput*](pymatgen.io.md#pymatgen.io.gaussian.GaussianOutput)) – Pymatgen GaussianOutput object



* **Returns**

    QuasiRRHO class instantiated from a Gaussian Output



* **Return type**

    QuasiRRHO



#### _classmethod_ from_qc_output(output: [QCOutput](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput), \*\*kwargs)

* **Parameters**

    **output** ([*QCOutput*](pymatgen.io.qchem.md#pymatgen.io.qchem.outputs.QCOutput)) – Pymatgen QCOutput object



* **Returns**

    QuasiRRHO class instantiated from a QChem Output



* **Return type**

    QuasiRRHO



### get_avg_mom_inertia(mol)
Calculate the average moment of inertia of a molecule


* **Parameters**

    **mol** ([*Molecule*](pymatgen.core.md#pymatgen.core.structure.Molecule)) – Pymatgen Molecule



* **Returns**

    average moment of inertia, eigenvalues of the inertia tensor



* **Return type**

    int, list


## pymatgen.analysis.reaction_calculator module

This module provides classes that define a chemical reaction.


### _class_ BalancedReaction(reactants_coeffs, products_coeffs)
Bases: `MSONable`

An object representing a complete chemical reaction.

Reactants and products to be specified as dict of {Composition: coeff}.


* **Parameters**


    * **reactants_coeffs** (*dict**[*[*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)*, **float**]*) – Reactants as dict of {Composition: amt}.


    * **products_coeffs** (*dict**[*[*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)*, **float**]*) – Products as dict of {Composition: amt}.



#### TOLERANCE(_ = 1e-0_ )

#### _classmethod_ _str_from_comp(coeffs, compositions, reduce=False)

#### _classmethod_ _str_from_formulas(coeffs, formulas)

#### _property_ all_comp()
List of all compositions in the reaction.


#### as_dict()

* **Returns**

    A dictionary representation of BalancedReaction.



#### as_entry(energies)
Returns a ComputedEntry representation of the reaction.


#### calculate_energy(energies)
Calculates the energy of the reaction.


* **Parameters**

    **(****{Composition** (*energies*) – float}): Energy for each composition.
    E.g ., {comp1: energy1, comp2: energy2}.



* **Returns**

    reaction energy as a float.



#### _property_ coeffs()
Final coefficients of the calculated reaction.


#### _property_ elements()
List of elements in the reaction.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – from as_dict().



* **Returns**

    A BalancedReaction object.



#### _static_ from_str(rxn_str)
Generates a balanced reaction from a string. The reaction must
already be balanced.


* **Parameters**

    **rxn_string** (*str*) – The reaction string. For example, “4 Li + O2 -> 2Li2O”



* **Returns**

    BalancedReaction



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_coeff(comp)
Returns coefficient for a particular composition.


#### get_el_amount(element)
Returns the amount of the element in the reaction.


* **Parameters**

    **element** (*Element/Species*) – Element in the reaction



* **Returns**

    Amount of that element in the reaction.



#### normalize_to(comp, factor=1)
Normalizes the reaction to one of the compositions.
By default, normalizes such that the composition given has a
coefficient of 1. Another factor can be specified.


* **Parameters**


    * **comp** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – Composition to normalize to


    * **factor** (*float*) – Factor to normalize to. Defaults to 1.



#### normalize_to_element(element, factor=1)
Normalizes the reaction to one of the elements.
By default, normalizes such that the amount of the element is 1.
Another factor can be specified.


* **Parameters**


    * **element** (*Element/Species*) – Element to normalize to.


    * **factor** (*float*) – Factor to normalize to. Defaults to 1.



#### _property_ normalized_repr()
A normalized representation of the reaction. All factors are converted
to lowest common factors.


#### normalized_repr_and_factor()
Normalized representation for a reaction
For example, `4 Li + 2 O -> 2Li2O` becomes `2 Li + O -> Li2O`.


#### _property_ products()
List of products.


#### _property_ reactants()
List of reactants.


### _class_ ComputedReaction(reactant_entries, product_entries)
Bases: `Reaction`

Convenience class to generate a reaction from ComputedEntry objects, with
some additional attributes, such as a reaction energy based on computed
energies.


* **Parameters**


    * **reactant_entries** (*[*[*ComputedEntry*](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)*]*) – List of reactant_entries.


    * **product_entries** (*[*[*ComputedEntry*](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedEntry)*]*) – List of product_entries.



#### _property_ all_entries()
Equivalent of all_comp but returns entries, in the same order as the
coefficients.


#### as_dict()

* **Returns**

    A dictionary representation of ComputedReaction.



#### _property_ calculated_reaction_energy()
Returns (float):
The calculated reaction energy.


#### _property_ calculated_reaction_energy_uncertainty()
Calculates the uncertainty in the reaction energy based on the uncertainty in the
energies of the products and reactants.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – from as_dict().



* **Returns**

    A ComputedReaction object.



### _class_ Reaction(reactants, products)
Bases: `BalancedReaction`

A more flexible class representing a Reaction. The reaction amounts will
be automatically balanced. Reactants and products can swap sides so that
all coefficients are positive, however this class will find the solution
with the minimum number of swaps and coefficients of 0. Normalizes so that
the *FIRST* product (or products, if underdetermined) has a coefficient of one.

Reactants and products to be specified as list of
pymatgen.core.structure.Composition. e.g., [comp1, comp2].


* **Parameters**


    * **reactants** (*[*[*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)*]*) – List of reactants.


    * **products** (*[*[*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)*]*) – List of products.



#### _balance_coeffs(comp_matrix, max_num_constraints)

#### as_dict()

* **Returns**

    A dictionary representation of Reaction.



#### copy()
Returns a copy of the Reaction object.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – from as_dict().



* **Returns**

    A Reaction object.



### _exception_ ReactionError(msg)
Bases: `Exception`

Exception class for Reactions. Allows more information in exception
messages to cover situations not covered by standard exception classes.

Create a ReactionError.


* **Parameters**

    **msg** (*str*) – More information about the ReactionError.


## pymatgen.analysis.structure_analyzer module

This module provides classes to perform topological analyses of structures.


### _class_ OxideType(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), relative_cutoff=1.1)
Bases: `object`

Separate class for determining oxide type.


* **Parameters**


    * **structure** – Input structure.


    * **relative_cutoff** – Relative_cutoff \* act. cutoff stipulates the max.
    distance two O atoms must be from each other. Default value is
    1.1. At most 1.1 is recommended, nothing larger, otherwise the
    script cannot distinguish between superoxides and peroxides.



#### parse_oxide()
Determines if an oxide is a peroxide/superoxide/ozonide/normal oxide.


* **Returns**

    Type of oxide
    ozonide/peroxide/superoxide/hydroxide/None.
    nbonds (int): Number of peroxide/superoxide/hydroxide bonds in structure.



* **Return type**

    oxide_type (str)



### _class_ RelaxationAnalyzer(initial_structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), final_structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Bases: `object`

This class analyzes the relaxation in a calculation.

Please note that the input and final structures should have the same
ordering of sites. This is typically the case for most computational codes.


* **Parameters**


    * **initial_structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Initial input structure to
    calculation.


    * **final_structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Final output structure from
    calculation.



* **Raises**

    **ValueError** – If initial and final structures have different formulas.



#### get_percentage_bond_dist_changes(max_radius: float = 3.0)
Returns the percentage bond distance changes for each site up to a
maximum radius for nearest neighbors.


* **Parameters**

    **max_radius** (*float*) – Maximum radius to search for nearest
    neighbors. This radius is applied to the initial structure,
    not the final structure.



* **Returns**

    Bond distance changes in the form {index1: {index2: 0.011, …}}.

        For economy of representation, the index1 is always less than index2, i.e., since bonding
        between site1 and site_n is the same as bonding between site_n and site1, there is no
        reason to duplicate the information or computation.




* **Return type**

    dict[int, dict[int, float]]



#### get_percentage_lattice_parameter_changes()
Returns the percentage lattice parameter changes.


* **Returns**

    Percent changes in lattice parameter, e.g.,

        {‘a’: 0.012, ‘b’: 0.021, ‘c’: -0.031} implies a change of 1.2%,
        2.1% and -3.1% in the a, b and c lattice parameters respectively.




* **Return type**

    dict[str, float]



#### get_percentage_volume_change()
Returns the percentage volume change.


* **Returns**

    Volume change in percent. 0.055 means a 5.5% increase.



* **Return type**

    float



### _class_ VoronoiAnalyzer(cutoff=5.0, qhull_options='Qbb Qc Qz')
Bases: `object`

Performs a statistical analysis of Voronoi polyhedra around each site.
Each Voronoi polyhedron is described using Schaefli notation.
That is a set of indices {c_i} where c_i is the number of faces with i
number of vertices. E.g. for a bcc crystal, there is only one polyhedron
notation of which is [0,6,0,8,0,0,…].
In perfect crystals, these also corresponds to the Wigner-Seitz cells.
For distorted-crystals, liquids or amorphous structures, rather than one-type,
there is a statistical distribution of polyhedra.
See ref: Microstructure and its relaxation in Fe-B amorphous system
simulated by molecular dynamics,

> Stepanyuk et al., J. Non-cryst. Solids (1993), 159, 80-87.


* **Parameters**


    * **cutoff** (*float*) – cutoff distance to search for neighbors of a given atom
    (default = 5.0)


    * **qhull_options** (*str*) – options to pass to qhull (optional).



#### analyze(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), n=0)
Performs Voronoi analysis and returns the polyhedra around atom n
in Schlaefli notation.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure to analyze


    * **n** (*int*) – index of the center atom in structure



* **Returns**

    <c3,c4,c6,c6,c7,c8,c9,c10>

        where c_i denotes number of facets with i vertices.




* **Return type**

    voronoi index of n



#### analyze_structures(structures, step_freq=10, most_frequent_polyhedra=15)
Perform Voronoi analysis on a list of Structures.
Note that this might take a significant amount of time depending on the
size and number of structures.


* **Parameters**


    * **structures** (*list*) – list of Structures


    * **(****float** (*cutoff*) – cutoff distance around an atom to search for
    neighbors


    * **step_freq** (*int*) – perform analysis every step_freq steps


    * **qhull_options** (*str*) – options to pass to qhull


    * **most_frequent_polyhedra** (*int*) – this many unique polyhedra with
    highest frequencies is stored.



* **Returns**

    A list of tuples in the form (voronoi_index,frequency)



#### _static_ plot_vor_analysis(voronoi_ensemble: list[tuple[str, float]])
Plot the Voronoi analysis.


* **Parameters**

    **voronoi_ensemble** (*list**[**tuple**[**str**, **float**]**]*) – List of tuples containing labels and
    values for Voronoi analysis.



* **Returns**

    Matplotlib Axes object with the plotted Voronoi analysis.



* **Return type**

    plt.Axes



### _class_ VoronoiConnectivity(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), cutoff=10)
Bases: `object`

Computes the solid angles swept out by the shared face of the voronoi
polyhedron between two sites.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure


    * **cutoff** (*float*) –



#### _property_ connectivity_array()
Provides connectivity array.


* **Returns**

    An array of shape [atom_i, atom_j, image_j]. atom_i is
    the index of the atom in the input structure. Since the second
    atom can be outside of the unit cell, it must be described
    by both an atom index and an image index. Array data is the
    solid angle of polygon between atom_i and image_j of atom_j



* **Return type**

    connectivity



#### get_connections()
Returns a list of site pairs that are Voronoi Neighbors, along
with their real-space distances.


#### get_sitej(site_index, image_index)
Assuming there is some value in the connectivity array at indices
(1, 3, 12). sitei can be obtained directly from the input structure
(structure[1]). sitej can be obtained by passing 3, 12 to this function.


* **Parameters**


    * **site_index** (*int*) – index of the site (3 in the example)


    * **image_index** (*int*) – index of the image (12 in the example)



#### _property_ max_connectivity()
Returns the 2d array [site_i, site_j] that represents the maximum connectivity of
site i to any periodic image of site j.


### average_coordination_number(structures, freq=10)
Calculates the ensemble averaged Voronoi coordination numbers
of a list of Structures using VoronoiNN.
Typically used for analyzing the output of a Molecular Dynamics run.


* **Parameters**


    * **structures** (*list*) – list of Structures.


    * **freq** (*int*) – sampling frequency of coordination number [every freq steps].



* **Returns**

    Dictionary of elements as keys and average coordination numbers as values.



### contains_peroxide(structure, relative_cutoff=1.1)
Determines if a structure contains peroxide anions.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure.


    * **relative_cutoff** – The peroxide bond distance is 1.49 Angstrom.
    Relative_cutoff \* 1.49 stipulates the maximum distance two O
    atoms must be to each other to be considered a peroxide.



* **Returns**

    Boolean indicating if structure contains a peroxide anion.



### get_max_bond_lengths(structure, el_radius_updates=None)
Provides max bond length estimates for a structure based on the JMol
table and algorithms.


* **Parameters**


    * **structure** – (structure)


    * **el_radius_updates** – (dict) symbol->float to update atom_ic radii



* **Returns**

    The two elements are ordered by Z.



* **Return type**

    dict[(Element1, Element2)], float]



### oxide_type(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), relative_cutoff: float = 1.1, return_nbonds: bool = False)
Determines if an oxide is a peroxide/superoxide/ozonide/normal oxide.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure.


    * **relative_cutoff** (*float*) – Relative_cutoff \* act. cutoff stipulates the
    max distance two O atoms must be from each other.


    * **return_nbonds** (*bool*) – Should number of bonds be requested?



### solid_angle(center, coords)
Helper method to calculate the solid angle of a set of coords from the center.


* **Parameters**


    * **center** (*3x1 array*) – Center to measure solid angle from.


    * **coords** (*Nx3 array*) – List of coords to determine solid angle.



* **Returns**

    The solid angle.



* **Return type**

    float



### sulfide_type(structure)
Determines if a structure is a sulfide/polysulfide/sulfate.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Input structure.



* **Returns**

    (str) sulfide/polysulfide or None if structure is a sulfate.


## pymatgen.analysis.structure_matcher module

This module provides classes to perform fitting of structures.


### _class_ AbstractComparator()
Bases: `MSONable`

Abstract Comparator class. A Comparator defines how sites are compared in
a structure.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### _abstract_ are_equal(sp1, sp2)
Defines how the species of two sites are considered equal. For
example, one can consider sites to have the same species only when
the species are exactly the same, i.e., Fe2+ matches Fe2+ but not
Fe3+. Or one can define that only the element matters,
and all oxidation state information are ignored.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    Boolean indicating whether species are considered equal.



#### as_dict()
MSONable dict


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    Comparator.



#### _abstract_ get_hash(composition)
Defines a hash to group structures. This allows structures to be
grouped efficiently for comparison. The hash must be invariant under
supercell creation. (e.g. composition is not a good hash, but
fractional_composition might be). Reduced formula is not a good formula,
due to weird behavior with fractional occupancy.

Composition is used here instead of structure because for anonymous
matches it is much quicker to apply a substitution to a composition
object than a structure object.


* **Parameters**

    **composition** ([*Composition*](pymatgen.core.md#pymatgen.core.composition.Composition)) – composition of the structure



* **Returns**

    A hashable object. Examples can be string formulas, integers etc.



### _class_ ElementComparator()
Bases: `AbstractComparator`

A Comparator that matches elements. i.e. oxidation states are
ignored.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### are_equal(sp1, sp2)
True if element:amounts are exactly the same, i.e.,
oxidation state is not considered.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    Boolean indicating whether species are the same based on element
    and amounts.



#### get_hash(composition)
Returns: Fractional element composition.


### _class_ FrameworkComparator()
Bases: `AbstractComparator`

A Comparator that matches sites, regardless of species.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### are_equal(sp1, sp2)
True if there are atoms on both sites.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    True always



#### get_hash(composition)
No hash possible.


### _class_ OccupancyComparator()
Bases: `AbstractComparator`

A Comparator that matches occupancies on sites,
irrespective of the species of those sites.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### are_equal(sp1, sp2)

* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    True if sets of occupancies (amt) are equal on both sites.



* **Return type**

    bool



#### get_hash(composition)

* **Parameters**

    **composition** – Composition.



* **Returns**


    1. Difficult to define sensible hash




### _class_ OrderDisorderElementComparator()
Bases: `AbstractComparator`

A Comparator that matches sites, given some overlap in the element
composition.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### are_equal(sp1, sp2)
True if there is some overlap in composition between the species.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    True always



#### get_hash(composition)
Returns: Fractional composition.


### _class_ SpeciesComparator()
Bases: `AbstractComparator`

A Comparator that matches species exactly. The default used in StructureMatcher.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### are_equal(sp1, sp2)
True if species are exactly the same, i.e., Fe2+ == Fe2+ but not Fe3+.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    Boolean indicating whether species are equal.



#### get_hash(composition: [Composition](pymatgen.core.md#pymatgen.core.composition.Composition))
Returns: Fractional composition.


### _class_ SpinComparator()
Bases: `AbstractComparator`

A Comparator that matches magnetic structures to their inverse spins.
This comparator is primarily used to filter magnetically ordered
structures with opposite spins, which are equivalent.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### are_equal(sp1, sp2)
True if species are exactly the same, i.e., Fe2+ == Fe2+ but not
Fe3+. and the spins are reversed. i.e., spin up maps to spin down,
and vice versa.


* **Parameters**


    * **sp1** – First species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.


    * **sp2** – Second species. A dict of {specie/element: amt} as per the
    definition in Site and PeriodicSite.



* **Returns**

    Boolean indicating whether species are equal.



#### get_hash(composition)
Returns: Fractional composition.


### _class_ StructureMatcher(ltol: float = 0.2, stol: float = 0.3, angle_tol: float = 5, primitive_cell: bool = True, scale: bool = True, attempt_supercell: bool = False, allow_subset: bool = False, comparator: AbstractComparator | None = None, supercell_size: Literal['num_sites', 'num_atoms', 'volume'] = 'num_sites', ignored_species: Sequence[SpeciesLike] = ())
Bases: `MSONable`

Class to match structures by similarity.

Algorithm:


1. Given two structures: s1 and s2


2. Optional: Reduce to primitive cells.


3. If the number of sites do not match, return False


4. Reduce to s1 and s2 to Niggli Cells


5. Optional: Scale s1 and s2 to same volume.


6. Optional: Remove oxidation states associated with sites


7. Find all possible lattice vectors for s2 within shell of ltol.


8. For s1, translate an atom in the smallest set to the origin


9. For s2: find all valid lattices from permutations of the list
of lattice vectors (invalid if: det(Lattice Matrix) < half
volume of original s2 lattice)


10. For each valid lattice:


    1. If the lattice angles of are within tolerance of s1,
basis change s2 into new lattice.


    2. For each atom in the smallest set of s2:

> i. Translate to origin and compare fractional sites in
> structure within a fractional tolerance.
> ii. If true:

> > ia. Convert both lattices to Cartesian and place
> > both structures on an average lattice
> > ib. Compute and return the average and max rms
> > displacement between the two structures normalized
> > by the average free length per atom

> > if fit function called:

> >     if normalized max rms displacement is less than
> >     stol. Return True

> > if get_rms_dist function called:

> >     if normalized average rms displacement is less
> >     than the stored rms displacement, store and
> >     continue. (This function will search all possible
> >     lattices for the smallest average rms displacement
> >     between the two structures)


* **Parameters**


    * **ltol** (*float*) – Fractional length tolerance. Default is 0.2.


    * **stol** (*float*) – Site tolerance. Defined as the fraction of the
    average free length per atom := ( V / Nsites ) \*\* (1/3)
    Default is 0.3.


    * **angle_tol** (*float*) – Angle tolerance in degrees. Default is 5 degrees.


    * **primitive_cell** (*bool*) – If true: input structures will be reduced to
    primitive cells prior to matching. Default to True.


    * **scale** (*bool*) – Input structures are scaled to equivalent volume if
    true; For exact matching, set to False.


    * **attempt_supercell** (*bool*) – If set to True and number of sites in
    cells differ after a primitive cell reduction (divisible by an
    integer) attempts to generate a supercell transformation of the
    smaller cell which is equivalent to the larger structure.


    * **allow_subset** (*bool*) – Allow one structure to match to the subset of
    another structure. Eg. Matching of an ordered structure onto a
    disordered one, or matching a delithiated to a lithiated
    structure. This option cannot be combined with
    attempt_supercell, or with structure grouping.


    * **comparator** (*Comparator*) – A comparator object implementing an equals
    method that declares equivalency of sites. Default is
    SpeciesComparator, which implies rigid species
    mapping, i.e., Fe2+ only matches Fe2+ and not Fe3+.

    Other comparators are provided, e.g., ElementComparator which
    matches only the elements and not the species.

    The reason why a comparator object is used instead of
    supplying a comparison function is that it is not possible to
    pickle a function, which makes it otherwise difficult to use
    StructureMatcher with Python’s multiprocessing.



    * **supercell_size** (*str** or **list*) – Method to use for determining the
    size of a supercell (if applicable). Possible values are
    ‘num_sites’, ‘num_atoms’, ‘volume’, or an element or list of elements
    present in both structures.


    * **ignored_species** (*list*) – A list of ions to be ignored in matching.
    Useful for matching structures that have similar frameworks
    except for certain ions, e.g., Li-ion intercalation frameworks.
    This is more useful than allow_subset because it allows better
    control over what species are ignored in the matching.



#### _anonymous_match(struct1: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), struct2: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), fu: int, s1_supercell=True, use_rms=False, break_on_match=False, single_match=False)
Tries all permutations of matching struct1 to struct2.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – First structure


    * **struct2** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Second structure


    * **fu** (*int*) – Factor of unit cell of struct1 to match to struct2


    * **s1_supercell** (*bool*) – whether to create the supercell of struct1 (vs struct2)


    * **use_rms** (*bool*) – Whether to minimize the rms of the matching


    * **break_on_match** (*bool*) – Whether to break search on first match


    * **single_match** (*bool*) – Whether to return only the best match



* **Returns**

    List of (mapping, match)



#### _classmethod_ _cart_dists(s1, s2, avg_lattice, mask, normalization, lll_frac_tol=None)
Finds a matching in Cartesian space. Finds an additional
fractional translation vector to minimize RMS distance.


* **Parameters**


    * **s1** – numpy array of fractional coordinates.


    * **s2** – numpy array of fractional coordinates. len(s1) >= len(s2)


    * **avg_lattice** – Lattice on which to calculate distances


    * **mask** – numpy array of booleans. mask[i, j] = True indicates
    that s2[i] cannot be matched to s1[j]


    * **normalization** (*float*) – inverse normalization length


    * **lll_frac_tol** (*float*) – tolerance for Lenstra-Lenstra-Lovász lattice basis reduction algorithm



* **Returns**

    Distances from s2 to s1, normalized by (V/atom) ^ 1/3
    Fractional translation vector to apply to s2.
    Mapping from s1 to s2, i.e. with numpy slicing, s1[mapping] => s2



#### _classmethod_ _cmp_fstruct(s1, s2, frac_tol, mask)
Returns true if a matching exists between s2 and s2
under frac_tol. s2 should be a subset of s1.


#### _get_lattices(target_lattice, s, supercell_size=1)
Yields lattices for s with lengths and angles close to the lattice of target_s. If
supercell_size is specified, the returned lattice will have that number of primitive
cells in it.


* **Parameters**


    * **target_lattice** ([*Lattice*](pymatgen.core.md#pymatgen.core.lattice.Lattice)) – target lattice.


    * **s** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – input structure.


    * **supercell_size** (*int*) – Number of primitive cells in returned lattice



#### _get_mask(struct1, struct2, fu, s1_supercell)
Returns mask for matching struct2 to struct1. If struct1 has sites
a b c, and fu = 2, assumes supercells of struct2 will be ordered
aabbcc (rather than abcabc).


* **Returns**

    mask, struct1 translation indices, struct2 translation index



#### _classmethod_ _get_reduced_structure(struct: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), primitive_cell: bool = True, niggli: bool = True)
Helper method to find a reduced structure.


#### _get_supercell_size(s1, s2)
Returns the supercell size, and whether the supercell should be applied to s1.
If fu == 1, s1_supercell is returned as true, to avoid ambiguity.


#### _get_supercells(struct1, struct2, fu, s1_supercell)
Computes all supercells of one structure close to the lattice of the
other
if s1_supercell is True, it makes the supercells of struct1, otherwise
it makes them of s2.

yields: s1, s2, supercell_matrix, average_lattice, supercell_matrix


#### _match(struct1, struct2, fu, s1_supercell=True, use_rms=False, break_on_match=False)
Matches one struct onto the other.


#### _preprocess(struct1, struct2, niggli=True, skip_structure_reduction: bool = False)
Rescales, finds the reduced structures (primitive and niggli),
and finds fu, the supercell size to make struct1 comparable to
s2.
If skip_structure_reduction is True, skip to get reduced structures (by primitive transformation and
niggli reduction). This option is useful for fitting a set of structures several times.


#### _process_species(structures)

#### _strict_match(struct1: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), struct2: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), fu: int, s1_supercell: bool = True, use_rms: bool = False, break_on_match: bool = False)
Matches struct2 onto struct1 (which should contain all sites in
struct2).


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure to match onto


    * **struct2** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – structure to match


    * **fu** (*int*) – size of supercell to create


    * **s1_supercell** (*bool*) – whether to create the supercell of struct1 (vs struct2)


    * **use_rms** (*bool*) – whether to minimize the rms of the matching


    * **break_on_match** (*bool*) – whether to stop search at first match



* **Returns**

    (rms, max_dist, mask, cost, mapping)

        if a match is found, else None




* **Return type**

    tuple[float, float, np.ndarray, float, Mapping]



#### as_dict()
MSONable dict


#### fit(struct1: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), struct2: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), symmetric: bool = False, skip_structure_reduction: bool = False)
Fit two structures.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 2nd structure


    * **symmetric** (*bool*) – Defaults to False
    If True, check the equality both ways.
    This only impacts a small percentage of structures


    * **skip_structure_reduction** (*bool*) – Defaults to False
    If True, skip to get a primitive structure and perform Niggli reduction for struct1 and struct2



* **Returns**

    True if the structures are equivalent



* **Return type**

    bool



#### fit_anonymous(struct1: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), struct2: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), niggli: bool = True, skip_structure_reduction: bool = False)
Performs an anonymous fitting, which allows distinct species in one structure to map
to another. E.g., to compare if the Li2O and Na2O structures are similar.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 2nd structure


    * **niggli** (*bool*) – If true, perform Niggli reduction for struct1 and struct2


    * **skip_structure_reduction** (*bool*) – Defaults to False
    If True, skip to get a primitive structure and perform Niggli reduction for struct1 and struct2



* **Returns**

    Whether a species mapping can map struct1 to struct2



* **Return type**

    bool



#### _classmethod_ from_dict(d)

* **Parameters**

    **d** – Dict representation



* **Returns**

    StructureMatcher



#### get_all_anonymous_mappings(struct1, struct2, niggli=True, include_dist=False)
Performs an anonymous fitting, which allows distinct species in one
structure to map to another. Returns a dictionary of species
substitutions that are within tolerance.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 2nd structure


    * **niggli** (*bool*) – Find niggli cell in preprocessing


    * **include_dist** (*bool*) – Return the maximin distance with each mapping



* **Returns**

    list of species mappings that map struct1 to struct2.



#### get_best_electronegativity_anonymous_mapping(struct1: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), struct2: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure))
Performs an anonymous fitting, which allows distinct species in one
structure to map to another. E.g., to compare if the Li2O and Na2O
structures are similar. If multiple substitutions are within tolerance
this will return the one which minimizes the difference in
electronegativity between the matches species.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 2nd structure



* **Returns**

    Mapping of struct1 species to struct2 species



* **Return type**

    min_mapping (dict)



#### get_mapping(superset, subset)
Calculate the mapping from superset to subset.


* **Parameters**


    * **superset** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure containing at least the sites in
    subset (within the structure matching tolerance)


    * **subset** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure containing some of the sites in
    superset (within the structure matching tolerance)



* **Returns**

    numpy array such that superset.sites[mapping] is within matching
    tolerance of subset.sites or None if no such mapping is possible



#### get_rms_anonymous(struct1, struct2)
Performs an anonymous fitting, which allows distinct species in one
structure to map to another. E.g., to compare if the Li2O and Na2O
structures are similar.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 2nd structure



* **Returns**

    (min_rms, min_mapping)
    min_rms is the minimum rms distance, and min_mapping is the
    corresponding minimal species mapping that would map
    struct1 to struct2. (None, None) is returned if the minimax_rms
    exceeds the threshold.



#### get_rms_dist(struct1, struct2)
Calculate RMS displacement between two structures.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 1st structure


    * **struct2** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – 2nd structure



* **Returns**

    rms displacement normalized by (Vol / nsites) \*\* (1/3)
    and maximum distance between paired sites. If no matching
    lattice is found None is returned.



#### get_s2_like_s1(struct1, struct2, include_ignored_species=True)
Performs transformations on struct2 to put it in a basis similar to
struct1 (without changing any of the inter-site distances).


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Reference structure


    * **struct2** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure to transform.


    * **include_ignored_species** (*bool*) – Defaults to True,
    the ignored_species is also transformed to the struct1
    lattice orientation, though obviously there is no direct
    matching to existing sites.



* **Returns**

    A structure object similar to struct1, obtained by making a
    supercell, sorting, and translating struct2.



#### get_supercell_matrix(supercell, struct)
Returns the matrix for transforming struct to supercell. This
can be used for very distorted ‘supercells’ where the primitive cell
is impossible to find.


#### get_transformation(struct1, struct2)
Returns the supercell transformation, fractional translation vector,
and a mapping to transform struct2 to be similar to struct1.


* **Parameters**


    * **struct1** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Reference structure


    * **struct2** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure to transform.



* **Returns**

    supercell matrix
    vector (numpy.ndarray(3)): fractional translation vector
    mapping (list(int or None)):

    > The first len(struct1) items of the mapping vector are the
    > indices of struct1’s corresponding sites in struct2 (or None
    > if there is no corresponding site), and the other items are
    > the remaining site indices of struct2.




* **Return type**

    supercell (numpy.ndarray(3, 3))



#### group_structures(s_list, anonymous=False)
Given a list of structures, use fit to group
them by structural equality.


* **Parameters**


    * **s_list** (*[*[*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)*]*) – List of structures to be grouped


    * **anonymous** (*bool*) – Whether to use anonymous mode.



* **Returns**

    A list of lists of matched structures
    Assumption: if s1 == s2 but s1 != s3, than s2 and s3 will be put
    in different groups without comparison.


## pymatgen.analysis.surface_analysis module

This module defines tools to analyze surface and adsorption related
quantities as well as related plots. If you use this module, please
consider citing the following works:

```default
R. Tran, Z. Xu, B. Radhakrishnan, D. Winston, W. Sun, K. A. Persson,
S. P. Ong, "Surface Energies of Elemental Crystals", Scientific
Data, 2016, 3:160080, doi: 10.1038/sdata.2016.80.

and

Kang, S., Mo, Y., Ong, S. P., & Ceder, G. (2014). Nanoscale
stabilization of sodium oxides: Implications for Na-O2 batteries.
Nano Letters, 14(2), 1016-1020. https://doi.org/10.1021/nl404557w

and

Montoya, J. H., & Persson, K. A. (2017). A high-throughput framework
    for determining adsorption energies on solid surfaces. Npj
    Computational Materials, 3(1), 14.
    https://doi.org/10.1038/s41524-017-0017-z
```

Todo:
- Still assumes individual elements have their own chempots

> in a molecular adsorbate instead of considering a single
> chempot for a single molecular adsorbate. E.g. for an OH
> adsorbate, the surface energy is a function of delu_O and
> delu_H instead of delu_OH


* Need a method to automatically get chempot range when

    dealing with non-stoichiometric slabs


* Simplify the input for SurfaceEnergyPlotter such that the

    user does not need to generate a dict


### _class_ NanoscaleStability(se_analyzers, symprec=1e-05)
Bases: `object`

A class for analyzing the stability of nanoparticles of different
polymorphs with respect to size. The Wulff shape will be the model for the
nanoparticle. Stability will be determined by an energetic competition between the
weighted surface energy (surface energy of the Wulff shape) and the bulk energy. A
future release will include a 2D phase diagram (e.g. wrt size vs chempot for adsorbed
or non-stoichiometric surfaces). Based on the following work:

Kang, S., Mo, Y., Ong, S. P., & Ceder, G. (2014). Nanoscale

    stabilization of sodium oxides: Implications for Na-O2
    batteries. Nano Letters, 14(2), 1016-1020.
    [https://doi.org/10.1021/nl404557w](https://doi.org/10.1021/nl404557w)


#### se_analyzers()
Each item corresponds to a different polymorph.


* **Type**

    list[SurfaceEnergyPlotter]



#### symprec()
Tolerance for symmetry finding. See WulffShape.


* **Type**

    float


Analyzes the nanoscale stability of different polymorphs.


#### _static_ bulk_gform(bulk_entry)
Returns the formation energy of the bulk.


* **Parameters**

    **bulk_entry** ([*ComputedStructureEntry*](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)) – Entry of the corresponding bulk.



* **Returns**

    bulk formation energy (in eV)



* **Return type**

    float



#### plot_all_stability_map(max_r, increments=50, delu_dict=None, delu_default=0, plt=None, labels=None, from_sphere_area=False, e_units='keV', r_units='nanometers', normalize=False, scale_per_atom=False)
Returns the plot of the formation energy of a particles

    of different polymorphs against its effect radius.


* **Parameters**


    * **max_r** (*float*) – The maximum radius of the particle to plot up to.


    * **increments** (*int*) – Number of plot points


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **plt** (*pyplot*) – Plot


    * **labels** (*list*) – List of labels for each plot, corresponds to the
    list of se_analyzers


    * **from_sphere_area** (*bool*) – There are two ways to calculate the bulk
    formation energy. Either by treating the volume and thus surface
    area of the particle as a perfect sphere, or as a Wulff shape.



#### plot_one_stability_map(analyzer, max_r, delu_dict=None, label='', increments=50, delu_default=0, plt=None, from_sphere_area=False, e_units='keV', r_units='nanometers', normalize=False, scale_per_atom=False)
Returns the plot of the formation energy of a particle against its

    effect radius.


* **Parameters**


    * **analyzer** (*SurfaceEnergyPlotter*) – Analyzer associated with the
    first polymorph


    * **max_r** (*float*) – The maximum radius of the particle to plot up to.


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **label** (*str*) – Label of the plot for legend


    * **increments** (*int*) – Number of plot points


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **plt** (*pyplot*) – Plot


    * **from_sphere_area** (*bool*) – There are two ways to calculate the bulk
    formation energy. Either by treating the volume and thus surface
    area of the particle as a perfect sphere, or as a Wulff shape.


    * **r_units** (*str*) – Can be nanometers or Angstrom


    * **e_units** (*str*) – Can be keV or eV


    * **normalize** (*str*) – Whether or not to normalize energy by volume



#### scaled_wulff(wulffshape, r)
Scales the Wulff shape with an effective radius r. Note that the resulting

    Wulff does not necessarily have the same effective radius as the one
    provided. The Wulff shape is scaled by its surface energies where first
    the surface energies are scale by the minimum surface energy and then
    multiplied by the given effective radius.


* **Parameters**


    * **wulffshape** (*WulffShape*) – Initial, unscaled WulffShape


    * **r** (*float*) – Arbitrary effective radius of the WulffShape



* **Returns**

    WulffShape (scaled by r)



#### solve_equilibrium_point(analyzer1, analyzer2, delu_dict=None, delu_default=0, units='nanometers')
Gives the radial size of two particles where equilibrium is reached

    between both particles. NOTE: the solution here is not the same
    as the solution visualized in the plot because solving for r
    requires that both the total surface area and volume of the
    particles are functions of r.


* **Parameters**


    * **analyzer1** (*SurfaceEnergyPlotter*) – Analyzer associated with the
    first polymorph


    * **analyzer2** (*SurfaceEnergyPlotter*) – Analyzer associated with the
    second polymorph


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **units** (*str*) – Can be nanometers or Angstrom



* **Returns**

    Particle radius in nm



#### wulff_gform_and_r(wulffshape, bulk_entry, r, from_sphere_area=False, r_units='nanometers', e_units='keV', normalize=False, scale_per_atom=False)
Calculates the formation energy of the particle with arbitrary radius r.


* **Parameters**


    * **wulffshape** (*WulffShape*) – Initial, unscaled WulffShape


    * **bulk_entry** ([*ComputedStructureEntry*](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)) – Entry of the corresponding bulk.


    * **r** (*float** (**Ang**)*) – Arbitrary effective radius of the WulffShape


    * **from_sphere_area** (*bool*) – There are two ways to calculate the bulk
    formation energy. Either by treating the volume and thus surface
    area of the particle as a perfect sphere, or as a Wulff shape.


    * **r_units** (*str*) – Can be nanometers or Angstrom


    * **e_units** (*str*) – Can be keV or eV


    * **normalize** (*bool*) – Whether or not to normalize energy by volume


    * **scale_per_atom** (*True*) – Whether or not to normalize by number of
    atoms in the particle



* **Returns**

    particle formation energy (float in keV), effective radius



### _class_ SlabEntry(structure, energy, miller_index, correction=0.0, parameters=None, data=None, entry_id=None, label=None, adsorbates=None, clean_entry=None, marker=None, color=None)
Bases: [`ComputedStructureEntry`](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)

A ComputedStructureEntry object encompassing all data relevant to a

    slab for analyzing surface thermodynamics.


#### miller_index()
Miller index of plane parallel to surface.


* **Type**

    tuple



#### label()
Brief description for this slab.


* **Type**

    str



#### adsorbates()
List of ComputedStructureEntry for the types of adsorbates.


* **Type**

    list



#### clean_entry()
SlabEntry for the corresponding clean slab for an adsorbed slab.


* **Type**

    SlabEntry



#### ads_entries_dict()
Dictionary where the key is the reduced composition of the
adsorbate entry and value is the entry itself.


* **Type**

    dict


Make a SlabEntry containing all relevant surface thermodynamics data.


* **Parameters**


    * **structure** ([*Slab*](pymatgen.core.md#pymatgen.core.surface.Slab)) – The primary slab associated with this entry.


    * **energy** (*float*) – Energy from total energy calculation


    * **miller_index** (*tuple**(**h**, **k**, **l**)*) – Miller index of plane parallel
    to surface


    * **correction** (*float*) – See ComputedSlabEntry


    * **parameters** (*dict*) – See ComputedSlabEntry


    * **data** (*dict*) – See ComputedSlabEntry


    * **entry_id** (*str*) – See ComputedSlabEntry


    * **data** – See ComputedSlabEntry


    * **entry_id** – See ComputedSlabEntry


    * **label** (*str*) – Any particular label for this slab, e.g. “Tasker 2”,
    “non-stoichiometric”, “reconstructed”


    * **adsorbates** (*[*[*ComputedStructureEntry*](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)*]*) – List of reference entries
    for the adsorbates on the slab, can be an isolated molecule
    (e.g. O2 for O or O2 adsorption), a bulk structure (eg. fcc
    Cu for Cu adsorption) or anything.


    * **clean_entry** ([*ComputedStructureEntry*](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)) – If the SlabEntry is for an
    adsorbed slab, this is the corresponding SlabEntry for the
    clean slab


    * **marker** (*str*) – Custom marker for gamma plots (”–” and “-” are typical)


    * **color** (*str** or **rgba*) – Custom color for gamma plots



#### _property_ Nads_in_slab()
Returns the TOTAL number of adsorbates in the slab on BOTH sides.


#### _property_ Nsurfs_ads_in_slab()
Returns the TOTAL number of adsorbed surfaces in the slab.


#### _abc_impl(_ = <_abc._abc_data object_ )

#### as_dict()
Returns dict which contains Slab Entry data.


#### _property_ cleaned_up_slab()
Returns a slab with the adsorbates removed.


#### _property_ create_slab_label()
Returns a label (str) for this particular slab based on composition, coverage and Miller index.


#### _classmethod_ from_computed_structure_entry(entry, miller_index, label=None, adsorbates=None, clean_entry=None, \*\*kwargs)
Returns SlabEntry from a ComputedStructureEntry.


#### _classmethod_ from_dict(dct)
Returns a SlabEntry by reading in an dictionary.


#### _property_ get_monolayer()
Returns the primitive unit surface area density of the
adsorbate.


#### _property_ get_unit_primitive_area()
Returns the surface area of the adsorbed system per
unit area of the primitive slab system.


#### gibbs_binding_energy(eads=False)
Returns the adsorption energy or Gibbs binding energy of an adsorbate on a surface.


* **Parameters**

    **eads** (*bool*) – Whether to calculate the adsorption energy
    (True) or the binding energy (False) which is just
    adsorption energy normalized by number of adsorbates.



#### _property_ surface_area()
Calculates the surface area of the slab.


#### surface_energy(ucell_entry, ref_entries=None)
Calculates the surface energy of this SlabEntry.


* **Parameters**


    * **ucell_entry** (*entry*) – An entry object for the bulk


    * **(****list** (*ref_entries*) – [entry]): A list of entries for each type
    of element to be used as a reservoir for non-stoichiometric
    systems. The length of this list MUST be n-1 where n is the
    number of different elements in the bulk entry. The chempot
    of the element ref_entry that is not in the list will be
    treated as a variable.


Returns (Add (Sympy class)): Surface energy


### _class_ SurfaceEnergyPlotter(all_slab_entries, ucell_entry, ref_entries=None)
Bases: `object`

A class used for generating plots to analyze the thermodynamics of surfaces
of a material. Produces stability maps of different slab configurations,
phases diagrams of two parameters to determine stability of configurations
(future release), and Wulff shapes.


#### all_slab_entries()
Either a list of SlabEntry objects (note for a list, the
SlabEntry must have the adsorbates and clean_entry parameter plugged in) or a Nested
dictionary containing a list of entries for slab calculations as
items and the corresponding Miller index of the slab as the key.
To account for adsorption, each value is a sub-dictionary with the
entry of a clean slab calculation as the sub-key and a list of
entries for adsorption calculations as the sub-value. The sub-value
can contain different adsorption configurations such as a different
site or a different coverage, however, ordinarily only the most stable
configuration for a particular coverage will be considered as the
function of the adsorbed surface energy has an intercept dependent on
the adsorption energy (ie an adsorption site with a higher adsorption
energy will always provide a higher surface energy than a site with a
lower adsorption energy). An example parameter is provided:
{(h1,k1,l1): {clean_entry1: [ads_entry1, ads_entry2, …], clean_entry2: […], …}, (h2,k2,l2): {…}}
where clean_entry1 can be a pristine surface and clean_entry2 can be a
reconstructed surface while ads_entry1 can be adsorption at site 1 with
a 2x2 coverage while ads_entry2 can have a 3x3 coverage. If adsorption
entries are present (i.e. if all_slab_entries[(h,k,l)][clean_entry1]), we
consider adsorption in all plots and analysis for this particular facet.


* **Type**

    dict | list



#### color_dict()
Dictionary of colors (r,g,b,a) when plotting surface energy stability.
The keys are individual surface entries where clean surfaces have a solid color while
the corresponding adsorbed surface will be transparent.


* **Type**

    dict



#### ucell_entry()
ComputedStructureEntry of the bulk reference for
this particular material.


* **Type**

    [ComputedStructureEntry](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)



#### ref_entries()
List of ComputedStructureEntries to be used for calculating chemical potential.


* **Type**

    list



#### facet_color_dict()
Randomly generated dictionary of colors associated with each facet.


* **Type**

    dict


Object for plotting surface energy in different ways for clean and

    adsorbed surfaces.


* **Parameters**


    * **all_slab_entries** (*dict** or **list*) – Dictionary or list containing
    all entries for slab calculations. See attributes.


    * **ucell_entry** ([*ComputedStructureEntry*](pymatgen.entries.md#pymatgen.entries.computed_entries.ComputedStructureEntry)) – ComputedStructureEntry
    of the bulk reference for this particular material.


    * **ref_entries** (*[**ComputedStructureEntries**]*) – A list of entries for
    each type of element to be used as a reservoir for
    non-stoichiometric systems. The length of this list MUST be
    n-1 where n is the number of different elements in the bulk
    entry. The bulk energy term in the grand surface potential can
    be defined by a summation of the chemical potentials for each
    element in the system. As the bulk energy is already provided,
    one can solve for one of the chemical potentials as a function
    of the other chemical potentials and bulk energy. i.e. there
    are n-1 variables (chempots). e.g. if your ucell_entry is for
    LiFePO4 than your ref_entries should have an entry for Li, Fe,
    and P if you want to use the chempot of O as the variable.



#### BE_vs_clean_SE(delu_dict, delu_default=0, plot_eads=False, annotate_monolayer=True, JPERM2=False)
For each facet, plot the clean surface energy against the most

    stable binding energy.


* **Parameters**


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **plot_eads** (*bool*) – Option to plot the adsorption energy (binding
    energy multiplied by number of adsorbates) instead.


    * **annotate_monolayer** (*bool*) – Whether or not to label each data point
    with its monolayer (adsorbate density per unit primiitve area)


    * **JPERM2** (*bool*) – Whether to plot surface energy in /m^2 (True) or
    eV/A^2 (False)



* **Returns**

    Plot of clean surface energy vs binding energy for

        all facets.




* **Return type**

    (Plot)



#### area_frac_vs_chempot_plot(ref_delu: Symbol, chempot_range: list[float], delu_dict: dict[Symbol, float] | None = None, delu_default: float = 0, increments: int = 10, no_clean: bool = False, no_doped: bool = False)
1D plot. Plots the change in the area contribution
of each facet as a function of chemical potential.


* **Parameters**


    * **ref_delu** (*Symbol*) – The free variable chempot with the format:
    Symbol(“delu_el”) where el is the name of the element.


    * **chempot_range** (*list**[**float**]*) – Min/max range of chemical potential to plot along.


    * **delu_dict** (*dict**[**Symbol**, **float**]*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials.


    * **increments** (*int*) – Number of data points between min/max or point
    of intersection. Defaults to 10 points.


    * **no_clean** (*bool*) – Some parameter, description missing.


    * **no_doped** (*bool*) – Some parameter, description missing.



* **Returns**

    Plot of area frac on the Wulff shape for each facet vs chemical potential.



* **Return type**

    plt.Axes



#### _static_ chempot_plot_addons(ax, xrange, ref_el, pad=2.4, rect=None, ylim=None)
Helper function to a chempot plot look nicer.


* **Parameters**


    * **plt** (*Plot*) –


    * **xrange** (*list*) – xlim parameter


    * **ref_el** (*str*) – Element of the referenced chempot.


    * **axes** (*axes*) –


    * **pad** (*float*) –


    * **rect** (*list*) – For tight layout


    * **ylim** (*ylim parameter*) –


return (Plot): Modified plot with addons.
return (Plot): Modified plot with addons.


#### chempot_vs_gamma(ref_delu, chempot_range, miller_index=(), delu_dict=None, delu_default=0, JPERM2=False, show_unstable=False, ylim=None, plt=None, no_clean=False, no_doped=False, use_entry_labels=False, no_label=False)
Plots the surface energy as a function of chemical potential.

    Each facet will be associated with its own distinct colors.
    Dashed lines will represent stoichiometries different from that
    of the mpid’s compound. Transparent lines indicates adsorption.


* **Parameters**


    * **ref_delu** (*sympy Symbol*) – The range stability of each slab is based
    on the chempot range of this chempot. Should be a sympy Symbol
    object of the format: Symbol(“delu_el”) where el is the name of
    the element


    * **chempot_range** (*[**max_chempot**, **min_chempot**]*) – Range to consider the
    stability of the slabs.


    * **miller_index** (*list*) – Miller index for a specific facet to get a
    dictionary for.


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **JPERM2** (*bool*) – Whether to plot surface energy in /m^2 (True) or
    eV/A^2 (False)


    * **show_unstable** (*bool*) – Whether or not to show parts of the surface
    energy plot outside the region of stability.


    * **ylim** (*[**ymax**, **ymin**]*) – Range of y axis


    * **no_doped** (*bool*) – Whether to plot for the clean slabs only.


    * **no_clean** (*bool*) – Whether to plot for the doped slabs only.


    * **use_entry_labels** (*bool*) – If True, will label each slab configuration
    according to their given label in the SlabEntry object.


    * **no_label** (*bool*) – Option to turn off labels.



* **Returns**

    Plot of surface energy vs chempot for all entries.



* **Return type**

    (Plot)



#### chempot_vs_gamma_plot_one(ax: plt.Axes, entry: SlabEntry, ref_delu: Symbol, chempot_range: list[float], delu_dict: dict[Symbol, float] | None = None, delu_default: float = 0, label: str = '', JPERM2: bool = False)
Helper function to help plot the surface energy of a
single SlabEntry as a function of chemical potential.


* **Parameters**


    * **ax** (*plt.Axes*) – Matplotlib Axes instance for plotting.


    * **entry** – Entry of the slab whose surface energy we want
    to plot. (Add appropriate description for type)


    * **ref_delu** (*Symbol*) – The range stability of each slab is based
    on the chempot range of this chempot.


    * **chempot_range** (*list**[**float**]*) – Range to consider the stability of the slabs.


    * **delu_dict** (*dict**[**Symbol**, **float**]*) – Dictionary of the chemical potentials.


    * **delu_default** (*float*) – Default value for all unset chemical potentials.


    * **label** (*str*) – Label of the slab for the legend.


    * **JPERM2** (*bool*) – Whether to plot surface energy in /m^2 (True) or
    eV/A^2 (False).



* **Returns**

    Plot of surface energy vs chemical potential for one entry.



* **Return type**

    plt.Axes



#### color_palette_dict(alpha=0.35)
Helper function to assign each facet a unique color using a dictionary.


* **Parameters**

    **alpha** (*float*) – Degree of transparency


return (dict): Dictionary of colors (r,g,b,a) when plotting surface

    energy stability. The keys are individual surface entries where
    clean surfaces have a solid color while the corresponding adsorbed
    surface will be transparent.


#### get_stable_entry_at_u(miller_index, delu_dict=None, delu_default=0, no_doped=False, no_clean=False)
Returns the entry corresponding to the most stable slab for a particular

    facet at a specific chempot. We assume that surface energy is constant
    so all free variables must be set with delu_dict, otherwise they are
    assumed to be equal to delu_default.


* **Parameters**


    * **miller_index** (*(**h**,**k**,**l**)*) – The facet to find the most stable slab in


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **no_doped** (*bool*) – Consider stability of clean slabs only.


    * **no_clean** (*bool*) – Consider stability of doped slabs only.



* **Returns**

    SlabEntry, surface_energy (float)



#### get_surface_equilibrium(slab_entries, delu_dict=None)
Takes in a list of SlabEntries and calculates the chemical potentials

    at which all slabs in the list coexists simultaneously. Useful for
    building surface phase diagrams. Note that to solve for x equations
    (x slab_entries), there must be x free variables (chemical potentials).
    Adjust delu_dict as need be to get the correct number of free variables.


* **Parameters**


    * **slab_entries** (*array*) – The coefficients of the first equation


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.



* **Returns**

    Array containing a solution to x equations with x

        variables (x-1 chemical potential and 1 surface energy)




* **Return type**

    (array)



#### monolayer_vs_BE(plot_eads=False)
Plots the binding energy as a function of monolayers (ML), i.e.

    the fractional area adsorbate density for all facets. For each
    facet at a specific monolayer, only plot the lowest binding energy.


* **Parameters**

    **plot_eads** (*bool*) – Option to plot the adsorption energy (binding
    energy multiplied by number of adsorbates) instead.



* **Returns**

    Plot of binding energy vs monolayer for all facets.



* **Return type**

    (Plot)



#### set_all_variables(delu_dict, delu_default)
Sets all chemical potential values and returns a dictionary where

    the key is a sympy Symbol and the value is a float (chempot).


* **Parameters**


    * **entry** (*SlabEntry*) – Computed structure entry of the slab


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials



* **Returns**

    Dictionary of set chemical potential values



#### stable_u_range_dict(chempot_range, ref_delu, no_doped=True, no_clean=False, delu_dict=None, miller_index=(), dmu_at_0=False, return_se_dict=False)
Creates a dictionary where each entry is a key pointing to a
chemical potential range where the surface of that entry is stable.
Does so by enumerating through all possible solutions (intersect)
for surface energies of a specific facet.


* **Parameters**


    * **chempot_range** (*[**max_chempot**, **min_chempot**]*) – Range to consider the
    stability of the slabs.


    * **ref_delu** (*sympy Symbol*) – The range stability of each slab is based
    on the chempot range of this chempot. Should be a sympy Symbol
    object of the format: Symbol(“delu_el”) where el is the name of
    the element


    * **no_doped** (*bool*) – Consider stability of clean slabs only.


    * **no_clean** (*bool*) – Consider stability of doped slabs only.


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **miller_index** (*list*) – Miller index for a specific facet to get a
    dictionary for.


    * **dmu_at_0** (*bool*) – If True, if the surface energies corresponding to
    the chemical potential range is between a negative and positive
    value, the value is a list of three chemical potentials with the
    one in the center corresponding a surface energy of 0. Uselful
    in identifying unphysical ranges of surface energies and their
    chemical potential range.


    * **return_se_dict** (*bool*) – Whether or not to return the corresponding
    dictionary of surface energies



#### surface_chempot_range_map(elements, miller_index, ranges, incr=50, no_doped=False, no_clean=False, delu_dict=None, ax=None, annotate=True, show_unphyiscal_only=False, fontsize=10)
Adapted from the get_chempot_range_map() method in the PhaseDiagram

    class. Plot the chemical potential range map based on surface
    energy stability. Currently works only for 2-component PDs. At
    the moment uses a brute force method by enumerating through the
    range of the first element chempot with a specified increment
    and determines the chempot rangeo fht e second element for each
    SlabEntry. Future implementation will determine the chempot range
    map first by solving systems of equations up to 3 instead of 2.


* **Parameters**


    * **elements** (*list*) – Sequence of elements to be considered as independent
    variables. E.g., if you want to show the stability ranges of
    all Li-Co-O phases wrt to duLi and duO, you will supply
    [Element(“Li”), Element(“O”)]


    * **miller_index** (*[**h**, **k**, **l**]*) – Miller index of the surface we are interested in


    * **ranges** (*[**[**range1**]**, **[**range2**]**]*) – List of chempot ranges (max and min values)
    for the first and second element.


    * **incr** (*int*) – Number of points to sample along the range of the first chempot


    * **no_doped** (*bool*) – Whether or not to include doped systems.


    * **no_clean** (*bool*) – Whether or not to include clean systems.


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **ax** (*plt.Axes*) – Axes object to plot on. If None, will create a new plot.


    * **annotate** (*bool*) – Whether to annotate each “phase” with the label of
    the entry. If no label, uses the reduced formula


    * **show_unphyiscal_only** (*bool*) – Whether to only show the shaded region where
    surface energy is negative. Useful for drawing other chempot range maps.


    * **fontsize** (*int*) – Font size of the annotation



#### wulff_from_chempot(delu_dict=None, delu_default=0, symprec=1e-05, no_clean=False, no_doped=False)
Method to get the Wulff shape at a specific chemical potential.


* **Parameters**


    * **delu_dict** (*dict*) – Dictionary of the chemical potentials to be set as
    constant. Note the key should be a sympy Symbol object of the
    format: Symbol(“delu_el”) where el is the name of the element.


    * **delu_default** (*float*) – Default value for all unset chemical potentials


    * **symprec** (*float*) – See WulffShape.


    * **no_doped** (*bool*) – Consider stability of clean slabs only.


    * **no_clean** (*bool*) – Consider stability of doped slabs only.



* **Returns**

    The WulffShape at u_ref and u_ads.



* **Return type**

    (WulffShape)



### _class_ WorkFunctionAnalyzer(structure: [Structure](pymatgen.core.md#pymatgen.core.structure.Structure), locpot_along_c, efermi, shift=0, blength=3.5)
Bases: `object`

A class used for calculating the work function from a slab model and
visualizing the behavior of the local potential along the slab.


#### efermi()
The Fermi energy.


* **Type**

    float



#### locpot_along_c()
Local potential in eV along points along the c axis.


* **Type**

    list



#### vacuum_locpot()
The maximum local potential along the c direction for the slab model,
i.e. the potential at the vacuum.


* **Type**

    float



#### work_function()
The minimum energy needed to move an electron from the surface to infinity.
Defined as the difference between the potential at the vacuum and the Fermi energy.


* **Type**

    float



#### slab()
The slab structure model.


* **Type**

    [Slab](pymatgen.core.md#pymatgen.core.surface.Slab)



#### along_c()
Points along the c direction with same increments as the locpot in the c axis.


* **Type**

    list



#### ave_locpot()
Mean of the minimum and maximum (vacuum) locpot along c.


* **Type**

    float



#### sorted_sites()
List of sites from the slab sorted along the c direction.


* **Type**

    list



#### ave_bulk_p()
The average locpot of the slab region along the c direction.


* **Type**

    float


Initializes the WorkFunctionAnalyzer class.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)) – Structure object modelling the surface


    * **locpot_along_c** (*list*) – Local potential along the c direction


    * **outcar** (*MSONable*) – Outcar vasp output object


    * **shift** (*float*) – Parameter to translate the slab (and
    therefore the vacuum) of the slab structure, thereby
    translating the plot along the x axis.


    * **blength** (*float** (**Ang**)*) – The longest bond length in the material.
    Used to handle pbc for noncontiguous slab layers



#### _classmethod_ from_files(poscar_filename, locpot_filename, outcar_filename, shift=0, blength=3.5)
Initializes a WorkFunctionAnalyzer from POSCAR, LOCPOT, and OUTCAR files.


* **Parameters**


    * **poscar_filename** (*str*) – The path to the POSCAR file.


    * **locpot_filename** (*str*) – The path to the LOCPOT file.


    * **outcar_filename** (*str*) – The path to the OUTCAR file.


    * **shift** (*float*) – The shift value. Defaults to 0.


    * **blength** (*float*) – The longest bond length in the material.
    Used to handle pbc for noncontiguous slab layers. Defaults to 3.5.



* **Returns**

    A WorkFunctionAnalyzer instance.



* **Return type**

    WorkFunctionAnalyzer



#### get_labels(plt, label_fontsize=10)
Handles the optional labelling of the plot with relevant quantities


* **Parameters**


    * **plt** (*plt*) – Plot of the locpot vs c axis


    * **label_fontsize** (*float*) – Fontsize of labels


Returns Labelled plt.


#### get_locpot_along_slab_plot(label_energies=True, plt=None, label_fontsize=10)
Returns a plot of the local potential (eV) vs the

    position along the c axis of the slab model (Ang).


* **Parameters**


    * **label_energies** (*bool*) – Whether to label relevant energy
    quantities such as the work function, Fermi energy,
    vacuum locpot, bulk-like locpot


    * **plt** (*plt*) – Matplotlib pyplot object


    * **label_fontsize** (*float*) – Fontsize of labels


Returns plt of the locpot vs c axis


#### is_converged(min_points_frac=0.015, tol: float = 0.0025)
A well converged work function should have a flat electrostatic

    potential within some distance (min_point) about where the peak
    electrostatic potential is found along the c direction of the
    slab. This is dependent on the size of the slab.


* **Parameters**


    * **min_point** (*fractional coordinates*) – The number of data points
    +/- the point of where the electrostatic potential is at
    its peak along the c direction.


    * **tol** (*float*) – If the electrostatic potential stays the same
    within this tolerance, within the min_points, it is converged.


Returns a bool (whether or not the work function is converged)


### entry_dict_from_list(all_slab_entries)
Converts a list of SlabEntry to an appropriate dictionary. It is
assumed that if there is no adsorbate, then it is a clean SlabEntry
and that adsorbed SlabEntry has the clean_entry parameter set.


* **Parameters**

    **all_slab_entries** (*list*) – List of SlabEntry objects



* **Returns**

    Dictionary of SlabEntry with the Miller index as the main

        key to a dictionary with a clean SlabEntry as the key to a
        list of adsorbed SlabEntry.




* **Return type**

    (dict)



### sub_chempots(gamma_dict, chempots)
Uses dot product of numpy array to sub chemical potentials

    into the surface grand potential. This is much faster
    than using the subs function in sympy.


* **Parameters**


    * **gamma_dict** (*dict*) – Surface grand potential equation
    as a coefficient dictionary


    * **chempots** (*dict*) – Dictionary assigning each chemical
    potential (key) in gamma a value



* **Returns**

    Surface energy as a float


## pymatgen.analysis.thermochemistry module

A module to perform experimental thermochemical data analysis.


### _class_ ThermoData(data_type, cpdname, phaseinfo, formula, value, ref='', method='', temp_range=(298, 298), uncertainty=None)
Bases: `object`

A object container for an experimental Thermochemical Data.


* **Parameters**


    * **data_type** – The thermochemical data type. Should be one of the
    following: fH - Formation enthalpy, S - Entropy,
    A, B, C, D, E, F, G, H - variables for use in the various
    equations for generating formation enthalpies or Cp at
    various temperatures.


    * **cpdname** (*str*) – A name for the compound. For example, hematite for
    Fe2O3.


    * **phaseinfo** (*str*) – Denoting the phase. For example, “solid”, “liquid”,
    “gas” or “tetragonal”.


    * **formula** (*str*) – A proper string formula, e.g., Fe2O3


    * **value** (*float*) – The value of the data.


    * **ref** (*str*) – A reference, if any, for the data.


    * **method** (*str*) – The method by which the data was determined,
    if available.


    * **temp_range** (*[**float**, **float**]*) – Temperature range of validity for the
    data in Kelvin. Defaults to 298 K only.


    * **uncertainty** (*float*) – An uncertainty for the data, if available.



#### as_dict()
Returns: MSONable dict.


#### _classmethod_ from_dict(d)

* **Parameters**

    **d** (*dict*) – Dict representation.



* **Returns**

    ThermoData


## pymatgen.analysis.transition_state module

Some reimplementation of Henkelman’s Transition State Analysis utilities,
which are originally in Perl. Additional features beyond those offered by
Henkelman’s utilities will be added.

This allows the usage and customization in Python.


### _class_ NEBAnalysis(r, energies, forces, structures, spline_options=None)
Bases: `MSONable`

An NEBAnalysis class.

Initializes an NEBAnalysis from the cumulative root mean squared distances
between structures, the energies, the forces, the structures and the
interpolation_order for the analysis.


* **Parameters**


    * **r** – Root mean square distances between structures


    * **energies** – Energies of each structure along reaction coordinate


    * **forces** – Tangent forces along the reaction coordinate.


    * **structures** (*[*[*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)*]*) – List of Structures along reaction
    coordinate.


    * **spline_options** (*dict*) – Options for cubic spline. For example,
    {“saddle_point”: “zero_slope”} forces the slope at the saddle to
    be zero.



#### as_dict()
Dict representation of NEBAnalysis.


* **Returns**

    JSON-serializable dict representation.



#### _classmethod_ from_dir(root_dir, relaxation_dirs=None, \*\*kwargs)
Initializes a NEBAnalysis object from a directory of a NEB run.
Note that OUTCARs must be present in all image directories. For the
terminal OUTCARs from relaxation calculations, you can specify the
locations using relaxation_dir. If these are not specified, the code
will attempt to look for the OUTCARs in 00 and 0n directories,
followed by subdirs “start”, “end” or “initial”, “final” in the
root_dir. These are just some typical conventions used
preferentially in Shyue Ping’s MAVRL research group. For the
non-terminal points, the CONTCAR is read to obtain structures. For
terminal points, the POSCAR is used. The image directories are
assumed to be the only directories that can be resolved to integers.
E.g., “00”, “01”, “02”, “03”, “04”, “05”, “06”. The minimum
sub-directory structure that can be parsed is of the following form (
a 5-image example is shown):

00:
- POSCAR
- OUTCAR
01, 02, 03, 04, 05:
- CONTCAR
- OUTCAR
06:
- POSCAR
- OUTCAR


* **Parameters**


    * **root_dir** (*str*) – Path to the root directory of the NEB calculation.


    * **relaxation_dirs** (*tuple*) – This specifies the starting and ending
    relaxation directories from which the OUTCARs are read for the
    terminal points for the energies.



* **Returns**

    NEBAnalysis object.



#### _classmethod_ from_outcars(outcars, structures, \*\*kwargs)
Initializes an NEBAnalysis from Outcar and Structure objects. Use
the static constructors, e.g., from_dir instead if you
prefer to have these automatically generated from a directory of NEB
calculations.


* **Parameters**


    * **outcars** (*[*[*Outcar*](pymatgen.io.vasp.md#pymatgen.io.vasp.outputs.Outcar)*]*) – List of Outcar objects. Note that these have
    to be ordered from start to end along reaction coordinates.


    * **structures** (*[*[*Structure*](pymatgen.core.md#pymatgen.core.structure.Structure)*]*) – List of Structures along reaction
    coordinate. Must be same length as outcar.


    * **interpolation_order** (*int*) – Order of polynomial to use to
    interpolate between images. Same format as order parameter in
    scipy.interplotate.PiecewisePolynomial.



#### get_extrema(normalize_rxn_coordinate=True)
Returns the positions of the extrema along the MEP. Both local
minimums and maximums are returned.


* **Parameters**

    **normalize_rxn_coordinate** (*bool*) – Whether to normalize the
    reaction coordinate to between 0 and 1. Defaults to True.



* **Returns**

    (min_extrema, max_extrema), where the extrema are given as
    [(x1, y1), (x2, y2), …].



#### get_plot(normalize_rxn_coordinate: bool = True, label_barrier: bool = True)
Returns the NEB plot. Uses Henkelman’s approach of spline fitting
each section of the reaction path based on tangent force and energies.


* **Parameters**


    * **normalize_rxn_coordinate** (*bool*) – Whether to normalize the
    reaction coordinate to between 0 and 1. Defaults to True.


    * **label_barrier** (*bool*) – Whether to label the maximum barrier.



* **Returns**

    matplotlib axes object.



* **Return type**

    plt.Axes



#### setup_spline(spline_options=None)
Setup of the options for the spline interpolation.


* **Parameters**

    **spline_options** (*dict*) – Options for cubic spline. For example,
    {“saddle_point”: “zero_slope”} forces the slope at the saddle to
    be zero.



### combine_neb_plots(neb_analyses, arranged_neb_analyses=False, reverse_plot=False)
neb_analyses: a list of NEBAnalysis objects.

arranged_neb_analyses: The code connects two end points with the
smallest-energy difference. If all end points have very close energies, it’s
likely to result in an inaccurate connection. Manually arrange neb_analyses
if the combined plot is not as expected compared with all individual plots.
E.g., if there are two NEBAnalysis objects to combine, arrange in such a
way that the end-point energy of the first NEBAnalysis object is the
start-point energy of the second NEBAnalysis object.
Note that the barrier labeled in y-axis in the combined plot might be
different from that in the individual plot due to the reference energy used.
reverse_plot: reverse the plot or percolation direction.
return: a NEBAnalysis object

## pymatgen.analysis.wulff module

This module define a WulffShape class to generate the Wulff shape from
a lattice, a list of indices and their corresponding surface energies,
and the total area and volume of the Wulff shape, the weighted surface energy,
the anisotropy and shape_factor can also be calculated.
In support of plotting from a given view in terms of miller index.

The lattice is from the conventional unit cell, and (hkil) for hexagonal
lattices.

If you use this code extensively, consider citing the following:

Tran, R.; Xu, Z.; Radhakrishnan, B.; Winston, D.; Persson, K. A.; Ong, S. P.
(2016). Surface energies of elemental crystals. Scientific Data.


### _class_ WulffFacet(normal, e_surf, normal_pt, dual_pt, index, m_ind_orig, miller)
Bases: `object`

Helper container for each Wulff plane.


* **Parameters**


    * **normal** –


    * **e_surf** –


    * **normal_pt** –


    * **dual_pt** –


    * **index** –


    * **m_ind_orig** –


    * **miller** –



### _class_ WulffShape(lattice: [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice), miller_list, e_surf_list, symprec=1e-05)
Bases: `object`

Generate Wulff Shape from list of miller index and surface energies,
with given conventional unit cell.
surface energy (Jm^2) is the length of normal.

Wulff shape is the convex hull.
Based on:
[http://scipy.github.io/devdocs/generated/scipy.spatial.ConvexHull.html](http://scipy.github.io/devdocs/generated/scipy.spatial.ConvexHull.html)

Process:


    1. get Wulff simplices


    2. label with color


    3. get wulff_area and other properties


#### debug()
Whether to print debug information.


* **Type**

    bool



#### alpha()
Transparency of the Wulff shape.


* **Type**

    float



#### color_set()
colors to use for facets.


* **Type**

    list



#### grid_off()
Whether to turn off the grid.


* **Type**

    bool



#### axis_off()
Whether to turn off the axis.


* **Type**

    bool



#### show_area()
Whether to show the area of each facet.


* **Type**

    bool



#### off_color()
Color of facets not on the Wulff shape.


* **Type**

    str



#### structure()
Input conventional unit cell (with H) from lattice.


* **Type**

    [Structure](pymatgen.core.md#pymatgen.core.structure.Structure)



#### miller_list()
input Miller indices, for hcp in the form of hkil.


* **Type**

    list



#### hkl_list()
Modified Miller indices in the same order as input_miller.


* **Type**

    list



#### e_surf_list()
input surface energies in the same order as input_miller.


* **Type**

    list



#### lattice()
Input lattice for the conventional unit cell.


* **Type**

    [Lattice](pymatgen.core.md#pymatgen.core.lattice.Lattice)



#### facets()
WulffFacet objects considering symmetry.


* **Type**

    list



#### dual_cv_simp()
Simplices from the dual convex hull (dual_pt).


* **Type**

    list



#### wulff_pt_list()
Wulff points.


* **Type**

    list



#### wulff_cv_simp()
Simplices from the convex hull of wulff_pt_list.


* **Type**

    list



#### on_wulff()
List for all input_miller, True if on the Wulff shape.


* **Type**

    list



#### color_area()
List for all input_miller, total area on the Wulff shape, off_wulff = 0.


* **Type**

    list



#### miller_area()
Dictionary of Miller indices and their corresponding areas.


* **Type**

    dict



* **Parameters**


    * **lattice** – Lattice object of the conventional unit cell


    * **miller_list** (*[**(**hkl*) – list of hkl or hkil for hcp


    * **e_surf_list** (*[**float**]*) – list of corresponding surface energies


    * **symprec** (*float*) – for reciprocal lattice operation, default is 1e-5.



#### _get_all_miller_e()
From self: get miller_list(unique_miller), e_surf_list and symmetry operations(symm_ops)
according to lattice apply symm_ops to get all the miller index, then get normal, get
all the facets functions for Wulff shape calculation:

```
|normal|
```

 = 1, e_surf is plane’s
distance to (0, 0, 0), normal[0]x + normal[1]y + normal[2]z = e_surf.


* **Returns**

    [WulffFacet]



#### _get_azimuth_elev(miller_index)

* **Parameters**

    **miller_index** – viewing direction.



* **Returns**

    azim, elev for plotting



#### _get_colors(color_set, alpha, off_color, custom_colors=None)
Assign colors according to the surface energies of on_wulff facets.


* **Returns**

    color_list, color_proxy, color_proxy_on_wulff, miller_on_wulff,
    e_surf_on_wulff_list



* **Return type**

    tuple



#### _get_cross_pt_dual_simp(dual_simp)


```
|normal|
```

 = 1, e_surf is plane’s distance to (0, 0, 0),
plane function:

> normal[0]x + normal[1]y + normal[2]z = e_surf.

from self:

    normal_e_m to get the plane functions
    dual_simp: (i, j, k) simplices from the dual convex hull

    > i, j, k: plane index(same order in normal_e_m)


#### _get_simpx_plane()
Locate the plane for simpx of on wulff_cv, by comparing the center of
the simpx triangle with the plane functions.


#### _property_ anisotropy()
Returns:
(float) Coefficient of Variation from weighted surface energy
The ideal sphere is 0.


#### _property_ area_fraction_dict()
Returns:
(dict): {hkl: area_hkl/total area on wulff}.


#### _property_ effective_radius()
Radius of the Wulffshape when the
Wulffshape is approximated as a sphere.


* **Returns**

    (float) radius.



#### get_line_in_facet(facet)
Returns the sorted pts in a facet used to draw a line.


#### get_plot(color_set='PuBu', grid_off=True, axis_off=True, show_area=False, alpha=1, off_color='red', direction=None, bar_pos=(0.75, 0.15, 0.05, 0.65), bar_on=False, units_in_JPERM2=True, legend_on=True, aspect_ratio=(8, 8), custom_colors=None)
Get the Wulff shape plot.


* **Parameters**


    * **color_set** – default is ‘PuBu’


    * **grid_off** (*bool*) – default is True


    * **axis_off** (*bool*) – default is True


    * **show_area** (*bool*) – default is False


    * **alpha** (*float*) – chosen from 0 to 1 (float), default is 1


    * **off_color** – Default color for facets not present on the Wulff shape.


    * **direction** – default is (1, 1, 1)


    * **bar_pos** – default is [0.75, 0.15, 0.05, 0.65]


    * **bar_on** (*bool*) – default is False


    * **legend_on** (*bool*) – default is True


    * **aspect_ratio** – default is (8, 8)


    * **(****{****(****h** (*custom_colors*) – [r,g,b,alpha]}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **k** – [r,g,b,alpha]}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **l}** – [r,g,b,alpha]}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **units_in_JPERM2** (*bool*) – Units of surface energy, defaults to
    Joules per square meter (True)



* **Returns**

    (matplotlib.pyplot)



#### get_plotly(color_set='PuBu', off_color='red', alpha=1, custom_colors=None, units_in_JPERM2=True)
Get the Wulff shape as a plotly Figure object.


* **Parameters**


    * **color_set** – default is ‘PuBu’


    * **alpha** (*float*) – chosen from 0 to 1 (float), default is 1


    * **off_color** – Default color for facets not present on the Wulff shape.


    * **(****{****(****h** (*custom_colors*) – [r,g,b,alpha}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **k** – [r,g,b,alpha}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **l}** – [r,g,b,alpha}): Customize color of each
    facet with a dictionary. The key is the corresponding Miller
    index and value is the color. Undefined facets will use default
    color site. Note: If you decide to set your own colors, it
    probably won’t make any sense to have the color bar on.


    * **units_in_JPERM2** (*bool*) – Units of surface energy, defaults to
    Joules per square meter (True)



* **Returns**

    (plotly.graph_objs.Figure)



#### _property_ miller_area_dict()
area_hkl on wulff}.


* **Type**

    Returns {hkl



#### _property_ miller_energy_dict()
surface energy_hkl}.


* **Type**

    Returns {hkl



#### _property_ shape_factor()
This is useful for determining the critical nucleus size.
A large shape factor indicates great anisotropy.
See Ballufi, R. W., Allen, S. M. & Carter, W. C. Kinetics

> of Materials. (John Wiley & Sons, 2005), p.461.


* **Returns**

    (float) Shape factor.



#### show(\*args, \*\*kwargs)
Show the Wulff plot.


* **Parameters**


    * **\*args** – Passed to get_plot.


    * **\*\*kwargs** – Passed to get_plot.



#### _property_ surface_area()
Total surface area of Wulff shape.


#### _property_ tot_corner_sites()
Returns the number of vertices in the convex hull.
Useful for identifying catalytically active sites.


#### _property_ tot_edges()
Returns the number of edges in the convex hull.
Useful for identifying catalytically active sites.


#### _property_ total_surface_energy()
Total surface energy of the Wulff shape.


* **Returns**

    (float) sum(surface_energy_hkl \* area_hkl)



#### _property_ volume()
Volume of the Wulff shape.


#### _property_ weighted_surface_energy()
Returns:
sum(surface_energy_hkl \* area_hkl)/ sum(area_hkl).


### get_tri_area(pts)
Given a list of coords for 3 points,
Compute the area of this triangle.


* **Parameters**

    **pts** – [a, b, c] three points



### hkl_tuple_to_str(hkl)
Prepare for display on plots “(hkl)” for surfaces


* **Parameters**

    **hkl** – in the form of [h, k, l] or (h, k, l).


## pymatgen.analysis.xps module

This is a module for XPS analysis. It is modelled after the Galore package ([https://github.com/SMTG-UCL/galore](https://github.com/SMTG-UCL/galore)), but
with some modifications for easier analysis from pymatgen itself. Please cite the following original work if you use
this:

```default
Adam J. Jackson, Alex M. Ganose, Anna Regoutz, Russell G. Egdell, David O. Scanlon (2018). Galore: Broadening and
weighting for simulation of photoelectron spectroscopy. Journal of Open Source Software, 3(26), 773,
doi: 10.21105/joss.007733
```

You may wish to look at the optional dependency galore for more functionality such as plotting and other cross-sections.
Note that the atomic_subshell_photoionization_cross_sections.csv has been reparsed from the original compilation:

```default
Yeh, J. J.; Lindau, I. Atomic Subshell Photoionization Cross Sections and Asymmetry Parameters: 1 ⩽ Z ⩽ 103.
Atomic Data and Nuclear Data Tables 1985, 32 (1), 1-155. https://doi.org/10.1016/0092-640X(85)90016-6.
```

This version contains all detailed information for all orbitals.


### _class_ XPS(x: ArrayLike, y: ArrayLike, \*args, \*\*kwargs)
Bases: [`Spectrum`](pymatgen.core.md#pymatgen.core.spectrum.Spectrum)

Class representing an X-ray photoelectron spectra.


* **Parameters**


    * **x** (*ndarray*) – A ndarray of N values.


    * **y** (*ndarray*) – A ndarray of N x k values. The first dimension must be
    the same as that of x. Each of the k values are interpreted as separate.


    * **\*args** – All subclasses should provide args other than x and y
    when calling super, e.g., super().__init__(
    x, y, arg1, arg2, kwarg1=val1, ..). This guarantees the +, -,

    ```
    *
    ```

    ,
    etc. operators work properly.


    * **\*\*kwargs** – Same as that for

    ```
    *
    ```

    args.



#### XLABEL(_ = 'Binding Energy (eV)_ )

#### YLABEL(_ = 'Intensity_ )

#### _classmethod_ from_dos(dos: [CompleteDos](pymatgen.electronic_structure.md#pymatgen.electronic_structure.dos.CompleteDos))

* **Parameters**


    * **dos** – CompleteDos object with project element-orbital DOS. Can be obtained from Vasprun.get_complete_dos.


    * **sigma** – Smearing for Gaussian.



* **Returns**

    XPS



### _load_cross_sections(fname)