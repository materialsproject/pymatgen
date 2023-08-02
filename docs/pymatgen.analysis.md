---
layout: default
title: pymatgen.analysis.md
nav_exclude: true
---

# pymatgen.analysis namespace

## Subpackages


* [pymatgen.analysis.chemenv package](pymatgen.analysis.chemenv.md)


    * [Subpackages](pymatgen.analysis.chemenv.md#subpackages)


        * [pymatgen.analysis.chemenv.connectivity package](pymatgen.analysis.chemenv.connectivity.md)




                * [pymatgen.analysis.chemenv.connectivity.connected_components module](pymatgen.analysis.chemenv.connectivity.connected_components.md)


                    * [`ConnectedComponent`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent)


                        * [`ConnectedComponent.as_dict()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.as_dict)


                        * [`ConnectedComponent.compute_periodicity()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.compute_periodicity)


                        * [`ConnectedComponent.compute_periodicity_all_simple_paths_algorithm()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.compute_periodicity_all_simple_paths_algorithm)


                        * [`ConnectedComponent.compute_periodicity_cycle_basis()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.compute_periodicity_cycle_basis)


                        * [`ConnectedComponent.coordination_sequence()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.coordination_sequence)


                        * [`ConnectedComponent.description()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.description)


                        * [`ConnectedComponent.elastic_centered_graph()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.elastic_centered_graph)


                        * [`ConnectedComponent.from_dict()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.from_dict)


                        * [`ConnectedComponent.from_graph()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.from_graph)


                        * [`ConnectedComponent.graph`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.graph)


                        * [`ConnectedComponent.is_0d`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.is_0d)


                        * [`ConnectedComponent.is_1d`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.is_1d)


                        * [`ConnectedComponent.is_2d`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.is_2d)


                        * [`ConnectedComponent.is_3d`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.is_3d)


                        * [`ConnectedComponent.is_periodic`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.is_periodic)


                        * [`ConnectedComponent.make_supergraph()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.make_supergraph)


                        * [`ConnectedComponent.periodicity`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.periodicity)


                        * [`ConnectedComponent.periodicity_vectors`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.periodicity_vectors)


                        * [`ConnectedComponent.show_graph()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.ConnectedComponent.show_graph)


                    * [`draw_network()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.draw_network)


                    * [`make_supergraph()`](pymatgen.analysis.chemenv.connectivity.connected_components.md#pymatgen.analysis.chemenv.connectivity.connected_components.make_supergraph)


                * [pymatgen.analysis.chemenv.connectivity.connectivity_finder module](pymatgen.analysis.chemenv.connectivity.connectivity_finder.md)


                    * [`ConnectivityFinder`](pymatgen.analysis.chemenv.connectivity.connectivity_finder.md#pymatgen.analysis.chemenv.connectivity.connectivity_finder.ConnectivityFinder)


                        * [`ConnectivityFinder.get_structure_connectivity()`](pymatgen.analysis.chemenv.connectivity.connectivity_finder.md#pymatgen.analysis.chemenv.connectivity.connectivity_finder.ConnectivityFinder.get_structure_connectivity)


                        * [`ConnectivityFinder.setup_parameters()`](pymatgen.analysis.chemenv.connectivity.connectivity_finder.md#pymatgen.analysis.chemenv.connectivity.connectivity_finder.ConnectivityFinder.setup_parameters)


                * [pymatgen.analysis.chemenv.connectivity.environment_nodes module](pymatgen.analysis.chemenv.connectivity.environment_nodes.md)


                    * [`AbstractEnvironmentNode`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode)


                        * [`AbstractEnvironmentNode.ATOM`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.ATOM)


                        * [`AbstractEnvironmentNode.CE_NNBCES_NBCES_LIGANDS`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.CE_NNBCES_NBCES_LIGANDS)


                        * [`AbstractEnvironmentNode.COORDINATION_ENVIRONMENT`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.COORDINATION_ENVIRONMENT)


                        * [`AbstractEnvironmentNode.DEFAULT_EXTENSIONS`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.DEFAULT_EXTENSIONS)


                        * [`AbstractEnvironmentNode.LIGANDS_ARRANGEMENT`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.LIGANDS_ARRANGEMENT)


                        * [`AbstractEnvironmentNode.NEIGHBORING_CES`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NEIGHBORING_CES)


                        * [`AbstractEnvironmentNode.NEIGHBORING_COORDINATION_ENVIRONMENTS`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NEIGHBORING_COORDINATION_ENVIRONMENTS)


                        * [`AbstractEnvironmentNode.NEIGHBORS_LIGANDS_ARRANGEMENT`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NEIGHBORS_LIGANDS_ARRANGEMENT)


                        * [`AbstractEnvironmentNode.NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_CE`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_CE)


                        * [`AbstractEnvironmentNode.NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_COORDINATION_ENVIRONMENT`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NUMBER_OF_LIGANDS_FOR_EACH_NEIGHBORING_COORDINATION_ENVIRONMENT)


                        * [`AbstractEnvironmentNode.NUMBER_OF_NEIGHBORING_CES`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NUMBER_OF_NEIGHBORING_CES)


                        * [`AbstractEnvironmentNode.NUMBER_OF_NEIGHBORING_COORDINATION_ENVIRONMENTS`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.NUMBER_OF_NEIGHBORING_COORDINATION_ENVIRONMENTS)


                        * [`AbstractEnvironmentNode.atom_symbol`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.atom_symbol)


                        * [`AbstractEnvironmentNode.ce`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.ce)


                        * [`AbstractEnvironmentNode.ce_symbol`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.ce_symbol)


                        * [`AbstractEnvironmentNode.coordination_environment`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.coordination_environment)


                        * [`AbstractEnvironmentNode.everything_equal()`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.everything_equal)


                        * [`AbstractEnvironmentNode.isite`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.isite)


                        * [`AbstractEnvironmentNode.mp_symbol`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.AbstractEnvironmentNode.mp_symbol)


                    * [`EnvironmentNode`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.EnvironmentNode)


                        * [`EnvironmentNode.coordination_environment`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.EnvironmentNode.coordination_environment)


                        * [`EnvironmentNode.everything_equal()`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.EnvironmentNode.everything_equal)


                    * [`get_environment_node()`](pymatgen.analysis.chemenv.connectivity.environment_nodes.md#pymatgen.analysis.chemenv.connectivity.environment_nodes.get_environment_node)


                * [pymatgen.analysis.chemenv.connectivity.structure_connectivity module](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md)


                    * [`StructureConnectivity`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity)


                        * [`StructureConnectivity.add_bonds()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.add_bonds)


                        * [`StructureConnectivity.add_sites()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.add_sites)


                        * [`StructureConnectivity.as_dict()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.as_dict)


                        * [`StructureConnectivity.environment_subgraph()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.environment_subgraph)


                        * [`StructureConnectivity.from_dict()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.from_dict)


                        * [`StructureConnectivity.get_connected_components()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.get_connected_components)


                        * [`StructureConnectivity.print_links()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.print_links)


                        * [`StructureConnectivity.setup_atom_environment_subgraph()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.setup_atom_environment_subgraph)


                        * [`StructureConnectivity.setup_atom_environments_subgraph()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.setup_atom_environments_subgraph)


                        * [`StructureConnectivity.setup_connectivity_description()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.setup_connectivity_description)


                        * [`StructureConnectivity.setup_environment_subgraph()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.setup_environment_subgraph)


                        * [`StructureConnectivity.setup_environments_subgraph()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.StructureConnectivity.setup_environments_subgraph)


                    * [`get_delta_image()`](pymatgen.analysis.chemenv.connectivity.structure_connectivity.md#pymatgen.analysis.chemenv.connectivity.structure_connectivity.get_delta_image)


        * [pymatgen.analysis.chemenv.coordination_environments package](pymatgen.analysis.chemenv.coordination_environments.md)


            * [Subpackages](pymatgen.analysis.chemenv.coordination_environments.md#subpackages)


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


        * [pymatgen.analysis.chemenv.utils package](pymatgen.analysis.chemenv.utils.md)




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


* [pymatgen.analysis.diffraction package](pymatgen.analysis.diffraction.md)




        * [pymatgen.analysis.diffraction.core module](pymatgen.analysis.diffraction.core.md)


            * [`AbstractDiffractionPatternCalculator`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator)


                * [`AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL)


                * [`AbstractDiffractionPatternCalculator.TWO_THETA_TOL`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.TWO_THETA_TOL)


                * [`AbstractDiffractionPatternCalculator.get_pattern()`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.get_pattern)


                * [`AbstractDiffractionPatternCalculator.get_plot()`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.get_plot)


                * [`AbstractDiffractionPatternCalculator.plot_structures()`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.plot_structures)


                * [`AbstractDiffractionPatternCalculator.show_plot()`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.AbstractDiffractionPatternCalculator.show_plot)


            * [`DiffractionPattern`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.DiffractionPattern)


                * [`DiffractionPattern.XLABEL`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.DiffractionPattern.XLABEL)


                * [`DiffractionPattern.YLABEL`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.DiffractionPattern.YLABEL)


            * [`get_unique_families()`](pymatgen.analysis.diffraction.core.md#pymatgen.analysis.diffraction.core.get_unique_families)


        * [pymatgen.analysis.diffraction.neutron module](pymatgen.analysis.diffraction.neutron.md)


            * [`NDCalculator`](pymatgen.analysis.diffraction.neutron.md#pymatgen.analysis.diffraction.neutron.NDCalculator)


                * [`NDCalculator.get_pattern()`](pymatgen.analysis.diffraction.neutron.md#pymatgen.analysis.diffraction.neutron.NDCalculator.get_pattern)


        * [pymatgen.analysis.diffraction.tem module](pymatgen.analysis.diffraction.tem.md)


            * [`TEMCalculator`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator)


                * [`TEMCalculator.bragg_angles()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.bragg_angles)


                * [`TEMCalculator.cell_intensity()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.cell_intensity)


                * [`TEMCalculator.cell_scattering_factors()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.cell_scattering_factors)


                * [`TEMCalculator.electron_scattering_factors()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.electron_scattering_factors)


                * [`TEMCalculator.generate_points()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.generate_points)


                * [`TEMCalculator.get_first_point()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_first_point)


                * [`TEMCalculator.get_interplanar_angle()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_interplanar_angle)


                * [`TEMCalculator.get_interplanar_spacings()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_interplanar_spacings)


                * [`TEMCalculator.get_pattern()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_pattern)


                * [`TEMCalculator.get_plot_2d()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_plot_2d)


                * [`TEMCalculator.get_plot_2d_concise()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_plot_2d_concise)


                * [`TEMCalculator.get_plot_coeffs()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_plot_coeffs)


                * [`TEMCalculator.get_positions()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_positions)


                * [`TEMCalculator.get_s2()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.get_s2)


                * [`TEMCalculator.is_parallel()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.is_parallel)


                * [`TEMCalculator.normalized_cell_intensity()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.normalized_cell_intensity)


                * [`TEMCalculator.tem_dots()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.tem_dots)


                * [`TEMCalculator.wavelength_rel()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.wavelength_rel)


                * [`TEMCalculator.x_ray_factors()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.x_ray_factors)


                * [`TEMCalculator.zone_axis_filter()`](pymatgen.analysis.diffraction.tem.md#pymatgen.analysis.diffraction.tem.TEMCalculator.zone_axis_filter)


        * [pymatgen.analysis.diffraction.xrd module](pymatgen.analysis.diffraction.xrd.md)


            * [`XRDCalculator`](pymatgen.analysis.diffraction.xrd.md#pymatgen.analysis.diffraction.xrd.XRDCalculator)


                * [`XRDCalculator.AVAILABLE_RADIATION`](pymatgen.analysis.diffraction.xrd.md#pymatgen.analysis.diffraction.xrd.XRDCalculator.AVAILABLE_RADIATION)


                * [`XRDCalculator.get_pattern()`](pymatgen.analysis.diffraction.xrd.md#pymatgen.analysis.diffraction.xrd.XRDCalculator.get_pattern)


* [pymatgen.analysis.elasticity package](pymatgen.analysis.elasticity.md)




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


* [pymatgen.analysis.ferroelectricity package](pymatgen.analysis.ferroelectricity.md)




        * [pymatgen.analysis.ferroelectricity.polarization module](pymatgen.analysis.ferroelectricity.polarization.md)


            * [`EnergyTrend`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.EnergyTrend)


                * [`EnergyTrend.endpoints_minima()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.EnergyTrend.endpoints_minima)


                * [`EnergyTrend.max_spline_jump()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.EnergyTrend.max_spline_jump)


                * [`EnergyTrend.smoothness()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.EnergyTrend.smoothness)


                * [`EnergyTrend.spline()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.EnergyTrend.spline)


            * [`Polarization`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.Polarization)


                * [`Polarization.from_outcars_and_structures()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.from_outcars_and_structures)


                * [`Polarization.get_lattice_quanta()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.get_lattice_quanta)


                * [`Polarization.get_pelecs_and_pions()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.get_pelecs_and_pions)


                * [`Polarization.get_polarization_change()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.get_polarization_change)


                * [`Polarization.get_polarization_change_norm()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.get_polarization_change_norm)


                * [`Polarization.get_same_branch_polarization_data()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.get_same_branch_polarization_data)


                * [`Polarization.max_spline_jumps()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.max_spline_jumps)


                * [`Polarization.same_branch_splines()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.same_branch_splines)


                * [`Polarization.smoothness()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.Polarization.smoothness)


            * [`PolarizationLattice`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.PolarizationLattice)


                * [`PolarizationLattice.get_nearest_site()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.PolarizationLattice.get_nearest_site)


            * [`calc_ionic()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.calc_ionic)


            * [`get_total_ionic_dipole()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.get_total_ionic_dipole)


            * [`zval_dict_from_potcar()`](pymatgen.analysis.ferroelectricity.polarization.md#pymatgen.analysis.ferroelectricity.polarization.zval_dict_from_potcar)


* [pymatgen.analysis.gb package](pymatgen.analysis.gb.md)




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


* [pymatgen.analysis.interfaces package](pymatgen.analysis.interfaces.md)




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


* [pymatgen.analysis.magnetism package](pymatgen.analysis.magnetism.md)




        * [pymatgen.analysis.magnetism.analyzer module](pymatgen.analysis.magnetism.analyzer.md)


            * [`CollinearMagneticStructureAnalyzer`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer)


                * [`CollinearMagneticStructureAnalyzer.get_exchange_group_info()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_exchange_group_info)


                * [`CollinearMagneticStructureAnalyzer.get_ferromagnetic_structure()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_ferromagnetic_structure)


                * [`CollinearMagneticStructureAnalyzer.get_nonmagnetic_structure()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_nonmagnetic_structure)


                * [`CollinearMagneticStructureAnalyzer.get_structure_with_only_magnetic_atoms()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_structure_with_only_magnetic_atoms)


                * [`CollinearMagneticStructureAnalyzer.get_structure_with_spin()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.get_structure_with_spin)


                * [`CollinearMagneticStructureAnalyzer.is_magnetic`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.is_magnetic)


                * [`CollinearMagneticStructureAnalyzer.magmoms`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.magmoms)


                * [`CollinearMagneticStructureAnalyzer.magnetic_species_and_magmoms`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.magnetic_species_and_magmoms)


                * [`CollinearMagneticStructureAnalyzer.matches_ordering()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.matches_ordering)


                * [`CollinearMagneticStructureAnalyzer.number_of_magnetic_sites`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.number_of_magnetic_sites)


                * [`CollinearMagneticStructureAnalyzer.number_of_unique_magnetic_sites()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.number_of_unique_magnetic_sites)


                * [`CollinearMagneticStructureAnalyzer.ordering`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.ordering)


                * [`CollinearMagneticStructureAnalyzer.types_of_magnetic_specie`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.types_of_magnetic_specie)


                * [`CollinearMagneticStructureAnalyzer.types_of_magnetic_species`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.CollinearMagneticStructureAnalyzer.types_of_magnetic_species)


            * [`MagneticDeformation`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation)


                * [`MagneticDeformation.deformation`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation.deformation)


                * [`MagneticDeformation.type`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.MagneticDeformation.type)


            * [`MagneticStructureEnumerator`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator)


                * [`MagneticStructureEnumerator.available_strategies`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.MagneticStructureEnumerator.available_strategies)


            * [`Ordering`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering)


                * [`Ordering.AFM`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering.AFM)


                * [`Ordering.FM`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering.FM)


                * [`Ordering.FiM`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering.FiM)


                * [`Ordering.NM`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering.NM)


                * [`Ordering.Unknown`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.Ordering.Unknown)


            * [`OverwriteMagmomMode`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode)


                * [`OverwriteMagmomMode.none`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.none)


                * [`OverwriteMagmomMode.normalize`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.normalize)


                * [`OverwriteMagmomMode.replace_all`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.replace_all)


                * [`OverwriteMagmomMode.respect_sign`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.respect_sign)


                * [`OverwriteMagmomMode.respect_zero`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.OverwriteMagmomMode.respect_zero)


            * [`magnetic_deformation()`](pymatgen.analysis.magnetism.analyzer.md#pymatgen.analysis.magnetism.analyzer.magnetic_deformation)


        * [pymatgen.analysis.magnetism.heisenberg module](pymatgen.analysis.magnetism.heisenberg.md)


            * [`HeisenbergMapper`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper)


                * [`HeisenbergMapper.estimate_exchange()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.estimate_exchange)


                * [`HeisenbergMapper.get_exchange()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_exchange)


                * [`HeisenbergMapper.get_heisenberg_model()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_heisenberg_model)


                * [`HeisenbergMapper.get_interaction_graph()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_interaction_graph)


                * [`HeisenbergMapper.get_low_energy_orderings()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_low_energy_orderings)


                * [`HeisenbergMapper.get_mft_temperature()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergMapper.get_mft_temperature)


            * [`HeisenbergModel`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel)


                * [`HeisenbergModel.as_dict()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel.as_dict)


                * [`HeisenbergModel.from_dict()`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergModel.from_dict)


            * [`HeisenbergScreener`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener)


                * [`HeisenbergScreener.screened_structures`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener.screened_structures)


                * [`HeisenbergScreener.screened_energies`](pymatgen.analysis.magnetism.heisenberg.md#pymatgen.analysis.magnetism.heisenberg.HeisenbergScreener.screened_energies)


        * [pymatgen.analysis.magnetism.jahnteller module](pymatgen.analysis.magnetism.jahnteller.md)


            * [`JahnTellerAnalyzer`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer)


                * [`JahnTellerAnalyzer.get_analysis()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_analysis)


                * [`JahnTellerAnalyzer.get_analysis_and_structure()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_analysis_and_structure)


                * [`JahnTellerAnalyzer.get_magnitude_of_effect_from_species()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_magnitude_of_effect_from_species)


                * [`JahnTellerAnalyzer.get_magnitude_of_effect_from_spin_config()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.get_magnitude_of_effect_from_spin_config)


                * [`JahnTellerAnalyzer.is_jahn_teller_active()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.is_jahn_teller_active)


                * [`JahnTellerAnalyzer.mu_so()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.mu_so)


                * [`JahnTellerAnalyzer.tag_structure()`](pymatgen.analysis.magnetism.jahnteller.md#pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer.tag_structure)


* [pymatgen.analysis.solar package](pymatgen.analysis.solar.md)




        * [pymatgen.analysis.solar.slme module](pymatgen.analysis.solar.slme.md)


            * [`absorption_coefficient()`](pymatgen.analysis.solar.slme.md#pymatgen.analysis.solar.slme.absorption_coefficient)


            * [`get_dir_indir_gap()`](pymatgen.analysis.solar.slme.md#pymatgen.analysis.solar.slme.get_dir_indir_gap)


            * [`matrix_eigvals()`](pymatgen.analysis.solar.slme.md#pymatgen.analysis.solar.slme.matrix_eigvals)


            * [`optics()`](pymatgen.analysis.solar.slme.md#pymatgen.analysis.solar.slme.optics)


            * [`parse_dielectric_data()`](pymatgen.analysis.solar.slme.md#pymatgen.analysis.solar.slme.parse_dielectric_data)


            * [`slme()`](pymatgen.analysis.solar.slme.md#pymatgen.analysis.solar.slme.slme)


            * [`to_matrix()`](pymatgen.analysis.solar.slme.md#pymatgen.analysis.solar.slme.to_matrix)


* [pymatgen.analysis.structure_prediction package](pymatgen.analysis.structure_prediction.md)




        * [pymatgen.analysis.structure_prediction.dopant_predictor module](pymatgen.analysis.structure_prediction.dopant_predictor.md)


            * [`get_dopants_from_shannon_radii()`](pymatgen.analysis.structure_prediction.dopant_predictor.md#pymatgen.analysis.structure_prediction.dopant_predictor.get_dopants_from_shannon_radii)


            * [`get_dopants_from_substitution_probabilities()`](pymatgen.analysis.structure_prediction.dopant_predictor.md#pymatgen.analysis.structure_prediction.dopant_predictor.get_dopants_from_substitution_probabilities)


        * [pymatgen.analysis.structure_prediction.substitution_probability module](pymatgen.analysis.structure_prediction.substitution_probability.md)


            * [`SubstitutionPredictor`](pymatgen.analysis.structure_prediction.substitution_probability.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionPredictor)


                * [`SubstitutionPredictor.composition_prediction()`](pymatgen.analysis.structure_prediction.substitution_probability.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionPredictor.composition_prediction)


                * [`SubstitutionPredictor.list_prediction()`](pymatgen.analysis.structure_prediction.substitution_probability.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionPredictor.list_prediction)


            * [`SubstitutionProbability`](pymatgen.analysis.structure_prediction.substitution_probability.md#pymatgen.analysis.structure_prediction.substitution_probability.SubstitutionProbability)


        * [pymatgen.analysis.structure_prediction.substitutor module](pymatgen.analysis.structure_prediction.substitutor.md)


            * [`Substitutor`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor)


                * [`Substitutor.as_dict()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.as_dict)


                * [`Substitutor.from_dict()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.from_dict)


                * [`Substitutor.get_allowed_species()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.get_allowed_species)


                * [`Substitutor.pred_from_comp()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.pred_from_comp)


                * [`Substitutor.pred_from_list()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.pred_from_list)


                * [`Substitutor.pred_from_structures()`](pymatgen.analysis.structure_prediction.substitutor.md#pymatgen.analysis.structure_prediction.substitutor.Substitutor.pred_from_structures)


        * [pymatgen.analysis.structure_prediction.volume_predictor module](pymatgen.analysis.structure_prediction.volume_predictor.md)


            * [`DLSVolumePredictor`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.DLSVolumePredictor)


                * [`DLSVolumePredictor.get_predicted_structure()`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.DLSVolumePredictor.get_predicted_structure)


                * [`DLSVolumePredictor.predict()`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.DLSVolumePredictor.predict)


            * [`RLSVolumePredictor`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.RLSVolumePredictor)


                * [`RLSVolumePredictor.get_predicted_structure()`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.RLSVolumePredictor.get_predicted_structure)


                * [`RLSVolumePredictor.predict()`](pymatgen.analysis.structure_prediction.volume_predictor.md#pymatgen.analysis.structure_prediction.volume_predictor.RLSVolumePredictor.predict)


* [pymatgen.analysis.topological package](pymatgen.analysis.topological.md)




        * [pymatgen.analysis.topological.spillage module](pymatgen.analysis.topological.spillage.md)


            * [`SOCSpillage`](pymatgen.analysis.topological.spillage.md#pymatgen.analysis.topological.spillage.SOCSpillage)


                * [`SOCSpillage.isclose()`](pymatgen.analysis.topological.spillage.md#pymatgen.analysis.topological.spillage.SOCSpillage.isclose)


                * [`SOCSpillage.orth()`](pymatgen.analysis.topological.spillage.md#pymatgen.analysis.topological.spillage.SOCSpillage.orth)


                * [`SOCSpillage.overlap_so_spinpol()`](pymatgen.analysis.topological.spillage.md#pymatgen.analysis.topological.spillage.SOCSpillage.overlap_so_spinpol)


* [pymatgen.analysis.xas package](pymatgen.analysis.xas.md)




        * [pymatgen.analysis.xas.spectrum module](pymatgen.analysis.xas.spectrum.md)


            * [`XAS`](pymatgen.analysis.xas.spectrum.md#pymatgen.analysis.xas.spectrum.XAS)


                * [`XAS.XLABEL`](pymatgen.analysis.xas.spectrum.md#pymatgen.analysis.xas.spectrum.XAS.XLABEL)


                * [`XAS.YLABEL`](pymatgen.analysis.xas.spectrum.md#pymatgen.analysis.xas.spectrum.XAS.YLABEL)


                * [`XAS.stitch()`](pymatgen.analysis.xas.spectrum.md#pymatgen.analysis.xas.spectrum.XAS.stitch)


            * [`site_weighted_spectrum()`](pymatgen.analysis.xas.spectrum.md#pymatgen.analysis.xas.spectrum.site_weighted_spectrum)




* [pymatgen.analysis.adsorption module](pymatgen.analysis.adsorption.md)


    * [`AdsorbateSiteFinder`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder)


        * [`AdsorbateSiteFinder.add_adsorbate()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.add_adsorbate)


        * [`AdsorbateSiteFinder.adsorb_both_surfaces()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.adsorb_both_surfaces)


        * [`AdsorbateSiteFinder.assign_selective_dynamics()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.assign_selective_dynamics)


        * [`AdsorbateSiteFinder.assign_site_properties()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.assign_site_properties)


        * [`AdsorbateSiteFinder.ensemble_center()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.ensemble_center)


        * [`AdsorbateSiteFinder.find_adsorption_sites()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.find_adsorption_sites)


        * [`AdsorbateSiteFinder.find_surface_sites_by_height()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.find_surface_sites_by_height)


        * [`AdsorbateSiteFinder.from_bulk_and_miller()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.from_bulk_and_miller)


        * [`AdsorbateSiteFinder.generate_adsorption_structures()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.generate_adsorption_structures)


        * [`AdsorbateSiteFinder.generate_substitution_structures()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.generate_substitution_structures)


        * [`AdsorbateSiteFinder.get_extended_surface_mesh()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.get_extended_surface_mesh)


        * [`AdsorbateSiteFinder.near_reduce()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.near_reduce)


        * [`AdsorbateSiteFinder.subsurface_sites()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.subsurface_sites)


        * [`AdsorbateSiteFinder.surface_sites`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.surface_sites)


        * [`AdsorbateSiteFinder.symm_reduce()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.AdsorbateSiteFinder.symm_reduce)


    * [`get_mi_vec()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.get_mi_vec)


    * [`get_rot()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.get_rot)


    * [`plot_slab()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.plot_slab)


    * [`put_coord_inside()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.put_coord_inside)


    * [`reorient_z()`](pymatgen.analysis.adsorption.md#pymatgen.analysis.adsorption.reorient_z)


* [pymatgen.analysis.bond_dissociation module](pymatgen.analysis.bond_dissociation.md)


    * [`BondDissociationEnergies`](pymatgen.analysis.bond_dissociation.md#pymatgen.analysis.bond_dissociation.BondDissociationEnergies)


        * [`BondDissociationEnergies.build_new_entry()`](pymatgen.analysis.bond_dissociation.md#pymatgen.analysis.bond_dissociation.BondDissociationEnergies.build_new_entry)


        * [`BondDissociationEnergies.filter_fragment_entries()`](pymatgen.analysis.bond_dissociation.md#pymatgen.analysis.bond_dissociation.BondDissociationEnergies.filter_fragment_entries)


        * [`BondDissociationEnergies.fragment_and_process()`](pymatgen.analysis.bond_dissociation.md#pymatgen.analysis.bond_dissociation.BondDissociationEnergies.fragment_and_process)


        * [`BondDissociationEnergies.search_fragment_entries()`](pymatgen.analysis.bond_dissociation.md#pymatgen.analysis.bond_dissociation.BondDissociationEnergies.search_fragment_entries)


* [pymatgen.analysis.bond_valence module](pymatgen.analysis.bond_valence.md)


    * [`BVAnalyzer`](pymatgen.analysis.bond_valence.md#pymatgen.analysis.bond_valence.BVAnalyzer)


        * [`BVAnalyzer.CHARGE_NEUTRALITY_TOLERANCE`](pymatgen.analysis.bond_valence.md#pymatgen.analysis.bond_valence.BVAnalyzer.CHARGE_NEUTRALITY_TOLERANCE)


        * [`BVAnalyzer.get_oxi_state_decorated_structure()`](pymatgen.analysis.bond_valence.md#pymatgen.analysis.bond_valence.BVAnalyzer.get_oxi_state_decorated_structure)


        * [`BVAnalyzer.get_valences()`](pymatgen.analysis.bond_valence.md#pymatgen.analysis.bond_valence.BVAnalyzer.get_valences)


    * [`add_oxidation_state_by_site_fraction()`](pymatgen.analysis.bond_valence.md#pymatgen.analysis.bond_valence.add_oxidation_state_by_site_fraction)


    * [`calculate_bv_sum()`](pymatgen.analysis.bond_valence.md#pymatgen.analysis.bond_valence.calculate_bv_sum)


    * [`calculate_bv_sum_unordered()`](pymatgen.analysis.bond_valence.md#pymatgen.analysis.bond_valence.calculate_bv_sum_unordered)


    * [`get_z_ordered_elmap()`](pymatgen.analysis.bond_valence.md#pymatgen.analysis.bond_valence.get_z_ordered_elmap)


* [pymatgen.analysis.chempot_diagram module](pymatgen.analysis.chempot_diagram.md)


    * [`ChemicalPotentialDiagram`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram)


        * [`ChemicalPotentialDiagram.border_hyperplanes`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram.border_hyperplanes)


        * [`ChemicalPotentialDiagram.chemical_system`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram.chemical_system)


        * [`ChemicalPotentialDiagram.domains`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram.domains)


        * [`ChemicalPotentialDiagram.el_refs`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram.el_refs)


        * [`ChemicalPotentialDiagram.entry_dict`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram.entry_dict)


        * [`ChemicalPotentialDiagram.get_plot()`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram.get_plot)


        * [`ChemicalPotentialDiagram.hyperplane_entries`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram.hyperplane_entries)


        * [`ChemicalPotentialDiagram.hyperplanes`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram.hyperplanes)


        * [`ChemicalPotentialDiagram.lims`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.ChemicalPotentialDiagram.lims)


    * [`get_2d_orthonormal_vector()`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.get_2d_orthonormal_vector)


    * [`get_centroid_2d()`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.get_centroid_2d)


    * [`simple_pca()`](pymatgen.analysis.chempot_diagram.md#pymatgen.analysis.chempot_diagram.simple_pca)


* [pymatgen.analysis.cost module](pymatgen.analysis.cost.md)


    * [`CostAnalyzer`](pymatgen.analysis.cost.md#pymatgen.analysis.cost.CostAnalyzer)


        * [`CostAnalyzer.get_cost_per_kg()`](pymatgen.analysis.cost.md#pymatgen.analysis.cost.CostAnalyzer.get_cost_per_kg)


        * [`CostAnalyzer.get_cost_per_mol()`](pymatgen.analysis.cost.md#pymatgen.analysis.cost.CostAnalyzer.get_cost_per_mol)


        * [`CostAnalyzer.get_lowest_decomposition()`](pymatgen.analysis.cost.md#pymatgen.analysis.cost.CostAnalyzer.get_lowest_decomposition)


    * [`CostDB`](pymatgen.analysis.cost.md#pymatgen.analysis.cost.CostDB)


        * [`CostDB.get_entries()`](pymatgen.analysis.cost.md#pymatgen.analysis.cost.CostDB.get_entries)


    * [`CostDBCSV`](pymatgen.analysis.cost.md#pymatgen.analysis.cost.CostDBCSV)


        * [`CostDBCSV.get_entries()`](pymatgen.analysis.cost.md#pymatgen.analysis.cost.CostDBCSV.get_entries)


    * [`CostEntry`](pymatgen.analysis.cost.md#pymatgen.analysis.cost.CostEntry)


* [pymatgen.analysis.dimensionality module](pymatgen.analysis.dimensionality.md)


    * [`calculate_dimensionality_of_site()`](pymatgen.analysis.dimensionality.md#pymatgen.analysis.dimensionality.calculate_dimensionality_of_site)


    * [`find_clusters()`](pymatgen.analysis.dimensionality.md#pymatgen.analysis.dimensionality.find_clusters)


    * [`find_connected_atoms()`](pymatgen.analysis.dimensionality.md#pymatgen.analysis.dimensionality.find_connected_atoms)


    * [`get_dimensionality_cheon()`](pymatgen.analysis.dimensionality.md#pymatgen.analysis.dimensionality.get_dimensionality_cheon)


    * [`get_dimensionality_gorai()`](pymatgen.analysis.dimensionality.md#pymatgen.analysis.dimensionality.get_dimensionality_gorai)


    * [`get_dimensionality_larsen()`](pymatgen.analysis.dimensionality.md#pymatgen.analysis.dimensionality.get_dimensionality_larsen)


    * [`get_structure_components()`](pymatgen.analysis.dimensionality.md#pymatgen.analysis.dimensionality.get_structure_components)


    * [`zero_d_graph_to_molecule_graph()`](pymatgen.analysis.dimensionality.md#pymatgen.analysis.dimensionality.zero_d_graph_to_molecule_graph)


* [pymatgen.analysis.disorder module](pymatgen.analysis.disorder.md)


    * [`get_warren_cowley_parameters()`](pymatgen.analysis.disorder.md#pymatgen.analysis.disorder.get_warren_cowley_parameters)


* [pymatgen.analysis.energy_models module](pymatgen.analysis.energy_models.md)


    * [`EnergyModel`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.EnergyModel)


        * [`EnergyModel.from_dict()`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.EnergyModel.from_dict)


        * [`EnergyModel.get_energy()`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.EnergyModel.get_energy)


    * [`EwaldElectrostaticModel`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.EwaldElectrostaticModel)


        * [`EwaldElectrostaticModel.as_dict()`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.EwaldElectrostaticModel.as_dict)


        * [`EwaldElectrostaticModel.get_energy()`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.EwaldElectrostaticModel.get_energy)


    * [`IsingModel`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.IsingModel)


        * [`IsingModel.as_dict()`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.IsingModel.as_dict)


        * [`IsingModel.get_energy()`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.IsingModel.get_energy)


    * [`NsitesModel`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.NsitesModel)


        * [`NsitesModel.as_dict()`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.NsitesModel.as_dict)


        * [`NsitesModel.get_energy()`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.NsitesModel.get_energy)


    * [`SymmetryModel`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.SymmetryModel)


        * [`SymmetryModel.as_dict()`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.SymmetryModel.as_dict)


        * [`SymmetryModel.get_energy()`](pymatgen.analysis.energy_models.md#pymatgen.analysis.energy_models.SymmetryModel.get_energy)


* [pymatgen.analysis.eos module](pymatgen.analysis.eos.md)


    * [`Birch`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.Birch)


    * [`BirchMurnaghan`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.BirchMurnaghan)


    * [`DeltaFactor`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.DeltaFactor)


        * [`DeltaFactor.fit()`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.DeltaFactor.fit)


    * [`EOS`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOS)


        * [`EOS.MODELS`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOS.MODELS)


        * [`EOS.fit()`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOS.fit)


    * [`EOSBase`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase)


        * [`EOSBase.b0`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase.b0)


        * [`EOSBase.b0_GPa`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase.b0_GPa)


        * [`EOSBase.b1`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase.b1)


        * [`EOSBase.e0`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase.e0)


        * [`EOSBase.fit()`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase.fit)


        * [`EOSBase.func()`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase.func)


        * [`EOSBase.plot()`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase.plot)


        * [`EOSBase.plot_ax()`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase.plot_ax)


        * [`EOSBase.results`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase.results)


        * [`EOSBase.v0`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSBase.v0)


    * [`EOSError`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.EOSError)


    * [`Murnaghan`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.Murnaghan)


    * [`NumericalEOS`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.NumericalEOS)


        * [`NumericalEOS.fit()`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.NumericalEOS.fit)


    * [`PolynomialEOS`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.PolynomialEOS)


        * [`PolynomialEOS.fit()`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.PolynomialEOS.fit)


    * [`PourierTarantola`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.PourierTarantola)


    * [`Vinet`](pymatgen.analysis.eos.md#pymatgen.analysis.eos.Vinet)


* [pymatgen.analysis.ewald module](pymatgen.analysis.ewald.md)


    * [`EwaldMinimizer`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer)


        * [`EwaldMinimizer.ALGO_BEST_FIRST`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.ALGO_BEST_FIRST)


        * [`EwaldMinimizer.ALGO_COMPLETE`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.ALGO_COMPLETE)


        * [`EwaldMinimizer.ALGO_FAST`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.ALGO_FAST)


        * [`EwaldMinimizer.ALGO_TIME_LIMIT`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.ALGO_TIME_LIMIT)


        * [`EwaldMinimizer.add_m_list()`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.add_m_list)


        * [`EwaldMinimizer.best_case()`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.best_case)


        * [`EwaldMinimizer.best_m_list`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.best_m_list)


        * [`EwaldMinimizer.get_next_index()`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.get_next_index)


        * [`EwaldMinimizer.minimize_matrix()`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.minimize_matrix)


        * [`EwaldMinimizer.minimized_sum`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.minimized_sum)


        * [`EwaldMinimizer.output_lists`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldMinimizer.output_lists)


    * [`EwaldSummation`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation)


        * [`EwaldSummation.CONV_FACT`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.CONV_FACT)


        * [`EwaldSummation.as_dict()`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.as_dict)


        * [`EwaldSummation.compute_partial_energy()`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.compute_partial_energy)


        * [`EwaldSummation.compute_sub_structure()`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.compute_sub_structure)


        * [`EwaldSummation.eta`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.eta)


        * [`EwaldSummation.forces`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.forces)


        * [`EwaldSummation.from_dict()`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.from_dict)


        * [`EwaldSummation.get_site_energy()`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.get_site_energy)


        * [`EwaldSummation.point_energy`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.point_energy)


        * [`EwaldSummation.point_energy_matrix`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.point_energy_matrix)


        * [`EwaldSummation.real_space_energy`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.real_space_energy)


        * [`EwaldSummation.real_space_energy_matrix`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.real_space_energy_matrix)


        * [`EwaldSummation.reciprocal_space_energy`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.reciprocal_space_energy)


        * [`EwaldSummation.reciprocal_space_energy_matrix`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.reciprocal_space_energy_matrix)


        * [`EwaldSummation.total_energy`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.total_energy)


        * [`EwaldSummation.total_energy_matrix`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.EwaldSummation.total_energy_matrix)


    * [`compute_average_oxidation_state()`](pymatgen.analysis.ewald.md#pymatgen.analysis.ewald.compute_average_oxidation_state)


* [pymatgen.analysis.excitation module](pymatgen.analysis.excitation.md)


    * [`ExcitationSpectrum`](pymatgen.analysis.excitation.md#pymatgen.analysis.excitation.ExcitationSpectrum)


        * [`ExcitationSpectrum.XLABEL`](pymatgen.analysis.excitation.md#pymatgen.analysis.excitation.ExcitationSpectrum.XLABEL)


        * [`ExcitationSpectrum.YLABEL`](pymatgen.analysis.excitation.md#pymatgen.analysis.excitation.ExcitationSpectrum.YLABEL)


* [pymatgen.analysis.fragmenter module](pymatgen.analysis.fragmenter.md)


    * [`Fragmenter`](pymatgen.analysis.fragmenter.md#pymatgen.analysis.fragmenter.Fragmenter)


    * [`open_ring()`](pymatgen.analysis.fragmenter.md#pymatgen.analysis.fragmenter.open_ring)


* [pymatgen.analysis.functional_groups module](pymatgen.analysis.functional_groups.md)


    * [`FunctionalGroupExtractor`](pymatgen.analysis.functional_groups.md#pymatgen.analysis.functional_groups.FunctionalGroupExtractor)


        * [`FunctionalGroupExtractor.categorize_functional_groups()`](pymatgen.analysis.functional_groups.md#pymatgen.analysis.functional_groups.FunctionalGroupExtractor.categorize_functional_groups)


        * [`FunctionalGroupExtractor.get_all_functional_groups()`](pymatgen.analysis.functional_groups.md#pymatgen.analysis.functional_groups.FunctionalGroupExtractor.get_all_functional_groups)


        * [`FunctionalGroupExtractor.get_basic_functional_groups()`](pymatgen.analysis.functional_groups.md#pymatgen.analysis.functional_groups.FunctionalGroupExtractor.get_basic_functional_groups)


        * [`FunctionalGroupExtractor.get_heteroatoms()`](pymatgen.analysis.functional_groups.md#pymatgen.analysis.functional_groups.FunctionalGroupExtractor.get_heteroatoms)


        * [`FunctionalGroupExtractor.get_special_carbon()`](pymatgen.analysis.functional_groups.md#pymatgen.analysis.functional_groups.FunctionalGroupExtractor.get_special_carbon)


        * [`FunctionalGroupExtractor.link_marked_atoms()`](pymatgen.analysis.functional_groups.md#pymatgen.analysis.functional_groups.FunctionalGroupExtractor.link_marked_atoms)


* [pymatgen.analysis.graphs module](pymatgen.analysis.graphs.md)


    * [`ConnectedSite`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.ConnectedSite)


        * [`ConnectedSite.dist`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.ConnectedSite.dist)


        * [`ConnectedSite.index`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.ConnectedSite.index)


        * [`ConnectedSite.jimage`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.ConnectedSite.jimage)


        * [`ConnectedSite.site`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.ConnectedSite.site)


        * [`ConnectedSite.weight`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.ConnectedSite.weight)


    * [`MolGraphSplitError`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MolGraphSplitError)


    * [`MoleculeGraph`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph)


        * [`MoleculeGraph.add_edge()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.add_edge)


        * [`MoleculeGraph.alter_edge()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.alter_edge)


        * [`MoleculeGraph.as_dict()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.as_dict)


        * [`MoleculeGraph.break_edge()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.break_edge)


        * [`MoleculeGraph.build_unique_fragments()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.build_unique_fragments)


        * [`MoleculeGraph.diff()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.diff)


        * [`MoleculeGraph.draw_graph_to_file()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.draw_graph_to_file)


        * [`MoleculeGraph.edge_weight_name`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.edge_weight_name)


        * [`MoleculeGraph.edge_weight_unit`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.edge_weight_unit)


        * [`MoleculeGraph.find_rings()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.find_rings)


        * [`MoleculeGraph.from_dict()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.from_dict)


        * [`MoleculeGraph.get_connected_sites()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.get_connected_sites)


        * [`MoleculeGraph.get_coordination_of_site()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.get_coordination_of_site)


        * [`MoleculeGraph.get_disconnected_fragments()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.get_disconnected_fragments)


        * [`MoleculeGraph.insert_node()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.insert_node)


        * [`MoleculeGraph.isomorphic_to()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.isomorphic_to)


        * [`MoleculeGraph.name`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.name)


        * [`MoleculeGraph.remove_nodes()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.remove_nodes)


        * [`MoleculeGraph.replace_group()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.replace_group)


        * [`MoleculeGraph.set_node_attributes()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.set_node_attributes)


        * [`MoleculeGraph.sort()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.sort)


        * [`MoleculeGraph.split_molecule_subgraphs()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.split_molecule_subgraphs)


        * [`MoleculeGraph.substitute_group()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.substitute_group)


        * [`MoleculeGraph.with_edges()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.with_edges)


        * [`MoleculeGraph.with_empty_graph()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.with_empty_graph)


        * [`MoleculeGraph.with_local_env_strategy()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.MoleculeGraph.with_local_env_strategy)


    * [`StructureGraph`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph)


        * [`StructureGraph.add_edge()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.add_edge)


        * [`StructureGraph.alter_edge()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.alter_edge)


        * [`StructureGraph.as_dict()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.as_dict)


        * [`StructureGraph.break_edge()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.break_edge)


        * [`StructureGraph.diff()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.diff)


        * [`StructureGraph.draw_graph_to_file()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.draw_graph_to_file)


        * [`StructureGraph.edge_weight_name`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.edge_weight_name)


        * [`StructureGraph.edge_weight_unit`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.edge_weight_unit)


        * [`StructureGraph.from_dict()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.from_dict)


        * [`StructureGraph.get_connected_sites()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.get_connected_sites)


        * [`StructureGraph.get_coordination_of_site()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.get_coordination_of_site)


        * [`StructureGraph.get_subgraphs_as_molecules()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.get_subgraphs_as_molecules)


        * [`StructureGraph.insert_node()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.insert_node)


        * [`StructureGraph.name`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.name)


        * [`StructureGraph.remove_nodes()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.remove_nodes)


        * [`StructureGraph.set_node_attributes()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.set_node_attributes)


        * [`StructureGraph.sort()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.sort)


        * [`StructureGraph.substitute_group()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.substitute_group)


        * [`StructureGraph.types_and_weights_of_connections`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.types_and_weights_of_connections)


        * [`StructureGraph.types_of_coordination_environments()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.types_of_coordination_environments)


        * [`StructureGraph.weight_statistics`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.weight_statistics)


        * [`StructureGraph.with_edges()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.with_edges)


        * [`StructureGraph.with_empty_graph()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.with_empty_graph)


        * [`StructureGraph.with_local_env_strategy()`](pymatgen.analysis.graphs.md#pymatgen.analysis.graphs.StructureGraph.with_local_env_strategy)


* [pymatgen.analysis.hhi module](pymatgen.analysis.hhi.md)


* [pymatgen.analysis.interface module](pymatgen.analysis.interface.md)


* [pymatgen.analysis.interface_reactions module](pymatgen.analysis.interface_reactions.md)


    * [`GrandPotentialInterfacialReactivity`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.GrandPotentialInterfacialReactivity)


        * [`GrandPotentialInterfacialReactivity.get_no_mixing_energy()`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.GrandPotentialInterfacialReactivity.get_no_mixing_energy)


    * [`InterfacialReactivity`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.InterfacialReactivity)


        * [`InterfacialReactivity.EV_TO_KJ_PER_MOL`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.InterfacialReactivity.EV_TO_KJ_PER_MOL)


        * [`InterfacialReactivity.get_chempot_correction()`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.InterfacialReactivity.get_chempot_correction)


        * [`InterfacialReactivity.get_critical_original_kink_ratio()`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.InterfacialReactivity.get_critical_original_kink_ratio)


        * [`InterfacialReactivity.get_dataframe()`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.InterfacialReactivity.get_dataframe)


        * [`InterfacialReactivity.get_kinks()`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.InterfacialReactivity.get_kinks)


        * [`InterfacialReactivity.labels`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.InterfacialReactivity.labels)


        * [`InterfacialReactivity.minimum`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.InterfacialReactivity.minimum)


        * [`InterfacialReactivity.plot()`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.InterfacialReactivity.plot)


        * [`InterfacialReactivity.products`](pymatgen.analysis.interface_reactions.md#pymatgen.analysis.interface_reactions.InterfacialReactivity.products)


* [pymatgen.analysis.local_env module](pymatgen.analysis.local_env.md)


    * [`BrunnerNN_real`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_real)


        * [`BrunnerNN_real.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_real.get_nn_info)


        * [`BrunnerNN_real.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_real.molecules_allowed)


        * [`BrunnerNN_real.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_real.structures_allowed)


    * [`BrunnerNN_reciprocal`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_reciprocal)


        * [`BrunnerNN_reciprocal.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_reciprocal.get_nn_info)


        * [`BrunnerNN_reciprocal.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_reciprocal.molecules_allowed)


        * [`BrunnerNN_reciprocal.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_reciprocal.structures_allowed)


    * [`BrunnerNN_relative`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_relative)


        * [`BrunnerNN_relative.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_relative.get_nn_info)


        * [`BrunnerNN_relative.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_relative.molecules_allowed)


        * [`BrunnerNN_relative.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.BrunnerNN_relative.structures_allowed)


    * [`CovalentBondNN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CovalentBondNN)


        * [`CovalentBondNN.extend_structure_molecules`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CovalentBondNN.extend_structure_molecules)


        * [`CovalentBondNN.get_bonded_structure()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CovalentBondNN.get_bonded_structure)


        * [`CovalentBondNN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CovalentBondNN.get_nn_info)


        * [`CovalentBondNN.get_nn_shell_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CovalentBondNN.get_nn_shell_info)


        * [`CovalentBondNN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CovalentBondNN.molecules_allowed)


        * [`CovalentBondNN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CovalentBondNN.structures_allowed)


    * [`Critic2NN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.Critic2NN)


        * [`Critic2NN.extend_structure_molecules`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.Critic2NN.extend_structure_molecules)


        * [`Critic2NN.get_bonded_structure()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.Critic2NN.get_bonded_structure)


        * [`Critic2NN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.Critic2NN.get_nn_info)


        * [`Critic2NN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.Critic2NN.molecules_allowed)


        * [`Critic2NN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.Critic2NN.structures_allowed)


    * [`CrystalNN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN)


        * [`CrystalNN.NNData`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.NNData)


            * [`CrystalNN.NNData.all_nninfo`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.NNData.all_nninfo)


            * [`CrystalNN.NNData.cn_nninfo`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.NNData.cn_nninfo)


            * [`CrystalNN.NNData.cn_weights`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.NNData.cn_weights)


        * [`CrystalNN.get_cn()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.get_cn)


        * [`CrystalNN.get_cn_dict()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.get_cn_dict)


        * [`CrystalNN.get_nn_data()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.get_nn_data)


        * [`CrystalNN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.get_nn_info)


        * [`CrystalNN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.molecules_allowed)


        * [`CrystalNN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.structures_allowed)


        * [`CrystalNN.transform_to_length()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CrystalNN.transform_to_length)


    * [`CutOffDictNN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CutOffDictNN)


        * [`CutOffDictNN.extend_structure_molecules`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CutOffDictNN.extend_structure_molecules)


        * [`CutOffDictNN.from_preset()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CutOffDictNN.from_preset)


        * [`CutOffDictNN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CutOffDictNN.get_nn_info)


        * [`CutOffDictNN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CutOffDictNN.molecules_allowed)


        * [`CutOffDictNN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.CutOffDictNN.structures_allowed)


    * [`EconNN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.EconNN)


        * [`EconNN.extend_structure_molecules`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.EconNN.extend_structure_molecules)


        * [`EconNN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.EconNN.get_nn_info)


        * [`EconNN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.EconNN.molecules_allowed)


        * [`EconNN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.EconNN.structures_allowed)


    * [`IsayevNN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.IsayevNN)


        * [`IsayevNN.get_all_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.IsayevNN.get_all_nn_info)


        * [`IsayevNN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.IsayevNN.get_nn_info)


    * [`JmolNN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.JmolNN)


        * [`JmolNN.extend_structure_molecules`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.JmolNN.extend_structure_molecules)


        * [`JmolNN.get_max_bond_distance()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.JmolNN.get_max_bond_distance)


        * [`JmolNN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.JmolNN.get_nn_info)


        * [`JmolNN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.JmolNN.molecules_allowed)


        * [`JmolNN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.JmolNN.structures_allowed)


    * [`LocalStructOrderParams`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.LocalStructOrderParams)


        * [`LocalStructOrderParams.compute_trigonometric_terms()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.LocalStructOrderParams.compute_trigonometric_terms)


        * [`LocalStructOrderParams.get_order_parameters()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.LocalStructOrderParams.get_order_parameters)


        * [`LocalStructOrderParams.get_parameters()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.LocalStructOrderParams.get_parameters)


        * [`LocalStructOrderParams.get_q2()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.LocalStructOrderParams.get_q2)


        * [`LocalStructOrderParams.get_q4()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.LocalStructOrderParams.get_q4)


        * [`LocalStructOrderParams.get_q6()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.LocalStructOrderParams.get_q6)


        * [`LocalStructOrderParams.get_type()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.LocalStructOrderParams.get_type)


        * [`LocalStructOrderParams.last_nneigh`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.LocalStructOrderParams.last_nneigh)


        * [`LocalStructOrderParams.num_ops`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.LocalStructOrderParams.num_ops)


    * [`MinimumDistanceNN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumDistanceNN)


        * [`MinimumDistanceNN.extend_structure_molecules`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumDistanceNN.extend_structure_molecules)


        * [`MinimumDistanceNN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumDistanceNN.get_nn_info)


        * [`MinimumDistanceNN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumDistanceNN.molecules_allowed)


        * [`MinimumDistanceNN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumDistanceNN.structures_allowed)


    * [`MinimumOKeeffeNN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumOKeeffeNN)


        * [`MinimumOKeeffeNN.extend_structure_molecules`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumOKeeffeNN.extend_structure_molecules)


        * [`MinimumOKeeffeNN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumOKeeffeNN.get_nn_info)


        * [`MinimumOKeeffeNN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumOKeeffeNN.molecules_allowed)


        * [`MinimumOKeeffeNN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumOKeeffeNN.structures_allowed)


    * [`MinimumVIRENN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumVIRENN)


        * [`MinimumVIRENN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumVIRENN.get_nn_info)


        * [`MinimumVIRENN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumVIRENN.molecules_allowed)


        * [`MinimumVIRENN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.MinimumVIRENN.structures_allowed)


    * [`NearNeighbors`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors)


        * [`NearNeighbors.extend_structure_molecules`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.extend_structure_molecules)


        * [`NearNeighbors.get_all_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.get_all_nn_info)


        * [`NearNeighbors.get_bonded_structure()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.get_bonded_structure)


        * [`NearNeighbors.get_cn()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.get_cn)


        * [`NearNeighbors.get_cn_dict()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.get_cn_dict)


        * [`NearNeighbors.get_local_order_parameters()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.get_local_order_parameters)


        * [`NearNeighbors.get_nn()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.get_nn)


        * [`NearNeighbors.get_nn_images()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.get_nn_images)


        * [`NearNeighbors.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.get_nn_info)


        * [`NearNeighbors.get_nn_shell_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.get_nn_shell_info)


        * [`NearNeighbors.get_weights_of_nn_sites()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.get_weights_of_nn_sites)


        * [`NearNeighbors.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.molecules_allowed)


        * [`NearNeighbors.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.NearNeighbors.structures_allowed)


    * [`OpenBabelNN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.OpenBabelNN)


        * [`OpenBabelNN.extend_structure_molecules`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.OpenBabelNN.extend_structure_molecules)


        * [`OpenBabelNN.get_bonded_structure()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.OpenBabelNN.get_bonded_structure)


        * [`OpenBabelNN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.OpenBabelNN.get_nn_info)


        * [`OpenBabelNN.get_nn_shell_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.OpenBabelNN.get_nn_shell_info)


        * [`OpenBabelNN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.OpenBabelNN.molecules_allowed)


        * [`OpenBabelNN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.OpenBabelNN.structures_allowed)


    * [`ValenceIonicRadiusEvaluator`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.ValenceIonicRadiusEvaluator)


        * [`ValenceIonicRadiusEvaluator.radii`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.ValenceIonicRadiusEvaluator.radii)


        * [`ValenceIonicRadiusEvaluator.structure`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.ValenceIonicRadiusEvaluator.structure)


        * [`ValenceIonicRadiusEvaluator.valences`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.ValenceIonicRadiusEvaluator.valences)


    * [`VoronoiNN`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.VoronoiNN)


        * [`VoronoiNN.get_all_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.VoronoiNN.get_all_nn_info)


        * [`VoronoiNN.get_all_voronoi_polyhedra()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.VoronoiNN.get_all_voronoi_polyhedra)


        * [`VoronoiNN.get_nn_info()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.VoronoiNN.get_nn_info)


        * [`VoronoiNN.get_voronoi_polyhedra()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.VoronoiNN.get_voronoi_polyhedra)


        * [`VoronoiNN.molecules_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.VoronoiNN.molecules_allowed)


        * [`VoronoiNN.structures_allowed`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.VoronoiNN.structures_allowed)


    * [`get_neighbors_of_site_with_index()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.get_neighbors_of_site_with_index)


    * [`get_okeeffe_distance_prediction()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.get_okeeffe_distance_prediction)


    * [`get_okeeffe_params()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.get_okeeffe_params)


    * [`gramschmidt()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.gramschmidt)


    * [`metal_edge_extender()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.metal_edge_extender)


    * [`oxygen_edge_extender()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.oxygen_edge_extender)


    * [`site_is_of_motif_type()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.site_is_of_motif_type)


    * [`solid_angle()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.solid_angle)


    * [`vol_tetra()`](pymatgen.analysis.local_env.md#pymatgen.analysis.local_env.vol_tetra)


* [pymatgen.analysis.molecule_matcher module](pymatgen.analysis.molecule_matcher.md)


    * [`AbstractMolAtomMapper`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.AbstractMolAtomMapper)


        * [`AbstractMolAtomMapper.from_dict()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.AbstractMolAtomMapper.from_dict)


        * [`AbstractMolAtomMapper.get_molecule_hash()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.AbstractMolAtomMapper.get_molecule_hash)


        * [`AbstractMolAtomMapper.uniform_labels()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.AbstractMolAtomMapper.uniform_labels)


    * [`BruteForceOrderMatcher`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.BruteForceOrderMatcher)


        * [`BruteForceOrderMatcher.fit()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.BruteForceOrderMatcher.fit)


        * [`BruteForceOrderMatcher.match()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.BruteForceOrderMatcher.match)


        * [`BruteForceOrderMatcher.permutations()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.BruteForceOrderMatcher.permutations)


    * [`GeneticOrderMatcher`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.GeneticOrderMatcher)


        * [`GeneticOrderMatcher.fit()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.GeneticOrderMatcher.fit)


        * [`GeneticOrderMatcher.match()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.GeneticOrderMatcher.match)


        * [`GeneticOrderMatcher.permutations()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.GeneticOrderMatcher.permutations)


    * [`HungarianOrderMatcher`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.HungarianOrderMatcher)


        * [`HungarianOrderMatcher.fit()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.HungarianOrderMatcher.fit)


        * [`HungarianOrderMatcher.get_principal_axis()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.HungarianOrderMatcher.get_principal_axis)


        * [`HungarianOrderMatcher.match()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.HungarianOrderMatcher.match)


        * [`HungarianOrderMatcher.permutations()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.HungarianOrderMatcher.permutations)


        * [`HungarianOrderMatcher.rotation_matrix_vectors()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.HungarianOrderMatcher.rotation_matrix_vectors)


    * [`InchiMolAtomMapper`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.InchiMolAtomMapper)


        * [`InchiMolAtomMapper.as_dict()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.InchiMolAtomMapper.as_dict)


        * [`InchiMolAtomMapper.from_dict()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.InchiMolAtomMapper.from_dict)


        * [`InchiMolAtomMapper.get_molecule_hash()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.InchiMolAtomMapper.get_molecule_hash)


        * [`InchiMolAtomMapper.uniform_labels()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.InchiMolAtomMapper.uniform_labels)


    * [`IsomorphismMolAtomMapper`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.IsomorphismMolAtomMapper)


        * [`IsomorphismMolAtomMapper.as_dict()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.IsomorphismMolAtomMapper.as_dict)


        * [`IsomorphismMolAtomMapper.from_dict()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.IsomorphismMolAtomMapper.from_dict)


        * [`IsomorphismMolAtomMapper.get_molecule_hash()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.IsomorphismMolAtomMapper.get_molecule_hash)


        * [`IsomorphismMolAtomMapper.uniform_labels()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.IsomorphismMolAtomMapper.uniform_labels)


    * [`KabschMatcher`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.KabschMatcher)


        * [`KabschMatcher.fit()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.KabschMatcher.fit)


        * [`KabschMatcher.kabsch()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.KabschMatcher.kabsch)


        * [`KabschMatcher.match()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.KabschMatcher.match)


    * [`MoleculeMatcher`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.MoleculeMatcher)


        * [`MoleculeMatcher.as_dict()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.MoleculeMatcher.as_dict)


        * [`MoleculeMatcher.fit()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.MoleculeMatcher.fit)


        * [`MoleculeMatcher.from_dict()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.MoleculeMatcher.from_dict)


        * [`MoleculeMatcher.get_rmsd()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.MoleculeMatcher.get_rmsd)


        * [`MoleculeMatcher.group_molecules()`](pymatgen.analysis.molecule_matcher.md#pymatgen.analysis.molecule_matcher.MoleculeMatcher.group_molecules)


* [pymatgen.analysis.molecule_structure_comparator module](pymatgen.analysis.molecule_structure_comparator.md)


    * [`CovalentRadius`](pymatgen.analysis.molecule_structure_comparator.md#pymatgen.analysis.molecule_structure_comparator.CovalentRadius)


        * [`CovalentRadius.radius`](pymatgen.analysis.molecule_structure_comparator.md#pymatgen.analysis.molecule_structure_comparator.CovalentRadius.radius)


    * [`MoleculeStructureComparator`](pymatgen.analysis.molecule_structure_comparator.md#pymatgen.analysis.molecule_structure_comparator.MoleculeStructureComparator)


        * [`MoleculeStructureComparator.are_equal()`](pymatgen.analysis.molecule_structure_comparator.md#pymatgen.analysis.molecule_structure_comparator.MoleculeStructureComparator.are_equal)


        * [`MoleculeStructureComparator.as_dict()`](pymatgen.analysis.molecule_structure_comparator.md#pymatgen.analysis.molecule_structure_comparator.MoleculeStructureComparator.as_dict)


        * [`MoleculeStructureComparator.from_dict()`](pymatgen.analysis.molecule_structure_comparator.md#pymatgen.analysis.molecule_structure_comparator.MoleculeStructureComparator.from_dict)


        * [`MoleculeStructureComparator.get_13_bonds()`](pymatgen.analysis.molecule_structure_comparator.md#pymatgen.analysis.molecule_structure_comparator.MoleculeStructureComparator.get_13_bonds)


        * [`MoleculeStructureComparator.halogen_list`](pymatgen.analysis.molecule_structure_comparator.md#pymatgen.analysis.molecule_structure_comparator.MoleculeStructureComparator.halogen_list)


        * [`MoleculeStructureComparator.ionic_element_list`](pymatgen.analysis.molecule_structure_comparator.md#pymatgen.analysis.molecule_structure_comparator.MoleculeStructureComparator.ionic_element_list)


* [pymatgen.analysis.nmr module](pymatgen.analysis.nmr.md)


    * [`ChemicalShielding`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding)


        * [`ChemicalShielding.HaeberlenNotation`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.HaeberlenNotation)


            * [`ChemicalShielding.HaeberlenNotation.delta_sigma_iso`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.HaeberlenNotation.delta_sigma_iso)


            * [`ChemicalShielding.HaeberlenNotation.eta`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.HaeberlenNotation.eta)


            * [`ChemicalShielding.HaeberlenNotation.sigma_iso`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.HaeberlenNotation.sigma_iso)


            * [`ChemicalShielding.HaeberlenNotation.zeta`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.HaeberlenNotation.zeta)


        * [`ChemicalShielding.MarylandNotation`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.MarylandNotation)


            * [`ChemicalShielding.MarylandNotation.kappa`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.MarylandNotation.kappa)


            * [`ChemicalShielding.MarylandNotation.omega`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.MarylandNotation.omega)


            * [`ChemicalShielding.MarylandNotation.sigma_iso`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.MarylandNotation.sigma_iso)


        * [`ChemicalShielding.MehringNotation`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.MehringNotation)


            * [`ChemicalShielding.MehringNotation.sigma_11`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.MehringNotation.sigma_11)


            * [`ChemicalShielding.MehringNotation.sigma_22`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.MehringNotation.sigma_22)


            * [`ChemicalShielding.MehringNotation.sigma_33`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.MehringNotation.sigma_33)


            * [`ChemicalShielding.MehringNotation.sigma_iso`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.MehringNotation.sigma_iso)


        * [`ChemicalShielding.from_maryland_notation()`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.from_maryland_notation)


        * [`ChemicalShielding.haeberlen_values`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.haeberlen_values)


        * [`ChemicalShielding.maryland_values`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.maryland_values)


        * [`ChemicalShielding.mehring_values`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.mehring_values)


        * [`ChemicalShielding.principal_axis_system`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ChemicalShielding.principal_axis_system)


    * [`ElectricFieldGradient`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ElectricFieldGradient)


        * [`ElectricFieldGradient.V_xx`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ElectricFieldGradient.V_xx)


        * [`ElectricFieldGradient.V_yy`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ElectricFieldGradient.V_yy)


        * [`ElectricFieldGradient.V_zz`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ElectricFieldGradient.V_zz)


        * [`ElectricFieldGradient.asymmetry`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ElectricFieldGradient.asymmetry)


        * [`ElectricFieldGradient.coupling_constant()`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ElectricFieldGradient.coupling_constant)


        * [`ElectricFieldGradient.principal_axis_system`](pymatgen.analysis.nmr.md#pymatgen.analysis.nmr.ElectricFieldGradient.principal_axis_system)


* [pymatgen.analysis.path_finder module](pymatgen.analysis.path_finder.md)


    * [`ChgcarPotential`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.ChgcarPotential)


    * [`FreeVolumePotential`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.FreeVolumePotential)


    * [`MixedPotential`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.MixedPotential)


    * [`NEBPathfinder`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.NEBPathfinder)


        * [`NEBPathfinder.images`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.NEBPathfinder.images)


        * [`NEBPathfinder.interpolate()`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.NEBPathfinder.interpolate)


        * [`NEBPathfinder.plot_images()`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.NEBPathfinder.plot_images)


        * [`NEBPathfinder.string_relax()`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.NEBPathfinder.string_relax)


    * [`StaticPotential`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.StaticPotential)


        * [`StaticPotential.gaussian_smear()`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.StaticPotential.gaussian_smear)


        * [`StaticPotential.get_v()`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.StaticPotential.get_v)


        * [`StaticPotential.normalize()`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.StaticPotential.normalize)


        * [`StaticPotential.rescale_field()`](pymatgen.analysis.path_finder.md#pymatgen.analysis.path_finder.StaticPotential.rescale_field)


* [pymatgen.analysis.phase_diagram module](pymatgen.analysis.phase_diagram.md)


    * [`CompoundPhaseDiagram`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.CompoundPhaseDiagram)


        * [`CompoundPhaseDiagram.amount_tol`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.CompoundPhaseDiagram.amount_tol)


        * [`CompoundPhaseDiagram.as_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.CompoundPhaseDiagram.as_dict)


        * [`CompoundPhaseDiagram.from_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.CompoundPhaseDiagram.from_dict)


        * [`CompoundPhaseDiagram.transform_entries()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.CompoundPhaseDiagram.transform_entries)


    * [`GrandPotPDEntry`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.GrandPotPDEntry)


        * [`GrandPotPDEntry.as_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.GrandPotPDEntry.as_dict)


        * [`GrandPotPDEntry.chemical_energy`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.GrandPotPDEntry.chemical_energy)


        * [`GrandPotPDEntry.composition`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.GrandPotPDEntry.composition)


        * [`GrandPotPDEntry.energy`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.GrandPotPDEntry.energy)


        * [`GrandPotPDEntry.from_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.GrandPotPDEntry.from_dict)


    * [`GrandPotentialPhaseDiagram`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.GrandPotentialPhaseDiagram)


        * [`GrandPotentialPhaseDiagram.as_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.GrandPotentialPhaseDiagram.as_dict)


        * [`GrandPotentialPhaseDiagram.from_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.GrandPotentialPhaseDiagram.from_dict)


    * [`PDEntry`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry)


        * [`PDEntry.composition`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry.composition)


        * [`PDEntry.energy`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry.energy)


        * [`PDEntry.name`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry.name)


        * [`PDEntry.attribute`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry.attribute)


        * [`PDEntry.as_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry.as_dict)


        * [`PDEntry.energy`](pymatgen.analysis.phase_diagram.md#id0)


        * [`PDEntry.from_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDEntry.from_dict)


    * [`PDPlotter`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDPlotter)


        * [`PDPlotter.get_chempot_range_map_plot()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDPlotter.get_chempot_range_map_plot)


        * [`PDPlotter.get_contour_pd_plot()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDPlotter.get_contour_pd_plot)


        * [`PDPlotter.get_plot()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDPlotter.get_plot)


        * [`PDPlotter.pd_plot_data`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDPlotter.pd_plot_data)


        * [`PDPlotter.plot_chempot_range_map()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDPlotter.plot_chempot_range_map)


        * [`PDPlotter.plot_element_profile()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDPlotter.plot_element_profile)


        * [`PDPlotter.show()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDPlotter.show)


        * [`PDPlotter.write_image()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PDPlotter.write_image)


    * [`PatchedPhaseDiagram`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram)


        * [`PatchedPhaseDiagram.all_entries`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.all_entries)


        * [`PatchedPhaseDiagram.min_entries`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.min_entries)


        * [`PatchedPhaseDiagram.el_refs`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.el_refs)


        * [`PatchedPhaseDiagram.elements`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.elements)


        * [`PatchedPhaseDiagram.as_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.as_dict)


        * [`PatchedPhaseDiagram.from_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.from_dict)


        * [`PatchedPhaseDiagram.get_all_chempots()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_all_chempots)


        * [`PatchedPhaseDiagram.get_chempot_range_map()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_chempot_range_map)


        * [`PatchedPhaseDiagram.get_chempot_range_stability_phase()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_chempot_range_stability_phase)


        * [`PatchedPhaseDiagram.get_composition_chempots()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_composition_chempots)


        * [`PatchedPhaseDiagram.get_critical_compositions()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_critical_compositions)


        * [`PatchedPhaseDiagram.get_decomp_and_e_above_hull()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_decomp_and_e_above_hull)


        * [`PatchedPhaseDiagram.get_decomposition()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_decomposition)


        * [`PatchedPhaseDiagram.get_element_profile()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_element_profile)


        * [`PatchedPhaseDiagram.get_equilibrium_reaction_energy()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_equilibrium_reaction_energy)


        * [`PatchedPhaseDiagram.get_pd_for_entry()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_pd_for_entry)


        * [`PatchedPhaseDiagram.get_transition_chempots()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.get_transition_chempots)


        * [`PatchedPhaseDiagram.getmu_vertices_stability_phase()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PatchedPhaseDiagram.getmu_vertices_stability_phase)


    * [`PhaseDiagram`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram)


        * [`PhaseDiagram.dim`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.dim)


        * [`PhaseDiagram.elements`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.elements)


        * [`PhaseDiagram.el_refs`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.el_refs)


        * [`PhaseDiagram.all_entries`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.all_entries)


        * [`PhaseDiagram.qhull_entries`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.qhull_entries)


        * [`PhaseDiagram.qhull_data`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.qhull_data)


        * [`PhaseDiagram.facets`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.facets)


        * [`PhaseDiagram.simplices`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.simplices)


        * [`PhaseDiagram.all_entries_hulldata`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.all_entries_hulldata)


        * [`PhaseDiagram.as_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.as_dict)


        * [`PhaseDiagram.formation_energy_tol`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.formation_energy_tol)


        * [`PhaseDiagram.from_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.from_dict)


        * [`PhaseDiagram.get_all_chempots()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_all_chempots)


        * [`PhaseDiagram.get_chempot_range_map()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_chempot_range_map)


        * [`PhaseDiagram.get_chempot_range_stability_phase()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_chempot_range_stability_phase)


        * [`PhaseDiagram.get_composition_chempots()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_composition_chempots)


        * [`PhaseDiagram.get_critical_compositions()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_critical_compositions)


        * [`PhaseDiagram.get_decomp_and_e_above_hull()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_decomp_and_e_above_hull)


        * [`PhaseDiagram.get_decomp_and_hull_energy_per_atom()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_decomp_and_hull_energy_per_atom)


        * [`PhaseDiagram.get_decomp_and_phase_separation_energy()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_decomp_and_phase_separation_energy)


        * [`PhaseDiagram.get_decomposition()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_decomposition)


        * [`PhaseDiagram.get_e_above_hull()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_e_above_hull)


        * [`PhaseDiagram.get_element_profile()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_element_profile)


        * [`PhaseDiagram.get_equilibrium_reaction_energy()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_equilibrium_reaction_energy)


        * [`PhaseDiagram.get_form_energy()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_form_energy)


        * [`PhaseDiagram.get_form_energy_per_atom()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_form_energy_per_atom)


        * [`PhaseDiagram.get_hull_energy()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_hull_energy)


        * [`PhaseDiagram.get_hull_energy_per_atom()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_hull_energy_per_atom)


        * [`PhaseDiagram.get_phase_separation_energy()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_phase_separation_energy)


        * [`PhaseDiagram.get_plot()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_plot)


        * [`PhaseDiagram.get_reference_energy_per_atom()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_reference_energy_per_atom)


        * [`PhaseDiagram.get_transition_chempots()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.get_transition_chempots)


        * [`PhaseDiagram.getmu_vertices_stability_phase()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.getmu_vertices_stability_phase)


        * [`PhaseDiagram.numerical_tol`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.numerical_tol)


        * [`PhaseDiagram.pd_coords()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.pd_coords)


        * [`PhaseDiagram.stable_entries`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.stable_entries)


        * [`PhaseDiagram.unstable_entries`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagram.unstable_entries)


    * [`PhaseDiagramError`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.PhaseDiagramError)


    * [`ReactionDiagram`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.ReactionDiagram)


        * [`ReactionDiagram.get_compound_pd()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.ReactionDiagram.get_compound_pd)


    * [`TransformedPDEntry`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.TransformedPDEntry)


        * [`TransformedPDEntry.amount_tol`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.TransformedPDEntry.amount_tol)


        * [`TransformedPDEntry.as_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.TransformedPDEntry.as_dict)


        * [`TransformedPDEntry.composition`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.TransformedPDEntry.composition)


        * [`TransformedPDEntry.from_dict()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.TransformedPDEntry.from_dict)


    * [`TransformedPDEntryError`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.TransformedPDEntryError)


    * [`get_facets()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.get_facets)


    * [`order_phase_diagram()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.order_phase_diagram)


    * [`tet_coord()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.tet_coord)


    * [`triangular_coord()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.triangular_coord)


    * [`uniquelines()`](pymatgen.analysis.phase_diagram.md#pymatgen.analysis.phase_diagram.uniquelines)


* [pymatgen.analysis.piezo module](pymatgen.analysis.piezo.md)


    * [`PiezoTensor`](pymatgen.analysis.piezo.md#pymatgen.analysis.piezo.PiezoTensor)


        * [`PiezoTensor.from_vasp_voigt()`](pymatgen.analysis.piezo.md#pymatgen.analysis.piezo.PiezoTensor.from_vasp_voigt)


* [pymatgen.analysis.piezo_sensitivity module](pymatgen.analysis.piezo_sensitivity.md)


    * [`BornEffectiveCharge`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.BornEffectiveCharge)


        * [`BornEffectiveCharge.get_BEC_operations()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.BornEffectiveCharge.get_BEC_operations)


        * [`BornEffectiveCharge.get_rand_BEC()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.BornEffectiveCharge.get_rand_BEC)


    * [`ForceConstantMatrix`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.ForceConstantMatrix)


        * [`ForceConstantMatrix.get_FCM_operations()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.ForceConstantMatrix.get_FCM_operations)


        * [`ForceConstantMatrix.get_asum_FCM()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.ForceConstantMatrix.get_asum_FCM)


        * [`ForceConstantMatrix.get_rand_FCM()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.ForceConstantMatrix.get_rand_FCM)


        * [`ForceConstantMatrix.get_stable_FCM()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.ForceConstantMatrix.get_stable_FCM)


        * [`ForceConstantMatrix.get_symmetrized_FCM()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.ForceConstantMatrix.get_symmetrized_FCM)


        * [`ForceConstantMatrix.get_unstable_FCM()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.ForceConstantMatrix.get_unstable_FCM)


    * [`InternalStrainTensor`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.InternalStrainTensor)


        * [`InternalStrainTensor.get_IST_operations()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.InternalStrainTensor.get_IST_operations)


        * [`InternalStrainTensor.get_rand_IST()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.InternalStrainTensor.get_rand_IST)


    * [`get_piezo()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.get_piezo)


    * [`rand_piezo()`](pymatgen.analysis.piezo_sensitivity.md#pymatgen.analysis.piezo_sensitivity.rand_piezo)


* [pymatgen.analysis.pourbaix_diagram module](pymatgen.analysis.pourbaix_diagram.md)


    * [`IonEntry`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.IonEntry)


        * [`IonEntry.name`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.IonEntry.name)


        * [`IonEntry.as_dict()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.IonEntry.as_dict)


        * [`IonEntry.from_dict()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.IonEntry.from_dict)


    * [`MultiEntry`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.MultiEntry)


        * [`MultiEntry.as_dict()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.MultiEntry.as_dict)


        * [`MultiEntry.from_dict()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.MultiEntry.from_dict)


        * [`MultiEntry.name`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.MultiEntry.name)


    * [`PourbaixDiagram`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram)


        * [`PourbaixDiagram.all_entries`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.all_entries)


        * [`PourbaixDiagram.as_dict()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.as_dict)


        * [`PourbaixDiagram.find_stable_entry()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.find_stable_entry)


        * [`PourbaixDiagram.from_dict()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.from_dict)


        * [`PourbaixDiagram.get_decomposition_energy()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.get_decomposition_energy)


        * [`PourbaixDiagram.get_hull_energy()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.get_hull_energy)


        * [`PourbaixDiagram.get_pourbaix_domains()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.get_pourbaix_domains)


        * [`PourbaixDiagram.get_stable_entry()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.get_stable_entry)


        * [`PourbaixDiagram.process_multientry()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.process_multientry)


        * [`PourbaixDiagram.stable_entries`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.stable_entries)


        * [`PourbaixDiagram.unprocessed_entries`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.unprocessed_entries)


        * [`PourbaixDiagram.unstable_entries`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixDiagram.unstable_entries)


    * [`PourbaixEntry`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry)


        * [`PourbaixEntry.as_dict()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.as_dict)


        * [`PourbaixEntry.composition`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.composition)


        * [`PourbaixEntry.conc_term`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.conc_term)


        * [`PourbaixEntry.energy`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.energy)


        * [`PourbaixEntry.energy_at_conditions()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.energy_at_conditions)


        * [`PourbaixEntry.energy_per_atom`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.energy_per_atom)


        * [`PourbaixEntry.from_dict()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.from_dict)


        * [`PourbaixEntry.get_element_fraction()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.get_element_fraction)


        * [`PourbaixEntry.nH2O`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.nH2O)


        * [`PourbaixEntry.nPhi`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.nPhi)


        * [`PourbaixEntry.name`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.name)


        * [`PourbaixEntry.normalization_factor`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.normalization_factor)


        * [`PourbaixEntry.normalized_energy`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.normalized_energy)


        * [`PourbaixEntry.normalized_energy_at_conditions()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.normalized_energy_at_conditions)


        * [`PourbaixEntry.npH`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.npH)


        * [`PourbaixEntry.num_atoms`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.num_atoms)


        * [`PourbaixEntry.to_pretty_string()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixEntry.to_pretty_string)


    * [`PourbaixPlotter`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixPlotter)


        * [`PourbaixPlotter.domain_vertices()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixPlotter.domain_vertices)


        * [`PourbaixPlotter.get_pourbaix_plot()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixPlotter.get_pourbaix_plot)


        * [`PourbaixPlotter.plot_entry_stability()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixPlotter.plot_entry_stability)


        * [`PourbaixPlotter.show()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.PourbaixPlotter.show)


    * [`generate_entry_label()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.generate_entry_label)


    * [`ion_or_solid_comp_object()`](pymatgen.analysis.pourbaix_diagram.md#pymatgen.analysis.pourbaix_diagram.ion_or_solid_comp_object)


* [pymatgen.analysis.prototypes module](pymatgen.analysis.prototypes.md)


    * [`AflowPrototypeMatcher`](pymatgen.analysis.prototypes.md#pymatgen.analysis.prototypes.AflowPrototypeMatcher)


        * [`AflowPrototypeMatcher.get_prototypes()`](pymatgen.analysis.prototypes.md#pymatgen.analysis.prototypes.AflowPrototypeMatcher.get_prototypes)


* [pymatgen.analysis.quasiharmonic module](pymatgen.analysis.quasiharmonic.md)


    * [`QuasiharmonicDebyeApprox`](pymatgen.analysis.quasiharmonic.md#pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox)


        * [`QuasiharmonicDebyeApprox.debye_integral()`](pymatgen.analysis.quasiharmonic.md#pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox.debye_integral)


        * [`QuasiharmonicDebyeApprox.debye_temperature()`](pymatgen.analysis.quasiharmonic.md#pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox.debye_temperature)


        * [`QuasiharmonicDebyeApprox.get_summary_dict()`](pymatgen.analysis.quasiharmonic.md#pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox.get_summary_dict)


        * [`QuasiharmonicDebyeApprox.gruneisen_parameter()`](pymatgen.analysis.quasiharmonic.md#pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox.gruneisen_parameter)


        * [`QuasiharmonicDebyeApprox.optimize_gibbs_free_energy()`](pymatgen.analysis.quasiharmonic.md#pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox.optimize_gibbs_free_energy)


        * [`QuasiharmonicDebyeApprox.optimizer()`](pymatgen.analysis.quasiharmonic.md#pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox.optimizer)


        * [`QuasiharmonicDebyeApprox.thermal_conductivity()`](pymatgen.analysis.quasiharmonic.md#pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox.thermal_conductivity)


        * [`QuasiharmonicDebyeApprox.vibrational_free_energy()`](pymatgen.analysis.quasiharmonic.md#pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox.vibrational_free_energy)


        * [`QuasiharmonicDebyeApprox.vibrational_internal_energy()`](pymatgen.analysis.quasiharmonic.md#pymatgen.analysis.quasiharmonic.QuasiharmonicDebyeApprox.vibrational_internal_energy)


* [pymatgen.analysis.reaction_calculator module](pymatgen.analysis.reaction_calculator.md)


    * [`BalancedReaction`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction)


        * [`BalancedReaction.TOLERANCE`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.TOLERANCE)


        * [`BalancedReaction.all_comp`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.all_comp)


        * [`BalancedReaction.as_dict()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.as_dict)


        * [`BalancedReaction.as_entry()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.as_entry)


        * [`BalancedReaction.calculate_energy()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.calculate_energy)


        * [`BalancedReaction.coeffs`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.coeffs)


        * [`BalancedReaction.elements`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.elements)


        * [`BalancedReaction.from_dict()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.from_dict)


        * [`BalancedReaction.from_str()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.from_str)


        * [`BalancedReaction.from_string()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.from_string)


        * [`BalancedReaction.get_coeff()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.get_coeff)


        * [`BalancedReaction.get_el_amount()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.get_el_amount)


        * [`BalancedReaction.normalize_to()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.normalize_to)


        * [`BalancedReaction.normalize_to_element()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.normalize_to_element)


        * [`BalancedReaction.normalized_repr`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.normalized_repr)


        * [`BalancedReaction.normalized_repr_and_factor()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.normalized_repr_and_factor)


        * [`BalancedReaction.products`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.products)


        * [`BalancedReaction.reactants`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.BalancedReaction.reactants)


    * [`ComputedReaction`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.ComputedReaction)


        * [`ComputedReaction.all_entries`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.ComputedReaction.all_entries)


        * [`ComputedReaction.as_dict()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.ComputedReaction.as_dict)


        * [`ComputedReaction.calculated_reaction_energy`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.ComputedReaction.calculated_reaction_energy)


        * [`ComputedReaction.calculated_reaction_energy_uncertainty`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.ComputedReaction.calculated_reaction_energy_uncertainty)


        * [`ComputedReaction.from_dict()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.ComputedReaction.from_dict)


    * [`Reaction`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.Reaction)


        * [`Reaction.as_dict()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.Reaction.as_dict)


        * [`Reaction.copy()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.Reaction.copy)


        * [`Reaction.from_dict()`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.Reaction.from_dict)


    * [`ReactionError`](pymatgen.analysis.reaction_calculator.md#pymatgen.analysis.reaction_calculator.ReactionError)


* [pymatgen.analysis.structure_analyzer module](pymatgen.analysis.structure_analyzer.md)


    * [`OxideType`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.OxideType)


        * [`OxideType.parse_oxide()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.OxideType.parse_oxide)


    * [`RelaxationAnalyzer`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.RelaxationAnalyzer)


        * [`RelaxationAnalyzer.get_percentage_bond_dist_changes()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.RelaxationAnalyzer.get_percentage_bond_dist_changes)


        * [`RelaxationAnalyzer.get_percentage_lattice_parameter_changes()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.RelaxationAnalyzer.get_percentage_lattice_parameter_changes)


        * [`RelaxationAnalyzer.get_percentage_volume_change()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.RelaxationAnalyzer.get_percentage_volume_change)


    * [`VoronoiAnalyzer`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.VoronoiAnalyzer)


        * [`VoronoiAnalyzer.analyze()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.VoronoiAnalyzer.analyze)


        * [`VoronoiAnalyzer.analyze_structures()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.VoronoiAnalyzer.analyze_structures)


        * [`VoronoiAnalyzer.plot_vor_analysis()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.VoronoiAnalyzer.plot_vor_analysis)


    * [`VoronoiConnectivity`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.VoronoiConnectivity)


        * [`VoronoiConnectivity.connectivity_array`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.VoronoiConnectivity.connectivity_array)


        * [`VoronoiConnectivity.get_connections()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.VoronoiConnectivity.get_connections)


        * [`VoronoiConnectivity.get_sitej()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.VoronoiConnectivity.get_sitej)


        * [`VoronoiConnectivity.max_connectivity`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.VoronoiConnectivity.max_connectivity)


    * [`average_coordination_number()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.average_coordination_number)


    * [`contains_peroxide()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.contains_peroxide)


    * [`get_max_bond_lengths()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.get_max_bond_lengths)


    * [`oxide_type()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.oxide_type)


    * [`solid_angle()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.solid_angle)


    * [`sulfide_type()`](pymatgen.analysis.structure_analyzer.md#pymatgen.analysis.structure_analyzer.sulfide_type)


* [pymatgen.analysis.structure_matcher module](pymatgen.analysis.structure_matcher.md)


    * [`AbstractComparator`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.AbstractComparator)


        * [`AbstractComparator.are_equal()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.AbstractComparator.are_equal)


        * [`AbstractComparator.as_dict()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.AbstractComparator.as_dict)


        * [`AbstractComparator.from_dict()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.AbstractComparator.from_dict)


        * [`AbstractComparator.get_hash()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.AbstractComparator.get_hash)


    * [`ElementComparator`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.ElementComparator)


        * [`ElementComparator.are_equal()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.ElementComparator.are_equal)


        * [`ElementComparator.get_hash()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.ElementComparator.get_hash)


    * [`FrameworkComparator`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.FrameworkComparator)


        * [`FrameworkComparator.are_equal()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.FrameworkComparator.are_equal)


        * [`FrameworkComparator.get_hash()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.FrameworkComparator.get_hash)


    * [`OccupancyComparator`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.OccupancyComparator)


        * [`OccupancyComparator.are_equal()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.OccupancyComparator.are_equal)


        * [`OccupancyComparator.get_hash()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.OccupancyComparator.get_hash)


    * [`OrderDisorderElementComparator`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.OrderDisorderElementComparator)


        * [`OrderDisorderElementComparator.are_equal()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.OrderDisorderElementComparator.are_equal)


        * [`OrderDisorderElementComparator.get_hash()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.OrderDisorderElementComparator.get_hash)


    * [`SpeciesComparator`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.SpeciesComparator)


        * [`SpeciesComparator.are_equal()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.SpeciesComparator.are_equal)


        * [`SpeciesComparator.get_hash()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.SpeciesComparator.get_hash)


    * [`SpinComparator`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.SpinComparator)


        * [`SpinComparator.are_equal()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.SpinComparator.are_equal)


        * [`SpinComparator.get_hash()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.SpinComparator.get_hash)


    * [`StructureMatcher`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher)


        * [`StructureMatcher.as_dict()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.as_dict)


        * [`StructureMatcher.fit()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.fit)


        * [`StructureMatcher.fit_anonymous()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.fit_anonymous)


        * [`StructureMatcher.from_dict()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.from_dict)


        * [`StructureMatcher.get_all_anonymous_mappings()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.get_all_anonymous_mappings)


        * [`StructureMatcher.get_best_electronegativity_anonymous_mapping()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.get_best_electronegativity_anonymous_mapping)


        * [`StructureMatcher.get_mapping()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.get_mapping)


        * [`StructureMatcher.get_rms_anonymous()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.get_rms_anonymous)


        * [`StructureMatcher.get_rms_dist()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.get_rms_dist)


        * [`StructureMatcher.get_s2_like_s1()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.get_s2_like_s1)


        * [`StructureMatcher.get_supercell_matrix()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.get_supercell_matrix)


        * [`StructureMatcher.get_transformation()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.get_transformation)


        * [`StructureMatcher.group_structures()`](pymatgen.analysis.structure_matcher.md#pymatgen.analysis.structure_matcher.StructureMatcher.group_structures)


* [pymatgen.analysis.surface_analysis module](pymatgen.analysis.surface_analysis.md)


    * [`NanoscaleStability`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.NanoscaleStability)


        * [`NanoscaleStability.se_analyzers`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.NanoscaleStability.se_analyzers)


        * [`NanoscaleStability.symprec`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.NanoscaleStability.symprec)


        * [`NanoscaleStability.bulk_gform()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.NanoscaleStability.bulk_gform)


        * [`NanoscaleStability.plot_all_stability_map()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.NanoscaleStability.plot_all_stability_map)


        * [`NanoscaleStability.plot_one_stability_map()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.NanoscaleStability.plot_one_stability_map)


        * [`NanoscaleStability.scaled_wulff()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.NanoscaleStability.scaled_wulff)


        * [`NanoscaleStability.solve_equilibrium_point()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.NanoscaleStability.solve_equilibrium_point)


        * [`NanoscaleStability.wulff_gform_and_r()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.NanoscaleStability.wulff_gform_and_r)


    * [`SlabEntry`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry)


        * [`SlabEntry.miller_index`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.miller_index)


        * [`SlabEntry.label`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.label)


        * [`SlabEntry.adsorbates`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.adsorbates)


        * [`SlabEntry.Nads_in_slab`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.Nads_in_slab)


        * [`SlabEntry.Nsurfs_ads_in_slab`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.Nsurfs_ads_in_slab)


        * [`SlabEntry.as_dict()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.as_dict)


        * [`SlabEntry.cleaned_up_slab`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.cleaned_up_slab)


        * [`SlabEntry.create_slab_label`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.create_slab_label)


        * [`SlabEntry.from_computed_structure_entry()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.from_computed_structure_entry)


        * [`SlabEntry.from_dict()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.from_dict)


        * [`SlabEntry.get_monolayer`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.get_monolayer)


        * [`SlabEntry.get_unit_primitive_area`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.get_unit_primitive_area)


        * [`SlabEntry.gibbs_binding_energy()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.gibbs_binding_energy)


        * [`SlabEntry.surface_area`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.surface_area)


        * [`SlabEntry.surface_energy()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SlabEntry.surface_energy)


    * [`SurfaceEnergyPlotter`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter)


        * [`SurfaceEnergyPlotter.all_slab_entries`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.all_slab_entries)


        * [`SurfaceEnergyPlotter.ucell_entry`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.ucell_entry)


        * [`SurfaceEnergyPlotter.ref_entries`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.ref_entries)


        * [`SurfaceEnergyPlotter.color_dict`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.color_dict)


        * [`SurfaceEnergyPlotter.BE_vs_clean_SE()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.BE_vs_clean_SE)


        * [`SurfaceEnergyPlotter.area_frac_vs_chempot_plot()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.area_frac_vs_chempot_plot)


        * [`SurfaceEnergyPlotter.chempot_plot_addons()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.chempot_plot_addons)


        * [`SurfaceEnergyPlotter.chempot_vs_gamma()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.chempot_vs_gamma)


        * [`SurfaceEnergyPlotter.chempot_vs_gamma_plot_one()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.chempot_vs_gamma_plot_one)


        * [`SurfaceEnergyPlotter.color_palette_dict()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.color_palette_dict)


        * [`SurfaceEnergyPlotter.get_stable_entry_at_u()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.get_stable_entry_at_u)


        * [`SurfaceEnergyPlotter.get_surface_equilibrium()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.get_surface_equilibrium)


        * [`SurfaceEnergyPlotter.monolayer_vs_BE()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.monolayer_vs_BE)


        * [`SurfaceEnergyPlotter.set_all_variables()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.set_all_variables)


        * [`SurfaceEnergyPlotter.stable_u_range_dict()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.stable_u_range_dict)


        * [`SurfaceEnergyPlotter.surface_chempot_range_map()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.surface_chempot_range_map)


        * [`SurfaceEnergyPlotter.wulff_from_chempot()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.SurfaceEnergyPlotter.wulff_from_chempot)


    * [`WorkFunctionAnalyzer`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer)


        * [`WorkFunctionAnalyzer.efermi`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.efermi)


        * [`WorkFunctionAnalyzer.locpot_along_c`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.locpot_along_c)


        * [`WorkFunctionAnalyzer.vacuum_locpot`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.vacuum_locpot)


        * [`WorkFunctionAnalyzer.work_function`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.work_function)


        * [`WorkFunctionAnalyzer.slab`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.slab)


        * [`WorkFunctionAnalyzer.along_c`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.along_c)


        * [`WorkFunctionAnalyzer.ave_locpot`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.ave_locpot)


        * [`WorkFunctionAnalyzer.sorted_sites`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.sorted_sites)


        * [`WorkFunctionAnalyzer.ave_bulk_p`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.ave_bulk_p)


        * [`WorkFunctionAnalyzer.from_files()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.from_files)


        * [`WorkFunctionAnalyzer.get_labels()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.get_labels)


        * [`WorkFunctionAnalyzer.get_locpot_along_slab_plot()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.get_locpot_along_slab_plot)


        * [`WorkFunctionAnalyzer.is_converged()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.WorkFunctionAnalyzer.is_converged)


    * [`entry_dict_from_list()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.entry_dict_from_list)


    * [`sub_chempots()`](pymatgen.analysis.surface_analysis.md#pymatgen.analysis.surface_analysis.sub_chempots)


* [pymatgen.analysis.thermochemistry module](pymatgen.analysis.thermochemistry.md)


    * [`ThermoData`](pymatgen.analysis.thermochemistry.md#pymatgen.analysis.thermochemistry.ThermoData)


        * [`ThermoData.as_dict()`](pymatgen.analysis.thermochemistry.md#pymatgen.analysis.thermochemistry.ThermoData.as_dict)


        * [`ThermoData.from_dict()`](pymatgen.analysis.thermochemistry.md#pymatgen.analysis.thermochemistry.ThermoData.from_dict)


* [pymatgen.analysis.transition_state module](pymatgen.analysis.transition_state.md)


    * [`NEBAnalysis`](pymatgen.analysis.transition_state.md#pymatgen.analysis.transition_state.NEBAnalysis)


        * [`NEBAnalysis.as_dict()`](pymatgen.analysis.transition_state.md#pymatgen.analysis.transition_state.NEBAnalysis.as_dict)


        * [`NEBAnalysis.from_dir()`](pymatgen.analysis.transition_state.md#pymatgen.analysis.transition_state.NEBAnalysis.from_dir)


        * [`NEBAnalysis.from_outcars()`](pymatgen.analysis.transition_state.md#pymatgen.analysis.transition_state.NEBAnalysis.from_outcars)


        * [`NEBAnalysis.get_extrema()`](pymatgen.analysis.transition_state.md#pymatgen.analysis.transition_state.NEBAnalysis.get_extrema)


        * [`NEBAnalysis.get_plot()`](pymatgen.analysis.transition_state.md#pymatgen.analysis.transition_state.NEBAnalysis.get_plot)


        * [`NEBAnalysis.setup_spline()`](pymatgen.analysis.transition_state.md#pymatgen.analysis.transition_state.NEBAnalysis.setup_spline)


    * [`combine_neb_plots()`](pymatgen.analysis.transition_state.md#pymatgen.analysis.transition_state.combine_neb_plots)


* [pymatgen.analysis.wulff module](pymatgen.analysis.wulff.md)


    * [`WulffFacet`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffFacet)


    * [`WulffShape`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape)


        * [`WulffShape.debug`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.debug)


        * [`WulffShape.alpha`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.alpha)


        * [`WulffShape.transparency`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.transparency)


        * [`WulffShape.color_set`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.color_set)


        * [`WulffShape.grid_off`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.grid_off)


        * [`WulffShape.axis_off`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.axis_off)


        * [`WulffShape.show_area`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.show_area)


        * [`WulffShape.off_color`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.off_color)


        * [`WulffShape.structure`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.structure)


        * [`WulffShape.miller_list`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.miller_list)


        * [`WulffShape.hkl_list`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.hkl_list)


        * [`WulffShape.e_surf_list`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.e_surf_list)


        * [`WulffShape.lattice`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.lattice)


        * [`WulffShape.facets`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.facets)


        * [`WulffShape.dual_cv_simp`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.dual_cv_simp)


        * [`WulffShape.wulff_pt_list`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.wulff_pt_list)


        * [`WulffShape.wulff_cv_simp`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.wulff_cv_simp)


        * [`WulffShape.on_wulff`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.on_wulff)


        * [`WulffShape.color_area`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.color_area)


        * [`WulffShape.miller_area`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.miller_area)


        * [`WulffShape.anisotropy`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.anisotropy)


        * [`WulffShape.area_fraction_dict`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.area_fraction_dict)


        * [`WulffShape.effective_radius`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.effective_radius)


        * [`WulffShape.get_line_in_facet()`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.get_line_in_facet)


        * [`WulffShape.get_plot()`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.get_plot)


        * [`WulffShape.get_plotly()`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.get_plotly)


        * [`WulffShape.miller_area_dict`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.miller_area_dict)


        * [`WulffShape.miller_energy_dict`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.miller_energy_dict)


        * [`WulffShape.shape_factor`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.shape_factor)


        * [`WulffShape.show()`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.show)


        * [`WulffShape.surface_area`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.surface_area)


        * [`WulffShape.tot_corner_sites`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.tot_corner_sites)


        * [`WulffShape.tot_edges`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.tot_edges)


        * [`WulffShape.total_surface_energy`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.total_surface_energy)


        * [`WulffShape.volume`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.volume)


        * [`WulffShape.weighted_surface_energy`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.WulffShape.weighted_surface_energy)


    * [`get_tri_area()`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.get_tri_area)


    * [`hkl_tuple_to_str()`](pymatgen.analysis.wulff.md#pymatgen.analysis.wulff.hkl_tuple_to_str)


* [pymatgen.analysis.xps module](pymatgen.analysis.xps.md)


    * [`XPS`](pymatgen.analysis.xps.md#pymatgen.analysis.xps.XPS)


        * [`XPS.XLABEL`](pymatgen.analysis.xps.md#pymatgen.analysis.xps.XPS.XLABEL)


        * [`XPS.YLABEL`](pymatgen.analysis.xps.md#pymatgen.analysis.xps.XPS.YLABEL)


        * [`XPS.from_dos()`](pymatgen.analysis.xps.md#pymatgen.analysis.xps.XPS.from_dos)