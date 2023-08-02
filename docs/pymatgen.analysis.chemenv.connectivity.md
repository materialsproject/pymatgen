---
layout: default
title: pymatgen.analysis.chemenv.connectivity.md
nav_exclude: true
---

# pymatgen.analysis.chemenv.connectivity package

Package for analyzing connectivity.



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