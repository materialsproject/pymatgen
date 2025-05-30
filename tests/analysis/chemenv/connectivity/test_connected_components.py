from __future__ import annotations

import copy

import networkx as nx
import numpy as np
import orjson
import pytest
from numpy.testing import assert_allclose

from pymatgen.analysis.chemenv.connectivity.connected_components import ConnectedComponent
from pymatgen.analysis.chemenv.connectivity.connectivity_finder import ConnectivityFinder
from pymatgen.analysis.chemenv.connectivity.environment_nodes import EnvironmentNode
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import (
    LightStructureEnvironments,
    StructureEnvironments,
)
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.util.testing import TEST_FILES_DIR, MatSciTest

__author__ = "waroquiers"


class TestConnectedComponent(MatSciTest):
    def test_init(self):
        # Generic connected component not using EnvironmentNodes
        # (as_dict won't work on such a ConnectedComponent instance)
        cc = ConnectedComponent(
            environments=["a", "b", "c", "d", "e", "f"],
            links=[
                ("a", "b", 0),
                ("a", "c", 0),
                ("b", "c", 0),
                ("a", "a", 0),
                ("a", "b", 1),
                ("c", "d"),
                ("c", "d"),
                ("c", "d"),
                ("d", "e"),
                ("e", "f"),
            ],
            environments_data={"a": {"var1": 2, "var2": 3}, "b": {"var1": 3}},
            links_data={
                ("c", "b"): {"bcinfo": 2},
                ("a", "b", 1): {"ab1info": 4},
                ("d", "c"): {"dcinfo": 8},
            },
        )
        assert isinstance(cc.graph, nx.MultiGraph)
        nodes = list(cc.graph.nodes())
        assert set(nodes) == {"a", "b", "c", "d", "e", "f"}
        edges = cc.graph.edges()
        assert len(edges) == 10
        assert cc.graph["a"]["b"] == {0: {}, 1: {"ab1info": 4}}
        assert cc.graph["b"]["a"] == {0: {}, 1: {"ab1info": 4}}
        assert cc.graph["c"]["b"] == {0: {"bcinfo": 2}}
        assert cc.graph["b"]["c"] == {0: {"bcinfo": 2}}
        assert len(cc.graph["c"]["d"]) == 3
        assert cc.graph["c"]["d"][0] == {"dcinfo": 8}
        assert cc.graph["c"]["d"][1] == {"dcinfo": 8}
        assert cc.graph["c"]["d"][2] == {"dcinfo": 8}

        assert cc.graph.nodes(data=True)["a"] == {"var1": 2, "var2": 3}
        assert cc.graph.nodes(data=True)["b"] == {"var1": 3}
        assert cc.graph.nodes(data=True)["c"] == {}

        # Make a deep copy of the graph to actually check that the objects are different instances
        graph = copy.deepcopy(cc.graph)
        assert isinstance(graph, nx.MultiGraph)  # Check that it is indeed the same type of graph

        cc2 = ConnectedComponent(graph=graph)
        assert set(graph.nodes()) == set(cc2.graph.nodes())
        assert set(graph.edges()) == set(cc2.graph.edges())
        assert len(cc2.graph) == 6

    def test_serialization(self):
        lattice = Lattice.hexagonal(a=2.0, c=2.5)
        en1 = EnvironmentNode(
            central_site=PeriodicSite("Si", coords=np.array([0.0, 0.0, 0.0]), lattice=lattice),
            i_central_site=3,
            ce_symbol="T:4",
        )
        en2 = EnvironmentNode(
            central_site=PeriodicSite("Ag", coords=np.array([0.0, 0.0, 0.5]), lattice=lattice),
            i_central_site=5,
            ce_symbol="T:4",
        )
        en3 = EnvironmentNode(
            central_site=PeriodicSite("Ag", coords=np.array([0.0, 0.5, 0.5]), lattice=lattice),
            i_central_site=8,
            ce_symbol="O:6",
        )

        graph = nx.MultiGraph()
        graph.add_nodes_from([en1, en2, en3])

        graph.add_edge(
            en1,
            en2,
            start=en1.isite,
            end=en2.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            en1,
            en3,
            start=en1.isite,
            end=en2.isite,
            delta=(0, 0, 0),
            ligands=[(10, (0, 0, 1), (0, 0, 1)), (11, (0, 0, 1), (0, 0, 1))],
        )

        cc = ConnectedComponent(graph=graph)
        ref_sorted_edges = [[en1, en2], [en1, en3]]
        sorted_edges = sorted(sorted(edge) for edge in cc.graph.edges())
        assert sorted_edges == ref_sorted_edges

        cc_from_dict = ConnectedComponent.from_dict(cc.as_dict())
        cc_from_json = ConnectedComponent.from_dict(orjson.loads(orjson.dumps(cc.as_dict()).decode()))
        loaded_cc_list = [cc_from_dict, cc_from_json]
        json_str = self.assert_msonable(cc)
        cc_from_json = ConnectedComponent.from_dict(orjson.loads(json_str))
        loaded_cc_list.append(cc_from_json)
        for loaded_cc in loaded_cc_list:
            assert loaded_cc.graph.number_of_nodes() == 3
            assert loaded_cc.graph.number_of_edges() == 2
            assert set(cc.graph.nodes()) == set(loaded_cc.graph.nodes())
            assert sorted_edges == sorted(sorted(edge) for edge in loaded_cc.graph.edges())

            for edge in sorted_edges:
                assert cc.graph[edge[0]][edge[1]] == loaded_cc.graph[edge[0]][edge[1]]

            for node in loaded_cc.graph.nodes():
                assert isinstance(node.central_site, PeriodicSite)

    def test_serialization_private_methods(self):
        # Testing _edge_key_to_edge_dict_key
        key = ConnectedComponent._edge_key_to_edge_dict_key(3)
        assert key == "3"
        with pytest.raises(
            RuntimeError,
            match=r"Cannot pass an edge key which is a str representation of an int\x2E",
        ):
            key = ConnectedComponent._edge_key_to_edge_dict_key("5")
        key = ConnectedComponent._edge_key_to_edge_dict_key("some-key")
        assert key == "some-key"
        with pytest.raises(ValueError, match=r"Edge key should be either a str or an int\x2E"):
            key = ConnectedComponent._edge_key_to_edge_dict_key(0.2)

    def test_periodicity(self):
        env_node1 = EnvironmentNode(central_site="Si", i_central_site=3, ce_symbol="T:4")
        env_node2 = EnvironmentNode(central_site="Ag", i_central_site=5, ce_symbol="T:4")
        env_node3 = EnvironmentNode(central_site="Ag", i_central_site=8, ce_symbol="O:6")
        env_node4 = EnvironmentNode(central_site="Fe", i_central_site=23, ce_symbol="C:8")

        graph = nx.MultiGraph()
        graph.add_nodes_from([env_node1, env_node2, env_node3])
        graph.add_edge(
            env_node1,
            env_node2,
            start=env_node1.isite,
            end=env_node2.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node1,
            env_node3,
            start=env_node1.isite,
            end=env_node2.isite,
            delta=(0, 0, 0),
            ligands=[(10, (0, 0, 1), (0, 0, 1)), (11, (0, 0, 1), (0, 0, 1))],
        )
        cc = ConnectedComponent(graph=graph)
        assert cc.is_0d
        assert not cc.is_1d
        assert not cc.is_2d
        assert not cc.is_3d
        assert not cc.is_periodic
        assert cc.periodicity == "0D"

        graph = nx.MultiGraph()
        graph.add_nodes_from([env_node1, env_node2, env_node3])
        graph.add_edge(
            env_node1,
            env_node2,
            start=env_node1.isite,
            end=env_node2.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node1,
            env_node3,
            start=env_node1.isite,
            end=env_node3.isite,
            delta=(0, 0, 0),
            ligands=[(10, (0, 0, 1), (0, 0, 1)), (11, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node2,
            env_node3,
            start=env_node2.isite,
            end=env_node3.isite,
            delta=(0, 0, 1),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        cc = ConnectedComponent(graph=graph)
        assert not cc.is_0d
        assert cc.is_1d
        assert not cc.is_2d
        assert not cc.is_3d
        assert cc.is_periodic
        assert cc.periodicity == "1D"

        graph = nx.MultiGraph()
        graph.add_nodes_from([env_node1, env_node2, env_node3])
        graph.add_edge(
            env_node1,
            env_node2,
            start=env_node1.isite,
            end=env_node2.isite,
            delta=(0, 0, 1),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node1,
            env_node3,
            start=env_node1.isite,
            end=env_node3.isite,
            delta=(0, 0, 0),
            ligands=[(10, (0, 0, 1), (0, 0, 1)), (11, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node2,
            env_node3,
            start=env_node2.isite,
            end=env_node3.isite,
            delta=(0, 0, -1),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        cc = ConnectedComponent(graph=graph)
        assert cc.periodicity == "0D"

        # Test errors when computing periodicity
        graph = nx.MultiGraph()
        graph.add_nodes_from([env_node1, env_node2, env_node3])
        graph.add_edge(
            env_node1,
            env_node1,
            start=env_node1.isite,
            end=env_node1.isite,
            delta=(0, 0, 1),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node1,
            env_node1,
            start=env_node1.isite,
            end=env_node1.isite,
            delta=(0, 0, 1),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        cc = ConnectedComponent(graph=graph)
        with pytest.raises(
            ValueError,
            match=r"There should not be self loops with the same \x28or opposite\x29 delta image\x2E",
        ):
            cc.compute_periodicity_all_simple_paths_algorithm()

        graph = nx.MultiGraph()
        graph.add_nodes_from([env_node1, env_node2, env_node3])
        graph.add_edge(
            env_node1,
            env_node1,
            start=env_node1.isite,
            end=env_node1.isite,
            delta=(3, 2, -1),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node1,
            env_node1,
            start=env_node1.isite,
            end=env_node1.isite,
            delta=(-3, -2, 1),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        cc = ConnectedComponent(graph=graph)
        with pytest.raises(
            ValueError,
            match=r"There should not be self loops with the same \x28or opposite\x29 delta image\x2E",
        ):
            cc.compute_periodicity_all_simple_paths_algorithm()

        graph = nx.MultiGraph()
        graph.add_nodes_from([env_node1, env_node2, env_node3])
        graph.add_edge(
            env_node1,
            env_node1,
            start=env_node1.isite,
            end=env_node1.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        cc = ConnectedComponent(graph=graph)
        with pytest.raises(
            ValueError,
            match=r"There should not be self loops with delta image = \x280, 0, 0\x29\x2E",
        ):
            cc.compute_periodicity_all_simple_paths_algorithm()

        # Test a 2d periodicity
        graph = nx.MultiGraph()
        graph.add_nodes_from([env_node1, env_node2, env_node3, env_node4])
        graph.add_edge(
            env_node1,
            env_node2,
            start=env_node1.isite,
            end=env_node2.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node1,
            env_node3,
            start=env_node1.isite,
            end=env_node3.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node4,
            env_node2,
            start=env_node4.isite,
            end=env_node2.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node3,
            env_node4,
            start=env_node4.isite,
            end=env_node3.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node3,
            env_node4,
            start=env_node4.isite,
            end=env_node3.isite,
            delta=(0, -1, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node3,
            env_node2,
            start=env_node2.isite,
            end=env_node3.isite,
            delta=(-1, -1, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        cc = ConnectedComponent(graph=graph)
        assert not cc.is_0d
        assert not cc.is_1d
        assert cc.is_2d
        assert not cc.is_3d
        assert cc.is_periodic
        assert cc.periodicity == "2D"
        assert_allclose(cc.periodicity_vectors, [[0, 1, 0], [1, 1, 0]])
        assert isinstance(cc.periodicity_vectors, list)
        assert cc.periodicity_vectors[0].dtype is np.dtype(np.int64)

        # Test a 3d periodicity
        graph = nx.MultiGraph()
        graph.add_nodes_from([env_node1, env_node2, env_node3, env_node4])
        graph.add_edge(
            env_node1,
            env_node2,
            start=env_node1.isite,
            end=env_node2.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node1,
            env_node3,
            start=env_node1.isite,
            end=env_node3.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node4,
            env_node2,
            start=env_node4.isite,
            end=env_node2.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node3,
            env_node4,
            start=env_node4.isite,
            end=env_node3.isite,
            delta=(0, 0, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node3,
            env_node4,
            start=env_node4.isite,
            end=env_node3.isite,
            delta=(0, -1, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node3,
            env_node2,
            start=env_node2.isite,
            end=env_node3.isite,
            delta=(-1, -1, 0),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        graph.add_edge(
            env_node3,
            env_node3,
            start=env_node3.isite,
            end=env_node3.isite,
            delta=(-1, -1, -1),
            ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))],
        )
        cc = ConnectedComponent(graph=graph)
        assert not cc.is_0d
        assert not cc.is_1d
        assert not cc.is_2d
        assert cc.is_3d
        assert cc.is_periodic
        assert cc.periodicity == "3D"
        assert_allclose(cc.periodicity_vectors, [[0, 1, 0], [1, 1, 0], [1, 1, 1]])
        assert isinstance(cc.periodicity_vectors, list)
        assert cc.periodicity_vectors[0].dtype is np.dtype(np.int64)

    def test_real_systems(self):
        # Initialize geometry and connectivity finders
        strategy = SimplestChemenvStrategy()
        lgf = LocalGeometryFinder()
        cf = ConnectivityFinder()

        # Connectivity of LiFePO4
        struct = self.get_structure("LiFePO4")
        lgf.setup_structure(structure=struct)
        se = lgf.compute_structure_environments(only_atoms=["Li", "Fe", "P"], maximum_distance_factor=1.2)
        lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)
        # Make sure the initial structure and environments are correct
        for site_idx in range(4):
            assert lse.structure[site_idx].specie.symbol == "Li"
            assert lse.coordination_environments[site_idx][0]["ce_symbol"] == "O:6"
        for site_idx in range(4, 8):
            assert lse.structure[site_idx].specie.symbol == "Fe"
            assert lse.coordination_environments[site_idx][0]["ce_symbol"] == "O:6"
        for site_idx in range(8, 12):
            assert lse.structure[site_idx].specie.symbol == "P"
            assert lse.coordination_environments[site_idx][0]["ce_symbol"] == "T:4"
        # Get the connectivity including all environments and check results
        sc = cf.get_structure_connectivity(lse)
        assert len(sc.environment_subgraphs) == 0  # Connected component not computed by default
        ccs = sc.get_connected_components()  # by default, will use all the environments (O:6 and T:4 here)
        assert list(sc.environment_subgraphs) == ["O:6-T:4"]  # Now the default components are there
        assert len(sc.environment_subgraphs) == 1
        assert len(ccs) == 1
        cc = ccs[0]
        assert len(cc) == 12
        assert cc.periodicity == "3D"
        expected_desc = """Connected component with environment nodes :
Node #0 Li (O:6)
Node #1 Li (O:6)
Node #2 Li (O:6)
Node #3 Li (O:6)
Node #4 Fe (O:6)
Node #5 Fe (O:6)
Node #6 Fe (O:6)
Node #7 Fe (O:6)
Node #8 P (T:4)
Node #9 P (T:4)
Node #10 P (T:4)
Node #11 P (T:4)"""
        assert cc.description() == expected_desc

        expected_full_desc = """Connected component with environment nodes :
Node #0 Li (O:6), connected to :
  - Node #1 Li (O:6) with delta image cells
     (-1 0 1)
     (-1 1 1)
  - Node #4 Fe (O:6) with delta image cells
     (-1 1 1)
     (0 1 1)
  - Node #5 Fe (O:6) with delta image cells
     (0 0 1)
  - Node #6 Fe (O:6) with delta image cells
     (-1 1 0)
  - Node #7 Fe (O:6) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #8 P (T:4) with delta image cells
     (-1 0 1)
     (0 0 1)
  - Node #11 P (T:4) with delta image cells
     (-1 1 0)
     (0 1 0)
Node #1 Li (O:6), connected to :
  - Node #0 Li (O:6) with delta image cells
     (1 -1 -1)
     (1 0 -1)
  - Node #4 Fe (O:6) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #5 Fe (O:6) with delta image cells
     (1 0 0)
  - Node #6 Fe (O:6) with delta image cells
     (0 0 -1)
  - Node #7 Fe (O:6) with delta image cells
     (0 0 -1)
     (1 0 -1)
  - Node #8 P (T:4) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #11 P (T:4) with delta image cells
     (0 0 -1)
     (1 0 -1)
Node #2 Li (O:6), connected to :
  - Node #3 Li (O:6) with delta image cells
     (0 0 0)
     (0 1 0)
  - Node #4 Fe (O:6) with delta image cells
     (0 1 0)
  - Node #5 Fe (O:6) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #6 Fe (O:6) with delta image cells
     (-1 1 0)
     (0 1 0)
  - Node #7 Fe (O:6) with delta image cells
     (0 0 0)
  - Node #9 P (T:4) with delta image cells
     (0 1 0)
     (1 1 0)
  - Node #10 P (T:4) with delta image cells
     (-1 0 0)
     (0 0 0)
Node #3 Li (O:6), connected to :
  - Node #2 Li (O:6) with delta image cells
     (0 -1 0)
     (0 0 0)
  - Node #4 Fe (O:6) with delta image cells
     (0 0 0)
  - Node #5 Fe (O:6) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #6 Fe (O:6) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #7 Fe (O:6) with delta image cells
     (0 0 0)
  - Node #9 P (T:4) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #10 P (T:4) with delta image cells
     (-1 0 0)
     (0 0 0)
Node #4 Fe (O:6), connected to :
  - Node #0 Li (O:6) with delta image cells
     (0 -1 -1)
     (1 -1 -1)
  - Node #1 Li (O:6) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #2 Li (O:6) with delta image cells
     (0 -1 0)
  - Node #3 Li (O:6) with delta image cells
     (0 0 0)
  - Node #5 Fe (O:6) with delta image cells
     (0 -1 0)
     (0 0 0)
     (1 -1 0)
     (1 0 0)
  - Node #8 P (T:4) with delta image cells
     (0 -1 0)
     (0 0 0)
  - Node #9 P (T:4) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #11 P (T:4) with delta image cells
     (0 0 -1)
Node #5 Fe (O:6), connected to :
  - Node #0 Li (O:6) with delta image cells
     (0 0 -1)
  - Node #1 Li (O:6) with delta image cells
     (-1 0 0)
  - Node #2 Li (O:6) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #3 Li (O:6) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #4 Fe (O:6) with delta image cells
     (-1 0 0)
     (-1 1 0)
     (0 0 0)
     (0 1 0)
  - Node #8 P (T:4) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #9 P (T:4) with delta image cells
     (0 0 0)
     (0 1 0)
  - Node #10 P (T:4) with delta image cells
     (-1 0 0)
Node #6 Fe (O:6), connected to :
  - Node #0 Li (O:6) with delta image cells
     (1 -1 0)
  - Node #1 Li (O:6) with delta image cells
     (0 0 1)
  - Node #2 Li (O:6) with delta image cells
     (0 -1 0)
     (1 -1 0)
  - Node #3 Li (O:6) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #7 Fe (O:6) with delta image cells
     (0 -1 0)
     (0 0 0)
     (1 -1 0)
     (1 0 0)
  - Node #9 P (T:4) with delta image cells
     (1 0 0)
  - Node #10 P (T:4) with delta image cells
     (0 -1 0)
     (0 0 0)
  - Node #11 P (T:4) with delta image cells
     (0 0 0)
     (1 0 0)
Node #7 Fe (O:6), connected to :
  - Node #0 Li (O:6) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #1 Li (O:6) with delta image cells
     (-1 0 1)
     (0 0 1)
  - Node #2 Li (O:6) with delta image cells
     (0 0 0)
  - Node #3 Li (O:6) with delta image cells
     (0 0 0)
  - Node #6 Fe (O:6) with delta image cells
     (-1 0 0)
     (-1 1 0)
     (0 0 0)
     (0 1 0)
  - Node #8 P (T:4) with delta image cells
     (0 0 1)
  - Node #10 P (T:4) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #11 P (T:4) with delta image cells
     (0 0 0)
     (0 1 0)
Node #8 P (T:4), connected to :
  - Node #0 Li (O:6) with delta image cells
     (0 0 -1)
     (1 0 -1)
  - Node #1 Li (O:6) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #4 Fe (O:6) with delta image cells
     (0 0 0)
     (0 1 0)
  - Node #5 Fe (O:6) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #7 Fe (O:6) with delta image cells
     (0 0 -1)
Node #9 P (T:4), connected to :
  - Node #2 Li (O:6) with delta image cells
     (-1 -1 0)
     (0 -1 0)
  - Node #3 Li (O:6) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #4 Fe (O:6) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #5 Fe (O:6) with delta image cells
     (0 -1 0)
     (0 0 0)
  - Node #6 Fe (O:6) with delta image cells
     (-1 0 0)
Node #10 P (T:4), connected to :
  - Node #2 Li (O:6) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #3 Li (O:6) with delta image cells
     (0 0 0)
     (1 0 0)
  - Node #5 Fe (O:6) with delta image cells
     (1 0 0)
  - Node #6 Fe (O:6) with delta image cells
     (0 0 0)
     (0 1 0)
  - Node #7 Fe (O:6) with delta image cells
     (0 0 0)
     (1 0 0)
Node #11 P (T:4), connected to :
  - Node #0 Li (O:6) with delta image cells
     (0 -1 0)
     (1 -1 0)
  - Node #1 Li (O:6) with delta image cells
     (-1 0 1)
     (0 0 1)
  - Node #4 Fe (O:6) with delta image cells
     (0 0 1)
  - Node #6 Fe (O:6) with delta image cells
     (-1 0 0)
     (0 0 0)
  - Node #7 Fe (O:6) with delta image cells
     (0 -1 0)
     (0 0 0)"""
        assert cc.description(full=True) == expected_full_desc
        # Get the connectivity for T:4 and O:6 separately and check results
        # Only tetrahedral
        sc.setup_environment_subgraph(environments_symbols=["T:4"])
        assert list(sc.environment_subgraphs) == ["O:6-T:4", "T:4"]
        ccs = sc.get_connected_components()
        assert len(ccs) == 4
        for cc in ccs:
            assert cc.periodicity == "0D"
        # Only octahedral
        sc.setup_environment_subgraph(environments_symbols=["O:6"])
        assert list(sc.environment_subgraphs) == ["O:6-T:4", "T:4", "O:6"]
        ccs = sc.get_connected_components()
        assert len(ccs) == 1
        cc = ccs[0]
        assert cc.periodicity == "3D"
        # Only Manganese octahedral
        sc.setup_environment_subgraph(environments_symbols=["O:6"], only_atoms=["Fe"])
        assert list(sc.environment_subgraphs) == [
            "O:6-T:4",
            "T:4",
            "O:6",
            "O:6#Fe",
        ]
        ccs = sc.get_connected_components()
        assert len(ccs) == 2
        for cc in ccs:
            assert cc.periodicity == "2D"

        # Connectivity of Li4Fe3Mn1(PO4)4
        struct = Structure.from_file(f"{TEST_FILES_DIR}/cif/Li4Fe3Mn1(PO4)4.cif")
        lgf.setup_structure(structure=struct)
        se = lgf.compute_structure_environments(only_atoms=["Li", "Fe", "Mn", "P"], maximum_distance_factor=1.2)
        lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)
        # Make sure the initial structure and environments are correct
        for site_idx in range(4):
            assert lse.structure[site_idx].specie.symbol == "Li"
            assert lse.coordination_environments[site_idx][0]["ce_symbol"] == "O:6"
        for site_idx in range(4, 5):
            assert lse.structure[site_idx].specie.symbol == "Mn"
            assert lse.coordination_environments[site_idx][0]["ce_symbol"] == "O:6"
        for site_idx in range(5, 8):
            assert lse.structure[site_idx].specie.symbol == "Fe"
            assert lse.coordination_environments[site_idx][0]["ce_symbol"] == "O:6"
        for site_idx in range(8, 12):
            assert lse.structure[site_idx].specie.symbol == "P"
            assert lse.coordination_environments[site_idx][0]["ce_symbol"] == "T:4"
        # Get the connectivity including all environments and check results
        sc = cf.get_structure_connectivity(lse)
        assert len(sc.environment_subgraphs) == 0  # Connected component not computed by default
        ccs = sc.get_connected_components()  # by default, will use all the environments (O:6 and T:4 here)
        # Now connected components for the defaults are there :
        assert list(sc.environment_subgraphs) == ["O:6-T:4"]
        assert len(sc.environment_subgraphs) == 1
        assert len(ccs) == 1
        cc = ccs[0]
        assert cc.periodicity == "3D"
        # Get the connectivity for Li octahedral only
        ccs = sc.get_connected_components(environments_symbols=["O:6"], only_atoms=["Li"])
        assert list(sc.environment_subgraphs) == ["O:6-T:4", "O:6#Li"]
        assert len(ccs) == 2
        for cc in ccs:
            assert cc.periodicity == "1D"
            assert len(cc) == 2
        # Sort connected components as they might
        # come in a different order depending on
        # the algorithm used to get them.
        sorted_ccs = sorted(ccs, key=lambda x: min(x.graph.nodes()))
        expected_full_desc_0 = """Connected component with environment nodes :
Node #0 Li (O:6), connected to :
  - Node #1 Li (O:6) with delta image cells
     (1 -1 1)
     (1 0 1)
Node #1 Li (O:6), connected to :
  - Node #0 Li (O:6) with delta image cells
     (-1 0 -1)
     (-1 1 -1)"""
        assert sorted_ccs[0].description(full=True) == expected_full_desc_0

        expected_full_desc_1 = """Connected component with environment nodes :
Node #2 Li (O:6), connected to :
  - Node #3 Li (O:6) with delta image cells
     (0 -1 0)
     (0 0 0)
Node #3 Li (O:6), connected to :
  - Node #2 Li (O:6) with delta image cells
     (0 0 0)
     (0 1 0)"""
        assert sorted_ccs[1].description(full=True) == expected_full_desc_1
        # Get the connectivity for Mn octahedral only
        ccs = sc.get_connected_components(environments_symbols=["O:6"], only_atoms=["Mn"])
        assert list(sc.environment_subgraphs) == ["O:6-T:4", "O:6#Li", "O:6#Mn"]
        assert len(ccs) == 1
        assert ccs[0].periodicity == "0D"
        # Get the connectivity for Fe octahedral only
        ccs = sc.get_connected_components(environments_symbols=["O:6"], only_atoms=["Fe"])
        assert list(sc.environment_subgraphs) == [
            "O:6-T:4",
            "O:6#Li",
            "O:6#Mn",
            "O:6#Fe",
        ]
        assert len(ccs) == 2
        ccs_periodicities = {cc.periodicity for cc in ccs}
        assert ccs_periodicities == {"0D", "2D"}

    def test_coordination_sequences(self):
        BaTiO3_se_fpath = f"{TEST_FILES_DIR}/analysis/chemenv/structure_environments/se_mp-5020.json"
        with open(BaTiO3_se_fpath, "rb") as file:
            dct = orjson.loads(file.read())
        struct_envs = StructureEnvironments.from_dict(dct)
        lse = LightStructureEnvironments.from_structure_environments(
            strategy=SimplestChemenvStrategy(), structure_environments=struct_envs
        )
        cf = ConnectivityFinder()
        sc = cf.get_structure_connectivity(light_structure_environments=lse)
        ccs_oct = sc.get_connected_components(environments_symbols=["O:6"])
        ccs_all = sc.get_connected_components(environments_symbols=["O:6", "C:12"])
        assert len(ccs_oct) == 1
        assert len(ccs_all) == 1
        cc_oct = ccs_oct[0]
        cc_all = ccs_all[0]
        cc_oct_node = next(iter(cc_oct.graph.nodes()))
        cseq = cc_oct.coordination_sequence(source_node=cc_oct_node, path_size=6)
        assert cseq == {1: 6, 2: 18, 3: 38, 4: 66, 5: 102, 6: 146}
        cc_all_oct_node = next(n for n in cc_all.graph.nodes() if n.coordination_environment == "O:6")
        cc_all_cuboct_node = next(n for n in cc_all.graph.nodes() if n.coordination_environment == "C:12")
        cseq = cc_all.coordination_sequence(source_node=cc_all_oct_node, path_size=6)
        assert cseq == {1: 14, 2: 74, 3: 218, 4: 442, 5: 746, 6: 1130}
        cseq = cc_all.coordination_sequence(source_node=cc_all_cuboct_node, path_size=6)
        assert cseq == {1: 26, 2: 122, 3: 298, 4: 554, 5: 890, 6: 1306}
        cseq = cc_all.coordination_sequence(source_node=cc_all_cuboct_node, path_size=6, include_source=True)
        assert cseq == {0: 1, 1: 26, 2: 122, 3: 298, 4: 554, 5: 890, 6: 1306}
        cseq = cc_all.coordination_sequence(source_node=cc_all_oct_node, path_size=4, coordination="env:number")
        assert cseq == {
            1: {"O:6": 6, "C:12": 8},
            2: {"O:6": 26, "C:12": 48},
            3: {"O:6": 90, "C:12": 128},
            4: {"O:6": 194, "C:12": 248},
        }
        cseq = cc_all.coordination_sequence(source_node=cc_all_cuboct_node, path_size=4, coordination="env:number")
        assert cseq == {
            1: {"O:6": 8, "C:12": 18},
            2: {"O:6": 48, "C:12": 74},
            3: {"O:6": 128, "C:12": 170},
            4: {"O:6": 248, "C:12": 306},
        }
        cseq = cc_all.coordination_sequence(
            source_node=cc_all_cuboct_node,
            path_size=4,
            coordination="env:number",
            include_source=True,
        )
        assert cseq == {
            0: {"C:12": 1},
            1: {"O:6": 8, "C:12": 18},
            2: {"O:6": 48, "C:12": 74},
            3: {"O:6": 128, "C:12": 170},
            4: {"O:6": 248, "C:12": 306},
        }
