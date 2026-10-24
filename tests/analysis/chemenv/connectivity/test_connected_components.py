from __future__ import annotations

import copy

import matplotlib as mpl
import matplotlib.pyplot as plt
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
from pymatgen.analysis.chemenv.utils.chemenv_errors import ChemenvError
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

__author__ = "waroquiers"

mpl.use("Agg")


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


SE_DIR = f"{TEST_FILES_DIR}/analysis/chemenv/structure_environments"


def _structure_connectivity(fname: str):
    with open(f"{SE_DIR}/{fname}", "rb") as file:
        struct_envs = StructureEnvironments.from_dict(orjson.loads(file.read()))
    lse = LightStructureEnvironments.from_structure_environments(
        strategy=SimplestChemenvStrategy(), structure_environments=struct_envs
    )
    return ConnectivityFinder().get_structure_connectivity(light_structure_environments=lse)


@pytest.fixture(scope="module")
def BaTiO3_sc():
    return _structure_connectivity("se_mp-5020.json")


@pytest.fixture(scope="module")
def K2NaNb2Fe7Si8H4O31_sc():
    return _structure_connectivity("se_mp-743972.json")


def _chain_graph():
    """Infinite 1D chain A-B-A-B... along a: A-B in the cell, B-A across the cell boundary."""
    lattice = Lattice.cubic(4.0)
    node_a = EnvironmentNode(PeriodicSite("Si", [0.0, 0.0, 0.0], lattice), i_central_site=0, ce_symbol="T:4")
    node_b = EnvironmentNode(PeriodicSite("Si", [0.5, 0.0, 0.0], lattice), i_central_site=1, ce_symbol="T:4")
    graph = nx.MultiGraph()
    graph.add_nodes_from([node_a, node_b])
    graph.add_edge(node_a, node_b, start=0, end=1, delta=(0, 0, 0), ligands=[])
    graph.add_edge(node_a, node_b, start=1, end=0, delta=(1, 0, 0), ligands=[])
    return graph, node_a, node_b


class TestConnectedComponentMethods(MatSciTest):
    def test_init_errors_and_links_data(self):
        with pytest.raises(ChemenvError, match="Trying to add edge with some unexistent node"):
            ConnectedComponent(environments=["a", "b"], links=[("a", "z")])

        cc = ConnectedComponent(
            environments=["a", "b", "c"],
            links=[("a", "b", 0), ("b", "c", 0), ("a", "c")],
            links_data={("b", "a", 0): {"info": 1}, ("c", "b"): {"info": 2}},
        )
        assert sorted(cc.graph.nodes()) == ["a", "b", "c"]
        assert cc.graph["a"]["b"][0] == {"info": 1}
        assert cc.graph["b"]["c"][0] == {"info": 2}
        assert cc.graph["a"]["c"][0] == {}

    def test_edgedictkey_to_edgekey(self):
        assert ConnectedComponent._edgedictkey_to_edgekey(4) == 4
        assert ConnectedComponent._edgedictkey_to_edgekey("4") == 4
        assert ConnectedComponent._edgedictkey_to_edgekey("some-key") == "some-key"
        with pytest.raises(TypeError, match="Edge key in a dict of dicts representation"):
            ConnectedComponent._edgedictkey_to_edgekey(0.5)

    def test_order_periodicity_vectors(self):
        vectors = [np.array([1, 1, 0]), np.array([0, -1, 1]), np.array([0, 1, 1])]
        ordered = ConnectedComponent._order_vectors(vectors)
        assert [vec.tolist() for vec in ordered] == [[0, 1, -1], [0, 1, 1], [1, 1, 0]]

        cc = ConnectedComponent(graph=_chain_graph()[0])
        cc._periodicity_vectors = [np.array(vec) for vec in np.eye(4, 3, dtype=int)]
        with pytest.raises(ValueError, match="Number of periodicity vectors is larger than 3"):
            cc._order_periodicity_vectors()

    def test_coordination_sequence_errors(self):
        graph, node_a, _node_b = _chain_graph()
        cc = ConnectedComponent(graph=graph)
        # Infinite chain: always 2 new neighbors per shell
        assert cc.coordination_sequence(node_a, path_size=3) == {1: 2, 2: 2, 3: 2}
        assert cc.coordination_sequence(node_a, path_size=2, coordination="env:number", include_source=True) == {
            0: {"T:4": 1},
            1: {"T:4": 2},
            2: {"T:4": 2},
        }
        other_node = EnvironmentNode(central_site="Fe", i_central_site=7, ce_symbol="O:6")
        with pytest.raises(ValueError, match="Node not in Connected Component"):
            cc.coordination_sequence(other_node)
        for include_source in (True, False):
            with pytest.raises(ValueError, match="Coordination type 'wrong' is not valid"):
                cc.coordination_sequence(node_a, coordination="wrong", include_source=include_source)

    def test_compute_periodicity_invalid_algorithm(self):
        cc = ConnectedComponent(graph=_chain_graph()[0])
        with pytest.raises(ValueError, match="algorithm='wrong' is not allowed to compute periodicity"):
            cc.compute_periodicity(algorithm="wrong")

    @pytest.mark.parametrize(
        ("fname", "envs", "periodicity", "pvecs_simple_paths", "pvecs_cycle_basis"),
        [
            ("se_mp-5020.json", ["O:6"], "3D", [[0, 0, 1], [0, 1, 0], [1, 0, 0]], [[0, 0, 1], [0, 1, 0], [1, 0, 0]]),
            (
                "se_mp-5020.json",
                ["O:6", "C:12"],
                "3D",
                [[0, 0, 1], [0, 1, 0], [1, 0, 0]],
                [[0, 0, 1], [0, 1, 0], [1, 0, 0]],
            ),
            (
                "se_mp-7000.json",
                ["T:4"],
                "3D",
                [[0, 0, 1], [0, 1, -1], [1, 0, 0]],
                [[0, 1, -1], [1, 0, 1], [1, 1, -1]],
            ),
            ("se_mp-743972.json", ["BS_1:10"], "1D", [[1, 0, 0]], [[1, 0, 0]]),
            ("se_mp-743972.json", ["T:4"], "1D", [[1, 0, 0]], [[1, 0, 0]]),
        ],
    )
    def test_compute_periodicity_algorithms(self, fname, envs, periodicity, pvecs_simple_paths, pvecs_cycle_basis):
        sc = _structure_connectivity(fname)
        ccs = sc.get_connected_components(environments_symbols=envs)
        for cc in ccs:
            cc.compute_periodicity(algorithm="all_simple_paths")
            assert cc.periodicity == periodicity
            assert [vec.tolist() for vec in cc.periodicity_vectors] == pvecs_simple_paths
            cc_cb = ConnectedComponent(graph=cc.graph)
            cc_cb.compute_periodicity(algorithm="cycle_basis")
            assert cc_cb.periodicity == periodicity
            assert [vec.tolist() for vec in cc_cb.periodicity_vectors] == pvecs_cycle_basis

    def test_compute_periodicity_cycle_basis_parallel_edges(self):
        # Two nodes linked by parallel edges only (no cycle in the simple graph)
        graph, node_a, node_b = _chain_graph()
        cc = ConnectedComponent(graph=graph)
        cc.compute_periodicity(algorithm="cycle_basis")
        assert cc.periodicity == "1D"
        assert [vec.tolist() for vec in cc.periodicity_vectors] == [[1, 0, 0]]

        # Add more parallel edges along b and c to get a 3D net
        graph.add_edge(node_a, node_b, start=0, end=1, delta=(0, 1, 0), ligands=[])
        graph.add_edge(node_a, node_b, start=0, end=1, delta=(0, 0, 1), ligands=[])
        cc = ConnectedComponent(graph=graph)
        cc.compute_periodicity(algorithm="cycle_basis")
        assert cc.is_3d
        assert [vec.tolist() for vec in cc.periodicity_vectors] == [[0, 0, 1], [0, 1, 0], [1, 0, 0]]

    def test_elastic_centered_graph(self, K2NaNb2Fe7Si8H4O31_sc):
        ccs = K2NaNb2Fe7Si8H4O31_sc.get_connected_components(environments_symbols=["T:4"])
        assert len(ccs) == 2
        for cc in ccs:
            original_deltas = sorted(tuple(data["delta"]) for *_, data in cc.graph.edges(data=True))
            centered = cc.elastic_centered_graph()
            assert set(centered.nodes()) == set(cc.graph.nodes())
            assert centered.number_of_edges() == cc.graph.number_of_edges() == 4
            deltas = sorted(tuple(int(abs(dd)) for dd in data["delta"]) for *_, data in centered.edges(data=True))
            # 4 tetrahedra forming a chain: 3 links inside the cell and one link to the next cell along a
            assert deltas == [(0, 0, 0), (0, 0, 0), (0, 0, 0), (1, 0, 0)]
            # The original graph is not modified
            assert sorted(tuple(data["delta"]) for *_, data in cc.graph.edges(data=True)) == original_deltas
            # In the input graph, more than one link crosses the cell boundary
            assert sum(delta != (0, 0, 0) for delta in original_deltas) > 1

    def test_elastic_centered_graph_errors(self):
        graph, node_a, node_b = _chain_graph()
        graph.add_edge(node_a, node_b, start=0, end=1, delta=(0, 0, 0), ligands=[])
        with pytest.raises(ValueError, match="Should not have more than one 000 delta"):
            ConnectedComponent(graph=graph).elastic_centered_graph()

        graph = nx.MultiGraph()
        graph.add_nodes_from([node_a, node_b])
        with pytest.raises(RuntimeError, match="Could not find a centered graph"):
            ConnectedComponent(graph=graph).elastic_centered_graph()

    @pytest.mark.parametrize(
        ("multiplicity", "expected_edges"),
        [
            (1, [(0, 0, (1, 0, 0))]),
            (2, [(0, 1, (0, 0, 0)), (1, 0, (1, 0, 0))]),
            (3, [(0, 1, (0, 0, 0)), (1, 2, (0, 0, 0)), (2, 0, (1, 0, 0))]),
            ([3], [(0, 1, (0, 0, 0)), (1, 2, (0, 0, 0)), (2, 0, (1, 0, 0))]),
        ],
    )
    def test_make_supergraph(self, K2NaNb2Fe7Si8H4O31_sc, multiplicity, expected_edges):
        # Single BS_1:10 environment linked to its own image along a (1D chain)
        ccs = K2NaNb2Fe7Si8H4O31_sc.get_connected_components(environments_symbols=["BS_1:10"])
        assert len(ccs) == 1
        cc = ccs[0]
        assert cc.periodicity == "1D"
        super_graph = cc.make_supergraph(multiplicity)
        mult = multiplicity if isinstance(multiplicity, int) else multiplicity[0]
        assert super_graph.number_of_nodes() == mult
        edges = sorted(
            (int(data["start"]), int(data["end"]), tuple(int(dd) for dd in data["delta"]))
            for *_, data in super_graph.edges(data=True)
        )
        assert edges == expected_edges
        for *_, data in super_graph.edges(data=True):
            assert len(data["ligands"]) == 1

    @pytest.mark.parametrize("sign", [1, -1])
    @pytest.mark.parametrize("multiplicity", [2, 3])
    def test_make_supergraph_chain(self, sign, multiplicity):
        # A-B chain: the periodic link may be stored with either orientation of its delta image
        graph, node_a, node_b = _chain_graph()
        graph.remove_edge(node_a, node_b, key=1)
        graph.add_edge(node_a, node_b, start=0, end=1, delta=(-sign, 0, 0), ligands=[])
        cc = ConnectedComponent(graph=graph)
        assert [vec.tolist() for vec in cc.periodicity_vectors] == [[1, 0, 0]]
        super_graph = cc.make_supergraph(multiplicity)
        # The supergraph is a closed ring of 2 * multiplicity environments
        assert super_graph.number_of_nodes() == 2 * multiplicity
        assert super_graph.number_of_edges() == 2 * multiplicity
        assert nx.is_connected(super_graph)
        assert all(degree == 2 for _, degree in super_graph.degree())
        n_periodic = sum(any(data["delta"]) for *_, data in super_graph.edges(data=True))
        assert n_periodic == 1

    @pytest.mark.parametrize("prop", ["is_0d", "is_1d", "is_2d", "is_3d", "periodicity_vectors", "periodicity"])
    def test_lazy_periodicity(self, prop):
        cc = ConnectedComponent(graph=_chain_graph()[0])
        assert cc._periodicity_vectors is None
        expected = {
            "is_0d": False,
            "is_1d": True,
            "is_2d": False,
            "is_3d": False,
            "periodicity": "1D",
        }
        value = getattr(cc, prop)
        if prop == "periodicity_vectors":
            assert [vec.tolist() for vec in value] == [[1, 0, 0]]
        else:
            assert value == expected[prop]
        assert cc._periodicity_vectors is not None

    def test_make_supergraph_not_implemented(self, BaTiO3_sc):
        cc = BaTiO3_sc.get_connected_components(environments_symbols=["O:6"])[0]
        with pytest.raises(NotImplementedError, match="make_supergraph not yet implemented"):
            cc.make_supergraph([2, 2])

    @pytest.mark.parametrize(
        ("fname", "envs"),
        [("se_mp-5020.json", ["O:6", "C:12"]), ("se_mp-743972.json", ["T:4"]), ("se_mp-743972.json", ["BS_1:10"])],
    )
    def test_show_graph(self, fname, envs, tmp_path):
        cc = _structure_connectivity(fname).get_connected_components(environments_symbols=envs)[0]
        _ = cc.periodicity  # Periodicity vectors are used to highlight periodic links in 1D components
        graph = copy.deepcopy(cc.graph)
        save_file = tmp_path / "graph.png"
        with np.errstate(invalid="ignore"):  # Self-loops have zero-length arcs
            cc.show_graph(graph=graph, save_file=str(save_file))
        assert save_file.is_file()
        ax = plt.gca()
        # One circle per node and one arrow per edge
        assert len(ax.patches) == graph.number_of_nodes() + graph.number_of_edges()
        assert all("patch" in data for _, data in graph.nodes(data=True))
        plt.close("all")

        cc.show_graph(graph=nx.MultiGraph(cc.graph), drawing_type="draw_random")
        assert len(plt.gca().collections) > 0
        plt.close("all")
