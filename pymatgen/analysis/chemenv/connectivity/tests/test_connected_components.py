#!/usr/bin/env python


__author__ = 'waroquiers'

import networkx as nx
from pymatgen.analysis.chemenv.connectivity.connected_components import ConnectedComponent
from pymatgen.analysis.chemenv.connectivity.environment_nodes import EnvironmentNode
from pymatgen.core.sites import PeriodicSite
from pymatgen.core.lattice import Lattice
from pymatgen.util.testing import PymatgenTest

import bson
import json
import pytest
import numpy as np
import copy


class ConnectedComponentTest(PymatgenTest):

    def test_init(self):
        # Generic connected component not using EnvironmentNodes
        # (as_dict won't work on such a ConnectedComponent instance)
        cc = ConnectedComponent(environments=['a', 'b', 'c', 'd', 'e', 'f'],
                                links=[('a', 'b', 0), ('a', 'c', 0), ('b', 'c', 0),
                                       ('a', 'a', 0), ('a', 'b', 1),
                                       ('c', 'd'), ('c', 'd'), ('c', 'd'), ('d', 'e'), ('e', 'f')],
                                environments_data={'a': {'var1': 2, 'var2': 3},
                                                   'b': {'var1': 3}},
                                links_data={('c', 'b'): {'bcinfo': 2},
                                            ('a', 'b', 1): {'ab1info': 4},
                                            ('d', 'c'): {'dcinfo': 8}})
        assert isinstance(cc.graph, nx.MultiGraph)
        nodes = list(cc.graph.nodes())
        assert set(nodes) == {'a', 'b', 'c', 'd', 'e', 'f'}
        edges = cc.graph.edges()
        assert len(edges) == 10
        assert cc.graph['a']['b'] == {0: {}, 1: {'ab1info': 4}}
        assert cc.graph['b']['a'] == {0: {}, 1: {'ab1info': 4}}
        assert cc.graph['c']['b'] == {0: {'bcinfo': 2}}
        assert cc.graph['b']['c'] == {0: {'bcinfo': 2}}
        assert len(cc.graph['c']['d']) == 3
        assert cc.graph['c']['d'][0] == {'dcinfo': 8}
        assert cc.graph['c']['d'][1] == {'dcinfo': 8}
        assert cc.graph['c']['d'][2] == {'dcinfo': 8}

        assert cc.graph.nodes(data=True)['a'] == {'var1': 2, 'var2': 3}
        assert cc.graph.nodes(data=True)['b'] == {'var1': 3}
        assert cc.graph.nodes(data=True)['c'] == {}

        # Make a deep copy of the graph to actually check that the objects are different instances
        mygraph = copy.deepcopy(cc.graph)
        assert isinstance(mygraph, nx.MultiGraph)  # Check that it is indeed the same type of graph

        cc2 = ConnectedComponent(graph=mygraph)
        assert set(list(mygraph.nodes())) == set(list(cc2.graph.nodes()))
        assert set(list(mygraph.edges())) == set(list(cc2.graph.edges()))
        assert len(cc2.graph) == 6

    def test_serialization(self):
        lat = Lattice.hexagonal(a=2.0, c=2.5)
        en1 = EnvironmentNode(central_site=PeriodicSite('Si',
                                                        coords=np.array([0.0, 0.0, 0.0]),
                                                        lattice=lat),
                              i_central_site=3, ce_symbol='T:4')
        en2 = EnvironmentNode(central_site=PeriodicSite('Ag',
                                                        coords=np.array([0.0, 0.0, 0.5]),
                                                        lattice=lat),
                              i_central_site=5, ce_symbol='T:4')
        en3 = EnvironmentNode(central_site=PeriodicSite('Ag',
                                                        coords=np.array([0.0, 0.5, 0.5]),
                                                        lattice=lat),
                              i_central_site=8, ce_symbol='O:6')

        graph = nx.MultiGraph()
        graph.add_nodes_from([en1, en2, en3])

        graph.add_edge(en1, en2, start=en1.isite, end=en2.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en1, en3, start=en1.isite, end=en2.isite,
                       delta=(0, 0, 0), ligands=[(10, (0, 0, 1), (0, 0, 1)), (11, (0, 0, 1), (0, 0, 1))])

        cc = ConnectedComponent(graph=graph)
        ref_sorted_edges = [[en1, en2], [en1, en3]]
        sorted_edges = sorted([sorted(e) for e in cc.graph.edges()])
        assert sorted_edges == ref_sorted_edges

        ccfromdict = ConnectedComponent.from_dict(cc.as_dict())
        ccfromjson = ConnectedComponent.from_dict(json.loads(json.dumps(cc.as_dict())))
        bson_data = bson.BSON.encode(cc.as_dict())
        ccfrombson = ConnectedComponent.from_dict(bson_data.decode())
        for loaded_cc in [ccfromdict, ccfromjson, ccfrombson]:
            assert loaded_cc.graph.number_of_nodes() == 3
            assert loaded_cc.graph.number_of_edges() == 2
            assert set(list(cc.graph.nodes())) == set(list(loaded_cc.graph.nodes()))
            assert sorted_edges == sorted([sorted(e) for e in loaded_cc.graph.edges()])

            for ii, e in enumerate(sorted_edges):
                assert cc.graph[e[0]][e[1]] == loaded_cc.graph[e[0]][e[1]]

            for node in loaded_cc.graph.nodes():
                assert isinstance(node.central_site, PeriodicSite)

    def test_periodicity(self):
        en1 = EnvironmentNode(central_site='Si',
                              i_central_site=3, ce_symbol='T:4')
        en2 = EnvironmentNode(central_site='Ag',
                              i_central_site=5, ce_symbol='T:4')
        en3 = EnvironmentNode(central_site='Ag',
                              i_central_site=8, ce_symbol='O:6')
        en4 = EnvironmentNode(central_site='Fe',
                              i_central_site=23, ce_symbol='C:8')

        graph = nx.MultiGraph()
        graph.add_nodes_from([en1, en2, en3])
        graph.add_edge(en1, en2, start=en1.isite, end=en2.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en1, en3, start=en1.isite, end=en2.isite,
                       delta=(0, 0, 0), ligands=[(10, (0, 0, 1), (0, 0, 1)), (11, (0, 0, 1), (0, 0, 1))])
        cc = ConnectedComponent(graph=graph)
        assert cc.is_0d
        assert not cc.is_1d
        assert not cc.is_2d
        assert not cc.is_3d
        assert not cc.is_periodic
        assert cc.periodicity == '0D'

        graph = nx.MultiGraph()
        graph.add_nodes_from([en1, en2, en3])
        graph.add_edge(en1, en2, start=en1.isite, end=en2.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en1, en3, start=en1.isite, end=en3.isite,
                       delta=(0, 0, 0), ligands=[(10, (0, 0, 1), (0, 0, 1)), (11, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en2, en3, start=en2.isite, end=en3.isite,
                       delta=(0, 0, 1), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        cc = ConnectedComponent(graph=graph)
        assert not cc.is_0d
        assert cc.is_1d
        assert not cc.is_2d
        assert not cc.is_3d
        assert cc.is_periodic
        assert cc.periodicity == '1D'

        graph = nx.MultiGraph()
        graph.add_nodes_from([en1, en2, en3])
        graph.add_edge(en1, en2, start=en1.isite, end=en2.isite,
                       delta=(0, 0, 1), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en1, en3, start=en1.isite, end=en3.isite,
                       delta=(0, 0, 0), ligands=[(10, (0, 0, 1), (0, 0, 1)), (11, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en2, en3, start=en2.isite, end=en3.isite,
                       delta=(0, 0, -1), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        cc = ConnectedComponent(graph=graph)
        assert cc.periodicity == '0D'

        # Test errors when computing periodicity
        graph = nx.MultiGraph()
        graph.add_nodes_from([en1, en2, en3])
        graph.add_edge(en1, en1, start=en1.isite, end=en1.isite,
                       delta=(0, 0, 1), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en1, en1, start=en1.isite, end=en1.isite,
                       delta=(0, 0, 1), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        cc = ConnectedComponent(graph=graph)
        with pytest.raises(ValueError, match=r'There should not be self loops with the same '
                                             r'\x28or opposite\x29 delta image\x2E'):
            cc.compute_periodicity_all_simple_paths_algorithm()

        graph = nx.MultiGraph()
        graph.add_nodes_from([en1, en2, en3])
        graph.add_edge(en1, en1, start=en1.isite, end=en1.isite,
                       delta=(3, 2, -1), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en1, en1, start=en1.isite, end=en1.isite,
                       delta=(-3, -2, 1), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        cc = ConnectedComponent(graph=graph)
        with pytest.raises(ValueError, match=r'There should not be self loops with the same '
                                             r'\x28or opposite\x29 delta image\x2E'):
            cc.compute_periodicity_all_simple_paths_algorithm()

        graph = nx.MultiGraph()
        graph.add_nodes_from([en1, en2, en3])
        graph.add_edge(en1, en1, start=en1.isite, end=en1.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        cc = ConnectedComponent(graph=graph)
        with pytest.raises(ValueError, match=r'There should not be self loops with delta image = '
                                             r'\x280, 0, 0\x29\x2E'):
            cc.compute_periodicity_all_simple_paths_algorithm()

        # Test a 2d periodicity
        graph = nx.MultiGraph()
        graph.add_nodes_from([en1, en2, en3, en4])
        graph.add_edge(en1, en2, start=en1.isite, end=en2.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en1, en3, start=en1.isite, end=en3.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en4, en2, start=en4.isite, end=en2.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en3, en4, start=en4.isite, end=en3.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en3, en4, start=en4.isite, end=en3.isite,
                       delta=(0, -1, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en3, en2, start=en2.isite, end=en3.isite,
                       delta=(-1, -1, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        cc = ConnectedComponent(graph=graph)
        assert not cc.is_0d
        assert not cc.is_1d
        assert cc.is_2d
        assert not cc.is_3d
        assert cc.is_periodic
        assert cc.periodicity == '2D'
        assert np.allclose(cc.periodicity_vectors, [np.array([0, 1, 0]), np.array([1, 1, 0])])
        assert type(cc.periodicity_vectors) is list
        assert cc.periodicity_vectors[0].dtype is np.dtype(int)

        # Test a 3d periodicity
        graph = nx.MultiGraph()
        graph.add_nodes_from([en1, en2, en3, en4])
        graph.add_edge(en1, en2, start=en1.isite, end=en2.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en1, en3, start=en1.isite, end=en3.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en4, en2, start=en4.isite, end=en2.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en3, en4, start=en4.isite, end=en3.isite,
                       delta=(0, 0, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en3, en4, start=en4.isite, end=en3.isite,
                       delta=(0, -1, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en3, en2, start=en2.isite, end=en3.isite,
                       delta=(-1, -1, 0), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        graph.add_edge(en3, en3, start=en3.isite, end=en3.isite,
                       delta=(-1, -1, -1), ligands=[(2, (0, 0, 1), (0, 0, 1)), (1, (0, 0, 1), (0, 0, 1))])
        cc = ConnectedComponent(graph=graph)
        assert not cc.is_0d
        assert not cc.is_1d
        assert not cc.is_2d
        assert cc.is_3d
        assert cc.is_periodic
        assert cc.periodicity == '3D'
        assert np.allclose(cc.periodicity_vectors, [np.array([0, 1, 0]), np.array([1, 1, 0]), np.array([1, 1, 1])])
        assert type(cc.periodicity_vectors) is list
        assert cc.periodicity_vectors[0].dtype is np.dtype(int)
