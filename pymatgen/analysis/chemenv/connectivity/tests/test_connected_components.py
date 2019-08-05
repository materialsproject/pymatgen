#!/usr/bin/env python


__author__ = 'waroquiers'

# import unittest2 as unittest
import networkx as nx
from pymatgen.analysis.chemenv.connectivity.connected_components import ConnectedComponent
from pymatgen.analysis.chemenv.connectivity.environment_nodes import get_environment_node, EnvironmentNode
from pymatgen.util.testing import PymatgenTest


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
        self.assertIsInstance(cc.graph, nx.MultiGraph)
        nodes = list(cc.graph.nodes())
        self.assertEqual(set(nodes), {'a', 'b', 'c', 'd', 'e', 'f'})
        edges = cc.graph.edges()
        self.assertEqual(len(edges), 10)
        self.assertEqual(cc.graph['a']['b'], {0: {}, 1: {'ab1info': 4}})
        self.assertEqual(cc.graph['b']['a'], {0: {}, 1: {'ab1info': 4}})
        self.assertEqual(cc.graph['c']['b'], {0: {'bcinfo': 2}})
        self.assertEqual(cc.graph['b']['c'], {0: {'bcinfo': 2}})
        self.assertEqual(len(cc.graph['c']['d']), 3)
        self.assertEqual(cc.graph['c']['d'][0], {'dcinfo': 8})
        self.assertEqual(cc.graph['c']['d'][1], {'dcinfo': 8})
        self.assertEqual(cc.graph['c']['d'][2], {'dcinfo': 8})

        self.assertEqual(cc.graph.nodes(data=True)['a'], {'var1': 2, 'var2': 3})
        self.assertEqual(cc.graph.nodes(data=True)['b'], {'var1': 3})
        self.assertEqual(cc.graph.nodes(data=True)['c'], {})

        import copy
        mygraph = copy.deepcopy(cc.graph)
        self.assertIsInstance(mygraph, nx.MultiGraph)  # Check that it is indeed the same type of graph

        cc2 = ConnectedComponent(graph=mygraph)
        self.assertEqual(set(list(mygraph.nodes())), set(list(cc2.graph.nodes())))
        self.assertEqual(set(list(mygraph.edges())), set(list(cc2.graph.edges())))
        self.assertEqual(cc2.graph, mygraph)

    def test_as_dict(self):
        pass
        # cc = ConnectedComponent()
        # cc.as_dict()
        # pass


if __name__ == "__main__":
    import unittest
    unittest.main()
