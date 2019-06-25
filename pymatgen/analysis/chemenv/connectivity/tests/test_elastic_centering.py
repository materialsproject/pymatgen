#!/usr/bin/env python


__author__ = 'waroquiers'

# import unittest2 as unittest
import networkx as nx
import os
import shutil
from pymatgen.analysis.chemenv.connectivity.environment_nodes import get_environment_node
from pymatgen.analysis.chemenv.connectivity.connected_components import ConnectedComponent
from pymatgen.util.testing import PymatgenTest

json_files_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..",
                              'test_files', "chemenv", "json_test_files")


class FakeSite(object):
    def __hash__(self):
        return 1


class VoronoiContainerTest(PymatgenTest):

    @classmethod
    def setUpClass(cls):
        os.makedirs('tmp_dir')

    def test_elastic_centering(self):
        fake_site = FakeSite()
        graph = nx.MultiGraph()
        env_nodes = {}
        for isite in range(1, 7):
            env_nodes[isite] = get_environment_node(fake_site, isite, 'O:6')
            graph.add_node(env_nodes[isite])



        edges = [(1, 3, (-1, 0, 0)),
                 (3, 6, (1, 0, 0)),
                 (6, 5, (-1, 0, 0)),
                 (6, 2, (0, 0, 0)),
                 (2, 4, (0, 0, 0)),
                 (4, 5, (0, 0, 0)),
                 (5, 1, (1, 1, 0)),
                 ]

        # Node positions for visualisation of the graph if necessary
        # nodes_pos = {env_nodes[1]: [0.1, 0.1],
        #              env_nodes[2]: [0.3, 0.85],
        #              env_nodes[3]: [0.7, 0.35],
        #              env_nodes[4]: [0.5, 0.7],
        #              env_nodes[5]: [0.7, 0.85],
        #              env_nodes[6]: [0.1, 0.6],
        #              }



        for in1, in2, delta in edges:
            n1 = env_nodes[in1]
            n2 = env_nodes[in2]
            graph.add_edge(n1, n2, start=n1.isite, end=n2.isite, delta=delta)

        connected_component = ConnectedComponent.from_graph(graph)

        elastic_centered_graph = connected_component.elastic_centered_graph(start_node=env_nodes[1])

        for in1, in2 in [(1, 3), (1, 5), (3, 6), (5, 4), (2, 6), (6, 5)]:
            n1 = env_nodes[in1]
            n2 = env_nodes[in2]
            edges_n1_n2 = elastic_centered_graph[n1][n2]
            self.assertEqual(len(edges_n1_n2), 1)
            self.assertTupleEqual(edges_n1_n2[0]['delta'], (0, 0, 0))

    @classmethod
    def tearDownClass(cls):
        #Remove the directory in which the temporary files have been created
        shutil.rmtree('tmp_dir')

if __name__ == "__main__":
    fake_site = FakeSite()
    graph = nx.MultiGraph()
    env_nodes = {}
    for isite in range(1, 7):
        env_nodes[isite] = get_environment_node(fake_site, isite, 'O:6')
        graph.add_node(env_nodes[isite])

    edges = [(1, 3, (-1, 0, 0)),
             (3, 6, (1, 0, 0)),
             (6, 5, (-1, 0, 0)),
             (6, 2, (0, 0, 0)),
             (2, 4, (0, 0, 0)),
             (4, 5, (0, 0, 0)),
             (5, 1, (1, 1, 0)),
             ]

    # Node positions for visualisation of the graph if necessary
    # nodes_pos = {env_nodes[1]: [0.1, 0.1],
    #              env_nodes[2]: [0.3, 0.85],
    #              env_nodes[3]: [0.7, 0.35],
    #              env_nodes[4]: [0.5, 0.7],
    #              env_nodes[5]: [0.7, 0.85],
    #              env_nodes[6]: [0.1, 0.6],
    #              }



    for in1, in2, delta in edges:
        n1 = env_nodes[in1]
        n2 = env_nodes[in2]
        graph.add_edge(n1, n2, start=n1.isite, end=n2.isite, delta=delta)

    from networkx.algorithms.centrality import closeness_centrality, load_centrality, betweenness_centrality, current_flow_closeness_centrality, current_flow_betweenness_centrality, eigenvector_centrality_numpy, katz_centrality_numpy, harmonic_centrality
    from networkx.algorithms.distance_measures import eccentricity
    from networkx.algorithms.centrality import edge_betweenness_centrality, edge_current_flow_betweenness_centrality

    loadc = load_centrality(graph)
    closenessc = closeness_centrality(graph)
    betweennessc = betweenness_centrality(graph)
    cflowclosenessc = current_flow_closeness_centrality(graph)
    cflowbetweennessc = current_flow_betweenness_centrality(graph)
    eigvectorc = eigenvector_centrality_numpy(graph)
    harmc = harmonic_centrality(graph)
    ecc = eccentricity(graph)
    edge_betweennessc = edge_betweenness_centrality(graph)
    edge_currentflowc = edge_current_flow_betweenness_centrality(graph)
    print(ecc)
    for inode in range(1, 7):
        print(inode, loadc[env_nodes[inode]], closenessc[env_nodes[inode]], betweennessc[env_nodes[inode]],
              cflowclosenessc[env_nodes[inode]], cflowbetweennessc[env_nodes[inode]], eigvectorc[env_nodes[inode]],
              harmc[env_nodes[inode]])
    for inode in range(1, 7):
        print(inode, ecc[env_nodes[inode]])
    print(edge_betweennessc)
    for edge in graph.edges():
        print(edge[0].isite, edge[1].isite, edge_betweennessc[edge])
    from networkx.algorithms.traversal import bfs_tree
    tree = bfs_tree(graph, source=env_nodes[5])
    print(tree)
    from networkx.drawing import draw_networkx

    draw_networkx(graph, labels={env_nodes[ii]: str(ii) for ii in range(1, 7)})
    import matplotlib.pyplot as plt

    from networkx.drawing import draw_networkx
    plt.figure(2)
    draw_networkx(tree, labels={env_nodes[ii]: str(ii) for ii in range(1, 7)})
    import matplotlib.pyplot as plt
    plt.show()
    # unittest.main()