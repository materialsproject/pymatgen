import networkx as nx
import numpy as np

__author__ = 'waroquiers'


def get_all_simple_paths_edges(graph, source, target, cutoff=None, data=True):
    edge_paths = []
    if not graph.is_multigraph():
        for path in nx.all_simple_paths(graph, source, target, cutoff=cutoff):
            edge_paths.append([(path[ii], path[ii+1]) for ii in range(len(path) - 1)])
        return edge_paths

    node_paths = []
    for path in nx.all_simple_paths(graph, source, target, cutoff=cutoff):
        exists = False
        for path2 in node_paths:
            if len(path2) == len(path) and np.all(path == path2):
                exists = True
                break
        if exists:
            continue
        node_paths.append(path)
        current_edge_paths = [[]]
        for (node1, node2) in [(node1, path[inode1 + 1]) for inode1, node1 in enumerate(path[:-1])]:
            new_edge_paths = []
            for key, edge_data in graph[node1][node2].items():
                for tmp_edge_path in current_edge_paths:
                    if data:
                        new_path = [(n1, n2, k, d) for (n1, n2, k, d) in tmp_edge_path]
                        new_path.append((node1, node2, key, edge_data))
                        new_edge_paths.append(new_path)
                    else:
                        new_path = [(n1, n2, k) for (n1, n2, k) in tmp_edge_path]
                        new_path.append((node1, node2, key))
                        new_edge_paths.append(new_path)
            current_edge_paths = new_edge_paths
        edge_paths.extend(current_edge_paths)
    return edge_paths