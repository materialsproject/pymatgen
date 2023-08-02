---
layout: default
title: pymatgen.util.graph_hashing.md
nav_exclude: true
---

# pymatgen.util.graph_hashing module

Copyright (C) 2004-2022, NetworkX Developers
Aric Hagberg <[hagberg@lanl.gov](mailto:hagberg@lanl.gov)>
Dan Schult <[dschult@colgate.edu](mailto:dschult@colgate.edu)>
Pieter Swart <[swart@lanl.gov](mailto:swart@lanl.gov)>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

>
> * Redistributions of source code must retain the above copyright
> notice, this list of conditions and the following disclaimer.


> * Redistributions in binary form must reproduce the above
> copyright notice, this list of conditions and the following
> disclaimer in the documentation and/or other materials provided
> with the distribution.


> * Neither the name of the NetworkX Developers nor the names of its
> contributors may be used to endorse or promote products derived
> from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
“AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


---

Functions for hashing graphs to strings.
Isomorphic graphs should be assigned identical hashes.
For now, only Weisfeiler-Lehman hashing is implemented.


### pymatgen.util.graph_hashing.weisfeiler_lehman_graph_hash(G, edge_attr=None, node_attr=None, iterations=3, digest_size=16)
Return Weisfeiler Lehman (WL) graph hash.

The function iteratively aggregates and hashes neighborhoods of each node.
After each node’s neighbors are hashed to obtain updated node labels,
a hashed histogram of resulting labels is returned as the final hash.

Hashes are identical for isomorphic graphs and strong guarantees that
non-isomorphic graphs will get different hashes. See

```
[1]_
```

 for details.

If no node or edge attributes are provided, the degree of each node
is used as its initial label.
Otherwise, node and/or edge labels are used to compute the hash.


* **Parameters**


    * **G** – graph
    The graph to be hashed.
    Can have node and/or edge attributes. Can also have no attributes.


    * **edge_attr** – string, default=None
    The key in edge attribute dictionary to be used for hashing.
    If None, edge labels are ignored.


    * **node_attr** – string, default=None
    The key in node attribute dictionary to be used for hashing.
    If None, and no edge_attr given, use the degrees of the nodes as labels.


    * **iterations** – int, default=3
    Number of neighbor aggregations to perform.
    Should be larger for larger graphs.


    * **digest_size** – int, default=16
    Size (in bits) of blake2b hash digest to use for hashing node labels.



* **Returns**

    string

        Hexadecimal string corresponding to hash of the input graph.




* **Return type**

    h


### Notes

To return the WL hashes of each subgraph of a graph, use
weisfeiler_lehman_subgraph_hashes

Similarity between hashes does not imply similarity between graphs.

### References

Kurt Mehlhorn, and Karsten M. Borgwardt. Weisfeiler Lehman
Graph Kernels. Journal of Machine Learning Research. 2011.
[http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf](http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf)


### pymatgen.util.graph_hashing.weisfeiler_lehman_subgraph_hashes(G, edge_attr=None, node_attr=None, iterations=3, digest_size=16)
Return a dictionary of subgraph hashes by node.

The dictionary is keyed by node to a list of hashes in increasingly
sized induced subgraphs containing the nodes within 2\*k edges
of the key node for increasing integer k until all nodes are included.

The function iteratively aggregates and hashes neighborhoods of each node.
This is achieved for each step by replacing for each node its label from
the previous iteration with its hashed 1-hop neighborhood aggregate.
The new node label is then appended to a list of node labels for each
node.

To aggregate neighborhoods at each step for a node $n$, all labels of
nodes adjacent to $n$ are concatenated. If the edge_attr parameter is set,
labels for each neighboring node are prefixed with the value of this attribute
along the connecting edge from this neighbor to node $n$. The resulting string
is then hashed to compress this information into a fixed digest size.

Thus, at the $i$th iteration nodes within $2i$ distance influence any given
hashed node label. We can therefore say that at depth $i$ for node $n$
we have a hash for a subgraph induced by the $2i$-hop neighborhood of $n$.

Can be used to to create general Weisfeiler-Lehman graph kernels, or
generate features for graphs or nodes, for example to generate ‘words’ in a
graph as seen in the ‘graph2vec’ algorithm.
See

```
[1]_
```

 &

```
[2]_
```

 respectively for details.

Hashes are identical for isomorphic subgraphs and there exist strong
guarantees that non-isomorphic graphs will get different hashes.
See

```
[1]_
```

 for details.

If no node or edge attributes are provided, the degree of each node
is used as its initial label.
Otherwise, node and/or edge labels are used to compute the hash.


* **Parameters**


    * **G** – graph
    The graph to be hashed.
    Can have node and/or edge attributes. Can also have no attributes.


    * **edge_attr** – string, default=None
    The key in edge attribute dictionary to be used for hashing.
    If None, edge labels are ignored.


    * **node_attr** – string, default=None
    The key in node attribute dictionary to be used for hashing.
    If None, and no edge_attr given, use the degrees of the nodes as labels.


    * **iterations** – int, default=3
    Number of neighbor aggregations to perform.
    Should be larger for larger graphs.


    * **digest_size** – int, default=16
    Size (in bits) of blake2b hash digest to use for hashing node labels.
    The default size is 16 bits



* **Returns**

    dict

        A dictionary with each key given by a node in G, and each value given
        by the subgraph hashes in order of depth from the key node.




* **Return type**

    node_subgraph_hashes


### Notes

To hash the full graph when subgraph hashes are not needed, use
weisfeiler_lehman_graph_hash for efficiency.

Similarity between hashes does not imply similarity between graphs.

### References

Kurt Mehlhorn, and Karsten M. Borgwardt. Weisfeiler Lehman
Graph Kernels. Journal of Machine Learning Research. 2011.
[http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf](http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf)
.. [2] Annamalai Narayanan, Mahinthan Chandramohan, Rajasekar Venkatesan,
Lihui Chen, Yang Liu and Shantanu Jaiswa. graph2vec: Learning
Distributed Representations of Graphs. arXiv. 2017
[https://arxiv.org/pdf/1707.05005.pdf](https://arxiv.org/pdf/1707.05005.pdf)