"""
Copyright (C) 2004-2022, NetworkX Developers
Aric Hagberg <hagberg@lanl.gov>
Dan Schult <dschult@colgate.edu>
Pieter Swart <swart@lanl.gov>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.

  * Neither the name of the NetworkX Developers nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from __future__ import annotations

import networkx as nx
import pytest


def test_graph_hash():
    pytest.importorskip("emmet.core.graph_hashing")

    from emmet.core.graph_hashing import weisfeiler_lehman_graph_hash

    G1 = nx.Graph()
    edges = [(1, 2, {"label": "A"}), (2, 3, {"label": "A"}), (3, 1, {"label": "A"}), (1, 4, {"label": "B"})]
    G1.add_edges_from(edges)
    G2 = nx.Graph()
    edges = [(5, 6, {"label": "B"}), (6, 7, {"label": "A"}), (7, 5, {"label": "A"}), (7, 8, {"label": "A"})]
    G2.add_edges_from(edges)

    assert weisfeiler_lehman_graph_hash(G1) == weisfeiler_lehman_graph_hash(G2)
    assert weisfeiler_lehman_graph_hash(G1, edge_attr="label") == "c653d85538bcf041d88c011f4f905f10"
    assert weisfeiler_lehman_graph_hash(G2, edge_attr="label") == "3dcd84af1ca855d0eff3c978d88e7ec7"


def test_subgraph_hashes():
    pytest.importorskip("emmet.core.graph_hashing")

    from emmet.core.graph_hashing import weisfeiler_lehman_subgraph_hashes

    G1 = nx.Graph()
    G1.add_edges_from([(1, 2), (2, 3), (2, 4), (3, 5), (4, 6), (5, 7), (6, 7)])
    G2 = nx.Graph()
    G2.add_edges_from([(1, 3), (2, 3), (1, 6), (1, 5), (4, 6)])

    g1_hashes = weisfeiler_lehman_subgraph_hashes(G1, iterations=3, digest_size=8)
    g2_hashes = weisfeiler_lehman_subgraph_hashes(G2, iterations=3, digest_size=8)

    assert g1_hashes[1] == ["a93b64973cfc8897", "db1b43ae35a1878f", "57872a7d2059c1c0"]
    assert g2_hashes[5] == ["a93b64973cfc8897", "db1b43ae35a1878f", "1716d2a4012fa4bc"]
