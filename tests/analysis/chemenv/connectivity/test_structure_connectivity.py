from __future__ import annotations

import copy
from types import SimpleNamespace

import networkx as nx
import numpy as np
import orjson
import pytest

from pymatgen.analysis.chemenv.connectivity.connectivity_finder import ConnectivityFinder
from pymatgen.analysis.chemenv.connectivity.structure_connectivity import StructureConnectivity, get_delta_image
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
    MultiWeightsChemenvStrategy,
    SimplestChemenvStrategy,
)
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import (
    LightStructureEnvironments,
    StructureEnvironments,
)
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

__author__ = "waroquiers"

SE_DIR = f"{TEST_FILES_DIR}/analysis/chemenv/structure_environments"


def _load_se(fname: str) -> StructureEnvironments:
    with open(f"{SE_DIR}/{fname}", "rb") as file:
        return StructureEnvironments.from_dict(orjson.loads(file.read()))


@pytest.fixture(scope="module")
def BaTiO3_se():
    return _load_se("se_mp-5020.json")


@pytest.fixture(scope="module")
def BaTiO3_lse(BaTiO3_se):
    return LightStructureEnvironments.from_structure_environments(
        strategy=SimplestChemenvStrategy(), structure_environments=BaTiO3_se
    )


@pytest.fixture(scope="module")
def K2NaNb2Fe7Si8H4O31_se():
    # mp-743972: contains sites without any coordination environment (empty neighbors sets)
    # as well as sites whose environment is a mix of several environments.
    return _load_se("se_mp-743972.json")


class TestStructureConnectivity(MatSciTest):
    def test_serialization(self):
        BaTiO3_se_fpath = f"{TEST_FILES_DIR}/analysis/chemenv/structure_environments/se_mp-5020.json"
        with open(BaTiO3_se_fpath, "rb") as file:
            dd = orjson.loads(file.read())
        struct_envs = StructureEnvironments.from_dict(dd)
        lse = LightStructureEnvironments.from_structure_environments(
            strategy=SimplestChemenvStrategy(), structure_environments=struct_envs
        )
        cf = ConnectivityFinder()
        sc = cf.get_structure_connectivity(light_structure_environments=lse)
        sc_from_dict = StructureConnectivity.from_dict(sc.as_dict())
        assert sc.light_structure_environments == sc_from_dict.light_structure_environments
        assert set(sc._graph.nodes()) == set(sc_from_dict._graph.nodes())
        assert set(sc._graph.edges()) == set(sc_from_dict._graph.edges())

        sc_from_json = StructureConnectivity.from_dict(orjson.loads(orjson.dumps(sc.as_dict()).decode()))
        assert sc.light_structure_environments == sc_from_json.light_structure_environments
        assert set(sc._graph.nodes()) == set(sc_from_json._graph.nodes())
        assert set(sc._graph.edges()) == set(sc_from_json._graph.edges())

        json_str = self.assert_msonable(sc)
        sc_from_json = StructureConnectivity.from_dict(orjson.loads(json_str))
        assert sc.light_structure_environments == sc_from_json.light_structure_environments
        assert set(sc._graph.nodes()) == set(sc_from_json._graph.nodes())
        assert set(sc._graph.edges()) == set(sc_from_json._graph.edges())

    def test_environment_subgraph(self, BaTiO3_lse):
        sc = ConnectivityFinder().get_structure_connectivity(BaTiO3_lse)
        # 5 sites (Ba, Ti, 3 O), each of the 18 Ba-O and Ti-O bonds stored once
        assert sc._graph.number_of_nodes() == 5
        assert sc._graph.number_of_edges() == 18

        # Default subgraph uses all identified environments (Ba C:12 and Ti O:6)
        subgraph = sc.environment_subgraph()
        assert sorted(node.isite for node in subgraph.nodes()) == [0, 1]
        assert list(sc.environment_subgraphs) == ["C:12-O:6"]
        # Second call returns the already computed subgraph
        assert sc.environment_subgraph() is subgraph

        # Explicit environment symbols
        oct_subgraph = sc.environment_subgraph(environments_symbols=["O:6"])
        assert [node.isite for node in oct_subgraph.nodes()] == [1]
        # Corner-sharing octahedra: 3 self-loops along a, b and c (opposite images are merged)
        assert oct_subgraph.number_of_edges() == 3
        deltas = sorted(tuple(int(dd) for dd in np.abs(data["delta"])) for *_, data in oct_subgraph.edges(data=True))
        assert deltas == [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
        for *_, data in oct_subgraph.edges(data=True):
            assert len(data["ligands"]) == 1

        # Cached subgraph is reused when asking again for the same environments
        sc.setup_environment_subgraph(environments_symbols=["O:6"])
        assert sc._environment_subgraph is oct_subgraph
        assert list(sc.environment_subgraphs) == ["C:12-O:6", "O:6"]

        # only_atoms on an ordered structure
        sc.setup_environment_subgraph(environments_symbols=["O:6", "C:12"], only_atoms=["Ba"])
        assert [node.isite for node in sc._environment_subgraph.nodes()] == [0]
        assert "C:12-O:6#Ba" in sc.environment_subgraphs

    def test_environment_subgraph_disordered(self, BaTiO3_lse):
        lse = copy.deepcopy(BaTiO3_lse)
        lse.structure.replace(1, {"Ti": 0.5, "Zr": 0.5})
        assert not lse.structure.is_ordered
        sc = ConnectivityFinder().get_structure_connectivity(lse)
        sc.setup_environment_subgraph(environments_symbols=["O:6", "C:12"], only_atoms=["Zr"])
        assert [node.isite for node in sc._environment_subgraph.nodes()] == [1]
        assert sc._environment_subgraph.number_of_edges() == 3
        sc.setup_environment_subgraph(environments_symbols=["O:6"], only_atoms=["Ba"])
        assert sc._environment_subgraph.number_of_nodes() == 0

    def test_print_links(self, BaTiO3_lse, capsys):
        sc = ConnectivityFinder().get_structure_connectivity(BaTiO3_lse)
        sc.print_links()
        out = capsys.readouterr().out.splitlines()
        assert out[0] == "Links in graph :"
        assert "0  is connected with : " in out
        assert "1  is connected with : " in out
        # Ba (C:12) shares faces/edges/corners with Ti (O:6) and other Ba
        links = [line for line in out if line.startswith("  - ")]
        assert len(links) == 2 * sc.environment_subgraph().number_of_edges() - sum(
            1 for n1, n2 in sc.environment_subgraph().edges() if n1 == n2
        )
        # Links listed from both ends, with opposite delta images
        assert any(line.startswith("  - 0 by") for line in links)
        assert any(line.startswith("  - 1 by") for line in links)

    @pytest.mark.parametrize(
        ("isite1", "isite2", "data1", "data2", "expected"),
        [
            (0, 1, {"start": 0, "delta": (1, 0, 0)}, {"start": 1, "delta": (0, 1, 0)}, (1, -1, 0)),
            (0, 1, {"start": 0, "delta": (1, 0, 0)}, {"start": 5, "delta": (0, 1, 0)}, (1, 1, 0)),
            (0, 1, {"start": 5, "delta": (1, 0, 0)}, {"start": 1, "delta": (0, 1, 0)}, (-1, -1, 0)),
            (0, 1, {"start": 5, "delta": (1, 0, 0)}, {"start": 5, "delta": (0, 1, 0)}, (-1, 1, 0)),
        ],
    )
    def test_get_delta_image(self, isite1, isite2, data1, data2, expected):
        assert tuple(get_delta_image(isite1, isite2, data1, data2)) == expected

    def test_add_bonds(self):
        sc = StructureConnectivity(light_structure_environment=None)
        sc._graph.add_nodes_from(range(3))

        def nbs(*idx_images):
            return SimpleNamespace(
                neighb_indices_and_images=[{"index": idx, "image_cell": np.array(image)} for idx, image in idx_images]
            )

        # Same neighbor in the same cell twice: added only once
        sc.add_bonds(0, nbs((1, (0, 0, 0)), (1, (0, 0, 0))))
        assert sc._graph.number_of_edges(0, 1) == 1
        # Self-images: (1, 0, 0) and (-1, 0, 0) describe the same bond
        sc.add_bonds(2, nbs((2, (1, 0, 0)), (2, (-1, 0, 0)), (2, (0, 1, 0))))
        assert sc._graph.number_of_edges(2, 2) == 2
        # Periodic bond seen from both ends is only stored once
        sc.add_bonds(1, nbs((2, (1, 0, 0))))
        sc.add_bonds(2, nbs((1, (-1, 0, 0))))
        sc.add_bonds(1, nbs((2, (1, 0, 0))))
        assert sc._graph.number_of_edges(1, 2) == 1
        # A different image is a different bond
        sc.add_bonds(2, nbs((1, (0, 0, 1))))
        assert sc._graph.number_of_edges(1, 2) == 2
        edge_data = sorted(
            (data["start"], data["end"], tuple(data["delta"])) for *_, data in sc._graph.edges(data=True)
        )
        assert edge_data == [
            (0, 1, (0, 0, 0)),
            (1, 2, (1, 0, 0)),
            (2, 1, (0, 0, 1)),
            (2, 2, (0, 1, 0)),
            (2, 2, (1, 0, 0)),
        ]

        # Inconsistent edge data (edge not starting nor ending on the site)
        sc._graph.add_edge(0, 2, start=5, end=6, delta=(0, 0, 0))
        with pytest.raises(ValueError, match="SHOULD NOT HAPPEN"):
            sc.add_bonds(0, nbs((2, (1, 0, 0))))

    def test_not_implemented(self):
        sc = StructureConnectivity(light_structure_environment=None)
        assert sc.setup_connectivity_description() is None
        with pytest.raises(NotImplementedError):
            sc.setup_atom_environment_subgraph("Fe:O:6")
        with pytest.raises(NotImplementedError):
            sc.setup_environments_subgraph(["O:6"])
        with pytest.raises(NotImplementedError):
            sc.setup_atom_environments_subgraph(["Fe:O:6"])

    def test_init_with_graphs(self):
        graph = nx.MultiGraph()
        graph.add_edge(0, 1, start=0, end=1, delta=(0, 0, 0))
        subgraphs = {"O:6": nx.MultiGraph()}
        sc = StructureConnectivity(None, connectivity_graph=graph, environment_subgraphs=subgraphs)
        assert sc._graph is graph
        assert sc.environment_subgraphs is subgraphs


class TestConnectivityFinder:
    def test_invalid_option(self):
        with pytest.raises(ValueError, match="Option 'WRONG' for multiple_environments_choice is not allowed"):
            ConnectivityFinder(multiple_environments_choice="WRONG")

    def test_multiple_environments(self, BaTiO3_se, K2NaNb2Fe7Si8H4O31_se):
        strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
        lse = LightStructureEnvironments.from_structure_environments(
            strategy=strategy, structure_environments=BaTiO3_se
        )
        # Ba site is a mix of C:12, I:12 and AC:12 with this strategy
        assert len(lse.neighbors_sets[0]) == 3
        with pytest.raises(ValueError, match="Local environment of site 0 is a mix and nothing is asked about it"):
            ConnectivityFinder().get_structure_connectivity(lse)

        sc = ConnectivityFinder(multiple_environments_choice="TAKE_HIGHEST_FRACTION").get_structure_connectivity(lse)
        assert sc._graph.number_of_nodes() == 5
        assert sc._graph.number_of_edges() == 18  # Same as with the cuboctahedral C:12 environment only
        ccs = sc.get_connected_components()
        assert len(ccs) == 1
        assert ccs[0].periodicity == "3D"

        # Sites with an empty list of neighbors sets are skipped
        lse = LightStructureEnvironments.from_structure_environments(
            strategy=strategy, structure_environments=K2NaNb2Fe7Si8H4O31_se
        )
        assert lse.neighbors_sets[0] == []
        sc = ConnectivityFinder(multiple_environments_choice="TAKE_HIGHEST_FRACTION").get_structure_connectivity(lse)
        assert sc._graph.number_of_nodes() == 55
        assert sc._graph.number_of_edges() == 96
        assert sc._graph.degree(0) == 0
