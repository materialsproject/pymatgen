from __future__ import annotations

import warnings

import networkx as nx
import pytest
from monty.serialization import loadfn

from pymatgen.analysis.dimensionality import (
    calculate_dimensionality_of_site,
    get_dimensionality_cheon,
    get_dimensionality_gorai,
    get_dimensionality_larsen,
    get_structure_components,
    zero_d_graph_to_molecule_graph,
)
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core.structure import Structure
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestLarsenDimensionality(PymatgenTest):
    def setUp(self):
        cnn = CrystalNN()
        self.lifepo = cnn.get_bonded_structure(self.get_structure("LiFePO4"))
        self.graphite = cnn.get_bonded_structure(self.get_structure("Graphite"))
        self.cscl = cnn.get_bonded_structure(self.get_structure("CsCl"))

        tricky_structure = Structure(
            [5.79, 0.0, 0.0, 0, 5.79, 0.0, 0.0, 0.0, 5.79],
            ["B", "C", "C", "C", "C", "N", "N", "N", "N", "Ag"],
            [
                [0.0, 0.0, 0.0],
                [0.842, 0.842, 0.842],
                [0.158, 0.842, 0.158],
                [0.158, 0.158, 0.842],
                [0.842, 0.158, 0.158],
                [0.726, 0.726, 0.726],
                [0.274, 0.726, 0.274],
                [0.274, 0.274, 0.726],
                [0.726, 0.274, 0.274],
                [0.5, 0.5, 0.5],
            ],
        )
        self.tricky_structure = cnn.get_bonded_structure(tricky_structure)

        mol_structure = Structure(
            [[-2.316, 2.316, 2.160], [2.316, -2.316, 2.160], [2.316, 2.316, -2.160]],
            ["H", "C", "N"],
            [[0.752, 0.752, 0], [0.004, 0.004, 0.0], [0.272, 0.272, 0.0]],
        )
        self.mol_structure = cnn.get_bonded_structure(mol_structure)
        warnings.simplefilter("ignore")

    def tearDown(self) -> None:
        warnings.simplefilter("default")

    def test_get_dimensionality(self):
        assert get_dimensionality_larsen(self.lifepo) == 3
        assert get_dimensionality_larsen(self.graphite) == 2
        assert get_dimensionality_larsen(self.cscl) == 3

    def test_tricky_structure(self):
        """
        Test for a tricky structure that other dimensionality finders say is
        2D but is actually an interpenetrated 3D structure.
        """
        assert get_dimensionality_larsen(self.tricky_structure) == 3

    def test_get_structure_components(self):
        # test components are returned correctly with the right keys
        components = get_structure_components(self.tricky_structure)
        assert len(components) == 1
        assert components[0]["dimensionality"] == 3
        assert isinstance(components[0]["structure_graph"], StructureGraph)
        assert len(components[0]["structure_graph"]) == 10

        # test 2D structure and get orientation information
        components = get_structure_components(self.graphite, inc_orientation=True)
        assert len(components) == 2
        assert components[0]["dimensionality"] == 2
        assert isinstance(components[0]["structure_graph"], StructureGraph)
        assert len(components[0]["structure_graph"]) == 2
        assert components[0]["orientation"] == (0, 0, 1)

        # test getting component graphs
        assert list(components[0]["structure_graph"].graph.nodes()) == [0, 1]

    def test_calculate_dimensionality_of_site(self):
        dim = calculate_dimensionality_of_site(self.tricky_structure, 0)
        assert dim == 3

        # test vertices returned correctly
        dim, vertices = calculate_dimensionality_of_site(self.cscl, 0, inc_vertices=True)
        assert dim == 3
        assert len(vertices) == 4

    def test_zero_d_to_molecule_graph(self):
        comp_graphs = [
            self.mol_structure.graph.subgraph(c) for c in nx.weakly_connected_components(self.mol_structure.graph)
        ]

        mol_graph = zero_d_graph_to_molecule_graph(self.mol_structure, comp_graphs[0])

        assert mol_graph.get_connected_sites(0)[0].index == 1
        assert mol_graph.get_connected_sites(1)[1].index == 2
        assert len(mol_graph.molecule) == 3

        # test catching non zero dimensionality graphs
        comp_graphs = [self.graphite.graph.subgraph(c) for c in nx.weakly_connected_components(self.graphite.graph)]
        with pytest.raises(ValueError, match="Graph component is not zero-dimensional"):
            zero_d_graph_to_molecule_graph(self.graphite, comp_graphs[0])

        # test for a troublesome structure
        s = loadfn(f"{TEST_FILES_DIR}/PH7CN3O3F.json.gz")
        bs = CrystalNN().get_bonded_structure(s)
        comp_graphs = [bs.graph.subgraph(c) for c in nx.weakly_connected_components(bs.graph)]
        mol_graph = zero_d_graph_to_molecule_graph(bs, comp_graphs[0])
        assert len(mol_graph.molecule) == 12


class TestCheonDimensionality(PymatgenTest):
    def test_get_dimensionality(self):
        struct = self.get_structure("LiFePO4")
        assert get_dimensionality_cheon(struct) == "intercalated ion"

        struct = self.get_structure("Graphite")
        assert get_dimensionality_cheon(struct) == "2D"

    def test_get_dimensionality_with_bonds(self):
        struct = self.get_structure("CsCl")
        assert get_dimensionality_cheon(struct) == "intercalated ion"
        assert get_dimensionality_cheon(struct, ldict={"Cs": 3.7, "Cl": 3}) == "3D"

    def test_tricky_structure(self):
        tricky_structure = Structure(
            [5.79, 0.0, 0.0, 0, 5.79, 0.0, 0.0, 0.0, 5.79],
            ["B", "C", "C", "C", "C", "N", "N", "N", "N", "Ag"],
            [
                [0.0, 0.0, 0.0],
                [0.842, 0.842, 0.842],
                [0.158, 0.842, 0.158],
                [0.158, 0.158, 0.842],
                [0.842, 0.158, 0.158],
                [0.726, 0.726, 0.726],
                [0.274, 0.726, 0.274],
                [0.274, 0.274, 0.726],
                [0.726, 0.274, 0.274],
                [0.5, 0.5, 0.5],
            ],
        )

        # cheon dimensionality gets wrong structure using default parameters
        assert get_dimensionality_cheon(tricky_structure) == "2D"
        # cheon dimensionality gets tricky structure right using a
        # bigger supercell
        assert get_dimensionality_cheon(tricky_structure, larger_cell=True) == "3D"


class TestGoraiDimensionality(PymatgenTest):
    def test_get_dimensionality(self):
        struct = self.get_structure("LiFePO4")
        assert get_dimensionality_gorai(struct) == 3

        struct = self.get_structure("Graphite")
        assert get_dimensionality_gorai(struct) == 2

    def test_get_dimensionality_with_bonds(self):
        struct = self.get_structure("CsCl")
        assert get_dimensionality_gorai(struct) == 1
        assert get_dimensionality_gorai(struct, bonds={("Cs", "Cl"): 3.7}) == 3
