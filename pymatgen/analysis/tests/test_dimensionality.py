import os
import networkx as nx
import warnings
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.dimensionality import (
    get_dimensionality_gorai, get_dimensionality_cheon,
    get_dimensionality_larsen, calculate_dimensionality_of_site,
    get_structure_components, zero_d_graph_to_molecule_graph)
from pymatgen.util.testing import PymatgenTest
import unittest
from monty.serialization import loadfn

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class LarsenDimensionalityTest(PymatgenTest):

    def setUp(self):
        cnn = CrystalNN()
        self.lifepo = cnn.get_bonded_structure(self.get_structure('LiFePO4'))
        self.graphite = cnn.get_bonded_structure(self.get_structure('Graphite'))
        self.cscl = cnn.get_bonded_structure(self.get_structure('CsCl'))

        tricky_structure = Structure(
            [5.79, 0., 0., 0, 5.79, 0., 0., 0., 5.79],
            ['B', 'C', 'C', 'C', 'C', 'N', 'N', 'N', 'N', 'Ag'],
            [[0.0, 0.0, 0.0], [0.842, 0.842, 0.842], [0.158, 0.842, 0.158],
             [0.158, 0.158, 0.842], [0.842, 0.158, 0.158],
             [0.726, 0.726, 0.726], [0.274, 0.726, 0.274],
             [0.274, 0.274, 0.726], [0.726, 0.274, 0.274], [0.5, 0.5, 0.5]])
        self.tricky_structure = cnn.get_bonded_structure(tricky_structure)

        mol_structure = Structure(
            [[-2.316, 2.316, 2.160], [2.316, -2.316, 2.160],
             [2.316, 2.316, -2.160]],
            ['H', 'C', 'N'],
            [[0.752, 0.752, 0.000], [0.004, 0.004, 0.], [0.272, 0.272, 0.]])
        self.mol_structure = cnn.get_bonded_structure(mol_structure)
        warnings.simplefilter("ignore")

    def tearDown(self) -> None:
        warnings.simplefilter("default")

    def test_get_dimensionality(self):
        self.assertEqual(get_dimensionality_larsen(self.lifepo), 3)
        self.assertEqual(get_dimensionality_larsen(self.graphite), 2)
        self.assertEqual(get_dimensionality_larsen(self.cscl), 3)

    def test_tricky_structure(self):
        """
        Test for a tricky structure that other dimensionality finders say is
        2D but is actually an interpenetrated 3D structure.
        """
        self.assertEqual(get_dimensionality_larsen(self.tricky_structure), 3)

    def test_get_structure_components(self):
        # test components are returned correctly with the right keys
        components = get_structure_components(self.tricky_structure)
        self.assertEqual(len(components), 1)
        self.assertEqual(components[0]['dimensionality'], 3)
        self.assertTrue(
            isinstance(components[0]['structure_graph'], StructureGraph))
        self.assertEqual(components[0]['structure_graph'].structure.num_sites,
                         10)

        # test 2D structure and get orientation information
        components = get_structure_components(
            self.graphite, inc_orientation=True)
        self.assertEqual(len(components), 2)
        self.assertEqual(components[0]['dimensionality'], 2)
        self.assertTrue(
            isinstance(components[0]['structure_graph'], StructureGraph))
        self.assertEqual(components[0]['structure_graph'].structure.num_sites,
                         2)
        self.assertEqual(components[0]['orientation'], (0, 0, 1))

        # test getting component graphs
        self.assertEqual(list(components[0]['structure_graph'].graph.nodes()),
                         [0, 1])

    def test_calculate_dimensionality_of_site(self):
        dimen = calculate_dimensionality_of_site(self.tricky_structure, 0)
        self.assertEqual(dimen, 3)

        # test vertices returned correctly
        dimen, vertices = calculate_dimensionality_of_site(
            self.cscl, 0, inc_vertices=True)
        self.assertEqual(dimen, 3)
        self.assertEqual(len(vertices), 4)
        self.assertEqual(vertices[0], (-1, 1, 0))

    def test_zero_d_to_molecule_graph(self):
        comp_graphs = [self.mol_structure.graph.subgraph(c) for c in
                       nx.weakly_connected_components(self.mol_structure.graph)]

        mol_graph = zero_d_graph_to_molecule_graph(self.mol_structure,
                                                   comp_graphs[0])

        self.assertEqual(mol_graph.get_connected_sites(0)[0].index, 1)
        self.assertEqual(mol_graph.get_connected_sites(1)[1].index, 2)
        self.assertEqual(mol_graph.molecule.num_sites, 3)

        # test catching non zero dimensionality graphs
        comp_graphs = [self.graphite.graph.subgraph(c) for c in
                       nx.weakly_connected_components(self.graphite.graph)]
        self.assertRaises(ValueError, zero_d_graph_to_molecule_graph,
                          self.graphite, comp_graphs[0])

        # test for a troublesome structure
        s = loadfn(os.path.join(test_dir, "PH7CN3O3F.json.gz"))
        bs = CrystalNN().get_bonded_structure(s)
        comp_graphs = [bs.graph.subgraph(c) for c in
                       nx.weakly_connected_components(bs.graph)]
        mol_graph = zero_d_graph_to_molecule_graph(bs, comp_graphs[0])
        self.assertEqual(mol_graph.molecule.num_sites, 12)


class CheonDimensionalityTest(PymatgenTest):

    def test_get_dimensionality(self):
        s = self.get_structure('LiFePO4')
        self.assertEqual(get_dimensionality_cheon(s), 'intercalated ion')

        s = self.get_structure('Graphite')
        self.assertEqual(get_dimensionality_cheon(s), '2D')

    def test_get_dimensionality_with_bonds(self):
        s = self.get_structure('CsCl')
        self.assertEqual(get_dimensionality_cheon(s), 'intercalated ion')
        self.assertEqual(
            get_dimensionality_cheon(s, ldict={"Cs": 3.7, "Cl": 3}), '3D')

    def test_tricky_structure(self):
        tricky_structure = Structure(
            [5.79, 0., 0., 0, 5.79, 0., 0., 0., 5.79],
            ['B', 'C', 'C', 'C', 'C', 'N', 'N', 'N', 'N', 'Ag'],
            [[0.0, 0.0, 0.0], [0.842, 0.842, 0.842], [0.158, 0.842, 0.158],
             [0.158, 0.158, 0.842], [0.842, 0.158, 0.158],
             [0.726, 0.726, 0.726], [0.274, 0.726, 0.274],
             [0.274, 0.274, 0.726], [0.726, 0.274, 0.274], [0.5, 0.5, 0.5]])

        # cheon dimensionality gets wrong structure using default parameters
        self.assertEqual(get_dimensionality_cheon(tricky_structure), '2D')
        # cheon dimensionality gets tricky structure right using a
        # bigger supercell
        self.assertEqual(get_dimensionality_cheon(tricky_structure, larger_cell=True), '3D')


class GoraiDimensionalityTest(PymatgenTest):

    def test_get_dimensionality(self):
        s = self.get_structure('LiFePO4')
        self.assertEqual(get_dimensionality_gorai(s), 3)

        s = self.get_structure('Graphite')
        self.assertEqual(get_dimensionality_gorai(s), 2)

    def test_get_dimensionality_with_bonds(self):
        s = self.get_structure('CsCl')
        self.assertEqual(get_dimensionality_gorai(s), 1)
        self.assertEqual(get_dimensionality_gorai(s, bonds={("Cs", "Cl"): 3.7}),
                         3)


if __name__ == "__main__":
    unittest.main()
