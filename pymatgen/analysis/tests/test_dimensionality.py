
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.dimensionality import (
    get_dimensionality_gorai, get_dimensionality_cheon,
    get_dimensionality_larsen, calculate_dimensionality_of_site,
    get_structure_component_info)
from pymatgen.util.testing import PymatgenTest


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

    def test_get_structure_component_info(self):
        # test components are returned correctly with the right keys
        components = get_structure_component_info(self.tricky_structure)
        self.assertEqual(len(components), 1)
        self.assertEqual(components[0]['dimensionality'], 3)
        self.assertTrue(isinstance(components[0]['structure'], Structure))
        self.assertEqual(components[0]['structure'].num_sites, 10)

        # test 2D structure and get orientation information
        components = get_structure_component_info(
            self.graphite, inc_orientation=True)
        self.assertEqual(len(components), 2)
        self.assertEqual(components[0]['dimensionality'], 2)
        self.assertTrue(isinstance(components[0]['structure'], Structure))
        self.assertEqual(components[0]['structure'].num_sites, 2)
        self.assertEqual(components[0]['orientation'], (0, 0, 1))

    def test_calculate_dimensionality_of_site(self):
        dimen = calculate_dimensionality_of_site(self.tricky_structure, 0)
        self.assertEqual(dimen, 3)

        # test vertices returned correctly
        dimen, vertices = calculate_dimensionality_of_site(
            self.cscl, 0, inc_vertices=True)
        self.assertEqual(dimen, 3)
        self.assertEqual(len(vertices), 4)
        self.assertEqual(vertices[0], (-1, 1, 0))


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

        # cheon dimensionality gets wrong structure
        self.assertEqual(get_dimensionality_cheon(tricky_structure), '2D')


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
