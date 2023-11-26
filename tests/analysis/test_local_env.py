from __future__ import annotations

import unittest
from math import pi
from shutil import which
from typing import get_args

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.graphs import MoleculeGraph, StructureGraph
from pymatgen.analysis.local_env import (
    BrunnerNN_real,
    BrunnerNN_reciprocal,
    BrunnerNN_relative,
    CovalentBondNN,
    Critic2NN,
    CrystalNN,
    CutOffDictNN,
    EconNN,
    IsayevNN,
    JmolNN,
    LocalStructOrderParams,
    MinimumDistanceNN,
    MinimumOKeeffeNN,
    MinimumVIRENN,
    NearNeighbors,
    OpenBabelNN,
    ValenceIonicRadiusEvaluator,
    VoronoiNN,
    get_neighbors_of_site_with_index,
    metal_edge_extender,
    on_disorder_options,
    oxygen_edge_extender,
    site_is_of_motif_type,
    solid_angle,
)
from pymatgen.core import Element, Lattice, Molecule, Structure
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

test_dir = f"{TEST_FILES_DIR}/fragmenter_files"


class TestValenceIonicRadiusEvaluator(PymatgenTest):
    def setUp(self):
        """Setup MgO rocksalt structure for testing Vacancy."""
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [
            [0, 0, 0],
            [0.5, 0.5, 0],
            [0.5, 0, 0.5],
            [0, 0.5, 0.5],
            [0.5, 0, 0],
            [0, 0.5, 0],
            [0, 0, 0.5],
            [0.5, 0.5, 0.5],
        ]
        self._mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, validate_proximity=True, to_unit_cell=True)
        self._mgo_val_rad_evaluator = ValenceIonicRadiusEvaluator(self._mgo_uc)

    def test_valences_ionic_structure(self):
        valence_dict = self._mgo_val_rad_evaluator.valences
        for val in list(valence_dict.values()):
            assert val in {2, -2}

    def test_radii_ionic_structure(self):
        radii_dict = self._mgo_val_rad_evaluator.radii
        for rad in list(radii_dict.values()):
            assert rad in {0.86, 1.26}


class TestVoronoiNN(PymatgenTest):
    def setUp(self):
        self.struct = self.get_structure("LiFePO4")
        self.nn = VoronoiNN(targets=[Element("O")])
        self.s_sic = self.get_structure("Si")
        self.s_sic["Si"] = {"Si": 0.5, "C": 0.5}
        self.nn_sic = VoronoiNN()

    def test_get_voronoi_polyhedra(self):
        assert len(self.nn.get_voronoi_polyhedra(self.struct, 0).items()) == 8

    def test_get_cn(self):
        site_0_coord_num = self.nn.get_cn(self.struct, 0, use_weights=True, on_disorder="take_max_species")
        assert site_0_coord_num == approx(5.809265748999465, abs=1e-7)

        site_0_coord_num = self.nn_sic.get_cn(self.s_sic, 0, use_weights=True, on_disorder="take_max_species")
        assert site_0_coord_num == approx(4.5381161643940668, abs=1e-7)

    def test_get_coordinated_sites(self):
        assert len(self.nn.get_nn(self.struct, 0)) == 8

    def test_volume(self):
        self.nn.targets = None
        volume = 0
        for n in range(len(self.struct)):
            for nn in self.nn.get_voronoi_polyhedra(self.struct, n).values():
                volume += nn["volume"]
        assert self.struct.volume == approx(volume)

    def test_solid_angle(self):
        self.nn.targets = None
        for n in range(len(self.struct)):
            angle = 0
            for nn in self.nn.get_voronoi_polyhedra(self.struct, n).values():
                angle += nn["solid_angle"]
            assert 4 * np.pi == approx(angle)
        assert solid_angle([0, 0, 0], [[1, 0, 0], [-1, 0, 0], [0, 1, 0]]) == pi

    def test_nn_shell(self):
        # First, make a SC lattice. Make my math easier
        struct = Structure([[1, 0, 0], [0, 1, 0], [0, 0, 1]], ["Cu"], [[0, 0, 0]])

        # Get the 1NN shell
        self.nn.targets = None
        nns = self.nn.get_nn_shell_info(struct, 0, 1)
        assert len(nns) == 6

        # Test the 2nd NN shell
        nns = self.nn.get_nn_shell_info(struct, 0, 2)
        assert len(nns) == 18
        assert_allclose([1] * 6, [x["weight"] for x in nns if max(np.abs(x["image"])) == 2])
        assert_allclose([2] * 12, [x["weight"] for x in nns if max(np.abs(x["image"])) == 1])

        # Test the 3rd NN shell
        nns = self.nn.get_nn_shell_info(struct, 0, 3)
        for nn in nns:
            #  Check that the coordinates were set correctly
            assert_allclose(nn["site"].frac_coords, nn["image"])

        # Test with a structure that has unequal faces
        cscl = Structure(
            Lattice([[4.209, 0, 0], [0, 4.209, 0], [0, 0, 4.209]]),
            ["Cl1-", "Cs1+"],
            [[2.1045, 2.1045, 2.1045], [0, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.nn.weight = "area"
        nns = self.nn.get_nn_shell_info(cscl, 0, 1)
        assert len(nns) == 14
        assert np.isclose([x["weight"] for x in nns], 0.125 / 0.32476).sum() == 6  # Square faces
        assert np.isclose([x["weight"] for x in nns], 1).sum() == 8

        nns = self.nn.get_nn_shell_info(cscl, 0, 2)
        # Weight of getting back on to own site
        #  Square-square hop: 6*5 options times (0.125/0.32476)^2 weight each
        #  Hex-hex hop: 8*7 options times 1 weight each
        assert approx(np.sum([x["weight"] for x in nns if x["site_index"] == 0]), abs=1e-3) == 60.4444

    def test_adj_neighbors(self):
        # Make a simple cubic structure
        struct = Structure([[1, 0, 0], [0, 1, 0], [0, 0, 1]], ["Cu"], [[0, 0, 0]])

        # Compute the NNs with adjacency
        self.nn.targets = None
        neighbors = self.nn.get_voronoi_polyhedra(struct, 0)

        # Each neighbor has 4 adjacent neighbors, all orthogonal
        for nn_info in neighbors.values():
            assert len(nn_info["adj_neighbors"]) == 4

            for adj_key in nn_info["adj_neighbors"]:
                assert np.dot(nn_info["normal"], neighbors[adj_key]["normal"]) == approx(0)

    def test_all_at_once(self):
        # Get all of the sites for LiFePO4
        all_sites = self.nn.get_all_voronoi_polyhedra(self.struct)

        # Make sure they are the same as the single-atom ones
        for i, site in enumerate(all_sites):
            # Compute the tessellation using only one site
            by_one = self.nn.get_voronoi_polyhedra(self.struct, i)

            # Match the coordinates the of the neighbors, as site matching does not seem to work?
            all_coords = np.sort([x["site"].coords for x in site.values()], axis=0)
            by_one_coords = np.sort([x["site"].coords for x in by_one.values()], axis=0)

            assert_allclose(all_coords, by_one_coords)

        # Test the nn_info operation
        all_nn_info = self.nn.get_all_nn_info(self.struct)
        for i, info in enumerate(all_nn_info):
            # Compute using the by-one method
            by_one = self.nn.get_nn_info(self.struct, i)

            # Get the weights
            all_weights = sorted(x["weight"] for x in info)
            by_one_weights = sorted(x["weight"] for x in by_one)

            assert_allclose(all_weights, by_one_weights)

    def test_Cs2O(self):
        """A problematic structure in the Materials Project."""
        struct = Structure(
            [
                [4.358219, 0.192833, 6.406960],
                [2.114414, 3.815824, 6.406960],
                [0.311360, 0.192833, 7.742498],
            ],
            ["O", "Cs", "Cs"],
            [[0, 0, 0], [0.264318, 0.264318, 0.264318], [0.735682, 0.735682, 0.735682]],
            coords_are_cartesian=False,
        )

        # Compute the voronoi tessellation
        result = VoronoiNN().get_all_voronoi_polyhedra(struct)
        assert len(result) == 3

    def test_filtered(self):
        nn = VoronoiNN(weight="area")

        # Make a bcc crystal
        bcc = Structure(
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            ["Cu", "Cu"],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            coords_are_cartesian=False,
        )

        # Compute the weight of the little face
        big_face_area = np.sqrt(3) * 3 / 2 * (2 / 4 / 4)
        small_face_area = 0.125
        little_weight = small_face_area / big_face_area

        # Run one test where you get the small neighbors
        nn.tol = little_weight * 0.99
        nns = nn.get_nn_info(bcc, 0)
        assert len(nns) == 14

        # Run a second test where we screen out little faces
        nn.tol = little_weight * 1.01
        nns = nn.get_nn_info(bcc, 0)
        assert len(nns) == 8

        # Make sure it works for the `get_all` operation
        all_nns = nn.get_all_nn_info(bcc * [2, 2, 2])
        assert [
            8,
        ] * 16 == [len(x) for x in all_nns]

    def tearDown(self):
        del self.struct
        del self.nn


class TestJmolNN(PymatgenTest):
    def setUp(self):
        self.jmol = JmolNN()
        self.jmol_update = JmolNN(el_radius_updates={"Li": 1})

    def test_get_nn(self):
        struct = self.get_structure("LiFePO4")

        # Test the default near-neighbor finder.
        n_sites_checked = 0

        for site_idx, site in enumerate(struct):
            if site.specie == Element("Li"):
                assert self.jmol.get_cn(struct, site_idx) == 0
                n_sites_checked += 1
            elif site.specie == Element("Fe"):
                assert self.jmol.get_cn(struct, site_idx) == 6
                n_sites_checked += 1
            elif site.specie == Element("P"):
                assert self.jmol.get_cn(struct, site_idx) == 4
                n_sites_checked += 1
        assert n_sites_checked == 12

        # Test a user override that would cause Li to show up as 2-coordinated
        assert self.jmol_update.get_cn(struct, 0) == 2

        # Verify get_nn function works
        assert len(self.jmol_update.get_nn(struct, 0)) == 2


class TestIsayevNN(PymatgenTest):
    def test_get_nn(self):
        inn = IsayevNN()
        struct = self.get_structure("LiFePO4")

        assert inn.get_cn(struct, 0) == 2
        assert inn.get_cn(struct, 5) == 6
        assert inn.get_cn(struct, 10) == 4
        assert len(inn.get_nn(struct, 0)) == 2


class TestOpenBabelNN(PymatgenTest):
    def setUp(self):
        pytest.importorskip("openbabel")
        self.benzene = Molecule.from_file(f"{TEST_FILES_DIR}/benzene.xyz")
        self.acetylene = Molecule.from_file(f"{TEST_FILES_DIR}/acetylene.xyz")

    def test_nn_orders(self):
        strategy = OpenBabelNN()
        acetylene = strategy.get_nn_info(self.acetylene, 0)
        assert acetylene[0]["weight"] == 3
        assert acetylene[1]["weight"] == 1

        # Currently, benzene bonds register either as double or single,
        # not aromatic
        # Instead of searching for aromatic bonds, we check that bonds are
        # detected in the same way from both sides
        assert strategy.get_nn_info(self.benzene, 0)[0]["weight"] == strategy.get_nn_info(self.benzene, 1)[0]["weight"]

    def test_nn_length(self):
        strategy = OpenBabelNN(order=False)

        benzene_bonds = strategy.get_nn_info(self.benzene, 0)

        c_bonds = [b for b in benzene_bonds if str(b["site"].specie) == "C"]
        h_bonds = [b for b in benzene_bonds if str(b["site"].specie) == "H"]

        assert c_bonds[0]["weight"] == approx(1.41, abs=1e-2)
        assert h_bonds[0]["weight"] == approx(1.02, abs=1e-2)

        assert strategy.get_nn_info(self.acetylene, 0)[0]["weight"] == approx(1.19, abs=1e-2)


class TestCovalentBondNN(PymatgenTest):
    def setUp(self):
        self.benzene = Molecule.from_file(f"{TEST_FILES_DIR}/benzene.xyz")
        self.acetylene = Molecule.from_file(f"{TEST_FILES_DIR}/acetylene.xyz")

    def test_nn_orders(self):
        strategy = CovalentBondNN()

        acetylene = strategy.get_nn_info(self.acetylene, 0)
        assert acetylene[0]["weight"] == 3
        assert acetylene[1]["weight"] == 1

        benzene = strategy.get_nn_info(self.benzene, 0)
        assert benzene[0]["weight"] == approx(1.6596, abs=1e-4)

    def test_nn_length(self):
        strategy = CovalentBondNN(order=False)

        benzene_bonds = strategy.get_nn_info(self.benzene, 0)

        c_bonds = [b for b in benzene_bonds if str(b["site"].specie) == "C"]
        h_bonds = [b for b in benzene_bonds if str(b["site"].specie) == "H"]

        assert c_bonds[0]["weight"] == approx(1.41, abs=1e-2)
        assert h_bonds[0]["weight"] == approx(1.02, abs=1e-2)

        acetylene = strategy.get_nn_info(self.acetylene, 0)
        assert acetylene[0]["weight"] == approx(1.19, abs=1e-2)

    def test_bonded_structure(self):
        strategy = CovalentBondNN()

        benzene = strategy.get_bonded_structure(self.benzene)
        assert len(benzene.find_rings()) == 1

        acetylene = strategy.get_bonded_structure(self.acetylene)
        assert len(acetylene.graph.nodes) == 4


class TestMiniDistNN(PymatgenTest):
    def setUp(self):
        self.diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264], [0, 0, 2.528]]),
            ["C0+", "C0+"],
            [[2.554, 1.806, 4.423], [0.365, 0.258, 0.632]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.nacl = Structure(
            Lattice([[3.485, 0, 2.012], [1.162, 3.286, 2.012], [0, 0, 4.025]]),
            ["Na1+", "Cl1-"],
            [[0, 0, 0], [2.324, 1.643, 4.025]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.cscl = Structure(
            Lattice([[4.209, 0, 0], [0, 4.209, 0], [0, 0, 4.209]]),
            ["Cl1-", "Cs1+"],
            [[2.105, 2.105, 2.105], [0, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.mos2 = Structure(
            Lattice([[3.19, 0, 0], [-1.595, 2.763, 0], [0, 0, 17.44]]),
            ["Mo", "S", "S"],
            [[-1e-06, 1.842, 3.72], [1.595, 0.92, 5.29], [1.595, 0.92, 2.155]],
            coords_are_cartesian=True,
        )
        self.lifepo4 = self.get_structure("LiFePO4")
        self.lifepo4.add_oxidation_state_by_guess()

    def test_all_nn_classes(self):
        assert MinimumDistanceNN(cutoff=5, get_all_sites=True).get_cn(self.cscl, 0) == 14
        assert MinimumDistanceNN().get_cn(self.diamond, 0) == 4
        assert MinimumDistanceNN().get_cn(self.nacl, 0) == 6
        assert MinimumDistanceNN().get_cn(self.lifepo4, 0) == 6
        assert MinimumDistanceNN(tol=0.01).get_cn(self.cscl, 0) == 8
        assert MinimumDistanceNN(tol=0.1).get_cn(self.mos2, 0) == 6

        for image in MinimumDistanceNN(tol=0.1).get_nn_images(self.mos2, 0):
            assert image in [(0, 0, 0), (0, 1, 0), (-1, 0, 0), (0, 0, 0), (0, 1, 0), (-1, 0, 0)]

        okeeffe = MinimumOKeeffeNN(tol=0.01)
        assert okeeffe.get_cn(self.diamond, 0) == 4
        assert okeeffe.get_cn(self.nacl, 0) == 6
        assert okeeffe.get_cn(self.cscl, 0) == 8
        assert okeeffe.get_cn(self.lifepo4, 0) == 2

        virenn = MinimumVIRENN(tol=0.01)
        assert virenn.get_cn(self.diamond, 0) == 4
        assert virenn.get_cn(self.nacl, 0) == 6
        assert virenn.get_cn(self.cscl, 0) == 8
        assert virenn.get_cn(self.lifepo4, 0) == 2

        brunner_recip = BrunnerNN_reciprocal(tol=0.01)
        assert brunner_recip.get_cn(self.diamond, 0) == 4
        assert brunner_recip.get_cn(self.nacl, 0) == 6
        assert brunner_recip.get_cn(self.cscl, 0) == 14
        assert brunner_recip.get_cn(self.lifepo4, 0) == 6

        brunner_rel = BrunnerNN_relative(tol=0.01)
        assert brunner_rel.get_cn(self.diamond, 0) == 4
        assert brunner_rel.get_cn(self.nacl, 0) == 6
        assert brunner_rel.get_cn(self.cscl, 0) == 14
        assert brunner_rel.get_cn(self.lifepo4, 0) == 6

        brunner_real = BrunnerNN_real(tol=0.01)
        assert brunner_real.get_cn(self.diamond, 0) == 4
        assert brunner_real.get_cn(self.nacl, 0) == 6
        assert brunner_real.get_cn(self.cscl, 0) == 14
        assert brunner_real.get_cn(self.lifepo4, 0) == 30

        econn = EconNN()
        assert econn.get_cn(self.diamond, 0) == 4
        assert econn.get_cn(self.nacl, 0) == 6
        assert econn.get_cn(self.cscl, 0) == 14
        assert econn.get_cn(self.lifepo4, 0) == 6

        voroinn = VoronoiNN(tol=0.5)
        assert voroinn.get_cn(self.diamond, 0) == 4
        assert voroinn.get_cn(self.nacl, 0) == 6
        assert voroinn.get_cn(self.cscl, 0) == 8
        assert voroinn.get_cn(self.lifepo4, 0) == 6

        crystalnn = CrystalNN()
        assert crystalnn.get_cn(self.diamond, 0) == 4
        assert crystalnn.get_cn(self.nacl, 0) == 6
        assert crystalnn.get_cn(self.cscl, 0) == 8
        assert crystalnn.get_cn(self.lifepo4, 0) == 6

    def test_get_local_order_params(self):
        nn = MinimumDistanceNN()
        ops = nn.get_local_order_parameters(self.diamond, 0)
        assert ops["tetrahedral"] == approx(0.9999934389036574)

        ops = nn.get_local_order_parameters(self.nacl, 0)
        assert ops["octahedral"] == approx(0.9999995266669)


class TestMotifIdentification(PymatgenTest):
    def setUp(self):
        self.silicon = Structure(
            Lattice.cubic(5.47),
            ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"],
            [
                [0.000000, 0.000000, 0.500000],
                [0.750000, 0.750000, 0.750000],
                [0.000000, 0.500000, 1.000000],
                [0.750000, 0.250000, 0.250000],
                [0.500000, 0.000000, 1.000000],
                [0.250000, 0.750000, 0.250000],
                [0.500000, 0.500000, 0.500000],
                [0.250000, 0.250000, 0.750000],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=False,
            site_properties=None,
        )
        self.diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264], [0, 0, 2.528]]),
            ["C0+", "C0+"],
            [[2.554, 1.806, 4.423], [0.365, 0.258, 0.632]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.nacl = Structure(
            Lattice([[3.485, 0, 2.012], [1.162, 3.286, 2.012], [0, 0, 4.025]]),
            ["Na1+", "Cl1-"],
            [[0, 0, 0], [2.324, 1.643, 4.025]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.cscl = Structure(
            Lattice([[4.209, 0, 0], [0, 4.209, 0], [0, 0, 4.209]]),
            ["Cl1-", "Cs1+"],
            [[2.105, 2.105, 2.105], [0, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.square_pyramid = Structure(
            Lattice([[100, 0, 0], [0, 100, 0], [0, 0, 100]]),
            ["C", "C", "C", "C", "C", "C"],
            [[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.trigonal_bipyramid = Structure(
            Lattice([[100, 0, 0], [0, 100, 0], [0, 0, 100]]),
            ["P", "Cl", "Cl", "Cl", "Cl", "Cl"],
            [
                [0, 0, 0],
                [0, 0, 2.14],
                [0, 2.02, 0],
                [1.74937, -1.01, 0],
                [-1.74937, -1.01, 0],
                [0, 0, -2.14],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )

    def test_site_is_of_motif_type(self):
        for i in range(len(self.diamond)):
            assert site_is_of_motif_type(self.diamond, i) == "tetrahedral"
        for i in range(len(self.nacl)):
            assert site_is_of_motif_type(self.nacl, i) == "octahedral"
        for i in range(len(self.cscl)):
            assert site_is_of_motif_type(self.cscl, i) == "bcc"
        assert site_is_of_motif_type(self.square_pyramid, 0) == "square pyramidal"
        for i in range(1, len(self.square_pyramid)):
            assert site_is_of_motif_type(self.square_pyramid, i) == "unrecognized"
        assert site_is_of_motif_type(self.trigonal_bipyramid, 0) == "trigonal bipyramidal"
        for i in range(1, len(self.trigonal_bipyramid)):
            assert site_is_of_motif_type(self.trigonal_bipyramid, i) == "unrecognized"

    def test_get_neighbors_of_site_with_index(self):
        assert len(get_neighbors_of_site_with_index(self.diamond, 0)) == 4
        assert len(get_neighbors_of_site_with_index(self.nacl, 0)) == 6
        assert len(get_neighbors_of_site_with_index(self.cscl, 0)) == 8
        assert len(get_neighbors_of_site_with_index(self.diamond, 0, delta=0.01)) == 4
        assert len(get_neighbors_of_site_with_index(self.diamond, 0, cutoff=6)) == 4
        assert len(get_neighbors_of_site_with_index(self.diamond, 0, approach="voronoi")) == 4
        assert len(get_neighbors_of_site_with_index(self.diamond, 0, approach="min_OKeeffe")) == 4
        assert len(get_neighbors_of_site_with_index(self.diamond, 0, approach="min_VIRE")) == 4

    def tearDown(self):
        del self.silicon
        del self.diamond
        del self.nacl
        del self.cscl


class TestNearNeighbor(PymatgenTest):
    def setUp(self):
        self.diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264], [0, 0, 2.528]]),
            ["C0+", "C0+"],
            [[2.554, 1.806, 4.423], [0.365, 0.258, 0.632]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )

    def set_nn_info(self):
        # check conformance
        # implicitly assumes that all NearNeighbors subclasses
        # will correctly identify bonds in diamond, if it
        # can't there are probably bigger problems
        subclasses = NearNeighbors.__subclasses__()
        for subclass in subclasses:
            # Critic2NN has external dependency, is tested separately
            if "Critic2" not in str(subclass):
                nn_info = subclass().get_nn_info(self.diamond, 0)
                assert nn_info[0]["site_index"] == 1
                assert nn_info[0]["image"][0] == 1

    def test_on_disorder_options(self):
        assert get_args(on_disorder_options) == (
            "take_majority_strict",
            "take_majority_drop",
            "take_max_species",
            "error",
        )


class TestLocalStructOrderParams(PymatgenTest):
    def setUp(self):
        self.single_bond = Structure(
            Lattice.cubic(10),
            ["H", "H", "H"],
            [[1, 0, 0], [0, 0, 0], [6, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.linear = Structure(
            Lattice.cubic(10),
            ["H", "H", "H"],
            [[1, 0, 0], [0, 0, 0], [2, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.bent45 = Structure(
            Lattice.cubic(10),
            ["H", "H", "H"],
            [[0, 0, 0], [0.707, 0.707, 0], [0.707, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.cubic = Structure(
            Lattice.cubic(1),
            ["H"],
            [[0, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=False,
            site_properties=None,
        )
        self.bcc = Structure(
            Lattice.cubic(1),
            ["H", "H"],
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=False,
            site_properties=None,
        )
        self.fcc = Structure(
            Lattice.cubic(1),
            ["H", "H", "H", "H"],
            [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=False,
            site_properties=None,
        )
        self.hcp = Structure(
            Lattice.hexagonal(1, 1.633),
            ["H", "H"],
            [[0.3333, 0.6667, 0.25], [0.6667, 0.3333, 0.75]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=False,
            site_properties=None,
        )
        self.diamond = Structure(
            Lattice.cubic(1),
            ["H", "H", "H", "H", "H", "H", "H", "H"],
            [
                [0, 0, 0.5],
                [0.75, 0.75, 0.75],
                [0, 0.5, 0],
                [0.75, 0.25, 0.25],
                [0.5, 0, 0],
                [0.25, 0.75, 0.25],
                [0.5, 0.5, 0.5],
                [0.25, 0.25, 0.75],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=False,
            site_properties=None,
        )
        self.trigonal_off_plane = Structure(
            Lattice.cubic(100),
            ["H", "H", "H", "H"],
            [
                [0.50, 0.50, 0.50],
                [0.25, 0.75, 0.25],
                [0.25, 0.25, 0.75],
                [0.75, 0.25, 0.25],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.regular_triangle = Structure(
            Lattice.cubic(30),
            ["H", "H", "H", "H"],
            [[15, 15.28867, 15.65], [14.5, 15, 15], [15.5, 15, 15], [15, 15.866, 15]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.trigonal_planar = Structure(
            Lattice.cubic(30),
            ["H", "H", "H", "H"],
            [[15, 15.28867, 15], [14.5, 15, 15], [15.5, 15, 15], [15, 15.866, 15]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.square_planar = Structure(
            Lattice.cubic(30),
            ["H", "H", "H", "H", "H"],
            [
                [15, 15, 15],
                [14.75, 14.75, 15],
                [14.75, 15.25, 15],
                [15.25, 14.75, 15],
                [15.25, 15.25, 15],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.square = Structure(
            Lattice.cubic(30),
            ["H", "H", "H", "H", "H"],
            [
                [15, 15, 15.707],
                [14.75, 14.75, 15],
                [14.75, 15.25, 15],
                [15.25, 14.75, 15],
                [15.25, 15.25, 15],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.T_shape = Structure(
            Lattice.cubic(30),
            ["H", "H", "H", "H"],
            [[15, 15, 15], [15, 15, 15.5], [15, 15.5, 15], [15, 14.5, 15]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.square_pyramid = Structure(
            Lattice.cubic(30),
            ["H", "H", "H", "H", "H", "H"],
            [
                [15, 15, 15],
                [15, 15, 15.3535],
                [14.75, 14.75, 15],
                [14.75, 15.25, 15],
                [15.25, 14.75, 15],
                [15.25, 15.25, 15],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.pentagonal_planar = Structure(
            Lattice.cubic(30),
            ["Xe", "F", "F", "F", "F", "F"],
            [
                [0, -1.6237, 0],
                [1.17969, 0, 0],
                [-1.17969, 0, 0],
                [1.90877, -2.24389, 0],
                [-1.90877, -2.24389, 0],
                [0, -3.6307, 0],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.pentagonal_pyramid = Structure(
            Lattice.cubic(30),
            ["Xe", "F", "F", "F", "F", "F", "F"],
            [
                [0, -1.6237, 0],
                [0, -1.6237, 1.17969],
                [1.17969, 0, 0],
                [-1.17969, 0, 0],
                [1.90877, -2.24389, 0],
                [-1.90877, -2.24389, 0],
                [0, -3.6307, 0],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.pentagonal_bipyramid = Structure(
            Lattice.cubic(30),
            ["Xe", "F", "F", "F", "F", "F", "F", "F"],
            [
                [0, -1.6237, 0],
                [0, -1.6237, -1.17969],
                [0, -1.6237, 1.17969],
                [1.17969, 0, 0],
                [-1.17969, 0, 0],
                [1.90877, -2.24389, 0],
                [-1.90877, -2.24389, 0],
                [0, -3.6307, 0],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.hexagonal_planar = Structure(
            Lattice.cubic(30),
            ["H", "C", "C", "C", "C", "C", "C"],
            [
                [0, 0, 0],
                [0.71, 1.2298, 0],
                [-0.71, 1.2298, 0],
                [0.71, -1.2298, 0],
                [-0.71, -1.2298, 0],
                [1.4199, 0, 0],
                [-1.4199, 0, 0],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.hexagonal_pyramid = Structure(
            Lattice.cubic(30),
            ["H", "Li", "C", "C", "C", "C", "C", "C"],
            [
                [0, 0, 0],
                [0, 0, 1.675],
                [0.71, 1.2298, 0],
                [-0.71, 1.2298, 0],
                [0.71, -1.2298, 0],
                [-0.71, -1.2298, 0],
                [1.4199, 0, 0],
                [-1.4199, 0, 0],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.hexagonal_bipyramid = Structure(
            Lattice.cubic(30),
            ["H", "Li", "Li", "C", "C", "C", "C", "C", "C"],
            [
                [0, 0, 0],
                [0, 0, 1.675],
                [0, 0, -1.675],
                [0.71, 1.2298, 0],
                [-0.71, 1.2298, 0],
                [0.71, -1.2298, 0],
                [-0.71, -1.2298, 0],
                [1.4199, 0, 0],
                [-1.4199, 0, 0],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.trigonal_pyramid = Structure(
            Lattice.cubic(30),
            ["P", "Cl", "Cl", "Cl", "Cl"],
            [
                [0, 0, 0],
                [0, 0, 2.14],
                [0, 2.02, 0],
                [1.74937, -1.01, 0],
                [-1.74937, -1.01, 0],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.trigonal_bipyramidal = Structure(
            Lattice.cubic(30),
            ["P", "Cl", "Cl", "Cl", "Cl", "Cl"],
            [
                [0, 0, 0],
                [0, 0, 2.14],
                [0, 2.02, 0],
                [1.74937, -1.01, 0],
                [-1.74937, -1.01, 0],
                [0, 0, -2.14],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.cuboctahedron = Structure(
            Lattice.cubic(30),
            ["H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"],
            [
                [15, 15, 15],
                [15, 14.5, 14.5],
                [15, 14.5, 15.5],
                [15, 15.5, 14.5],
                [15, 15.5, 15.5],
                [14.5, 15, 14.5],
                [14.5, 15, 15.5],
                [15.5, 15, 14.5],
                [15.5, 15, 15.5],
                [14.5, 14.5, 15],
                [14.5, 15.5, 15],
                [15.5, 14.5, 15],
                [15.5, 15.5, 15],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.see_saw_rect = Structure(
            Lattice.cubic(30),
            ["H", "H", "H", "H", "H"],
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, -1.0],
                [-1.0, 0.0, 0.0],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )
        self.sq_face_capped_trig_pris = Structure(
            Lattice.cubic(30),
            ["H", "H", "H", "H", "H", "H", "H", "H"],
            [
                [0, 0, 0],
                [-0.6546536707079771, -0.37796447300922725, 0.6546536707079771],
                [0.6546536707079771, -0.37796447300922725, 0.6546536707079771],
                [0.0, 0.7559289460184545, 0.6546536707079771],
                [-0.6546536707079771, -0.37796447300922725, -0.6546536707079771],
                [0.6546536707079771, -0.37796447300922725, -0.6546536707079771],
                [0.0, 0.7559289460184545, -0.6546536707079771],
                [0.0, -1.0, 0.0],
            ],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None,
        )

    def test_init(self):
        assert LocalStructOrderParams(["cn"], parameters=None, cutoff=0.99) is not None

        parameters = [{"norm": 2}]
        lostops = LocalStructOrderParams(["cn"], parameters=parameters)
        tmp = lostops.get_parameters(0)
        parameters[0]["norm"] = 3
        assert tmp == lostops.get_parameters(0)

    def test_get_order_parameters(self):
        # Set up everything.
        op_types = [
            "cn",
            "bent",
            "bent",
            "tet",
            "oct",
            "bcc",
            "q2",
            "q4",
            "q6",
            "reg_tri",
            "sq",
            "sq_pyr_legacy",
            "tri_bipyr",
            "sgl_bd",
            "tri_plan",
            "sq_plan",
            "pent_plan",
            "sq_pyr",
            "tri_pyr",
            "pent_pyr",
            "hex_pyr",
            "pent_bipyr",
            "hex_bipyr",
            "T",
            "cuboct",
            "see_saw_rect",
            "hex_plan_max",
            "tet_max",
            "oct_max",
            "tri_plan_max",
            "sq_plan_max",
            "pent_plan_max",
            "cuboct_max",
            "tet_max",
            "sq_face_cap_trig_pris",
        ]
        op_params = [None for i in range(len(op_types))]
        op_params[1] = {"TA": 1, "IGW_TA": 1.0 / 0.0667}
        op_params[2] = {"TA": 45.0 / 180, "IGW_TA": 1.0 / 0.0667}
        op_params[33] = {
            "TA": 0.6081734479693927,
            "IGW_TA": 18.33,
            "fac_AA": 1.5,
            "exp_cos_AA": 2,
        }
        ops_044 = LocalStructOrderParams(op_types, parameters=op_params, cutoff=0.44)
        ops_071 = LocalStructOrderParams(op_types, parameters=op_params, cutoff=0.71)
        ops_087 = LocalStructOrderParams(op_types, parameters=op_params, cutoff=0.87)
        ops_099 = LocalStructOrderParams(op_types, parameters=op_params, cutoff=0.99)
        ops_101 = LocalStructOrderParams(op_types, parameters=op_params, cutoff=1.01)
        ops_501 = LocalStructOrderParams(op_types, parameters=op_params, cutoff=5.01)
        _ = LocalStructOrderParams(op_types, parameters=op_params)

        # Single bond.
        op_vals = ops_101.get_order_parameters(self.single_bond, 0)
        assert op_vals[13] == approx(1)
        op_vals = ops_501.get_order_parameters(self.single_bond, 0)
        assert op_vals[13] == approx(0.7999999)
        op_vals = ops_101.get_order_parameters(self.linear, 0)
        assert op_vals[13] == approx(0.0)

        # Linear motif.
        op_vals = ops_101.get_order_parameters(self.linear, 0)
        assert op_vals[1] == approx(1)

        # 45 degrees-bent motif.
        op_vals = ops_101.get_order_parameters(self.bent45, 0)
        assert op_vals[2] == approx(1)

        # T-shape motif.
        op_vals = ops_101.get_order_parameters(self.T_shape, 0, indices_neighs=[1, 2, 3])
        assert op_vals[23] == approx(1)

        # Cubic structure.
        op_vals = ops_099.get_order_parameters(self.cubic, 0)
        assert op_vals[0] == approx(0.0)
        assert op_vals[3:9] == [None] * 6
        op_vals = ops_101.get_order_parameters(self.cubic, 0)
        assert op_vals[0] == approx(6.0)
        assert op_vals[3:9] == approx([0.023, 1, 0.333, 0, 0.763, 0.353], abs=1e-3)
        assert op_vals[28] == 1

        # Bcc structure.
        op_vals = ops_087.get_order_parameters(self.bcc, 0)
        assert op_vals[0] == approx(8.0)
        assert op_vals[3:9] == approx([0.2, 0.145, 1, 0, 0.509, 0.628], abs=1e-3)

        # Fcc structure.
        op_vals = ops_071.get_order_parameters(self.fcc, 0)
        assert op_vals[0] == approx(12.0)
        assert op_vals[3:9] == approx([0.036, 0.078, -0.002, 0, 0.19, 0.574], abs=1e-3)

        # Hcp structure.
        op_vals = ops_101.get_order_parameters(self.hcp, 0)
        assert op_vals[0] == approx(12.0)
        assert op_vals[3:9] == approx([0.033, 0.082, -0.018, 0, 0.097, 0.484], abs=1e-3)

        # Diamond structure.
        op_vals = ops_044.get_order_parameters(self.diamond, 0)
        assert op_vals[0] == approx(4.0)
        assert op_vals[3:9] == approx([1, 0.037, 0.75, 0, 0.509, 0.628], abs=1e-3)
        assert op_vals[27] == 1

        # Trigonal off-plane molecule.
        op_vals = ops_044.get_order_parameters(self.trigonal_off_plane, 0)
        assert op_vals[0] == approx(3.0)
        assert op_vals[3] == 1
        assert op_vals[33] == 1

        # Trigonal-planar motif.
        op_vals = ops_101.get_order_parameters(self.trigonal_planar, 0)
        assert int(op_vals[0] + 0.5) == 3
        assert op_vals[14] == approx(1)
        assert op_vals[29] == approx(1)

        # Regular triangle motif.
        op_vals = ops_101.get_order_parameters(self.regular_triangle, 0)
        assert op_vals[9] == approx(0.9999999)

        # Square-planar motif.
        op_vals = ops_101.get_order_parameters(self.square_planar, 0)
        assert op_vals[15] == approx(1)
        assert op_vals[30] == approx(1)

        # Square motif.
        op_vals = ops_101.get_order_parameters(self.square, 0)
        assert op_vals[10] == approx(1)

        # Pentagonal planar.
        op_vals = ops_101.get_order_parameters(self.pentagonal_planar.sites, 0, indices_neighs=[1, 2, 3, 4, 5])
        assert op_vals[12] == approx(0.1260699690)
        assert op_vals[16] == approx(1)
        assert op_vals[31] == approx(1)

        # Trigonal pyramid motif.
        op_vals = ops_101.get_order_parameters(self.trigonal_pyramid, 0, indices_neighs=[1, 2, 3, 4])
        assert op_vals[18] == approx(1)

        # Square pyramid motif.
        op_vals = ops_101.get_order_parameters(self.square_pyramid, 0)
        assert op_vals[11] == approx(1)
        assert op_vals[12] == approx(2 / 3)
        assert op_vals[17] == approx(1)

        # Pentagonal pyramid motif.
        op_vals = ops_101.get_order_parameters(self.pentagonal_pyramid, 0, indices_neighs=[1, 2, 3, 4, 5, 6])
        assert op_vals[19] == approx(1)

        # Hexagonal pyramid motif.
        op_vals = ops_101.get_order_parameters(self.hexagonal_pyramid, 0, indices_neighs=[1, 2, 3, 4, 5, 6, 7])
        assert op_vals[20] == approx(1)

        # Trigonal bipyramidal.
        op_vals = ops_101.get_order_parameters(self.trigonal_bipyramidal.sites, 0, indices_neighs=[1, 2, 3, 4, 5])
        assert op_vals[12] == approx(1)

        # Pentagonal bipyramidal.
        op_vals = ops_101.get_order_parameters(self.pentagonal_bipyramid.sites, 0, indices_neighs=[1, 2, 3, 4, 5, 6, 7])
        assert op_vals[21] == approx(1)

        # Hexagonal bipyramid motif.
        op_vals = ops_101.get_order_parameters(self.hexagonal_bipyramid, 0, indices_neighs=[1, 2, 3, 4, 5, 6, 7, 8])
        assert op_vals[22] == approx(1)

        # Cuboctahedral motif.
        op_vals = ops_101.get_order_parameters(self.cuboctahedron, 0, indices_neighs=list(range(1, 13)))
        assert op_vals[24] == approx(1)
        assert op_vals[32] == approx(1)

        # See-saw motif.
        op_vals = ops_101.get_order_parameters(self.see_saw_rect, 0, indices_neighs=list(range(1, 5)))
        assert op_vals[25] == approx(1)

        # Hexagonal planar motif.
        op_vals = ops_101.get_order_parameters(self.hexagonal_planar, 0, indices_neighs=[1, 2, 3, 4, 5, 6])
        assert op_vals[26] == approx(1)

        # Square face capped trigonal prism.
        op_vals = ops_101.get_order_parameters(self.sq_face_capped_trig_pris, 0, indices_neighs=list(range(1, 8)))
        assert op_vals[34] == approx(1)

        # Test providing explicit neighbor lists.
        op_vals = ops_101.get_order_parameters(self.bcc, 0, indices_neighs=[1])
        assert op_vals[0] is not None
        assert op_vals[3] is None
        with pytest.raises(ValueError, match="Neighbor site index beyond maximum!"):
            ops_101.get_order_parameters(self.bcc, 0, indices_neighs=[2])


class TestCrystalNN(PymatgenTest):
    def setUp(self):
        self.lifepo4 = self.get_structure("LiFePO4")
        self.lifepo4.add_oxidation_state_by_guess()
        self.he_bcc = self.get_structure("He_BCC")
        self.he_bcc.add_oxidation_state_by_guess()

        self.disordered_struct = Structure(
            Lattice.cubic(3), [{"Fe": 0.4, "C": 0.3, "Mn": 0.3}, "O"], [[0, 0, 0], [0.5, 0.5, 0.5]]
        )
        self.disordered_struct_with_majority = Structure(
            Lattice.cubic(3), [{"Fe": 0.6, "C": 0.4}, "O"], [[0, 0, 0], [0.5, 0.5, 0.5]]
        )

    def test_sanity(self):
        cnn = CrystalNN()
        expected_msg = "The weighted_cn parameter and use_weights parameter should match"
        with pytest.raises(ValueError, match=expected_msg):
            cnn.get_cn(self.lifepo4, 0, use_weights=True)

        cnn = CrystalNN(weighted_cn=True)
        with pytest.raises(ValueError, match=expected_msg):
            cnn.get_cn(self.lifepo4, 0, use_weights=False)

    def test_discrete_cn(self):
        cnn = CrystalNN()
        cn_array = []
        expected_array = 8 * [6] + 20 * [4]
        for idx, _ in enumerate(self.lifepo4):
            cn_array.append(cnn.get_cn(self.lifepo4, idx))

        assert cn_array == expected_array

    def test_weighted_cn(self):
        cnn = CrystalNN(weighted_cn=True)
        cn_array = []

        # fmt: off
        expected_array = [
            5.863, 5.8716, 5.863, 5.8716, 5.7182, 5.7182, 5.719, 5.7181, 3.991, 3.991, 3.991,
            3.9907, 3.5997, 3.525, 3.4133, 3.4714, 3.4727, 3.4133, 3.525, 3.5997, 3.5997, 3.525,
            3.4122, 3.4738, 3.4728, 3.4109, 3.5259, 3.5997,
        ]
        # fmt: on
        for idx, _ in enumerate(self.lifepo4):
            cn_array.append(cnn.get_cn(self.lifepo4, idx, use_weights=True))

        assert_allclose(expected_array, cn_array, 2)

    def test_weighted_cn_no_oxid(self):
        cnn = CrystalNN(weighted_cn=True)
        cn_array = []
        # fmt: off
        expected_array = [
            5.8962, 5.8996, 5.8962, 5.8996, 5.7195, 5.7195, 5.7202, 5.7194, 4.0012, 4.0012,
            4.0012, 4.0009, 3.3897, 3.2589, 3.1218, 3.1914, 3.1914, 3.1218, 3.2589, 3.3897,
            3.3897, 3.2589, 3.1207, 3.1924, 3.1915, 3.1207, 3.2598, 3.3897,
        ]
        # fmt: on
        s = self.lifepo4.copy()
        s.remove_oxidation_states()
        for idx in range(len(s)):
            cn_array.append(cnn.get_cn(s, idx, use_weights=True))

        assert_allclose(expected_array, cn_array, 2)

    def test_fixed_length(self):
        cnn = CrystalNN(fingerprint_length=30)
        nn_data = cnn.get_nn_data(self.lifepo4, 0)
        assert len(nn_data.cn_weights) == 30
        assert len(nn_data.cn_nninfo) == 30

    def test_cation_anion(self):
        cnn = CrystalNN(weighted_cn=True, cation_anion=True)
        assert cnn.get_cn(self.lifepo4, 0, use_weights=True) == approx(5.8630, abs=1e-2)

    def test_x_diff_weight(self):
        cnn = CrystalNN(weighted_cn=True, x_diff_weight=0)
        assert cnn.get_cn(self.lifepo4, 0, use_weights=True) == approx(5.8630, abs=1e-2)

    def test_noble_gas_material(self):
        cnn = CrystalNN()

        assert cnn.get_cn(self.he_bcc, 0, use_weights=False) == 0

        cnn = CrystalNN(distance_cutoffs=(1.25, 5))
        assert cnn.get_cn(self.he_bcc, 0, use_weights=False) == 8

    def test_shifted_sites(self):
        cnn = CrystalNN()

        sites = [[0.0, 0.2, 0.2], [0, 0, 0]]
        struct = Structure([7, 0, 0, 0, 7, 0, 0, 0, 7], ["I"] * len(sites), sites)
        bonded_struct = cnn.get_bonded_structure(struct)

        sites_shifted = [[1.0, 0.2, 0.2], [0, 0, 0]]
        struct_shifted = Structure([7, 0, 0, 0, 7, 0, 0, 0, 7], ["I"] * len(sites_shifted), sites_shifted)
        bonded_struct_shifted = cnn.get_bonded_structure(struct_shifted)

        assert len(bonded_struct.get_connected_sites(0)) == len(bonded_struct_shifted.get_connected_sites(0))

    def test_get_cn(self):
        cnn = CrystalNN()

        site_0_coord_num = cnn.get_cn(self.disordered_struct, 0, on_disorder="take_max_species")
        site_0_coord_num_strict_majority = cnn.get_cn(
            self.disordered_struct_with_majority, 0, on_disorder="take_majority_strict"
        )
        assert site_0_coord_num == 8
        assert site_0_coord_num == site_0_coord_num_strict_majority

        with pytest.raises(
            ValueError, match="Site 0 has no majority species, the max species is Fe with occupancy 0.4"
        ):
            cnn.get_cn(self.disordered_struct, 0, on_disorder="take_majority_strict")
        with pytest.raises(
            ValueError, match="enerating StructureGraphs for disordered Structures is unsupported. Pass on_disorder="
        ):
            cnn.get_cn(self.disordered_struct, 0, on_disorder="error")

    def test_get_bonded_structure(self):
        cnn = CrystalNN()

        structure_graph = cnn.get_bonded_structure(self.disordered_struct, on_disorder="take_max_species")
        structure_graph_strict_majority = cnn.get_bonded_structure(
            self.disordered_struct_with_majority, on_disorder="take_majority_strict"
        )
        structure_graph_drop_majority = cnn.get_bonded_structure(
            self.disordered_struct_with_majority, on_disorder="take_majority_drop"
        )

        assert isinstance(structure_graph, StructureGraph)
        assert len(structure_graph) == 2
        assert structure_graph == structure_graph_strict_majority == structure_graph_drop_majority

        with pytest.raises(
            ValueError, match="Site 0 has no majority species, the max species is Fe with occupancy 0.4"
        ):
            cnn.get_bonded_structure(self.disordered_struct, 0, on_disorder="take_majority_strict")

        expected_msg = "Generating StructureGraphs for disordered Structures is unsupported. Pass on_disorder="
        with pytest.raises(ValueError, match=expected_msg):
            cnn.get_bonded_structure(self.disordered_struct, 0, on_disorder="error")


class TestCutOffDictNN(PymatgenTest):
    def setUp(self):
        self.diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264], [0, 0, 2.528]]),
            ["C", "C"],
            [[2.554, 1.806, 4.423], [0.365, 0.258, 0.632]],
            coords_are_cartesian=True,
        )

    def test_cn(self):
        nn = CutOffDictNN({("C", "C"): 2})
        assert nn.get_cn(self.diamond, 0) == 4

        nn_null = CutOffDictNN()
        assert nn_null.get_cn(self.diamond, 0) == 0

    def test_from_preset(self):
        nn = CutOffDictNN.from_preset("vesta_2019")
        assert nn.get_cn(self.diamond, 0) == 4

        # test error thrown on unknown preset
        with pytest.raises(ValueError, match="Unknown preset='test'"):
            CutOffDictNN.from_preset("test")


@unittest.skipIf(not which("critic2"), "critic2 executable not present")
class TestCritic2NN(PymatgenTest):
    def setUp(self):
        self.diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264], [0, 0, 2.528]]),
            ["C", "C"],
            [[2.554, 1.806, 4.423], [0.365, 0.258, 0.632]],
            coords_are_cartesian=True,
        )

    def test_cn(self):
        Critic2NN()
        # assert nn.get_cn(self.diamond, 0) == 4


class TestMetalEdgeExtender(PymatgenTest):
    def setUp(self):
        self.LiEC = Molecule.from_file(f"{test_dir}/LiEC.xyz")
        self.phsh = Molecule.from_file(f"{test_dir}/phsh.xyz")
        self.phsh_graph = MoleculeGraph.with_edges(
            molecule=self.phsh,
            edges={
                (0, 1): None,
                (0, 2): None,
                (0, 3): None,
                (0, 4): None,
                (4, 5): None,
                (4, 6): None,
                (4, 18): None,
                (6, 7): None,
                (6, 16): None,
                (7, 8): None,
                (7, 9): None,
                (9, 10): None,
                (9, 11): None,
                (11, 14): None,
                (12, 13): None,
                (12, 25): None,
                (14, 15): None,
                (14, 16): None,
                (16, 17): None,
                (18, 19): None,
                (18, 20): None,
                (18, 21): None,
                (21, 22): None,
                (21, 23): None,
                (21, 24): None,
            },
        )
        self.LiEC_graph = MoleculeGraph.with_edges(
            molecule=self.LiEC,
            edges={
                (0, 2): None,
                (0, 1): None,
                (1, 3): None,
                (1, 4): None,
                (2, 7): None,
                (2, 5): None,
                (2, 8): None,
                (3, 6): None,
                (4, 5): None,
                (5, 9): None,
                (5, 10): None,
            },
        )

        # potassium + 7 H2O. 4 at ~2.5 Ang and 3 more within 4.25 Ang
        uncharged_K_cluster = Molecule.from_file(f"{test_dir}/water_cluster_K.xyz")
        K_sites = [s.coords for s in uncharged_K_cluster]
        K_species = [s.species for s in uncharged_K_cluster]
        charged_K_cluster = Molecule(K_species, K_sites, charge=1)
        self.water_cluster_K = MoleculeGraph.with_empty_graph(charged_K_cluster)
        assert len(self.water_cluster_K.graph.edges) == 0

        # Mg + 6 H2O at 1.94 Ang from Mg
        uncharged_Mg_cluster = Molecule.from_file(f"{test_dir}/water_cluster_Mg.xyz")
        Mg_sites = [s.coords for s in uncharged_Mg_cluster]
        Mg_species = [s.species for s in uncharged_Mg_cluster]
        charged_Mg_cluster = Molecule(Mg_species, Mg_sites, charge=2)
        self.water_cluster_Mg = MoleculeGraph.with_empty_graph(charged_Mg_cluster)

    def test_metal_edge_extender(self):
        assert len(self.LiEC_graph.graph.edges) == 11
        extended_mol_graph = metal_edge_extender(self.LiEC_graph)
        assert len(extended_mol_graph.graph.edges) == 12

    def test_oxygen_edge_extender(self):
        assert len(self.phsh_graph.graph.edges) == 25
        phsh_fixed_graph = oxygen_edge_extender(self.phsh_graph)
        assert len(phsh_fixed_graph.graph.edges) == 26

    def test_custom_metals(self):
        extended_mol_graph = metal_edge_extender(self.LiEC_graph, metals={"K"})
        assert len(extended_mol_graph.graph.edges) == 11

        # empty metals should exit cleanly with no change to graph
        mol_graph = metal_edge_extender(self.water_cluster_K, metals={}, cutoff=2.5)
        assert len(mol_graph.graph.edges) == 0

        mol_graph = metal_edge_extender(self.water_cluster_K, metals={"K"}, cutoff=2.5)
        assert len(mol_graph.graph.edges) == 4

        extended_graph = metal_edge_extender(self.water_cluster_K, metals={"K"}, cutoff=4.5)
        assert len(extended_graph.graph.edges) == 7

        # if None, should auto-detect Li
        extended_mol_graph = metal_edge_extender(self.LiEC_graph, metals=None)
        assert len(extended_mol_graph.graph.edges) == 12

    def test_custom_coordinators(self):
        # leave out Oxygen, graph should not change
        extended_mol_graph = metal_edge_extender(self.LiEC_graph, coordinators={"N", "F", "S", "Cl"})
        assert len(extended_mol_graph.graph.edges) == 11
        # empty coordinators should exit cleanly with no change
        extended_mol_graph = metal_edge_extender(self.LiEC_graph, coordinators={})
        assert len(extended_mol_graph.graph.edges) == 11

    def test_custom_cutoff(self):
        short_mol_graph = metal_edge_extender(self.LiEC_graph, cutoff=0.5)
        assert len(short_mol_graph.graph.edges) == 11

        # with a cutoff of 1.5, no edges should be found.
        # test that the 2nd pass analysis (auto increasing cutoff to 2.5) picks
        # up the six coordination bonds
        short_mol_graph = metal_edge_extender(self.water_cluster_Mg, cutoff=1.5)
        assert len(short_mol_graph.graph.edges) == 6
