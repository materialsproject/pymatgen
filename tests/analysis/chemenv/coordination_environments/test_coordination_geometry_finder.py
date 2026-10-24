from __future__ import annotations

import copy
import os
import random

import numpy as np
import orjson
import pytest
from numpy.testing import assert_allclose
from pytest import approx
from scipy.spatial.transform import Rotation

from pymatgen.analysis.chemenv.coordination_environments import coordination_geometry_finder as cgf_module
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import (
    SimpleAbundanceChemenvStrategy,
    SimplestChemenvStrategy,
)
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import (
    AllCoordinationGeometries,
    CoordinationGeometry,
    ExplicitPermutationsAlgorithm,
)
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import (
    AbstractGeometry,
    LocalGeometryFinder,
    find_rotation,
    find_scaling_factor,
    symmetry_measure,
)
from pymatgen.core import Species
from pymatgen.core.structure import Lattice, Structure
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

__author__ = "waroquiers"

json_dir = f"{TEST_FILES_DIR}/analysis/chemenv/json"


class TestCoordinationGeometryFinder(MatSciTest):
    def setup_method(self):
        self.lgf = LocalGeometryFinder()
        self.lgf.setup_parameters(
            centering_type="standard",
            structure_refinement=self.lgf.STRUCTURE_REFINEMENT_NONE,
        )

    # self.strategies = [SimplestChemenvStrategy(), SimpleAbundanceChemenvStrategy()]

    def test_abstract_geometry(self):
        cg_ts3 = self.lgf.allcg["TS:3"]
        cg_tet = self.lgf.allcg["T:4"]
        abstract_geom = AbstractGeometry.from_cg(cg=cg_ts3, centering_type="central_site")
        assert_allclose(abstract_geom.centre, [0.0, 0.0, 0.0])
        abstract_geom = AbstractGeometry.from_cg(cg=cg_ts3, centering_type="centroid")
        assert_allclose(abstract_geom.centre, [0.0, 0.0, 0.33333333333])
        with pytest.raises(
            ValueError,
            match="The center is the central site, no calculation of the centroid, "
            "variable include_central_site_in_centroid should be set to False",
        ):
            AbstractGeometry.from_cg(
                cg=cg_ts3,
                centering_type="central_site",
                include_central_site_in_centroid=True,
            )
        abstract_geom = AbstractGeometry.from_cg(
            cg=cg_ts3, centering_type="centroid", include_central_site_in_centroid=True
        )
        assert_allclose(abstract_geom.centre, [0.0, 0.0, 0.25])

        # WHY ARE WE TESTING STRINGS????
        # assert (
        #     str(abstract_geom) == "\nAbstract Geometry with 3 points :\n"
        #     "  [-1.    0.   -0.25]\n  [ 1.    0.   -0.25]\n  [ 0.   0.   0.75]\n"
        #     "Points are referenced to the centroid (calculated with the central site) :\n"
        #     "  [ 0.   0.   0.25]\n"
        # )

        symm_dict = symmetry_measure([[0.0, 0.0, 0.0]], [1.1, 2.2, 3.3])
        assert symm_dict["symmetry_measure"] == approx(0.0)
        assert symm_dict["scaling_factor"] is None
        assert symm_dict["rotation_matrix"] is None

        tio2_struct = self.get_structure("TiO2")

        envs = self.lgf.compute_coordination_environments(structure=tio2_struct, indices=[0])
        assert envs[0][0]["csm"] == approx(1.5309987846957258)
        assert envs[0][0]["ce_fraction"] == approx(1.0)
        assert envs[0][0]["ce_symbol"] == "O:6"
        assert sorted(envs[0][0]["permutation"]) == sorted([0, 4, 1, 5, 2, 3])

        self.lgf.setup_random_structure(coordination=5)
        assert len(self.lgf.structure) == 6

        self.lgf.setup_random_indices_local_geometry(coordination=5)
        assert self.lgf.icentral_site == 0
        assert len(self.lgf.indices) == 5

        self.lgf.setup_ordered_indices_local_geometry(coordination=5)
        assert self.lgf.icentral_site == 0
        assert self.lgf.indices == list(range(1, 6))

        self.lgf.setup_explicit_indices_local_geometry(explicit_indices=[3, 5, 2, 0, 1, 4])
        assert self.lgf.icentral_site == 0
        assert self.lgf.indices == [4, 6, 3, 1, 2, 5]

        LiFePO4_struct = self.get_structure("LiFePO4")
        site_idx = 10
        envs_LiFePO4 = self.lgf.compute_coordination_environments(structure=LiFePO4_struct, indices=[site_idx])
        assert envs_LiFePO4[site_idx][0]["csm"] == approx(0.140355832317)
        nbs_coords = [
            np.array([6.16700437, -4.55194317, -5.89031356]),
            np.array([4.71588167, -4.54248093, -3.75553856]),
            np.array([6.88012571, -5.79877503, -3.73177541]),
            np.array([6.90041188, -3.32797839, -3.71812416]),
        ]

        # test to check that one can pass voronoi_distance_cutoff
        struct = Structure(
            Lattice.cubic(25),
            ["O", "C", "O"],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 1.17], [0.0, 0.0, 2.34]],
        )
        self.lgf.setup_structure(structure=struct)
        self.lgf.compute_structure_environments(voronoi_distance_cutoff=25)

        self.lgf.setup_structure(LiFePO4_struct)
        self.lgf.setup_local_geometry(site_idx, coords=nbs_coords)

        perfect_tet = AbstractGeometry.from_cg(
            cg=cg_tet, centering_type="centroid", include_central_site_in_centroid=False
        )
        points_perfect_tet = perfect_tet.points_wcs_ctwcc()
        result = self.lgf.coordination_geometry_symmetry_measures_fallback_random(
            coordination_geometry=cg_tet, n_random=5, points_perfect=points_perfect_tet
        )
        (
            permutations_symmetry_measures,
            _permutations,
            _algos,
            _local2perfect_maps,
            _perfect2local_maps,
        ) = result
        for perm_csm_dict in permutations_symmetry_measures:
            assert perm_csm_dict["symmetry_measure"] == approx(0.140355832317)

    def _strategy_test(self, strategy):
        files = []
        for _dirpath, _dirnames, filenames in os.walk(json_dir):
            files.extend(filenames)
            break

        for json_file in files:
            with self.subTest(json_file=json_file):
                with open(f"{json_dir}/{json_file}", "rb") as file:
                    dct = orjson.loads(file.read())

                atom_indices = dct["atom_indices"]
                expected_geoms = dct["expected_geoms"]

                struct = Structure.from_dict(dct["structure"])

                struct = self.lgf.setup_structure(struct)
                se = self.lgf.compute_structure_environments_detailed_voronoi(
                    only_indices=atom_indices, maximum_distance_factor=1.5
                )

                # All strategies should get the correct environment with their default parameters
                strategy.set_structure_environments(se)
                for ienv, isite in enumerate(atom_indices):
                    ce = strategy.get_site_coordination_environment(struct[isite])
                    try:
                        coord_env = ce[0]
                    except TypeError:
                        coord_env = ce
                    # Check that the environment found is the expected one
                    assert coord_env == expected_geoms[ienv]

    @pytest.mark.xfail(reason="TODO: need someone to fix this")
    def test_simplest_chemenv_strategy(self):
        strategy = SimplestChemenvStrategy()
        self._strategy_test(strategy)

    @pytest.mark.xfail(reason="TODO: need someone to fix this")
    def test_simple_abundance_chemenv_strategy(self):
        strategy = SimpleAbundanceChemenvStrategy()
        self._strategy_test(strategy)

    def test_perfect_environments(self):
        allcg = AllCoordinationGeometries()
        indices_CN = {
            1: [0],
            2: [1, 0],
            3: [1, 0, 2],
            4: [2, 0, 3, 1],
            5: [2, 3, 1, 0, 4],
            6: [0, 2, 3, 1, 5, 4],
            7: [2, 6, 0, 3, 4, 5, 1],
            8: [1, 2, 6, 3, 7, 0, 4, 5],
            9: [5, 2, 6, 0, 4, 7, 3, 8, 1],
            10: [8, 5, 6, 3, 0, 7, 2, 4, 9, 1],
            11: [7, 6, 4, 1, 2, 5, 0, 8, 9, 10, 3],
            12: [5, 8, 9, 0, 3, 1, 4, 2, 6, 11, 10, 7],
            13: [4, 11, 5, 12, 1, 2, 8, 3, 0, 6, 9, 7, 10],
            20: [8, 12, 11, 0, 14, 10, 13, 6, 18, 1, 9, 17, 3, 19, 5, 7, 15, 2, 16, 4],
        }

        for coordination in range(1, 21):
            for mp_symbol in allcg.get_implemented_geometries(coordination=coordination, returned="mp_symbol"):
                cg = allcg.get_geometry_from_mp_symbol(mp_symbol=mp_symbol)
                self.lgf.allcg = AllCoordinationGeometries(only_symbols=[mp_symbol])
                self.lgf.setup_test_perfect_environment(
                    mp_symbol,
                    randomness=False,
                    indices=indices_CN[coordination],
                    random_translation="NONE",
                    random_rotation="NONE",
                    random_scale="NONE",
                )
                se = self.lgf.compute_structure_environments(
                    only_indices=[0],
                    maximum_distance_factor=1.01 * cg.distfactor_max,
                    min_cn=cg.coordination_number,
                    max_cn=cg.coordination_number,
                    only_symbols=[mp_symbol],
                )
                assert abs(se.get_csm(0, mp_symbol)["symmetry_measure"] - 0.0) < 1e-8, (
                    f"Failed to get perfect environment with {mp_symbol=}"
                )

    def test_disable_hints(self):
        allcg = AllCoordinationGeometries()
        mp_symbol = "SH:13"
        mp_symbols = ["SH:13", "HP:12"]
        cg = allcg.get_geometry_from_mp_symbol(mp_symbol=mp_symbol)
        cg_points = cg.points
        cg_points[-1] = [0.9 * cc for cc in cg_points[-1]]
        self.lgf.allcg = AllCoordinationGeometries(only_symbols=[mp_symbol])
        self.lgf.setup_test_perfect_environment(
            mp_symbol,
            randomness=False,
            indices=[4, 11, 5, 12, 1, 2, 8, 3, 0, 6, 9, 7, 10],
            random_translation="NONE",
            random_rotation="NONE",
            random_scale="NONE",
            points=cg_points,
        )
        se_nohints = self.lgf.compute_structure_environments(
            only_indices=[0],
            maximum_distance_factor=1.02 * cg.distfactor_max,
            min_cn=12,
            max_cn=13,
            only_symbols=mp_symbols,
            get_from_hints=False,
        )
        se_hints = self.lgf.compute_structure_environments(
            only_indices=[0],
            maximum_distance_factor=1.02 * cg.distfactor_max,
            min_cn=12,
            max_cn=13,
            only_symbols=mp_symbols,
            get_from_hints=True,
        )
        with pytest.raises(KeyError, match="12"):
            se_nohints.ce_list[0][12]
        assert se_hints.ce_list[0][13][0] == se_nohints.ce_list[0][13][0]
        assert set(se_nohints.ce_list[0]).issubset(set(se_hints.ce_list[0]))


@pytest.fixture(scope="module")
def nacl() -> Structure:
    """Rock-salt NaCl: every Na sits in a perfect octahedron of Cl."""
    return Structure.from_spacegroup("Fm-3m", Lattice.cubic(5.69), ["Na", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])


@pytest.fixture
def lgf() -> LocalGeometryFinder:
    finder = LocalGeometryFinder()
    finder.setup_parameters(
        centering_type="standard", structure_refinement=LocalGeometryFinder.STRUCTURE_REFINEMENT_NONE
    )
    return finder


@pytest.fixture
def seeded_rng(monkeypatch):
    """Make the unseeded RNGs used by the finder deterministic."""
    real_default_rng = np.random.default_rng
    monkeypatch.setattr(np.random, "default_rng", lambda *_args, **_kwargs: real_default_rng(42))
    monkeypatch.setattr(cgf_module, "shuffle", random.Random(42).shuffle)


SQUARE = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]]
OCTAHEDRON = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]


class TestAbstractGeometry:
    @pytest.mark.parametrize(
        ("bare_coords", "centering_type", "include_central", "central_site", "match"),
        [
            (SQUARE, "standard", True, [0, 0, 0], "include_central_site_in_centroid should be set to False"),
            (SQUARE, "standard", False, None, "the central site should be given"),
            (OCTAHEDRON, "standard", True, None, "no central site is given"),
            (OCTAHEDRON, "central_site", False, None, "the central site should be given"),
            (OCTAHEDRON, "centroid", True, None, "no central site is given"),
            (OCTAHEDRON, "unknown", False, [0, 0, 0], "Unknown centering_type 'unknown'"),
        ],
    )
    def test_invalid_centering(self, bare_coords, centering_type, include_central, central_site, match):
        with pytest.raises(ValueError, match=match):
            AbstractGeometry(
                central_site=central_site,
                bare_coords=bare_coords,
                centering_type=centering_type,
                include_central_site_in_centroid=include_central,
            )

    @pytest.mark.parametrize(
        ("bare_coords", "centering_type", "include_central", "expected_centre", "expected_str"),
        [
            # CN < 5 with standard centering: the reference point is the central site
            (SQUARE, "standard", False, [0.7, 0.0, 0.0], "to the centroid (calculated without the central site)"),
            # CN >= 5 with standard centering: the reference point is the centroid, with or without the central site
            (OCTAHEDRON, "standard", False, [0.0, 0.0, 0.0], "to the centroid (calculated without the central site)"),
            (OCTAHEDRON, "standard", True, [0.1, 0.0, 0.0], "to the centroid (calculated with the central site)"),
            (OCTAHEDRON, "central_site", False, [0.7, 0.0, 0.0], "Points are referenced to the central site"),
            (OCTAHEDRON, "centroid", False, [0.0, 0.0, 0.0], "(calculated without the central site) :"),
            (OCTAHEDRON, "centroid", True, [0.1, 0.0, 0.0], "(calculated with the central site) :"),
        ],
    )
    def test_centering(self, bare_coords, centering_type, include_central, expected_centre, expected_str):
        geom = AbstractGeometry(
            central_site=[0.7, 0.0, 0.0],
            bare_coords=bare_coords,
            centering_type=centering_type,
            include_central_site_in_centroid=include_central,
        )
        assert_allclose(geom.centre, expected_centre, atol=1e-12)
        assert_allclose(geom.coords, np.array(bare_coords) - expected_centre, atol=1e-12)
        assert_allclose(geom.central_site, np.array([0.7, 0.0, 0.0]) - expected_centre, atol=1e-12)
        assert geom.cn == geom.coordination_number == len(bare_coords)
        string = str(geom)
        assert string.startswith(f"\nAbstract Geometry with {len(bare_coords)} points :")
        assert expected_str in string

    def test_no_central_site(self):
        geom = AbstractGeometry(bare_coords=[[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0], [2.0, 2.0, 2.0], [1, 1, 1]])
        assert_allclose(geom.centre, [1.0, 1.0, 1.0])
        # Without a central site, the "central site" is the origin seen from the centre
        assert_allclose(geom.central_site, [-1.0, -1.0, -1.0])
        assert geom.bare_centre is None
        assert_allclose(geom.points_wcs_csc(), geom.bare_points_with_centre)

    def test_points_with_permutation(self):
        geom = AbstractGeometry(central_site=[0.5, 0.0, 0.0], bare_coords=OCTAHEDRON, centering_type="centroid")
        perm = [5, 4, 3, 2, 1, 0]
        octahedron = np.array(OCTAHEDRON)
        assert_allclose(geom.centroid_without_centre, [0.0, 0.0, 0.0])
        assert_allclose(geom.centroid_with_centre, [0.5 / 7, 0.0, 0.0])
        for method, centre in (
            (geom.points_wocs_csc, [0.5, 0, 0]),
            (geom.points_wocs_ctwcc, [0.5 / 7, 0, 0]),
            (geom.points_wocs_ctwocc, [0.0, 0, 0]),
        ):
            assert_allclose(method(), octahedron - centre)
            assert_allclose(method(permutation=perm), octahedron[perm] - centre)
        for method, centre in (
            (geom.points_wcs_csc, [0.5, 0, 0]),
            (geom.points_wcs_ctwcc, [0.5 / 7, 0, 0]),
            (geom.points_wcs_ctwocc, [0.0, 0, 0]),
        ):
            pts = method(permutation=perm)
            assert pts.shape == (7, 3)
            assert_allclose(pts[0], np.array([0.5, 0, 0]) - centre)
            assert_allclose(pts[1:], octahedron[perm] - centre)
            assert_allclose(method()[1:], octahedron - centre)


class TestSymmetryMeasure:
    def test_perfect_rotated_scaled(self):
        rot = Rotation.from_euler("xyz", [10, 20, 30], degrees=True).as_matrix()
        perfect = np.array(OCTAHEDRON)
        distorted = 1.7 * perfect @ rot.T
        sm = symmetry_measure(points_distorted=distorted, points_perfect=perfect)
        assert sm["symmetry_measure"] == approx(0.0, abs=1e-12)
        assert sm["scaling_factor"] == approx(1 / 1.7)
        assert_allclose(sm["rotation_matrix"], rot.T, atol=1e-10)

    def test_known_distortion(self):
        # A 2 x 1 rectangle vs. a unit square: the optimal scaling factor is 6/10 and the CSM is exactly 10
        perfect = np.array(SQUARE)
        distorted = perfect * [2.0, 1.0, 1.0]
        rot = find_rotation(points_distorted=distorted, points_perfect=perfect)
        scale, rotated, pts_perfect = find_scaling_factor(points_distorted=distorted, points_perfect=perfect, rot=rot)
        assert scale == approx(0.6)
        assert_allclose(rotated, distorted @ rot.T)
        assert pts_perfect is perfect
        sm = symmetry_measure(points_distorted=distorted, points_perfect=perfect)
        assert sm["symmetry_measure"] == approx(10.0)
        assert sm["scaling_factor"] == approx(0.6)


class TestLocalGeometryFinderSetup:
    def test_setup_parameters(self, lgf):
        lgf.setup_parameters(
            centering_type="centroid",
            include_central_site_in_centroid=True,
            bva_distance_scale_factor=LocalGeometryFinder.BVA_DISTANCE_SCALE_FACTORS["GGA_relaxed"],
            structure_refinement=LocalGeometryFinder.STRUCTURE_REFINEMENT_SYMMETRIZED,
            spg_analyzer_options={"symprec": 0.1, "angle_tolerance": 2},
        )
        assert lgf.centering_type == "centroid"
        assert lgf.include_central_site_in_centroid
        assert lgf.bva_distance_scale_factor == approx(1.015)
        assert lgf.spg_analyzer_options == {"symprec": 0.1, "angle_tolerance": 2}
        lgf.setup_parameter("centering_type", "central_site")
        assert lgf.centering_type == "central_site"

    @pytest.mark.parametrize(
        "refinement",
        [
            LocalGeometryFinder.STRUCTURE_REFINEMENT_NONE,
            LocalGeometryFinder.STRUCTURE_REFINEMENT_REFINED,
            LocalGeometryFinder.STRUCTURE_REFINEMENT_SYMMETRIZED,
        ],
    )
    # spglib itself emits this when called through SpacegroupAnalyzer
    @pytest.mark.filterwarnings("ignore:Set OLD_ERROR_HANDLING to false:DeprecationWarning")
    def test_setup_structure_refinement(self, lgf, nacl, refinement):
        primitive = nacl.get_primitive_structure()
        assert len(primitive) == 2
        lgf.setup_parameters(structure_refinement=refinement)
        lgf.setup_structure(primitive)
        assert lgf.initial_structure == primitive
        structure = lgf.get_structure()
        if refinement == LocalGeometryFinder.STRUCTURE_REFINEMENT_NONE:
            assert structure == primitive
            assert lgf.spg_analyzer is None
        else:
            # refinement gives back the conventional cubic cell
            assert len(structure) == 8
            assert lgf.spg_analyzer.get_space_group_symbol() == "Fm-3m"
        if refinement == LocalGeometryFinder.STRUCTURE_REFINEMENT_SYMMETRIZED:
            assert len(lgf.symmetrized_structure.equivalent_indices) == 2
        else:
            assert lgf.symmetrized_structure is None

    def test_set_structure(self, lgf):
        lgf.set_structure(
            lattice=Lattice.cubic(4.0),
            species=["Si", "Si"],
            coords=[[0, 0, 0], [1.0, 1.0, 1.0]],
            coords_are_cartesian=True,
        )
        assert_allclose(lgf.structure.frac_coords[1], [0.25, 0.25, 0.25])
        assert lgf.structure.charge == 0

    def test_setup_random_structure(self, lgf, seeded_rng):
        lgf.setup_random_structure(coordination=4)
        assert lgf.structure.composition.reduced_formula == "Si"
        assert len(lgf.structure) == 5
        assert sorted(lgf.indices) == [1, 2, 3, 4]

    @pytest.mark.parametrize(
        ("symbol", "symbol_type"),
        [("O:6", "mp_symbol"), ("O:6", "MP"), ("OC-6", "IUPAC"), (None, "CoordinationGeometry")],
    )
    def test_setup_test_perfect_environment_symbol_types(self, lgf, symbol, symbol_type):
        if symbol_type == "CoordinationGeometry":
            symbol = lgf.allcg["O:6"]
        lgf.setup_test_perfect_environment(symbol, symbol_type=symbol_type, indices="ORDERED")
        assert lgf.structure.formula == "Cu1 O6"
        assert lgf.perfect_geometry.cn == 6
        assert_allclose(lgf.local_geometry.bare_coords, lgf.allcg["O:6"].points)

    def test_setup_test_perfect_environment_wrong_type(self, lgf):
        with pytest.raises(ValueError, match="Wrong mp_symbol to setup coordination geometry"):
            lgf.setup_test_perfect_environment("O:6", symbol_type="name")

    def test_setup_test_perfect_environment_random(self, lgf, seeded_rng):
        lgf.setup_test_perfect_environment(
            "T:4",
            randomness=True,
            max_random_dist=0.05,
            indices="RANDOM",
            random_translation="RANDOM",
            random_rotation="RANDOM",
            random_scale="RANDOM",
        )
        assert lgf.structure.formula == "Cu1 O4"
        # The randomly displaced, scaled, rotated and translated tetrahedron is only slightly distorted
        result = lgf.get_coordination_symmetry_measures()
        assert set(result) == {"T:4", "S:4", "SY:4", "SS:4"}
        assert 0 < result["T:4"]["csm"] < 1.0
        assert min(result, key=lambda symb: result[symb]["csm"]) == "T:4"

    def test_setup_test_perfect_environment_explicit_transformations(self, lgf):
        rotation = Rotation.from_euler("z", 30, degrees=True).as_matrix()
        lgf.setup_test_perfect_environment(
            "T:4",
            indices=[3, 2, 1, 0],
            random_translation=[1.0, 2.0, 3.0],
            random_rotation=rotation.tolist(),
            random_scale=2.0,
        )
        points = np.array(lgf.allcg["T:4"].points)[[3, 2, 1, 0]]
        assert_allclose(lgf.structure.cart_coords[0], [1.0, 2.0, 3.0])
        assert_allclose(lgf.structure.cart_coords[1:], 2.0 * points @ rotation.T + [1.0, 2.0, 3.0])
        assert lgf.get_coordination_symmetry_measures()["T:4"]["csm"] == approx(0.0, abs=1e-10)


class TestLocalGeometryFinderMeasures:
    def test_cn1_without_geometries(self):
        finder = LocalGeometryFinder(only_symbols=["O:6"])
        finder.setup_parameters(centering_type="standard", structure_refinement="none")
        finder.setup_structure(Structure(Lattice.cubic(10), ["Cu", "O"], [[0, 0, 0], [0.2, 0, 0]]))
        finder.setup_local_geometry(0, coords=[[2.0, 0.0, 0.0]])
        assert finder.get_coordination_symmetry_measures() == {}

    def test_cn1(self, lgf):
        lgf.setup_test_perfect_environment("S:1", indices="ORDERED")
        result = lgf.get_coordination_symmetry_measures(all_csms=False)
        assert result["S:1"]["csm"] == 0.0
        assert "csm_wcs_ctwcc" not in result["S:1"]
        result = lgf.get_coordination_symmetry_measures(all_csms=True)
        assert result["S:1"]["csm_wcs_ctwcc"] == 0.0

    def test_not_only_minimum(self, lgf):
        lgf.setup_test_perfect_environment("T:4", indices="ORDERED")
        result = lgf.get_coordination_symmetry_measures(only_minimum=False)
        csms = result["T:4"]["csm"]
        assert len(csms) == len(result["T:4"]["indices"]) == len(lgf.allcg["T:4"].algorithms[0].permutations)
        assert min(sm["symmetry_measure"] for sm in csms) == approx(0.0, abs=1e-10)
        for l2p, p2l in zip(result["T:4"]["local2perfect_map"], result["T:4"]["perfect2local_map"], strict=True):
            assert all(l2p[p2l[ii]] == ii for ii in range(4))

    @pytest.mark.parametrize("optimization", [None, 2])
    def test_explicit_permutations_standard(self, lgf, optimization):
        lgf.setup_test_perfect_environment("T:4", indices=[2, 0, 3, 1])
        cg = lgf.allcg["T:4"]
        points_perfect = AbstractGeometry.from_cg(cg).points_wcs_ctwcc()
        csms, perms, algos, l2p_maps, p2l_maps = lgf.coordination_geometry_symmetry_measures_standard(
            cg, cg.algorithms[0], points_perfect=points_perfect, optimization=optimization
        )
        n_perms = len(cg.algorithms[0].permutations)
        assert len(csms) == len(perms) == len(algos) == len(l2p_maps) == len(p2l_maps) == n_perms
        assert set(algos) == {str(cg.algorithms[0])}
        assert min(sm["symmetry_measure"] for sm in csms) == approx(0.0, abs=1e-10)

    def test_optim_mixed_algorithms(self):
        finder = LocalGeometryFinder(only_symbols=["O:6"])
        finder.setup_parameters(centering_type="standard", structure_refinement="none")
        oct_dict = finder.allcg["O:6"].as_dict()
        oct_dict.update(mp_symbol="XO:6", _algorithms=[ExplicitPermutationsAlgorithm([list(range(6))]).as_dict()])
        finder.allcg.cg_list.append(CoordinationGeometry.from_dict(oct_dict))
        finder.setup_test_perfect_environment("O:6", indices="ORDERED")
        with pytest.raises(ValueError, match="All algorithms should be EXPLICIT_PERMUTATIONS or SEPARATION_PLANE"):
            finder.get_coordination_symmetry_measures_optim(nb_set=None)

    def test_permutations_safe_override(self):
        finder = LocalGeometryFinder(permutations_safe_override=True, only_symbols=["O:6"])
        assert finder.allcg["O:6"].permutations_safe_override
        finder.setup_parameters(centering_type="standard", structure_refinement="none")
        finder.setup_test_perfect_environment("O:6", indices="ORDERED")
        with pytest.raises(ValueError, match="No permutations safe override anymore"):
            finder.get_coordination_symmetry_measures()

    @pytest.mark.parametrize("plane_safe_permutations", [False, True])
    def test_separation_plane(self, plane_safe_permutations):
        finder = LocalGeometryFinder(
            only_symbols=["O:6", "T:6", "PP:6"], plane_safe_permutations=plane_safe_permutations
        )
        finder.setup_parameters(centering_type="standard", structure_refinement="none")
        finder.setup_test_perfect_environment("O:6", indices=[0, 2, 3, 1, 5, 4])
        result = finder.get_coordination_symmetry_measures()
        assert set(result) == {"O:6", "T:6", "PP:6"}
        assert result["O:6"]["csm"] == approx(0.0, abs=1e-10)
        assert result["T:6"]["csm"] > 1.0
        assert result["O:6"]["algo"] == "SEPARATION_PLANE"
        assert result["O:6"]["csm_wocs_ctwocc"] == approx(0.0, abs=1e-10)

        cg = finder.allcg["O:6"]
        points_perfect = AbstractGeometry.from_cg(cg).points_wcs_ctwcc()
        csms, perms, *_ = finder.coordination_geometry_symmetry_measures(
            cg, tested_permutations=True, points_perfect=points_perfect
        )
        csms_all, perms_all, *_ = finder.coordination_geometry_symmetry_measures(cg, points_perfect=points_perfect)
        # equivalent permutations are skipped, hence fewer permutations are tested
        assert 0 < len(perms) < len(perms_all)
        assert min(sm["symmetry_measure"] for sm in csms) == approx(0.0, abs=1e-10)
        assert min(sm["symmetry_measure"] for sm in csms_all) == approx(0.0, abs=1e-10)

        csms_t, perms_t, sep_perms = finder.coordination_geometry_symmetry_measures_separation_plane(
            cg, cg.algorithms[0], testing=True, points_perfect=points_perfect
        )
        assert len(csms_t) == len(perms_t) == len(sep_perms) > 0

    def test_fallback_random(self, lgf, seeded_rng):
        lgf.setup_test_perfect_environment("O:6", indices="ORDERED")
        cg = lgf.allcg["O:6"]
        points_perfect = AbstractGeometry.from_cg(cg).points_wcs_ctwcc()
        with pytest.warns(DeprecationWarning, match="NRANDOM is deprecated, use n_random instead"):
            csms, perms, algos, l2p_maps, p2l_maps = lgf.coordination_geometry_symmetry_measures_fallback_random(
                cg, points_perfect=points_perfect, NRANDOM=3
            )
        assert len(csms) == len(perms) == len(algos) == len(l2p_maps) == len(p2l_maps) == 3
        for perm, l2p in zip(perms, l2p_maps, strict=True):
            assert sorted(perm) == list(range(6))
            assert all(l2p[pp] == ii for ii, pp in enumerate(perm))

    def test_separation_plane_optim_invalid(self, lgf):
        lgf.setup_test_perfect_environment("O:6", indices="ORDERED")
        cg = lgf.allcg["O:6"]
        with pytest.raises(ValueError, match="Optimization should be 1 or 2"):
            lgf.coordination_geometry_symmetry_measures_separation_plane_optim(
                cg, cg.algorithms[0], nb_set=None, optimization=3
            )


class TestComputeStructureEnvironments:
    @pytest.mark.parametrize("optimization", [None, 1, 2])
    @pytest.mark.parametrize("plane_safe_permutations", [False, True])
    def test_optimization_levels(self, nacl, optimization, plane_safe_permutations):
        finder = LocalGeometryFinder(plane_safe_permutations=plane_safe_permutations)
        finder.setup_parameters(centering_type="standard", structure_refinement="none")
        finder.setup_structure(nacl)
        se = finder.compute_structure_environments(
            only_indices=[0], maximum_distance_factor=1.2, min_cn=6, max_cn=6, optimization=optimization
        )
        csm = se.get_csm(0, "O:6")
        assert csm["symmetry_measure"] == approx(0.0, abs=1e-10)
        assert csm["other_symmetry_measures"]["scaling_factor_wcs_csc"] == approx(5.69 / 2)
        assert se.get_csm(0, "T:6")["symmetry_measure"] > 1.0

    # Geometries whose separation plane algorithms have unordered planes (PP:6), asymmetric point groups
    # (SA:8) or ordered point groups on both sides of the plane (TT_1:9)
    @pytest.mark.parametrize(
        ("mp_symbol", "indices"),
        [("PP:6", [0, 2, 3, 1, 5, 4]), ("SA:8", [1, 2, 6, 3, 7, 0, 4, 5]), ("TT_1:9", [5, 2, 6, 0, 4, 7, 3, 8, 1])],
    )
    @pytest.mark.parametrize("optimization", [None, 1])
    def test_perfect_environments_optimization(self, lgf, mp_symbol, indices, optimization):
        cg = lgf.allcg[mp_symbol]
        lgf.allcg = AllCoordinationGeometries(only_symbols=[mp_symbol])
        lgf.setup_test_perfect_environment(mp_symbol, indices=indices)
        se = lgf.compute_structure_environments(
            only_indices=[0],
            maximum_distance_factor=1.01 * cg.distfactor_max,
            min_cn=cg.coordination_number,
            max_cn=cg.coordination_number,
            only_symbols=[mp_symbol],
            optimization=optimization,
        )
        assert se.get_csm(0, mp_symbol)["symmetry_measure"] == approx(0.0, abs=1e-8)

    def test_site_selection(self, lgf, nacl):
        lgf.setup_structure(nacl)
        kwargs = {"maximum_distance_factor": 1.2, "min_cn": 6, "max_cn": 6, "only_symbols": ["O:6"]}
        se = lgf.compute_structure_environments(only_atoms=["Cl"], **kwargs)
        assert [isite for isite, ce in enumerate(se.ce_list) if ce is not None] == [4, 5, 6, 7]
        se = lgf.compute_structure_environments(excluded_atoms=["Cl"], only_indices=[1, 2, 5], **kwargs)
        assert [isite for isite, ce in enumerate(se.ce_list) if ce is not None] == [1, 2]
        assert lgf.valences == "undefined"

    def test_valences_from_oxidation_states(self, lgf, nacl):
        lgf.setup_structure(nacl.copy().add_oxidation_state_by_element({"Na": 1, "Cl": -1}))
        se = lgf.compute_structure_environments(maximum_distance_factor=1.2, min_cn=6, max_cn=6, only_symbols=["O:6"])
        assert lgf.valences == [1] * 4 + [-1] * 4
        # only cations by default
        assert [isite for isite, ce in enumerate(se.ce_list) if ce is not None] == [0, 1, 2, 3]

    def test_valences_species_without_oxidation_state(self, lgf, nacl):
        structure = nacl.copy()
        structure.replace_species({"Na": Species("Na", None), "Cl": Species("Cl", None)})
        lgf.setup_structure(structure)
        lgf.compute_structure_environments(only_indices=[0], maximum_distance_factor=1.2, min_cn=6, max_cn=6)
        assert lgf.valences == "undefined"

    def test_recompute_and_timelimit(self, lgf, nacl):
        lgf.setup_structure(nacl)
        kwargs = {"maximum_distance_factor": 1.2, "only_symbols": ["O:6"]}
        # O:6 has no neighbors sets hints, so get_from_hints does not add anything
        se = lgf.compute_structure_environments(only_indices=[0], recompute={"cns": [6]}, get_from_hints=True, **kwargs)
        assert se.get_csm(0, "O:6")["symmetry_measure"] == approx(0.0, abs=1e-10)
        assert all(
            src["origin"] != "nb_set_hints"
            for nb_sets in se.neighbors_sets[0].values()
            for nb_set in nb_sets
            for src in nb_set.sources
        )
        # With no time left, only the first site is computed
        se = lgf.compute_structure_environments(timelimit=0.0, **kwargs)
        assert se.ce_list[0] is not None
        assert all(ce is None for ce in se.ce_list[1:])

    def test_initial_structure_environments(self, lgf, nacl):
        lgf.setup_structure(nacl)
        kwargs = {"maximum_distance_factor": 1.2, "min_cn": 6, "max_cn": 6, "only_symbols": ["O:6"]}
        se = lgf.compute_structure_environments(only_indices=[0], **kwargs)
        se2 = lgf.compute_structure_environments(
            only_indices=[0], initial_structure_environments=se, info={"origin": "test"}, **kwargs
        )
        assert se2 is se
        assert se.info["origin"] == "test"
        assert se.get_csm(0, "O:6")["symmetry_measure"] == approx(0.0, abs=1e-10)

        with pytest.raises(ValueError, match="Detailed Voronoi is not the same in initial_structure_environments"):
            lgf.compute_structure_environments(
                initial_structure_environments=se, voronoi_normalized_distance_tolerance=0.2, **kwargs
            )

        # A Voronoi container that only differs within tolerance is accepted and reused
        se_close = copy.deepcopy(se)
        se_close.voronoi.voronoi_list2[0][0]["distance"] += 1e-12
        se3 = lgf.compute_structure_environments(only_indices=[0], initial_structure_environments=se_close, **kwargs)
        assert se3 is se_close
        assert lgf.detailed_voronoi is se_close.voronoi

        lgf.setup_structure(nacl.copy().perturb(0.01, min_distance=0.01, seed=0))
        with pytest.raises(ValueError, match="Structure is not the same in initial_structure_environments"):
            lgf.compute_structure_environments(initial_structure_environments=se, **kwargs)

    def test_hints_out_of_cn_range(self, lgf):
        mp_symbol = "SH:13"
        cg = lgf.allcg[mp_symbol]
        points = [list(pp) for pp in cg.points]
        points[-1] = [0.9 * cc for cc in points[-1]]
        lgf.allcg = AllCoordinationGeometries(only_symbols=[mp_symbol])
        lgf.setup_test_perfect_environment(mp_symbol, indices=[4, 11, 5, 12, 1, 2, 8, 3, 0, 6, 9, 7, 10], points=points)
        kwargs = {
            "only_indices": [0],
            "maximum_distance_factor": 1.02 * cg.distfactor_max / 0.9,
            "max_cn": 13,
            "only_symbols": [mp_symbol],
            "get_from_hints": True,
        }
        # The single-cap hint suggests a CN=12 neighbors set ...
        se = lgf.compute_structure_environments(min_cn=12, **kwargs)
        assert se.neighbors_sets[0][12][0].sources[0]["origin"] == "nb_set_hints"
        assert se.get_csm(0, mp_symbol)["symmetry_measure"] > 0
        # ... which is discarded when it is below min_cn
        se = lgf.compute_structure_environments(min_cn=13, **kwargs)
        assert 13 in se.neighbors_sets[0]
        assert 12 not in se.neighbors_sets[0]

    def test_hints_already_present(self, lgf):
        mp_symbol = "SH:13"
        points = [list(pp) for pp in lgf.allcg[mp_symbol].points]
        # Moving the cap away: the uncapped CN=12 neighbors set is already found from the Voronoi analysis
        points[-1] = [1.1 * cc for cc in points[-1]]
        lgf.allcg = AllCoordinationGeometries(only_symbols=[mp_symbol])
        lgf.setup_test_perfect_environment(mp_symbol, indices="ORDERED", points=points)
        se = lgf.compute_structure_environments(
            only_indices=[0],
            maximum_distance_factor=1.3,
            min_cn=12,
            max_cn=13,
            only_symbols=[mp_symbol],
            get_from_hints=True,
        )
        assert len(se.neighbors_sets[0][12]) == 1
        assert se.neighbors_sets[0][12][0].sources[0]["origin"] == "dist_ang_ac_voronoi"


# the bond valence analysis calls spglib, which emits this itself
@pytest.mark.filterwarnings("ignore:Set OLD_ERROR_HANDLING to false:DeprecationWarning")
class TestComputeCoordinationEnvironments:
    def test_valences(self, lgf, nacl):
        envs = lgf.compute_coordination_environments(nacl, indices=[0], valences=[1] * 4 + [-1] * 4)
        assert envs[0][0]["ce_symbol"] == "O:6"
        assert envs[0][0]["csm"] == approx(0.0, abs=1e-10)
        envs = lgf.compute_coordination_environments(nacl, indices=[4], valences="undefined", only_cations=False)
        assert envs[4][0]["ce_symbol"] == "O:6"
        with pytest.raises(ValueError, match=r"Valences \(2\) do not match the number of sites in the structure \(8\)"):
            lgf.compute_coordination_environments(nacl, valences=[1, -1])

    def test_bond_valence_failure(self, lgf):
        # Valences cannot be assigned to a metal: chemenv falls back to undefined valences
        cu = Structure(Lattice.cubic(3.6), ["Cu"] * 4, [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
        envs = lgf.compute_coordination_environments(cu, indices=[0])
        assert lgf.valences == "undefined"
        assert envs[0][0]["ce_symbol"] == "C:12"
        assert envs[0][0]["csm"] == approx(0.0, abs=1e-10)
