from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pymatgen.analysis.chemenv.coordination_environments.voronoi import DetailedVoronoiContainer
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

__author__ = "waroquiers"

img_files_dir = f"{TEST_FILES_DIR}/analysis/chemenv/images"


class TestVoronoiContainer(MatSciTest):
    def test_voronoi(self):
        # Define a cubic lattice and a list of species (to be used for the fake structures)
        cubic_lattice = Lattice.cubic(10)
        species = ["Cu", "O", "O", "O", "O", "O", "O"]
        valences = "undefined"

        # First fake structure
        coords = [[5, 5, 5]]
        order_and_coords = [
            (1, [4, 5, 5]),
            (2, [6.01, 5, 5]),
            (3, [5, 3.98, 5]),
            (4, [5, 6.03, 5]),
            (5, [5, 5, 3.96]),
            (6, [5, 5, 6.05]),
        ]
        rng = np.random.default_rng()
        rng.shuffle(order_and_coords)
        arr_sorted = np.argsort([oc[0] for oc in order_and_coords]) + 1
        coords.extend([oc[1] for oc in order_and_coords])
        fake_structure = Structure(cubic_lattice, species, coords, coords_are_cartesian=True)

        # First fake structure with a given normalized_distance_tolerance of 0.0100001
        detailed_voronoi_container = DetailedVoronoiContainer(
            structure=fake_structure,
            valences=valences,
            normalized_distance_tolerance=0.0100001,
            isites=[0],
        )
        assert len(detailed_voronoi_container.voronoi_list2[0]) == 6
        neighbors = detailed_voronoi_container.neighbors(0, 1, 0.5)
        assert len(neighbors) == 6
        neighbors = detailed_voronoi_container.neighbors(0, 1.02, 0.5)
        assert len(neighbors) == 6
        neighbors = detailed_voronoi_container.neighbors(0, 1.026, 0.5)
        assert len(neighbors) == 6
        neighbors = detailed_voronoi_container.neighbors(0, 1.5, 0.5)
        assert len(neighbors) == 6

        # First fake structure with a given normalized_distance_tolerance of 0.001
        detailed_voronoi_container = DetailedVoronoiContainer(
            structure=fake_structure,
            valences=valences,
            normalized_distance_tolerance=0.001,
            isites=[0],
        )
        assert len(detailed_voronoi_container.voronoi_list2[0]) == 6
        neighbors = detailed_voronoi_container.neighbors(0, 1, 0.5)
        assert len(neighbors) == 1
        assert neighbors[0]["site"] == fake_structure[arr_sorted[0]]
        neighbors = detailed_voronoi_container.neighbors(0, 1.02, 0.5)
        nbs = [nb["site"] for nb in neighbors]
        assert len(neighbors) == 3
        assert fake_structure[arr_sorted[0]] in nbs
        assert fake_structure[arr_sorted[1]] in nbs
        assert fake_structure[arr_sorted[2]] in nbs
        neighbors = detailed_voronoi_container.neighbors(0, 1.026, 0.5)
        nbs = [nb["site"] for nb in neighbors]
        assert len(neighbors) == 3
        assert fake_structure[arr_sorted[0]] in nbs
        assert fake_structure[arr_sorted[1]] in nbs
        assert fake_structure[arr_sorted[2]] in nbs
        neighbors = detailed_voronoi_container.neighbors(0, 1.5, 0.5)
        assert len(neighbors) == 6

        # Second fake structure
        coords2 = [[5, 5, 5]]
        order_and_coords = [
            (1, [4, 5, 5]),
            (2, [6.01, 5, 5]),
            (3, [5, 3.98, 5]),
            (4, [5, 6.07, 5]),
            (5, [5, 5, 3.92]),
            (6, [5, 5, 6.09]),
        ]
        rng.shuffle(order_and_coords)
        arr_sorted = np.argsort([oc[0] for oc in order_and_coords]) + 1
        coords2.extend([oc[1] for oc in order_and_coords])
        fake_structure2 = Structure(cubic_lattice, species, coords2, coords_are_cartesian=True)

        # Second fake structure with a given normalized_distance_tolerance of 0.0100001
        detailed_voronoi_container = DetailedVoronoiContainer(
            structure=fake_structure2,
            valences=valences,
            normalized_distance_tolerance=0.0100001,
            isites=[0],
        )
        assert len(detailed_voronoi_container.voronoi_list2[0]) == 6
        neighbors = detailed_voronoi_container.neighbors(0, 1, 0.5)
        nbs = [nb["site"] for nb in neighbors]
        assert len(neighbors) == 3
        assert fake_structure2[arr_sorted[0]] in nbs
        assert fake_structure2[arr_sorted[1]] in nbs
        assert fake_structure2[arr_sorted[2]] in nbs
        neighbors = detailed_voronoi_container.neighbors(0, 1.02, 0.5)
        nbs = [nb["site"] for nb in neighbors]
        assert len(neighbors) == 3
        assert fake_structure2[arr_sorted[0]] in nbs
        assert fake_structure2[arr_sorted[1]] in nbs
        assert fake_structure2[arr_sorted[2]] in nbs
        neighbors = detailed_voronoi_container.neighbors(0, 1.026, 0.5)
        nbs = [nb["site"] for nb in neighbors]
        assert len(neighbors) == 3
        assert fake_structure2[arr_sorted[0]] in nbs
        assert fake_structure2[arr_sorted[1]] in nbs
        assert fake_structure2[arr_sorted[2]] in nbs
        neighbors = detailed_voronoi_container.neighbors(0, 1.5, 0.5)
        assert len(neighbors) == 6

        species = ["Cu", "Cu", "O", "O", "O", "Cu", "O"]
        valences = [2, 2, -2, -2, -2, 2, -2]

        # Third fake structure (test of the only_anion_cation_bonds)
        coords = [
            [5, 5, 5],
            [6.01, 5, 5],
            [5, 5, 3.96],
            [4, 5, 5],
            [5, 6.03, 5],
            [5, 3.98, 5],
            [5, 5, 6.05],
        ]
        fake_structure3 = Structure(cubic_lattice, species, coords, coords_are_cartesian=True)
        detailed_voronoi_container = DetailedVoronoiContainer(
            structure=fake_structure3,
            valences=valences,
            normalized_distance_tolerance=0.0100001,
            isites=[0],
            additional_conditions=[DetailedVoronoiContainer.AC.ONLY_ACB],
        )
        assert len(detailed_voronoi_container.voronoi_list2[0]) == 6
        neighbors = detailed_voronoi_container.neighbors(0, 1.01, 0.5)
        nbs = [nb["site"] for nb in neighbors]
        assert len(neighbors) == 6
        assert fake_structure3[1] in nbs
        assert fake_structure3[2] in nbs
        assert fake_structure3[3] in nbs
        assert fake_structure3[4] in nbs
        assert fake_structure3[5] in nbs
        assert fake_structure3[6] in nbs

        # Test of the as_dict() and from_dict() methods as well as __eq__ method
        other_detailed_voronoi_container = DetailedVoronoiContainer.from_dict(detailed_voronoi_container.as_dict())
        assert detailed_voronoi_container == other_detailed_voronoi_container
        assert detailed_voronoi_container.is_close_to(other_detailed_voronoi_container)

    def test_get_vertices_dist_ang_indices(self):
        cubic_lattice = Lattice.cubic(10)
        species = ["Cu", "O", "O", "O", "O", "O", "O"]
        valences = "undefined"

        # First fake structure
        coords = [
            [5, 5, 5],
            [6.01, 5, 5],
            [5, 5, 3.96],
            [4, 5, 5],
            [5, 6.03, 5],
            [5, 3.98, 5],
            [5, 5, 6.05],
        ]
        fake_structure = Structure(cubic_lattice, species, coords, coords_are_cartesian=True)

        # First fake structure with a given normalized_distance_tolerance of 0.0100001
        detailed_voronoi_container = DetailedVoronoiContainer(
            structure=fake_structure,
            valences=valences,
            normalized_distance_tolerance=0.0100001,
            isites=[0],
        )
        fake_parameter_indices = [
            *[(ii, jj) for ii in range(2, 5) for jj in range(7, 14)],
            *[(ii, jj) for ii in range(5, 7) for jj in range(10, 14)],
        ]

        points = detailed_voronoi_container._get_vertices_dist_ang_indices(fake_parameter_indices)
        assert points[0] == (2, 7)
        assert points[1] == (4, 7)
        assert points[2] == (4, 10)
        assert points[3] == (6, 10)
        assert points[4] == (6, 13)
        assert points[5] == (2, 13)


mpl.use("Agg")

CSCL_A = 4.12
# CsCl (B2): 8 Cl nearest neighbors at a*sqrt(3)/2 and 6 Cs next-nearest neighbors at a.
# Both shells share a face with the Voronoi cell (truncated octahedron).
NN_DIST = CSCL_A * np.sqrt(3) / 2
NNN_NORM_DIST = 2 / np.sqrt(3)
# Ratio of the solid angles of the square (Cs) and hexagonal (Cl) faces of the truncated octahedron.
NNN_NORM_ANG = 0.3601016153951434


@pytest.fixture(scope="module")
def cscl():
    return Structure.from_spacegroup("Pm-3m", Lattice.cubic(CSCL_A), ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])


@pytest.fixture(scope="module")
def dvc(cscl):
    return DetailedVoronoiContainer(structure=cscl, valences=[1, -1], isites=[0])


@pytest.fixture(scope="module")
def dvc_limited(cscl):
    return DetailedVoronoiContainer(
        structure=cscl, valences=[1, -1], isites=[0], maximum_distance_factor=1.1, minimum_angle_factor=0.9
    )


def test_voronoi_list(dvc):
    assert dvc.voronoi_list2[1] is None
    nbs = dvc.voronoi_list2[0]
    assert len(nbs) == 14
    cl_nbs = [nb for nb in nbs if nb["site"].specie.symbol == "Cl"]
    cs_nbs = [nb for nb in nbs if nb["site"].specie.symbol == "Cs"]
    assert len(cl_nbs) == 8
    assert len(cs_nbs) == 6
    assert {nb["index"] for nb in cl_nbs} == {1}
    assert {nb["index"] for nb in cs_nbs} == {0}
    for nb in cl_nbs:
        assert nb["distance"] == pytest.approx(NN_DIST)
        assert nb["normalized_distance"] == pytest.approx(1)
        assert nb["normalized_angle"] == pytest.approx(1)
    for nb in cs_nbs:
        assert nb["distance"] == pytest.approx(CSCL_A)
        assert nb["normalized_distance"] == pytest.approx(NNN_NORM_DIST)
        assert nb["normalized_angle"] == pytest.approx(NNN_NORM_ANG)
    # Full solid angle is 4 pi
    assert sum(nb["angle"] for nb in nbs) == pytest.approx(4 * np.pi)


def test_distance_and_angle_groups(dvc):
    nd_groups = dvc.neighbors_normalized_distances[0]
    assert np.array([(grp["min"], grp["max"]) for grp in nd_groups], dtype=float) == pytest.approx(
        np.array([(1, 1), (NNN_NORM_DIST, NNN_NORM_DIST)])
    )
    assert [len(grp["nb_indices"]) for grp in nd_groups] == [8, 14]
    assert [len(grp["dnb_indices"]) for grp in nd_groups] == [8, 6]
    assert nd_groups[0]["next"] == pytest.approx(NNN_NORM_DIST)
    # Tail "next" is the default Voronoi cutoff relative to the smallest distance
    assert nd_groups[1]["next"] == pytest.approx(DetailedVoronoiContainer.default_voronoi_cutoff / NN_DIST)

    d_groups = dvc.neighbors_distances[0]
    assert [grp["min"] for grp in d_groups] == pytest.approx([NN_DIST, CSCL_A])
    assert d_groups[0]["next"] == pytest.approx(CSCL_A)
    assert d_groups[1]["next"] == pytest.approx(DetailedVoronoiContainer.default_voronoi_cutoff)

    na_groups = dvc.neighbors_normalized_angles[0]
    assert np.array([(grp["max"], grp["min"]) for grp in na_groups], dtype=float) == pytest.approx(
        np.array([(1, 1), (NNN_NORM_ANG, NNN_NORM_ANG)])
    )
    assert [len(grp["nb_indices"]) for grp in na_groups] == [8, 14]
    assert [len(grp["dnb_indices"]) for grp in na_groups] == [8, 6]
    assert na_groups[0]["next"] == pytest.approx(NNN_NORM_ANG)
    assert na_groups[1]["next"] == 0

    a_groups = dvc.neighbors_angles[0]
    assert a_groups[1]["max"] / a_groups[0]["max"] == pytest.approx(NNN_NORM_ANG)
    assert a_groups[0]["next"] == pytest.approx(a_groups[1]["max"])
    assert a_groups[1]["next"] == 0


def test_distance_and_angle_groups_with_limits(dvc_limited):
    # Second shell lies beyond maximum_distance_factor and below minimum_angle_factor
    assert len(dvc_limited.voronoi_list2[0]) == 14
    (nd_group,) = dvc_limited.neighbors_normalized_distances[0]
    assert (nd_group["min"], nd_group["max"], nd_group["next"]) == pytest.approx((1, 1, 1.1))
    assert len(nd_group["nb_indices"]) == 8
    (d_group,) = dvc_limited.neighbors_distances[0]
    assert d_group["next"] == pytest.approx(1.1 * NN_DIST)

    (na_group,) = dvc_limited.neighbors_normalized_angles[0]
    assert (na_group["max"], na_group["min"], na_group["next"]) == pytest.approx((1, 1, 0.9))
    assert len(na_group["nb_indices"]) == 8
    (a_group,) = dvc_limited.neighbors_angles[0]
    assert a_group["next"] == pytest.approx(0.9 * a_group["max"])


@pytest.mark.parametrize(
    ("distfactor", "angfactor", "n_cl", "n_cs"),
    [(1.0, 0.5, 8, 0), (1.1, 0.3, 8, 0), (1.2, 0.5, 8, 0), (1.2, 0.3, 8, 6), (2.0, 0.0, 8, 6)],
)
def test_neighbors(dvc, distfactor, angfactor, n_cl, n_cs):
    nbs = dvc.neighbors(0, distfactor, angfactor)
    symbols = [nb["site"].specie.symbol for nb in nbs]
    assert symbols.count("Cl") == n_cl
    assert symbols.count("Cs") == n_cs


@pytest.mark.parametrize(("distfactor", "angfactor"), [(0.9, 0.5), (1.5, 1.1)])
def test_neighbors_out_of_range(dvc, distfactor, angfactor):
    with pytest.raises(ValueError, match="Distance or angle parameter not found"):
        dvc.neighbors(0, distfactor, angfactor)


@pytest.mark.parametrize(
    ("plot_type", "distance_bounds", "dist_limits", "angle_bounds"),
    [
        (None, [0, 1 - 1 / NNN_NORM_DIST, 1], [0, 1], [1, 0, 1 - NNN_NORM_ANG]),
        (
            {"distance_parameter": ("initial_normalized", None), "angle_parameter": ("initial_normalized", None)},
            [1, NNN_NORM_DIST, 2],
            [1, 2],
            [0, 1, NNN_NORM_ANG],
        ),
        (
            {"distance_parameter": ("initial_inverse3_opposite", None), "angle_parameter": ("initial_opposite", None)},
            [0, 1 - 1 / NNN_NORM_DIST**3, 1],
            [0, 1],
            [1, 0, 1 - NNN_NORM_ANG],
        ),
    ],
)
def test_voronoi_parameters_bounds_and_limits(dvc, plot_type, distance_bounds, dist_limits, angle_bounds):
    bounds = dvc.voronoi_parameters_bounds_and_limits(0, plot_type, 2.0)
    assert bounds["distance_bounds"] == pytest.approx(distance_bounds)
    assert bounds["distance_limits"] == pytest.approx(dist_limits)
    assert bounds["angle_bounds"] == pytest.approx(angle_bounds)
    assert bounds["angle_limits"] == [0, 1]


@pytest.mark.parametrize(
    ("plot_type", "match"),
    [
        ({"distance_parameter": ("foo", None), "angle_parameter": ("initial_opposite", None)}, "for the distance"),
        ({"distance_parameter": ("initial_normalized", None), "angle_parameter": ("foo", None)}, "for the angle"),
    ],
)
def test_voronoi_parameters_bounds_and_limits_invalid(dvc, plot_type, match):
    with pytest.raises(NotImplementedError, match=match):
        dvc.voronoi_parameters_bounds_and_limits(0, plot_type, 2.0)


def test_site_without_voronoi_returns_none(dvc):
    assert dvc.voronoi_parameters_bounds_and_limits(1, None, 2.0) is None
    assert dvc.neighbors_surfaces(1) is None
    assert dvc.neighbors_surfaces_bounded(1) is None
    assert dvc.maps_and_surfaces(1) is None
    assert dvc.maps_and_surfaces_bounded(1) is None


def test_neighbors_surfaces(dvc):
    surfaces = dvc.neighbors_surfaces(0)
    dist_plateaus = [1 - 1 / NNN_NORM_DIST, 1 / NNN_NORM_DIST]
    ang_plateaus = [1, 1 - NNN_NORM_ANG]
    expected = np.zeros((3, 3))
    expected[:2, :2] = np.outer(dist_plateaus, ang_plateaus)
    assert surfaces == pytest.approx(expected)


@pytest.mark.parametrize(
    ("options", "expected_row"),
    [
        (None, [0.239734, 0.158893, 0]),
        (
            {
                "type": "standard_diamond",
                "distance_bounds": {"lower": 1.2, "upper": 1.8, "delta": 0.1},
                "angle_bounds": {"lower": 0.1, "upper": 0.8, "delta": 0.1},
            },
            [0.119881, 0.07741, 0],
        ),
        (
            {
                "type": "standard_spline",
                "distance_bounds": {"lower": 1.2, "upper": 1.8},
                "angle_bounds": {"lower": 0.1, "upper": 0.8},
                "lower_points": [(1.2, 0.8), (1.5, 0.3), (1.8, 0.1)],
                "upper_points": [(1.2, 0.8), (1.6, 0.6), (1.8, 0.1)],
                "degree": 2,
            },
            [0.171319, 0.125705, 0],
        ),
        (
            # Angle window [0.1, 0.3] lies entirely below the second angle plateau [0.36, 1]
            {
                "type": "standard_elliptic",
                "distance_bounds": {"lower": 1.2, "upper": 1.8},
                "angle_bounds": {"lower": 0.1, "upper": 0.3},
            },
            [0.068496, 0, 0],
        ),
    ],
)
def test_neighbors_surfaces_bounded(dvc, options, expected_row):
    surfaces = dvc.neighbors_surfaces_bounded(0, surface_calculation_options=options)
    assert surfaces.shape == (3, 3)
    # Only the second distance plateau (1.155 -> 2.0) overlaps the [1.2, 1.8] distance window
    assert surfaces[0] == pytest.approx(np.zeros(3))
    assert surfaces[1] == pytest.approx(expected_row, abs=1e-6)
    assert surfaces[2] == pytest.approx(np.zeros(3))


def test_neighbors_surfaces_bounded_invalid_type(dvc):
    options = {
        "type": "foo",
        "distance_bounds": {"lower": 1.2, "upper": 1.8},
        "angle_bounds": {"lower": 0.1, "upper": 0.8},
    }
    with pytest.raises(ValueError, match="Type 'foo' for the surface calculation"):
        dvc.neighbors_surfaces_bounded(0, surface_calculation_options=options)


def test_precompute_conditions(dvc):
    voronoi = dvc.to_bson_voronoi_list2()[0]
    is_cl = [nb["site"].specie.symbol == "Cl" for nb in dvc.voronoi_list2[0]]

    conditions = dvc._precompute_additional_conditions(0, voronoi, dvc.valences)
    assert set(conditions) == {dvc.AC.NONE, dvc.AC.ONLY_ACB}
    assert all(conditions[dvc.AC.NONE])
    assert conditions[dvc.AC.ONLY_ACB] == is_cl

    angle_conditions = dvc._precompute_angle_conditions(0, voronoi)
    assert len(angle_conditions) == 2
    assert angle_conditions[0] == is_cl
    assert all(angle_conditions[1])


def test_as_dict_from_dict(dvc, dvc_limited):
    for container in (dvc, dvc_limited):
        dct = container.as_dict()
        assert dct["bson_nb_voro_list2"][1] is None
        assert len(dct["bson_nb_voro_list2"][0]) == 14
        new = DetailedVoronoiContainer.from_dict(dct)
        assert new == container
        assert new.is_close_to(container)
        assert container.is_close_to(new)
        assert new.maximum_distance_factor == container.maximum_distance_factor
        assert new.minimum_angle_factor == container.minimum_angle_factor
        assert new.neighbors_normalized_distances[0] == container.neighbors_normalized_distances[0]
        assert new.neighbors_normalized_angles[0] == container.neighbors_normalized_angles[0]


def test_is_close_to(cscl, dvc):
    assert dvc.is_close_to(dvc)

    # Different tolerances/valences
    other = DetailedVoronoiContainer(structure=cscl, valences=[1, -1], isites=[0], normalized_angle_tolerance=0.1)
    assert not dvc.is_close_to(other)
    other = DetailedVoronoiContainer(structure=cscl, valences="undefined", isites=[0])
    assert not dvc.is_close_to(other)

    # Different set of computed sites
    other = DetailedVoronoiContainer(structure=cscl, valences=[1, -1], isites=[0, 1])
    assert not dvc.is_close_to(other)
    assert not other.is_close_to(dvc)

    # Slightly strained structure: neighbors match but distances and normalized quantities differ
    strained_struct = cscl.copy()
    strained_struct.apply_strain([0.001, 0, 0])
    strained = DetailedVoronoiContainer(structure=strained_struct, valences=[1, -1], isites=[0])
    assert not dvc.is_close_to(strained)

    # Different number of sites
    other = DetailedVoronoiContainer(structure=cscl * (1, 1, 2), valences=[1, -1], isites=[0])
    assert not dvc.is_close_to(other)


def _perturbed_container(dvc, perturb):
    voronoi_list2 = [None if vl is None else [dict(nb) for nb in vl] for vl in dvc.voronoi_list2]
    perturb(voronoi_list2[0])
    return DetailedVoronoiContainer(
        structure=dvc.structure, voronoi_list2=voronoi_list2, valences=dvc.valences, isites=[0]
    )


@pytest.mark.parametrize("key", ["distance", "angle", "normalized_distance", "normalized_angle"])
def test_is_close_to_numeric_tolerance(dvc, key):
    def perturb(nbs):
        nbs[3][key] += 1e-6

    other = _perturbed_container(dvc, perturb)
    assert not dvc.is_close_to(other)
    assert dvc.is_close_to(other, atol=1e-5)


@pytest.mark.parametrize(
    "perturb",
    [
        pytest.param(lambda nbs: nbs.pop(), id="missing-neighbor"),
        pytest.param(lambda nbs: nbs[0].update(index=1 - nbs[0]["index"]), id="wrong-index"),
    ],
)
def test_is_close_to_neighbor_mismatch(dvc, perturb):
    assert not dvc.is_close_to(_perturbed_container(dvc, perturb))


def test_infinite_vertex():
    # Slab without neighbors below the bottom layer within the cutoff: the Voronoi cell is open
    slab = Structure(Lattice.orthorhombic(3, 3, 30), ["Cu", "Cu"], [[0, 0, 0], [0, 0, 2]], coords_are_cartesian=True)
    with pytest.raises(RuntimeError, match="infinite vertex"):
        DetailedVoronoiContainer(structure=slab, voronoi_cutoff=3.5, isites=[0])


def test_eq(dvc):
    assert dvc != "not a voronoi container"
    assert dvc.__eq__(1) is NotImplemented


@pytest.mark.parametrize(
    ("normalized", "step_function", "x_end", "y_start"),
    [
        (True, {"type": "step_function", "scale": 0.0001}, 1.1 * NNN_NORM_DIST, 0),
        (False, {"type": "step_function", "scale": 0.0001}, 1.1 * CSCL_A, 0),
        # Default normal_cdf maps the first shell (normalized distance 1) to x = 0, i.e. half of its 8 neighbors
        (True, None, 1.1 * (1 - 1 / NNN_NORM_DIST**3), 4),
    ],
)
def test_get_rdf_figure(dvc, normalized, step_function, x_end, y_start):
    fig = dvc.get_rdf_figure(0, normalized=normalized, figsize=(4, 3), step_function=step_function)
    xx, yy = fig.axes[0].lines[0].get_data()
    assert xx[-1] == pytest.approx(x_end)
    assert yy[0] == pytest.approx(y_start)
    assert yy[-1] == pytest.approx(14)
    if step_function is not None:
        # Steps of 8 (Cl shell) then 6 (Cs shell)
        assert list(yy) == [0, 0, 8, 8, 14, 14]
    plt.close(fig)


@pytest.mark.parametrize(
    ("normalized", "step_function"),
    [(True, None), (True, {"type": "normal_cdf", "scale": 0.0001}), (False, None)],
)
def test_get_sadf_figure(dvc, normalized, step_function):
    fig = dvc.get_sadf_figure(0, normalized=normalized, step_function=step_function)
    xx, yy = fig.axes[0].lines[0].get_data()
    if normalized:
        assert xx[-1] == pytest.approx(1.1 * NNN_NORM_ANG**-0.1)
    assert yy[0] == pytest.approx(0)
    assert yy[-1] == pytest.approx(14)
    plt.close(fig)


@pytest.mark.parametrize("method", ["get_rdf_figure", "get_sadf_figure"])
def test_figure_invalid_step_function(dvc, method):
    with pytest.raises(ValueError, match="Step function of type 'foo' is not allowed"):
        getattr(dvc, method)(0, step_function={"type": "foo"})
    plt.close("all")
